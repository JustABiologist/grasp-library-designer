"""One target RNA → Golden Gate oligos for the binder gene.

The combinatorial GRASP library and ARELF junctions live on the library path.
One-shot co-designs cut sites, overhangs, and coding sequence for a single
binder, then wraps the fragments as orderable Type IIS oligos.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import numpy as np
import pandas as pd

from .assembly_interfaces import (
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    extract_order_payload,
    resolve_assembly_interfaces,
    reverse_complement,
)
from .binder import describe_binder, normalize_target_rna
from .codon_validation import verify_cds_for_organism
from .dna import clean_dna, gc_fraction, mask_matches
from .gga_split import (
    assembly_plan_frame,
    co_design_gene_assembly,
    oligo_flank_overhead,
    resolve_wrap_enzyme,
    wrap_geometry,
)
from .ligation_fidelity import LigationFidelityCalculator
from .optimizer import synthesis_qc
from .synthesis_vendors import ligation_protocol_for_level


ENTRY_CLONING_ENZYME = "BsaI"
MODULE_RELEASE_ENZYME = "BpiI / BbsI"


def sanitize_rna_name(target_rna: str) -> str:
    return normalize_target_rna(target_rna)


def validate_order_fragment_in_silico(
    sequence: str, interfaces: Dict[str, Any]
) -> Dict[str, Any]:
    """Check order-fragment geometry without claiming wet-lab/vector validation."""
    sequence = clean_dna(sequence)
    payload = extract_order_payload(sequence, interfaces)
    entry = interfaces["level_minus1_entry"]
    entry_insert = (
        entry[FIVE_PRIME_CODING_SITE]
        + payload
        + entry[THREE_PRIME_CODING_SITE]
    )
    context_5 = entry.get("completion_context_5p")
    context_3 = entry.get("completion_context_3p")
    release_site = entry.get("release_recognition_site")
    context_available = bool(context_5 and context_3 and release_site)
    cloned_context = (
        clean_dna(context_5) + entry_insert + clean_dna(context_3)
        if context_available
        else ""
    )
    context_validated = bool(
        context_available
        and cloned_context.startswith(clean_dna(release_site))
        and cloned_context.endswith(reverse_complement(clean_dna(release_site)))
    )
    if context_available and not context_validated:
        raise ValueError("configured entry-vector context does not complete release sites")
    return {
        "assembly_interface_profile": interfaces["profile_name"],
        "terminal_site_notation": interfaces["notation"],
        "coding_strand_direction": interfaces["coding_strand_direction"],
        "entry_insert_5to3": entry_insert,
        "entry_five_prime_end_overhang": entry[FIVE_PRIME_END],
        "entry_three_prime_end_overhang": entry[THREE_PRIME_END],
        "entry_five_prime_assembled_coding_site": entry_insert[:4],
        "entry_three_prime_assembled_coding_site": entry_insert[-4:],
        "cloned_entry_context_5to3": cloned_context,
        "module_release_payload_5to3": payload,
        "module_release_five_prime_assembled_coding_site": payload[:4],
        "module_release_three_prime_assembled_coding_site": payload[-4:],
        "order_fragment_requirements_checked": True,
        "entry_interface_requirements_checked": True,
        "entry_vector_context_in_silico_validated": context_validated,
        "entry_vector_sequence_provided": bool(entry.get("vector_sequence")),
        "entry_vector_sequence_in_silico_validated": False,
        "module_release_requirements_checked": True,
    }


def validate_pagm1311_order_fragment(sequence: str) -> Dict[str, Any]:
    """Backward-compatible deposited-preset validator.

    New exports use the explicit ``*_requirements_checked`` and
    ``*_in_silico_validated`` fields returned by the generic validator.
    """
    result = validate_order_fragment_in_silico(
        sequence, resolve_assembly_interfaces(preset="deposited_grasp")
    )
    # Preserve historical sequence-field spellings for callers, not claims.
    result["entry_insertion_overhang_5"] = result[
        "entry_five_prime_assembled_coding_site"
    ]
    result["entry_insertion_overhang_3"] = result[
        "entry_three_prime_assembled_coding_site"
    ]
    result["cloned_pagm1311_context_5to3"] = result["cloned_entry_context_5to3"]
    result["bpii_release_payload_5to3"] = result["module_release_payload_5to3"]
    result["bpii_release_oh5"] = result[
        "module_release_five_prime_assembled_coding_site"
    ]
    result["bpii_release_oh3"] = result[
        "module_release_three_prime_assembled_coding_site"
    ]
    return result


def _fidelity_calculator(
    config: Dict[str, Any],
    wrap_enzyme: str,
    fidelity: Optional[LigationFidelityCalculator],
) -> LigationFidelityCalculator:
    if fidelity is not None:
        return fidelity
    geometry = wrap_geometry(wrap_enzyme)
    lig = dict(config.get("ligation") or {})
    table_name = str(lig.get("table_name") or "")
    potapov = "potapov" in table_name.lower() or "t4 ligase" in table_name.lower()
    if potapov:
        return LigationFidelityCalculator(
            temperature=lig.get("temperature"),
            hours=lig.get("hours"),
            ligation_table=lig.get("ligation_table"),
        )
    protocol = ligation_protocol_for_level(
        config, str(geometry.get("fidelity_level") or "level_minus1")
    )
    return LigationFidelityCalculator(
        temperature=protocol.get("temperature"),
        hours=protocol.get("hours"),
        ligation_table=protocol.get("ligation_table"),
    )


def _destination_overhangs(config: Dict[str, Any]) -> tuple[str, str]:
    oneshot = dict(config.get("oneshot") or {})
    dest5 = oneshot.get("destination_5prime_overhang")
    dest3 = oneshot.get("destination_3prime_overhang")
    if dest5 and dest3:
        return str(dest5), str(dest3)
    try:
        interfaces = resolve_assembly_interfaces(config)
        outer = interfaces["level0"]["acceptor_outer"]
        return outer[FIVE_PRIME_END], outer[THREE_PRIME_END]
    except Exception:
        return "CTCA", "CTCG"


def _oligo_table(
    design,
    *,
    target_rna: str,
    config: Dict[str, Any],
    codon_data: Dict,
    plan: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    aa = design.aa_sequence
    translation = verify_cds_for_organism(
        design.cds,
        aa,
        genetic_code=int(config.get("genetic_code", 1)),
        codon_data=codon_data,
    )
    cds_qc = synthesis_qc(
        design.cds,
        config,
        sequence_kind="cds",
        rescue_notes=design.rescue_notes,
    )
    for order, (payload, oligo) in enumerate(
        zip(design.payloads, design.oligos), start=1
    ):
        frag_id = f"F{order}"
        plan_row = plan.loc[plan["fragment_id"] == frag_id].iloc[0]
        oligo_qc = synthesis_qc(
            oligo,
            config,
            forbidden_scan=design.cds,
            sequence_kind="oligo",
            rescue_notes=design.rescue_notes if order == 1 else (),
        )
        rows.append(
            {
                "order_fragment_id": f"{target_rna}_{frag_id}",
                "fragment_id": frag_id,
                "part_id": frag_id,
                "optimized_part_id": f"{target_rna}_{frag_id}",
                "assembly_order": order,
                "order_quantity": 1,
                "aa_start_0based": int(plan_row.aa_start_0based),
                "aa_end_0based": int(plan_row.aa_end_0based),
                "aa_sequence": str(plan_row.aa_sequence),
                "optimized_cds": design.cds[
                    int(plan_row.cds_start) : int(plan_row.cds_end)
                ],
                "payload_5to3": payload,
                "oligo_sequence_5to3": oligo,
                "order_sequence_5to3": oligo,
                "oh5_coding_site_5to3": str(plan_row.oh5_coding_site_5to3),
                "oh3_coding_site_5to3": str(plan_row.oh3_coding_site_5to3),
                "five_prime_end_overhang": str(plan_row.five_prime_end_overhang),
                "three_prime_end_overhang": str(plan_row.three_prime_end_overhang),
                "wrap_enzyme": design.wrap_enzyme,
                "oligo_length": len(oligo),
                "oligo_gc": gc_fraction(oligo),
                "oligo_gc_pct": oligo_qc["gc_pct"],
                "oligo_gc_status": oligo_qc["gc_status"],
                "oligo_homopolymers": oligo_qc["homopolymer_detail"],
                "oligo_warnings": oligo_qc["warnings"],
                "oligo_failures": oligo_qc["failures"],
                "codon_score": design.codon_score,
                "ligation_fidelity_set": design.ligation_fidelity,
                "translation_verified": bool(translation["ok"]),
                "genetic_code_ok": bool(translation["genetic_code_ok"]),
                "codon_table_ok": bool(translation["codon_table_ok"]),
                "genetic_code": int(config.get("genetic_code", 1)),
                "mask_verified": mask_matches(design.cds, design.coding_mask),
                "cds_gc_pct": cds_qc["gc_pct"],
                "cds_gc_status": cds_qc["gc_status"],
                "cds_local_gc": cds_qc["local_gc_detail"],
                "cds_longest_homopolymer": cds_qc["longest_homopolymer"],
                "cds_homopolymers": cds_qc["homopolymer_detail"],
                "cds_repeats": cds_qc["repeat_detail"],
                "cds_warnings": cds_qc["warnings"],
                "cds_failures": cds_qc["failures"],
                "blacklist_tested": cds_qc["blacklist_tested"],
                "blacklist_hits": cds_qc["blacklist_hits"],
                "rescue_codon_count": cds_qc["rescue_codon_count"],
                "rescue_codon_detail": cds_qc["rescue_codon_detail"],
                "qc_status": (
                    "FAIL"
                    if (
                        not cds_qc["hard_constraints_passed"]
                        or not oligo_qc["hard_constraints_passed"]
                        or not bool(translation["ok"])
                    )
                    else "WARNING"
                    if (cds_qc["warnings"] or oligo_qc["warnings"])
                    else "PASS"
                ),
                "hard_constraints_passed": (
                    cds_qc["hard_constraints_passed"]
                    and oligo_qc["hard_constraints_passed"]
                    and bool(translation["ok"])
                ),
                "qc_passed": (
                    cds_qc["passed"]
                    and oligo_qc["passed"]
                    and bool(translation["ok"])
                ),
                "vendor_acceptance_confirmed": False,
                "sequence_type": "double-stranded DNA synthesis fragment",
            }
        )
    return pd.DataFrame(rows)


def run_oneshot_design(
    *,
    target_rna: str,
    codon_data: Dict,
    config: Dict[str, Any],
    output_dir: Path,
    seed: int = 42,
    n_fragments: Optional[int] = None,
    max_fragment_cds: Optional[int] = None,
    oligo_prefix: Optional[str] = None,
    oligo_suffix: Optional[str] = None,
    fidelity=None,
    log: Callable[[str], None] = print,
) -> Dict[str, Any]:
    """Design Golden Gate oligos that assemble into one binder gene.

    Input is a target RNA. Output is a set of Type IIS oligos whose payloads
    ligate to the codon-optimized binder CDS. Cut positions, junction overhangs,
    and synonymous coding sequence are chosen together.
    """
    if oligo_prefix is not None or oligo_suffix is not None:
        raise ValueError(
            "Custom one-shot flanks are chosen via wrap_enzyme and "
            "destination overhangs, not oligo_prefix/oligo_suffix"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    info = describe_binder(target_rna)
    rna = info["target_rna"]
    oneshot_cfg = dict(config.get("oneshot") or {})
    wrap_enzyme = resolve_wrap_enzyme(
        oneshot_cfg.get("wrap_enzyme") or config.get("wrap_enzyme") or "BsaI"
    )
    dest5, dest3 = _destination_overhangs(config)
    requested_fragments = n_fragments if n_fragments is not None else oneshot_cfg.get(
        "n_fragments"
    )
    if requested_fragments in (0, "auto", "Auto"):
        requested_fragments = None

    local_config = dict(config)
    forbidden = dict(local_config.get("forbidden_sites") or {})
    site = str(wrap_geometry(wrap_enzyme)["recognition_site"])
    forbidden.setdefault(wrap_enzyme, site)
    local_config["forbidden_sites"] = forbidden

    calc = _fidelity_calculator(local_config, wrap_enzyme, fidelity)
    max_oligo = oneshot_cfg.get("max_oligo_length")
    if max_fragment_cds is not None:
        max_oligo = int(max_fragment_cds) + oligo_flank_overhead(
            wrap_geometry(wrap_enzyme)
        )
    old_py_state = random.getstate()
    old_np_state = np.random.get_state()
    random.seed(seed)
    np.random.seed(seed)
    try:
        design = co_design_gene_assembly(
            info["aa_sequence"],
            codon_data,
            local_config,
            n_fragments=None if requested_fragments is None else int(requested_fragments),
            destination_5prime=dest5,
            destination_3prime=dest3,
            wrap_enzyme=wrap_enzyme,
            fidelity=calc,
            max_oligo_length=max_oligo,
        )
    finally:
        random.setstate(old_py_state)
        np.random.set_state(old_np_state)

    plan = assembly_plan_frame(design)
    oligos = _oligo_table(
        design,
        target_rna=rna,
        config=local_config,
        codon_data=codon_data,
        plan=plan,
    )
    translation = verify_cds_for_organism(
        design.cds,
        info["aa_sequence"],
        genetic_code=int(local_config.get("genetic_code", 1)),
        codon_data=codon_data,
    )
    if not translation.get("ok"):
        raise AssertionError("Designed CDS does not translate to the binder protein")

    assembled = {
        "assembled_cds": design.cds,
        "expected_protein": info["aa_sequence"],
        "observed_protein": info["aa_sequence"],
        "translation_verified": True,
        "genetic_code_ok": bool(translation["genetic_code_ok"]),
        "codon_table_ok": bool(translation["codon_table_ok"]),
        "genetic_code": int(local_config.get("genetic_code", 1)),
        "n_fragments": len(design.oligos),
        "ligation_fidelity": design.ligation_fidelity,
        "codon_score": design.codon_score,
    }

    plan_path = output_dir / f"assembly_plan_{rna}.csv"
    order_csv = output_dir / f"oneshot_{rna}_orderable_fragments.csv"
    legacy_csv = output_dir / f"oneshot_{rna}_oligos.csv"
    order_fasta = output_dir / f"oneshot_{rna}_orderable_fragments.fasta"
    legacy_fasta = output_dir / f"oneshot_{rna}_oligos.fasta"
    gene_fasta = output_dir / f"oneshot_{rna}_gene.fasta"

    plan.to_csv(plan_path, index=False)
    oligos.to_csv(order_csv, index=False)
    oligos.to_csv(legacy_csv, index=False)
    fasta_text = "".join(
        f">{row.order_fragment_id}|{row.wrap_enzyme}|{row.oh5_coding_site_5to3}.."
        f"{row.oh3_coding_site_5to3}|{row.oligo_length}bp\n{row.order_sequence_5to3}\n"
        for row in oligos.itertuples(index=False)
    )
    order_fasta.write_text(fasta_text)
    legacy_fasta.write_text(fasta_text)
    gene_fasta.write_text(
        f">GRASP_{rna}|binder_CDS|{info['aa_length']}aa|{len(design.cds)}nt\n"
        f"{design.cds}\n"
    )

    summary = {
        "target_rna": rna,
        "architecture": f"{len(rna)}S",
        "ppr_code": info["ppr_code"],
        "aa_length": info["aa_length"],
        "cds_length": len(design.cds),
        "n_fragments": len(design.oligos),
        "n_unique_order_fragments": len(oligos),
        "wrap_enzyme": design.wrap_enzyme,
        "destination_5prime_overhang": design.destination_5prime,
        "destination_3prime_overhang": design.destination_3prime,
        "destination_3prime_coding_site": design.destination_3prime_coding,
        "junction_overhangs": ";".join(design.overhangs),
        "cut_aa": ";".join(str(cut) for cut in design.cuts),
        "ligation_fidelity": design.ligation_fidelity,
        "codon_score": design.codon_score,
        "translation_verified": True,
        "mask_verified": bool(oligos["mask_verified"].all()),
        "fragments_clean_qc": bool(oligos["qc_passed"].all()),
        "fragments_hard_constraints_passed": bool(
            oligos["hard_constraints_passed"].all()
        ),
        "vendor_acceptance_confirmed": False,
        "organism": local_config.get("selected_organism")
        or local_config.get("selected_organism_label"),
        "standalone_expression_cassette": False,
    }
    pd.DataFrame([summary]).to_csv(
        output_dir / f"oneshot_{rna}_summary.csv", index=False
    )

    log(
        f"✓ One-shot {len(rna)}S: {len(oligos)} {design.wrap_enzyme} oligos; "
        f"fidelity {design.ligation_fidelity:.4f}; translation checked"
    )
    return {
        "target_rna": rna,
        "binder": info,
        "aa_sequence": info["aa_sequence"],
        "cds": design.cds,
        "binding_tract_cds": design.cds,
        "design": design,
        "gga_plan": plan,
        "assembly_plan": plan,
        "oligos": oligos,
        "orderable_fragments": oligos,
        "assembled": assembled,
        "output_dir": output_dir,
        "oligo_csv": legacy_csv,
        "oligo_fasta": legacy_fasta,
        "order_csv": order_csv,
        "order_fasta": order_fasta,
        "gene_fasta": gene_fasta,
        "selected_overhangs": {
            f"J{i + 1}": oh for i, oh in enumerate(design.overhangs)
        },
        "summary": summary,
    }


def subset_for_target(*_args, **_kwargs):
    raise RuntimeError(
        "subset_for_target is obsolete; run_oneshot_design now co-designs "
        "Golden Gate oligos for the continuous binder gene"
    )
