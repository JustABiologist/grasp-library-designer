"""One-shot GRASP design: RNA → binder protein → free GGA cuts → oligos.

No combinatorial library modules. The full binder CDS is codon/synthesis
optimized, then cut sites are chosen so the overhangs already present in that
DNA form a high-fidelity Golden Gate set.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, Optional

import pandas as pd
from Bio.Seq import Seq

from .binder import describe_binder, normalize_target_rna
from .dna import clean_dna, translate_dna
from .gga_split import plan_gga_from_optimized_cds, suggest_fragment_count
from .ligation_fidelity import LigationFidelityCalculator
from .optimizer import optimize_coding_sequence, synthesis_qc


def sanitize_rna_name(target_rna: str) -> str:
    return normalize_target_rna(target_rna)


DEFAULT_OLIGO_PREFIX = "ACATCTC"
DEFAULT_OLIGO_SUFFIX = "TTGTCTTC"


def run_oneshot_design(
    *,
    target_rna: str,
    codon_data: Dict,
    config: Dict[str, Any],
    output_dir: Path,
    seed: int = 42,
    n_fragments: Optional[int] = None,
    max_fragment_cds: Optional[int] = None,
    oligo_prefix: str = DEFAULT_OLIGO_PREFIX,
    oligo_suffix: str = DEFAULT_OLIGO_SUFFIX,
    fidelity: Optional[LigationFidelityCalculator] = None,
    log: Callable[[str], None] = print,
) -> Dict[str, Any]:
    """
    Protein-first one-shot pipeline (no library parts):

    1. RNA → PPR code → continuous binder AA
    2. Codon + synthesis optimize the full CDS
    3. Choose codon-aligned cut sites whose DNA overhangs are high-fidelity
    4. Export GGA oligos (flanks + fragment DNA)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    info = describe_binder(target_rna)
    rna = info["target_rna"]
    aa = info["aa_sequence"]
    log(
        f"▶ One-shot · RNA {rna} → {info['aa_length']} aa binder "
        f"(PPR {info['ppr_code']})"
    )

    synth = config.get("synthesis", {})
    if max_fragment_cds is None:
        max_oligo = int(synth.get("max_oligo_length", 300))
        flank = len(oligo_prefix) + len(oligo_suffix)
        max_fragment_cds = max(120, max_oligo - flank - 10)

    log(
        f"  Annealing full CDS ({info['cds_length']} nt, "
        f"{config['optimizer']['iterations_per_part']:,} iters)…"
    )
    full_mask = "N" * (3 * len(aa))
    cds, objective = optimize_coding_sequence(
        aa_sequence=aa,
        coding_mask=full_mask,
        codon_data=codon_data,
        config=config,
    )
    cds = clean_dna(cds)
    if translate_dna(cds, config.get("genetic_code", 1)) != aa:
        raise AssertionError("Full-CDS anneal changed the protein sequence")
    log(f"  Full CDS objective={objective:.3f}")

    lig = config.get("ligation", {})
    calc = fidelity or LigationFidelityCalculator(
        temperature=lig.get("temperature", 25),
        hours=lig.get("hours", 18),
        ligation_table=lig.get("ligation_table"),
        min_efficiency=lig.get("min_efficiency", 0.25),
        min_fidelity=lig.get("min_fidelity", 0.9),
    )

    if n_fragments is None:
        n_fragments = suggest_fragment_count(
            len(cds), max_fragment_cds=max_fragment_cds
        )

    plan = plan_gga_from_optimized_cds(
        cds,
        aa,
        n_fragments=n_fragments,
        max_fragment_cds=max_fragment_cds,
        fidelity=calc,
        seed=seed,
    )
    # Rebuild fragment DNA slices properly (codon-complete upstream;
    # downstream includes 5′ overhang for Type IIS sticky end)
    cuts = [int(x) for x in plan["aa_end_0based"].tolist()[:-1]]
    ohs = [str(x) for x in plan["oh3"].tolist()[:-1]]
    fid = float(plan["ligation_fidelity_set"].iloc[0])
    log(f"  GGA · {len(plan)} fragments · overhang set fidelity {fid:.4f}")
    log(f"  Overhangs: {', '.join(ohs)}")

    rows = []
    n = len(plan)
    for i, prow in enumerate(plan.itertuples(index=False)):
        aa0, aa1 = int(prow.aa_start_0based), int(prow.aa_end_0based)
        if i == 0:
            frag_cds = cds[0 : 3 * aa1]
        else:
            frag_cds = cds[3 * aa0 - 4 : 3 * aa1]
        full_oligo = clean_dna(oligo_prefix) + frag_cds + clean_dna(oligo_suffix)
        cds_qc = synthesis_qc(frag_cds, config)
        oligo_qc = synthesis_qc(full_oligo, config, forbidden_scan=frag_cds)
        rows.append(
            {
                "fragment_id": prow.fragment_id,
                "assembly_order": int(prow.assembly_order),
                "aa_start_0based": aa0,
                "aa_end_0based": aa1,
                "aa_length": aa1 - aa0,
                "oh5": prow.oh5,
                "oh3": prow.oh3,
                "fragment_cds": frag_cds,
                "oligo_prefix": clean_dna(oligo_prefix),
                "oligo_suffix": clean_dna(oligo_suffix),
                "oligo_sequence_5to3": full_oligo,
                "oligo_length": len(full_oligo),
                "oligo_gc": oligo_qc["gc_fraction"],
                "oligo_warnings": oligo_qc["warnings"],
                "oligo_failures": oligo_qc["failures"],
                "qc_passed": cds_qc["passed"] and oligo_qc["passed"],
                "ligation_fidelity_set": fid,
            }
        )
        log(
            f"    {prow.fragment_id}: aa {aa0}–{aa1} · "
            f"{len(frag_cds)} nt CDS · oligo {len(full_oligo)} nt · "
            f"oh5={prow.oh5 or '—'} oh3={prow.oh3 or '—'}"
        )

    oligos = pd.DataFrame(rows)
    plan.to_csv(output_dir / f"gga_plan_{rna}.csv", index=False)
    oligos.to_csv(output_dir / f"oneshot_{rna}_oligos.csv", index=False)
    with open(output_dir / f"full_cds_{rna}.fasta", "w") as handle:
        handle.write(f">GRASP_oneshot_{rna}_full_CDS\n{cds}\n")
    with open(output_dir / f"oneshot_{rna}_oligos.fasta", "w") as handle:
        for row in oligos.itertuples(index=False):
            handle.write(
                f">{rna}|{row.fragment_id}|order{row.assembly_order}\n"
                f"{row.oligo_sequence_5to3}\n"
            )
    with open(output_dir / f"oneshot_{rna}_assembled.fasta", "w") as handle:
        handle.write(f">GRASP_oneshot_{rna}\n{cds}\n")

    # Reassemble from oligos: F1 full ORF slice; Fi drop 4-nt 5′ overhang
    bits = []
    for i, row in enumerate(oligos.itertuples(index=False)):
        bits.append(row.fragment_cds if i == 0 else row.fragment_cds[4:])
    assembled = clean_dna("".join(bits))
    translation_ok = (
        len(assembled) % 3 == 0
        and translate_dna(assembled, config.get("genetic_code", 1)) == aa
        and assembled == cds
    )

    summary = {
        "target_rna": rna,
        "ppr_code": info["ppr_code"],
        "aa_length": info["aa_length"],
        "n_fragments": n,
        "ligation_fidelity": fid,
        "translation_verified": translation_ok,
        "fragments_qc": bool(oligos["qc_passed"].all()),
        "full_cds_objective": objective,
    }
    pd.DataFrame([summary]).to_csv(
        output_dir / f"oneshot_{rna}_summary.csv", index=False
    )

    log(
        f"✓ One-shot done · translation_verified={translation_ok} · "
        f"fragment_qc={summary['fragments_qc']} · "
        f"oligos → oneshot_{rna}_oligos.csv"
    )

    return {
        "target_rna": rna,
        "binder": info,
        "aa_sequence": aa,
        "full_cds": cds,
        "gga_plan": plan,
        "oligos": oligos,
        "assembled": {
            "assembled_cds": assembled,
            "expected_protein": aa,
            "observed_protein": (
                translate_dna(assembled, config.get("genetic_code", 1))
                if len(assembled) % 3 == 0
                else ""
            ),
            "translation_verified": translation_ok,
            "ligation_fidelity": fid,
            "binder_info": info,
        },
        "output_dir": output_dir,
        "oligo_csv": output_dir / f"oneshot_{rna}_oligos.csv",
        "oligo_fasta": output_dir / f"oneshot_{rna}_oligos.fasta",
        "assembled_fasta": output_dir / f"oneshot_{rna}_assembled.fasta",
        "selected_overhangs": {f"J{i+1}": oh for i, oh in enumerate(ohs)},
        "pareto_front": None,
        "summary": summary,
    }


def subset_for_target(*_args, **_kwargs):
    raise RuntimeError(
        "subset_for_target is obsolete. One-shot design no longer uses library "
        "modules — use run_oneshot_design (RNA → protein → free GGA cuts)."
    )
