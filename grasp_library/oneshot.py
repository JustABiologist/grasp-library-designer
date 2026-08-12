"""Target-specific GRASP order fragments and in-silico assembly checks."""

from __future__ import annotations

import copy
import random
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import numpy as np
import pandas as pd

from .binder import describe_binder, normalize_target_rna
from .arelf import materialize_arelf_parts
from .assembly_interfaces import (
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    extract_order_payload,
    resolve_assembly_interfaces,
    reverse_complement,
)
from .dna import clean_dna
from .import_grasp import (
    build_parts_table,
    compile_target_gap,
    load_grasp_records,
)
from .optimizer import optimize_library, simulate_assembled_cds
from .paths import bundled_profile_genbank


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
        "module_release_oh5": payload[:4],
        "module_release_oh3": payload[-4:],
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
    result["bpii_release_oh5"] = result["module_release_oh5"]
    result["bpii_release_oh3"] = result["module_release_oh3"]
    return result


def _validate_level0_groups(
    plan: pd.DataFrame,
    order_rows: pd.DataFrame,
    interfaces: Dict[str, Any],
) -> pd.DataFrame:
    """Join every five-part block and check configured Level 0 requirements."""
    by_id = order_rows.set_index("optimized_part_id")
    rows = []
    for group, group_plan in plan.groupby("assembly_group", sort=False):
        group_plan = group_plan.sort_values("assembly_order")
        payloads = []
        part_ids = []
        for row in group_plan.itertuples(index=False):
            data = by_id.loc[str(row.optimized_part_id)]
            payloads.append(str(data.module_release_payload_5to3))
            part_ids.append(str(row.part_id))

        if len(payloads) != 5:
            raise AssertionError(f"{group}: expected five Level -1 modules")
        outer = interfaces["level0"].get("acceptor_outer")
        if outer is not None and payloads[0][:4] != outer[FIVE_PRIME_CODING_SITE]:
            raise AssertionError(f"{group}: first module misses the Level 0 outer interface")
        if outer is not None and payloads[-1][-4:] != outer[THREE_PRIME_CODING_SITE]:
            raise AssertionError(f"{group}: last module misses the Level 0 outer interface")
        for left, right in zip(payloads, payloads[1:]):
            if left[-4:] != right[:4]:
                raise AssertionError(
                    f"{group}: BpiI junction mismatch {left[-4:]} != {right[:4]}"
                )

        assembled = payloads[0] + "".join(payload[4:] for payload in payloads[1:])
        block_insert = assembled[4:-4] if outer is not None else assembled
        rows.append(
            {
                "assembly_group": group,
                "level_minus1_parts": ";".join(part_ids),
                "n_parts": len(part_ids),
                "level0_outer_left_5to3": assembled[:4] if outer else "",
                "level0_outer_right_5to3": assembled[-4:] if outer else "",
                "ppr_block_5to3": block_insert,
                "ppr_block_length": len(block_insert),
                "level0_interface_requirements_checked": True,
                "level0_module_chain_in_silico_validated": True,
                "level0_vector_sequence_provided": bool(
                    interfaces["level0"].get("vector_sequence")
                ),
                "level0_vector_sequence_in_silico_validated": False,
            }
        )
    return pd.DataFrame(rows)


def _validate_ppr_block_chain(
    plan: pd.DataFrame,
    level0: pd.DataFrame,
    library: pd.DataFrame,
    parts_full: pd.DataFrame,
    assembled_cds: str,
    interfaces_profile: Dict[str, Any],
    architecture: str,
) -> Dict[str, Any]:
    """Join the PPR blocks and validate only that in-silico block chain.

    The physical Level 1 module begins at the first external coding overhang,
    whereas the in-frame translation model also needs the few leading bases
    that precede that overhang in the deposited first-part coding window.  The
    comparison below proves that the exported two-stage assembly reconstructs
    the same binding tract used for translation QC.
    """
    group_order = (
        plan.groupby("assembly_group", sort=False)["assembly_slot"]
        .min()
        .sort_values()
        .index.tolist()
    )
    by_group = level0.set_index("assembly_group")
    if set(group_order) != set(by_group.index):
        raise AssertionError("Level 0 products do not match the assembly plan")
    layout = interfaces_profile["architectures"][architecture]
    if group_order != layout["blocks"]:
        raise AssertionError(
            f"{architecture}: block order {group_order} does not match {layout['blocks']}"
        )
    blocks = [str(by_group.loc[group, "ppr_block_5to3"]) for group in group_order]
    ppr_outer = interfaces_profile["level0"]["ppr_outer"]
    if blocks[0][:4] != ppr_outer[FIVE_PRIME_CODING_SITE]:
        raise AssertionError("first PPR block misses the configured PPR N interface")
    if blocks[-1][-4:] != ppr_outer[THREE_PRIME_CODING_SITE]:
        raise AssertionError("last PPR block misses the configured PPR C interface")
    for left, right, join_name in zip(blocks, blocks[1:], layout["joins"]):
        junction = interfaces_profile["junctions"][join_name]
        upstream = junction["upstream_three_prime_end_overhang"]
        downstream = junction["downstream_five_prime_end_overhang"]
        if reverse_complement(upstream) != downstream:
            raise AssertionError(f"{join_name}: physical sticky ends are incompatible")
        expected = junction["assembled_coding_site"]
        if left[-4:] != expected or right[:4] != expected:
            raise AssertionError(
                f"{join_name}: PPR block junction mismatch "
                f"{left[-4:]}/{right[:4]} != {expected}"
            )
    module_core = blocks[0] + "".join(block[4:] for block in blocks[1:])

    first_plan = plan.sort_values("assembly_slot").iloc[0]
    first_part = str(first_plan["part_id"])
    first_optimized = str(
        library.set_index("optimized_part_id").loc[
            str(first_plan["optimized_part_id"]), "optimized_cds"
        ]
    )
    oh5 = int(parts_full.set_index("part_id").loc[first_part, "oh5_mask_start"])
    coding_frame_prefix = clean_dna(first_optimized)[:oh5]
    reconstructed = coding_frame_prefix + module_core
    if reconstructed != clean_dna(assembled_cds):
        raise AssertionError(
            "Exported pAGM9121 blocks do not reconstruct the translated GRASP tract"
        )

    interface_chain = [blocks[0][:4]] + [block[-4:] for block in blocks]
    final = interfaces_profile["final_cassette"]
    return {
        "level0_group_order": ";".join(group_order),
        "ppr_interface_chain": ";".join(interface_chain),
        "ppr_block_chain_5to3": module_core,
        "coding_frame_prefix_5to3": coding_frame_prefix,
        "binding_tract_cds_5to3": reconstructed,
        "ppr_block_chain_in_silico_validated": True,
        "final_cassette_vector_id": final["vector_id"],
        "final_cassette_five_prime_end_overhang": final[FIVE_PRIME_END],
        "final_cassette_three_prime_end_overhang": final[THREE_PRIME_END],
        "final_cassette_requirements_checked": True,
        "final_cassette_vector_sequence_provided": bool(
            interfaces_profile["final_cassette"].get("vector_sequence")
        ),
        "final_cassette_vector_sequence_in_silico_validated": False,
        "standalone_expression_cassette": False,
    }


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
    """Design target-specific dsDNA fragments through the configured entry vector.

    ``n_fragments`` and free custom flanks belonged to the former generic
    splitter and are rejected because they cannot guarantee interface
    compatibility. GRASP fixes five Level -1 modules per PPR block.
    """
    _ = fidelity  # retained only for source-level API compatibility
    if n_fragments is not None:
        raise ValueError(
            "GRASP one-shot uses fixed five-part Level 0 blocks; "
            "n_fragments is no longer configurable"
        )
    if max_fragment_cds is not None:
        raise ValueError("max_fragment_cds is not used by the GRASP module topology")
    if oligo_prefix is not None or oligo_suffix is not None:
        raise ValueError(
            "Custom one-shot flanks must be supplied through assembly_interfaces"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    info = describe_binder(target_rna)
    rna = info["target_rna"]
    architecture = f"{len(rna)}S"
    plan = compile_target_gap(
        rna,
        architecture=architecture,
        nterm_overhang=str(config.get("nterm_overhang", "AGGT")),
    )

    records = load_grasp_records([bundled_profile_genbank() / "GRASP_-1.gb"])
    parts_full = build_parts_table(records)
    local_config = copy.deepcopy(config)
    interfaces = resolve_assembly_interfaces(local_config)
    if (
        str(
            local_config.get("overhang_redesign", {}).get(
                "cut_mode", "native_fixed"
            )
        ).lower()
        == "movable_arelf"
    ):
        # One-shot does not run a Pareto search, but it must honor the dashboard's
        # configured inter-block overhangs and ARELF offsets.
        parts_full = materialize_arelf_parts(
            parts_full,
            {},
            config=local_config,
            codon_data=codon_data,
        )
    available = parts_full.set_index("part_id", drop=False)
    missing = sorted(set(plan["part_id"]) - set(available.index))
    if missing:
        raise ValueError(f"Deposited GRASP profile lacks selected modules: {missing}")

    # A repeated module is cloned once and reused in every listed assembly slot.
    unique_ids = list(dict.fromkeys(plan["part_id"].astype(str)))
    selected_parts = available.loc[unique_ids].reset_index(drop=True)

    local_config.setdefault("optimizer", {})["orthogonal_versions_per_part"] = 1
    local_config.setdefault("forbidden_sites", {}).update(
        {"BsaI": "GGTCTC", "BpiI": "GAAGAC"}
    )
    old_py_state = random.getstate()
    old_np_state = np.random.get_state()
    random.seed(seed)
    np.random.seed(seed)
    try:
        library = optimize_library(selected_parts, codon_data, local_config, log=log)
    finally:
        random.setstate(old_py_state)
        np.random.set_state(old_np_state)

    library = library[library["version"] == 1].copy().reset_index(drop=True)
    validation = [
        validate_order_fragment_in_silico(seq, interfaces)
        for seq in library["oligo_sequence_5to3"]
    ]
    validation_df = pd.DataFrame(validation)
    orderable = pd.concat([library.reset_index(drop=True), validation_df], axis=1)
    orderable.insert(0, "order_fragment_id", [f"GRASP_{pid}" for pid in orderable.part_id])
    orderable["sequence_type"] = "double-stranded DNA synthesis fragment"
    orderable["assembly_interface_profile"] = interfaces["profile_name"]
    orderable["entry_vector"] = interfaces["level_minus1_entry"]["vector_id"]
    orderable["entry_cloning_enzyme"] = ENTRY_CLONING_ENZYME
    orderable["module_release_enzyme"] = MODULE_RELEASE_ENZYME
    orderable["level0_acceptor"] = interfaces["level0"]["acceptor_id"]
    orderable["order_quantity"] = orderable["optimized_part_id"].map(
        plan["optimized_part_id"].value_counts()
    )
    orderable["used_in_groups"] = orderable["optimized_part_id"].map(
        plan.groupby("optimized_part_id", sort=False)["assembly_group"]
        .agg(lambda values: ";".join(dict.fromkeys(values)))
    )
    orderable["oligo_sequence_5to3"] = orderable["oligo_sequence_5to3"].map(clean_dna)
    orderable["order_sequence_5to3"] = orderable["oligo_sequence_5to3"]

    level0 = _validate_level0_groups(plan, orderable, interfaces)
    assembled = simulate_assembled_cds(
        plan,
        library,
        parts_full,
        genetic_code=int(local_config.get("genetic_code", 1)),
        config=local_config,
        codon_data=codon_data,
    )
    if not assembled.get("translation_verified"):
        raise AssertionError(
            "Target-specific GRASP modules did not reconstruct the expected binder"
        )
    ppr_chain = _validate_ppr_block_chain(
        plan,
        level0,
        library,
        parts_full,
        assembled["assembled_cds"],
        interfaces,
        architecture,
    )

    plan = plan.merge(
        orderable[
            [
                "optimized_part_id",
                "order_fragment_id",
                "module_release_oh5",
                "module_release_oh3",
            ]
        ],
        on="optimized_part_id",
        how="left",
        validate="many_to_one",
    )

    plan_path = output_dir / f"assembly_plan_{rna}.csv"
    order_csv = output_dir / f"oneshot_{rna}_orderable_fragments.csv"
    legacy_csv = output_dir / f"oneshot_{rna}_oligos.csv"
    order_fasta = output_dir / f"oneshot_{rna}_orderable_fragments.fasta"
    legacy_fasta = output_dir / f"oneshot_{rna}_oligos.fasta"
    level0_path = output_dir / f"oneshot_{rna}_ppr_blocks.csv"
    ppr_chain_path = output_dir / f"oneshot_{rna}_ppr_block_chain.csv"
    binding_tract_fasta = output_dir / f"oneshot_{rna}_binding_tract_context.fasta"

    plan.to_csv(plan_path, index=False)
    orderable.to_csv(order_csv, index=False)
    orderable.to_csv(legacy_csv, index=False)
    level0.to_csv(level0_path, index=False)
    pd.DataFrame([ppr_chain]).to_csv(ppr_chain_path, index=False)
    fasta_text = "".join(
        f">{row.order_fragment_id}|dsDNA|qty={row.order_quantity}|{row.entry_vector}_BsaI\n"
        f"{row.order_sequence_5to3}\n"
        for row in orderable.itertuples(index=False)
    )
    order_fasta.write_text(fasta_text)
    legacy_fasta.write_text(fasta_text)
    binding_tract_fasta.write_text(
        f">GRASP_{rna}|binding_tract_CDS_context|not_expression_cassette\n"
        f"{assembled['assembled_cds']}\n"
    )

    summary = {
        "target_rna": rna,
        "architecture": architecture,
        "ppr_code": info["ppr_code"],
        "aa_length": info["aa_length"],
        "n_assembly_slots": len(plan),
        "n_unique_order_fragments": len(orderable),
        "n_level0_assemblies": len(level0),
        "assembly_interface_profile": interfaces["profile_name"],
        "terminal_site_notation": interfaces["notation"],
        "coding_strand_direction": interfaces["coding_strand_direction"],
        "entry_vector": interfaces["level_minus1_entry"]["vector_id"],
        "entry_five_prime_end_overhang": interfaces["level_minus1_entry"][
            FIVE_PRIME_END
        ],
        "entry_three_prime_end_overhang": interfaces["level_minus1_entry"][
            THREE_PRIME_END
        ],
        "entry_five_prime_assembled_coding_site": interfaces["level_minus1_entry"][
            FIVE_PRIME_CODING_SITE
        ],
        "entry_three_prime_assembled_coding_site": interfaces["level_minus1_entry"][
            THREE_PRIME_CODING_SITE
        ],
        "entry_cloning_enzyme": ENTRY_CLONING_ENZYME,
        "module_release_enzyme": MODULE_RELEASE_ENZYME,
        "level0_acceptor": interfaces["level0"]["acceptor_id"],
        "ppr_five_prime_end_overhang": interfaces["level0"]["ppr_outer"][
            FIVE_PRIME_END
        ],
        "ppr_three_prime_end_overhang": interfaces["level0"]["ppr_outer"][
            THREE_PRIME_END
        ],
        "final_cassette_vector": interfaces["final_cassette"]["vector_id"],
        "final_cassette_five_prime_end_overhang": interfaces["final_cassette"][
            FIVE_PRIME_END
        ],
        "final_cassette_three_prime_end_overhang": interfaces["final_cassette"][
            THREE_PRIME_END
        ],
        "translation_verified": bool(assembled["translation_verified"]),
        "order_fragment_requirements_checked": bool(
            orderable["order_fragment_requirements_checked"].all()
        ),
        "entry_interface_requirements_checked": bool(
            orderable["entry_interface_requirements_checked"].all()
        ),
        "entry_vector_context_in_silico_validated": bool(
            orderable["entry_vector_context_in_silico_validated"].all()
        ),
        "entry_vector_sequence_in_silico_validated": bool(
            orderable["entry_vector_sequence_in_silico_validated"].all()
        ),
        "module_release_requirements_checked": bool(
            orderable["module_release_requirements_checked"].all()
        ),
        "level0_module_chain_in_silico_validated": bool(
            level0["level0_module_chain_in_silico_validated"].all()
        ),
        "level0_vector_sequence_in_silico_validated": bool(
            level0["level0_vector_sequence_in_silico_validated"].all()
        ),
        "ppr_block_chain_in_silico_validated": bool(
            ppr_chain["ppr_block_chain_in_silico_validated"]
        ),
        "final_cassette_requirements_checked": bool(
            ppr_chain["final_cassette_requirements_checked"]
        ),
        "final_cassette_vector_sequence_in_silico_validated": bool(
            ppr_chain["final_cassette_vector_sequence_in_silico_validated"]
        ),
        "standalone_expression_cassette": False,
        "fragments_clean_qc": bool(orderable["qc_passed"].all()),
        "fragments_hard_constraints_passed": bool(
            orderable["hard_constraints_passed"].all()
        ),
        "vendor_acceptance_confirmed": False,
        "organism": local_config.get("selected_organism")
        or local_config.get("selected_organism_label"),
    }
    pd.DataFrame([summary]).to_csv(
        output_dir / f"oneshot_{rna}_summary.csv", index=False
    )

    log(
        f"✓ One-shot {architecture}: {len(orderable)} unique dsDNA fragments; "
        f"{len(plan)} assembly uses; order geometry, PPR block chain, "
        f"and translation checked in silico"
    )
    return {
        "target_rna": rna,
        "binder": info,
        "aa_sequence": info["aa_sequence"],
        "binding_tract_cds": assembled["assembled_cds"],
        "gga_plan": plan,
        "assembly_plan": plan,
        "oligos": orderable,
        "orderable_fragments": orderable,
        "level0_assemblies": level0,
        "ppr_block_chain": ppr_chain,
        "assembled": assembled,
        "output_dir": output_dir,
        "oligo_csv": legacy_csv,
        "oligo_fasta": legacy_fasta,
        "order_csv": order_csv,
        "order_fasta": order_fasta,
        "level0_csv": level0_path,
        "ppr_block_chain_csv": ppr_chain_path,
        "binding_tract_fasta": binding_tract_fasta,
        "selected_overhangs": {},
        "pareto_front": None,
        "summary": summary,
    }


def subset_for_target(*_args, **_kwargs):
    raise RuntimeError(
        "subset_for_target is obsolete; run_oneshot_design now creates the "
        "target-specific pAGM1311 Level -1 subset directly"
    )
