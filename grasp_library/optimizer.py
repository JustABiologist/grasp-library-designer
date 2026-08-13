"""Masked codon / synthesis optimizer used by library and one-shot notebooks."""

from __future__ import annotations

import json
import math
import random
from collections import Counter
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .codon_validation import cds_matches_organism, verify_cds_for_organism
from .dna import (
    clean_dna,
    clean_mask,
    contains_forbidden_site,
    format_forbidden_hits,
    format_gc_windows,
    format_homopolymer_runs,
    format_pct,
    format_repeat_hits,
    gc_fraction,
    gc_window_violations,
    kmer_counts,
    local_gc_penalty,
    longest_homopolymer,
    mask_matches,
    notable_homopolymers,
    repeated_kmer_hits,
    repeated_kmer_penalty,
    reverse_complement,
    translate_dna,
)


def codon_matches_mask(codon: str, codon_mask: str) -> bool:
    return all(
        mask_base == "N" or codon_base == mask_base
        for codon_base, mask_base in zip(codon, codon_mask)
    )


def build_allowed_codons(
    aa_sequence: str,
    coding_mask: str,
    codon_data: Mapping[str, Sequence[dict]],
    minimum_relative_adaptiveness: float = 0.20,
) -> List[List[dict]]:
    aa_sequence = str(aa_sequence).upper().replace(" ", "")
    coding_mask = clean_mask(coding_mask)

    if len(coding_mask) != 3 * len(aa_sequence):
        raise ValueError(
            f"Mask length {len(coding_mask)} does not match "
            f"{len(aa_sequence)} amino acids. Expected: "
            f"{3 * len(aa_sequence)}."
        )

    allowed = []
    for index, amino_acid in enumerate(aa_sequence):
        if amino_acid not in codon_data:
            raise ValueError(f"No codons available for amino acid {amino_acid!r}.")

        codon_mask = coding_mask[3 * index : 3 * index + 3]
        candidates = [
            item
            for item in codon_data[amino_acid]
            if codon_matches_mask(item["codon"], codon_mask)
        ]
        if not candidates:
            raise ValueError(
                f"No synonymous codons for amino acid {amino_acid} "
                f"at position {index + 1} match the mask {codon_mask}."
            )

        preferred = [
            item
            for item in candidates
            if item["relative_adaptiveness"] >= minimum_relative_adaptiveness
        ]
        # Mask compatibility takes priority over codon frequency.
        allowed.append(preferred if preferred else candidates)
    return allowed


def codon_score(
    sequence: str,
    aa_sequence: str,
    codon_data: Mapping[str, Sequence[dict]],
    epsilon: float = 1e-12,
) -> float:
    """Mean log relative adaptiveness (0 best; negative worse)."""
    sequence = clean_dna(sequence)
    aa_sequence = str(aa_sequence).upper()

    lookup = {}
    for aa, entries in codon_data.items():
        for entry in entries:
            lookup[(aa, entry["codon"])] = entry["relative_adaptiveness"]

    scores = []
    for index, aa in enumerate(aa_sequence):
        codon = sequence[3 * index : 3 * index + 3]
        adaptiveness = lookup.get((aa, codon), 0.0)
        scores.append(math.log(max(adaptiveness, epsilon)))
    return float(np.mean(scores)) if scores else 0.0


def global_gc_penalty(sequence: str, gc_min: float, gc_max: float) -> float:
    value = gc_fraction(sequence)
    if value < gc_min:
        return (gc_min - value) ** 2
    if value > gc_max:
        return (value - gc_max) ** 2
    return 0.0


def library_similarity_penalty(
    sequence: str,
    external_kmers: Optional[Counter],
    k: int = 12,
) -> float:
    if not external_kmers:
        return 0.0
    counts = kmer_counts(sequence, k)
    return float(
        sum(count * external_kmers.get(kmer, 0) for kmer, count in counts.items())
    )


def sequence_objective(
    sequence: str,
    aa_sequence: str,
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    external_kmers: Optional[Counter] = None,
) -> float:
    synthesis = config["synthesis"]
    weights = config["weights"]

    if contains_forbidden_site(sequence, config["forbidden_sites"]):
        return -1e12

    score = 0.0
    score += weights["codon"] * codon_score(sequence, aa_sequence, codon_data)
    score -= weights["global_gc"] * global_gc_penalty(
        sequence, synthesis["global_gc_min"], synthesis["global_gc_max"]
    )
    score -= weights["local_gc"] * local_gc_penalty(
        sequence,
        synthesis["window_size"],
        synthesis["window_gc_min"],
        synthesis["window_gc_max"],
    )
    excess_homopolymer = max(
        0, longest_homopolymer(sequence) - synthesis["max_homopolymer"]
    )
    score -= weights["homopolymer"] * excess_homopolymer**2
    score -= weights["internal_repeat"] * repeated_kmer_penalty(
        sequence, k=synthesis["repeat_k"]
    )
    score -= weights["library_similarity"] * library_similarity_penalty(
        sequence, external_kmers or Counter(), k=synthesis["repeat_k"]
    )
    return float(score)


def weighted_initial_sequence(allowed_codons: Sequence[Sequence[dict]]) -> str:
    codons = []
    for candidates in allowed_codons:
        probabilities = np.array(
            [item["probability"] for item in candidates], dtype=float
        )
        probabilities /= probabilities.sum()
        selected = np.random.choice(len(candidates), p=probabilities)
        codons.append(candidates[selected]["codon"])
    return "".join(codons)


def greedy_coding_sequence(allowed_codons: Sequence[Sequence[dict]]) -> str:
    """Fast encode: best relative-adaptiveness codon per position (mask-safe)."""
    return "".join(
        max(candidates, key=lambda c: c["relative_adaptiveness"])["codon"]
        for candidates in allowed_codons
    )


def _codon_indices_for_span(start: int, length: int, n_codons: int) -> List[int]:
    last = (start + max(length, 1) - 1) // 3
    return [index for index in range(start // 3, last + 1) if 0 <= index < n_codons]


def repair_forbidden_sites(
    sequence: str,
    allowed_codons: Sequence[Sequence[dict]],
    forbidden_sites: Mapping[str, str],
    *,
    max_passes: int = 80,
) -> Optional[str]:
    """Recode overlapping synonymous codons until blacklist sites are gone."""
    sequence = clean_dna(sequence)
    if not forbidden_sites:
        return sequence
    if not contains_forbidden_site(sequence, forbidden_sites):
        return sequence

    adaptiveness = {
        (index, item["codon"]): item["relative_adaptiveness"]
        for index, candidates in enumerate(allowed_codons)
        for item in candidates
    }
    codons = [sequence[i : i + 3] for i in range(0, len(sequence), 3)]
    for _ in range(max_passes):
        hits = contains_forbidden_site("".join(codons), forbidden_sites)
        if not hits:
            return "".join(codons)
        repaired = False
        for hit in hits:
            site_len = len(str(hit["site"]))
            for index in _codon_indices_for_span(
                int(hit["start_0based"]), site_len, len(codons)
            ):
                current = codons[index]
                alternatives = sorted(
                    (
                        item["codon"]
                        for item in allowed_codons[index]
                        if item["codon"] != current
                    ),
                    key=lambda codon: adaptiveness.get((index, codon), 0.0),
                    reverse=True,
                )
                for alt in alternatives:
                    trial = list(codons)
                    trial[index] = alt
                    if len(contains_forbidden_site("".join(trial), forbidden_sites)) < len(
                        hits
                    ):
                        codons = trial
                        repaired = True
                        break
                if repaired:
                    break
            if repaired:
                break
        if not repaired:
            return None
    return None


def initial_sequence_without_forbidden_sites(
    allowed_codons: Sequence[Sequence[dict]],
    forbidden_sites: Mapping[str, str],
    *,
    random_tries: int = 40,
) -> str:
    """Greedy encode, then recode around blacklist hits; sample if needed."""
    seed = greedy_coding_sequence(allowed_codons)
    repaired = repair_forbidden_sites(seed, allowed_codons, forbidden_sites)
    if repaired is not None:
        return repaired
    for _ in range(random_tries):
        candidate = weighted_initial_sequence(allowed_codons)
        repaired = repair_forbidden_sites(
            candidate, allowed_codons, forbidden_sites
        )
        if repaired is not None:
            return repaired
    hits = contains_forbidden_site(seed, forbidden_sites)
    detail = format_forbidden_hits(hits) or "unknown site"
    raise RuntimeError(
        "No synonymous CDS avoids the cut-site blacklist "
        f"({detail}). The coding mask may lock a forbidden site."
    )


def optimize_coding_sequence(
    aa_sequence: str,
    coding_mask: str,
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    external_kmers: Optional[Counter] = None,
    iterations: Optional[int] = None,
):
    aa_sequence = str(aa_sequence).upper().replace(" ", "")
    coding_mask = clean_mask(coding_mask)
    minimum_adaptiveness = config["codon_optimization"][
        "minimum_relative_adaptiveness"
    ]
    allowed_codons = build_allowed_codons(
        aa_sequence,
        coding_mask,
        codon_data,
        minimum_relative_adaptiveness=minimum_adaptiveness,
    )

    if iterations is None:
        iterations = config["optimizer"]["iterations_per_part"]

    forbidden_sites = dict(config.get("forbidden_sites") or {})
    current_sequence = initial_sequence_without_forbidden_sites(
        allowed_codons, forbidden_sites
    )
    current_codons = [
        current_sequence[i : i + 3] for i in range(0, len(current_sequence), 3)
    ]
    current_score = sequence_objective(
        current_sequence, aa_sequence, codon_data, config, external_kmers=external_kmers
    )

    # iterations <= 0 → greedy / repaired only (used during Pareto search)
    if iterations <= 0:
        return current_sequence, current_score

    best_sequence = current_sequence
    best_score = current_score

    mutable_positions = [
        index
        for index, candidates in enumerate(allowed_codons)
        if len(candidates) > 1
    ]
    if not mutable_positions:
        return best_sequence, best_score

    t_initial = config["optimizer"]["initial_temperature"]
    t_final = config["optimizer"]["final_temperature"]

    for iteration in range(iterations):
        fraction = iteration / max(1, iterations - 1)
        temperature = t_initial * (t_final / t_initial) ** fraction
        position = random.choice(mutable_positions)
        old_codon = current_codons[position]
        alternatives = [
            item["codon"]
            for item in allowed_codons[position]
            if item["codon"] != old_codon
        ]
        if not alternatives:
            continue

        new_codon = random.choice(alternatives)
        candidate_codons = current_codons.copy()
        candidate_codons[position] = new_codon
        candidate_sequence = "".join(candidate_codons)
        candidate_score = sequence_objective(
            candidate_sequence,
            aa_sequence,
            codon_data,
            config,
            external_kmers=external_kmers,
        )
        delta = candidate_score - current_score
        accept = delta >= 0 or random.random() < math.exp(
            delta / max(temperature, 1e-12)
        )
        if accept:
            current_codons = candidate_codons
            current_sequence = candidate_sequence
            current_score = candidate_score
            if current_score > best_score:
                best_sequence = current_sequence
                best_score = current_score

    if not cds_matches_organism(
        best_sequence,
        aa_sequence,
        genetic_code=int(config["genetic_code"]),
        codon_data=codon_data,
    ):
        raise AssertionError(
            "Optimization failed organism codon-table translation check "
            f"(genetic code {config['genetic_code']})."
        )
    if not mask_matches(best_sequence, coding_mask):
        raise AssertionError("Optimized sequence violates the coding mask.")

    return best_sequence, best_score


def synthesis_qc(
    sequence: str,
    config: Mapping,
    *,
    forbidden_scan: Optional[str] = None,
    sequence_kind: str = "generic",
):
    """QC a DNA string with explicit clean/warning/failure semantics.

    ``passed`` means a clean pass: no soft guideline warnings and no hard
    failures. ``hard_constraints_passed`` remains true for a sequence that has
    only warnings. Length limits are checked only when ``sequence_kind`` is
    ``"oligo"`` or ``"gene"``; this avoids applying gene-product limits to an
    individual GRASP module CDS.

    This is a transparent heuristic, not confirmation that a synthesis vendor
    will accept an order. Non-machine-checkable vendor rules are returned for
    manual review.
    """
    synthesis = config["synthesis"]
    sequence = clean_dna(sequence)
    scan = clean_dna(forbidden_scan) if forbidden_scan is not None else sequence

    valid_kinds = {"generic", "cds", "oligo", "gene"}
    if sequence_kind not in valid_kinds:
        raise ValueError(
            f"Unknown sequence_kind {sequence_kind!r}; expected one of {sorted(valid_kinds)}"
        )

    gc = gc_fraction(sequence)
    homopolymer = longest_homopolymer(sequence)
    repeat_penalty = repeated_kmer_penalty(sequence, synthesis["repeat_k"])
    forbidden_sites = dict(config.get("forbidden_sites") or {})
    forbidden_hits = contains_forbidden_site(scan, forbidden_sites)
    local_windows = gc_window_violations(
        sequence,
        synthesis["window_size"],
        synthesis["window_gc_min"],
        synthesis["window_gc_max"],
    )
    local_penalty = local_gc_penalty(
        sequence,
        synthesis["window_size"],
        synthesis["window_gc_min"],
        synthesis["window_gc_max"],
    )
    over_homopolymers = [
        run
        for run in notable_homopolymers(sequence, synthesis["max_homopolymer"])
        if run["length"] > synthesis["max_homopolymer"]
    ]
    homopolymer_detail = format_homopolymer_runs(
        notable_homopolymers(sequence, synthesis["max_homopolymer"])
    )
    repeat_hits = repeated_kmer_hits(sequence, synthesis["repeat_k"])
    forbidden_detail = format_forbidden_hits(forbidden_hits)
    local_gc_detail = format_gc_windows(local_windows)
    gc_min = float(synthesis["global_gc_min"])
    gc_max = float(synthesis["global_gc_max"])
    if gc > gc_max:
        gc_status = f"HIGH {format_pct(gc)} (max {format_pct(gc_max)})"
    elif gc < gc_min:
        gc_status = f"LOW {format_pct(gc)} (min {format_pct(gc_min)})"
    else:
        gc_status = f"OK {format_pct(gc)}"

    warnings = []
    failures = []
    if gc > gc_max:
        warnings.append(f"Global GC {gc_status}")
    elif gc < gc_min:
        warnings.append(f"Global GC {gc_status}")
    if local_windows:
        warnings.append(f"Local GC outside target range: {local_gc_detail}")
    if over_homopolymers:
        warnings.append(
            f"Homopolymer too long: {format_homopolymer_runs(over_homopolymers)}"
        )
    if repeat_penalty > synthesis["max_repeat_count"]:
        detail = format_repeat_hits(repeat_hits)
        warnings.append(
            f"Elevated internal DNA repeats: {detail}" if detail
            else "Elevated internal DNA repeats"
        )
    if forbidden_hits:
        failures.append(
            f"Forbidden restriction site in CDS: {forbidden_detail}"
        )

    if sequence_kind == "oligo":
        minimum = synthesis.get("min_oligo_length")
        maximum = synthesis.get("max_oligo_length")
        if minimum is not None and len(sequence) < int(minimum):
            failures.append(f"Oligo shorter than configured minimum ({minimum} bp)")
        if maximum is not None and len(sequence) > int(maximum):
            failures.append(f"Oligo longer than configured maximum ({maximum} bp)")
    elif sequence_kind == "gene":
        minimum = synthesis.get("min_gene_length")
        maximum = synthesis.get("max_gene_length")
        if minimum is not None and len(sequence) < int(minimum):
            failures.append(f"Gene shorter than configured minimum ({minimum} bp)")
        if maximum is not None and len(sequence) > int(maximum):
            failures.append(f"Gene longer than configured maximum ({maximum} bp)")

    vendor_meta = config.get("synthesis_vendor_meta", {})
    machine_hard = vendor_meta.get("machine_hard_constraints", {}) or {}
    hard_max_homopolymer = machine_hard.get("max_homopolymer")
    if (
        hard_max_homopolymer is not None
        and homopolymer > int(hard_max_homopolymer)
    ):
        failures.append(
            "Homopolymer exceeds vendor hard maximum "
            f"({hard_max_homopolymer} nt)"
            + (f": {homopolymer_detail}" if homopolymer_detail else "")
        )

    # Preserve order while avoiding duplicate diagnostics.
    warnings = list(dict.fromkeys(warnings))
    failures = list(dict.fromkeys(failures))
    hard_constraints_passed = len(failures) == 0
    passed = hard_constraints_passed and len(warnings) == 0
    if failures:
        status = "FAIL"
    elif warnings:
        status = "WARNING"
    else:
        status = "PASS"

    manual_rules = vendor_meta.get(
        "manual_review_rules", vendor_meta.get("hard_rules", [])
    ) or []

    return {
        "length": len(sequence),
        "gc_fraction": gc,
        "gc_pct": round(gc * 100, 1),
        "gc_status": gc_status,
        "longest_homopolymer": homopolymer,
        "homopolymer_detail": homopolymer_detail,
        "homopolymer_count": len(over_homopolymers),
        "repeat_penalty": repeat_penalty,
        "repeat_detail": format_repeat_hits(repeat_hits),
        "local_gc_penalty": local_penalty,
        "local_gc_detail": local_gc_detail,
        "forbidden_hits": json.dumps(forbidden_hits),
        "forbidden_detail": forbidden_detail,
        "blacklist_tested": ", ".join(forbidden_sites) or "none",
        "blacklist_hits": forbidden_detail or "none",
        "warnings": "; ".join(warnings),
        "failures": "; ".join(failures),
        "warning_count": len(warnings),
        "failure_count": len(failures),
        "status": status,
        "hard_constraints_passed": hard_constraints_passed,
        "passed": passed,
        "vendor_acceptance_confirmed": False,
        "manual_vendor_rules": "; ".join(str(rule) for rule in manual_rules),
    }


def optimize_library(
    parts: pd.DataFrame,
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    *,
    log=print,
) -> pd.DataFrame:
    """Anneal every part; return oligo table."""
    results = []
    library_kmers: Counter = Counter()
    number_of_versions = config["optimizer"]["orthogonal_versions_per_part"]
    repeat_k = config["synthesis"]["repeat_k"]
    parts = parts.reset_index(drop=True)

    for part_index, row in parts.iterrows():
        aa_sequence = str(row["aa_sequence"]).upper().replace(" ", "")
        coding_mask = clean_mask(row["coding_mask"])
        log(f"[{part_index + 1}/{len(parts)}] Optimizing {row['part_id']}")

        part_versions = []
        for version in range(1, number_of_versions + 1):
            sequence, objective = optimize_coding_sequence(
                aa_sequence=aa_sequence,
                coding_mask=coding_mask,
                codon_data=codon_data,
                config=config,
                external_kmers=library_kmers,
            )
            if {"oh5_mask_start", "oh3_mask_start"} <= set(parts.columns):
                from .assembly_interfaces import order_fragment_arms, resolve_assembly_interfaces
                from .import_grasp import build_configured_order_fragment

                interfaces = resolve_assembly_interfaces(config)
                full_oligo = build_configured_order_fragment(
                    sequence,
                    part_id=str(row["part_id"]),
                    oh5_mask_start=int(row["oh5_mask_start"]),
                    oh3_mask_start=int(row["oh3_mask_start"]),
                    interfaces=interfaces,
                )
                actual_prefix, actual_suffix = order_fragment_arms(interfaces)
            else:
                raise ValueError(
                    "parts table lacks oh5_mask_start/oh3_mask_start; re-import "
                    "the bundled GRASP profile before generating order fragments"
                )
            cds_qc = synthesis_qc(sequence, config, sequence_kind="cds")
            oligo_qc = synthesis_qc(
                full_oligo,
                config,
                forbidden_scan=sequence,
                sequence_kind="oligo",
            )
            translation = verify_cds_for_organism(
                sequence,
                aa_sequence,
                genetic_code=int(config["genetic_code"]),
                codon_data=codon_data,
            )

            result = {
                "part_id": row["part_id"],
                "version": version,
                "optimized_part_id": f"{row['part_id']}_v{version}",
                "aa_sequence": aa_sequence,
                "coding_mask": coding_mask,
                "optimized_cds": sequence,
                "oligo_prefix": actual_prefix,
                "oligo_suffix": actual_suffix,
                "oligo_sequence_5to3": full_oligo,
                "objective_score": objective,
                "codon_score": codon_score(sequence, aa_sequence, codon_data),
                "translation_verified": bool(translation["ok"]),
                "genetic_code_ok": bool(translation["genetic_code_ok"]),
                "codon_table_ok": bool(translation["codon_table_ok"]),
                "genetic_code": int(config["genetic_code"]),
                "mask_verified": mask_matches(sequence, coding_mask),
                "cds_gc": cds_qc["gc_fraction"],
                "cds_gc_pct": cds_qc["gc_pct"],
                "cds_gc_status": cds_qc["gc_status"],
                "cds_longest_homopolymer": cds_qc["longest_homopolymer"],
                "cds_homopolymers": cds_qc["homopolymer_detail"],
                "cds_local_gc": cds_qc["local_gc_detail"],
                "cds_repeat_penalty": cds_qc["repeat_penalty"],
                "cds_repeats": cds_qc["repeat_detail"],
                "cds_warnings": cds_qc["warnings"],
                "cds_failures": cds_qc["failures"],
                "blacklist_tested": cds_qc["blacklist_tested"],
                "blacklist_hits": cds_qc["blacklist_hits"],
                "oligo_length": oligo_qc["length"],
                "oligo_gc": oligo_qc["gc_fraction"],
                "oligo_gc_pct": oligo_qc["gc_pct"],
                "oligo_gc_status": oligo_qc["gc_status"],
                "oligo_homopolymers": oligo_qc["homopolymer_detail"],
                "oligo_warnings": oligo_qc["warnings"],
                "oligo_failures": oligo_qc["failures"],
                "qc_status": (
                    "FAIL"
                    if (
                        not cds_qc["hard_constraints_passed"]
                        or not oligo_qc["hard_constraints_passed"]
                        or not bool(translation["ok"])
                        or not mask_matches(sequence, coding_mask)
                    )
                    else "WARNING"
                    if (cds_qc["warnings"] or oligo_qc["warnings"])
                    else "PASS"
                ),
                "hard_constraints_passed": (
                    cds_qc["hard_constraints_passed"]
                    and oligo_qc["hard_constraints_passed"]
                    and bool(translation["ok"])
                    and mask_matches(sequence, coding_mask)
                ),
                "vendor_acceptance_confirmed": False,
                "manual_vendor_rules": oligo_qc["manual_vendor_rules"],
                "qc_passed": (
                    cds_qc["passed"]
                    and oligo_qc["passed"]
                    and bool(translation["ok"])
                    and mask_matches(sequence, coding_mask)
                ),
            }
            for metadata_column in (
                "oh5_coding_site_5to3",
                "oh3_coding_site_5to3",
                "five_prime_end_overhang",
                "three_prime_end_overhang",
                "overhang_notation",
            ):
                if metadata_column in row.index:
                    result[metadata_column] = row[metadata_column]
            part_versions.append(result)
            results.append(result)
            library_kmers.update(kmer_counts(sequence, repeat_k))

        sequences = [item["optimized_cds"] for item in part_versions]
        if len(set(sequences)) != len(sequences):
            log(f"  Warning: not all versions for {row['part_id']} were DNA-distinct.")

    from .binder import annotate_module_roles

    return annotate_module_roles(pd.DataFrame(results))


def simulate_assembled_cds(
    assembly_plan: pd.DataFrame,
    optimized_library: pd.DataFrame,
    parts_full: Optional[pd.DataFrame] = None,
    *,
    genetic_code: int = 1,
    config: Optional[Mapping] = None,
    codon_data: Optional[Mapping[str, Sequence[dict]]] = None,
) -> Dict:
    """Stitch optimized module CDS with shared 4-nt overhangs deduplicated.

    Uses overhang-bounded inserts from ``parts_full`` coordinates:
    keep ``cds[oh5 : oh3+4]`` (plus any 5′-side lead bases before ``oh5``
    on the first module). Downstream modules drop the shared 5′ overhang.
    This trims GenBank window trail past the 3′ overhang (e.g. C modules)
    so native 9S ORFs stay in-frame and match the GAP binder protein.
    """
    from .binder import rna_to_binder_aa

    lib_cols = ["optimized_part_id", "optimized_cds"]
    if "aa_sequence" not in assembly_plan.columns:
        lib_cols.append("aa_sequence")
    if "part_id" not in assembly_plan.columns:
        lib_cols.append("part_id")

    selected = assembly_plan.merge(
        optimized_library[lib_cols],
        on="optimized_part_id",
        how="left",
        validate="many_to_one",
    )
    # Normalize possible merge suffixes
    if "part_id" not in selected.columns and "part_id_x" in selected.columns:
        selected = selected.rename(columns={"part_id_x": "part_id"})
    if "aa_sequence" not in selected.columns and "aa_sequence_x" in selected.columns:
        selected = selected.rename(columns={"aa_sequence_x": "aa_sequence"})
    if "aa_sequence" not in selected.columns and "aa_sequence_y" in selected.columns:
        selected = selected.rename(columns={"aa_sequence_y": "aa_sequence"})

    if selected["optimized_cds"].isna().any():
        missing = selected.loc[
            selected["optimized_cds"].isna(), "optimized_part_id"
        ].tolist()
        raise ValueError(f"Missing optimized sequences for: {missing}")

    sort_cols = [
        c
        for c in ("assembly_group", "assembly_order", "assembly_slot")
        if c in selected.columns
    ]
    selected = selected.sort_values(sort_cols)

    chunks = []
    full_index = (
        parts_full.set_index("part_id") if parts_full is not None else None
    )
    used_overhang_bounds = False
    previous_right_site = None
    previous_three_prime_end = None
    previous_group = None
    coding_junctions_checked = 0
    directional_terminal_pairs_checked = 0

    for i, row in enumerate(selected.itertuples(index=False)):
        cds = clean_dna(row.optimized_cds)
        part_id = str(row.part_id)
        if (
            full_index is not None
            and part_id in full_index.index
            and "oh5_mask_start" in full_index.columns
            and "oh3_mask_start" in full_index.columns
        ):
            meta = full_index.loc[part_id]
            oh5 = int(meta.oh5_mask_start)
            oh3 = int(meta.oh3_mask_start)
            if not (0 <= oh5 < oh3 + 4 <= len(cds)):
                raise ValueError(
                    f"{part_id}: overhang window oh5={oh5} oh3={oh3} "
                    f"incompatible with CDS length {len(cds)}"
                )
            # Bounded insert: 5′ overhang … 3′ overhang (inclusive)
            insert = cds[oh5 : oh3 + 4]
            left_site = insert[:4]
            right_site = insert[-4:]
            declared_left = clean_dna(
                meta.get("oh5_coding_site_5to3", meta.get("oh5", left_site))
            )
            declared_right = clean_dna(
                meta.get("oh3_coding_site_5to3", meta.get("oh3", right_site))
            )
            if left_site != declared_left or right_site != declared_right:
                raise ValueError(
                    f"{part_id}: optimized coding sites {left_site}/{right_site} "
                    f"do not match declared sites {declared_left}/{declared_right}"
                )
            if previous_right_site is not None:
                coding_junctions_checked += 1
                if previous_right_site != left_site:
                    raise ValueError(
                        f"{part_id}: assembled coding-strand junction mismatch "
                        f"{previous_right_site} != {left_site}"
                    )

            current_group = getattr(row, "assembly_group", None)
            current_five_prime_end = meta.get("five_prime_end_overhang")
            crossing_blocks = (
                previous_group is not None
                and current_group is not None
                and str(previous_group) != str(current_group)
            )
            if crossing_blocks and previous_three_prime_end is not None and pd.notna(
                current_five_prime_end
            ):
                upstream = clean_dna(previous_three_prime_end)
                downstream = clean_dna(current_five_prime_end)
                directional_terminal_pairs_checked += 1
                if reverse_complement(upstream) != downstream:
                    raise ValueError(
                        f"{previous_group}->{current_group}: directional terminal "
                        f"mismatch reverse_complement({upstream}) != {downstream}"
                    )
            if i == 0:
                chunks.append(cds[:oh5] + insert)
            else:
                chunks.append(insert[4:])
            previous_right_site = right_site
            previous_three_prime_end = meta.get("three_prime_end_overhang")
            previous_group = current_group
            used_overhang_bounds = True
        else:
            # Fallback when parts_full coordinates are unavailable
            chunks.append(cds)

    assembled_cds = "".join(chunks)

    expected_protein = ""
    expected_source = ""
    if "target_rna" in selected.columns and len(selected):
        rna = str(selected["target_rna"].iloc[0]).strip()
        if rna and rna.upper() != "NAN":
            expected_protein = rna_to_binder_aa(rna)
            expected_source = "binder_scaffold"

    result = {
        "assembled_cds": assembled_cds,
        "expected_protein": expected_protein,
        "observed_protein": "",
        "translation_verified": False,
        "stitch_warning": "",
        "stitch_mode": (
            "overhang_bounded" if used_overhang_bounds else "full_cds_concat"
        ),
        "expected_source": expected_source,
        "coding_junctions_checked": coding_junctions_checked,
        "directional_terminal_pairs_checked": directional_terminal_pairs_checked,
        "coding_junctions_verified": bool(
            used_overhang_bounds and coding_junctions_checked == max(0, len(selected) - 1)
        ),
    }
    if len(assembled_cds) % 3 != 0:
        result["stitch_warning"] = (
            f"Assembled CDS length {len(assembled_cds)} is not divisible by 3; "
            "skipped ORF translation check."
        )
        if config is not None:
            qc = synthesis_qc(assembled_cds, config, sequence_kind="gene")
            qc["passed"] = False
            result.update(qc)
        return result

    observed_protein = translate_dna(assembled_cds, genetic_code)
    result["observed_protein"] = observed_protein
    result["genetic_code"] = int(genetic_code)

    if not expected_protein:
        # No independent scaffold — treat ORF translation as the reference.
        expected_protein = observed_protein
        expected_source = "assembled_translation"
        result["expected_protein"] = expected_protein
        result["expected_source"] = expected_source

    if codon_data is not None:
        tx = verify_cds_for_organism(
            assembled_cds,
            expected_protein,
            genetic_code=int(genetic_code),
            codon_data=codon_data,
        )
        result["observed_protein"] = tx["observed_protein"]
        result["translation_verified"] = bool(tx["ok"])
        result["genetic_code_ok"] = bool(tx["genetic_code_ok"])
        result["codon_table_ok"] = bool(tx["codon_table_ok"])
    else:
        result["translation_verified"] = observed_protein == expected_protein

    if config is not None:
        result.update(synthesis_qc(assembled_cds, config, sequence_kind="gene"))
    return result
