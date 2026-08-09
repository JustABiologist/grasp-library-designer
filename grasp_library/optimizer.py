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
    gc_fraction,
    kmer_counts,
    local_gc_penalty,
    longest_homopolymer,
    mask_matches,
    repeated_kmer_penalty,
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

    # iterations <= 0 → greedy only (used during Pareto search)
    if iterations <= 0:
        sequence = greedy_coding_sequence(allowed_codons)
        score = sequence_objective(
            sequence, aa_sequence, codon_data, config, external_kmers=external_kmers
        )
        return sequence, score

    current_sequence = greedy_coding_sequence(allowed_codons)
    current_codons = [
        current_sequence[i : i + 3] for i in range(0, len(current_sequence), 3)
    ]
    current_score = sequence_objective(
        current_sequence, aa_sequence, codon_data, config, external_kmers=external_kmers
    )

    for _ in range(30):
        if current_score > -1e11:
            break
        current_sequence = weighted_initial_sequence(allowed_codons)
        current_codons = [
            current_sequence[i : i + 3] for i in range(0, len(current_sequence), 3)
        ]
        current_score = sequence_objective(
            current_sequence,
            aa_sequence,
            codon_data,
            config,
            external_kmers=external_kmers,
        )
    else:
        raise RuntimeError(
            "No initial sequence without forbidden sites found. "
            "The coding mask may be too restrictive."
        )

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


def synthesis_qc(sequence: str, config: Mapping, *, forbidden_scan: Optional[str] = None):
    """QC a DNA string. Forbidden sites scanned on `forbidden_scan` if given."""
    synthesis = config["synthesis"]
    sequence = clean_dna(sequence)
    scan = clean_dna(forbidden_scan) if forbidden_scan is not None else sequence

    gc = gc_fraction(sequence)
    homopolymer = longest_homopolymer(sequence)
    repeat_penalty = repeated_kmer_penalty(sequence, synthesis["repeat_k"])
    forbidden_hits = contains_forbidden_site(scan, config["forbidden_sites"])
    local_penalty = local_gc_penalty(
        sequence,
        synthesis["window_size"],
        synthesis["window_gc_min"],
        synthesis["window_gc_max"],
    )

    warnings = []
    failures = []
    if not (synthesis["global_gc_min"] <= gc <= synthesis["global_gc_max"]):
        warnings.append("Global GC outside target range")
    if local_penalty > 0:
        warnings.append("Local GC outside target range")
    if homopolymer > synthesis["max_homopolymer"]:
        warnings.append("Homopolymer too long")
    if repeat_penalty > synthesis["max_repeat_count"]:
        warnings.append("Elevated internal DNA repeats")
    if forbidden_hits:
        failures.append("Forbidden restriction site in CDS")

    return {
        "length": len(sequence),
        "gc_fraction": gc,
        "longest_homopolymer": homopolymer,
        "repeat_penalty": repeat_penalty,
        "local_gc_penalty": local_penalty,
        "forbidden_hits": json.dumps(forbidden_hits),
        "warnings": "; ".join(warnings),
        "failures": "; ".join(failures),
        "passed": len(failures) == 0,
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
            full_oligo = (
                clean_dna(row["oligo_prefix"])
                + sequence
                + clean_dna(row["oligo_suffix"])
            )
            cds_qc = synthesis_qc(sequence, config)
            oligo_qc = synthesis_qc(full_oligo, config, forbidden_scan=sequence)
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
                "oligo_prefix": clean_dna(row["oligo_prefix"]),
                "oligo_suffix": clean_dna(row["oligo_suffix"]),
                "oligo_sequence_5to3": full_oligo,
                "objective_score": objective,
                "codon_score": codon_score(sequence, aa_sequence, codon_data),
                "translation_verified": bool(translation["ok"]),
                "genetic_code_ok": bool(translation["genetic_code_ok"]),
                "codon_table_ok": bool(translation["codon_table_ok"]),
                "genetic_code": int(config["genetic_code"]),
                "mask_verified": mask_matches(sequence, coding_mask),
                "cds_gc": cds_qc["gc_fraction"],
                "cds_longest_homopolymer": cds_qc["longest_homopolymer"],
                "cds_repeat_penalty": cds_qc["repeat_penalty"],
                "oligo_length": oligo_qc["length"],
                "oligo_gc": oligo_qc["gc_fraction"],
                "oligo_warnings": oligo_qc["warnings"],
                "oligo_failures": oligo_qc["failures"],
                "qc_passed": (
                    cds_qc["passed"]
                    and oligo_qc["passed"]
                    and bool(translation["ok"])
                    and mask_matches(sequence, coding_mask)
                ),
            }
            part_versions.append(result)
            results.append(result)
            library_kmers.update(kmer_counts(sequence, repeat_k))

        sequences = [item["optimized_cds"] for item in part_versions]
        if len(set(sequences)) != len(sequences):
            log(f"  Warning: not all versions for {row['part_id']} were DNA-distinct.")

    return pd.DataFrame(results)


def simulate_assembled_cds(
    assembly_plan: pd.DataFrame,
    optimized_library: pd.DataFrame,
    parts_full: Optional[pd.DataFrame] = None,
    *,
    genetic_code: int = 1,
    config: Optional[Mapping] = None,
    codon_data: Optional[Mapping[str, Sequence[dict]]] = None,
) -> Dict:
    """Stitch optimized module CDS with shared 4-nt overhangs deduplicated."""
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
    expected_chunks = []
    full_index = (
        parts_full.set_index("part_id") if parts_full is not None else None
    )

    for i, row in enumerate(selected.itertuples(index=False)):
        cds = clean_dna(row.optimized_cds)
        aa = str(row.aa_sequence).upper().replace(" ", "")
        part_id = str(row.part_id)
        if full_index is not None and part_id in full_index.index:
            meta = full_index.loc[part_id]
            oh5 = int(meta.oh5_mask_start)
            if i == 0:
                chunks.append(cds)
                expected_chunks.append(aa)
            else:
                chunks.append(cds[oh5 + 4 :])
                drop_aa = (oh5 + 4) // 3
                expected_chunks.append(aa[drop_aa:])
        else:
            chunks.append(cds)
            expected_chunks.append(aa)

    assembled_cds = "".join(chunks)
    expected_protein = "".join(expected_chunks)
    result = {
        "assembled_cds": assembled_cds,
        "expected_protein": expected_protein,
        "observed_protein": "",
        "translation_verified": False,
        "stitch_warning": "",
    }
    if len(assembled_cds) % 3 != 0:
        # Shared 4-nt overhangs are not always codon-aligned at every junction
        # (e.g. D modules with oh5 at index 0). Per-module translation is still
        # verified during anneal; assembled ORF QC is best-effort.
        result["stitch_warning"] = (
            f"Assembled CDS length {len(assembled_cds)} is not divisible by 3; "
            "skipped ORF translation check."
        )
        if config is not None:
            qc = synthesis_qc(assembled_cds, config)
            qc["passed"] = False
            result.update(qc)
        return result

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
        result["genetic_code"] = int(genetic_code)
    else:
        observed_protein = translate_dna(assembled_cds, genetic_code)
        result["observed_protein"] = observed_protein
        result["translation_verified"] = observed_protein == expected_protein
        result["genetic_code"] = int(genetic_code)
    if config is not None:
        result.update(synthesis_qc(assembled_cds, config))
    return result
