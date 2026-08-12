"""Three design objectives: ligation fidelity, codon optimality, synthesis."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .dna import (
    clean_dna,
    contains_forbidden_site,
    gc_fraction,
    local_gc_penalty,
    longest_homopolymer,
    repeated_kmer_penalty,
)
from .ligation_fidelity import LigationFidelityCalculator


@dataclass(frozen=True)
class ObjectiveScores:
    """All three objectives are maximized."""

    ligation_fidelity: float
    codon_optimality: float
    synthesis: float

    def as_tuple(self) -> tuple[float, float, float]:
        return (
            self.ligation_fidelity,
            self.codon_optimality,
            self.synthesis,
        )

    def as_dict(self) -> dict:
        return asdict(self)


def codon_optimality_score(
    sequence: str,
    aa_sequence: str,
    codon_data: Mapping[str, Sequence[dict]],
    epsilon: float = 1e-12,
    *,
    codon_indices: Optional[Sequence[int]] = None,
) -> float:
    """Mean log relative adaptiveness (0 theoretical best, negative worse)."""
    sequence = clean_dna(sequence)
    aa_sequence = str(aa_sequence).upper()

    lookup = {}
    for aa, entries in codon_data.items():
        for entry in entries:
            lookup[(aa, entry["codon"])] = entry["relative_adaptiveness"]

    if codon_indices is None:
        codon_indices = range(len(aa_sequence))

    scores = []
    for index in codon_indices:
        aa = aa_sequence[index]
        codon = sequence[3 * index : 3 * index + 3]
        adaptiveness = lookup.get((aa, codon), 0.0)
        scores.append(math.log(max(adaptiveness, epsilon)))
    return float(np.mean(scores)) if scores else 0.0


def synthesis_score(
    sequence: str,
    config: Mapping,
    external_kmers: Optional[Counter] = None,
    *,
    include_library_similarity: bool = True,
    include_global_gc: bool = True,
    include_internal_repeat: bool = True,
    forbidden_scan: Optional[str] = None,
    prefer_ideal: bool = False,
) -> float:
    """
    Higher is better. Soft synthesis penalties are subtracted; forbidden sites
    yield a large negative constant (hard failure).

    `forbidden_scan`: if set, Type IIS sites are checked only on this substring
    (use the CDS when scoring full oligos whose MoClo flanks intentionally
    contain BpiI/BsaI).

    `prefer_ideal`: when True, always apply continuous preferences (GC→0.5,
    soft homopolymer ramp) so Pareto ranking still moves under loose vendor
    hard limits.
    """
    synthesis = config["synthesis"]
    weights = config["weights"]
    sequence = clean_dna(sequence)
    scan = clean_dna(forbidden_scan) if forbidden_scan is not None else sequence

    if contains_forbidden_site(scan, config["forbidden_sites"]):
        return -1e6

    score = 0.0
    gc = gc_fraction(sequence)

    if include_global_gc:
        if gc < synthesis["global_gc_min"]:
            score -= weights["global_gc"] * (synthesis["global_gc_min"] - gc) ** 2
        elif gc > synthesis["global_gc_max"]:
            score -= weights["global_gc"] * (gc - synthesis["global_gc_max"]) ** 2

    if prefer_ideal:
        ideal_gc = float(synthesis.get("ideal_gc", 0.5))
        score -= float(weights.get("ideal_gc", 2.0)) * (gc - ideal_gc) ** 2

    window_size = min(int(synthesis["window_size"]), max(len(sequence), 1))
    score -= weights["local_gc"] * local_gc_penalty(
        sequence,
        window_size,
        synthesis["window_gc_min"],
        synthesis["window_gc_max"],
    )

    homo = longest_homopolymer(sequence)
    excess_homopolymer = max(0, homo - synthesis["max_homopolymer"])
    score -= weights["homopolymer"] * excess_homopolymer**2

    if prefer_ideal:
        soft_hp = int(synthesis.get("soft_homopolymer", 6))
        soft_excess = max(0, homo - soft_hp)
        score -= 0.35 * weights["homopolymer"] * soft_excess**2

    if include_internal_repeat:
        score -= weights["internal_repeat"] * repeated_kmer_penalty(
            sequence,
            k=synthesis["repeat_k"],
        )

    if include_library_similarity and external_kmers:
        from .dna import kmer_counts

        counts = kmer_counts(sequence, synthesis["repeat_k"])
        similarity = sum(
            count * external_kmers.get(kmer, 0) for kmer, count in counts.items()
        )
        score -= weights["library_similarity"] * similarity

    return float(score)


def build_oligos_from_cds(
    cds_by_part: Mapping[str, str],
    flanks_by_part: Mapping[str, tuple],
    config: Optional[Mapping] = None,
) -> Dict[str, str]:
    """Build coordinate-safe synthesis fragments for the configured entry vector."""
    from .import_grasp import build_configured_order_fragment

    oligos: Dict[str, str] = {}
    for part_id, cds in cds_by_part.items():
        if part_id not in flanks_by_part:
            continue
        flank_data = flanks_by_part[part_id]
        if len(flank_data) != 4:
            raise ValueError(
                f"{part_id}: missing overhang coordinates for order-fragment construction"
            )
        _, _, oh5, oh3 = flank_data
        oligos[part_id] = build_configured_order_fragment(
            cds,
            part_id=part_id,
            oh5_mask_start=int(oh5),
            oh3_mask_start=int(oh3),
            config=dict(config or {}),
        )
    return oligos


def oligo_synthesis_score(
    oligos_by_part: Mapping[str, str],
    cds_by_part: Mapping[str, str],
    config: Mapping,
    *,
    cross_part_similarity: bool = False,
) -> float:
    """
    Mean synthesis fitness over full oligos (prefix + CDS + suffix).

    Soft constraints (GC, homopolymers, repeats) use the whole oligo.
    Forbidden Type IIS sites are scanned on the CDS only — MoClo flanks
    intentionally carry BpiI/BsaI.
    """
    part_ids = [p for p in oligos_by_part if p in cds_by_part]
    if not part_ids:
        return 0.0

    from .dna import kmer_counts

    external: Optional[Counter] = Counter() if cross_part_similarity else None
    scores = []
    repeat_k = config["synthesis"]["repeat_k"]
    for part_id in part_ids:
        oligo = oligos_by_part[part_id]
        cds = cds_by_part[part_id]
        scores.append(
            synthesis_score(
                oligo,
                config,
                external_kmers=external,
                include_library_similarity=cross_part_similarity,
                forbidden_scan=cds,
                prefer_ideal=True,
            )
        )
        if external is not None:
            external.update(kmer_counts(oligo, repeat_k))
    return float(np.mean(scores))


def library_codon_score(
    sequences: Iterable[str],
    aa_sequences: Iterable[str],
    codon_data: Mapping[str, Sequence[dict]],
) -> float:
    values = [
        codon_optimality_score(seq, aa, codon_data)
        for seq, aa in zip(sequences, aa_sequences)
    ]
    return float(np.mean(values)) if values else 0.0


def library_synthesis_score(
    sequences: Iterable[str],
    config: Mapping,
) -> float:
    seqs = list(sequences)
    if not seqs:
        return 0.0
    external = Counter()
    scores = []
    from .dna import kmer_counts

    repeat_k = config["synthesis"]["repeat_k"]
    for seq in seqs:
        scores.append(synthesis_score(seq, config, external_kmers=external))
        external.update(kmer_counts(seq, repeat_k))
    return float(np.mean(scores))


def _junction_span(
    start_0based: int,
    overhang_len: int = 4,
    flank: int = 9,
    seq_len: Optional[int] = None,
) -> Tuple[int, int]:
    lo = max(0, start_0based - flank)
    hi = start_0based + overhang_len + flank
    if seq_len is not None:
        hi = min(seq_len, hi)
    return lo, hi


def extract_junction_windows(
    cds_by_part: Mapping[str, str],
    junction_map: pd.DataFrame,
    *,
    overhang_len: int = 4,
    flank: int = 9,
) -> List[str]:
    """DNA windows centered on each part's coded overhang (junction ± flank)."""
    windows: List[str] = []
    for row in junction_map.itertuples(index=False):
        part_id = str(row.part_id)
        if part_id not in cds_by_part:
            continue
        seq = clean_dna(cds_by_part[part_id])
        start = int(row.mask_start_0based)
        lo, hi = _junction_span(start, overhang_len, flank, len(seq))
        if hi > lo:
            windows.append(seq[lo:hi])
    return windows


def junction_codon_indices(
    junction_map: pd.DataFrame,
    part_id: str,
    *,
    overhang_len: int = 4,
) -> List[int]:
    """Codon indices whose triplets overlap a junction overhang on this part."""
    indices = set()
    for row in junction_map.itertuples(index=False):
        if str(row.part_id) != part_id:
            continue
        start = int(row.mask_start_0based)
        first = start // 3
        last = (start + overhang_len - 1) // 3
        indices.update(range(first, last + 1))
    return sorted(indices)


def junction_codon_score(
    cds_by_part: Mapping[str, str],
    aa_by_part: Mapping[str, str],
    junction_map: pd.DataFrame,
    codon_data: Mapping[str, Sequence[dict]],
    *,
    overhang_len: int = 4,
) -> float:
    """Mean log relative adaptiveness over codons that touch junction overhangs."""
    values = []
    for part_id, seq in cds_by_part.items():
        aa = aa_by_part[part_id]
        indices = junction_codon_indices(
            junction_map, part_id, overhang_len=overhang_len
        )
        indices = [i for i in indices if 0 <= i < len(aa)]
        if not indices:
            continue
        values.append(
            codon_optimality_score(
                seq, aa, codon_data, codon_indices=indices
            )
        )
    return float(np.mean(values)) if values else 0.0


def junction_synthesis_score(
    cds_by_part: Mapping[str, str],
    junction_map: pd.DataFrame,
    config: Mapping,
    *,
    overhang_len: int = 4,
    flank: int = 9,
) -> float:
    """
    Synthesis fitness for junction-local DNA only (overhang ± flank nt).

    Skips whole-library k-mer similarity and internal-repeat terms that dominate
    full-CDS scoring and drown out junction trade-offs.
    """
    windows = extract_junction_windows(
        cds_by_part,
        junction_map,
        overhang_len=overhang_len,
        flank=flank,
    )
    if not windows:
        return 0.0
    scores = [
        synthesis_score(
            window,
            config,
            external_kmers=None,
            include_library_similarity=False,
            include_global_gc=False,
            include_internal_repeat=False,
        )
        for window in windows
    ]
    return float(np.mean(scores))


def evaluate_design(
    *,
    overhangs: Sequence[str] | Mapping[str, str],
    cds_sequences: Sequence[str],
    aa_sequences: Sequence[str],
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    fidelity_calculator: Optional[LigationFidelityCalculator] = None,
    cds_by_part: Optional[Mapping[str, str]] = None,
    aa_by_part: Optional[Mapping[str, str]] = None,
    oligos_by_part: Optional[Mapping[str, str]] = None,
    flanks_by_part: Optional[Mapping[str, Tuple[str, str]]] = None,
    junction_map: Optional[pd.DataFrame] = None,
    junction_flank: int = 9,
) -> ObjectiveScores:
    """
    Score a complete design on all three Pareto objectives.

    - Ligation: if ``overhangs`` is a junction-name mapping, report the
      orientation-invariant score of one physical Level 0 reaction. Every PPR
      block uses the same junction set in a separate tube, so the objective does
      not multiply unrelated transformations/screening steps.
      A bare sequence retains the legacy single-reaction behavior for generic
      Golden Gate and one-shot designs.
    - Codon: mean log relative adaptiveness over full CDS sequences.
    - Synthesis: mean fitness over full oligos (prefix + CDS + suffix) when
      flanks/oligos are provided; otherwise full CDS. Forbidden sites are
      checked on CDS only so intentional MoClo sites in flanks do not fail QC.
    """
    _ = (junction_map, junction_flank)  # API compat; codon/synthesis use full sequences

    calc = fidelity_calculator or LigationFidelityCalculator(
        temperature=config.get("ligation", {}).get("temperature", 25),
        hours=config.get("ligation", {}).get("hours", 18),
    )

    if isinstance(overhangs, Mapping) and overhangs:
        from .assembly_interfaces import resolve_assembly_interfaces

        assembly_profile = resolve_assembly_interfaces(config)
        configured_outer = assembly_profile["level0"].get("acceptor_outer")
        external_overhangs = (
            (
                configured_outer["n_overhang_5p"],
                configured_outer["c_overhang_5p"],
            )
            if configured_outer is not None
            else None
        )
        fidelity = calc.grasp_first_stage_fidelity(
            overhangs,
            architecture=str(config.get("architecture", "9S")),
            external_overhangs=external_overhangs,
        )
    elif overhangs:
        fidelity = calc.set_fidelity(overhangs)
    else:
        fidelity = 1.0

    if cds_by_part is not None and aa_by_part is not None:
        part_ids = sorted(cds_by_part)
        codon = library_codon_score(
            [cds_by_part[p] for p in part_ids],
            [aa_by_part[p] for p in part_ids],
            codon_data,
        )
        if oligos_by_part is None and flanks_by_part is not None:
            oligos_by_part = build_oligos_from_cds(
                cds_by_part, flanks_by_part, config=config
            )
        if oligos_by_part:
            synthesis = oligo_synthesis_score(
                oligos_by_part,
                cds_by_part,
                config,
                cross_part_similarity=False,
            )
        else:
            synthesis = float(
                np.mean(
                    [
                        synthesis_score(cds_by_part[p], config, external_kmers=None)
                        for p in part_ids
                    ]
                )
            )
    else:
        codon = library_codon_score(cds_sequences, aa_sequences, codon_data)
        synthesis = library_synthesis_score(cds_sequences, config)

    return ObjectiveScores(
        ligation_fidelity=fidelity,
        codon_optimality=codon,
        synthesis=synthesis,
    )
