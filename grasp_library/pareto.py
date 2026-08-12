"""Pareto-front search over ligation fidelity, codon optimality, and synthesis."""

from __future__ import annotations

import math
import random
from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .arelf import format_cut_token, parse_cut_token, selection_overhangs
from .dna import (
    apply_overhang_to_mask,
    clean_dna,
    is_self_reverse_complement,
    reverse_complement,
)
from .ligation_fidelity import LigationFidelityCalculator, fidelity_calculator_for_level
from .objectives import ObjectiveScores, evaluate_design


@dataclass
class ParetoPoint:
    overhangs: Dict[str, str]
    scores: ObjectiveScores
    cds_by_part: Dict[str, str] = field(default_factory=dict)
    meta: Dict = field(default_factory=dict)

    def as_row(self) -> dict:
        row = {
            "overhangs": ";".join(f"{k}={v}" for k, v in sorted(self.overhangs.items())),
            "ligation_fidelity": self.scores.ligation_fidelity,
            "level_minus1_fidelity": self.scores.level_minus1_fidelity,
            "level0_fidelity": self.scores.level0_fidelity or self.scores.ligation_fidelity,
            "level1_fidelity": self.scores.level1_fidelity,
            "codon_optimality": self.scores.codon_optimality,
            "synthesis": self.scores.synthesis,
        }
        row.update(self.meta)
        return row


def dominates(a: ObjectiveScores, b: ObjectiveScores) -> bool:
    """True if a Pareto-dominates b (all objectives maximized)."""
    a_vals = a.as_tuple()
    b_vals = b.as_tuple()
    return all(x >= y for x, y in zip(a_vals, b_vals)) and any(
        x > y for x, y in zip(a_vals, b_vals)
    )


def pareto_front(points: Sequence[ParetoPoint]) -> List[ParetoPoint]:
    front: List[ParetoPoint] = []
    for candidate in points:
        if any(dominates(other.scores, candidate.scores) for other in points if other is not candidate):
            continue
        front.append(candidate)
    return front


def _valid_overhang_set(selection: Mapping[str, str]) -> bool:
    overhangs = list(selection_overhangs(selection).values())
    if any(is_self_reverse_complement(o) for o in overhangs):
        return False
    used = set()
    for overhang in overhangs:
        rev = reverse_complement(overhang)
        if overhang in used or rev in used:
            return False
        used.add(overhang)
        used.add(rev)
    return True


def score_overhang_set_ggassembler(
    selection: Mapping[str, str],
    fidelity_calculator: LigationFidelityCalculator,
    *,
    architecture: str = "9S",
    external_overhangs: Optional[Sequence[str]] = None,
) -> float:
    if not selection:
        return 1.0
    overhang_selection = selection_overhangs(selection)
    outer = tuple(
        external_overhangs
        or fidelity_calculator.GRASP_LEVEL0_EXTERNAL_OVERHANGS
    )
    if len(outer) != 2:
        raise ValueError("Level 0 external overhangs require 5′ and 3′ values")
    physical = {
        "J_level0_left": outer[0],
        **overhang_selection,
        "J_level0_right": outer[1],
    }
    if not _valid_overhang_set(physical):
        return 0.0
    try:
        complete = set(fidelity_calculator.GRASP_LEVEL0_INTERNAL_JUNCTIONS) <= set(
            overhang_selection
        )
        if complete:
            return fidelity_calculator.grasp_first_stage_fidelity(
                overhang_selection,
                architecture=architecture,
                external_overhangs=outer,
            )
        # Partial beam states cannot yet form a physical reaction.  Retain a
        # deterministic generic score only to rank them until all four named
        # internal junctions are present.
        return fidelity_calculator.set_fidelity(overhang_selection.values())
    except (KeyError, ValueError):
        return 0.0


def optimize_pareto_overhangs(
    candidates: pd.DataFrame,
    *,
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    optimize_cds_for_masks: Callable[[Dict[str, str]], Tuple[Dict[str, str], Dict[str, str]]],
    fidelity_calculator: Optional[LigationFidelityCalculator] = None,
    junction_map: Optional[pd.DataFrame] = None,
    flanks_by_part: Optional[Mapping[str, tuple]] = None,
    max_evaluations: int = 2_000,
    beam_width: int = 200,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Search overhang combinations and evaluate the three-objective Pareto front.

    `candidates` columns: junction, overhang [, fidelity optional]
    `optimize_cds_for_masks(overhang_selection) -> (cds_by_part, aa_by_part)`
    should apply selected overhangs into coding masks, codon-optimize each part
    for synthesis+codon (single-objective SA), and return resulting CDS/AA maps.

    Synthesis is scored on full oligos (prefix + CDS + suffix) when
    `flanks_by_part` is provided; codon on the full CDS.
    """
    random.seed(seed)

    calc = fidelity_calculator or fidelity_calculator_for_level(
        config,
        str(config.get("ligation", {}).get("redesign_level", "level0")),
    )
    from .assembly_interfaces import (
        FIVE_PRIME_CODING_SITE,
        THREE_PRIME_CODING_SITE,
        resolve_assembly_interfaces,
    )

    assembly_profile = resolve_assembly_interfaces(config)
    configured_outer = assembly_profile["level0"].get("acceptor_outer")
    external_overhangs = (
        (
            configured_outer[FIVE_PRIME_CODING_SITE],
            configured_outer[THREE_PRIME_CODING_SITE],
        )
        if configured_outer is not None
        else calc.GRASP_LEVEL0_EXTERNAL_OVERHANGS
    )

    table = candidates.copy()
    table["junction"] = table["junction"].astype(str)
    table["overhang"] = (
        table["overhang"].astype(str).str.upper().str.replace("U", "T", regex=False)
    )
    table = table[~table["overhang"].map(is_self_reverse_complement)].copy()

    if "cut_token" in table.columns:
        table["candidate_token"] = table["cut_token"].astype(str)
    elif "motif_offset_nt" in table.columns:
        table["candidate_token"] = [
            format_cut_token(overhang, int(offset))
            for overhang, offset in zip(table["overhang"], table["motif_offset_nt"])
        ]
    else:
        # Backwards-compatible fixed-cut candidate tables carry only the 4-mer.
        table["candidate_token"] = table["overhang"]
    grouped = {
        junction: list(dict.fromkeys(group["candidate_token"].tolist()))
        for junction, group in table.groupby("junction")
    }
    junctions = sorted(grouped, key=lambda j: len(grouped[j]))

    # Beam search on GGAssembler fidelity, then fully evaluate top states.
    states: List[Tuple[Dict[str, str], float]] = [({}, 1.0)]
    for junction in junctions:
        new_states: List[Tuple[Dict[str, str], float]] = []
        for selection, _ in states:
            used = set(selection_overhangs(selection).values()) | set(
                external_overhangs
            )
            used_rc = {reverse_complement(o) for o in used}
            for candidate_token in grouped[junction]:
                overhang, _ = parse_cut_token(candidate_token)
                if overhang in used or overhang in used_rc:
                    continue
                if reverse_complement(overhang) in used:
                    continue
                new_selection = dict(selection)
                new_selection[junction] = candidate_token
                fidelity = score_overhang_set_ggassembler(
                    new_selection,
                    calc,
                    architecture=str(config.get("architecture", "9S")),
                    external_overhangs=external_overhangs,
                )
                if fidelity <= 0:
                    continue
                new_states.append((new_selection, fidelity))
        if not new_states:
            raise RuntimeError(f"No compatible solution from junction {junction}.")
        new_states.sort(key=lambda item: item[1], reverse=True)
        states = new_states[:beam_width]

    # Cap evaluations; sample if the beam is large.
    if len(states) > max_evaluations:
        states = random.sample(states, max_evaluations)

    points: List[ParetoPoint] = []
    for selection, fidelity_guess in states:
        try:
            optimized = optimize_cds_for_masks(selection)
            if len(optimized) == 2:
                cds_by_part, aa_by_part = optimized
                point_flanks = flanks_by_part
                point_junction_map = junction_map
            elif len(optimized) == 4:
                (
                    cds_by_part,
                    aa_by_part,
                    point_flanks,
                    point_junction_map,
                ) = optimized
            else:
                raise ValueError(
                    "CDS optimizer must return (cds, aa) or "
                    "(cds, aa, dynamic_flanks, dynamic_junction_map)"
                )
        except Exception as exc:  # keep searching other sets
            points.append(
                ParetoPoint(
                    overhangs=selection,
                    scores=ObjectiveScores(0.0, -1e6, -1e6),
                    meta={"error": str(exc), "fidelity_beam": fidelity_guess},
                )
            )
            continue

        part_ids = sorted(cds_by_part)
        scores = evaluate_design(
            overhangs=selection_overhangs(selection),
            cds_sequences=[cds_by_part[p] for p in part_ids],
            aa_sequences=[aa_by_part[p] for p in part_ids],
            codon_data=codon_data,
            config=config,
            fidelity_calculator=calc,
            cds_by_part=cds_by_part,
            aa_by_part=aa_by_part,
            flanks_by_part=point_flanks,
            junction_map=point_junction_map,
        )
        points.append(
            ParetoPoint(
                overhangs=selection,
                scores=scores,
                cds_by_part=cds_by_part,
                meta={
                    "fidelity_beam": fidelity_guess,
                    "synthesis_scope": "full_oligo" if point_flanks else "full_cds",
                },
            )
        )

    front = pareto_front(
        [p for p in points if p.scores.ligation_fidelity > 0 and p.scores.synthesis > -1e5]
    )
    if not front:
        front = pareto_front(points)

    rows = [p.as_row() for p in sorted(front, key=lambda p: p.scores.ligation_fidelity, reverse=True)]
    return pd.DataFrame(rows)


def knee_point(front: pd.DataFrame) -> pd.Series:
    """
    Pick a compromise point: minimize distance to the utopia point after
    min-max normalizing the three maximized objectives.
    """
    if front.empty:
        raise ValueError("Empty Pareto front.")

    cols = ["ligation_fidelity", "codon_optimality", "synthesis"]
    values = front[cols].astype(float)
    mins = values.min()
    maxs = values.max()
    spans = (maxs - mins).replace(0, 1.0)
    norm = (values - mins) / spans
    utopia = norm.max()
    dist = ((norm - utopia) ** 2).sum(axis=1).pow(0.5)
    return front.loc[dist.idxmin()]
