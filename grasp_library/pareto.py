"""Pareto-front search over ligation fidelity, codon optimality, and synthesis."""

from __future__ import annotations

import math
import random
from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .dna import (
    apply_overhang_to_mask,
    clean_dna,
    is_self_reverse_complement,
    reverse_complement,
)
from .ligation_fidelity import LigationFidelityCalculator
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
    overhangs = list(selection.values())
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
) -> float:
    if not selection:
        return 1.0
    if not _valid_overhang_set(selection):
        return 0.0
    try:
        return fidelity_calculator.set_fidelity(selection.values())
    except ValueError:
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

    calc = fidelity_calculator or LigationFidelityCalculator(
        temperature=config.get("ligation", {}).get("temperature", 25),
        hours=config.get("ligation", {}).get("hours", 18),
    )

    table = candidates.copy()
    table["junction"] = table["junction"].astype(str)
    table["overhang"] = (
        table["overhang"].astype(str).str.upper().str.replace("U", "T", regex=False)
    )
    table = table[~table["overhang"].map(is_self_reverse_complement)].copy()

    grouped = {
        junction: group["overhang"].tolist()
        for junction, group in table.groupby("junction")
    }
    junctions = sorted(grouped, key=lambda j: len(grouped[j]))

    # Beam search on GGAssembler fidelity, then fully evaluate top states.
    states: List[Tuple[Dict[str, str], float]] = [({}, 1.0)]
    for junction in junctions:
        new_states: List[Tuple[Dict[str, str], float]] = []
        for selection, _ in states:
            used = set(selection.values())
            used_rc = {reverse_complement(o) for o in used}
            for overhang in grouped[junction]:
                if overhang in used or overhang in used_rc:
                    continue
                if reverse_complement(overhang) in used:
                    continue
                new_selection = dict(selection)
                new_selection[junction] = overhang
                fidelity = score_overhang_set_ggassembler(new_selection, calc)
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
            cds_by_part, aa_by_part = optimize_cds_for_masks(selection)
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
            overhangs=list(selection.values()),
            cds_sequences=[cds_by_part[p] for p in part_ids],
            aa_sequences=[aa_by_part[p] for p in part_ids],
            codon_data=codon_data,
            config=config,
            fidelity_calculator=calc,
            cds_by_part=cds_by_part,
            aa_by_part=aa_by_part,
            flanks_by_part=flanks_by_part,
            junction_map=junction_map,
        )
        points.append(
            ParetoPoint(
                overhangs=selection,
                scores=scores,
                cds_by_part=cds_by_part,
                meta={
                    "fidelity_beam": fidelity_guess,
                    "synthesis_scope": "full_oligo" if flanks_by_part else "full_cds",
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
