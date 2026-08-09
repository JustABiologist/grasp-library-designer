"""Runnable workflows (no ipywidgets). Safe to call from notebook cells or CLI."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

import pandas as pd

from .dna import apply_overhang_to_mask, clean_mask, clean_dna
from .ligation_fidelity import LigationFidelityCalculator
from .objectives import evaluate_design, build_oligos_from_cds
from .pareto import knee_point, optimize_pareto_overhangs


def flanks_from_parts(parts_df: pd.DataFrame) -> Dict[str, Tuple[str, str]]:
    return {
        str(row.part_id): (
            clean_dna(row.oligo_prefix),
            clean_dna(row.oligo_suffix),
        )
        for row in parts_df.itertuples(index=False)
    }


def write_overhangs_into_parts(
    parts_df: pd.DataFrame,
    selection: Dict[str, str],
    junction_map: pd.DataFrame,
) -> pd.DataFrame:
    updated = parts_df.copy()
    mask_by_id = {
        row.part_id: clean_mask(row.coding_mask)
        for row in updated.itertuples(index=False)
    }
    for row in junction_map.itertuples(index=False):
        overhang = selection[str(row.junction)]
        mask_by_id[row.part_id] = apply_overhang_to_mask(
            mask_by_id[row.part_id],
            int(row.mask_start_0based),
            overhang,
        )
    updated["coding_mask"] = updated["part_id"].map(mask_by_id)
    return updated


def parse_overhang_selection(text: str) -> Dict[str, str]:
    selection: Dict[str, str] = {}
    for chunk in str(text).split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        junction, overhang = chunk.split("=", 1)
        selection[junction.strip()] = overhang.strip().upper().replace("U", "T")
    return selection


def make_cds_optimizer(
    parts_df: pd.DataFrame,
    junction_map: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    optimize_coding_sequence: Callable,
    *,
    greedy: bool = True,
):
    def optimize_cds_for_masks(selection: Dict[str, str]):
        local_config = copy.deepcopy(config)
        if greedy:
            local_config["optimizer"]["iterations_per_part"] = 0
        updated_parts = write_overhangs_into_parts(parts_df, selection, junction_map)
        cds_by_part, aa_by_part = {}, {}
        for row in updated_parts.itertuples(index=False):
            aa = str(row.aa_sequence).upper().replace(" ", "")
            sequence, _ = optimize_coding_sequence(
                aa_sequence=aa,
                coding_mask=clean_mask(row.coding_mask),
                codon_data=codon_data,
                config=local_config,
            )
            cds_by_part[row.part_id] = sequence
            aa_by_part[row.part_id] = aa
        return cds_by_part, aa_by_part

    return optimize_cds_for_masks


def run_overhang_redesign(
    *,
    parts: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    optimize_coding_sequence: Callable,
    input_dir: Path,
    output_dir: Path,
    seed: int = 42,
    fidelity: Optional[LigationFidelityCalculator] = None,
    log: Callable[[str], None] = print,
) -> Tuple[pd.DataFrame, Dict[str, str], pd.DataFrame]:
    """Pareto overhang search → selected set → parts with updated masks."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not config.get("overhang_redesign", {}).get("enabled", True):
        log("Overhang redesign disabled — returning parts unchanged.")
        return pd.DataFrame(), {}, parts

    if parts is None or len(parts) == 0:
        raise ValueError("parts is empty — import GRASP / load parts.csv first")
    if not codon_data:
        raise ValueError("codon_data empty — Apply codon table first")

    candidate_file = input_dir / "overhang_candidates.csv"
    junction_map_file = input_dir / "junction_map.csv"
    overhang_candidates = pd.read_csv(candidate_file)
    junction_map = pd.read_csv(junction_map_file)

    missing = sorted(set(junction_map["part_id"]) - set(parts["part_id"]))
    if missing:
        raise ValueError(f"junction_map unknown part_ids: {missing}")

    lig = config["ligation"]
    calc = fidelity or LigationFidelityCalculator(
        temperature=lig["temperature"],
        hours=lig["hours"],
        ligation_table=lig.get("ligation_table"),
        min_efficiency=lig.get("min_efficiency", 0.25),
        min_fidelity=lig.get("min_fidelity", 0.9),
    )

    optimize_cds_for_masks = make_cds_optimizer(
        parts,
        junction_map,
        codon_data,
        config,
        optimize_coding_sequence,
        greedy=True,
    )

    log(
        f"Searching ≤{config['pareto']['max_evaluations']} overhang sets "
        f"(beam {config['pareto']['beam_width']}, greedy CDS scoring)…"
    )
    front = optimize_pareto_overhangs(
        overhang_candidates,
        codon_data=codon_data,
        config=config,
        optimize_cds_for_masks=optimize_cds_for_masks,
        fidelity_calculator=calc,
        junction_map=junction_map,
        flanks_by_part=flanks_from_parts(parts),
        max_evaluations=config["pareto"]["max_evaluations"],
        beam_width=config["pareto"]["beam_width"],
        seed=seed,
    )
    if front is None or front.empty:
        raise RuntimeError("Pareto front empty")

    mode = config.get("overhang_redesign", {}).get("selection", "knee")
    if mode == "knee":
        chosen = knee_point(front)
    elif mode == "max_fidelity":
        chosen = front.sort_values("ligation_fidelity", ascending=False).iloc[0]
    else:
        match = front[front["overhangs"] == mode]
        if match.empty:
            raise ValueError(f"No Pareto point matches selection={mode!r}")
        chosen = match.iloc[0]

    selected = parse_overhang_selection(chosen["overhangs"])
    log(
        f"Selected ({mode}): {chosen['overhangs']}\n"
        f"  fidelity={chosen['ligation_fidelity']:.4f}  "
        f"codon={chosen['codon_optimality']:.4f}  "
        f"synthesis={chosen['synthesis']:.4f}"
    )

    updated = write_overhangs_into_parts(parts, selected, junction_map)
    front.to_csv(output_dir / "pareto_front.csv", index=False)
    pd.DataFrame([chosen]).to_csv(output_dir / "selected_overhangs.csv", index=False)
    updated.to_csv(output_dir / "parts_with_redesigned_junctions.csv", index=False)
    log(f"Wrote masks → {output_dir / 'parts_with_redesigned_junctions.csv'}")
    return front, selected, updated


def run_library_optimize(
    *,
    parts: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    optimize_library: Callable,
    output_dir: Path,
    log: Callable[[str], None] = print,
) -> pd.DataFrame:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if parts is None or len(parts) == 0:
        raise ValueError("parts is empty")
    if not codon_data:
        raise ValueError("codon_data empty")

    n = len(parts)
    iters = config["optimizer"]["iterations_per_part"]
    versions = config["optimizer"]["orthogonal_versions_per_part"]
    log(f"Annealing {n} parts × {versions} version(s) × {iters:,} iters…")
    library = optimize_library(parts, codon_data, config)
    out = output_dir / "optimized_library.csv"
    library.to_csv(out, index=False)
    log(f"Done — {len(library)} sequences → {out}")
    return library


def _cds_maps_from_library(library: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str]]:
    cds_by_part: Dict[str, str] = {}
    aa_by_part: Dict[str, str] = {}
    for row in library.itertuples(index=False):
        pid = str(row.part_id)
        # Prefer version 1 if multiple orthogonal versions exist
        if pid in cds_by_part and getattr(row, "version", 1) != 1:
            continue
        cds_by_part[pid] = clean_dna(row.optimized_cds)
        aa_by_part[pid] = str(row.aa_sequence).upper().replace(" ", "")
    return cds_by_part, aa_by_part


def rescore_pareto_front_after_anneal(
    front: pd.DataFrame,
    *,
    parts: pd.DataFrame,
    junction_map: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    optimize_coding_sequence: Callable,
    selected_overhangs: Optional[Dict[str, str]] = None,
    optimized_library: Optional[pd.DataFrame] = None,
    fidelity: Optional[LigationFidelityCalculator] = None,
    deep_all: bool = False,
    output_dir: Optional[Path] = None,
    log: Callable[[str], None] = print,
) -> pd.DataFrame:
    """
    Re-score the Pareto front after junction redesign + library anneal.

    Codon = full CDS; synthesis = full oligo (prefix + CDS + suffix).
    Selected overhang set uses `optimized_library`; other members use greedy
    CDS by default (`deep_all=True` fully anneals each).
    """
    if front is None or front.empty:
        raise ValueError("Pareto front is empty")
    if optimized_library is None or len(optimized_library) == 0:
        raise ValueError(
            "optimized_library is required — run full CDS anneal before plotting"
        )

    lig = config.get("ligation", {})
    calc = fidelity or LigationFidelityCalculator(
        temperature=lig.get("temperature", 25),
        hours=lig.get("hours", 18),
        ligation_table=lig.get("ligation_table"),
        min_efficiency=lig.get("min_efficiency", 0.25),
        min_fidelity=lig.get("min_fidelity", 0.9),
    )
    flanks = flanks_from_parts(parts)

    selected_key = None
    if selected_overhangs:
        selected_key = ";".join(
            f"{k}={v}" for k, v in sorted(selected_overhangs.items())
        )

    greedy_optimize = make_cds_optimizer(
        parts,
        junction_map,
        codon_data,
        config,
        optimize_coding_sequence,
        greedy=True,
    )
    deep_optimize = None
    if deep_all:
        deep_optimize = make_cds_optimizer(
            parts,
            junction_map,
            codon_data,
            config,
            optimize_coding_sequence,
            greedy=False,
        )

    rows = []
    n = len(front)
    for i, row in enumerate(front.itertuples(index=False), start=1):
        overhangs_text = str(row.overhangs)
        selection = parse_overhang_selection(overhangs_text)
        oligos_by_part = None

        if selected_key is not None and overhangs_text == selected_key:
            log(f"[{i}/{n}] Scoring selected set from full CDS anneal…")
            cds_by_part, aa_by_part = _cds_maps_from_library(optimized_library)
            if "oligo_sequence_5to3" in optimized_library.columns:
                oligos_by_part = {}
                for lib_row in optimized_library.itertuples(index=False):
                    pid = str(lib_row.part_id)
                    if pid in oligos_by_part and getattr(lib_row, "version", 1) != 1:
                        continue
                    oligos_by_part[pid] = clean_dna(lib_row.oligo_sequence_5to3)
            source = "full_anneal_library"
        elif deep_all and deep_optimize is not None:
            log(f"[{i}/{n}] Full-anneal scoring overhang set…")
            cds_by_part, aa_by_part = deep_optimize(selection)
            source = "full_anneal"
        else:
            log(f"[{i}/{n}] Scoring overhang set (greedy CDS + full oligo)…")
            cds_by_part, aa_by_part = greedy_optimize(selection)
            source = "greedy_cds"

        if oligos_by_part is None:
            oligos_by_part = build_oligos_from_cds(cds_by_part, flanks)

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
            oligos_by_part=oligos_by_part,
            flanks_by_part=flanks,
            junction_map=junction_map,
        )
        rows.append(
            {
                "overhangs": overhangs_text,
                "ligation_fidelity": scores.ligation_fidelity,
                "codon_optimality": scores.codon_optimality,
                "synthesis": scores.synthesis,
                "fidelity_beam": getattr(row, "fidelity_beam", scores.ligation_fidelity),
                "synthesis_scope": "full_oligo",
                "score_source": source,
            }
        )

    rescored = pd.DataFrame(rows)
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / "pareto_front_after_anneal.csv"
        rescored.to_csv(path, index=False)
        log(f"Wrote post-anneal front → {path}")
    return rescored
