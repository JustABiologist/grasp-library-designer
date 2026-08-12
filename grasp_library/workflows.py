"""Runnable workflows (no ipywidgets). Safe to call from notebook cells or CLI."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional, Tuple

import pandas as pd

from .arelf import (
    build_arelf_candidates,
    dynamic_junction_map,
    materialize_arelf_parts,
    selection_overhangs,
)
from .dna import apply_overhang_to_mask, clean_mask, clean_dna
from .import_grasp import compile_target_gap, import_grasp_profile
from .ligation_fidelity import (
    LigationFidelityCalculator,
    fidelity_calculator_for_level,
)
from .objectives import evaluate_design, build_oligos_from_cds
from .optimizer import optimize_library as default_optimize_library
from .optimizer import optimize_coding_sequence as default_optimize_coding_sequence
from .optimizer import simulate_assembled_cds
from .pareto import knee_point, optimize_pareto_overhangs
from .plotting import plot_pareto_front


def flanks_from_parts(parts_df: pd.DataFrame) -> Dict[str, tuple]:
    required = {"oh5_mask_start", "oh3_mask_start"}
    if not required <= set(parts_df.columns):
        raise ValueError(
            "parts table lacks GRASP overhang coordinates; re-import the profile"
        )
    return {
        str(row.part_id): (
            clean_dna(row.oligo_prefix),
            clean_dna(row.oligo_suffix),
            int(row.oh5_mask_start),
            int(row.oh3_mask_start),
        )
        for row in parts_df.itertuples(index=False)
    }


def write_overhangs_into_parts(
    parts_df: pd.DataFrame,
    selection: Dict[str, str],
    junction_map: pd.DataFrame,
    *,
    config: Optional[Mapping] = None,
    codon_data: Optional[Mapping] = None,
) -> pd.DataFrame:
    cut_mode = str(
        (config or {}).get("overhang_redesign", {}).get("cut_mode", "native_fixed")
    ).strip().lower()
    if cut_mode == "movable_arelf":
        return materialize_arelf_parts(
            parts_df,
            selection,
            config=config,
            codon_data=codon_data,
        )

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
        updated_parts = write_overhangs_into_parts(
            parts_df,
            selection,
            junction_map,
            config=local_config,
            codon_data=codon_data,
        )
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
        if str(
            local_config.get("overhang_redesign", {}).get(
                "cut_mode", "native_fixed"
            )
        ).strip().lower() == "movable_arelf":
            return (
                cds_by_part,
                aa_by_part,
                flanks_from_parts(updated_parts),
                dynamic_junction_map(updated_parts),
            )
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

    if not config.get("overhang_redesign", {}).get("enabled", False):
        log("Overhang redesign disabled — returning parts unchanged.")
        return pd.DataFrame(), {}, parts

    if parts is None or len(parts) == 0:
        raise ValueError("parts is empty — import GRASP / load parts.csv first")
    if not codon_data:
        raise ValueError("codon_data empty — Apply codon table first")

    junction_map_file = input_dir / "junction_map.csv"
    junction_map = pd.read_csv(junction_map_file)
    redesign_cfg = config.get("overhang_redesign", {})
    movable_arelf = (
        str(redesign_cfg.get("cut_mode", "native_fixed")).strip().lower()
        == "movable_arelf"
    )
    if movable_arelf:
        offsets = redesign_cfg.get("allowed_arelf_offsets_nt", range(12))
        overhang_candidates = build_arelf_candidates(
            codon_data,
            offsets=offsets,
        )
        overhang_candidates.to_csv(
            output_dir / "arelf_overhang_candidates.csv", index=False
        )
    else:
        candidate_file = input_dir / "overhang_candidates.csv"
        overhang_candidates = pd.read_csv(candidate_file)

    missing = sorted(set(junction_map["part_id"]) - set(parts["part_id"]))
    if missing:
        raise ValueError(f"junction_map unknown part_ids: {missing}")

    lig = config["ligation"]
    calc = fidelity or fidelity_calculator_for_level(
        config,
        lig.get("redesign_level", "level0"),
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

    updated = write_overhangs_into_parts(
        parts,
        selected,
        junction_map,
        config=config,
        codon_data=codon_data,
    )
    front.to_csv(output_dir / "pareto_front.csv", index=False)
    pd.DataFrame([chosen]).to_csv(output_dir / "selected_overhangs.csv", index=False)
    updated.to_csv(output_dir / "parts_with_redesigned_junctions.csv", index=False)
    log(f"Wrote masks → {output_dir / 'parts_with_redesigned_junctions.csv'}")
    return front, selected, updated


def write_oligo_fasta(
    library: pd.DataFrame,
    path: Path,
    *,
    id_col: str = "optimized_part_id",
    seq_col: str = "oligo_sequence_5to3",
) -> Path:
    """Write optimized oligo sequences as FASTA (5′→3′)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if seq_col not in library.columns:
        raise ValueError(f"Library missing {seq_col!r} — cannot write oligo FASTA")
    with open(path, "w") as handle:
        for row in library.itertuples(index=False):
            seq = clean_dna(getattr(row, seq_col))
            oid = str(getattr(row, id_col, "oligo"))
            length = getattr(row, "oligo_length", len(seq))
            qc = getattr(row, "qc_passed", "")
            handle.write(f">{oid}|length={length}|qc={qc}\n{seq}\n")
    return path


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
    # Preserve dynamic cut coordinates with each optimized version.  These are
    # required to release order fragments and reconstruct assembled binders when
    # ARELF cut positions differ from the deposited GenBank coordinates.
    coordinate_columns = [
        column
        for column in (
            "oh5",
            "oh3",
            "oh5_coding_site_5to3",
            "oh3_coding_site_5to3",
            "five_prime_end_overhang",
            "three_prime_end_overhang",
            "overhang_notation",
            "oh5_mask_start",
            "oh3_mask_start",
            "oh5_junction",
            "oh3_junction",
            "oh5_arelf_offset_nt",
            "oh3_arelf_offset_nt",
            "cut_mode",
            "full_window_start_nt",
            "full_window_end_nt",
        )
        if column in parts.columns
    ]
    if coordinate_columns:
        metadata = parts[["part_id", *coordinate_columns]].drop_duplicates("part_id")
        stale = [column for column in coordinate_columns if column in library.columns]
        if stale:
            library = library.drop(columns=stale)
        library = library.merge(metadata, on="part_id", how="left", validate="many_to_one")
    out = output_dir / "optimized_library.csv"
    library.to_csv(out, index=False)
    fasta = write_oligo_fasta(library, output_dir / "optimized_grasp_oligos.fasta")
    log(f"Done — {len(library)} sequences → {out}")
    log(f"Oligo FASTA → {fasta}")
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


def _unpack_cds_optimization(result, fallback_flanks, fallback_junction_map):
    """Normalize fixed-cut and movable-cut optimizer callback results."""
    if len(result) == 2:
        cds_by_part, aa_by_part = result
        return cds_by_part, aa_by_part, fallback_flanks, fallback_junction_map
    if len(result) == 4:
        return result
    raise ValueError(
        "CDS optimizer must return (cds, aa) or "
        "(cds, aa, dynamic_flanks, dynamic_junction_map)"
    )


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
    Re-score every Pareto member with the *same* method after anneal.

    Default: greedy CDS + full-oligo synthesis for all rows (comparable).
    ``deep_all=True`` fully anneals every member (slow, still uniform).

    If ``optimized_library`` matches ``selected_overhangs``, also attach
    ``full_anneal_*`` diagnostic columns for that row only — those are not
    used for knee selection (avoids mixed score sources on the plot).
    """
    if front is None or front.empty:
        raise ValueError("Pareto front is empty")
    if optimized_library is None or len(optimized_library) == 0:
        raise ValueError(
            "optimized_library is required — run full CDS anneal before plotting"
        )

    lig = config.get("ligation", {})
    calc = fidelity or fidelity_calculator_for_level(
        config,
        lig.get("redesign_level", "level0"),
        min_efficiency=lig.get("min_efficiency", 0.25),
        min_fidelity=lig.get("min_fidelity", 0.9),
    )
    flanks = flanks_from_parts(parts)

    annealed_key = None
    if selected_overhangs:
        annealed_key = ";".join(
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

    # Optional diagnostic: full-anneal scores for the library that was actually run
    full_anneal_scores = None
    if annealed_key is not None:
        cds_by_part, aa_by_part = _cds_maps_from_library(optimized_library)
        annealed_flanks = (
            flanks_from_parts(optimized_library)
            if {"oh5_mask_start", "oh3_mask_start"} <= set(optimized_library.columns)
            else flanks
        )
        oligos_by_part = None
        if "oligo_sequence_5to3" in optimized_library.columns:
            oligos_by_part = {}
            for lib_row in optimized_library.itertuples(index=False):
                pid = str(lib_row.part_id)
                if pid in oligos_by_part and getattr(lib_row, "version", 1) != 1:
                    continue
                oligos_by_part[pid] = clean_dna(lib_row.oligo_sequence_5to3)
        if oligos_by_part is None:
            oligos_by_part = build_oligos_from_cds(
                cds_by_part, annealed_flanks, config=config
            )
        part_ids = sorted(cds_by_part)
        sel = parse_overhang_selection(annealed_key)
        full_anneal_scores = evaluate_design(
            overhangs=selection_overhangs(sel),
            cds_sequences=[cds_by_part[p] for p in part_ids],
            aa_sequences=[aa_by_part[p] for p in part_ids],
            codon_data=codon_data,
            config=config,
            fidelity_calculator=calc,
            cds_by_part=cds_by_part,
            aa_by_part=aa_by_part,
            oligos_by_part=oligos_by_part,
            flanks_by_part=annealed_flanks,
            junction_map=junction_map,
        )

    rows = []
    n = len(front)
    source = "full_anneal" if deep_all else "greedy_cds"
    log(
        f"Re-scoring {n} Pareto members uniformly ({source}, full-oligo synthesis)…"
    )
    for i, row in enumerate(front.itertuples(index=False), start=1):
        overhangs_text = str(row.overhangs)
        selection = parse_overhang_selection(overhangs_text)

        if deep_all and deep_optimize is not None:
            log(f"[{i}/{n}] Full-anneal scoring overhang set…")
            optimized = deep_optimize(selection)
        else:
            log(f"[{i}/{n}] Scoring overhang set (greedy CDS + full oligo)…")
            optimized = greedy_optimize(selection)

        cds_by_part, aa_by_part, point_flanks, point_junction_map = (
            _unpack_cds_optimization(optimized, flanks, junction_map)
        )

        oligos_by_part = build_oligos_from_cds(
            cds_by_part, point_flanks, config=config
        )
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
            oligos_by_part=oligos_by_part,
            flanks_by_part=point_flanks,
            junction_map=point_junction_map,
        )
        entry = {
            "overhangs": overhangs_text,
            "ligation_fidelity": scores.ligation_fidelity,
            "codon_optimality": scores.codon_optimality,
            "synthesis": scores.synthesis,
            "fidelity_beam": getattr(row, "fidelity_beam", scores.ligation_fidelity),
            "synthesis_scope": "full_oligo",
            "score_source": source,
            "was_annealed_library": bool(
                annealed_key is not None and overhangs_text == annealed_key
            ),
        }
        if (
            full_anneal_scores is not None
            and annealed_key is not None
            and overhangs_text == annealed_key
        ):
            entry["full_anneal_codon"] = full_anneal_scores.codon_optimality
            entry["full_anneal_synthesis"] = full_anneal_scores.synthesis
            entry["full_anneal_fidelity"] = full_anneal_scores.ligation_fidelity
        rows.append(entry)

    rescored = pd.DataFrame(rows)
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / "pareto_front_after_anneal.csv"
        rescored.to_csv(path, index=False)
        log(f"Wrote post-anneal front → {path}")
    return rescored


def select_from_rescored_front(
    front: pd.DataFrame,
    *,
    selection_mode: str = "knee",
) -> Tuple[pd.Series, Dict[str, str]]:
    """Pick best overhang set on a uniformly scored front."""
    if front is None or front.empty:
        raise ValueError("Empty rescored Pareto front")
    mode = selection_mode or "knee"
    if mode == "max_fidelity":
        chosen = front.sort_values("ligation_fidelity", ascending=False).iloc[0]
    elif mode == "knee":
        chosen = knee_point(front)
    else:
        match = front[front["overhangs"] == mode]
        if match.empty:
            raise ValueError(f"No Pareto point matches selection={mode!r}")
        chosen = match.iloc[0]
    selected = parse_overhang_selection(str(chosen["overhangs"]))
    return chosen, selected


def load_and_validate_parts(path: Path | str) -> pd.DataFrame:
    """Load parts.csv and check mask length vs AA (skips REPLACE placeholders)."""
    parts = pd.read_csv(path).fillna("")
    required = {
        "part_id",
        "aa_sequence",
        "coding_mask",
        "oligo_prefix",
        "oligo_suffix",
        "oh5_mask_start",
        "oh3_mask_start",
    }
    missing = required - set(parts.columns)
    if missing:
        raise ValueError(f"parts.csv missing columns: {missing}")
    if parts["part_id"].duplicated().any():
        duplicates = parts.loc[parts["part_id"].duplicated(), "part_id"].tolist()
        raise ValueError(f"Duplicate part_id values: {duplicates}")

    for row in parts.itertuples(index=False):
        aa = str(row.aa_sequence).upper().replace(" ", "")
        mask = clean_mask(row.coding_mask)
        if "REPLACE" in aa or "REPLACE" in mask:
            continue
        if len(mask) != 3 * len(aa):
            raise ValueError(
                f"{row.part_id}: mask length {len(mask)}, expected {3 * len(aa)}."
            )
    return parts


def ensure_grasp_imported(
    *,
    profile_genbank_dir: Path,
    input_dir: Path,
    force: bool = False,
    log: Callable[[str], None] = print,
) -> Dict[str, pd.DataFrame]:
    """Import GenBank profile into input/, or reuse cached CSVs."""
    input_dir = Path(input_dir)
    parts_path = input_dir / "parts.csv"
    cached_text = parts_path.read_text() if parts_path.exists() else ""
    cached_lines = cached_text.splitlines()
    cached_header = cached_lines[0] if cached_lines else ""
    full_path = input_dir / "parts_full.csv"
    full_lines = full_path.read_text().splitlines() if full_path.exists() else []
    full_header = full_lines[0] if full_lines else ""
    cached_ok = (
        parts_path.exists()
        and "REPLACE" not in cached_text[:500]
        and "oh5_mask_start" in cached_header
        and "oh3_mask_start" in cached_header
        and "five_prime_end_overhang" in cached_header
        and "three_prime_end_overhang" in cached_header
        and "five_prime_end_overhang" in full_header
        and "three_prime_end_overhang" in full_header
        and (input_dir / "junction_map.csv").exists()
        and (input_dir / "overhang_candidates.csv").exists()
    )
    if cached_ok and not force:
        log(f"Using cached import under {input_dir}")
        return {
            "parts": load_and_validate_parts(parts_path),
            "parts_full": pd.read_csv(input_dir / "parts_full.csv")
            if (input_dir / "parts_full.csv").exists()
            else pd.DataFrame(),
            "junction_map": pd.read_csv(input_dir / "junction_map.csv"),
            "overhang_candidates": pd.read_csv(input_dir / "overhang_candidates.csv"),
            "target_map": pd.read_csv(input_dir / "target_map.csv")
            if (input_dir / "target_map.csv").exists()
            else pd.DataFrame(),
        }

    log(f"Importing GRASP GenBank from {profile_genbank_dir} …")
    return import_grasp_profile(profile_genbank_dir, input_dir)


def export_optimized_library(
    library: pd.DataFrame,
    output_dir: Path,
    *,
    selected_overhangs: Optional[Mapping[str, str]] = None,
) -> Dict[str, Path]:
    """Write CSV / FASTA / Excel for the annealed combinatorial library."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = library.copy()
    if selected_overhangs and "selected_overhangs" not in out.columns:
        tag = ";".join(f"{k}={v}" for k, v in sorted(selected_overhangs.items()))
        out["selected_overhangs"] = tag

    csv_path = output_dir / "optimized_grasp_oligos.csv"
    fasta_path = output_dir / "optimized_grasp_oligos.fasta"
    xlsx_path = output_dir / "optimized_grasp_library.xlsx"

    out.to_csv(csv_path, index=False)
    write_oligo_fasta(out, fasta_path)

    qc_cols = [
        c
        for c in [
            "optimized_part_id",
            "translation_verified",
            "mask_verified",
            "codon_score",
            "cds_gc",
            "cds_repeat_penalty",
            "oligo_warnings",
            "oligo_failures",
            "qc_passed",
            "selected_overhangs",
        ]
        if c in out.columns
    ]
    with pd.ExcelWriter(xlsx_path) as writer:
        out.to_excel(writer, sheet_name="oligos", index=False)
        if qc_cols:
            out[qc_cols].to_excel(writer, sheet_name="QC", index=False)

    return {"csv": csv_path, "fasta": fasta_path, "xlsx": xlsx_path}


def compile_and_assemble_target(
    *,
    target_rna: str,
    optimized_library: pd.DataFrame,
    config: Mapping[str, Any],
    input_dir: Path,
    output_dir: Path,
    architecture: str = "9S",
    five_prime_fusion_site: str = "AGGT",
    codon_data: Optional[Dict] = None,
) -> Dict[str, Any]:
    """GAP-compile target RNA against the annealed library and stitch CDS."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plan = compile_target_gap(
        target_rna,
        architecture=architecture,
        five_prime_fusion_site=five_prime_fusion_site,
    )
    parts_full = None
    if {"part_id", "oh5_mask_start", "oh3_mask_start"} <= set(
        optimized_library.columns
    ):
        # Movable ARELF cuts travel with the optimized library.  Prefer those
        # coordinates over the fixed GenBank profile used as the import source.
        parts_full = optimized_library.drop_duplicates("part_id").copy()
    else:
        parts_full_path = input_dir / "parts_full.csv"
        if parts_full_path.exists():
            parts_full = pd.read_csv(parts_full_path)

    assembled = simulate_assembled_cds(
        plan,
        optimized_library,
        parts_full=parts_full,
        genetic_code=int(config.get("genetic_code", 1)),
        config=config,
        codon_data=codon_data,
    )
    rna = str(target_rna).upper().replace("T", "U")
    plan_path = output_dir / f"assembly_plan_{rna}.csv"
    fasta_path = output_dir / f"assembled_{rna}.fasta"
    oligo_fasta_path = output_dir / f"oligos_{rna}.fasta"
    oligo_csv_path = output_dir / f"oligos_{rna}.csv"
    plan.to_csv(plan_path, index=False)
    with open(fasta_path, "w") as handle:
        handle.write(f">GRASP_{rna}|assembled_CDS|length={len(assembled['assembled_cds'])}\n")
        handle.write(assembled["assembled_cds"] + "\n")

    # Optimized oligos for this target, in GAP assembly order
    lib_by_id = optimized_library.set_index("optimized_part_id", drop=False)
    oligo_rows = []
    for slot, row in enumerate(plan.itertuples(index=False), start=1):
        oid = str(row.optimized_part_id)
        if oid not in lib_by_id.index:
            raise ValueError(f"Missing optimized oligo for {oid}")
        lib_row = lib_by_id.loc[oid]
        if isinstance(lib_row, pd.DataFrame):
            lib_row = lib_row.iloc[0]
        oligo_rows.append(
            {
                "assembly_slot": slot,
                "part_id": str(row.part_id),
                "optimized_part_id": oid,
                "assembly_group": getattr(row, "assembly_group", ""),
                "assembly_order": getattr(row, "assembly_order", ""),
                "oligo_length": int(lib_row.oligo_length),
                "qc_passed": bool(lib_row.qc_passed),
                "oligo_sequence_5to3": clean_dna(lib_row.oligo_sequence_5to3),
            }
        )
    oligos = pd.DataFrame(oligo_rows)
    oligos.to_csv(oligo_csv_path, index=False)
    with open(oligo_fasta_path, "w") as handle:
        for row in oligos.itertuples(index=False):
            handle.write(
                f">{rna}|slot{row.assembly_slot}|{row.optimized_part_id}"
                f"|length={row.oligo_length}|qc={row.qc_passed}\n"
                f"{row.oligo_sequence_5to3}\n"
            )

    return {
        "assembly_plan": plan,
        "assembled": assembled,
        "oligos": oligos,
        "plan_csv": plan_path,
        "assembled_fasta": fasta_path,
        "oligo_fasta": oligo_fasta_path,
        "oligo_csv": oligo_csv_path,
        "target_rna": rna,
    }


def plot_library_pareto_after_anneal(
    *,
    front: pd.DataFrame,
    parts: pd.DataFrame,
    junction_map: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    optimized_library: pd.DataFrame,
    selected_overhangs: Optional[Dict[str, str]] = None,
    fidelity: Optional[LigationFidelityCalculator] = None,
    deep_all: bool = False,
    output_dir: Optional[Path] = None,
    log: Callable[[str], None] = print,
):
    """
    Uniformly rescore the front, re-select the best point, and plot it.

    The star is the post-anneal knee / max-fidelity choice on comparable scores —
    not the locked pre-anneal selection.
    """
    rescored = rescore_pareto_front_after_anneal(
        front,
        parts=parts,
        junction_map=junction_map,
        codon_data=codon_data,
        config=config,
        optimize_coding_sequence=default_optimize_coding_sequence,
        selected_overhangs=selected_overhangs,
        optimized_library=optimized_library,
        fidelity=fidelity,
        deep_all=deep_all,
        output_dir=output_dir,
        log=log,
    )
    mode = config.get("overhang_redesign", {}).get("selection", "knee")
    chosen, selected = select_from_rescored_front(rescored, selection_mode=mode)

    if output_dir is not None:
        output_dir = Path(output_dir)
        pd.DataFrame([chosen]).to_csv(
            output_dir / "selected_overhangs_after_anneal.csv", index=False
        )
        # Keep a simple dict dump for reload
        pd.DataFrame(
            [{"junction": k, "overhang": v} for k, v in sorted(selected.items())]
        ).to_csv(output_dir / "selected_overhangs_after_anneal_map.csv", index=False)

    pre_key = (
        ";".join(f"{k}={v}" for k, v in sorted(selected_overhangs.items()))
        if selected_overhangs
        else None
    )
    post_key = str(chosen["overhangs"])
    if pre_key and pre_key != post_key:
        log(
            "⚠ Post-anneal selection differs from the overhang set used for "
            "library anneal.\n"
            f"  annealed (pre): {pre_key}\n"
            f"  best now (post): {post_key}\n"
            "  Re-run Anneal library if you want oligos for the new choice."
        )
    else:
        log(f"Post-anneal selection ({mode}): {post_key}")

    save_path = Path(output_dir) / "pareto_front.png" if output_dir else None
    scope = "full anneal" if deep_all else "greedy CDS"
    fig, ax, plotted = plot_pareto_front(
        rescored,
        selected_overhangs=selected,
        selection_mode=mode,
        save_path=save_path,
        title=(
            f"Pareto front · post-anneal ({scope}, full oligo) · "
            f"selected = {mode}"
        ),
    )
    return {
        "front": rescored,
        "figure": fig,
        "axes": ax,
        "chosen": chosen,
        "selected_overhangs": selected,
        "pre_anneal_selected_overhangs": selected_overhangs,
        "selection_changed": bool(pre_key and pre_key != post_key),
    }


def run_library_redesign_and_anneal(
    *,
    parts: pd.DataFrame,
    codon_data: Dict,
    config: Dict[str, Any],
    input_dir: Path,
    output_dir: Path,
    seed: int = 42,
    fidelity: Optional[LigationFidelityCalculator] = None,
    log: Callable[[str], None] = print,
) -> Dict[str, Any]:
    """Overhang redesign (optional) → full CDS anneal for all modules."""
    front, selected, updated = run_overhang_redesign(
        parts=parts,
        codon_data=codon_data,
        config=config,
        optimize_coding_sequence=default_optimize_coding_sequence,
        input_dir=input_dir,
        output_dir=output_dir,
        seed=seed,
        fidelity=fidelity,
        log=log,
    )
    source = updated if updated is not None and len(updated) else parts
    library = run_library_optimize(
        parts=source,
        codon_data=codon_data,
        config=config,
        optimize_library=default_optimize_library,
        output_dir=output_dir,
        log=log,
    )
    if selected:
        tag = ";".join(f"{k}={v}" for k, v in sorted(selected.items()))
        library = library.copy()
        library["selected_overhangs"] = tag
        library.to_csv(Path(output_dir) / "optimized_library.csv", index=False)
    return {
        "pareto_front": front,
        "selected_overhangs": selected,
        "parts": source,
        "optimized_library": library,
    }
