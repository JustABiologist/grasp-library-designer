"""Pareto front visualization — 3D trade-off space with lab-notebook styling."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator

from .pareto import knee_point

# Match notebook_ui palette
INK = "#1c2b24"
INK_SOFT = "#3d5248"
PAPER = "#f4f1ea"
ACCENT = "#0f6b4c"
STAR = "#d4551a"
POINT = "#2f6f9f"
GRID = "#b7c2b8"


def _chosen_row(
    front: pd.DataFrame,
    selected_overhangs: Optional[Mapping[str, str]] = None,
    selection_mode: str = "knee",
) -> pd.Series:
    if selected_overhangs:
        key = ";".join(f"{k}={v}" for k, v in sorted(selected_overhangs.items()))
        match = front[front["overhangs"] == key]
        if not match.empty:
            return match.iloc[0]
    if selection_mode == "knee":
        return knee_point(front)
    return front.sort_values("ligation_fidelity", ascending=False).iloc[0]


def _pad(lo: float, hi: float, frac: float = 0.14, floor: float = 1e-9) -> Tuple[float, float]:
    span = max(hi - lo, floor)
    return lo - frac * span, hi + frac * span


def _fmt_fidelity(value: float, _pos=None) -> str:
    if abs(value - 1.0) < 5e-7:
        return "1.0"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _fmt_signed(value: float, _pos=None) -> str:
    if abs(value) < 1e-10:
        return "0"
    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.4f}".rstrip("0").rstrip(".")


def plot_pareto_front(
    front: pd.DataFrame,
    *,
    selected_overhangs: Optional[Mapping[str, str]] = None,
    selection_mode: str = "knee",
    save_path: Optional[Union[str, Path]] = None,
    title: str = "Pareto front",
):
    """
    3D scatter of ligation × junction-codon × junction-synthesis.

    Styled for the GRASP notebook (paper / forest accent). Stems drop to the
    fidelity–codon floor so depth reads clearly; selected point is a star.
    """
    if front is None or front.empty:
        raise ValueError("Empty Pareto front.")

    df = front.copy()
    chosen = _chosen_row(df, selected_overhangs, selection_mode)

    # Drop hard synthesis failures from the axis (keep them listed in the table)
    plot_df = df[df["synthesis"].astype(float) > -1e3].copy()
    if plot_df.empty:
        plot_df = df.copy()

    fid = plot_df["ligation_fidelity"].astype(float).to_numpy()
    codon = plot_df["codon_optimality"].astype(float).to_numpy()
    synth = plot_df["synthesis"].astype(float).to_numpy()

    cx = float(chosen["ligation_fidelity"])
    cy = float(chosen["codon_optimality"])
    cz = float(chosen["synthesis"])

    fig = plt.figure(figsize=(9.2, 6.8), facecolor=PAPER, dpi=140)
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.set_facecolor(PAPER)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.90)

    # Soft panes + thin grid (the “old plot” silhouette, cleaned up)
    pane = (0.96, 0.95, 0.92, 0.92)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(pane)
        axis.pane.set_edgecolor("#d5ddd6")
        axis.pane.set_alpha(0.95)
        axis._axinfo["grid"]["color"] = (0.72, 0.76, 0.73, 0.55)
        axis._axinfo["grid"]["linewidth"] = 0.7
        axis._axinfo["tick"]["inward_factor"] = 0.0
        axis._axinfo["tick"]["outward_factor"] = 0.35

    ax.view_init(elev=24, azim=-52)
    ax.dist = 11

    z_floor = float(np.min(synth))
    # Stems to the floor — classic readable 3D scatter
    for x, y, z in zip(fid, codon, synth):
        ax.plot(
            [x, x],
            [y, y],
            [z_floor, z],
            color=GRID,
            lw=1.0,
            alpha=0.85,
            zorder=1,
        )
        ax.scatter(
            [x],
            [y],
            [z_floor],
            s=18,
            c=GRID,
            alpha=0.7,
            depthshade=False,
            zorder=1,
        )

    # Front members
    ax.scatter(
        fid,
        codon,
        synth,
        s=110,
        c=POINT,
        depthshade=True,
        edgecolors=INK,
        linewidths=0.55,
        alpha=0.92,
        zorder=3,
        label="Pareto set",
    )
    # Soft halo under points
    ax.scatter(
        fid,
        codon,
        synth,
        s=260,
        c=POINT,
        alpha=0.12,
        depthshade=False,
        linewidths=0,
        zorder=2,
    )

    # Selected
    ax.plot(
        [cx, cx],
        [cy, cy],
        [z_floor, cz],
        color=STAR,
        lw=1.4,
        alpha=0.9,
        zorder=4,
    )
    ax.scatter(
        [cx],
        [cy],
        [cz],
        s=420,
        marker="*",
        c=STAR,
        edgecolors=INK,
        linewidths=0.6,
        depthshade=False,
        zorder=5,
        label="selected",
    )

    ax.set_xlabel("ligation fidelity", labelpad=10, color=INK_SOFT, fontsize=11)
    ax.set_ylabel("codon optimality", labelpad=12, color=INK_SOFT, fontsize=11)
    ax.set_zlabel("oligo synthesis", labelpad=10, color=INK_SOFT, fontsize=11)

    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_fidelity))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_signed))
    ax.zaxis.set_major_formatter(FuncFormatter(_fmt_signed))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.zaxis.set_major_locator(MaxNLocator(5))

    xlim = _pad(float(fid.min()), float(fid.max()), 0.2, 1e-7)
    ylim = _pad(float(codon.min()), float(codon.max()), 0.16, 1e-4)
    zlim = _pad(float(synth.min()), float(synth.max()), 0.22, 1e-4)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.label.set_color(INK_SOFT)
        axis.set_tick_params(colors=INK_SOFT, labelsize=9, pad=2)

    # Title block
    fig.text(
        0.03,
        0.955,
        title,
        fontsize=16,
        color=INK,
        fontname="DejaVu Serif",
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.03,
        0.915,
        "ligation fidelity  ·  full-CDS codon  ·  full-oligo synthesis (prefix + CDS + suffix)",
        fontsize=10,
        color=ACCENT,
        fontname="DejaVu Sans",
        va="top",
    )

    legend = ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, 0.96),
        frameon=True,
        fancybox=False,
        edgecolor="#c9d0c8",
        facecolor="white",
        framealpha=0.92,
        fontsize=9,
    )
    for text in legend.get_texts():
        text.set_color(INK)

    # Score chip for the selected point
    chip = (
        f"selected   fidelity {cx:.6f}   "
        f"codon {cy:.3f}   synthesis {cz:.4f}"
    )
    fig.text(
        0.97,
        0.04,
        chip,
        ha="right",
        va="bottom",
        fontsize=9,
        color=INK_SOFT,
        fontname="DejaVu Sans",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#c9d0c8",
            "alpha": 0.92,
        },
    )

    if save_path is not None:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            path,
            dpi=200,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
            edgecolor="none",
        )

    return fig, ax, chosen
