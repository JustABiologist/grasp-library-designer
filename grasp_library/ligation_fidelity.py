"""Wrapper around GGAssembler/dawdlib GGData.reaction_fidelity."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, Tuple

import sys

_ROOT = Path(__file__).resolve().parents[1]
_VENDOR = _ROOT / "third_party"
if str(_VENDOR) not in sys.path:
    sys.path.insert(0, str(_VENDOR))

from dawdlib_golden_gate import GGData  # noqa: E402


class LigationFidelityCalculator:
    """
    GGAssembler ligation-fidelity engine (Potapov matrices via GGData).

    `reaction_fidelity` returns (directional_fwd, directional_rev, bidirectional).
    By default we report the directional NEB-style product used in GGAssembler
    path selection (`reaction_fidelity(...)[0]`).
    """

    DEFAULT_TEMPERATURE = 25
    DEFAULT_HOURS = 18
    _RESOURCES = (
        Path(__file__).resolve().parents[1]
        / "third_party"
        / "dawdlib_golden_gate"
        / "resources"
    )

    def __init__(
        self,
        temperature: int = DEFAULT_TEMPERATURE,
        hours: int = DEFAULT_HOURS,
        ligation_table: str | Path | None = None,
        min_efficiency: float = GGData.MIN_EFFICIENCY,
        min_fidelity: float = GGData.MIN_FIDELITY,
    ) -> None:
        self.temperature = temperature
        self.hours = hours
        self.ligation_table = ligation_table
        # Potapov (hours, temp) keys are limited; bootstrap with a valid pair
        # then optionally swap in an enzyme-specific CSV.
        try:
            self.ggdata = GGData(
                temperature=temperature,
                hours=hours,
                min_efficiency=min_efficiency,
                min_fidelity=min_fidelity,
            )
        except ValueError:
            self.ggdata = GGData(
                temperature=self.DEFAULT_TEMPERATURE,
                hours=self.DEFAULT_HOURS,
                min_efficiency=min_efficiency,
                min_fidelity=min_fidelity,
            )
        if ligation_table is not None:
            table_path = self.resolve_table(ligation_table)
            self.ggdata.set_default_df(str(table_path))
            self.ggdata.init()

    @classmethod
    def resolve_table(cls, ligation_table: str | Path) -> Path:
        path = Path(ligation_table)
        if path.exists():
            return path
        candidate = cls._RESOURCES / Path(ligation_table).name
        if candidate.exists():
            return candidate
        raise FileNotFoundError(
            f"Ligation table not found: {ligation_table} (also tried {candidate})"
        )

    def pair_fidelity(self, overhang_a: str, overhang_b: str) -> float:
        return float(self.ggdata.overhangs_fidelity(overhang_a.upper(), overhang_b.upper()))

    def reaction_fidelity(
        self,
        overhangs: Sequence[str],
    ) -> Tuple[float, float, float]:
        cleaned = [str(o).upper().replace("U", "T") for o in overhangs]
        neb, fwd, rev = self.ggdata.reaction_fidelity(*cleaned)
        return float(neb), float(fwd), float(rev)

    def set_fidelity(self, overhangs: Iterable[str]) -> float:
        """Primary scalar used in Pareto: GGAssembler directional fidelity."""
        neb, _, _ = self.reaction_fidelity(list(overhangs))
        return neb

    def efficient_overhangs(self, filter_gc: bool = False) -> list[str]:
        return list(self.ggdata.filter_self_binding_gates(filter_gc=filter_gc))
