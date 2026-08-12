"""Wrapper around GGAssembler/dawdlib GGData.reaction_fidelity."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

# Prefer the installed package; fall back to the repo's third_party/ tree.
try:
    import dawdlib_golden_gate  # type: ignore
except ImportError:  # pragma: no cover - editable/source checkout
    _vendor = Path(__file__).resolve().parents[1] / "third_party"
    if _vendor.is_dir() and str(_vendor) not in sys.path:
        sys.path.insert(0, str(_vendor))
    import dawdlib_golden_gate  # type: ignore

from dawdlib_golden_gate import GGData  # noqa: E402


class LigationFidelityCalculator:
    """
    GGAssembler ligation-fidelity engine (Potapov matrices via GGData).

    `reaction_fidelity` returns (directional_fwd, directional_rev, bidirectional).
    `set_fidelity` reports the orientation-invariant geometric mean of the two
    directional products, matching the value reported by NEB's Ligase Fidelity
    Viewer.  These are model scores derived from ligation-count matrices, not a
    vendor acceptance or cloning guarantee.
    """

    # Every GRASP five-part Level -1 reaction is released from pAGM1311 and
    # inserted into pAGM9121 with fixed external BpiI overhangs CTCA/CGAG.  The
    # four internal junctions are shared by CDS1/CDS14/CDS19/CDS2.  The coding
    # interfaces AGGT/GTGA/CACG/CTTC/TTCG belong to the later BsaI/MoClo stage
    # and must not be mixed into this reaction score.
    GRASP_LEVEL0_INTERNAL_JUNCTIONS = (
        "J_ACTC",
        "J_AAGA",
        "J_GCAC",
        "J_TGAA",
    )
    GRASP_LEVEL0_EXTERNAL_OVERHANGS = ("CTCA", "CGAG")
    GRASP_ARCHITECTURE_GROUPS = {
        "9S": ("CDS1", "CDS2"),
        "14S": ("CDS1", "CDS14", "CDS2"),
        "19S": ("CDS1", "CDS14", "CDS19", "CDS2"),
    }

    DEFAULT_TEMPERATURE = 25
    DEFAULT_HOURS = 18
    _RESOURCES = Path(dawdlib_golden_gate.__file__).resolve().parent / "resources"

    def __init__(
        self,
        temperature: int | None = DEFAULT_TEMPERATURE,
        hours: int | None = DEFAULT_HOURS,
        ligation_table: str | Path | None = None,
        min_efficiency: float = GGData.MIN_EFFICIENCY,
        min_fidelity: float = GGData.MIN_FIDELITY,
    ) -> None:
        self.temperature = temperature
        self.hours = hours
        self.ligation_table = ligation_table

        init_kwargs = {
            "min_efficiency": min_efficiency,
            "min_fidelity": min_fidelity,
        }
        score_kwargs = {
            # Reaction scoring must retain every overhang in the empirical
            # matrix.  The user threshold is only a candidate-design filter;
            # applying it to fixed GRASP junctions otherwise causes KeyError or
            # silently changes the denominator of the measured reaction.
            "min_efficiency": -1.0,
            "min_fidelity": min_fidelity,
        }
        if ligation_table is None:
            # Let GGData reject unmeasured static temperature/time pairs.  In
            # particular, never relabel the 25 C/18 h data as a 16 C matrix.
            if temperature is None or hours is None:
                raise ValueError(
                    "Static ligation datasets require measured temperature "
                    "and duration values."
                )
            self.ggdata = GGData(
                temperature=temperature,
                hours=hours,
                **init_kwargs,
            )
            self._scoring_ggdata = GGData(
                temperature=temperature,
                hours=hours,
                **score_kwargs,
            )
        else:
            # GGData requires a built-in dataset at construction, but this
            # bootstrap is fully replaced before either object is exposed.
            self.ggdata = GGData(
                temperature=self.DEFAULT_TEMPERATURE,
                hours=self.DEFAULT_HOURS,
                **init_kwargs,
            )
            self._scoring_ggdata = GGData(
                temperature=self.DEFAULT_TEMPERATURE,
                hours=self.DEFAULT_HOURS,
                **score_kwargs,
            )
            table_path = self.resolve_table(ligation_table)
            for dataset in (self.ggdata, self._scoring_ggdata):
                dataset.set_default_df(str(table_path))
                dataset.init()

        self._validate_four_base_matrix(self._scoring_ggdata.lig_df)

    @staticmethod
    def _validate_four_base_matrix(matrix) -> None:
        """Reject tables that cannot represent GRASP's four-base junctions."""
        rows = [str(label).upper() for label in matrix.index]
        columns = [str(label).upper() for label in matrix.columns]
        alphabet = set("ACGT")
        valid_rows = all(len(label) == 4 and set(label) <= alphabet for label in rows)
        valid_columns = all(
            len(label) == 4 and set(label) <= alphabet for label in columns
        )
        values_are_valid = (
            not matrix.isna().any().any()
            and not (matrix < 0).any().any()
            and float(matrix.to_numpy().sum()) > 0.0
        )
        complement = str.maketrans("ACGT", "TGCA")
        has_watson_crick_support = all(
            float(matrix.loc[label, label.translate(complement)[::-1]]) > 0.0
            for label in rows
        )
        if (
            matrix.shape != (256, 256)
            or not valid_rows
            or not valid_columns
            or set(rows) != set(columns)
            or not values_are_valid
            or not has_watson_crick_support
        ):
            raise ValueError(
                "GRASP fidelity scoring requires a complete 256 x 256 "
                "four-base overhang matrix."
            )

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
        return float(
            self._scoring_ggdata.overhangs_fidelity(
                overhang_a.upper(), overhang_b.upper()
            )
        )

    def reaction_fidelity(
        self,
        overhangs: Sequence[str],
    ) -> Tuple[float, float, float]:
        cleaned = [str(o).upper().replace("U", "T") for o in overhangs]
        directional_fwd, directional_rev, bidirectional = (
            self._scoring_ggdata.reaction_fidelity(*cleaned)
        )
        return (
            float(directional_fwd),
            float(directional_rev),
            float(bidirectional),
        )

    def set_fidelity(self, overhangs: Iterable[str]) -> float:
        """NEB-style, orientation-invariant fidelity for one reaction."""
        directional_fwd, directional_rev, _ = self.reaction_fidelity(list(overhangs))
        return float(math.sqrt(directional_fwd * directional_rev))

    def grasp_level0_reaction_overhangs(
        self,
        junction_overhangs: Mapping[str, str],
        *,
        external_overhangs: Sequence[str] | None = None,
    ) -> list[str]:
        """Return the six overhangs in one configured Level 0 assembly tube."""
        missing = sorted(set(self.GRASP_LEVEL0_INTERNAL_JUNCTIONS) - set(junction_overhangs))
        if missing:
            raise ValueError(
                "Missing GRASP Level 0 junction overhangs: " + ", ".join(missing)
            )
        cleaned = {
            str(junction): str(overhang).upper().replace("U", "T")
            for junction, overhang in junction_overhangs.items()
        }
        outer = tuple(external_overhangs or self.GRASP_LEVEL0_EXTERNAL_OVERHANGS)
        if len(outer) != 2:
            raise ValueError(
                "Level 0 external_overhangs must contain 5′ and 3′ values"
            )
        left, right = (str(value).upper().replace("U", "T") for value in outer)
        return [left] + [cleaned[j] for j in self.GRASP_LEVEL0_INTERNAL_JUNCTIONS] + [right]

    def grasp_first_stage_fidelities(
        self,
        junction_overhangs: Mapping[str, str],
        *,
        architecture: str = "9S",
        external_overhangs: Sequence[str] | None = None,
    ) -> Dict[str, float]:
        """Score each physical pAGM1311-to-pAGM9121 assembly reaction."""
        architecture = str(architecture).strip().upper()
        if architecture not in self.GRASP_ARCHITECTURE_GROUPS:
            raise ValueError(f"Unsupported GRASP architecture {architecture!r}")
        reaction_overhangs = self.grasp_level0_reaction_overhangs(
            junction_overhangs, external_overhangs=external_overhangs
        )
        score = self.set_fidelity(reaction_overhangs)
        return {
            group: score for group in self.GRASP_ARCHITECTURE_GROUPS[architecture]
        }

    def grasp_first_stage_fidelity(
        self,
        junction_overhangs: Mapping[str, str],
        *,
        architecture: str = "9S",
        external_overhangs: Sequence[str] | None = None,
    ) -> float:
        """Representative fidelity of one physical Level 0 reaction.

        Every GRASP block uses the same six junction identities but is assembled,
        transformed, and screened in a separate tube.  Multiplying these scores
        is therefore not the fidelity of any experimental reaction.
        """
        scores = self.grasp_first_stage_fidelities(
            junction_overhangs,
            architecture=architecture,
            external_overhangs=external_overhangs,
        )
        return float(next(iter(scores.values())))

    def grasp_architecture_success_product_estimate(
        self,
        junction_overhangs: Mapping[str, str],
        *,
        architecture: str = "9S",
        external_overhangs: Sequence[str] | None = None,
    ) -> float:
        """Independence/no-screening estimate across all required Level 0 tubes."""
        scores = self.grasp_first_stage_fidelities(
            junction_overhangs,
            architecture=architecture,
            external_overhangs=external_overhangs,
        )
        product = 1.0
        for score in scores.values():
            product *= score
        return float(product)

    # Backwards-compatible names retained for callers of the initial 9S API.
    def grasp_9s_first_stage_fidelities(
        self, junction_overhangs: Mapping[str, str]
    ) -> Dict[str, float]:
        return self.grasp_first_stage_fidelities(junction_overhangs, architecture="9S")

    def grasp_9s_first_stage_fidelity(
        self, junction_overhangs: Mapping[str, str]
    ) -> float:
        return self.grasp_first_stage_fidelity(junction_overhangs, architecture="9S")

    def efficient_overhangs(self, filter_gc: bool = False) -> list[str]:
        return list(self.ggdata.filter_self_binding_gates(filter_gc=filter_gc))
