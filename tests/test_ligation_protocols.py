from __future__ import annotations

import hashlib
import math

import pandas as pd
import pytest

from grasp_library.control_panel import build_default_config
from grasp_library.ligation_fidelity import (
    LigationFidelityCalculator,
    fidelity_calculator_for_level,
)
from grasp_library.synthesis_vendors import (
    GRASP_LIGATION_BY_LEVEL,
    LEVEL0_LIGATION,
    LEVEL1_LIGATION,
    LEVEL_MINUS1_LIGATION,
    LIGATION_TABLES,
    apply_ligation_table_to_config,
    ligation_table_names,
    redesign_ligation_table_names,
)


GRASP_LEVEL0_OVERHANGS = ["CTCA", "ACTC", "AAGA", "GCAC", "TGAA", "CGAG"]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.parametrize(
    ("temperature", "hours", "expected"),
    [
        (25, 18, 0.8430188400019161),  # NEB viewer displays 84%
        (37, 18, 0.9049667570689219),  # NEB viewer displays 90%
        (25, 1, 0.8126849298956424),  # NEB viewer displays 81%
    ],
)
def test_potapov_score_matches_neb_viewer(
    temperature: int, hours: int, expected: float
) -> None:
    calc = LigationFidelityCalculator(temperature=temperature, hours=hours)

    assert calc.set_fidelity(GRASP_LEVEL0_OVERHANGS) == pytest.approx(expected)


def test_set_score_is_geometric_mean_and_orientation_invariant() -> None:
    calc = LigationFidelityCalculator(ligation_table="BbsI-HF.csv")
    forward, reverse, _ = calc.reaction_fidelity(GRASP_LEVEL0_OVERHANGS)
    reverse_complements = [
        _reverse_complement(overhang) for overhang in GRASP_LEVEL0_OVERHANGS
    ]

    assert calc.set_fidelity(GRASP_LEVEL0_OVERHANGS) == pytest.approx(
        math.sqrt(forward * reverse)
    )
    assert calc.set_fidelity(GRASP_LEVEL0_OVERHANGS) == pytest.approx(
        0.7516829347822497
    )
    assert calc.set_fidelity(reverse_complements) == pytest.approx(
        calc.set_fidelity(GRASP_LEVEL0_OVERHANGS)
    )


def test_fixed_grasp_overhangs_are_not_removed_by_candidate_filter() -> None:
    # TGAA is below the default efficiency cutoff in this dataset. It remains
    # part of the physical reaction score while staying out of design choices.
    calc = LigationFidelityCalculator(temperature=37, hours=18)

    assert "TGAA" not in calc.efficient_overhangs()
    assert calc.set_fidelity(GRASP_LEVEL0_OVERHANGS) == pytest.approx(
        0.9049667570689219
    )


def test_unsupported_static_temperature_is_not_silently_substituted() -> None:
    with pytest.raises(ValueError, match="No data was found"):
        LigationFidelityCalculator(temperature=16, hours=1)


def test_three_base_table_is_rejected_for_grasp_scoring() -> None:
    with pytest.raises(ValueError, match="256 x 256 four-base"):
        LigationFidelityCalculator(ligation_table="SapI.csv")


def test_dashboard_exposes_level_matched_four_base_protocols() -> None:
    names = ligation_table_names()
    redesign = redesign_ligation_table_names()

    assert names[:3] == [
        LEVEL_MINUS1_LIGATION,
        LEVEL0_LIGATION,
        LEVEL1_LIGATION,
    ]
    assert redesign[0] == LEVEL0_LIGATION
    assert LEVEL_MINUS1_LIGATION not in redesign
    assert LEVEL1_LIGATION not in redesign
    assert all("SapI" not in name for name in names)
    assert all("constant 37" not in name for name in names)
    assert all("constant 42" not in name for name in names)
    assert "T4 ligase only · 1 h · 37 °C (Potapov 2018)" in names


@pytest.mark.parametrize(
    ("name", "table", "restriction_enzyme", "level"),
    [
        (LEVEL_MINUS1_LIGATION, "BsaI-HFv2.csv", "BsaI-HFv2", "level_minus1"),
        (LEVEL0_LIGATION, "BbsI-HF.csv", "BbsI-HF", "level0"),
        (LEVEL1_LIGATION, "BsaI-HFv2.csv", "BsaI-HFv2", "level1"),
    ],
)
def test_cycling_proxy_metadata_is_not_relabelled_as_static_ligation(
    name: str, table: str, restriction_enzyme: str, level: str
) -> None:
    metadata = LIGATION_TABLES[name]

    assert metadata["table"] == table
    assert metadata["temperature"] is None
    assert metadata["hours"] is None
    assert metadata["assay_kind"] == "golden_gate_cycling"
    assert metadata["overhang_length"] == 4
    assert metadata["restriction_enzyme"] == restriction_enzyme
    assert metadata["cloning_level"] == level
    assert metadata["cycles"] == 30
    assert metadata["steps"] == [
        {"temperature_c": 37, "minutes": 5},
        {"temperature_c": 16, "minutes": 5},
    ]
    assert metadata["grasp_status"] == "proxy"
    assert metadata["grasp_reference_protocol"]["cycles"] == 26
    assert metadata["grasp_reference_protocol"]["steps"] == [
        {"temperature_c": 37, "minutes": 3},
        {"temperature_c": 16, "minutes": 4},
    ]
    assert "supplementary" in metadata["source_url"]


def test_selected_protocol_metadata_is_preserved_in_config() -> None:
    updated = apply_ligation_table_to_config({"ligation": {}}, LEVEL0_LIGATION)
    ligation = updated["ligation"]

    assert ligation["temperature"] is None
    assert ligation["hours"] is None
    assert ligation["ligation_table"] == "BbsI-HF.csv"
    assert ligation["protocol_metadata"]["assay_kind"] == "golden_gate_cycling"
    assert ligation["protocol_metadata"]["proxy_for"].startswith("GRASP Level 0")
    assert ligation["by_level"]["level_minus1"]["ligation_table"] == "BsaI-HFv2.csv"
    assert ligation["by_level"]["level0"]["ligation_table"] == "BbsI-HF.csv"
    assert ligation["by_level"]["level1"]["ligation_table"] == "BsaI-HFv2.csv"
    assert GRASP_LIGATION_BY_LEVEL["level_minus1"] == LEVEL_MINUS1_LIGATION


def test_selecting_bsai_level_does_not_retarget_level0_redesign() -> None:
    updated = apply_ligation_table_to_config({"ligation": {}}, LEVEL1_LIGATION)

    assert updated["ligation"]["ligation_table"] == "BbsI-HF.csv"
    assert updated["ligation"]["table_name"] == LEVEL0_LIGATION
    assert updated["ligation"]["by_level"]["level1"]["ligation_table"] == "BsaI-HFv2.csv"


def test_fidelity_calculator_for_level_uses_enzyme_matched_tables(tmp_path) -> None:
    cfg = build_default_config(tmp_path)
    level0 = fidelity_calculator_for_level(cfg, "level0")
    level1 = fidelity_calculator_for_level(cfg, "level1")
    entry = fidelity_calculator_for_level(cfg, "level_minus1")

    assert str(level0.ligation_table).endswith("BbsI-HF.csv")
    assert str(level1.ligation_table).endswith("BsaI-HFv2.csv")
    assert str(entry.ligation_table).endswith("BsaI-HFv2.csv")


@pytest.mark.parametrize(
    ("resource", "expected_sum", "expected_sha256"),
    [
        (
            "BbsI-HF.csv",
            189_282,
            "eed3c64c09ca372ca647ef2a7dc4c478670eefa1bdf4c78809d98923341b1dba",
        ),
        (
            "BsaI-HFv2.csv",
            203_364,
            "6cb5a033866cf8241f62d4049f4e29d2f155e561ef236c9428c5fe102d5e56c7",
        ),
        (
            "FileS_T4_01h_37C.csv",
            682_150,
            "f04a0fe179cb1fc1f92326357da09f8e76b7956a005ad59148678d4e74bc8cf3",
        ),
    ],
)
def test_primary_source_matrices_are_complete_and_canonical(
    resource: str, expected_sum: int, expected_sha256: str
) -> None:
    path = LigationFidelityCalculator.resolve_table(resource)
    matrix = pd.read_csv(path, index_col=0)
    expected_labels = [
        "".join((a, b, c, d))
        for a in "ACGT"
        for b in "ACGT"
        for c in "ACGT"
        for d in "ACGT"
    ]

    assert matrix.shape == (256, 256)
    assert list(matrix.columns) == expected_labels
    assert list(matrix.index) == [
        _reverse_complement(label) for label in expected_labels
    ]
    assert int(matrix.to_numpy().sum()) == expected_sum
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_sha256
