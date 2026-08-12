from pathlib import Path

import pytest

from grasp_library.import_grasp import compile_target_gap, load_grasp_records


@pytest.mark.parametrize(
    ("architecture", "target", "expected_groups", "boundary_parts"),
    [
        (
            "9S",
            "A" * 9,
            ["CDS1", "CDS2"],
            {0: "1A_5T_AGGT", 4: "1E_LN5T", 5: "2A_LN5T", 9: "2E_LN"},
        ),
        (
            "14S",
            "A" * 14,
            ["CDS1", "CDS14", "CDS2"],
            {
                0: "1A_5T_AGGT",
                4: "14E_LN5T",
                5: "14A_LN5T",
                9: "1E_LN5T",
                10: "2A_LN5T",
                14: "2E_LN",
            },
        ),
        (
            "19S",
            "A" * 19,
            ["CDS1", "CDS14", "CDS19", "CDS2"],
            {
                0: "1A_5T_AGGT",
                4: "14E_LN5T",
                5: "14A_LN5T",
                9: "19E_LN5T",
                10: "19A_LN5T",
                14: "1E_LN5T",
                15: "2A_LN5T",
                19: "2E_LN",
            },
        ),
    ],
)
def test_compile_architecture_uses_five_part_level_zero_blocks(
    architecture, target, expected_groups, boundary_parts
):
    plan = compile_target_gap(target, architecture=architecture)

    assert len(plan) == len(expected_groups) * 5
    assert plan["assembly_slot"].tolist() == list(range(1, len(plan) + 1))
    assert plan["assembly_group"].tolist() == [
        group for group in expected_groups for _ in range(5)
    ]
    assert plan["assembly_order"].tolist() == list(range(1, 6)) * len(
        expected_groups
    )
    for index, part_id in boundary_parts.items():
        assert plan.iloc[index]["part_id"] == part_id


def test_all_architecture_parts_exist_in_deposited_grasp_library():
    genbank = (
        Path(__file__).parents[1]
        / "grasp_library"
        / "data"
        / "profiles"
        / "grasp_nar2025"
        / "genbank"
        / "GRASP_-1.gb"
    )
    available = {record["part_id"] for record in load_grasp_records([genbank])}

    for length in (9, 14, 19):
        target = "ACGU" * (length // 4) + "ACG"[: length % 4]
        plan = compile_target_gap(target, architecture=f"{length}S")
        assert set(plan["part_id"]) <= available


@pytest.mark.parametrize(
    ("architecture", "target", "expected_message"),
    [
        ("9S", "A" * 14, "9-base"),
        ("14S", "A" * 9, "14-base"),
        ("19S", "A" * 14, "19-base"),
    ],
)
def test_architecture_target_length_mismatch_is_a_clear_value_error(
    architecture, target, expected_message
):
    with pytest.raises(ValueError, match=expected_message):
        compile_target_gap(target, architecture=architecture)


def test_architecture_and_dna_input_are_normalized():
    plan = compile_target_gap("atgcatgcatgcat", architecture=" 14s ")

    assert plan["architecture"].unique().tolist() == ["14S"]
    assert plan["target_rna"].unique().tolist() == ["AUGCAUGCAUGCAU"]


def test_unsupported_architecture_and_noncanonical_target_are_rejected():
    with pytest.raises(ValueError, match="Unsupported architecture"):
        compile_target_gap("A" * 9, architecture="10S")
    with pytest.raises(ValueError, match="noncanonical bases"):
        compile_target_gap("A" * 8 + "X", architecture="9S")
