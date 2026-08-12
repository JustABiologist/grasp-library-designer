from pathlib import Path

import pandas as pd
import pytest

from grasp_library.arelf import materialize_arelf_parts
from grasp_library.assembly_interfaces import resolve_assembly_interfaces
from grasp_library.dna import reverse_complement
from grasp_library.import_grasp import (
    build_parts_table,
    compile_target_gap,
    load_grasp_records,
)
from grasp_library.optimizer import simulate_assembled_cds


def _deposited_parts() -> pd.DataFrame:
    genbank = (
        Path(__file__).parents[1]
        / "grasp_library"
        / "data"
        / "profiles"
        / "grasp_nar2025"
        / "genbank"
        / "GRASP_-1.gb"
    )
    return build_parts_table(load_grasp_records([genbank]))


def test_deposited_entry_order_arm_uses_coding_site_not_sticky_end_label():
    profile = resolve_assembly_interfaces(preset="deposited_grasp")
    entry = profile["level_minus1_entry"]

    assert entry["three_prime_end_overhang"] == "ACAA"
    assert entry["three_prime_assembled_coding_site"] == "TTGT"


def test_deposited_cds1_cds2_site_and_directional_terminals_are_distinct():
    parts = _deposited_parts().set_index("part_id")
    cds1 = parts.loc["1E_LD5N"]
    cds2 = parts.loc["2A_LD5N"]

    # The two coding windows overlap at the same plus-strand sequence.
    assert cds1.oh3 == cds1.oh3_coding_site_5to3 == "CTTC"
    assert cds2.oh5 == cds2.oh5_coding_site_5to3 == "CTTC"

    # The physical ends are written separately and directionally, both 5′->3′.
    assert cds1.three_prime_end_overhang == "GAAG"
    assert cds2.five_prime_end_overhang == "CTTC"
    assert (
        reverse_complement(cds1.three_prime_end_overhang)
        == cds2.five_prime_end_overhang
    )


def test_movable_interblock_cut_keeps_coding_site_but_orients_terminals():
    config = {
        "overhang_redesign": {"cut_mode": "movable_arelf"},
        "assembly_interfaces": {
            "junctions": {
                "terminal_to_cds2": {
                    "upstream_three_prime_end_overhang": "GAAT",
                    "downstream_five_prime_end_overhang": "ATTC",
                    "assembled_coding_site": "ATTC",
                    "arelf_offset_nt": 11,
                }
            }
        },
    }
    parts = materialize_arelf_parts(
        _deposited_parts(), {}, config=config
    ).set_index("part_id")

    assert parts.loc["1E_LD5N", "oh3_coding_site_5to3"] == "ATTC"
    assert parts.loc["2A_LD5N", "oh5_coding_site_5to3"] == "ATTC"
    assert parts.loc["1E_LD5N", "three_prime_end_overhang"] == "GAAT"
    assert parts.loc["2A_LD5N", "five_prime_end_overhang"] == "ATTC"


@pytest.mark.parametrize(
    ("target", "architecture", "expected_block_pairs"),
    [
        ("UUACACGUG", "9S", 1),
        ("A" * 14, "14S", 2),
        ("G" * 19, "19S", 3),
    ],
)
def test_assembly_checks_both_coding_overlap_and_directional_block_pair(
    target, architecture, expected_block_pairs
):
    parts = _deposited_parts()
    plan = compile_target_gap(target, architecture=architecture)
    selected = parts.set_index("part_id").loc[
        list(dict.fromkeys(plan.part_id.astype(str)))
    ]
    library = pd.DataFrame(
        {
            "part_id": selected.index,
            "optimized_part_id": [f"{part_id}_v1" for part_id in selected.index],
            "aa_sequence": selected.aa_sequence.values,
            "optimized_cds": selected.native_cds.values,
        }
    )

    assembled = simulate_assembled_cds(plan, library, parts_full=parts)
    assert assembled["translation_verified"] is True
    assert assembled["coding_junctions_verified"] is True
    assert (
        assembled["directional_terminal_pairs_checked"] == expected_block_pairs
    )

    bad = parts.copy()
    bad.loc[
        bad.part_id.str.startswith("2A_"), "five_prime_end_overhang"
    ] = "AAAA"
    with pytest.raises(ValueError, match="directional terminal mismatch"):
        simulate_assembled_cds(plan, library, parts_full=bad)
