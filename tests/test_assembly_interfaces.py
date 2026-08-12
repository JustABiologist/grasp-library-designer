import copy

import pytest

from grasp_library.assembly_interfaces import (
    build_order_fragment,
    order_fragment_arms,
    resolve_assembly_interfaces,
)


def test_custom_default_uses_requested_directional_interfaces():
    profile = resolve_assembly_interfaces()

    assert profile["notation"] == "directional_terminal_5p"
    assert profile["level_minus1_entry"]["n_overhang_5p"] == "AACA"
    assert profile["level_minus1_entry"]["c_overhang_5p"] == "GGAG"
    assert profile["level0"]["ppr_outer"] == {
        "n_overhang_5p": "AGGT",
        "c_overhang_5p": "TTCG",
    }
    assert profile["final_cassette"]["n_overhang_5p"] == "GCCC"
    assert profile["final_cassette"]["c_overhang_5p"] == "GCGA"
    assert profile["level_minus1_entry"]["vector_sequence"] is None


def test_custom_order_fragment_has_two_inward_facing_bsai_sites():
    profile = resolve_assembly_interfaces()
    payload = "AGGTAAAATTCG"
    sequence = build_order_fragment(payload, profile)

    assert order_fragment_arms(profile) == (
        "TTTGGTCTCAAACA",
        "CTCCTGAGACCAAA",
    )
    assert sequence == "TTTGGTCTCAAACAAGGTAAAATTCGCTCCTGAGACCAAA"
    assert sequence.count("GGTCTC") == 1
    assert sequence.count("GAGACC") == 1


def test_deposited_grasp_preset_reproduces_entry_and_outer_interfaces():
    profile = resolve_assembly_interfaces(preset="deposited_grasp")

    assert order_fragment_arms(profile) == (
        "TTTGGTCTCAACAT",
        "TTGTTGAGACCAAA",
    )
    assert profile["level_minus1_entry"]["c_overhang_5p"] == "ACAA"
    assert profile["level0"]["acceptor_outer"] == {
        "n_overhang_5p": "CTCA",
        "c_overhang_5p": "CGAG",
    }
    assert profile["final_cassette"]["n_overhang_5p"] == "GGAG"
    assert profile["final_cassette"]["c_overhang_5p"] == "CGCT"


def test_dashboard_schema_is_normalized_and_not_ignored():
    config = {
        "assembly_interfaces": {
            "overhang_notation": "directional_terminal_5p",
            "level_minus1_entry": {
                "profile": "custom",
                "vector_name": "my entry",
                "n_terminal_overhang": "AAAA",
                "c_terminal_overhang": "CCCC",
            },
            "level0": {
                "acceptor_name": "my L0",
                "acceptor_n_terminal_overhang": "ATGC",
                "acceptor_c_terminal_overhang": "CGTA",
                "block_junctions": {
                    "cds1_to_cds2": {
                        "upstream_c": "CTTC",
                        "downstream_n": "GAAG",
                    }
                },
            },
            "level1": {
                "acceptor_name": "my L1",
                "n_terminal_overhang": "GCCC",
                "c_terminal_overhang": "GCGA",
            },
        }
    }
    profile = resolve_assembly_interfaces(config)

    assert profile["level_minus1_entry"]["vector_id"] == "my entry"
    assert profile["level_minus1_entry"]["n_overhang_5p"] == "AAAA"
    assert profile["level_minus1_entry"]["c_overhang_5p"] == "CCCC"
    assert profile["level0"]["acceptor_id"] == "my L0"
    assert profile["level0"]["acceptor_outer"]["n_overhang_5p"] == "ATGC"
    assert profile["final_cassette"]["vector_id"] == "my L1"


def test_directional_terminal_pair_is_reverse_complement_validated():
    profile = resolve_assembly_interfaces()
    bad = copy.deepcopy(profile)
    bad["junctions"]["terminal_to_cds2"]["downstream_n_5p"] = "AAAA"

    with pytest.raises(ValueError, match=r"reverse_complement\(CTTC\) != AAAA"):
        resolve_assembly_interfaces({"assembly_interfaces": bad})
