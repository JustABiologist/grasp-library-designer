import copy
from pathlib import Path

import pytest
from Bio import SeqIO

from grasp_library.assembly_interfaces import (
    CANONICAL_NOTATION,
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    build_order_fragment,
    order_fragment_arms,
    resolve_assembly_interfaces,
)


def test_custom_default_uses_requested_directional_interfaces():
    profile = resolve_assembly_interfaces()

    assert profile["notation"] == CANONICAL_NOTATION
    assert profile["coding_strand_direction"] == "5prime_N_to_3prime_C"
    assert profile["level_minus1_entry"][FIVE_PRIME_END] == "AACA"
    assert profile["level_minus1_entry"][THREE_PRIME_END] == "GGAG"
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "CTCC"
    assert profile["level0"]["ppr_outer"][FIVE_PRIME_END] == "AGGT"
    assert profile["level0"]["ppr_outer"][THREE_PRIME_END] == "TTCG"
    assert profile["final_cassette"][FIVE_PRIME_END] == "GCCC"
    assert profile["final_cassette"][THREE_PRIME_END] == "GCGA"
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
    assert profile["level_minus1_entry"][THREE_PRIME_END] == "ACAA"
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "TTGT"
    assert profile["level_minus1_entry"]["c_overhang_5p"] == "ACAA"
    assert profile["level0"]["acceptor_outer"][FIVE_PRIME_END] == "CTCA"
    assert profile["level0"]["acceptor_outer"][THREE_PRIME_END] == "CGAG"
    assert profile["final_cassette"][FIVE_PRIME_END] == "GGAG"
    assert profile["final_cassette"][THREE_PRIME_END] == "CGCT"


def test_deposited_entry_sites_reconstruct_annotated_clone_plus_strand():
    genbank = (
        Path(__file__).parents[1]
        / "grasp_library"
        / "data"
        / "profiles"
        / "grasp_nar2025"
        / "genbank"
        / "GRASP_-1.gb"
    )
    record = next(SeqIO.parse(genbank, "genbank"))
    inserts = [
        str(feature.extract(record.seq)).upper()
        for feature in record.features
        if feature.type == "misc_feature"
        and 40 < len(feature) < 200
        and str(feature.extract(record.seq)).upper().startswith("ACAT")
    ]
    assert len(inserts) == 1
    deposited_insert = inserts[0]
    assert deposited_insert.endswith("TTGT")

    profile = resolve_assembly_interfaces(preset="deposited_grasp")
    fragment = build_order_fragment(deposited_insert[4:-4], profile)
    order = profile["order_fragment"]
    insert_start = len(order["clamp_5p"] + order["recognition_site"] + order["spacer_5p"])
    insert_end = len(fragment) - len(
        order["spacer_3p"]
        + order["recognition_site"]
        + order["clamp_3p"]
    )
    assert fragment[insert_start:insert_end] == deposited_insert


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
    assert profile["level_minus1_entry"][FIVE_PRIME_END] == "AAAA"
    assert profile["level_minus1_entry"][THREE_PRIME_END] == "CCCC"
    assert profile["level0"]["acceptor_id"] == "my L0"
    assert profile["level0"]["acceptor_outer"][FIVE_PRIME_END] == "ATGC"
    assert profile["final_cassette"]["vector_id"] == "my L1"
    junction = profile["junctions"]["terminal_to_cds2"]
    assert junction["upstream_three_prime_end_overhang"] == "CTTC"
    assert junction["downstream_five_prime_end_overhang"] == "GAAG"
    assert junction["assembled_coding_site"] == "CTTC"


def test_directional_terminal_pair_is_reverse_complement_validated():
    profile = resolve_assembly_interfaces()
    bad = copy.deepcopy(profile)
    bad["junctions"]["terminal_to_cds2"][
        "downstream_five_prime_end_overhang"
    ] = "AAAA"

    with pytest.raises(ValueError, match=r"reverse_complement\(CTTC\) != AAAA"):
        resolve_assembly_interfaces({"assembly_interfaces": bad})


def test_legacy_entry_c_5p_is_mapped_from_opposite_strand_only():
    profile = resolve_assembly_interfaces(
        {
            "assembly_interfaces": {
                "notation": "directional_terminal_5p",
                "level_minus1_entry": {
                    "n_overhang_5p": "ACAT",
                    "c_overhang_5p": "ACAA",
                },
            }
        }
    )

    assert profile["level_minus1_entry"][FIVE_PRIME_END] == "ACAT"
    assert profile["level_minus1_entry"][THREE_PRIME_END] == "ACAA"
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "TTGT"
    assert order_fragment_arms(profile)[1].startswith("TTGT")


def test_entry_rejects_a_three_prime_site_copied_without_reverse_complementing():
    profile = resolve_assembly_interfaces()
    bad = copy.deepcopy(profile)
    bad["level_minus1_entry"][THREE_PRIME_CODING_SITE] = bad[
        "level_minus1_entry"
    ][THREE_PRIME_END]

    with pytest.raises(ValueError, match="3-prime/C-terminal-side overhang"):
        resolve_assembly_interfaces({"assembly_interfaces": bad})
