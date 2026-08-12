from pathlib import Path

import pandas as pd
import pytest
from Bio.Data import CodonTable

from grasp_library.arelf import (
    ARELF_MOTIF,
    INTERNAL_JUNCTIONS,
    NATIVE_INTERNAL_CHOICES,
    build_arelf_candidates,
    materialize_arelf_parts,
    parse_cut_token,
    selection_overhangs,
)
from grasp_library.import_grasp import (
    build_parts_table,
    compile_target_gap,
    load_grasp_records,
)
from grasp_library.optimizer import (
    build_allowed_codons,
    greedy_coding_sequence,
    simulate_assembled_cds,
)


@pytest.fixture(scope="module")
def codon_data():
    result = {}
    for codon, aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items():
        result.setdefault(aa, []).append(
            {"codon": codon, "relative_adaptiveness": 1.0, "probability": 1.0}
        )
    for entries in result.values():
        for entry in entries:
            entry["probability"] = 1.0 / len(entries)
    return result


@pytest.fixture(scope="module")
def grasp_parts():
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


def movable_config(**junctions):
    config = {
        "overhang_redesign": {
            "cut_mode": "movable_arelf",
            "allowed_arelf_offsets_nt": list(range(12)),
        }
    }
    if junctions:
        config["assembly_interfaces"] = {
            "junctions": junctions,
        }
    return config


def test_all_synonymous_arelf_candidates_keep_offset_identity(codon_data):
    candidates = build_arelf_candidates(codon_data)
    expected_counts = [8, 8, 24, 6, 4, 8, 4, 4, 12, 6, 4, 8]

    assert set(candidates["junction"]) == set(INTERNAL_JUNCTIONS)
    for junction in INTERNAL_JUNCTIONS:
        group = candidates[candidates["junction"] == junction]
        assert len(group) == 96
        assert group["cut_token"].nunique() == 96
        assert group.groupby("motif_offset_nt").size().tolist() == expected_counts
        assert group["native"].sum() == 1


def test_cut_tokens_are_backward_compatible_but_offsets_are_not_lost():
    assert parse_cut_token("actc") == ("ACTC", None)
    assert parse_cut_token("actc@8") == ("ACTC", 8)
    assert parse_cut_token("actc@ARELF+8") == ("ACTC", 8)
    selection = {"J_ACTC": "ACTC@8", "J_AAGA": "AAGA"}
    assert selection_overhangs(selection) == {
        "J_ACTC": "ACTC",
        "J_AAGA": "AAGA",
    }


def test_materialized_internal_cuts_are_inside_arelf_and_keep_roles(
    grasp_parts, codon_data
):
    selection = {
        "J_ACTC": "ACTG@8",
        "J_AAGA": "AAGG@2",
        "J_GCAC": "GCCC@0",
        "J_TGAA": "TGAG@5",
    }
    parts = materialize_arelf_parts(
        grasp_parts,
        selection,
        config=movable_config(),
        codon_data=codon_data,
    )

    assert len(parts) == 42
    assert parts["part_id"].nunique() == 42
    assert set(parts["cut_mode"]) == {"movable_arelf"}
    for row in parts.itertuples(index=False):
        assert len(row.coding_mask) == 3 * len(row.aa_sequence)
        assert row.coding_mask[row.oh5_mask_start : row.oh5_mask_start + 4] == row.oh5
        assert row.coding_mask[row.oh3_mask_start : row.oh3_mask_start + 4] == row.oh3
        full_motifs = [
            3 * index
            for index in range(len(row.full_aa_sequence))
            if row.full_aa_sequence.startswith(ARELF_MOTIF, index)
        ]
        for prefix in ("oh5", "oh3"):
            offset = getattr(row, f"{prefix}_arelf_offset_nt")
            if pd.isna(offset):
                continue
            absolute_cut = (
                int(getattr(row, f"{prefix}_mask_start"))
                + int(row.full_window_start_nt)
            )
            assert any(
                motif_start <= absolute_cut
                and absolute_cut + 4 <= motif_start + 15
                and absolute_cut - motif_start == int(offset)
                for motif_start in full_motifs
            )

    # Deposited 5′ prefixes and recognition-code residues survive.
    first = parts.set_index("part_id").loc["1A_5T_AGGT"]
    assert first.aa_sequence.startswith("QGGNSEEPRKSFDERPERGVVSWT")
    middle = parts.set_index("part_id").loc["B_LD5N"]
    assert "PERDVVS" in middle.aa_sequence
    assert "WNAM" in middle.aa_sequence


def test_plain_overhang_selection_uses_native_arelf_offsets(grasp_parts, codon_data):
    plain = {
        junction: choice[0]
        for junction, choice in NATIVE_INTERNAL_CHOICES.items()
    }
    parts = materialize_arelf_parts(
        grasp_parts,
        plain,
        config=movable_config(),
        codon_data=codon_data,
    )
    observed = {}
    for junction in INTERNAL_JUNCTIONS:
        matches = []
        for row in parts.itertuples(index=False):
            if row.oh5_junction == junction:
                matches.append(int(row.oh5_arelf_offset_nt))
            if row.oh3_junction == junction:
                matches.append(int(row.oh3_arelf_offset_nt))
        observed[junction] = set(matches)
    assert observed == {
        junction: {offset}
        for junction, (_, offset) in NATIVE_INTERNAL_CHOICES.items()
    }


def test_internal_offset_filter_does_not_disable_fixed_block_cuts(
    grasp_parts, codon_data
):
    config = movable_config()
    config["overhang_redesign"]["allowed_arelf_offsets_nt"] = [0, 2, 5, 8]
    parts = materialize_arelf_parts(
        grasp_parts,
        {},
        config=config,
        codon_data=codon_data,
    )
    assert len(parts) == 42


def test_configured_block_interfaces_are_paired_and_materialized(
    grasp_parts, codon_data
):
    custom = {
        "terminal_to_cds2": {
            "upstream_three_prime_end_overhang": "GAAT",
            "downstream_five_prime_end_overhang": "ATTC",
            "assembled_coding_site": "ATTC",
            "arelf_offset_nt": 11,
        }
    }
    parts = materialize_arelf_parts(
        grasp_parts,
        {},
        config=movable_config(**custom),
        codon_data=codon_data,
    ).set_index("part_id")
    assert parts.loc["1E_LD5N", "oh3"] == "ATTC"
    assert parts.loc["2A_LD5N", "oh5"] == "ATTC"

    custom["terminal_to_cds2"]["downstream_five_prime_end_overhang"] = "AAAA"
    with pytest.raises(ValueError, match="reverse_complement"):
        materialize_arelf_parts(
            grasp_parts,
            {},
            config=movable_config(**custom),
            codon_data=codon_data,
        )


@pytest.mark.parametrize("length", [9, 14, 19])
def test_movable_parts_reassemble_exact_binder_without_frameshifts(
    grasp_parts, codon_data, length
):
    selection = {
        "J_ACTC": "ACTG@8",
        "J_AAGA": "AAGG@2",
        "J_GCAC": "GCCC@0",
        "J_TGAA": "TGAG@5",
    }
    parts = materialize_arelf_parts(
        grasp_parts,
        selection,
        config=movable_config(),
        codon_data=codon_data,
    )
    library_rows = []
    for row in parts.itertuples(index=False):
        allowed = build_allowed_codons(
            row.aa_sequence,
            row.coding_mask,
            codon_data,
            minimum_relative_adaptiveness=0.0,
        )
        library_rows.append(
            {
                "part_id": row.part_id,
                "optimized_part_id": f"{row.part_id}_v1",
                "aa_sequence": row.aa_sequence,
                "optimized_cds": greedy_coding_sequence(allowed),
            }
        )
    target = ("ACGU" * 5)[:length]
    plan = compile_target_gap(target, architecture=f"{length}S")
    assembled = simulate_assembled_cds(
        plan,
        pd.DataFrame(library_rows),
        parts_full=parts,
    )

    assert assembled["translation_verified"] is True
    assert assembled["stitch_warning"] == ""
    assert assembled["observed_protein"].count(ARELF_MOTIF) == length
