import copy
from pathlib import Path

import pandas as pd
import pytest

from grasp_library.codon_tables import load_codon_usage
from grasp_library.control_panel import build_default_config
from grasp_library.import_grasp import build_pagm1311_order_fragment, import_grasp_profile
from grasp_library.oneshot import run_oneshot_design, validate_pagm1311_order_fragment
from grasp_library.paths import bundled_profile_genbank


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp("oneshot-input")
    imported = import_grasp_profile(bundled_profile_genbank(), root)
    source = Path(__file__).parents[1] / "grasp_library_project" / "input" / "codon_usage.csv"
    codon_file = root / "codon_usage.csv"
    codon_file.write_bytes(source.read_bytes())
    _, codon_data = load_codon_usage(codon_file)
    config = build_default_config(root)
    config["optimizer"]["iterations_per_part"] = 0
    return imported, codon_data, config


def test_deposited_part_reconstructs_orderable_pagm1311_fragment(inputs):
    imported, _, _ = inputs
    row = imported["parts_full"].set_index("part_id").loc["B_LD5N"]
    order_sequence = build_pagm1311_order_fragment(
        row.native_cds,
        part_id="B_LD5N",
        oh5_mask_start=row.oh5_mask_start,
        oh3_mask_start=row.oh3_mask_start,
    )
    result = validate_pagm1311_order_fragment(order_sequence)

    assert order_sequence.startswith("TTTGGTCTCAACAT")
    assert order_sequence.endswith("TTGTTGAGACCAAA")
    assert result["cloned_pagm1311_context_5to3"].startswith("GAAGAC")
    assert result["cloned_pagm1311_context_5to3"].endswith("GTCTTC")
    assert result["bpii_release_oh5"] == "ACTC"
    assert result["bpii_release_oh3"] == "AAGA"


def test_all_42_deposited_parts_reconstruct_valid_entry_fragments(inputs):
    imported, _, _ = inputs
    parts = imported["parts_full"]

    for row in parts.itertuples(index=False):
        order_sequence = build_pagm1311_order_fragment(
            row.native_cds,
            part_id=row.part_id,
            oh5_mask_start=row.oh5_mask_start,
            oh3_mask_start=row.oh3_mask_start,
        )
        result = validate_pagm1311_order_fragment(order_sequence)
        assert order_sequence.count("GGTCTC") == 1
        assert order_sequence.count("GAGACC") == 1
        assert result["entry_insertion_overhang_5"] == "ACAT"
        assert result["entry_insertion_overhang_3"] == "TTGT"
        assert result["entry_vector_context_in_silico_validated"] is True
        assert result["entry_vector_sequence_in_silico_validated"] is False
        expected_oh5 = "CTCA" if row.part_id.split("_", 1)[0].endswith("A") else row.oh5
        expected_oh3 = row.oh3
        if row.part_id.split("_", 1)[0].endswith("E"):
            expected_oh3 = "CGAG"
        assert result["bpii_release_oh5"] == expected_oh5
        assert result["bpii_release_oh3"] == expected_oh3


@pytest.mark.parametrize(
    ("target", "architecture", "interfaces"),
    [
        ("UUACACGUG", "9S", "AGGT;CTTC;TTCG"),
        ("A" * 14, "14S", "AGGT;GTGA;CTTC;TTCG"),
        ("G" * 19, "19S", "AGGT;GTGA;CACG;CTTC;TTCG"),
    ],
)
def test_oneshot_exports_custom_vector_requirements_and_ppr_chain(
    inputs, tmp_path, target, architecture, interfaces
):
    _, codon_data, config = inputs
    result = run_oneshot_design(
        target_rna=target,
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        seed=42,
        log=lambda *_: None,
    )

    assert result["summary"]["architecture"] == architecture
    assert result["summary"]["entry_vector"] == "Custom Level -1 entry vector"
    assert result["summary"]["level0_acceptor"] == "Custom Level 0 acceptor"
    assert result["summary"]["entry_five_prime_end_overhang"] == "AACA"
    assert result["summary"]["entry_three_prime_end_overhang"] == "GGAG"
    assert result["summary"]["entry_three_prime_assembled_coding_site"] == "CTCC"
    assert result["summary"][
        "final_cassette_five_prime_end_overhang"
    ] == "GCCC"
    assert result["summary"][
        "final_cassette_three_prime_end_overhang"
    ] == "GCGA"
    assert result["summary"]["translation_verified"] is True
    assert len(result["level0_assemblies"]) == len(interfaces.split(";")) - 1
    assert result["level0_assemblies"][
        "level0_module_chain_in_silico_validated"
    ].all()
    assert set(result["level0_assemblies"]["level0_outer_left_5to3"]) == {"CTCA"}
    assert set(result["level0_assemblies"]["level0_outer_right_5to3"]) == {"CGAG"}
    assert result["ppr_block_chain"]["ppr_interface_chain"] == interfaces
    assert result["ppr_block_chain"][
        "ppr_block_chain_in_silico_validated"
    ] is True
    assert result["summary"]["entry_interface_requirements_checked"] is True
    assert result["summary"]["entry_vector_sequence_in_silico_validated"] is False
    assert result["summary"]["standalone_expression_cassette"] is False
    assert result["orderable_fragments"]["order_sequence_5to3"].str.startswith(
        "TTTGGTCTCAAACA"
    ).all()
    assert result["orderable_fragments"]["order_sequence_5to3"].str.endswith(
        "CTCCTGAGACCAAA"
    ).all()
    assert result["orderable_fragments"]["order_sequence_5to3"].str.count(
        "GGTCTC"
    ).eq(1).all()
    assert result["orderable_fragments"]["order_sequence_5to3"].str.count(
        "GAGACC"
    ).eq(1).all()
    assert result["order_csv"].exists()
    assert result["order_fasta"].exists()
    assert result["ppr_block_chain_csv"].exists()
    assert result["binding_tract_fasta"].exists()


@pytest.mark.parametrize("target", ["A" * 9, "G" * 9, "C" * 9, "U" * 9, "ACGUACGUA"])
def test_difficult_targets_no_longer_crash_in_cut_search(inputs, tmp_path, target):
    _, codon_data, config = inputs
    result = run_oneshot_design(
        target_rna=target,
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path / target,
        seed=7,
        log=lambda *_: None,
    )
    assert result["summary"]["translation_verified"] is True


def test_one_shot_rejects_obsolete_arbitrary_fragmentation(inputs, tmp_path):
    _, codon_data, config = inputs
    with pytest.raises(ValueError, match="fixed five-part"):
        run_oneshot_design(
            target_rna="A" * 9,
            codon_data=codon_data,
            config=config,
            output_dir=tmp_path,
            n_fragments=4,
        )


def test_oneshot_deposited_grasp_preset_is_explicit(inputs, tmp_path):
    _, codon_data, default_config = inputs
    config = copy.deepcopy(default_config)
    config["assembly_interfaces"] = "deposited_grasp"
    result = run_oneshot_design(
        target_rna="UUACACGUG",
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        seed=42,
        log=lambda *_: None,
    )

    assert result["summary"]["entry_vector"] == "pAGM1311"
    assert result["summary"]["level0_acceptor"] == "pAGM9121"
    assert result["summary"]["entry_vector_context_in_silico_validated"] is True
    assert result["summary"]["entry_vector_sequence_in_silico_validated"] is False
    assert result["orderable_fragments"]["order_sequence_5to3"].str.startswith(
        "TTTGGTCTCAACAT"
    ).all()
    assert result["orderable_fragments"]["order_sequence_5to3"].str.endswith(
        "TTGTTGAGACCAAA"
    ).all()


def test_oneshot_honors_configured_level0_block_interface(inputs, tmp_path):
    _, codon_data, default_config = inputs
    config = copy.deepcopy(default_config)
    junction = config["assembly_interfaces"]["level0"]["block_junctions"][
        "cds1_to_cds2"
    ]
    junction.update(
        upstream_c="ATTC",
        downstream_n="GAAT",
        arelf_offset_nt=11,
    )
    result = run_oneshot_design(
        target_rna="UUACACGUG",
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        seed=42,
        log=lambda *_: None,
    )

    assert result["ppr_block_chain"]["ppr_interface_chain"] == "AGGT;ATTC;TTCG"
    assert result["summary"]["ppr_block_chain_in_silico_validated"] is True


def test_invalid_target_characters_are_not_silently_deleted():
    from grasp_library.binder import normalize_target_rna

    with pytest.raises(ValueError, match="non-ACGU"):
        normalize_target_rna("AAAAAAAAAX")
