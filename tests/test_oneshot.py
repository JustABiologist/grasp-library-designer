import copy
from pathlib import Path

import pandas as pd
import pytest

from grasp_library.codon_tables import load_codon_usage
from grasp_library.control_panel import build_default_config
from grasp_library.dna import reverse_complement, translate_dna
from grasp_library.gga_split import assemble_payloads_to_cds
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


@pytest.mark.parametrize("target", ["UUACACGUG", "A" * 14, "G" * 19])
def test_oneshot_oligos_assemble_into_the_binder_gene(inputs, tmp_path, target):
    _, codon_data, config = inputs
    result = run_oneshot_design(
        target_rna=target,
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path / target,
        seed=42,
        log=lambda *_: None,
    )

    oligos = result["oligos"]
    design = result["design"]
    info = result["binder"]
    assert result["summary"]["architecture"] == f"{len(target)}S"
    assert result["summary"]["translation_verified"] is True
    assert len(oligos) >= 2
    assert translate_dna(result["cds"]) == info["aa_sequence"]
    assembled = assemble_payloads_to_cds(
        design.payloads,
        destination_5prime=design.destination_5prime,
        destination_3prime_coding=design.destination_3prime_coding,
    )
    assert assembled == result["cds"]
    assert oligos["order_sequence_5to3"].str.contains("GGTCTC").all()
    assert oligos["order_sequence_5to3"].str.contains("GAGACC").all()
    assert oligos["hard_constraints_passed"].all()
    assert 0.0 <= result["summary"]["ligation_fidelity"] <= 1.0
    assert result["order_csv"].exists()
    assert result["order_fasta"].exists()
    assert result["gene_fasta"].exists()


def test_oneshot_junctions_are_synonymous_and_orthogonal(inputs, tmp_path):
    _, codon_data, config = inputs
    result = run_oneshot_design(
        target_rna="UUACACGUG",
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        seed=7,
        log=lambda *_: None,
    )
    design = result["design"]
    used = set()
    for cut, overhang in zip(design.cuts, design.overhangs):
        assert design.cds[3 * cut - 4 : 3 * cut] == overhang
        assert overhang != reverse_complement(overhang)
        assert overhang not in used
        assert reverse_complement(overhang) not in used
        used.add(overhang)
        used.add(reverse_complement(overhang))
    dest = {
        design.destination_5prime,
        reverse_complement(design.destination_5prime),
        design.destination_3prime_coding,
        reverse_complement(design.destination_3prime_coding),
    }
    assert used.isdisjoint(dest)


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


def test_oneshot_honors_fragment_count_and_destination_overhangs(inputs, tmp_path):
    _, codon_data, default_config = inputs
    config = copy.deepcopy(default_config)
    config["oneshot"] = {
        "n_fragments": 4,
        "destination_5prime_overhang": "AATG",
        "destination_3prime_overhang": "GCTT",
        "wrap_enzyme": "BsaI",
    }
    result = run_oneshot_design(
        target_rna="A" * 9,
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        n_fragments=4,
        seed=3,
        log=lambda *_: None,
    )

    assert len(result["oligos"]) == 4
    assert result["summary"]["destination_5prime_overhang"] == "AATG"
    assert result["summary"]["destination_3prime_overhang"] == "GCTT"
    first = result["oligos"].iloc[0].order_sequence_5to3
    last = result["oligos"].iloc[-1].order_sequence_5to3
    assert first.startswith("TTTGGTCTCAAATG")
    assert first.count("GGTCTC") == 1
    assert last.endswith("AAGCTGAGACCAAA")
    assert result["oligos"].iloc[0].oh5_coding_site_5to3 == "AATG"
    assert result["oligos"].iloc[-1].oh3_coding_site_5to3 == reverse_complement("GCTT")


def test_oneshot_n_fragments_kwarg_is_not_rejected(inputs, tmp_path):
    _, codon_data, config = inputs
    result = run_oneshot_design(
        target_rna="A" * 9,
        codon_data=codon_data,
        config=config,
        output_dir=tmp_path,
        n_fragments=4,
        seed=1,
        log=lambda *_: None,
    )
    assert len(result["oligos"]) == 4


def test_invalid_target_characters_are_not_silently_deleted():
    from grasp_library.binder import normalize_target_rna

    with pytest.raises(ValueError, match="non-ACGU"):
        normalize_target_rna("AAAAAAAAAX")
