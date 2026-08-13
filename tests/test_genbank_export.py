from pathlib import Path

from Bio import SeqIO

from grasp_library.assembly_interfaces import resolve_assembly_interfaces
from grasp_library.genbank_export import fragment_features, write_annotated_genbank
from grasp_library.import_grasp import (
    build_pagm1311_order_fragment,
    build_parts_table,
    load_grasp_records,
)
from grasp_library.paths import bundled_profile_genbank


def _parts():
    return build_parts_table(
        load_grasp_records([bundled_profile_genbank() / "GRASP_-1.gb"])
    ).set_index("part_id")


def _order_row(parts, part_id: str):
    row = parts.loc[part_id]
    sequence = build_pagm1311_order_fragment(
        row.native_cds,
        part_id=part_id,
        oh5_mask_start=int(row.oh5_mask_start),
        oh3_mask_start=int(row.oh3_mask_start),
    )
    data = row.to_dict()
    data["part_id"] = part_id
    data["optimized_part_id"] = f"{part_id}_v1"
    data["optimized_cds"] = row.native_cds
    data["oligo_sequence_5to3"] = sequence
    data["genetic_code"] = 1
    return data


def _label(feature) -> str:
    return feature.qualifiers["label"][0]


def test_b_module_gb_marks_cds_bsai_and_internal_overhangs():
    parts = _parts()
    row = _order_row(parts, "B_LD5N")
    sequence, features, meta = fragment_features(row)
    labels = [_label(feature) for feature in features]
    by_label = {_label(feature): feature for feature in features}

    assert sequence.startswith("TTTGGTCTCAACAT")
    assert sequence.endswith("TTGTTGAGACCAAA")
    assert sequence[meta["insert_start"] : meta["insert_start"] + 4] == "ACTC"
    assert sequence[meta["insert_end"] - 4 : meta["insert_end"]] == "AAGA"
    assert "BsaI recognition 5'" in labels
    assert "BsaI recognition 3'" in labels
    assert "BsaI 5' cut" in labels
    assert "Level -1 overhang 5' ACAT" in labels
    assert "Level -1 overhang 3' TTGT" in labels
    assert "module overhang 5' ACTC" in labels
    assert "module overhang 3' AAGA" in labels
    assert "B_LD5N CDS" in labels
    assert by_label["B_LD5N CDS"].type == "CDS"
    assert "Level 0 overhang 5'" not in " ".join(labels)
    assert by_label["BsaI recognition 3'"].location.strand == -1


def test_a_and_e_modules_mark_level0_outer_overhangs():
    parts = _parts()
    a_seq, a_features, a_meta = fragment_features(_order_row(parts, "1A_5T_AGGT"))
    e_seq, e_features, e_meta = fragment_features(_order_row(parts, "1E_LD5N"))
    a_labels = [_label(feature) for feature in a_features]
    e_labels = [_label(feature) for feature in e_features]

    assert "Level 0 overhang 5' CTCA" in a_labels
    assert a_seq[a_meta["payload_start"] : a_meta["payload_start"] + 4] == "CTCA"
    assert a_seq[a_meta["insert_start"] : a_meta["insert_start"] + 4] == "AGGT"
    assert "Level 0 overhang 3' CGAG" in e_labels
    assert e_seq[e_meta["payload_end"] - 4 : e_meta["payload_end"]] == "CGAG"
    assert e_seq[e_meta["insert_end"] - 4 : e_meta["insert_end"]] == "CTTC"


def test_write_annotated_genbank_roundtrips(tmp_path: Path):
    parts = _parts()
    import pandas as pd

    library = pd.DataFrame(
        [
            _order_row(parts, "B_LD5N"),
            _order_row(parts, "C_LN5T"),
            _order_row(parts, "1A_5T_AGGT"),
        ]
    )
    path = write_annotated_genbank(library, tmp_path / "library.gb")
    records = list(SeqIO.parse(path, "genbank"))
    assert [record.id for record in records] == [
        "B_LD5N_v1",
        "C_LN5T_v1",
        "1A_5T_AGGT_v1",
    ]
    assert all(record.annotations.get("molecule_type") == "DNA" for record in records)
    b_labels = [feature.qualifiers["label"][0] for feature in records[0].features]
    assert "module overhang 5' ACTC" in b_labels
    assert any(feature.type == "CDS" for feature in records[0].features)
    profile = resolve_assembly_interfaces()
    assert str(records[0].seq).startswith("TTTGGTCTCAACAT")
    assert profile["order_fragment"]["enzyme"] == "BsaI"
