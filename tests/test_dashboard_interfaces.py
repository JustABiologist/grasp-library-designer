from pathlib import Path

import pytest

from grasp_library.colab_forms import apply_form_settings
from grasp_library.control_panel import build_default_config
from grasp_library.assembly_interfaces import resolve_assembly_interfaces
from grasp_library.sample_codon_tables import sample_names
from grasp_library.synthesis_vendors import (
    enzyme_names,
    ligation_table_names,
    vendor_names,
)


def _form_kwargs():
    return {
        "organism": sample_names()[0],
        "genetic_code": 1,
        "target_rna": "UUACACGUG",
        "synthesis_vendor": vendor_names()[0],
        "assembly_enzyme": enzyme_names()[0],
        "ligation_table": ligation_table_names()[0],
        "prompt_upload_if_needed": False,
    }


def test_dashboard_defaults_expose_all_cloning_layers(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    assembly = cfg["assembly_interfaces"]

    assert assembly["overhang_notation"] == "directional_terminal_5p"
    assert assembly["level_minus1_entry"]["n_terminal_overhang"] == "AACA"
    assert assembly["level_minus1_entry"]["c_terminal_overhang"] == "GGAG"
    assert assembly["level0"]["block_junctions"]["cds1_to_cds2"] == {
        "upstream_c": "CTTC",
        "downstream_n": "GAAG",
        "arelf_offset_nt": 11,
    }
    assert assembly["level1"]["n_terminal_overhang"] == "GCCC"
    assert assembly["level1"]["c_terminal_overhang"] == "GCGA"
    assert cfg["overhang_redesign"]["cut_mode"] == "movable_arelf"
    assert cfg["overhang_redesign"]["allowed_arelf_offsets_nt"] == list(range(12))


def test_form_settings_apply_directional_pairs_and_arelf_offsets(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    applied = apply_form_settings(
        cfg,
        **_form_kwargs(),
        entry_n_overhang="AAAA",
        entry_c_overhang="CCCC",
        cds1_c_overhang="AGTC",
        cds2_n_overhang="GACT",
        cds1_cds2_arelf_offset_nt=7,
        level1_n_overhang="ACGT",
        level1_c_overhang="TGCA",
    )["config"]

    assembly = applied["assembly_interfaces"]
    assert assembly["level_minus1_entry"]["n_terminal_overhang"] == "AAAA"
    assert assembly["level_minus1_entry"]["c_terminal_overhang"] == "CCCC"
    assert assembly["level0"]["block_junctions"]["cds1_to_cds2"] == {
        "upstream_c": "AGTC",
        "downstream_n": "GACT",
        "arelf_offset_nt": 7,
    }
    assert assembly["level1"]["n_terminal_overhang"] == "ACGT"
    assert assembly["level1"]["c_terminal_overhang"] == "TGCA"


def test_form_settings_reject_mismatched_directional_junctions(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    with pytest.raises(ValueError, match="reverse complements"):
        apply_form_settings(
            cfg,
            **_form_kwargs(),
            cds1_c_overhang="CTTC",
            cds2_n_overhang="AAAA",
        )


def test_form_settings_reject_cut_outside_arelf(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    with pytest.raises(ValueError, match="between 0 and 11"):
        apply_form_settings(
            cfg,
            **_form_kwargs(),
            cds1_cds2_arelf_offset_nt=12,
        )


def test_form_can_select_deposited_grasp_vector_preset(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    applied = apply_form_settings(
        cfg,
        **_form_kwargs(),
        assembly_interface_preset="deposited_grasp",
    )["config"]
    profile = resolve_assembly_interfaces(applied)

    assert profile["level_minus1_entry"]["vector_id"] == "pAGM1311"
    assert profile["level_minus1_entry"]["n_overhang_5p"] == "ACAT"
    assert profile["level_minus1_entry"]["c_overhang_5p"] == "ACAA"
    assert profile["level0"]["acceptor_id"] == "pAGM9121"
