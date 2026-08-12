from pathlib import Path

import pytest

from grasp_library.colab_forms import apply_form_settings
from grasp_library.control_panel import GraspControlPanel, build_default_config
from grasp_library.assembly_interfaces import (
    CANONICAL_NOTATION,
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    resolve_assembly_interfaces,
)
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

    assert assembly["notation"] == CANONICAL_NOTATION
    assert assembly["overhang_sequence_notation"] == "5prime_to_3prime"
    assert assembly["terminal_side_convention"] == (
        "N-terminal side = 5prime; C-terminal side = 3prime"
    )
    # The old schema tag remains readable during migration.
    assert assembly["overhang_notation"] == "directional_terminal_5p"
    assert assembly["level_minus1_entry"]["n_terminal_overhang"] == "AACA"
    assert assembly["level_minus1_entry"]["c_terminal_overhang"] == "GGAG"
    assert assembly["level_minus1_entry"][FIVE_PRIME_END] == "AACA"
    assert assembly["level_minus1_entry"][THREE_PRIME_END] == "GGAG"
    assert assembly["level_minus1_entry"][FIVE_PRIME_CODING_SITE] == "AACA"
    assert assembly["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "CTCC"
    assert assembly["level0"]["block_junctions"]["cds1_to_cds2"] == {
        "upstream_c": "CTTC",
        "downstream_n": "GAAG",
        "arelf_offset_nt": 11,
    }
    assert assembly["level1"]["n_terminal_overhang"] == "GCCC"
    assert assembly["level1"]["c_terminal_overhang"] == "GCGA"
    junction = assembly["junctions"]["terminal_to_cds2"]
    assert junction["upstream_three_prime_end_overhang"] == "CTTC"
    assert junction["downstream_five_prime_end_overhang"] == "GAAG"
    assert junction["assembled_coding_site"] == "CTTC"
    assert cfg["overhang_redesign"]["cut_mode"] == "movable_arelf"
    assert cfg["overhang_redesign"]["allowed_arelf_offsets_nt"] == list(range(12))


def test_legacy_form_pairs_migrate_to_physical_terminal_ends_and_offsets(
    tmp_path: Path,
):
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
    assert assembly["junctions"]["terminal_to_cds2"] == {
        "upstream_three_prime_end_overhang": "AGTC",
        "downstream_five_prime_end_overhang": "GACT",
        "assembled_coding_site": "AGTC",
        "assembled_plus_site": "AGTC",
        "arelf_offset_nt": 7,
    }
    assert assembly["level1"]["n_terminal_overhang"] == "ACGT"
    assert assembly["level1"]["c_terminal_overhang"] == "TGCA"


def test_form_explicit_terminal_side_api_uses_physical_end_values(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    applied = apply_form_settings(
        cfg,
        **_form_kwargs(),
        entry_5prime_n_terminal_side_overhang="AAAA",
        entry_3prime_c_terminal_side_overhang="CCCC",
        cds1_to_cds2_upstream_3prime_c_terminal_side_overhang="AGTC",
        cds1_to_cds2_downstream_5prime_n_terminal_side_overhang="GACT",
        level1_5prime_n_terminal_side_overhang="ACGT",
        level1_3prime_c_terminal_side_overhang="TGCA",
    )["config"]

    assembly = applied["assembly_interfaces"]
    assert assembly["notation"] == CANONICAL_NOTATION
    assert assembly["level_minus1_entry"][FIVE_PRIME_END] == "AAAA"
    assert assembly["level_minus1_entry"][THREE_PRIME_END] == "CCCC"
    assert assembly["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "GGGG"
    junction = assembly["junctions"]["terminal_to_cds2"]
    assert junction["upstream_three_prime_end_overhang"] == "AGTC"
    assert junction["downstream_five_prime_end_overhang"] == "GACT"
    assert junction["assembled_coding_site"] == "AGTC"
    # Short historical aliases remain available to older consumers.
    assert assembly["level0"]["block_junctions"]["cds1_to_cds2"] == {
        "upstream_c": "AGTC",
        "downstream_n": "GACT",
        "arelf_offset_nt": 11,
    }


def test_dashboard_widget_labels_map_terminal_sides_to_sequence_ends(
    tmp_path: Path,
) -> None:
    panel = GraspControlPanel(build_default_config(tmp_path), input_dir=tmp_path)

    assert panel.entry_n.description == "Entry 5′ / N side"
    assert panel.entry_c.description == "Entry 3′ / C side"
    assert panel.cds1_c.description == "CDS1→2 3′ / C"
    assert panel.cds2_n.description == "CDS2←1 5′ / N"
    assert panel.entry_c.value == "GGAG"
    assert panel.cds1_c.value == "CTTC"
    assert panel.cds2_n.value == "GAAG"


def test_form_settings_reject_mismatched_physical_end_overhangs(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    with pytest.raises(ValueError, match="must be reverse complements"):
        apply_form_settings(
            cfg,
            **_form_kwargs(),
            cds1_to_cds2_upstream_3prime_c_terminal_side_overhang="CTTC",
            cds1_to_cds2_downstream_5prime_n_terminal_side_overhang="AAAA",
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
    assert profile["level_minus1_entry"][FIVE_PRIME_END] == "ACAT"
    assert profile["level_minus1_entry"][THREE_PRIME_END] == "ACAA"
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "TTGT"
    assert profile["level0"]["acceptor_id"] == "pAGM9121"
