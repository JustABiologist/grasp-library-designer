import json
from pathlib import Path

import pytest

from grasp_library.assembly_interfaces import (
    CANONICAL_NOTATION,
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    resolve_assembly_interfaces,
)
from grasp_library.colab_forms import apply_form_settings
from grasp_library.control_panel import GraspControlPanel, build_default_config
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


def _ends(profile):
    return {
        "level_minus1": (
            profile["level_minus1_entry"][FIVE_PRIME_END],
            profile["level_minus1_entry"][THREE_PRIME_END],
        ),
        "level0": (
            profile["level0"]["acceptor_outer"][FIVE_PRIME_END],
            profile["level0"]["acceptor_outer"][THREE_PRIME_END],
        ),
        "level1": (
            profile["final_cassette"][FIVE_PRIME_END],
            profile["final_cassette"][THREE_PRIME_END],
        ),
    }


def test_dashboard_defaults_are_the_deposited_grasp_toolbox_overhangs(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    profile = resolve_assembly_interfaces(cfg)

    assert profile["profile_name"] == "deposited_grasp"
    assert profile["notation"] == CANONICAL_NOTATION
    assert _ends(profile) == {
        "level_minus1": ("ACAT", "ACAA"),
        "level0": ("CTCA", "CTCG"),
        "level1": ("GGAG", "AGCG"),
    }
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "TTGT"
    assert profile["level0"]["acceptor_outer"][THREE_PRIME_CODING_SITE] == "CGAG"
    assert profile["final_cassette"][THREE_PRIME_CODING_SITE] == "CGCT"


def test_form_exposes_only_one_5prime_3prime_pair_per_level(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    applied = apply_form_settings(
        cfg,
        **_form_kwargs(),
        level_minus1_5prime_overhang="AAAA",
        level_minus1_3prime_overhang="CCCC",
        level0_5prime_overhang="ATGC",
        level0_3prime_overhang="CGTA",
        level1_5prime_overhang="ACGT",
        level1_3prime_overhang="TGCA",
    )["config"]
    profile = resolve_assembly_interfaces(applied)

    assert profile["profile_name"] == "custom"
    assert _ends(profile) == {
        "level_minus1": ("AAAA", "CCCC"),
        "level0": ("ATGC", "CGTA"),
        "level1": ("ACGT", "TGCA"),
    }
    assert profile["level_minus1_entry"][FIVE_PRIME_CODING_SITE] == "AAAA"
    assert profile["level_minus1_entry"][THREE_PRIME_CODING_SITE] == "GGGG"
    assert profile["level0"]["acceptor_outer"][THREE_PRIME_CODING_SITE] == "TACG"
    assert profile["final_cassette"][THREE_PRIME_CODING_SITE] == "TGCA"


def test_dashboard_overhang_section_has_exactly_six_plain_end_fields(tmp_path: Path):
    panel = GraspControlPanel(build_default_config(tmp_path), input_dir=tmp_path)
    fields = (
        panel.level_minus1_5prime,
        panel.level_minus1_3prime,
        panel.level0_5prime,
        panel.level0_3prime,
        panel.level1_5prime,
        panel.level1_3prime,
    )

    assert [field.description for field in fields] == [
        "5′ overhang",
        "3′ overhang",
        "5′ overhang",
        "3′ overhang",
        "5′ overhang",
        "3′ overhang",
    ]
    assert [field.value for field in fields] == [
        "ACAT",
        "ACAA",
        "CTCA",
        "CTCG",
        "GGAG",
        "AGCG",
    ]


def test_notebook_forms_expose_the_same_six_overhang_fields():
    root = Path(__file__).parents[1]
    notebooks = (
        root / "grasp_library_designer.ipynb",
        root / "grasp_oneshot_designer.ipynb",
        root / "grasp_library" / "notebooks" / "grasp_library_designer.ipynb",
        root / "grasp_library" / "notebooks" / "grasp_oneshot_designer.ipynb",
    )
    expected = {
        "level_minus1_5prime_overhang": "ACAT",
        "level_minus1_3prime_overhang": "ACAA",
        "level0_5prime_overhang": "CTCA",
        "level0_3prime_overhang": "CTCG",
        "level1_5prime_overhang": "GGAG",
        "level1_3prime_overhang": "AGCG",
    }

    for path in notebooks:
        notebook = json.loads(path.read_text())
        source = "".join(
            "".join(cell.get("source", [])) for cell in notebook["cells"]
        )
        observed = {}
        for line in source.splitlines():
            if "#@param" not in line or "_overhang =" not in line:
                continue
            name, value = line.split("=", 1)
            observed[name.strip()] = value.split("#", 1)[0].strip().strip('"')
        assert observed == expected, path
def test_form_rejects_invalid_level_overhang(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    with pytest.raises(ValueError, match="Level 0 3′ overhang"):
        apply_form_settings(
            cfg,
            **_form_kwargs(),
            level0_3prime_overhang="XYZ",
        )


def test_deposited_preset_rejects_edited_standard_overhangs(tmp_path: Path):
    cfg = build_default_config(tmp_path)
    with pytest.raises(ValueError, match="requires the deposited 5′/3′ overhangs"):
        apply_form_settings(
            cfg,
            **_form_kwargs(),
            assembly_interface_preset="deposited_grasp",
            level1_3prime_overhang="AAAA",
        )
