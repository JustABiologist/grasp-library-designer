from pathlib import Path

import pytest

from grasp_library.colab_forms import apply_form_settings
from grasp_library.control_panel import GraspControlPanel, build_default_config
from grasp_library.dna import contains_forbidden_site
from grasp_library.restriction_sites import (
    COMMON_RESTRICTION_SITES,
    DEFAULT_SITE_BLACKLIST,
    apply_site_blacklist_to_config,
    parse_site_blacklist,
    resolve_site_blacklist,
)
from grasp_library.sample_codon_tables import sample_names
from grasp_library.synthesis_vendors import (
    apply_enzyme_to_config,
    enzyme_names,
    redesign_ligation_table_names,
    vendor_names,
)


def test_catalog_has_one_hundred_common_enzymes() -> None:
    assert len(COMMON_RESTRICTION_SITES) == 100
    assert DEFAULT_SITE_BLACKLIST == ("SapI", "BsaI", "BpiI")
    for name in DEFAULT_SITE_BLACKLIST:
        assert name in COMMON_RESTRICTION_SITES


def test_default_config_forbids_sapi_bsai_and_bpii(tmp_path: Path) -> None:
    cfg = build_default_config(tmp_path)

    assert cfg["site_blacklist"] == ["SapI", "BsaI", "BpiI"]
    assert cfg["forbidden_sites"]["SapI"] == "GCTCTTC"
    assert cfg["forbidden_sites"]["BsaI"] == "GGTCTC"
    assert cfg["forbidden_sites"]["BpiI"] == "GAAGAC"


def test_parse_site_blacklist_accepts_labels_and_aliases() -> None:
    assert parse_site_blacklist("SapI, BsaI, BpiI") == ["SapI", "BsaI", "BpiI"]
    assert parse_site_blacklist("sapi / bspqi (GCTCTTC)") == ["SapI"]
    with pytest.raises(ValueError, match="Unknown restriction enzyme"):
        parse_site_blacklist("NotARealEnzyme")


def test_contains_forbidden_site_finds_sapi_and_reverse_complement() -> None:
    hits = contains_forbidden_site("AAAGCTCTTCAAA", {"SapI": "GCTCTTC"})
    assert hits[0]["enzyme"] == "SapI"
    assert hits[0]["start_0based"] == 3

    rc_hits = contains_forbidden_site("AAAGAAGAGCAAA", {"SapI": "GCTCTTC"})
    assert rc_hits[0]["site"] == "GAAGAGC"


def test_contains_forbidden_site_matches_iupac_motifs() -> None:
    hits = contains_forbidden_site("AAAGTCGACAAA", {"HincII": "GTYRAC"})
    assert hits[0]["enzyme"] == "HincII"
    assert hits[0]["site"] == "GTCGAC"


def test_apply_site_blacklist_merges_with_assembly_filter() -> None:
    cfg = apply_enzyme_to_config({}, "None (no enzyme filter)")
    updated = apply_site_blacklist_to_config(cfg, ["SapI", "EcoRI"])

    assert updated["site_blacklist"] == ["SapI", "EcoRI"]
    assert updated["forbidden_sites"] == {
        "SapI": "GCTCTTC",
        "EcoRI": "GAATTC",
    }

    cleared = apply_site_blacklist_to_config(
        apply_enzyme_to_config(updated, "GRASP default · BsaI + BpiI + BsmBI"),
        [],
    )
    assert cleared["site_blacklist"] == []
    assert "SapI" not in cleared["forbidden_sites"]
    assert set(cleared["forbidden_sites"]) == {"BsaI", "BpiI", "BsmBI"}


def test_dashboard_blacklist_defaults_to_sapi_bsai_bpii(tmp_path: Path) -> None:
    panel = GraspControlPanel(build_default_config(tmp_path), input_dir=tmp_path)

    assert tuple(panel.site_blacklist.value) == ("SapI", "BsaI", "BpiI")
    assert "EcoRI (GAATTC)" in [label for label, _ in panel.site_blacklist.options]


def test_form_settings_apply_default_blacklist(tmp_path: Path) -> None:
    applied = apply_form_settings(
        build_default_config(tmp_path),
        organism=sample_names()[0],
        genetic_code=1,
        target_rna="UUACACGUG",
        synthesis_vendor=vendor_names()[0],
        assembly_enzyme=enzyme_names()[0],
        ligation_table=redesign_ligation_table_names()[0],
        prompt_upload_if_needed=False,
    )["config"]

    assert applied["site_blacklist"] == ["SapI", "BsaI", "BpiI"]
    assert "SapI" in applied["forbidden_sites"]
    assert resolve_site_blacklist(DEFAULT_SITE_BLACKLIST)["SapI"] == "GCTCTTC"
