from __future__ import annotations

from unittest.mock import patch

import pytest

import pandas as pd

from grasp_library.ligation_fidelity import LigationFidelityCalculator
from grasp_library.optimizer import synthesis_qc
from grasp_library.pareto import score_overhang_set_ggassembler
from grasp_library.synthesis_vendors import (
    apply_enzyme_to_config,
    apply_vendor_to_config,
)
from grasp_library.workflows import qc_report_frame, write_qc_sheet


@pytest.fixture
def qc_config() -> dict:
    return {
        "synthesis": {
            "global_gc_min": 0.0,
            "global_gc_max": 1.0,
            "window_size": 10,
            "window_gc_min": 0.0,
            "window_gc_max": 1.0,
            "max_homopolymer": 3,
            "repeat_k": 8,
            "max_repeat_count": 999,
            "min_oligo_length": 20,
            "max_oligo_length": 300,
            "min_gene_length": 100,
            "max_gene_length": 5_000,
        },
        "forbidden_sites": {},
        "synthesis_vendor_meta": {
            "hard_rules": ["A manual vendor check"],
            "manual_review_rules": ["A manual vendor check"],
            "machine_hard_constraints": {"max_homopolymer": 13},
        },
    }


def test_homopolymer_qc_lists_every_overlong_stretch(qc_config: dict) -> None:
    result = synthesis_qc("ATGAAAAAAATGCTTTTTGC", qc_config, sequence_kind="cds")

    assert "A7@4" in result["homopolymer_detail"]
    assert "T5@14" in result["homopolymer_detail"]
    assert result["homopolymer_count"] == 2
    assert "Homopolymer too long: A7@4; T5@14" == result["warnings"]


def test_high_and_low_gc_are_labelled_with_percent(qc_config: dict) -> None:
    qc_config["synthesis"]["global_gc_min"] = 0.25
    qc_config["synthesis"]["global_gc_max"] = 0.65
    qc_config["synthesis"]["window_gc_min"] = 0.20
    qc_config["synthesis"]["window_gc_max"] = 0.80
    qc_config["synthesis"]["max_homopolymer"] = 20

    high = synthesis_qc("GCGCGCGCGCGC", qc_config, sequence_kind="cds")
    assert high["gc_status"].startswith("HIGH")
    assert "100%" in high["gc_status"]
    assert "Global GC HIGH" in high["warnings"]
    assert "HIGH" in high["local_gc_detail"]

    low = synthesis_qc("ATATATATATAT", qc_config, sequence_kind="cds")
    assert low["gc_status"].startswith("LOW")
    assert "Global GC LOW" in low["warnings"]


def test_blacklist_cutters_are_reported(qc_config: dict) -> None:
    qc_config["forbidden_sites"] = {"SapI": "GCTCTTC", "BsaI": "GGTCTC"}
    result = synthesis_qc("AAAGCTCTTCAAAGGTCTCAAA", qc_config, sequence_kind="cds")

    assert result["blacklist_tested"] == "SapI, BsaI"
    assert "SapI GCTCTTC @4" in result["blacklist_hits"]
    assert "BsaI GGTCTC @14" in result["blacklist_hits"]
    assert "Forbidden restriction site in CDS: SapI" in result["failures"]


def test_excel_qc_sheet_includes_detail_columns(tmp_path) -> None:
    library = pd.DataFrame(
        [
            {
                "optimized_part_id": "B_LD5N_v1",
                "part_id": "B_LD5N",
                "qc_status": "WARNING",
                "qc_passed": False,
                "translation_verified": True,
                "mask_verified": True,
                "codon_score": -0.2,
                "cds_gc_pct": 72.0,
                "cds_gc_status": "HIGH 72% (max 65%)",
                "cds_local_gc": "HIGH 88% @ 13–62",
                "cds_longest_homopolymer": 6,
                "cds_homopolymers": "A6@13; T5@40",
                "cds_repeats": "",
                "blacklist_tested": "SapI, BsaI, BpiI",
                "blacklist_hits": "none",
                "cds_warnings": "Global GC HIGH 72% (max 65%)",
                "cds_failures": "",
                "oligo_length": 180,
                "oligo_gc_pct": 55.0,
                "oligo_gc_status": "OK 55%",
                "oligo_homopolymers": "A6@27",
                "oligo_warnings": "Homopolymer too long: A6@27",
                "oligo_failures": "",
            }
        ]
    )
    xlsx = tmp_path / "qc.xlsx"
    with pd.ExcelWriter(xlsx) as writer:
        write_qc_sheet(writer, library)
    qc = pd.read_excel(xlsx, sheet_name="QC")

    assert list(qc.columns) == list(qc_report_frame(library).columns)
    assert qc.loc[0, "cds_homopolymers"] == "A6@13; T5@40"
    assert qc.loc[0, "cds_gc_status"] == "HIGH 72% (max 65%)"
    assert qc.loc[0, "blacklist_tested"] == "SapI, BsaI, BpiI"


def test_soft_synthesis_deviation_is_warning_not_clean_pass(qc_config: dict) -> None:
    result = synthesis_qc("ATGCAAAAATGC", qc_config, sequence_kind="cds")

    assert result["status"] == "WARNING"
    assert result["hard_constraints_passed"] is True
    assert result["passed"] is False
    assert "Homopolymer too long" in result["warnings"]
    assert result["failures"] == ""
    assert result["vendor_acceptance_confirmed"] is False
    assert result["manual_vendor_rules"] == "A manual vendor check"


def test_machine_checkable_vendor_hard_rule_fails(qc_config: dict) -> None:
    result = synthesis_qc("A" * 14, qc_config, sequence_kind="cds")

    assert result["status"] == "FAIL"
    assert result["hard_constraints_passed"] is False
    assert result["passed"] is False
    assert "vendor hard maximum (13 nt)" in result["failures"]


def test_orderable_oligo_length_bounds_are_hard(qc_config: dict) -> None:
    result = synthesis_qc("ATGCATGC", qc_config, sequence_kind="oligo")

    assert result["status"] == "FAIL"
    assert "Oligo shorter than configured minimum (20 bp)" in result["failures"]


def test_generic_cds_does_not_inherit_oligo_or_gene_length_bounds(
    qc_config: dict,
) -> None:
    result = synthesis_qc("ATGCATGC", qc_config, sequence_kind="cds")

    assert result["status"] == "PASS"
    assert result["passed"] is True


def test_unknown_sequence_kind_is_rejected(qc_config: dict) -> None:
    with pytest.raises(ValueError, match="Unknown sequence_kind"):
        synthesis_qc("ATGC", qc_config, sequence_kind="plasmid")


def test_twist_profile_exposes_machine_and_manual_hard_rules(
    qc_config: dict,
) -> None:
    updated = apply_vendor_to_config(
        qc_config, "Twist · Standard gene guidelines"
    )

    meta = updated["synthesis_vendor_meta"]
    assert meta["machine_hard_constraints"]["max_homopolymer"] == 13
    assert meta["manual_review_rules"] == ["Do not include CcdB toxin sequences"]


def test_enzyme_selector_is_explicitly_a_domestication_filter(
    qc_config: dict,
) -> None:
    updated = apply_enzyme_to_config(qc_config, "BpiI / BbsI (GAAGAC)")

    assert updated["forbidden_sites"] == {"BpiI": "GAAGAC"}
    assert updated["domestication_enzyme_filter"] == "BpiI / BbsI (GAAGAC)"
    assert updated["assembly_enzyme_semantics"] == "internal_site_filter_only"
    assert updated["assembly_flanks_modified"] is False


def test_grasp_9s_fidelity_reports_one_physical_level0_reaction() -> None:
    calc = LigationFidelityCalculator()
    junctions = {
        "J_ACTC": "ACTC",
        "J_AAGA": "AAGA",
        "J_GCAC": "GCAC",
        "J_TGAA": "TGAA",
    }

    with patch.object(calc, "set_fidelity", return_value=0.8) as scorer:
        combined = calc.grasp_9s_first_stage_fidelity(junctions)

    assert combined == pytest.approx(0.8)
    scorer.assert_called_once_with(["CTCA", "ACTC", "AAGA", "GCAC", "TGAA", "CGAG"])


def test_grasp_9s_fidelity_requires_named_junctions() -> None:
    calc = LigationFidelityCalculator()

    with pytest.raises(ValueError, match="J_TGAA"):
        calc.grasp_9s_first_stage_fidelity({"J_ACTC": "ACTC"})


def test_grasp_14s_fidelity_is_not_a_product_of_separate_tubes() -> None:
    calc = LigationFidelityCalculator()
    junctions = {
        "J_ACTC": "ACTC",
        "J_AAGA": "AAGA",
        "J_GCAC": "GCAC",
        "J_TGAA": "TGAA",
    }

    with patch.object(calc, "set_fidelity", return_value=0.9) as scorer:
        combined = calc.grasp_first_stage_fidelity(junctions, architecture="14S")

    assert combined == pytest.approx(0.9)
    scorer.assert_called_once()

    with patch.object(calc, "set_fidelity", return_value=0.9):
        estimate = calc.grasp_architecture_success_product_estimate(
            junctions, architecture="14S"
        )
    assert estimate == pytest.approx(0.9**3)


def test_redesign_rejects_collision_with_fixed_pagm9121_overhang() -> None:
    calc = LigationFidelityCalculator()
    selection = {
        "J_ACTC": "CTCA",  # collides with fixed pAGM9121 left overhang
        "J_AAGA": "AAGA",
        "J_GCAC": "GCAC",
        "J_TGAA": "TGAA",
    }

    assert score_overhang_set_ggassembler(selection, calc) == 0.0
