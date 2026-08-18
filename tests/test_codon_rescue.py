"""Cut-site depletion must outrank the codon-frequency floor.

The GRASP ``1E`` parts end inside the invariant ARELF motif, where the
``cds1_to_cds2`` coding site CTTC mask-locks the L and F codons to CTC/TTC.
On a strongly biased table (Chlamydomonas reinhardtii nuclear) the preceding
Glu codon greedily takes GAG, producing GAG|CTC|TTC -> SapI (GCTCTTC), and the
only escape (GAA) sits below the default minimum_relative_adaptiveness. The
optimizer must reach for that rare codon rather than fail, and must say so.
"""

from pathlib import Path

import pandas as pd
import pytest

from grasp_library import build_default_config, materialize_project
from grasp_library.codon_tables import load_codon_usage
from grasp_library.dna import clean_mask, contains_forbidden_site
from grasp_library.optimizer import (
    below_threshold_codons,
    build_allowed_codons,
    greedy_coding_sequence,
    optimize_coding_sequence,
    optimize_library,
    repair_forbidden_sites,
    synthesis_qc,
)
from grasp_library.sample_codon_tables import write_sample_codon_table

CHLAMY_NUCLEAR = "Chlamydomonas reinhardtii nuclear (Kazusa)"
SAPI = {"SapI": "GCTCTTC"}


@pytest.fixture(scope="module")
def chlamy_codon_data(tmp_path_factory):
    path = tmp_path_factory.mktemp("codon") / "chlamy.csv"
    write_sample_codon_table(CHLAMY_NUCLEAR, path)
    _, codon_data = load_codon_usage(path)
    return codon_data


@pytest.fixture(scope="module")
def project_config():
    project = materialize_project()
    return project, build_default_config(project / "input")


@pytest.fixture(scope="module")
def one_e_part(project_config):
    project, _ = project_config
    parts = pd.read_csv(Path(project) / "input" / "parts_full.csv")
    row = parts[parts["part_id"] == "1E_LD5N"].iloc[0]
    return str(row["aa_sequence"]), clean_mask(str(row["coding_mask"]))


def test_preferred_codons_alone_cannot_clear_sapi(chlamy_codon_data, one_e_part):
    """Guards the premise: without a rescue tier this case is unrepairable."""
    aa_sequence, coding_mask = one_e_part
    allowed = build_allowed_codons(aa_sequence, coding_mask, chlamy_codon_data, 0.20)
    seed = greedy_coding_sequence(allowed)

    assert contains_forbidden_site(seed, SAPI), "expected a SapI site in the greedy CDS"
    assert repair_forbidden_sites(seed, allowed, SAPI) is None


def test_rescue_tier_clears_sapi_and_is_reported(chlamy_codon_data, one_e_part):
    aa_sequence, coding_mask = one_e_part
    allowed = build_allowed_codons(aa_sequence, coding_mask, chlamy_codon_data, 0.20)
    rescue = build_allowed_codons(aa_sequence, coding_mask, chlamy_codon_data, 0.0)
    seed = greedy_coding_sequence(allowed)

    notes: list = []
    repaired = repair_forbidden_sites(
        seed, allowed, SAPI, rescue_codons=rescue, rescue_log=notes
    )

    assert repaired is not None
    assert not contains_forbidden_site(repaired, SAPI)
    assert len(notes) == 1
    note = notes[0]
    assert note["enzyme"] == "SapI"
    assert note["to_codon"] == "GAA"
    assert note["relative_adaptiveness"] < 0.20


def test_optimize_coding_sequence_logs_the_rescue(
    chlamy_codon_data, one_e_part, project_config
):
    _, config = project_config
    aa_sequence, coding_mask = one_e_part
    config = dict(config)
    config["forbidden_sites"] = dict(SAPI)

    notes: list = []
    sequence, _ = optimize_coding_sequence(
        aa_sequence=aa_sequence,
        coding_mask=coding_mask,
        codon_data=chlamy_codon_data,
        config=config,
        iterations=0,
        rescue_log=notes,
    )

    assert not contains_forbidden_site(sequence, SAPI)
    assert notes and notes[0]["amino_acid"] == "E"
    assert notes[0]["minimum_relative_adaptiveness"] == pytest.approx(0.20)


def test_rescue_codon_surfaces_as_a_qc_warning(
    chlamy_codon_data, one_e_part, project_config
):
    _, config = project_config
    aa_sequence, coding_mask = one_e_part
    config = dict(config)
    config["forbidden_sites"] = dict(SAPI)

    notes: list = []
    sequence, _ = optimize_coding_sequence(
        aa_sequence=aa_sequence,
        coding_mask=coding_mask,
        codon_data=chlamy_codon_data,
        config=config,
        iterations=0,
        rescue_log=notes,
    )
    qc = synthesis_qc(sequence, config, sequence_kind="cds", rescue_notes=notes)

    assert qc["status"] == "WARNING"
    assert qc["hard_constraints_passed"] is True
    assert qc["rescue_codon_count"] == 1
    assert "GAG→GAA" in qc["rescue_codon_detail"]
    assert "SapI" in qc["rescue_codon_detail"]
    assert "minimum relative adaptiveness" in qc["rescue_codon_reason"]


def test_rescue_is_scoped_to_the_single_blocking_codon(
    chlamy_codon_data, one_e_part, project_config
):
    """Exactly one codon may fall below the floor, and only the logged one."""
    _, config = project_config
    aa_sequence, coding_mask = one_e_part
    config = dict(config)
    config["forbidden_sites"] = dict(SAPI)
    threshold = config["codon_optimization"]["minimum_relative_adaptiveness"]

    notes: list = []
    sequence, _ = optimize_coding_sequence(
        aa_sequence=aa_sequence,
        coding_mask=coding_mask,
        codon_data=chlamy_codon_data,
        config=config,
        iterations=0,
        rescue_log=notes,
    )

    allowed = build_allowed_codons(
        aa_sequence, coding_mask, chlamy_codon_data, threshold
    )
    below = below_threshold_codons(sequence, allowed)

    assert below == [note["codon_index"] for note in notes]
    assert len(below) == 1
    # Every other codon in the CDS still respects the configured floor.
    assert len(below) == 1 and below[0] == 32


def test_rescue_scope_holds_across_the_whole_library(chlamy_codon_data, project_config):
    project, config = project_config
    config = dict(config)
    config["optimizer"] = {**config["optimizer"], "iterations_per_part": 0}
    parts = pd.read_csv(Path(project) / "input" / "parts_full.csv")
    threshold = config["codon_optimization"]["minimum_relative_adaptiveness"]

    library = optimize_library(parts, chlamy_codon_data, config, log=lambda _: None)

    assert len(library) == len(parts)
    for row in library.itertuples(index=False):
        allowed = build_allowed_codons(
            str(row.aa_sequence),
            clean_mask(str(row.coding_mask)),
            chlamy_codon_data,
            threshold,
        )
        below = below_threshold_codons(str(row.optimized_cds), allowed)
        # Unlogged below-threshold codons must never appear.
        assert len(below) == int(row.rescue_codon_count)
        if below:
            assert row.qc_status == "WARNING"
            assert "SapI" in row.rescue_codon_detail

    rescued = library[library["rescue_codon_count"] > 0]
    assert sorted(rescued["part_id"]) == [
        "1E_LD5N",
        "1E_LD5T",
        "1E_LN5N",
        "1E_LN5T",
    ]
    assert int(library["rescue_codon_count"].sum()) == 4


def test_unaffected_designs_use_no_rescue_codons(project_config):
    """A table where the preferred set already works must be untouched."""
    project, config = project_config
    _, codon_data = load_codon_usage(Path(project) / "input" / "codon_usage.csv")
    parts = pd.read_csv(Path(project) / "input" / "parts_full.csv")
    row = parts[parts["part_id"] == "1E_LD5N"].iloc[0]

    notes: list = []
    sequence, _ = optimize_coding_sequence(
        aa_sequence=str(row["aa_sequence"]),
        coding_mask=clean_mask(str(row["coding_mask"])),
        codon_data=codon_data,
        config=config,
        iterations=0,
        rescue_log=notes,
    )

    assert notes == []
    assert not contains_forbidden_site(sequence, config["forbidden_sites"])
    qc = synthesis_qc(sequence, config, sequence_kind="cds", rescue_notes=notes)
    assert qc["rescue_codon_count"] == 0
    assert qc["rescue_codon_detail"] == ""
