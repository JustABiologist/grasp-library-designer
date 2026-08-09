"""Load / normalize codon-usage tables with genetic-code-consistent amino acids."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio.Seq import Seq

from .codon_validation import (
    ValidationIssue,
    analyze_cut_site_aa_risks,
    format_issues,
    reconcile_codon_table_aas,
    validate_codon_table_aas,
)
from .sample_codon_tables import SAMPLE_CODON_TABLES


def load_codon_usage(
    path,
    genetic_code: int = 1,
    *,
    force_aa_from_genetic_code: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, list]]:
    """Load codon_usage.csv and return (table, codon_data).

    Amino acids are always reconciled to `genetic_code` when
    `force_aa_from_genetic_code` is True, so switching clades cannot leave
    stale AA labels on codons.
    """
    table = pd.read_csv(path)
    table.columns = [str(column).strip().lower() for column in table.columns]

    required = {"codon", "frequency"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"Codon table needs columns {required}. Missing: {missing}"
        )

    table["codon"] = (
        table["codon"]
        .astype(str)
        .str.upper()
        .str.replace("U", "T", regex=False)
        .str.replace(r"\s+", "", regex=True)
    )
    table["frequency"] = pd.to_numeric(table["frequency"], errors="raise")

    for codon in table["codon"]:
        if len(codon) != 3 or set(codon) - set("ACGT"):
            raise ValueError(f"Invalid codon: {codon}")

    if "aa" not in table.columns or force_aa_from_genetic_code:
        # Derive / overwrite AA from the active genetic code
        table["aa"] = [
            str(Seq(codon).translate(table=genetic_code)) for codon in table["codon"]
        ]
    else:
        table["aa"] = table["aa"].astype(str).str.upper()
        table = reconcile_codon_table_aas(table, genetic_code)

    table = table[table["aa"] != "*"].copy()

    if (table["frequency"] < 0).any():
        raise ValueError("Codon frequencies must not be negative.")
    if table.groupby("aa")["frequency"].sum().eq(0).any():
        raise ValueError("At least one amino acid has all-zero frequencies.")

    table["aa_total"] = table.groupby("aa")["frequency"].transform("sum")
    table["probability"] = table["frequency"] / table["aa_total"]
    table["aa_max"] = table.groupby("aa")["frequency"].transform("max")
    table["relative_adaptiveness"] = table["frequency"] / table["aa_max"]

    # Zero-frequency sense codons get adaptiveness 0 (kept for completeness)
    table.loc[table["aa_max"] == 0, "relative_adaptiveness"] = 0.0

    codon_data: Dict[str, list] = defaultdict(list)
    for row in table.itertuples(index=False):
        codon_data[row.aa].append(
            {
                "codon": row.codon,
                "frequency": float(row.frequency),
                "probability": float(row.probability),
                "relative_adaptiveness": float(row.relative_adaptiveness),
            }
        )

    # Hard check
    issues = validate_codon_table_aas(table, genetic_code)
    errors = [i for i in issues if i.level == "error"]
    if errors:
        raise ValueError(format_issues(errors))

    return table, dict(codon_data)


def apply_organism_codon_table(
    organism_name: str,
    codon_usage_file: Path,
    *,
    genetic_code: int | None = None,
) -> Tuple[pd.DataFrame, Dict[str, list], dict, List[ValidationIssue]]:
    """Write sample table (if built-in), load with AA reconciliation, return meta + issues."""
    meta = SAMPLE_CODON_TABLES[organism_name]
    code = int(genetic_code if genetic_code is not None else meta.get("genetic_code", 1))
    issues: List[ValidationIssue] = []

    if meta["frequencies"] is None:
        if not Path(codon_usage_file).exists():
            raise FileNotFoundError(codon_usage_file)
    else:
        from .sample_codon_tables import write_sample_codon_table

        write_sample_codon_table(organism_name, codon_usage_file)

    table, codon_data = load_codon_usage(
        codon_usage_file, genetic_code=code, force_aa_from_genetic_code=True
    )
    issues.extend(validate_codon_table_aas(table, code))
    if meta.get("notes"):
        issues.append(
            ValidationIssue("info", "organism_note", str(meta["notes"]))
        )
    return table, codon_data, {**meta, "genetic_code": code}, issues


def validate_parts_for_organism(
    parts: pd.DataFrame,
    codon_data: Dict[str, list],
    *,
    genetic_code: int,
    keep_cut_sites: bool,
    minimum_relative_adaptiveness: float = 0.20,
) -> List[ValidationIssue]:
    return analyze_cut_site_aa_risks(
        parts,
        codon_data,
        genetic_code=genetic_code,
        keep_cut_sites=keep_cut_sites,
        minimum_relative_adaptiveness=minimum_relative_adaptiveness,
    )
