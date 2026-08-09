"""Load / normalize codon-usage tables with genetic-code-consistent amino acids."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from Bio.Seq import Seq

from .codon_upload import frequencies_to_dataframe, save_codon_table, write_uploaded_codon_table
from .codon_validation import (
    ValidationIssue,
    analyze_cut_site_aa_risks,
    format_issues,
    reconcile_codon_table_aas,
    validate_codon_table_aas,
)
from .sample_codon_tables import (
    CUSTOM_FILE,
    FETCH_FROM_KAZUSA,
    SAMPLE_CODON_TABLES,
    UPLOAD_OWN_TABLE,
)


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


def write_frequencies_codon_table(
    frequencies: Dict[str, float],
    path,
    *,
    genetic_code: int = 1,
) -> pd.DataFrame:
    table = frequencies_to_dataframe(frequencies, genetic_code=genetic_code)
    save_codon_table(table, path)
    return table


def apply_organism_codon_table(
    organism_name: str,
    codon_usage_file: Path,
    *,
    genetic_code: int | None = None,
    kazusa_species_id: str | None = None,
    upload_bytes: Optional[Union[bytes, bytearray, memoryview, str]] = None,
    upload_filename: str | None = None,
) -> Tuple[pd.DataFrame, Dict[str, list], dict, List[ValidationIssue]]:
    """Write sample / fetched / uploaded table, load with AA reconciliation."""
    if organism_name not in SAMPLE_CODON_TABLES:
        raise KeyError(
            f"Unknown codon table {organism_name!r}. "
            f"Choose one of: {', '.join(SAMPLE_CODON_TABLES)}"
        )

    base_meta = SAMPLE_CODON_TABLES[organism_name]
    code = int(
        genetic_code if genetic_code is not None else base_meta.get("genetic_code", 1)
    )
    issues: List[ValidationIssue] = []
    meta = {**base_meta, "genetic_code": code}
    path = Path(codon_usage_file)
    path.parent.mkdir(parents=True, exist_ok=True)

    if organism_name == FETCH_FROM_KAZUSA:
        from .kazusa import fetch_kazusa_codon_table

        if not kazusa_species_id or not str(kazusa_species_id).strip():
            raise ValueError(
                "Choose Fetch from Kazusa and enter a species accession "
                "(e.g. 37762), or paste a Kazusa showcodon.cgi URL."
            )
        frequencies, fetched = fetch_kazusa_codon_table(str(kazusa_species_id))
        # Prefer page genetic code when the caller left the default / None
        if genetic_code is None and "genetic_code" in fetched:
            code = int(fetched["genetic_code"])
        write_frequencies_codon_table(frequencies, path, genetic_code=code)
        meta = {**meta, **fetched, "genetic_code": code}
        issues.append(
            ValidationIssue(
                "info",
                "kazusa_fetch",
                f"Fetched {meta.get('organism')} from Kazusa "
                f"(species {meta.get('species_id')}).",
            )
        )
    elif organism_name == UPLOAD_OWN_TABLE:
        if upload_bytes is None:
            if not path.exists():
                raise ValueError(
                    "Choose Upload your own codon table and provide a file "
                    "(CSV with codon,frequency or a Kazusa frequency block)."
                )
            # Re-use previously uploaded CSV already on disk
            meta = {
                **meta,
                "organism": "uploaded file",
                "source": f"Existing file {path}",
            }
        else:
            _, uploaded_meta = write_uploaded_codon_table(
                upload_bytes,
                path,
                genetic_code=code,
                filename=upload_filename,
            )
            meta = {**meta, **uploaded_meta, "genetic_code": code}
            issues.append(
                ValidationIssue(
                    "info",
                    "codon_upload",
                    f"Loaded uploaded codon table"
                    + (
                        f" ({upload_filename})"
                        if upload_filename
                        else ""
                    )
                    + ".",
                )
            )
    elif organism_name == CUSTOM_FILE or base_meta["frequencies"] is None:
        if not path.exists():
            raise FileNotFoundError(path)
    else:
        from .sample_codon_tables import write_sample_codon_table

        write_sample_codon_table(organism_name, path)

    table, codon_data = load_codon_usage(
        path, genetic_code=code, force_aa_from_genetic_code=True
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
