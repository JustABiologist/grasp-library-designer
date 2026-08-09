"""Validate codon tables and locked cut-site DNA against genetic code + organism usage."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq

from .dna import clean_dna, clean_mask


@dataclass
class ValidationIssue:
    level: str  # "error" | "warning" | "info"
    code: str
    message: str
    part_id: Optional[str] = None
    detail: Dict[str, Any] = field(default_factory=dict)


def translate_codon(codon: str, genetic_code: int = 1) -> str:
    codon = clean_dna(codon)
    if len(codon) != 3:
        raise ValueError(f"Not a codon: {codon!r}")
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    if codon in table.stop_codons:
        return "*"
    return table.forward_table[codon]


def validate_codon_table_aas(
    codon_table: pd.DataFrame,
    genetic_code: int = 1,
) -> List[ValidationIssue]:
    """Ensure every row's aa matches the selected NCBI genetic code."""
    issues: List[ValidationIssue] = []
    required = {"codon", "aa"}
    missing = required - set(codon_table.columns)
    if missing:
        issues.append(
            ValidationIssue(
                "error",
                "missing_columns",
                f"Codon table missing columns: {sorted(missing)}",
            )
        )
        return issues

    mismatches = []
    for row in codon_table.itertuples(index=False):
        codon = str(row.codon).upper().replace("U", "T")
        declared = str(row.aa).upper()
        try:
            expected = translate_codon(codon, genetic_code)
        except Exception as exc:
            issues.append(
                ValidationIssue(
                    "error",
                    "bad_codon",
                    f"Cannot translate {codon} with genetic code {genetic_code}: {exc}",
                )
            )
            continue
        if declared != expected:
            mismatches.append((codon, declared, expected))

    if mismatches:
        preview = ", ".join(f"{c}:{got}≠{exp}" for c, got, exp in mismatches[:8])
        more = f" (+{len(mismatches)-8} more)" if len(mismatches) > 8 else ""
        issues.append(
            ValidationIssue(
                "error",
                "aa_mismatch",
                f"{len(mismatches)} codon(s) disagree with genetic code {genetic_code}: "
                f"{preview}{more}. Re-derive aa from the genetic code before optimizing.",
                detail={"mismatches": mismatches},
            )
        )
    else:
        issues.append(
            ValidationIssue(
                "info",
                "aa_ok",
                f"All {len(codon_table)} codon→AA mappings match NCBI genetic code {genetic_code}.",
            )
        )
    return issues


def reconcile_codon_table_aas(
    codon_table: pd.DataFrame,
    genetic_code: int = 1,
) -> pd.DataFrame:
    """Return a copy with aa column forced to the selected genetic code."""
    out = codon_table.copy()
    out["codon"] = (
        out["codon"].astype(str).str.upper().str.replace("U", "T", regex=False)
    )
    out["aa"] = [translate_codon(c, genetic_code) for c in out["codon"]]
    return out


def codon_to_aa_map(
    codon_data: Mapping[str, Sequence[dict]],
) -> Dict[str, str]:
    """Map codon → AA from the active organism usage table."""
    mapping: Dict[str, str] = {}
    for aa, entries in codon_data.items():
        aa_key = str(aa).upper()
        for entry in entries:
            codon = clean_dna(entry["codon"])
            mapping[codon] = aa_key
    return mapping


def verify_cds_for_organism(
    cds: str,
    aa_sequence: str,
    *,
    genetic_code: int,
    codon_data: Mapping[str, Sequence[dict]],
) -> Dict[str, Any]:
    """
    Standard translation QC against the *selected organism codon table*.

    Passes only if:
      1. CDS length matches the protein
      2. NCBI genetic-code translation equals the expected AA
      3. Every codon is listed in ``codon_data`` for that expected AA
         (so Euglena designs are checked with Euglena codons, etc.)
    """
    from .dna import translate_dna

    cds = clean_dna(cds)
    aa = str(aa_sequence).upper().replace(" ", "").replace("*", "")
    result: Dict[str, Any] = {
        "ok": False,
        "genetic_code": int(genetic_code),
        "length_ok": len(cds) == 3 * len(aa),
        "genetic_code_ok": False,
        "codon_table_ok": False,
        "observed_protein": "",
        "bad_codons": [],
    }
    if not result["length_ok"]:
        return result

    observed = translate_dna(cds, int(genetic_code))
    result["observed_protein"] = observed
    result["genetic_code_ok"] = observed == aa

    by_codon = codon_to_aa_map(codon_data)
    bad: List[str] = []
    for i, expected_aa in enumerate(aa):
        codon = cds[3 * i : 3 * i + 3]
        table_aa = by_codon.get(codon)
        if table_aa != expected_aa:
            bad.append(f"{codon}@{i}:{table_aa or '?'}!={expected_aa}")
    result["bad_codons"] = bad[:12]
    result["codon_table_ok"] = len(bad) == 0
    result["ok"] = bool(result["genetic_code_ok"] and result["codon_table_ok"])
    return result


def cds_matches_organism(
    cds: str,
    aa_sequence: str,
    *,
    genetic_code: int,
    codon_data: Mapping[str, Sequence[dict]],
) -> bool:
    """True when CDS verifies against genetic code + organism codon table."""
    return bool(
        verify_cds_for_organism(
            cds,
            aa_sequence,
            genetic_code=genetic_code,
            codon_data=codon_data,
        )["ok"]
    )


def locked_codon_spans(coding_mask: str) -> List[tuple[int, str]]:
    """Return (aa_index, codon_mask) for positions with any fixed (non-N) base."""
    mask = clean_mask(coding_mask)
    spans = []
    for i in range(0, len(mask), 3):
        codon_mask = mask[i : i + 3]
        if len(codon_mask) == 3 and any(b != "N" for b in codon_mask):
            spans.append((i // 3, codon_mask))
    return spans


def analyze_cut_site_aa_risks(
    parts: pd.DataFrame,
    codon_data: Mapping[str, Sequence[dict]],
    *,
    genetic_code: int = 1,
    keep_cut_sites: bool = True,
    native_cds_col: str = "native_cds",
    minimum_relative_adaptiveness: float = 0.20,
) -> List[ValidationIssue]:
    """Warn when keeping fixed cut-site DNA can change or trap amino acids.

    For each part, inspect coding_mask positions with fixed bases (overhangs).
    If native CDS is available, decode the locked codon under `genetic_code` and
    compare to the declared aa_sequence. Also flag rare/absent codons in the
    new organism table.
    """
    issues: List[ValidationIssue] = []
    if not keep_cut_sites:
        issues.append(
            ValidationIssue(
                "info",
                "redesign_enabled",
                "Overhang redesign is ON — junction DNA may change (synonymously). "
                "Protein sequence must still be preserved.",
            )
        )
        # Still validate that redesign stays synonymous; native check below helps.
    else:
        issues.append(
            ValidationIssue(
                "warning",
                "keep_cut_sites",
                "Keeping GRASP cut-site / overhang DNA locks the junction codons. "
                "Under a new organism or genetic code those fixed triplets can encode "
                "different amino acids or become very rare — check warnings below.",
            )
        )

    if parts is None or len(parts) == 0:
        return issues

    has_native = native_cds_col in parts.columns
    aa_change_parts = []
    rare_parts = []
    impossible_parts = []

    adapt = {}
    for aa, entries in codon_data.items():
        for e in entries:
            adapt[(aa, e["codon"])] = float(e["relative_adaptiveness"])

    for row in parts.itertuples(index=False):
        part_id = row.part_id
        aa_seq = str(row.aa_sequence).upper().replace(" ", "")
        mask = clean_mask(row.coding_mask)
        if len(mask) != 3 * len(aa_seq):
            issues.append(
                ValidationIssue(
                    "error",
                    "mask_len",
                    f"Mask length {len(mask)} ≠ 3×{len(aa_seq)} AA",
                    part_id=part_id,
                )
            )
            continue

        native = None
        if has_native:
            native = clean_dna(getattr(row, native_cds_col))
            if len(native) != len(mask):
                native = None

        for aa_index, codon_mask in locked_codon_spans(mask):
            expected_aa = aa_seq[aa_index]
            # Reconstruct codon from native DNA if available, else only mask constraints
            if native is not None:
                codon = native[3 * aa_index : 3 * aa_index + 3]
                # Verify native matches fixed mask bases
                if not all(
                    m == "N" or m == b for m, b in zip(codon_mask, codon)
                ):
                    issues.append(
                        ValidationIssue(
                            "error",
                            "native_mask_conflict",
                            f"Native CDS codon {codon} conflicts with mask {codon_mask} "
                            f"at AA index {aa_index}",
                            part_id=part_id,
                        )
                    )
                    continue
                try:
                    decoded = translate_codon(codon, genetic_code)
                except Exception:
                    decoded = "?"
                if decoded != expected_aa:
                    aa_change_parts.append(
                        (part_id, aa_index, codon, expected_aa, decoded)
                    )
                # rarity under new organism
                rel = adapt.get((expected_aa, codon))
                if rel is None:
                    # codon does not encode expected_aa in codon_data (wrong code / missing)
                    if decoded == expected_aa:
                        impossible_parts.append(
                            (part_id, aa_index, codon, expected_aa)
                        )
                elif rel < minimum_relative_adaptiveness:
                    rare_parts.append((part_id, aa_index, codon, expected_aa, rel))
            else:
                # No native CDS: check whether any synonym matches the mask for expected AA
                options = [
                    e["codon"]
                    for e in codon_data.get(expected_aa, [])
                    if all(
                        m == "N" or m == b
                        for m, b in zip(codon_mask, e["codon"])
                    )
                ]
                if not options:
                    impossible_parts.append(
                        (part_id, aa_index, codon_mask, expected_aa)
                    )

    if aa_change_parts:
        preview = "; ".join(
            f"{p}@{i}:{codon} translates to {got} (part AA={exp}) under code {genetic_code}"
            for p, i, codon, exp, got in aa_change_parts[:6]
        )
        more = f" (+{len(aa_change_parts)-6} more)" if len(aa_change_parts) > 6 else ""
        issues.append(
            ValidationIssue(
                "error",
                "aa_would_change",
                f"Keeping cut-site DNA would CHANGE amino acids under genetic code "
                f"{genetic_code}: {preview}{more}. "
                f"Do not keep these overhangs — enable overhang redesign or update "
                f"the part AA/mask for this clade.",
                detail={"items": aa_change_parts},
            )
        )

    if impossible_parts:
        preview = ", ".join(
            f"{p}@{i}:{codon}/{aa}" for p, i, codon, aa in impossible_parts[:6]
        )
        issues.append(
            ValidationIssue(
                "error",
                "no_synonym_for_mask",
                f"Locked cut-site mask has no synonymous codon for the declared AA "
                f"in this organism table: {preview}",
                detail={"items": impossible_parts},
            )
        )

    if rare_parts and keep_cut_sites:
        preview = ", ".join(
            f"{p}@{i}:{codon}({aa},w={rel:.2f})"
            for p, i, codon, aa, rel in rare_parts[:6]
        )
        issues.append(
            ValidationIssue(
                "warning",
                "rare_locked_codon",
                f"{len(rare_parts)} locked cut-site codon(s) are rare in the selected "
                f"organism (relative adaptiveness < {minimum_relative_adaptiveness}): "
                f"{preview}. Consider enabling overhang redesign for synonym-compatible "
                f"junction DNA.",
                detail={"items": rare_parts},
            )
        )

    return issues


def format_issues(issues: Sequence[ValidationIssue]) -> str:
    lines = []
    for issue in issues:
        prefix = issue.level.upper()
        part = f" [{issue.part_id}]" if issue.part_id else ""
        lines.append(f"{prefix}{part}: {issue.message}")
    return "\n".join(lines)
