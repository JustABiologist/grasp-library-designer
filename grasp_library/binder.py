"""RNA → GRASP binder protein (PPR recognition code), no library parts required."""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

# Classic PPR code: (5th AA, last AA) of each repeat
_RNA_TO_CODE = {
    "A": ("T", "N"),
    "C": ("N", "N"),
    "G": ("T", "D"),
    "U": ("N", "D"),
    "T": ("N", "D"),
}
_CODE_TO_RNA = {pair: base for base, pair in _RNA_TO_CODE.items() if base != "T"}

_JOIN = {
    "A": ("", "B"),
    "B": ("A", "C"),
    "C": ("B", "D"),
    "D": ("C", "E"),
    "E": ("D", ""),
}

MODULE_ROLE_COLUMNS = (
    "assembly_role",
    "joins_upstream_role",
    "joins_downstream_role",
    "ppr_5th_aa",
    "ppr_last_aa",
    "target_rna_base",
)

# 5′-side solvating helix from Farley et al. GRASP 9S native assemblies
FIVE_PRIME_SOLVATING_HELIX = "QGGNSEEPRKSFDERPERGVVS"

# One ~31-aa PPR-like repeat; only the W?AM 5th and PER? last positions vary
REPEAT_TEMPLATE = "W{fifth}AMISGYAQNGRIDEARELFDKMPER{last}VVS"


def normalize_target_rna(sequence: str) -> str:
    rna = "".join(str(sequence).upper().replace("T", "U").split())
    if not rna:
        raise ValueError("Empty target RNA")
    invalid = sorted(set(rna) - set("ACGU"))
    if invalid:
        raise ValueError(
            "Target RNA contains non-ACGU characters: " + ", ".join(invalid)
        )
    if len(rna) not in (9, 14, 19):
        raise ValueError(
            f"Target RNA length must be 9, 14, or 19 (got {len(rna)}). "
            "GRASP binder scaffolds are defined for those sizes."
        )
    return rna


def rna_to_ppr_pairs(target_rna: str) -> List[Tuple[str, str]]:
    """Return [(5th, last), ...] recognition pairs for each RNA base."""
    rna = normalize_target_rna(target_rna)
    return [_RNA_TO_CODE[b] for b in rna]


def rna_to_binder_aa(target_rna: str) -> str:
    """
    Build the continuous GRASP binder protein for a target RNA.

    Does **not** use combinatorial library modules — only the PPR code and
    the validated GRASP repeat scaffold (matches oh-bounded native 9S assemblies).
    """
    pairs = rna_to_ppr_pairs(target_rna)
    parts = [FIVE_PRIME_SOLVATING_HELIX]
    for fifth, last in pairs:
        parts.append(REPEAT_TEMPLATE.format(fifth=fifth, last=last))
    return "".join(parts)


def ppr_code_string(target_rna: str) -> str:
    """Classic concatenated code, e.g. UUACACGUG → NDNDTNNNTNNNTDNDTD."""
    return "".join(f"{a}{b}" for a, b in rna_to_ppr_pairs(target_rna))


def describe_binder(target_rna: str) -> dict:
    rna = normalize_target_rna(target_rna)
    aa = rna_to_binder_aa(rna)
    pairs = rna_to_ppr_pairs(rna)
    return {
        "target_rna": rna,
        "n_bases": len(rna),
        "ppr_code": ppr_code_string(rna),
        "ppr_pairs": [f"{a}{b}" for a, b in pairs],
        "aa_sequence": aa,
        "aa_length": len(aa),
        "cds_length": 3 * len(aa),
        "n_repeats": len(pairs),
    }


def assembly_role_from_part_id(part_id: str) -> str:
    """Return the five-part Golden Gate slot A/B/C/D/E from a GRASP part id."""
    prefix = str(part_id).split("_", 1)[0]
    if prefix in _JOIN:
        return prefix
    if prefix.endswith("A"):
        return "A"
    if prefix.endswith("E"):
        return "E"
    return ""


def recognition_residues_from_part_id(part_id: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (last AA, 5th AA) encoded in a deposited GRASP part id."""
    match = re.search(r"_L([DN])5([NT])(?:_|$)", str(part_id))
    if match:
        return match.group(1), match.group(2)
    first = re.search(r"_5([NT])(?:_|$)", str(part_id))
    if first:
        return None, first.group(1)
    last = re.search(r"_L([DN])$", str(part_id))
    if last:
        return last.group(1), None
    return None, None


def target_rna_base_from_residues(
    fifth: Optional[str],
    last: Optional[str],
) -> str:
    """RNA base specified by a PPR (5th, last) pair; partial codes stay readable."""
    if fifth and last:
        return _CODE_TO_RNA.get((fifth, last), "")
    if fifth == "T":
        return "A or G"
    if fifth == "N":
        return "C or U"
    if last == "D":
        return "G or U"
    if last == "N":
        return "A or C"
    return ""


def describe_part_id(part_id: str) -> Dict[str, str]:
    """Label a library module by Golden Gate slot and PPR recognition code."""
    role = assembly_role_from_part_id(part_id)
    upstream, downstream = _JOIN.get(role, ("", ""))
    last, fifth = recognition_residues_from_part_id(part_id)
    return {
        "assembly_role": role,
        "joins_upstream_role": upstream,
        "joins_downstream_role": downstream,
        "ppr_5th_aa": fifth or "",
        "ppr_last_aa": last or "",
        "target_rna_base": target_rna_base_from_residues(fifth, last),
    }


def annotate_module_roles(df: pd.DataFrame) -> pd.DataFrame:
    """Add A–E slot and RNA-target columns next to part_id."""
    if df is None or len(df) == 0 or "part_id" not in df.columns:
        return df
    out = df.copy()
    extra = pd.DataFrame(
        [describe_part_id(str(part_id)) for part_id in out["part_id"]],
        index=out.index,
    )
    for column in MODULE_ROLE_COLUMNS:
        out[column] = extra[column]
    lead = [
        column
        for column in ("order_fragment_id", "optimized_part_id", "part_id")
        if column in out.columns
    ]
    rest = [
        column
        for column in out.columns
        if column not in lead and column not in MODULE_ROLE_COLUMNS
    ]
    return out[lead + list(MODULE_ROLE_COLUMNS) + rest]
