"""RNA → GRASP binder protein (PPR recognition code), no library parts required."""

from __future__ import annotations

from typing import List, Sequence, Tuple

# Classic PPR code: (5th AA, last AA) of each repeat
_RNA_TO_CODE = {
    "A": ("T", "N"),
    "C": ("N", "N"),
    "G": ("T", "D"),
    "U": ("N", "D"),
    "T": ("N", "D"),
}

# Solvating N-helix from Farley et al. GRASP 9S native assemblies
NTERM_HELIX = "QGGNSEEPRKSFDERPERGVVS"

# One ~31-aa PPR-like repeat; only the W?AM 5th and PER? last positions vary
REPEAT_TEMPLATE = "W{fifth}AMISGYAQNGRIDEARELFDKMPER{last}VVS"


def normalize_target_rna(sequence: str) -> str:
    rna = str(sequence).upper().replace("T", "U")
    rna = "".join(b for b in rna if b in "ACGU")
    if not rna:
        raise ValueError("Empty target RNA")
    if len(rna) not in (9, 14, 19):
        raise ValueError(
            f"Target RNA length must be 9, 14, or 19 (got {len(rna)}). "
            "GRASP binder scaffolds are defined for those sizes."
        )
    if any(b not in _RNA_TO_CODE for b in rna):
        raise ValueError(f"Non-ACGU base in target RNA: {sequence!r}")
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
    parts = [NTERM_HELIX]
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
