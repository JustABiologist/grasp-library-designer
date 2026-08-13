"""Shared DNA / mask helpers used by objectives and Pareto search."""

from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache
from typing import Dict, Iterable, List

from Bio.Seq import Seq

DNA_ALPHABET = set("ACGT")
MASK_ALPHABET = set("ACGTN")
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
_IUPAC_TO_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def clean_dna(sequence: str) -> str:
    sequence = str(sequence).upper().replace("U", "T")
    sequence = re.sub(r"\s+", "", sequence)
    invalid = set(sequence) - DNA_ALPHABET
    if invalid:
        raise ValueError(f"Invalid DNA characters: {invalid}")
    return sequence


def clean_mask(mask: str) -> str:
    mask = str(mask).upper().replace("U", "T")
    mask = re.sub(r"\s+", "", mask)
    invalid = set(mask) - MASK_ALPHABET
    if invalid:
        raise ValueError(f"Invalid mask characters: {invalid}")
    return mask


def reverse_complement(sequence: str) -> str:
    return str(Seq(clean_dna(sequence)).reverse_complement())


def translate_dna(sequence: str, genetic_code: int = 1) -> str:
    sequence = clean_dna(sequence)
    if len(sequence) % 3 != 0:
        raise ValueError("CDS length is not divisible by three.")
    return str(Seq(sequence).translate(table=genetic_code))


def gc_fraction(sequence: str) -> float:
    sequence = clean_dna(sequence)
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def mask_matches(sequence: str, mask: str) -> bool:
    sequence = clean_dna(sequence)
    mask = clean_mask(mask)
    if len(sequence) != len(mask):
        return False
    return all(
        mask_base == "N" or sequence_base == mask_base
        for sequence_base, mask_base in zip(sequence, mask)
    )


def longest_homopolymer(sequence: str) -> int:
    sequence = clean_dna(sequence)
    return max(
        (len(match.group(0)) for match in re.finditer(r"(A+|C+|G+|T+)", sequence)),
        default=0,
    )


def kmer_counts(sequence: str, k: int) -> Counter:
    sequence = clean_dna(sequence)
    if len(sequence) < k:
        return Counter()
    return Counter(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def repeated_kmer_penalty(sequence: str, k: int = 12) -> float:
    counts = kmer_counts(sequence, k)
    return float(sum((count - 1) ** 2 for count in counts.values() if count > 1))


def local_gc_penalty(
    sequence: str,
    window_size: int,
    gc_min: float,
    gc_max: float,
) -> float:
    sequence = clean_dna(sequence)
    if len(sequence) <= window_size:
        value = gc_fraction(sequence)
        if value < gc_min:
            return (gc_min - value) ** 2
        if value > gc_max:
            return (value - gc_max) ** 2
        return 0.0

    penalty = 0.0
    for start in range(len(sequence) - window_size + 1):
        window = sequence[start : start + window_size]
        value = gc_fraction(window)
        if value < gc_min:
            penalty += (gc_min - value) ** 2
        elif value > gc_max:
            penalty += (value - gc_max) ** 2
    return penalty


def clean_iupac_dna(sequence: str) -> str:
    sequence = str(sequence).upper().replace("U", "T")
    sequence = re.sub(r"\s+", "", sequence)
    invalid = set(sequence) - IUPAC_DNA
    if invalid:
        raise ValueError(f"Invalid IUPAC DNA characters: {invalid}")
    return sequence


def reverse_complement_iupac(sequence: str) -> str:
    return str(Seq(clean_iupac_dna(sequence)).reverse_complement())


@lru_cache(maxsize=256)
def _iupac_motif_regex(motif: str) -> re.Pattern[str]:
    return re.compile("".join(_IUPAC_TO_REGEX[base] for base in motif))


def contains_forbidden_site(
    sequence: str,
    forbidden_sites: Dict[str, str],
) -> List[dict]:
    sequence = clean_dna(sequence)
    hits = []
    for name, motif in forbidden_sites.items():
        motif = clean_iupac_dna(motif)
        queries = {motif, reverse_complement_iupac(motif)}
        for query in queries:
            if set(query) <= DNA_ALPHABET:
                start = sequence.find(query)
                while start != -1:
                    hits.append(
                        {
                            "enzyme": name,
                            "site": query,
                            "start_0based": start,
                        }
                    )
                    start = sequence.find(query, start + 1)
                continue
            for match in _iupac_motif_regex(query).finditer(sequence):
                hits.append(
                    {
                        "enzyme": name,
                        "site": match.group(0),
                        "start_0based": match.start(),
                    }
                )
    return hits


def is_self_reverse_complement(overhang: str) -> bool:
    overhang = clean_dna(overhang)
    return overhang == reverse_complement(overhang)


def apply_overhang_to_mask(mask: str, start_0based: int, overhang: str) -> str:
    """Write a fixed overhang into a coding mask."""
    mask = clean_mask(mask)
    overhang = clean_dna(overhang)
    if start_0based < 0 or start_0based + len(overhang) > len(mask):
        raise ValueError("Overhang lies outside the mask.")
    chars = list(mask)
    for offset, base in enumerate(overhang):
        chars[start_0based + offset] = base
    return "".join(chars)
