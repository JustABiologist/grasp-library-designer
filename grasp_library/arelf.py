"""Movable synonymous GRASP junctions anchored inside the ARELF motif.

GRASP Level -1 parts overlap at an invariant ``ARELF`` sequence.  A Type IIS
fusion boundary may start at any of the twelve nucleotide positions for which
the complete four-base overhang remains inside that 15-nt motif.  This module
builds those candidates and re-materializes the clipped Level -1 part windows
after a boundary position has been selected.
"""

from __future__ import annotations

import math
import re
from itertools import product
from typing import Dict, Mapping, Optional, Sequence, Tuple

import pandas as pd
from Bio.Data import CodonTable

from .binder import FIVE_PRIME_SOLVATING_HELIX
from .dna import clean_dna, reverse_complement

ARELF_MOTIF = "ARELF"
ARELF_NT_LENGTH = 15
ARELF_OFFSETS = tuple(range(12))

INTERNAL_JUNCTIONS = ("J_ACTC", "J_AAGA", "J_GCAC", "J_TGAA")
NATIVE_INTERNAL_CHOICES = {
    "J_ACTC": ("ACTC", 8),
    "J_AAGA": ("AAGA", 2),
    "J_GCAC": ("GCAC", 0),
    "J_TGAA": ("TGAA", 5),
}

DEFAULT_BLOCK_INTERFACES = {
    "cds1_to_cds14": {
        "upstream_three_prime_end_overhang": "TCAC",
        "downstream_five_prime_end_overhang": "GTGA",
        "assembled_coding_site": "GTGA",
        "arelf_offset_nt": 4,
    },
    "cds14_to_cds19": {
        "upstream_three_prime_end_overhang": "CGTG",
        "downstream_five_prime_end_overhang": "CACG",
        "assembled_coding_site": "CACG",
        "arelf_offset_nt": 1,
    },
    "cds1_to_cds2": {
        "upstream_three_prime_end_overhang": "GAAG",
        "downstream_five_prime_end_overhang": "CTTC",
        "assembled_coding_site": "CTTC",
        "arelf_offset_nt": 11,
    },
}

_LEFT_SUFFIX = "ARELFDKMPER{last}VVS"
_RIGHT_PREFIX = "W{fifth}AMISGYAQNGRIDEARELF"

_ROLE_BOUNDARIES = {
    "B": ("J_ACTC", "J_AAGA"),
    "C": ("J_AAGA", "J_GCAC"),
    "D": ("J_GCAC", "J_TGAA"),
    "1E": ("J_TGAA", "cds1_to_cds2"),
    "2A": ("cds1_to_cds2", "J_ACTC"),
    "14E": ("J_TGAA", "cds1_to_cds14"),
    "14A": ("cds1_to_cds14", "J_ACTC"),
    "19E": ("J_TGAA", "cds14_to_cds19"),
    "19A": ("cds14_to_cds19", "J_ACTC"),
}


def _dna4(value: str, *, label: str = "overhang") -> str:
    result = clean_dna(value)
    if len(result) != 4:
        raise ValueError(f"{label} must be exactly four DNA bases")
    return result


def format_cut_token(overhang: str, motif_offset_nt: int) -> str:
    """Serialize a movable cut without losing its motif-relative position."""
    overhang = _dna4(overhang)
    offset = int(motif_offset_nt)
    if offset not in ARELF_OFFSETS:
        raise ValueError("ARELF cut offset must be between 0 and 11")
    return f"{overhang}@{offset}"


def parse_cut_token(
    value: str,
    *,
    default_offset: Optional[int] = None,
) -> Tuple[str, Optional[int]]:
    """Return ``(overhang, offset)``; plain legacy overhangs remain valid."""
    text = str(value).strip().upper().replace("U", "T")
    match = re.fullmatch(r"([ACGT]{4})(?:@(?:ARELF\+)?(\d{1,2}))?", text)
    if not match:
        raise ValueError(
            f"Invalid cut choice {value!r}; expected OVERHANG or OVERHANG@OFFSET"
        )
    overhang = match.group(1)
    offset = int(match.group(2)) if match.group(2) is not None else default_offset
    if offset is not None and offset not in ARELF_OFFSETS:
        raise ValueError("ARELF cut offset must be between 0 and 11")
    return overhang, offset


def selection_overhangs(selection: Mapping[str, str]) -> Dict[str, str]:
    """Strip movable-cut position metadata for ligation calculations."""
    return {str(j): parse_cut_token(v)[0] for j, v in selection.items()}


def _codons_by_aa(
    codon_data: Optional[Mapping[str, Sequence[Mapping]]] = None,
) -> Dict[str, list[str]]:
    if codon_data is None:
        table = CodonTable.unambiguous_dna_by_id[1].forward_table
        return {
            aa: sorted(codon for codon, encoded in table.items() if encoded == aa)
            for aa in ARELF_MOTIF
        }
    result: Dict[str, list[str]] = {}
    for aa in ARELF_MOTIF:
        if aa not in codon_data:
            raise ValueError(f"Codon table has no entries for ARELF residue {aa!r}")
        codons = set()
        for item in codon_data[aa]:
            codon = clean_dna(item["codon"])
            if len(codon) != 3:
                raise ValueError(f"Invalid codon {codon!r} for ARELF residue {aa!r}")
            codons.add(codon)
        codons = sorted(codons)
        if not codons:
            raise ValueError(f"Codon table has no codons for ARELF residue {aa!r}")
        result[aa] = codons
    return result


def synonymous_arelf_overhangs(
    codon_data: Optional[Mapping[str, Sequence[Mapping]]] = None,
    *,
    offsets: Sequence[int] = ARELF_OFFSETS,
) -> Dict[int, set[str]]:
    """All four-mers obtainable by synonymous encoding at each ARELF offset."""
    codons = _codons_by_aa(codon_data)
    encoded = [
        "".join(items)
        for items in product(*(codons[aa] for aa in ARELF_MOTIF))
    ]
    result: Dict[int, set[str]] = {}
    for raw_offset in offsets:
        offset = int(raw_offset)
        if offset not in ARELF_OFFSETS:
            raise ValueError("ARELF candidate offsets must be between 0 and 11")
        result[offset] = {sequence[offset : offset + 4] for sequence in encoded}
    return result


def build_arelf_candidates(
    codon_data: Optional[Mapping[str, Sequence[Mapping]]] = None,
    *,
    junctions: Sequence[str] = INTERNAL_JUNCTIONS,
    offsets: Sequence[int] = ARELF_OFFSETS,
) -> pd.DataFrame:
    """Build runtime candidates whose identity retains overhang and cut offset."""
    by_offset = synonymous_arelf_overhangs(codon_data, offsets=offsets)
    rows = []
    for junction in junctions:
        native_overhang, native_offset = NATIVE_INTERNAL_CHOICES.get(
            str(junction), (None, None)
        )
        for offset in sorted(by_offset):
            for overhang in sorted(by_offset[offset]):
                rows.append(
                    {
                        "junction": str(junction),
                        "overhang": overhang,
                        "motif": ARELF_MOTIF,
                        "motif_offset_nt": offset,
                        "phase": offset % 3,
                        "cut_token": format_cut_token(overhang, offset),
                        "native": int(
                            overhang == native_overhang and offset == native_offset
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _block_interfaces(config: Optional[Mapping]) -> Dict[str, dict]:
    from .assembly_interfaces import resolve_assembly_interfaces

    canonical = resolve_assembly_interfaces(config)["junctions"]
    result: Dict[str, dict] = {}
    for name, defaults in DEFAULT_BLOCK_INTERFACES.items():
        canonical_name = "terminal_to_cds2" if name == "cds1_to_cds2" else name
        canonical_item = dict(canonical.get(canonical_name, {}))
        item = {
            **defaults,
            **canonical_item,
        }
        upstream = _dna4(
            item["upstream_three_prime_end_overhang"],
            label=f"{name} upstream 3′ overhang",
        )
        downstream = _dna4(
            item["downstream_five_prime_end_overhang"],
            label=f"{name} downstream 5′ overhang",
        )
        assembled = _dna4(
            item["assembled_coding_site"],
            label=f"{name} assembled coding-strand site",
        )
        if downstream != reverse_complement(upstream):
            raise ValueError(
                f"{name}: downstream directional overhang {downstream} must equal "
                f"reverse complement {reverse_complement(upstream)} of {upstream}"
            )
        offset = int(item["arelf_offset_nt"])
        if offset not in ARELF_OFFSETS:
            raise ValueError(f"{name}: ARELF cut offset must be between 0 and 11")
        result[name] = {
            "overhang": assembled,
            "offset": offset,
            "upstream_three_prime_end_overhang": upstream,
            "downstream_five_prime_end_overhang": downstream,
            "assembled_coding_site": assembled,
        }
    return result


def _recognition_residues(part_id: str) -> Tuple[Optional[str], Optional[str]]:
    """Return the left-last and right-fifth recognition residues for a part."""
    match = re.search(r"_L([DN])5([NT])(?:_|$)", part_id)
    if match:
        return match.group(1), match.group(2)
    first = re.search(r"_5([NT])(?:_|$)", part_id)
    if first:
        return None, first.group(1)
    last = re.search(r"_L([DN])$", part_id)
    if last:
        return last.group(1), None
    raise ValueError(f"Cannot parse GRASP recognition residues from {part_id!r}")


def _choice_for(
    name: str,
    selection: Mapping[str, str],
    interfaces: Mapping[str, Mapping],
    allowed_offsets: set[int],
    synonymous: Mapping[int, set[str]],
) -> Tuple[str, int]:
    if name in INTERNAL_JUNCTIONS:
        native_overhang, native_offset = NATIVE_INTERNAL_CHOICES[name]
        value = selection.get(name, native_overhang)
        overhang, offset = parse_cut_token(value, default_offset=native_offset)
    else:
        overhang = str(interfaces[name]["overhang"])
        offset = int(interfaces[name]["offset"])
    if name in INTERNAL_JUNCTIONS and offset not in allowed_offsets:
        raise ValueError(f"{name}: ARELF offset {offset} is disabled by configuration")
    if overhang not in synonymous[offset]:
        raise ValueError(
            f"{name}: {overhang}@{offset} cannot synonymously encode ARELF"
        )
    return overhang, offset


def _put(mask: list[str], start: int, overhang: str, *, part_id: str) -> None:
    for index, base in enumerate(overhang):
        position = start + index
        if not 0 <= position < len(mask):
            raise ValueError(
                f"{part_id}: overhang lies outside materialized coding window"
            )
        if mask[position] != "N" and mask[position] != base:
            raise ValueError(f"{part_id}: incompatible overlapping overhang constraints")
        mask[position] = base


def materialize_arelf_parts(
    parts_df: pd.DataFrame,
    selection: Mapping[str, str],
    *,
    config: Optional[Mapping] = None,
    codon_data: Optional[Mapping[str, Sequence[Mapping]]] = None,
) -> pd.DataFrame:
    """Rebuild all GRASP coding windows for movable ARELF fusion boundaries."""
    allowed = set(
        int(x)
        for x in (config or {})
        .get("overhang_redesign", {})
        .get("allowed_arelf_offsets_nt", ARELF_OFFSETS)
    )
    if not allowed or not allowed <= set(ARELF_OFFSETS):
        raise ValueError("Allowed ARELF offsets must be a non-empty subset of 0..11")
    synonymous = synonymous_arelf_overhangs(codon_data)
    interfaces = _block_interfaces(config)
    chosen = {
        name: _choice_for(name, selection, interfaces, allowed, synonymous)
        for name in (*INTERNAL_JUNCTIONS, *DEFAULT_BLOCK_INTERFACES)
    }

    rows = []
    for source in parts_df.to_dict(orient="records"):
        row = dict(source)
        part_id = str(row["part_id"])
        role = part_id.split("_", 1)[0]
        last, fifth = _recognition_residues(part_id)

        if role == "1A":
            current_aa = str(row["aa_sequence"]).upper().replace(" ", "")
            w_index = current_aa.find("W")
            five_prime_prefix = (
                current_aa[:w_index]
                if w_index >= 0
                else FIVE_PRIME_SOLVATING_HELIX
            )
            full_aa = five_prime_prefix + _RIGHT_PREFIX.format(fifth=fifth)
            right_anchor_nt = 3 * (len(five_prime_prefix) + 16)
            left = (_dna4(part_id.rsplit("_", 1)[-1]), 2, "J_5prime", None)
            overhang, offset = chosen["J_ACTC"]
            right = (overhang, right_anchor_nt + offset, "J_ACTC", offset)
        elif role == "2E":
            full_aa = _LEFT_SUFFIX.format(last=last)
            overhang, offset = chosen["J_TGAA"]
            left = (overhang, offset, "J_TGAA", offset)
            right = ("TTCG", 3 * len(full_aa) - 4, "J_3prime", None)
        else:
            if role not in _ROLE_BOUNDARIES or last is None or fifth is None:
                raise ValueError(f"Unsupported GRASP part role in {part_id!r}")
            full_aa = _LEFT_SUFFIX.format(last=last) + _RIGHT_PREFIX.format(
                fifth=fifth
            )
            left_name, right_name = _ROLE_BOUNDARIES[role]
            left_oh, left_offset = chosen[left_name]
            right_oh, right_offset = chosen[right_name]
            left = (left_oh, left_offset, left_name, left_offset)
            right_anchor_nt = 3 * len(_LEFT_SUFFIX.format(last=last)) + 3 * 16
            right = (
                right_oh,
                right_anchor_nt + right_offset,
                right_name,
                right_offset,
            )

        left_oh, left_full_start, left_name, left_offset = left
        right_oh, right_full_start, right_name, right_offset = right
        window_start = (left_full_start // 3) * 3
        window_end = int(math.ceil((right_full_start + 4) / 3.0) * 3)
        if not (0 <= window_start < window_end <= 3 * len(full_aa)):
            raise ValueError(f"{part_id}: invalid materialized coding window")
        aa_sequence = full_aa[window_start // 3 : window_end // 3]
        mask = ["N"] * (window_end - window_start)
        oh5 = left_full_start - window_start
        oh3 = right_full_start - window_start
        _put(mask, oh5, left_oh, part_id=part_id)
        _put(mask, oh3, right_oh, part_id=part_id)

        row.update(
            {
                "aa_sequence": aa_sequence,
                "coding_mask": "".join(mask),
                "oh5": left_oh,
                "oh3": right_oh,
                "oh5_coding_site_5to3": left_oh,
                "oh3_coding_site_5to3": right_oh,
                # Coding-site labels and physical end labels remain separate.
                "five_prime_end_overhang": (
                    interfaces[left_name]["downstream_five_prime_end_overhang"]
                    if left_name in interfaces
                    else left_oh
                ),
                "three_prime_end_overhang": (
                    interfaces[right_name]["upstream_three_prime_end_overhang"]
                    if right_name in interfaces
                    else reverse_complement(right_oh)
                ),
                "overhang_notation": "assembled_coding_strand_site",
                "oh5_mask_start": oh5,
                "oh3_mask_start": oh3,
                "oh5_junction": left_name,
                "oh3_junction": right_name,
                "oh5_arelf_offset_nt": left_offset,
                "oh3_arelf_offset_nt": right_offset,
                "cut_mode": "movable_arelf",
                "full_aa_sequence": full_aa,
                "full_window_start_nt": window_start,
                "full_window_end_nt": window_end,
            }
        )
        rows.append(row)
    from .binder import annotate_module_roles

    return annotate_module_roles(pd.DataFrame(rows))


def dynamic_junction_map(parts_df: pd.DataFrame) -> pd.DataFrame:
    """Describe internal junction coordinates after ARELF materialization."""
    rows = []
    for part in parts_df.itertuples(index=False):
        for end, prefix in (("5prime", "oh5"), ("3prime", "oh3")):
            junction = str(getattr(part, f"{prefix}_junction", ""))
            if junction not in INTERNAL_JUNCTIONS:
                continue
            rows.append(
                {
                    "junction": junction,
                    "part_id": str(part.part_id),
                    "mask_start_0based": int(getattr(part, f"{prefix}_mask_start")),
                    "native_overhang": str(getattr(part, prefix)),
                    "end": end,
                    "motif": ARELF_MOTIF,
                    "motif_offset_nt": int(
                        getattr(part, f"{prefix}_arelf_offset_nt")
                    ),
                }
            )
    return pd.DataFrame(rows)
