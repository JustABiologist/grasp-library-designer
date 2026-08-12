"""Import GRASP level −1 modules from deposited GenBank files (Farley et al., NAR 2025).

Source: https://github.com/farleykvdg/GRASP (GRASP_-1.gb / Assembly Planner Seqs)
Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.Data import CodonTable
from Bio.Seq import Seq

warnings.simplefilter("ignore", BiopythonParserWarning)

from .assembly_interfaces import (
    build_order_fragment,
    deposited_grasp_interface_preset,
    order_fragment_arms,
    resolve_assembly_interfaces,
    reverse_complement,
)

OVERHANG_SET = {
    "AGGT",
    "AATG",
    "ACTC",
    "AAGA",
    "GCAC",
    "TGAA",
    "CTTC",
    "TTCG",
    "GTGA",
    "CACG",
}

# pAGM1311 is the universal Level -1 entry vector used by GRASP.  An ordered
# dsDNA fragment is inserted through the vector's BsaI sites with ACAT/TTGT
# fusion sites.  The BsaI sites below are on the ordered fragment and disappear
# after cloning; the retained pAGM1311 backbone then completes flanking BpiI
# sites used to release the GRASP module into pAGM9121.
_DEPOSITED_INTERFACES = deposited_grasp_interface_preset()
_DEPOSITED_PREFIX, _DEPOSITED_SUFFIX = order_fragment_arms(_DEPOSITED_INTERFACES)
PAGM1311_ORDER_5P_ARM = _DEPOSITED_PREFIX[:-4]
PAGM1311_ORDER_3P_ARM = _DEPOSITED_SUFFIX[4:]
PAGM1311_INSERT_5P_FUSION = _DEPOSITED_INTERFACES["level_minus1_entry"]["n_overhang_5p"]
PAGM1311_INSERT_3P_FUSION = reverse_complement(
    _DEPOSITED_INTERFACES["level_minus1_entry"]["c_overhang_5p"]
)
PAGM1311_BACKBONE_5P_CONTEXT = _DEPOSITED_INTERFACES["level_minus1_entry"]["completion_context_5p"]
PAGM1311_BACKBONE_3P_CONTEXT = _DEPOSITED_INTERFACES["level_minus1_entry"]["completion_context_3p"]
PAGM9121_EXTERNAL_5P_OVERHANG = _DEPOSITED_INTERFACES["level0"]["acceptor_outer"]["n_overhang_5p"]
PAGM9121_EXTERNAL_3P_OVERHANG = _DEPOSITED_INTERFACES["level0"]["acceptor_outer"]["c_overhang_5p"]


def build_configured_order_fragment(
    optimized_cds: str,
    *,
    part_id: str,
    oh5_mask_start: int,
    oh3_mask_start: int,
    config: Optional[dict] = None,
    interfaces: Optional[dict] = None,
) -> str:
    """Build an order fragment for the configured entry/assembly profile."""
    profile = interfaces or resolve_assembly_interfaces(config)
    cds = re.sub(r"\s+", "", str(optimized_cds).upper().replace("U", "T"))
    if set(cds) - set("ACGT"):
        raise ValueError(f"{part_id}: invalid DNA in optimized CDS")
    oh5 = int(oh5_mask_start)
    oh3 = int(oh3_mask_start)
    if not (0 <= oh5 < oh3 + 4 <= len(cds)):
        raise ValueError(
            f"{part_id}: invalid overhang coordinates oh5={oh5}, oh3={oh3}, "
            f"CDS length={len(cds)}"
        )
    payload = cds[oh5 : oh3 + 4]
    outer = profile["level0"].get("acceptor_outer")
    role = str(part_id).split("_", 1)[0]
    if outer is not None and role.endswith("A"):
        payload = outer["n_overhang_5p"] + payload
    if outer is not None and role.endswith("E"):
        payload = payload + outer["c_overhang_5p"]
    return build_order_fragment(payload, profile)


def build_pagm1311_order_fragment(
    optimized_cds: str,
    *,
    part_id: str,
    oh5_mask_start: int,
    oh3_mask_start: int,
) -> str:
    """Build a dsDNA order sequence that survives both GRASP cloning stages.

    Only the overhang-bounded module is released by BpiI.  A/E parts also carry
    the fixed pAGM9121 external fusion site outside their coding overhang.  This
    coordinate-based construction prevents mutable in-frame padding from
    corrupting the pAGM1311 ACAT/TTGT entry-vector fusion sites.
    """
    return build_configured_order_fragment(
        optimized_cds,
        part_id=part_id,
        oh5_mask_start=oh5_mask_start,
        oh3_mask_start=oh3_mask_start,
        interfaces=_DEPOSITED_INTERFACES,
    )

# Coding junctions annotated in the deposited modules.  ACTC/AAGA/GCAC/TGAA
# are the four BpiI-released internal junctions used in every five-part Level 0
# assembly.  AGGT/CTTC/TTCG are later-stage BsaI/MoClo fusion interfaces and are
# deliberately not eligible for redesign.
JUNCTION_ORDER_9S = [
    ("J_Nterm", "AGGT"),  # 1A 5′ (AGGT preferred; AATG alternate)
    ("J_ACTC", "ACTC"),  # 1A/2A 3′ ↔ B 5′
    ("J_AAGA", "AAGA"),  # B 3′ ↔ C 5′
    ("J_GCAC", "GCAC"),  # C 3′ ↔ D 5′
    ("J_TGAA", "TGAA"),  # D 3′ ↔ 1E/2E 5′
    ("J_CTTC", "CTTC"),  # 1E 3′ ↔ 2A 5′
    ("J_Cterm", "TTCG"),  # 2E 3′
]

# Each binding architecture is first assembled as five-part Level 0 blocks,
# which are subsequently joined into the complete sPPR.  The block order mirrors
# GRASP_AP.py: pPR0_1, optionally pPR0_14 and pPR0_19, then pPR0_2.  Keep the
# existing CDS1/CDS2 labels for backwards compatibility with exported plans.
ARCHITECTURE_LAYOUTS = {
    "9S": {
        "target_length": 9,
        "assembly_groups": ("CDS1", "CDS2"),
    },
    "14S": {
        "target_length": 14,
        "assembly_groups": ("CDS1", "CDS14", "CDS2"),
    },
    "19S": {
        "target_length": 19,
        "assembly_groups": ("CDS1", "CDS14", "CDS19", "CDS2"),
    },
}

# GAP code_part ↔ filename suffix (GRASP_AP.partpicker)
CODE_TO_SUFFIX = {
    "T": "_5T",
    "N": "_5N",
    "DN": "_LD5N",
    "DT": "_LD5T",
    "NN": "_LN5N",
    "NT": "_LN5T",
    "D": "D",
    "N_last": "N",
}

STANDARD_CODONS = CodonTable.unambiguous_dna_by_id[1]


def _feature_label(feature) -> str:
    return feature.qualifiers.get("label", [""])[0]


def _parse_record(rec) -> dict:
    overhangs: List[Tuple[int, str]] = []
    cds = None
    binds: List[Tuple[int, str]] = []
    insert = None

    for feature in rec.features:
        start = int(feature.location.start)
        end = int(feature.location.end)
        seq = str(rec.seq[start:end]).upper()
        label = _feature_label(feature)

        if end - start == 4 and start < 220 and seq in OVERHANG_SET:
            overhangs.append((start, seq))

        if feature.type == "CDS" and end - start < 200:
            if cds is None or (end - start) < (cds[1] - cds[0]):
                cds = (start, end, seq)

        if feature.type == "protein_bind" and start < 220:
            binds.append((start, seq))

        if (
            feature.type == "misc_feature"
            and seq.startswith("ACAT")
            and 40 < (end - start) < 200
        ):
            if insert is None or (end - start) > (insert[1] - insert[0]):
                insert = (start, end, seq)

    overhangs = sorted(set(overhangs), key=lambda x: x[0])
    if len(overhangs) < 2:
        raise ValueError(f"{rec.name}: expected 2 overhangs, found {overhangs}")

    return {
        "name": rec.name,
        "seq": str(rec.seq).upper(),
        "overhangs": overhangs,
        "cds": cds,
        "binds": binds,
        "insert": insert,
    }


def _annotated_aa(cds_seq: str) -> str:
    trim = cds_seq[: len(cds_seq) // 3 * 3]
    aa = str(Seq(trim).translate(to_stop=False))
    return aa.replace("*", "")


def _choose_coding_window(parsed: dict) -> dict:
    """Pick an in-frame window that contains both overhangs and keeps the CDS AA."""
    seq = parsed["seq"]
    (o1, oh1), (o2, oh2) = parsed["overhangs"][0], parsed["overhangs"][-1]
    core_aa = _annotated_aa(parsed["cds"][2]) if parsed["cds"] else ""

    candidates = []
    for start in range(max(0, o1 - 9), o1 + 1):
        for end in range(o2 + 4, min(len(seq), o2 + 4 + 9)):
            if (end - start) % 3 != 0:
                continue
            dna = seq[start:end]
            aa = str(Seq(dna).translate(to_stop=False))
            if "*" in aa:
                continue
            if core_aa and core_aa not in aa and aa not in core_aa:
                # allow 2E-style short modules where GenBank CDS annotation is off-frame
                if len(aa) > 20:
                    continue
            extra = abs(len(aa) - len(core_aa)) if core_aa else len(aa)
            # prefer windows where annotated AA is a substring
            contains = 0 if (core_aa and core_aa in aa) else 2
            candidates.append(
                {
                    "start": start,
                    "end": end,
                    "dna": dna,
                    "aa": aa,
                    "oh5_pos": o1 - start,
                    "oh3_pos": o2 - start,
                    "oh5": oh1,
                    "oh3": oh2,
                    "score": (contains, extra, end - start),
                }
            )

    if not candidates:
        raise ValueError(f"{parsed['name']}: no in-frame window covering overhangs")

    best = sorted(candidates, key=lambda c: c["score"])[0]
    return best


def _flanks(parsed: dict, window_start: int, window_end: int) -> Tuple[str, str]:
    """Return order-ready BsaI arms for cloning the module into pAGM1311.

    The deposited plasmids contain the post-cloning product, so the BsaI sites
    used to insert the synthesized fragment are absent.  Reconstruct them around
    the annotated ACAT...TTGT insert rather than accidentally exporting one
    retained BpiI site from the entry-vector backbone.
    """
    seq = parsed["seq"]
    if parsed["insert"] is None:
        raise ValueError(f"{parsed['name']}: no annotated pAGM1311 insert")
    insert_start, insert_end, insert_seq = parsed["insert"]
    if not insert_seq.startswith(PAGM1311_INSERT_5P_FUSION):
        raise ValueError(f"{parsed['name']}: pAGM1311 insert does not start ACAT")
    if not insert_seq.endswith(PAGM1311_INSERT_3P_FUSION):
        raise ValueError(f"{parsed['name']}: pAGM1311 insert does not end TTGT")
    if not (insert_start <= window_start <= window_end <= insert_end):
        raise ValueError(f"{parsed['name']}: coding window lies outside entry insert")

    return PAGM1311_ORDER_5P_ARM, PAGM1311_ORDER_3P_ARM


def _part_id_from_name(name: str) -> str:
    # pPR-1_1A_5N_AGGT → 1A_5N_AGGT
    return re.sub(r"^pPR-1_", "", name)


def load_grasp_records(genbank_paths: Sequence[Path]) -> List[dict]:
    records = []
    for path in genbank_paths:
        for rec in SeqIO.parse(path, "genbank"):
            if not rec.name.startswith("pPR-1_"):
                continue
            parsed = _parse_record(rec)
            window = _choose_coding_window(parsed)
            prefix, suffix = _flanks(parsed, window["start"], window["end"])
            dna = window["dna"]
            mask_chars = ["N"] * len(dna)
            for pos, oh in (
                (window["oh5_pos"], window["oh5"]),
                (window["oh3_pos"], window["oh3"]),
            ):
                for i, base in enumerate(oh):
                    mask_chars[pos + i] = base
            records.append(
                {
                    "part_id": _part_id_from_name(parsed["name"]),
                    "source_name": parsed["name"],
                    "aa_sequence": window["aa"],
                    "coding_dna": dna,
                    "coding_mask": "".join(mask_chars),
                    "oligo_prefix": prefix,
                    "oligo_suffix": suffix,
                    "oh5": window["oh5"],
                    "oh3": window["oh3"],
                    "oh5_mask_start": window["oh5_pos"],
                    "oh3_mask_start": window["oh3_pos"],
                    "genome_start": window["start"],
                    "genome_end": window["end"],
                }
            )
    return records


def overhang_synonyms(
    coding_dna: str,
    aa_sequence: str,
    mask_start: int,
    overhang: str,
) -> List[str]:
    """Synonymous 4-mers at a fixed mask index that preserve the part AA."""
    overhang = overhang.upper().replace("U", "T")
    if coding_dna[mask_start : mask_start + 4] != overhang:
        raise ValueError("Overhang does not match coding DNA at mask_start")

    # Codon window covering the overhang
    codon_start = (mask_start // 3) * 3
    codon_end = ((mask_start + 4 + 2) // 3) * 3
    local_dna = coding_dna[codon_start:codon_end]
    local_aa = aa_sequence[codon_start // 3 : codon_end // 3]
    offset = mask_start - codon_start

    codon_options = []
    for aa in local_aa:
        opts = [c for c, a in STANDARD_CODONS.forward_table.items() if a == aa]
        if not opts:
            # stop or invalid
            return [overhang]
        codon_options.append(opts)

    from itertools import product

    syn = set()
    for codons in product(*codon_options):
        dna = "".join(codons)
        cand = dna[offset : offset + 4]
        if len(cand) == 4:
            syn.add(cand)
    return sorted(syn)


def build_parts_table(records: Sequence[dict]) -> pd.DataFrame:
    rows = [
        {
            "part_id": r["part_id"],
            "aa_sequence": r["aa_sequence"],
            "coding_mask": r["coding_mask"],
            "oligo_prefix": r["oligo_prefix"],
            "oligo_suffix": r["oligo_suffix"],
            "source_name": r["source_name"],
            "native_cds": r["coding_dna"],
            "oh5": r["oh5"],
            "oh3": r["oh3"],
            "oh5_mask_start": r["oh5_mask_start"],
            "oh3_mask_start": r["oh3_mask_start"],
        }
        for r in records
    ]
    return pd.DataFrame(rows).sort_values("part_id").reset_index(drop=True)


def build_junction_map_9s(records: Sequence[dict]) -> pd.DataFrame:
    """Map each 9S junction onto every part that carries that overhang."""
    by_id = {r["part_id"]: r for r in records}
    # Prefer AGGT 1A variants for N-term junction ownership; still list AATG parts.
    rows = []

    def add(junction: str, part_id: str, which: str):
        rec = by_id[part_id]
        start = rec["oh5_mask_start"] if which == "5" else rec["oh3_mask_start"]
        oh = rec["oh5"] if which == "5" else rec["oh3"]
        rows.append(
            {
                "junction": junction,
                "part_id": part_id,
                "mask_start_0based": int(start),
                "native_overhang": oh,
                "end": "5prime" if which == "5" else "3prime",
            }
        )

    for part_id, rec in by_id.items():
        role = part_id.split("_")[0]
        oh5, oh3 = rec["oh5"], rec["oh3"]

        # Prefer AGGT N-terminal fusion (MoClo standard); AATG variants stay in parts.csv
        # but are not part of the default redesign junction.
        if part_id.startswith("1A_") and oh5 == "AGGT":
            add("J_Nterm", part_id, "5")
        if part_id.startswith("2E_") and oh3 == "TTCG":
            add("J_Cterm", part_id, "3")

        # Every A extension feeds the shared B/C/D core.
        if role in {"1A", "2A", "14A", "19A"} and oh3 == "ACTC":
            add("J_ACTC", part_id, "3")
        if role == "B" and oh5 == "ACTC":
            add("J_ACTC", part_id, "5")
        if role == "B" and oh3 == "AAGA":
            add("J_AAGA", part_id, "3")
        if role == "C" and oh5 == "AAGA":
            add("J_AAGA", part_id, "5")
        if role == "C" and oh3 == "GCAC":
            add("J_GCAC", part_id, "3")
        if role == "D" and oh5 == "GCAC":
            add("J_GCAC", part_id, "5")
        if role == "D" and oh3 == "TGAA":
            add("J_TGAA", part_id, "3")
        # Every E extension receives the shared B/C/D core.
        if role in {"1E", "2E", "14E", "19E"} and oh5 == "TGAA":
            add("J_TGAA", part_id, "5")
        if role == "1E" and oh3 == "CTTC":
            add("J_CTTC", part_id, "3")
        if role == "2A" and oh5 == "CTTC":
            add("J_CTTC", part_id, "5")
    df = pd.DataFrame(rows).drop_duplicates()
    return df.sort_values(["junction", "part_id", "end"]).reset_index(drop=True)


def build_overhang_candidates(
    records: Sequence[dict],
    junction_map: pd.DataFrame,
    architecture: str = "9S",
) -> pd.DataFrame:
    by_id = {r["part_id"]: r for r in records}
    rows = []
    junctions = junction_map["junction"].unique()
    for junction in junctions:
        if architecture == "9S" and junction.startswith("J_"):
            # keep 14/19 helpers but mark them; still emit synonyms
            pass
        sub = junction_map[junction_map["junction"] == junction]
        # Intersection of synonym sets across all parts carrying this junction
        syn_sets = []
        native = None
        for row in sub.itertuples(index=False):
            rec = by_id[row.part_id]
            native = row.native_overhang
            syn = overhang_synonyms(
                rec["coding_dna"],
                rec["aa_sequence"],
                int(row.mask_start_0based),
                row.native_overhang,
            )
            syn_sets.append(set(syn))
        if not syn_sets:
            continue
        common = set.intersection(*syn_sets) if len(syn_sets) > 1 else syn_sets[0]
        if native:
            common.add(native)
        for oh in sorted(common):
            rows.append(
                {
                    "junction": junction,
                    "overhang": oh,
                    "native": int(oh == native),
                }
            )
    return pd.DataFrame(rows)


# --- Target / assembly planning (from GRASP_AP.py, Farley et al.) ---

def dna_to_ppr_code(sequence: str) -> str:
    sequence = sequence.upper().replace("U", "T")
    out = []
    for base in sequence:
        if base == "A":
            out.append("TN")
        elif base == "T":
            out.append("ND")
        elif base == "C":
            out.append("NN")
        elif base == "G":
            out.append("TD")
        else:
            raise ValueError(f"Noncanonical base: {base}")
    return "".join(out)


def codify_ppr(ppr_code: str) -> List[str]:
    """Split concatenated 5th/last code into GAP duplex list."""
    if len(ppr_code) > 14:
        trim = ppr_code[1:-1]
        duplex = [trim[i : i + 2] for i in range(0, len(trim), 2)]
        duplex.insert(0, ppr_code[0])
        duplex.append(ppr_code[-1])
    else:
        duplex = [ppr_code[i : i + 2] for i in range(0, len(ppr_code), 2)]
    return duplex


def part_prefixes(length: int) -> List[str]:
    if length == 10:  # 9-motif → 10 parts
        return [
            "1A",
            "B",
            "C",
            "D",
            "1E",
            "2A",
            "B",
            "C",
            "D",
            "2E_L",
        ]
    if length == 15:
        return [
            "1A",
            "B",
            "C",
            "D",
            "14E",
            "14A",
            "B",
            "C",
            "D",
            "1E",
            "2A",
            "B",
            "C",
            "D",
            "2E_L",
        ]
    if length == 20:
        return [
            "1A",
            "B",
            "C",
            "D",
            "14E",
            "14A",
            "B",
            "C",
            "D",
            "19E",
            "19A",
            "B",
            "C",
            "D",
            "1E",
            "2A",
            "B",
            "C",
            "D",
            "2E_L",
        ]
    raise ValueError(f"Unsupported duplex length {length}")


def pick_parts_for_target(
    target_rna: str,
    *,
    nterm_overhang: str = "AGGT",
) -> List[str]:
    """Return part_id list for a 9/14/19-base binding target (GAP logic)."""
    target_rna = target_rna.upper().replace("T", "U")
    dna = target_rna.replace("U", "T")
    n_motif = len(dna)
    if n_motif not in {9, 14, 19}:
        raise ValueError("Binding targets must be 9, 14, or 19 bases")

    ppr = dna_to_ppr_code(dna)
    duplex = codify_ppr(ppr)
    prefixes = part_prefixes(len(duplex))
    code_part = ["T", "N", "DN", "DT", "NN", "NT"]
    code_compute = ["_5T", "_5N", "_LD5N", "_LD5T", "_LN5N", "_LN5T"]

    part_ids = []
    for i, code in enumerate(duplex):
        prefix = prefixes[i]
        if i == len(duplex) - 1:
            # 2E_L{D|N}
            suffix = "D" if code == "D" else "N"
            # deposited files are 2E_LD / 2E_LN
            part_ids.append(f"2E_L{suffix}")
            continue
        if code not in code_part:
            raise ValueError(f"Unexpected GAP code {code!r} at slot {i}")
        suffix = code_compute[code_part.index(code)]
        if prefix == "1A":
            part_ids.append(f"1A{suffix}_{nterm_overhang}")
        else:
            part_ids.append(f"{prefix}{suffix}")
    return part_ids


def build_target_map_catalog(records: Sequence[dict]) -> pd.DataFrame:
    """Catalog rows for notebook display: architecture × slot × code → part_id."""
    rows = []
    # 9S slots
    slot_meta = [
        (1, "1A", "CDS1", 1, ["_5T", "_5N"]),
        (2, "B", "CDS1", 2, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (3, "C", "CDS1", 3, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (4, "D", "CDS1", 4, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (5, "1E", "CDS1", 5, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (6, "2A", "CDS2", 1, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (7, "B", "CDS2", 2, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (8, "C", "CDS2", 3, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (9, "D", "CDS2", 4, ["_LD5N", "_LD5T", "_LN5N", "_LN5T"]),
        (10, "2E", "CDS2", 5, ["_LD", "_LN"]),
    ]
    by_id = {r["part_id"] for r in records}

    # Map GAP first/last half-codes and middle codes onto RNA bases they help specify.
    # For notebook compatibility keep target_base as the RNA base preferentially
    # associated with each full TN/NN/TD/ND code; 1A/2E use half-codes.
    code_to_bases = {
        "_5T": ["A", "G"],  # 5th=T
        "_5N": ["C", "U"],  # 5th=N
        "_LD5N": ["U"],  # last D, 5th N → ND
        "_LD5T": ["G"],  # last D, 5th T → TD
        "_LN5N": ["C"],  # last N, 5th N → NN
        "_LN5T": ["A"],  # last N, 5th T → TN
        "_LD": ["U", "G"],  # last D
        "_LN": ["A", "C"],  # last N
    }

    for slot, role, group, order, suffixes in slot_meta:
        for suffix in suffixes:
            bases = code_to_bases[suffix]
            if role == "1A":
                for nt_oh in ("AGGT", "AATG"):
                    part_id = f"1A{suffix}_{nt_oh}"
                    if part_id not in by_id:
                        continue
                    for base in bases:
                        rows.append(
                            {
                                "architecture": "9S",
                                "target_position": slot,
                                "target_base": base,
                                "part_id": part_id,
                                "assembly_group": group,
                                "assembly_order": order,
                                "code_suffix": suffix,
                                "note": "half-code 5th AA; full motif depends on next part",
                            }
                        )
            elif role == "2E":
                part_id = f"2E{suffix}"
                if part_id not in by_id:
                    continue
                for base in bases:
                    rows.append(
                        {
                            "architecture": "9S",
                            "target_position": slot,
                            "target_base": base,
                            "part_id": part_id,
                            "assembly_group": group,
                            "assembly_order": order,
                            "code_suffix": suffix,
                            "note": "half-code last AA; full motif depends on previous part",
                        }
                    )
            else:
                part_id = f"{role}{suffix}"
                if part_id not in by_id:
                    continue
                for base in bases:
                    rows.append(
                        {
                            "architecture": "9S",
                            "target_position": slot,
                            "target_base": base,
                            "part_id": part_id,
                            "assembly_group": group,
                            "assembly_order": order,
                            "code_suffix": suffix,
                            "note": "GAP Last+5th code; selected from full target via compile_target_gap",
                        }
                    )
    return pd.DataFrame(rows)


def import_grasp_profile(
    genbank_dir: Path,
    output_dir: Path,
    *,
    prefer_multirecord: bool = True,
) -> Dict[str, pd.DataFrame]:
    genbank_dir = Path(genbank_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    multi = genbank_dir / "GRASP_-1.gb"
    if prefer_multirecord and multi.exists():
        paths = [multi]
    else:
        paths = sorted(genbank_dir.glob("pPR-1_*.gb"))

    records = load_grasp_records(paths)
    parts = build_parts_table(records)
    junction_map = build_junction_map_9s(records)
    # Only the four physical BpiI subassembly junctions may be redesigned.
    # MoClo interfaces AGGT/CTTC/TTCG remain native so the resulting Level 0
    # constructs stay compatible with standard flanking domains and vectors.
    redesignable = {"J_ACTC", "J_AAGA", "J_GCAC", "J_TGAA"}
    junction_9s = junction_map[
        junction_map["junction"].isin(redesignable)
    ].copy()
    candidates = build_overhang_candidates(records, junction_9s)
    target_map = build_target_map_catalog(records)

    parts_out = parts[
        [
            "part_id",
            "aa_sequence",
            "coding_mask",
            "oligo_prefix",
            "oligo_suffix",
            "oh5_mask_start",
            "oh3_mask_start",
        ]
    ]
    parts_out.to_csv(output_dir / "parts.csv", index=False)
    # Keep an extended sidecar for debugging / redesign
    parts.to_csv(output_dir / "parts_full.csv", index=False)
    junction_9s.to_csv(output_dir / "junction_map.csv", index=False)
    candidates.to_csv(output_dir / "overhang_candidates.csv", index=False)
    target_map.to_csv(output_dir / "target_map.csv", index=False)

    # Also write a README for the profile
    readme = output_dir.parent / "README.md"
    if output_dir.name == "input":
        readme = output_dir / "GRASP_SOURCE.md"
    else:
        readme = output_dir / "README.md"
    readme.write_text(
        """# GRASP sequences (Farley et al., Nucleic Acids Research 2025)

Imported from the deposited GenBank modules (`GRASP_-1.gb` / Assembly Planner).

- Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212
- Data: https://github.com/farleykvdg/GRASP

The CSV junction coordinates record the deposited design. Runtime redesign is
not tied to those indices: every synonymous four-base window fully inside the
invariant ARELF motif (offsets 0–11) is eligible.

Order fragments use configurable BsaI Level -1 entry interfaces, configurable
BpiI Level 0 block interfaces, and configurable final cassette interfaces. The
deposited GRASP preset remains available as pAGM1311 -> pAGM9121.

Deposited-preset order-strand geometry:
TTTGGTCTCAACAT{{pAGM1311 insert}}TTGTTGAGACCAAA

Deposited PPR block chains (not complete expression constructs):
- 9S: AGGT–CDS1–CTTC–CDS2–TTCG
- 14S: AGGT–CDS1–GTGA–CDS14–CTTC–CDS2–TTCG
- 19S: AGGT–CDS1–GTGA–CDS14–CACG–CDS19–CTTC–CDS2–TTCG
"""
    )

    return {
        "parts": parts_out,
        "parts_full": parts,
        "junction_map": junction_9s,
        "overhang_candidates": candidates,
        "target_map": target_map,
        "records": records,
    }


def compile_target_gap(
    target_rna: str,
    *,
    architecture: str = "9S",
    nterm_overhang: str = "AGGT",
) -> pd.DataFrame:
    """Compile a binding-site RNA into ordered GRASP part_ids (GAP algorithm).

    The returned assembly groups are five-part Level 0 blocks.  A 9S target
    produces the original ``CDS1``/``CDS2`` pair; 14S and 19S insert the
    corresponding ``CDS14`` and ``CDS19`` intermediate blocks described by the
    deposited GRASP assembly planner.
    """
    if not isinstance(architecture, str):
        raise ValueError("Architecture must be one of: 9S, 14S, 19S")
    architecture = architecture.strip().upper()
    if architecture not in ARCHITECTURE_LAYOUTS:
        raise ValueError(
            f"Unsupported architecture {architecture!r}; expected one of: "
            f"{', '.join(ARCHITECTURE_LAYOUTS)}"
        )

    if not isinstance(target_rna, str):
        raise ValueError("Binding target must be a DNA/RNA string")
    rna = target_rna.strip().upper().replace("T", "U")
    invalid = sorted(set(rna) - set("ACGU"))
    if invalid:
        raise ValueError(
            "Binding target contains noncanonical bases: " + ", ".join(invalid)
        )

    layout = ARCHITECTURE_LAYOUTS[architecture]
    expected_length = int(layout["target_length"])
    if len(rna) != expected_length:
        raise ValueError(
            f"Architecture {architecture} requires a {expected_length}-base "
            f"binding target; received {len(rna)} bases"
        )

    part_ids = pick_parts_for_target(rna, nterm_overhang=nterm_overhang)
    assembly_groups = tuple(layout["assembly_groups"])
    groups = [group for group in assembly_groups for _ in range(5)]
    orders = list(range(1, 6)) * len(assembly_groups)
    if len(part_ids) != len(groups):
        raise RuntimeError(
            f"Internal GRASP layout mismatch for {architecture}: selected "
            f"{len(part_ids)} parts for {len(groups)} assembly slots"
        )

    rows = []
    for i, part_id in enumerate(part_ids):
        rows.append(
            {
                "target_rna": rna,
                "architecture": architecture,
                "assembly_slot": i + 1,
                "part_id": part_id,
                "optimized_part_id": f"{part_id}_v1",
                "assembly_group": groups[i],
                "assembly_order": orders[i],
            }
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    profile = root / "grasp_library_project" / "profiles" / "grasp_nar2025" / "genbank"
    out = root / "grasp_library_project" / "input"
    result = import_grasp_profile(profile, out)
    print(f"parts: {len(result['parts'])}")
    print(f"junctions: {result['junction_map']['junction'].nunique()}")
    print(f"candidates: {len(result['overhang_candidates'])}")
    print(f"target_map rows: {len(result['target_map'])}")
    print(compile_target_gap("UUACACGUG"))
