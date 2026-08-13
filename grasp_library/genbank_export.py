"""Annotated GenBank export of GRASP order fragments for Geneious / SnapGene."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .assembly_interfaces import (
    FIVE_PRIME_CODING_SITE,
    THREE_PRIME_CODING_SITE,
    order_fragment_arms,
    resolve_assembly_interfaces,
    reverse_complement,
)
from .binder import assembly_role_from_part_id, describe_part_id
from .dna import clean_dna, translate_dna


def _qual(**items: str) -> dict[str, list[str]]:
    return {key: [value] for key, value in items.items() if value}


def _feature(
    start: int,
    end: int,
    ftype: str,
    label: str,
    *,
    strand: int = 1,
    **qualifiers: str,
) -> SeqFeature:
    if end <= start:
        raise ValueError(f"empty feature {label!r}: {start}:{end}")
    quals = _qual(label=label, **qualifiers)
    return SeqFeature(
        SimpleLocation(start, end, strand=strand),
        type=ftype,
        qualifiers=quals,
    )


def _row_get(row: Any, name: str, default: Any = None) -> Any:
    if isinstance(row, Mapping):
        return row.get(name, default)
    if hasattr(row, name):
        value = getattr(row, name)
        if value is not None and not (isinstance(value, float) and pd.isna(value)):
            return value
    return default


def _order_sequence(row: Any) -> str:
    for name in ("order_sequence_5to3", "oligo_sequence_5to3"):
        value = _row_get(row, name)
        if value:
            return clean_dna(value)
    raise ValueError("library row has no order/oligo DNA sequence")


def _sanitize_name(name: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_]", "_", str(name))[:16]
    return token or "GRASP_part"


def _codon_start_for_insert(oh5: int) -> int:
    leftover = int(oh5) % 3
    if leftover == 0:
        return 1
    if leftover == 1:
        return 3
    return 2


def _in_frame_translation(cds: str, oh5: int, oh3: int, genetic_code: int) -> str:
    start = ((int(oh5) + 2) // 3) * 3
    end = ((int(oh3) + 4) // 3) * 3
    window = cds[start:end]
    if len(window) < 3:
        return ""
    return translate_dna(window, genetic_code=genetic_code).rstrip("*")


def _add(features: list[SeqFeature], feature: SeqFeature | None) -> None:
    if feature is not None:
        features.append(feature)


def fragment_features(
    row: Any,
    *,
    config: Optional[Mapping] = None,
) -> tuple[str, list[SeqFeature], dict[str, Any]]:
    """Return order DNA, Geneious features, and coordinate metadata."""
    sequence = _order_sequence(row)
    part_id = str(_row_get(row, "part_id", "part"))
    cds = clean_dna(_row_get(row, "optimized_cds", "") or "")
    oh5 = int(_row_get(row, "oh5_mask_start", 0) or 0)
    oh3 = int(_row_get(row, "oh3_mask_start", max(len(cds) - 4, 0)) or 0)
    genetic_code = int(_row_get(row, "genetic_code", 1) or 1)
    role = str(_row_get(row, "assembly_role") or assembly_role_from_part_id(part_id))
    info = describe_part_id(part_id)
    interfaces = resolve_assembly_interfaces(config)
    prefix, suffix = order_fragment_arms(interfaces)
    order = interfaces["order_fragment"]
    entry = interfaces["level_minus1_entry"]
    outer = interfaces["level0"].get("acceptor_outer") or {}

    if not sequence.startswith(prefix) or not sequence.endswith(suffix):
        raise ValueError(f"{part_id}: order DNA does not match configured BsaI arms")

    payload_start = len(prefix)
    payload_end = len(sequence) - len(suffix)
    payload = sequence[payload_start:payload_end]
    insert = cds[oh5 : oh3 + 4] if cds else payload
    insert_in_payload = payload.find(insert) if insert else 0
    if insert_in_payload < 0:
        insert = payload
        insert_in_payload = 0
    insert_start = payload_start + insert_in_payload
    insert_end = insert_start + len(insert)

    rec5_start = len(order["clamp_5p"])
    rec5_end = rec5_start + len(order["recognition_site"])
    entry5_start = payload_start - 4
    entry5_end = payload_start
    entry3_start = payload_end
    entry3_end = payload_end + 4
    rec3_start = entry3_end + len(order["spacer_3p"])
    rec3_end = rec3_start + len(order["recognition_site"])

    oh5_seq = clean_dna(_row_get(row, "oh5", insert[:4] if insert else ""))
    oh3_seq = clean_dna(_row_get(row, "oh3", insert[-4:] if len(insert) >= 4 else ""))
    five_end = clean_dna(_row_get(row, "five_prime_end_overhang", oh5_seq))
    three_end = clean_dna(
        _row_get(row, "three_prime_end_overhang", reverse_complement(oh3_seq) if oh3_seq else "")
    )
    oh5_junc = str(_row_get(row, "oh5_junction", "") or "")
    oh3_junc = str(_row_get(row, "oh3_junction", "") or "")

    features: list[SeqFeature] = []
    _add(
        features,
        _feature(
            rec5_start,
            rec5_end,
            "misc_feature",
            f"{order['enzyme']} recognition 5'",
            note="Type IIS recognition on the order fragment; cut is downstream.",
        ),
    )
    _add(
        features,
        _feature(
            rec5_end,
            rec5_end + 1,
            "misc_feature",
            f"{order['enzyme']} 5' cut",
            note="Cut after GGTCTC-N; 4-base overhang follows.",
        ),
    )
    _add(
        features,
        _feature(
            entry5_start,
            entry5_end,
            "misc_feature",
            f"Level -1 overhang 5' {entry[FIVE_PRIME_CODING_SITE]}",
            note=(
                f"Physical 5' sticky end {entry.get('five_prime_end_overhang', '')} "
                f"(coding site {entry[FIVE_PRIME_CODING_SITE]})."
            ),
        ),
    )
    _add(
        features,
        _feature(
            payload_start,
            payload_end,
            "misc_feature",
            "BpiI-released payload",
            note="Module released from the Level -1 entry vector by BpiI / BbsI.",
        ),
    )

    if role == "A" and insert_in_payload >= 4:
        l0_5 = payload[:4]
        expected = str(outer.get(FIVE_PRIME_CODING_SITE, l0_5))
        _add(
            features,
            _feature(
                payload_start,
                payload_start + 4,
                "misc_feature",
                f"Level 0 overhang 5' {l0_5}",
                note=f"pAGM9121 outer 5' coding site (expected {expected}). A-parts only.",
            ),
        )
    if role == "E" and len(payload) >= len(insert) + 4:
        l0_3 = payload[-4:]
        expected = str(outer.get(THREE_PRIME_CODING_SITE, l0_3))
        _add(
            features,
            _feature(
                payload_end - 4,
                payload_end,
                "misc_feature",
                f"Level 0 overhang 3' {l0_3}",
                note=(
                    f"pAGM9121 outer 3' coding site (expected {expected}). "
                    f"Physical 3' sticky end is {reverse_complement(l0_3)}."
                ),
            ),
        )

    if len(insert) >= 4:
        _add(
            features,
            _feature(
                insert_start,
                insert_start + 4,
                "misc_feature",
                f"module overhang 5' {oh5_seq}",
                note=(
                    f"Coding-strand site {oh5_seq}"
                    + (f" ({oh5_junc})" if oh5_junc else "")
                    + (f"; physical 5' end {five_end}." if five_end else ".")
                ),
            ),
        )
        _add(
            features,
            _feature(
                insert_end - 4,
                insert_end,
                "misc_feature",
                f"module overhang 3' {oh3_seq}",
                note=(
                    f"Coding-strand site {oh3_seq}"
                    + (f" ({oh3_junc})" if oh3_junc else "")
                    + (
                        f"; physical 3' sticky end {three_end}."
                        if three_end
                        else "."
                    )
                ),
            ),
        )

    if len(insert) >= 3:
        codon_start = _codon_start_for_insert(oh5)
        translation = _in_frame_translation(cds, oh5, oh3, genetic_code)
        cds_note = (
            "Overhang-bounded coding insert on the order fragment. "
            "Golden Gate modules may start mid-codon; codon_start is set so "
            "the translation covers complete in-frame residues only."
        )
        if info.get("target_rna_base"):
            cds_note += f" PPR target RNA base: {info['target_rna_base']}."
        _add(
            features,
            _feature(
                insert_start,
                insert_end,
                "CDS",
                f"{part_id} CDS",
                gene=part_id,
                product=f"GRASP {role or 'module'} PPR fragment",
                codon_start=str(codon_start),
                transl_table=str(genetic_code),
                translation=translation,
                note=cds_note,
            ),
        )

    _add(
        features,
        _feature(
            entry3_start,
            entry3_end,
            "misc_feature",
            f"Level -1 overhang 3' {entry[THREE_PRIME_CODING_SITE]}",
            note=(
                f"Coding site {entry[THREE_PRIME_CODING_SITE]}; "
                f"physical 3' sticky end {entry.get('three_prime_end_overhang', '')}."
            ),
        ),
    )
    _add(
        features,
        _feature(
            rec3_start - 1,
            rec3_start,
            "misc_feature",
            f"{order['enzyme']} 3' cut",
            note="Reverse-strand Type IIS cut; 4-base overhang is upstream on the plus strand.",
        ),
    )
    _add(
        features,
        _feature(
            rec3_start,
            rec3_end,
            "misc_feature",
            f"{order['enzyme']} recognition 3'",
            strand=-1,
            note="Reverse complement of GGTCTC on the order-fragment plus strand.",
        ),
    )

    meta = {
        "part_id": part_id,
        "role": role,
        "payload_start": payload_start,
        "payload_end": payload_end,
        "insert_start": insert_start,
        "insert_end": insert_end,
        "oh5": oh5_seq,
        "oh3": oh3_seq,
    }
    return sequence, features, meta


def library_seqrecords(
    library: pd.DataFrame,
    *,
    config: Optional[Mapping] = None,
) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    used_names: set[str] = set()
    for row in library.itertuples(index=False):
        sequence, features, meta = fragment_features(row, config=config)
        part_id = meta["part_id"]
        rec_id = str(_row_get(row, "optimized_part_id", part_id))
        name = _sanitize_name(rec_id)
        if name in used_names:
            name = _sanitize_name(f"{name}{len(used_names)}")
        used_names.add(name)
        info = describe_part_id(part_id)
        description = (
            f"GRASP order fragment {part_id} slot {info['assembly_role'] or '?'} "
            f"joins {info['joins_upstream_role'] or 'L-1'}"
            f"→{info['assembly_role'] or '?'}"
            f"→{info['joins_downstream_role'] or 'L-1'}"
        )
        if info["target_rna_base"]:
            description += f"; PPR {info['ppr_5th_aa'] or '-'}5/{info['ppr_last_aa'] or '-'}last targets {info['target_rna_base']}"
        record = SeqRecord(
            Seq(sequence),
            id=rec_id,
            name=name,
            description=description,
            annotations={
                "molecule_type": "DNA",
                "topology": "linear",
                "data_file_division": "SYN",
            },
        )
        record.features = features
        records.append(record)
    return records


def write_annotated_genbank(
    library: pd.DataFrame,
    path: Path | str,
    *,
    config: Optional[Mapping] = None,
) -> Path:
    """Write one multi-record .gb file of annotated order fragments."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = library_seqrecords(library, config=config)
    if not records:
        raise ValueError("no order fragments to write as GenBank")
    SeqIO.write(records, path, "genbank")
    return path.resolve()
