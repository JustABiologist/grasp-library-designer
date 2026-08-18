"""Annotated GenBank export for Geneious / SnapGene.

Library path: ``write_annotated_genbank`` — pAGM1311 order fragments.
One-shot path: ``write_oneshot_genbank`` — assembled CDS, ligated insert, oligos.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

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
from .binder import (
    FIVE_PRIME_SOLVATING_HELIX,
    REPEAT_TEMPLATE,
    assembly_role_from_part_id,
    describe_part_id,
)
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


def _ppr_repeat_geometry() -> tuple[int, int, int]:
    marked = REPEAT_TEMPLATE.format(fifth="!", last=".")
    return marked.index("!"), marked.index("."), len(marked)


def _sort_features(features: list[SeqFeature]) -> list[SeqFeature]:
    return sorted(
        features,
        key=lambda feature: (
            int(feature.location.start),
            -(int(feature.location.end) - int(feature.location.start)),
        ),
    )


def _wrap_layout(sequence: str, geometry: Mapping[str, Any]) -> dict[str, Any]:
    clamp5 = clean_dna(str(geometry["clamp_5p"]))
    site = clean_dna(str(geometry["recognition_site"]))
    spacer = int(geometry["spacer_len"])
    clamp3 = clean_dna(str(geometry["clamp_3p"]))
    enzyme = str(geometry.get("enzyme") or "TypeIIS")
    prefix = clamp5 + site + ("A" * spacer)
    suffix = ("T" * spacer) + reverse_complement(site) + clamp3
    if not sequence.startswith(prefix) or not sequence.endswith(suffix):
        raise ValueError("one-shot oligo does not match wrap geometry")
    length = len(sequence)
    payload_start = len(prefix)
    payload_end = length - len(suffix)
    rec5_start = len(clamp5)
    rec5_end = rec5_start + len(site)
    rec3_start = payload_end + spacer
    rec3_end = rec3_start + len(site)
    return {
        "enzyme": enzyme,
        "prefix": prefix,
        "suffix": suffix,
        "clamp5": (0, len(clamp5)),
        "rec5": (rec5_start, rec5_end),
        "cut5": (payload_start - 1, payload_start),
        "payload": (payload_start, payload_end),
        "oh5": (payload_start, payload_start + 4),
        "oh3": (payload_end - 4, payload_end),
        "cut3": (payload_end, payload_end + 1),
        "rec3": (rec3_start, rec3_end),
        "clamp3": (rec3_end, length),
        "primer_f": (0, len(prefix)),
        "primer_r": (payload_end, length),
    }


def _oneshot_orf_features(
    cds: str,
    binder: Mapping[str, Any],
    *,
    offset: int = 0,
    genetic_code: int = 1,
) -> list[SeqFeature]:
    cds = clean_dna(cds)
    aa = str(binder.get("aa_sequence") or "")
    rna = str(binder.get("target_rna") or "")
    features: list[SeqFeature] = []
    _add(
        features,
        _feature(
            offset,
            offset + len(cds),
            "CDS",
            f"GRASP {len(rna)}S binder CDS" if rna else "binder CDS",
            gene=f"GRASP_{rna}" if rna else "GRASP_binder",
            product=f"GRASP {len(rna)}S PPR binder" if rna else "GRASP PPR binder",
            codon_start="1",
            transl_table=str(genetic_code),
            translation=aa,
            note=(
                "Continuous one-shot ORF: ATG-Met, native solvating helix, "
                "then one PPR-like repeat per target RNA base. Destination "
                "sticky ends are cloning adapters outside this feature."
            ),
        ),
    )
    aa_pos = 0
    if binder.get("include_start_codon", True) and cds.startswith("ATG"):
        _add(
            features,
            _feature(
                offset,
                offset + 3,
                "misc_feature",
                "start codon ATG",
                note="Initiating methionine. Not a destination overhang.",
            ),
        )
        aa_pos = 1
    helix_nt0 = offset + 3 * aa_pos
    helix_nt1 = helix_nt0 + 3 * len(FIVE_PRIME_SOLVATING_HELIX)
    _add(
        features,
        _feature(
            helix_nt0,
            helix_nt1,
            "misc_feature",
            "solvating helix",
            note=f"Native 5′ GRASP helix {FIVE_PRIME_SOLVATING_HELIX}.",
        ),
    )
    aa_pos += len(FIVE_PRIME_SOLVATING_HELIX)
    fifth_i, last_i, repeat_len = _ppr_repeat_geometry()
    pairs = list(binder.get("ppr_pairs") or [])
    for i, base in enumerate(rna):
        pair = pairs[i] if i < len(pairs) else ""
        nt0 = offset + 3 * aa_pos
        nt1 = nt0 + 3 * repeat_len
        fifth = pair[0] if pair else "?"
        last = pair[1] if len(pair) > 1 else "?"
        _add(
            features,
            _feature(
                nt0,
                nt1,
                "misc_feature",
                f"PPR{i + 1} {base} ({pair})",
                note=(
                    f"31-aa PPR-like repeat {i + 1} targeting {base}; "
                    f"5th/last code {fifth}/{last}."
                ),
            ),
        )
        _add(
            features,
            _feature(
                nt0 + 3 * fifth_i,
                nt0 + 3 * fifth_i + 3,
                "misc_feature",
                f"PPR{i + 1} 5th {fifth}",
                note=f"PPR 5th residue encoding target base {base}.",
            ),
        )
        _add(
            features,
            _feature(
                nt0 + 3 * last_i,
                nt0 + 3 * last_i + 3,
                "misc_feature",
                f"PPR{i + 1} last {last}",
                note=f"PPR last residue encoding target base {base}.",
            ),
        )
        aa_pos += repeat_len
    if aa and aa_pos != len(aa):
        raise AssertionError(
            f"one-shot architecture features covered {aa_pos} aa; protein is {len(aa)} aa"
        )
    return features


def _oneshot_junction_features(
    cds: str,
    cuts: Sequence[int],
    overhangs: Sequence[str],
    *,
    offset: int = 0,
) -> list[SeqFeature]:
    features: list[SeqFeature] = []
    cds = clean_dna(cds)
    for i, (cut, overhang) in enumerate(zip(cuts, overhangs), start=1):
        start = offset + 3 * int(cut) - 4
        end = offset + 3 * int(cut)
        oh = clean_dna(overhang)
        if cds[start - offset : end - offset] != oh:
            raise AssertionError(
                f"junction J{i} {oh} does not match CDS[{start - offset}:{end - offset}]"
            )
        _add(
            features,
            _feature(
                start,
                end,
                "misc_feature",
                f"junction overhang J{i} {oh}",
                note=f"Codon-aligned Golden Gate junction after residue {cut}.",
            ),
        )
    return features


def _oneshot_oligo_features(
    sequence: str,
    geometry: Mapping[str, Any],
    *,
    fragment_id: str,
    oh5: str,
    oh3: str,
    cds_slice: str,
    primers: Optional[Mapping[str, Any]],
    is_first: bool,
    is_last: bool,
) -> list[SeqFeature]:
    layout = _wrap_layout(sequence, geometry)
    enzyme = str(layout["enzyme"])
    payload_start, payload_end = layout["payload"]
    cds_slice = clean_dna(cds_slice)
    cds_start = payload_start + 4
    cds_end = cds_start + len(cds_slice)
    if sequence[cds_start:cds_end] != cds_slice:
        raise AssertionError(f"{fragment_id}: fragment CDS is not after the 5′ 4-nt overhang")
    primers = primers or {}
    forward = str(primers.get("forward") or layout["prefix"])
    reverse = str(primers.get("reverse") or reverse_complement(layout["suffix"]))
    tm_f = primers.get("tm_forward_c")
    tm_r = primers.get("tm_reverse_c")
    oh5_label = (
        f"destination overhang 5' {oh5}" if is_first else f"junction overhang 5' {oh5}"
    )
    oh3_label = (
        f"destination overhang 3' {oh3}" if is_last else f"junction overhang 3' {oh3}"
    )
    features: list[SeqFeature] = []
    _add(
        features,
        _feature(
            *layout["primer_f"],
            "primer_bind",
            "pool PCR forward",
            note=(
                f"Shared 5′ wrap arm / pool primer {forward}"
                + (f" (Tm {tm_f} °C)." if tm_f is not None else ".")
            ),
        ),
    )
    _add(
        features,
        _feature(
            *layout["clamp5"],
            "misc_feature",
            f"PCR adapter 5' {sequence[layout['clamp5'][0]:layout['clamp5'][1]]}",
            note="5′ clamp of the Type IIS wrap (part of the forward primer).",
        ),
    )
    _add(
        features,
        _feature(
            *layout["rec5"],
            "misc_feature",
            f"{enzyme} recognition 5'",
            note="Type IIS recognition on the order oligo; cut is downstream.",
        ),
    )
    _add(
        features,
        _feature(
            *layout["cut5"],
            "misc_feature",
            f"{enzyme} 5' cut",
            note="Plus-strand cut after the spacer N; 4-base overhang follows.",
        ),
    )
    _add(
        features,
        _feature(
            payload_start,
            payload_end,
            "misc_feature",
            f"{fragment_id} payload",
            note="DNA released by Type IIS digestion (overhangs + fragment CDS).",
        ),
    )
    _add(
        features,
        _feature(
            *layout["oh5"],
            "misc_feature",
            oh5_label,
            note=(
                f"Coding-strand 5′ sticky end {oh5}."
                + (
                    " Destination adapter 5′ of the ORF."
                    if is_first
                    else " Shared with the previous fragment 3′ overhang."
                )
            ),
        ),
    )
    _add(
        features,
        _feature(
            cds_start,
            cds_end,
            "misc_feature",
            f"{fragment_id} CDS slice",
            note="Unique coding DNA on this fragment (5′ 4-nt overlap excluded).",
        ),
    )
    _add(
        features,
        _feature(
            *layout["oh3"],
            "misc_feature",
            oh3_label,
            note=(
                f"Coding-strand 3′ site {oh3}; physical sticky end "
                f"{reverse_complement(oh3)}."
                + (
                    " Destination adapter 3′ of the ORF."
                    if is_last
                    else " Shared with the next fragment 5′ overhang."
                )
            ),
        ),
    )
    _add(
        features,
        _feature(
            *layout["cut3"],
            "misc_feature",
            f"{enzyme} 3' cut",
            note="Reverse-strand Type IIS cut; 4-base overhang is upstream on the plus strand.",
        ),
    )
    _add(
        features,
        _feature(
            *layout["rec3"],
            "misc_feature",
            f"{enzyme} recognition 3'",
            strand=-1,
            note=f"Reverse complement of {clean_dna(str(geometry['recognition_site']))} on the plus strand.",
        ),
    )
    _add(
        features,
        _feature(
            *layout["clamp3"],
            "misc_feature",
            f"PCR adapter 3' {sequence[layout['clamp3'][0]:layout['clamp3'][1]]}",
            note="3′ clamp of the Type IIS wrap (part of the reverse primer).",
        ),
    )
    _add(
        features,
        _feature(
            *layout["primer_r"],
            "primer_bind",
            "pool PCR reverse",
            strand=-1,
            note=(
                f"Binds the shared 3′ wrap arm; primer {reverse}"
                + (f" (Tm {tm_r} °C)." if tm_r is not None else ".")
            ),
        ),
    )
    return features


def _oneshot_record(
    sequence: str,
    *,
    rec_id: str,
    name: str,
    description: str,
    features: list[SeqFeature],
    comment: str = "",
) -> SeqRecord:
    annotations: dict[str, Any] = {
        "molecule_type": "DNA",
        "topology": "linear",
        "data_file_division": "SYN",
    }
    if comment:
        annotations["comment"] = comment
    record = SeqRecord(
        Seq(clean_dna(sequence)),
        id=rec_id,
        name=_sanitize_name(name),
        description=description,
        annotations=annotations,
    )
    record.features = _sort_features(features)
    return record


def oneshot_seqrecords(
    design: Any,
    binder: Mapping[str, Any],
    *,
    primers: Optional[Mapping[str, Any]] = None,
    genetic_code: int = 1,
) -> list[SeqRecord]:
    """Assembled CDS, ligated insert, and one annotated record per order oligo."""
    cds = clean_dna(design.cds)
    rna = str(binder.get("target_rna") or "")
    arch = f"{int(binder.get('n_bases') or len(rna) or 0)}S"
    dest5 = clean_dna(design.destination_5prime)
    dest3 = clean_dna(design.destination_3prime_coding)
    enzyme = str(design.wrap_enzyme)
    insert = dest5 + cds + dest3
    orf_features = _oneshot_orf_features(
        cds, binder, genetic_code=genetic_code
    )
    junction_features = _oneshot_junction_features(
        cds, design.cuts, design.overhangs
    )
    cds_comment = (
        f"GRASP one-shot {arch} binder ORF for {rna}. "
        f"Wrap enzyme {enzyme}. Destination adapters {dest5}/{dest3} sit outside "
        "this record; see the GG_insert record and each order oligo. "
        f"Internal junctions: {';'.join(design.overhangs) or 'none'}."
    )
    records = [
        _oneshot_record(
            cds,
            rec_id=f"GRASP_{rna}_CDS" if rna else "GRASP_CDS",
            name=f"G{arch}_CDS",
            description=f"Assembled GRASP {arch} binder CDS targeting {rna}",
            features=orf_features + junction_features,
            comment=cds_comment,
        )
    ]
    insert_features = _oneshot_orf_features(
        cds, binder, offset=len(dest5), genetic_code=genetic_code
    )
    insert_features.extend(
        _oneshot_junction_features(
            cds, design.cuts, design.overhangs, offset=len(dest5)
        )
    )
    _add(
        insert_features,
        _feature(
            0,
            len(dest5),
            "misc_feature",
            f"destination overhang 5' {dest5}",
            note="Cloning adapter 5′ of the ORF (not translated).",
        ),
    )
    _add(
        insert_features,
        _feature(
            len(dest5) + len(cds),
            len(insert),
            "misc_feature",
            f"destination overhang 3' {dest3}",
            note=(
                f"Coding-strand 3′ site {dest3}; physical sticky end "
                f"{reverse_complement(dest3)}."
            ),
        ),
    )
    records.append(
        _oneshot_record(
            insert,
            rec_id=f"GRASP_{rna}_insert" if rna else "GRASP_insert",
            name=f"G{arch}_GG",
            description=(
                f"Ligated Golden Gate insert ({dest5} + ORF + {dest3}) targeting {rna}"
            ),
            features=insert_features,
            comment=(
                "In-silico product after Type IIS digestion and ligation of the "
                f"order oligos. Plus-strand sticky-end sites are {dest5} and {dest3}."
            ),
        )
    )
    n = len(design.oligos)
    aa_bounds = [0, *list(design.cuts), len(design.aa_sequence)]
    for i, oligo in enumerate(design.oligos):
        frag_id = f"F{i + 1}"
        aa0, aa1 = int(aa_bounds[i]), int(aa_bounds[i + 1])
        cds_slice = cds[3 * aa0 : 3 * aa1]
        oh5 = dest5 if i == 0 else clean_dna(design.overhangs[i - 1])
        oh3 = dest3 if i == n - 1 else clean_dna(design.overhangs[i])
        features = _oneshot_oligo_features(
            clean_dna(oligo),
            design.geometry,
            fragment_id=frag_id,
            oh5=oh5,
            oh3=oh3,
            cds_slice=cds_slice,
            primers=primers,
            is_first=i == 0,
            is_last=i == n - 1,
        )
        records.append(
            _oneshot_record(
                oligo,
                rec_id=f"{rna}_{frag_id}" if rna else frag_id,
                name=f"G{arch}_{frag_id}",
                description=(
                    f"Order oligo {frag_id} ({enzyme}) {oh5}..{oh3} "
                    f"for GRASP {arch} {rna}"
                ),
                features=features,
                comment=(
                    f"Wrapped Type IIS oligo {frag_id}. Pool-amplify with the "
                    "shared PCR primers annotated as primer_bind features, digest "
                    f"with {enzyme}, and ligate on the labelled 4-nt overhangs."
                ),
            )
        )
    return records


def write_oneshot_genbank(
    path: Path | str,
    design: Any,
    binder: Mapping[str, Any],
    *,
    primers: Optional[Mapping[str, Any]] = None,
    genetic_code: int = 1,
) -> Path:
    """Write a multi-record .gb of the assembled CDS, ligated insert, and oligos.

    Library pAGM1311 export is ``write_annotated_genbank``; this is one-shot only.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = oneshot_seqrecords(
        design,
        binder,
        primers=primers,
        genetic_code=genetic_code,
    )
    if not records:
        raise ValueError("no one-shot sequences to write as GenBank")
    SeqIO.write(records, path, "genbank")
    return path.resolve()
