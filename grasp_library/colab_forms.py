"""Helpers to apply Colab Form / plain-dict settings into CONFIG."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Union

from .codon_tables import apply_organism_codon_table, load_codon_usage
from .sample_codon_tables import SAMPLE_CODON_TABLES, UPLOAD_OWN_TABLE
from .synthesis_vendors import (
    apply_enzyme_to_config,
    apply_ligation_table_to_config,
    apply_vendor_to_config,
)


def apply_form_settings(
    config: Dict[str, Any],
    *,
    organism: str,
    genetic_code: int,
    target_rna: str,
    nterm_overhang: str = "AGGT",
    architecture: str = "9S",
    synthesis_vendor: str,
    assembly_enzyme: str,
    ligation_table: str,
    overhang_redesign: bool = False,
    redesign_selection: str = "knee",
    assembly_interface_preset: str = "custom",
    entry_vector_name: str = "Custom Level -1 entry vector",
    entry_n_overhang: str = "AACA",
    entry_c_overhang: str = "GGAG",
    level0_acceptor_name: str = "Custom Level 0 acceptor",
    level0_acceptor_n_overhang: str = "CTCA",
    level0_acceptor_c_overhang: str = "CGAG",
    cds1_c_overhang: str = "CTTC",
    cds2_n_overhang: str = "GAAG",
    cds1_cds2_arelf_offset_nt: int = 11,
    cds1_cds14_c_overhang: str = "GTGA",
    cds14_n_overhang: str = "TCAC",
    cds1_cds14_arelf_offset_nt: int = 4,
    cds14_cds19_c_overhang: str = "CACG",
    cds19_n_overhang: str = "CGTG",
    cds14_cds19_arelf_offset_nt: int = 1,
    level1_acceptor_name: str = "Custom Level 1 acceptor",
    level1_n_overhang: str = "GCCC",
    level1_c_overhang: str = "GCGA",
    optimize_depth: int = 2000,
    n_fragments: Optional[int] = None,
    kazusa_species_id: str = "",
    upload_bytes: Optional[Union[bytes, bytearray, memoryview, str]] = None,
    upload_filename: Optional[str] = None,
    prompt_upload_if_needed: bool = True,
) -> Dict[str, Any]:
    """Apply Colab Form values; return dict with updated config + codon_data.

    Built-in organism tables own their NCBI genetic code (e.g. Chlamy
    chloroplast → 11). The form ``genetic_code`` is only used for custom /
    upload / Kazusa-fetch modes. Translation QC always checks the loaded
    organism codon table.
    """
    cfg = dict(config)
    cfg["target_rna"] = str(target_rna).strip().upper().replace("T", "U")
    cfg["nterm_overhang"] = nterm_overhang
    cfg["architecture"] = architecture

    from .dna import reverse_complement

    def _overhang(value: str, label: str) -> str:
        cleaned = str(value).strip().upper().replace("U", "T")
        if len(cleaned) != 4 or set(cleaned) - set("ACGT"):
            raise ValueError(f"{label} must be exactly four DNA bases (ACGT)")
        return cleaned

    def _pair(upstream: str, downstream: str, label: str):
        up = _overhang(upstream, f"{label} upstream C overhang")
        down = _overhang(downstream, f"{label} downstream N overhang")
        if reverse_complement(up) != down:
            raise ValueError(
                f"{label}: directional terminal overhangs must be reverse "
                f"complements ({up} pairs with {reverse_complement(up)}, not {down})"
            )
        return up, down

    def _offset(value: int, label: str) -> int:
        result = int(value)
        if not 0 <= result <= 11:
            raise ValueError(f"{label} must be between 0 and 11 within ARELF")
        return result

    cds1_c, cds2_n = _pair(cds1_c_overhang, cds2_n_overhang, "CDS1→CDS2")
    cds1_14_c, cds14_n = _pair(
        cds1_cds14_c_overhang, cds14_n_overhang, "CDS1→CDS14"
    )
    cds14_19_c, cds19_n = _pair(
        cds14_cds19_c_overhang, cds19_n_overhang, "CDS14→CDS19"
    )
    cfg["assembly_interfaces"] = {
        "overhang_notation": "directional_terminal_5p",
        "level_minus1_entry": {
            "profile": "custom",
            "vector_name": str(entry_vector_name).strip(),
            "enzyme": "BsaI",
            "n_terminal_overhang": _overhang(
                entry_n_overhang, "Entry N overhang"
            ),
            "c_terminal_overhang": _overhang(
                entry_c_overhang, "Entry C overhang"
            ),
        },
        "level0": {
            "acceptor_name": str(level0_acceptor_name).strip(),
            "release_enzyme": "BpiI / BbsI",
            "acceptor_n_terminal_overhang": _overhang(
                level0_acceptor_n_overhang, "Level 0 acceptor N overhang"
            ),
            "acceptor_c_terminal_overhang": _overhang(
                level0_acceptor_c_overhang, "Level 0 acceptor C overhang"
            ),
            "block_junctions": {
                "cds1_to_cds14": {
                    "upstream_c": cds1_14_c,
                    "downstream_n": cds14_n,
                    "arelf_offset_nt": _offset(
                        cds1_cds14_arelf_offset_nt, "CDS1→CDS14 ARELF offset"
                    ),
                },
                "cds14_to_cds19": {
                    "upstream_c": cds14_19_c,
                    "downstream_n": cds19_n,
                    "arelf_offset_nt": _offset(
                        cds14_cds19_arelf_offset_nt, "CDS14→CDS19 ARELF offset"
                    ),
                },
                "cds1_to_cds2": {
                    "upstream_c": cds1_c,
                    "downstream_n": cds2_n,
                    "arelf_offset_nt": _offset(
                        cds1_cds2_arelf_offset_nt, "CDS1→CDS2 ARELF offset"
                    ),
                },
            },
        },
        "level1": {
            "acceptor_name": str(level1_acceptor_name).strip(),
            "n_terminal_overhang": _overhang(
                level1_n_overhang, "Level 1 cassette N overhang"
            ),
            "c_terminal_overhang": _overhang(
                level1_c_overhang, "Level 1 cassette C overhang"
            ),
        },
    }
    if str(assembly_interface_preset) == "deposited_grasp":
        if nterm_overhang != "AGGT":
            raise ValueError(
                "The deposited GRASP preset requires the AGGT PPR N interface"
            )
        cfg["assembly_interfaces"] = {"preset": "deposited_grasp"}
    elif str(assembly_interface_preset) != "custom":
        raise ValueError(
            "assembly_interface_preset must be 'custom' or 'deposited_grasp'"
        )
    cfg["overhang_redesign"] = {
        "enabled": bool(overhang_redesign),
        "selection": redesign_selection,
        "cut_mode": "movable_arelf",
        "allowed_arelf_offsets_nt": list(range(12)),
    }
    cfg.setdefault("optimizer", {})
    cfg["optimizer"] = dict(cfg["optimizer"])
    cfg["optimizer"]["iterations_per_part"] = int(optimize_depth)
    if n_fragments is not None:
        raise ValueError(
            "GRASP uses fixed five-part BpiI assemblies; arbitrary one-shot "
            "fragment counts are not supported"
        )
    cfg.pop("oneshot_n_fragments", None)

    cfg = apply_vendor_to_config(cfg, synthesis_vendor)
    cfg = apply_enzyme_to_config(cfg, assembly_enzyme)
    cfg = apply_ligation_table_to_config(cfg, ligation_table)

    codon_path = Path(cfg["codon_usage_file"])
    codon_path.parent.mkdir(parents=True, exist_ok=True)
    cfg["selected_organism"] = organism
    if kazusa_species_id:
        cfg["kazusa_species_id"] = str(kazusa_species_id).strip()

    entry = SAMPLE_CODON_TABLES.get(organism, {})
    # Built-in frequency tables: organism genetic code is authoritative.
    if entry.get("frequencies") is not None:
        code: Optional[int] = int(entry.get("genetic_code", 1))
    else:
        code = int(genetic_code)

    resolved_upload = upload_bytes
    resolved_name = upload_filename
    if organism == UPLOAD_OWN_TABLE and resolved_upload is None and prompt_upload_if_needed:
        try:
            from .codon_upload import prompt_colab_codon_upload

            _, uploaded_meta = prompt_colab_codon_upload(
                codon_path, genetic_code=int(code)
            )
            resolved_name = uploaded_meta.get("filename")
            # File already written; apply_organism_codon_table reuses the path
        except RuntimeError as exc:
            if codon_path.exists():
                print(
                    f"▶ {exc} Using existing {codon_path} instead. "
                    "For drag-and-drop upload locally, use GraspControlPanel."
                )
            else:
                raise

    if organism in SAMPLE_CODON_TABLES:
        table, codon_data, meta, _issues = apply_organism_codon_table(
            organism,
            codon_path,
            genetic_code=code,
            kazusa_species_id=kazusa_species_id or cfg.get("kazusa_species_id"),
            upload_bytes=resolved_upload,
            upload_filename=resolved_name,
        )
        cfg["selected_organism_label"] = meta.get("organism", organism)
        if meta.get("species_id"):
            cfg["kazusa_species_id"] = meta["species_id"]
        code = int(meta.get("genetic_code", code))
    else:
        table, codon_data = load_codon_usage(codon_path, genetic_code=int(code))
        meta = {"organism": organism, "genetic_code": int(code)}

    cfg["genetic_code"] = int(code)

    return {
        "codon_table": table,
        "codon_data": codon_data,
        "config": cfg,
        "meta": meta,
    }
