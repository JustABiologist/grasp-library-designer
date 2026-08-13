"""Helpers to apply Colab Form / plain-dict settings into CONFIG."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Union

from .assembly_interfaces import (
    CANONICAL_NOTATION,
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    deposited_grasp_interface_preset,
)
from .codon_tables import apply_organism_codon_table, load_codon_usage
from .restriction_sites import DEFAULT_SITE_BLACKLIST, apply_site_blacklist_to_config
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
    architecture: str = "9S",
    synthesis_vendor: str,
    assembly_enzyme: str,
    ligation_table: str,
    site_blacklist: str = ",".join(DEFAULT_SITE_BLACKLIST),
    overhang_redesign: bool = False,
    redesign_plasmid_overhangs: bool = False,
    redesign_level0_junctions: Optional[bool] = None,
    redesign_selection: str = "knee",
    assembly_interface_preset: str = "auto",
    level_minus1_5prime_overhang: str = "ACAT",
    level_minus1_3prime_overhang: str = "ACAA",
    level0_5prime_overhang: str = "CTCA",
    level0_3prime_overhang: str = "CTCG",
    level1_5prime_overhang: str = "GGAG",
    level1_3prime_overhang: str = "AGCG",
    level_minus1_vector: Optional[str] = None,
    level0_vector: Optional[str] = None,
    level1_vector: Optional[str] = None,
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
    cfg["ppr_5prime_fusion_site"] = deposited_grasp_interface_preset()["level0"][
        "ppr_outer"
    ][FIVE_PRIME_CODING_SITE]
    cfg["architecture"] = architecture

    from .dna import reverse_complement

    def _overhang(value: str, label: str) -> str:
        cleaned = str(value).strip().upper().replace("U", "T")
        if len(cleaned) != 4 or set(cleaned) - set("ACGT"):
            raise ValueError(f"{label} must be exactly four DNA bases (ACGT)")
        return cleaned

    level_minus1_5prime = _overhang(
        level_minus1_5prime_overhang, "Level −1 5′ overhang"
    )
    level_minus1_3prime = _overhang(
        level_minus1_3prime_overhang, "Level −1 3′ overhang"
    )
    level0_5prime = _overhang(level0_5prime_overhang, "Level 0 5′ overhang")
    level0_3prime = _overhang(level0_3prime_overhang, "Level 0 3′ overhang")
    level1_5prime = _overhang(level1_5prime_overhang, "Level 1 5′ overhang")
    level1_3prime = _overhang(level1_3prime_overhang, "Level 1 3′ overhang")
    toolbox = deposited_grasp_interface_preset()
    entry_defaults = toolbox["level_minus1_entry"]
    level0_defaults = toolbox["level0"]["acceptor_outer"]
    level1_defaults = toolbox["final_cassette"]
    ppr_5prime = toolbox["level0"]["ppr_outer"][FIVE_PRIME_CODING_SITE]
    standard_values = (
        level_minus1_5prime == entry_defaults[FIVE_PRIME_END]
        and level_minus1_3prime == entry_defaults[THREE_PRIME_END]
        and level0_5prime == level0_defaults[FIVE_PRIME_END]
        and level0_3prime == level0_defaults[THREE_PRIME_END]
        and level1_5prime == level1_defaults[FIVE_PRIME_END]
        and level1_3prime == level1_defaults[THREE_PRIME_END]
    )
    preset_request = str(assembly_interface_preset).strip().lower()
    if preset_request == "auto":
        selected_profile = "deposited_grasp" if standard_values else "custom"
    elif preset_request == "deposited_grasp":
        if not standard_values:
            raise ValueError(
                "The deposited GRASP preset requires the deposited 5′/3′ overhangs"
            )
        selected_profile = "deposited_grasp"
    elif preset_request == "custom":
        selected_profile = "custom"
    else:
        raise ValueError(
            "assembly_interface_preset must be 'auto', 'custom', or 'deposited_grasp'"
        )

    level_minus1_vector = level_minus1_vector or (
        entry_defaults["vector_id"]
        if (level_minus1_5prime, level_minus1_3prime)
        == (
            entry_defaults[FIVE_PRIME_END],
            entry_defaults[THREE_PRIME_END],
        )
        else "custom_level_minus1_entry"
    )
    level0_vector = level0_vector or (
        toolbox["level0"]["acceptor_id"]
        if (level0_5prime, level0_3prime)
        == (
            level0_defaults[FIVE_PRIME_END],
            level0_defaults[THREE_PRIME_END],
        )
        else "custom_level0_acceptor"
    )
    level1_vector = level1_vector or (
        level1_defaults["vector_id"]
        if (level1_5prime, level1_3prime)
        == (
            level1_defaults[FIVE_PRIME_END],
            level1_defaults[THREE_PRIME_END],
        )
        else "custom_level1_acceptor"
    )
    offsets = {
        "terminal_to_cds2": 11,
        "cds1_to_cds14": 4,
        "cds14_to_cds19": 1,
    }
    junctions = {
        name: {
            "upstream_three_prime_end_overhang": values[
                "upstream_three_prime_end_overhang"
            ],
            "downstream_five_prime_end_overhang": values[
                "downstream_five_prime_end_overhang"
            ],
            "assembled_coding_site": values["assembled_coding_site"],
            "arelf_offset_nt": offsets[name],
        }
        for name, values in toolbox["junctions"].items()
    }
    cfg["assembly_interfaces"] = {
        "preset": selected_profile,
        "terminal_side_convention": "construct_ends_5prime_and_3prime",
        "overhang_sequence_notation": "5prime_to_3prime",
        "notation": CANONICAL_NOTATION,
        "coding_strand_direction": "5prime_to_3prime",
        "level_minus1_entry": {
            "profile": selected_profile,
            "vector_name": str(level_minus1_vector).strip(),
            "enzyme": "BsaI",
            FIVE_PRIME_END: level_minus1_5prime,
            THREE_PRIME_END: level_minus1_3prime,
            FIVE_PRIME_CODING_SITE: level_minus1_5prime,
            THREE_PRIME_CODING_SITE: reverse_complement(level_minus1_3prime),
        },
        "level0": {
            "acceptor_name": str(level0_vector).strip(),
            "release_enzyme": "BpiI / BbsI",
            "acceptor_outer": {
                FIVE_PRIME_END: level0_5prime,
                THREE_PRIME_END: level0_3prime,
                FIVE_PRIME_CODING_SITE: level0_5prime,
                THREE_PRIME_CODING_SITE: reverse_complement(level0_3prime),
            },
        },
        "junctions": junctions,
        "level1": {
            "acceptor_name": str(level1_vector).strip(),
            FIVE_PRIME_END: level1_5prime,
            THREE_PRIME_END: level1_3prime,
            FIVE_PRIME_CODING_SITE: level1_5prime,
            THREE_PRIME_CODING_SITE: reverse_complement(level1_3prime),
        },
    }
    level0_on = (
        bool(overhang_redesign)
        if redesign_level0_junctions is None
        else bool(redesign_level0_junctions)
    )
    cfg["overhang_redesign"] = {
        "plasmid_overhangs": bool(redesign_plasmid_overhangs),
        "level0_junctions": level0_on,
        "enabled": level0_on,
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
    cfg = apply_site_blacklist_to_config(cfg, site_blacklist)
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
