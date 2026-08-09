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
    overhang_redesign: bool = True,
    redesign_selection: str = "knee",
    optimize_depth: int = 2000,
    n_fragments: Optional[int] = None,
    kazusa_species_id: str = "",
    upload_bytes: Optional[Union[bytes, bytearray, memoryview, str]] = None,
    upload_filename: Optional[str] = None,
    prompt_upload_if_needed: bool = True,
) -> Dict[str, Any]:
    """Apply Colab Form values; return dict with updated config + codon_data."""
    cfg = dict(config)
    cfg["genetic_code"] = int(genetic_code)
    cfg["target_rna"] = str(target_rna).strip().upper().replace("T", "U")
    cfg["nterm_overhang"] = nterm_overhang
    cfg["architecture"] = architecture
    cfg["overhang_redesign"] = {
        "enabled": bool(overhang_redesign),
        "selection": redesign_selection,
    }
    cfg.setdefault("optimizer", {})
    cfg["optimizer"] = dict(cfg["optimizer"])
    cfg["optimizer"]["iterations_per_part"] = int(optimize_depth)
    if n_fragments is not None and int(n_fragments) > 0:
        cfg["oneshot_n_fragments"] = int(n_fragments)
    else:
        cfg.pop("oneshot_n_fragments", None)

    cfg = apply_vendor_to_config(cfg, synthesis_vendor)
    cfg = apply_enzyme_to_config(cfg, assembly_enzyme)
    cfg = apply_ligation_table_to_config(cfg, ligation_table)

    codon_path = Path(cfg["codon_usage_file"])
    codon_path.parent.mkdir(parents=True, exist_ok=True)
    cfg["selected_organism"] = organism
    if kazusa_species_id:
        cfg["kazusa_species_id"] = str(kazusa_species_id).strip()

    resolved_upload = upload_bytes
    resolved_name = upload_filename
    if organism == UPLOAD_OWN_TABLE and resolved_upload is None and prompt_upload_if_needed:
        try:
            from .codon_upload import prompt_colab_codon_upload

            _, uploaded_meta = prompt_colab_codon_upload(
                codon_path, genetic_code=int(genetic_code)
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
            genetic_code=int(genetic_code),
            kazusa_species_id=kazusa_species_id or cfg.get("kazusa_species_id"),
            upload_bytes=resolved_upload,
            upload_filename=resolved_name,
        )
        cfg["selected_organism_label"] = meta.get("organism", organism)
        if meta.get("species_id"):
            cfg["kazusa_species_id"] = meta["species_id"]
    else:
        table, codon_data = load_codon_usage(codon_path, genetic_code=int(genetic_code))
        meta = {"organism": organism}

    return {
        "codon_table": table,
        "codon_data": codon_data,
        "config": cfg,
        "meta": meta,
    }
