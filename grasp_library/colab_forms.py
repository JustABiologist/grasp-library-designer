"""Helpers to apply Colab Form / plain-dict settings into CONFIG."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .codon_tables import apply_organism_codon_table, load_codon_usage
from .sample_codon_tables import SAMPLE_CODON_TABLES
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
    if (
        organism in SAMPLE_CODON_TABLES
        and SAMPLE_CODON_TABLES[organism]["frequencies"] is not None
    ):
        apply_organism_codon_table(
            organism, codon_path, genetic_code=int(genetic_code)
        )
    table, codon_data = load_codon_usage(codon_path, genetic_code=int(genetic_code))
    return {"codon_table": table, "codon_data": codon_data, "config": cfg}
