"""Synthesis vendor constraint profiles for the GRASP library designer.

Twist Bioscience limits are taken from public FAQ / tech notes (2023–2025):
https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow
https://www.twistbioscience.com/products/genes/complex-genes

Twist scores sequences with a proprietary ML model; the values below are
published hard rules + design guidelines we can optimize against. They are
not a reimplementation of Twist’s scorer.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List

# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------

SYNTHESIS_VENDORS: Dict[str, Dict[str, Any]] = {
    "Twist · Express / Low complexity": {
        "vendor": "Twist Bioscience",
        "product": "Clonal Genes / Gene Fragments (aim for Low / Express)",
        "url": "https://www.twistbioscience.com/products/genes/gene-synthesis",
        "notes": (
            "Strictest practical targets for Express eligibility. "
            "Twist hard rule: avoid homopolymers ≥14 bp; no CcdB. "
            "Global GC 25–65%; keep 50 bp windows moderate; minimize repeats ≥12 bp."
        ),
        "synthesis": {
            "global_gc_min": 0.25,
            "global_gc_max": 0.65,
            "window_size": 50,
            "window_gc_min": 0.20,
            "window_gc_max": 0.80,
            "max_homopolymer": 3,
            "repeat_k": 12,
            "max_repeat_count": 1,
            "min_oligo_length": 20,
            "max_oligo_length": 300,  # oligo pool upper bound
            "min_gene_length": 300,
            "max_gene_length": 5000,
        },
        "hard_rules": [
            "Homopolymers must be < 14 bp (Twist hard rule)",
            "Do not include CcdB toxin sequences",
        ],
    },
    "Twist · Standard gene guidelines": {
        "vendor": "Twist Bioscience",
        "product": "Gene Fragments / Clonal Genes (Standard)",
        "url": "https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow",
        "notes": (
            "Published design guidelines: GC 25–65%, minimize homopolymers "
            "(≤13 nt tip), avoid repeats ≥20 bp, keep 50 bp GC windows from "
            "differing by >~50–52 points. Final acceptance is ML-scored."
        ),
        "synthesis": {
            "global_gc_min": 0.25,
            "global_gc_max": 0.65,
            "window_size": 50,
            "window_gc_min": 0.15,
            "window_gc_max": 0.85,
            "max_homopolymer": 3,
            "repeat_k": 16,
            "max_repeat_count": 1,
            "min_oligo_length": 20,
            "max_oligo_length": 300,
            "min_gene_length": 300,
            "max_gene_length": 5000,
        },
        "hard_rules": [
            "Homopolymers must be < 14 bp (Twist hard rule)",
            "Do not include CcdB toxin sequences",
        ],
    },
    "Twist · Complex Genes tolerant": {
        "vendor": "Twist Bioscience",
        "product": "Complex Clonal Genes",
        "url": "https://www.twistbioscience.com/products/genes/complex-genes",
        "notes": (
            "Upper envelope Twist can manufacture as High complexity: "
            "global GC 25–75%, local GC 10–90% (50 bp), homopolymers up to 30 bp, "
            "long repeats up to 200 bp. Expect longer TAT / higher price."
        ),
        "synthesis": {
            "global_gc_min": 0.25,
            "global_gc_max": 0.75,
            "window_size": 50,
            "window_gc_min": 0.10,
            "window_gc_max": 0.90,
            "max_homopolymer": 3,
            "repeat_k": 20,
            "max_repeat_count": 2,
            "min_oligo_length": 20,
            "max_oligo_length": 300,
            "min_gene_length": 300,
            "max_gene_length": 7000,
        },
        "hard_rules": [
            "Still avoid features outside Complex Genes envelope",
            "Do not include CcdB toxin sequences",
        ],
    },
    "IDT · gBlocks / eBlocks conservative": {
        "vendor": "IDT",
        "product": "gBlocks / eBlocks (conservative heuristic)",
        "url": "https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments",
        "notes": (
            "Conservative heuristic commonly used for short dsDNA fragments: "
            "moderate GC, short homopolymers, low repeat density. "
            "Verify against current IDT online checker before ordering."
        ),
        "synthesis": {
            "global_gc_min": 0.25,
            "global_gc_max": 0.70,
            "window_size": 40,
            "window_gc_min": 0.20,
            "window_gc_max": 0.80,
            "max_homopolymer": 3,
            "repeat_k": 12,
            "max_repeat_count": 2,
            "min_oligo_length": 40,
            "max_oligo_length": 500,
            "min_gene_length": 125,
            "max_gene_length": 3000,
        },
        "hard_rules": [],
    },
    "Generic · conservative (default)": {
        "vendor": "Generic",
        "product": "Vendor-agnostic conservative defaults",
        "url": "",
        "notes": "Balanced defaults used when no vendor is selected.",
        "synthesis": {
            "global_gc_min": 0.30,
            "global_gc_max": 0.65,
            "window_size": 40,
            "window_gc_min": 0.25,
            "window_gc_max": 0.75,
            "max_homopolymer": 3,
            "repeat_k": 12,
            "max_repeat_count": 2,
            "min_oligo_length": 40,
            "max_oligo_length": 500,
            "min_gene_length": 40,
            "max_gene_length": 5000,
        },
        "hard_rules": [],
    },
    "Custom (keep current CONFIG)": {
        "vendor": "Custom",
        "product": "User-defined",
        "url": "",
        "notes": "Does not overwrite synthesis parameters in CONFIG.",
        "synthesis": None,
        "hard_rules": [],
    },
}


ASSEMBLY_ENZYMES: Dict[str, Dict[str, str]] = {
    "BsaI (GGTCTC)": {"BsaI": "GGTCTC"},
    "BpiI / BbsI (GAAGAC)": {"BpiI": "GAAGAC"},
    "BsmBI / Esp3I (CGTCTC)": {"BsmBI": "CGTCTC"},
    "SapI / BspQI (GCTCTTC)": {"SapI": "GCTCTTC"},
    "GRASP default · BsaI + BpiI + BsmBI": {
        "BsaI": "GGTCTC",
        "BpiI": "GAAGAC",
        "BsmBI": "CGTCTC",
    },
    "None (no enzyme filter)": {},
}


LIGATION_TABLES: Dict[str, Dict[str, Any]] = {
    "T4 · 18 h · 25 °C (Potapov)": {"temperature": 25, "hours": 18, "table": None},
    "T4 · 18 h · 37 °C (Potapov)": {"temperature": 37, "hours": 18, "table": None},
    "T4 · 1 h · 25 °C (Potapov)": {"temperature": 25, "hours": 1, "table": None},
    "BsaI-HFv2 · constant 37 °C": {
        "temperature": 37,
        "hours": 18,
        "table": "BsaI-HFv2_T4_constant_37.csv",
    },
    "BsmBI-v2 · constant 42 °C": {
        "temperature": 42,
        "hours": 18,
        "table": "BsmBI-v2_T4_constant_42.csv",
    },
    "SapI": {
        "temperature": 25,
        "hours": 18,
        "table": "SapI.csv",
    },
}


def vendor_names() -> List[str]:
    return list(SYNTHESIS_VENDORS.keys())


def enzyme_names() -> List[str]:
    return list(ASSEMBLY_ENZYMES.keys())


def ligation_table_names() -> List[str]:
    return list(LIGATION_TABLES.keys())


def apply_vendor_to_config(config: Dict[str, Any], vendor_name: str) -> Dict[str, Any]:
    """Return a copy of config with synthesis block updated from the vendor profile."""
    updated = deepcopy(config)
    profile = SYNTHESIS_VENDORS[vendor_name]
    updated["synthesis_vendor"] = vendor_name
    updated["synthesis_vendor_meta"] = {
        "vendor": profile["vendor"],
        "product": profile["product"],
        "url": profile["url"],
        "notes": profile["notes"],
        "hard_rules": profile["hard_rules"],
    }
    if profile["synthesis"] is not None:
        # Keep unknown keys; overwrite known synthesis constraints
        synth = dict(updated.get("synthesis", {}))
        synth.update(profile["synthesis"])
        updated["synthesis"] = synth
    return updated


def apply_enzyme_to_config(config: Dict[str, Any], enzyme_name: str) -> Dict[str, Any]:
    updated = deepcopy(config)
    updated["assembly_enzyme"] = enzyme_name
    updated["forbidden_sites"] = dict(ASSEMBLY_ENZYMES[enzyme_name])
    return updated


def apply_ligation_table_to_config(
    config: Dict[str, Any], table_name: str
) -> Dict[str, Any]:
    updated = deepcopy(config)
    meta = LIGATION_TABLES[table_name]
    lig = dict(updated.get("ligation", {}))
    lig["temperature"] = meta["temperature"]
    lig["hours"] = meta["hours"]
    lig["table_name"] = table_name
    lig["ligation_table"] = meta["table"]
    updated["ligation"] = lig
    return updated


def twist_length_advice(length_bp: int) -> str:
    """Human-readable Twist product fit for a fragment length."""
    if length_bp < 20:
        return "Too short for Twist oligo pools (<20 bp)."
    if length_bp <= 300:
        return (
            f"{length_bp} bp fits Twist oligo pools (20–300 bp). "
            "Gene Fragments / Clonal Genes need ≥300 bp (pad or order as oligos)."
        )
    if length_bp <= 5000:
        return f"{length_bp} bp fits Twist Gene Fragments and Clonal Genes (0.3–5 kb)."
    if length_bp <= 7000:
        return f"{length_bp} bp may fit Complex Clonal Genes (up to ~7 kb); expect High complexity."
    return f"{length_bp} bp exceeds typical Twist clonal gene length (~7 kb)."
