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
from typing import Any, Dict, List, Mapping

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
        "machine_hard_constraints": {"max_homopolymer": 13},
        "manual_review_rules": ["Do not include CcdB toxin sequences"],
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
        "machine_hard_constraints": {"max_homopolymer": 13},
        "manual_review_rules": ["Do not include CcdB toxin sequences"],
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
        "machine_hard_constraints": {"max_homopolymer": 30},
        "manual_review_rules": [
            "Review features outside the Complex Genes envelope",
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
        "machine_hard_constraints": {},
        "manual_review_rules": [],
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
        "machine_hard_constraints": {},
        "manual_review_rules": [],
    },
    "Custom (keep current CONFIG)": {
        "vendor": "Custom",
        "product": "User-defined",
        "url": "",
        "notes": "Does not overwrite synthesis parameters in CONFIG.",
        "synthesis": None,
        "hard_rules": [],
        "machine_hard_constraints": {},
        "manual_review_rules": [],
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


_POTAPOV_DATA_URL = (
    "https://acs.figshare.com/articles/dataset/"
    "Comprehensive_Profiling_of_Four_Base_Overhang_Ligation_Fidelity_by_"
    "T4_DNA_Ligase_and_Application_to_DNA_Assembly/7267505"
)
_PRYOR_ARTICLE_URL = (
    "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592"
)


def _potapov_static(temperature: int, hours: int, table: str | None) -> Dict[str, Any]:
    return {
        "temperature": temperature,
        "hours": hours,
        "table": table,
        "assay_kind": "static_ligation",
        "overhang_length": 4,
        "ligase": "T4 DNA Ligase",
        "restriction_enzyme": None,
        "buffer": "1x T4 DNA Ligase Reaction Buffer",
        "cycles": None,
        "steps": [{"temperature_c": temperature, "hours": hours}],
        "terminal_step": None,
        "source_doi": "10.1021/acssynbio.8b00333",
        "source_url": _POTAPOV_DATA_URL,
        "grasp_status": "surrogate",
        "notes": (
            "T4-ligase-only endpoint matrix. Potapov et al. showed that the "
            "25 C/18 h dataset approximately predicts a typical 37/16 C "
            "Golden Gate cycling protocol; it is not a measured 16 C matrix."
        ),
    }


def _pryor_cycling(
    *, table: str, restriction_enzyme: str, supplement_url: str, proxy_for: str
) -> Dict[str, Any]:
    return {
        # Static temperature/time fields are deliberately unset: this is a
        # whole Golden Gate cycling assay, not a constant-temperature matrix.
        "temperature": None,
        "hours": None,
        "table": table,
        "assay_kind": "golden_gate_cycling",
        "overhang_length": 4,
        "ligase": "T4 DNA Ligase",
        "restriction_enzyme": restriction_enzyme,
        "buffer": "see Pryor et al. 2020 Methods for enzyme-specific formulation",
        "cycles": 30,
        "steps": [
            {"temperature_c": 37, "minutes": 5},
            {"temperature_c": 16, "minutes": 5},
        ],
        "terminal_step": {"temperature_c": 60, "minutes": 5},
        "source_doi": "10.1371/journal.pone.0238592",
        "source_url": supplement_url,
        "article_url": _PRYOR_ARTICLE_URL,
        "grasp_status": "proxy",
        "proxy_for": proxy_for,
        "grasp_reference_protocol": {
            "initial_step": {"temperature_c": 37, "seconds": 20},
            "cycles": 26,
            "steps": [
                {"temperature_c": 37, "minutes": 3},
                {"temperature_c": 16, "minutes": 4},
            ],
            "terminal_steps": [
                {"temperature_c": 50, "minutes": 5},
                {"temperature_c": 80, "minutes": 5},
            ],
            "source_url": (
                "https://academic.oup.com/nar/article/53/20/gkaf1169/8321212"
            ),
        },
        "notes": (
            "Measured endpoint matrix for the complete 30-cycle reaction. "
            "GRASP used 26 cycles of 3 min at 37 C and 4 min at 16 C with "
            "different enzyme amounts/formulation, so this is a proxy rather "
            "than an exact prediction of the GRASP experiment."
        ),
    }


# Pryor et al. 2020 measured whole Golden Gate 37↔16 °C cycling matrices.
# BbsI-HF is the published isoschizomer proxy for Thermo Fisher BpiI (GAAGAC).
LEVEL_MINUS1_LIGATION = (
    "GRASP Level −1 proxy · BsaI-HFv2 + T4 · 37↔16 °C cycling (Pryor 2020)"
)
LEVEL0_LIGATION = (
    "GRASP Level 0 proxy · BbsI-HF + T4 · 37↔16 °C cycling (Pryor 2020)"
)
LEVEL1_LIGATION = (
    "GRASP Level 1 proxy · BsaI-HFv2 + T4 · 37↔16 °C cycling (Pryor 2020)"
)

_BSAI_HFV2 = _pryor_cycling(
    table="BsaI-HFv2.csv",
    restriction_enzyme="BsaI-HFv2",
    supplement_url=(
        "https://journals.plos.org/plosone/article/file?id="
        "10.1371/journal.pone.0238592.s001&type=supplementary"
    ),
    proxy_for="GRASP BsaI-HFv2 + T4 DNA Ligase (Levels −1 and 1)",
)
_BBSI_HF = _pryor_cycling(
    table="BbsI-HF.csv",
    restriction_enzyme="BbsI-HF",
    supplement_url=(
        "https://journals.plos.org/plosone/article/file?id="
        "10.1371/journal.pone.0238592.s004&type=supplementary"
    ),
    proxy_for="GRASP Level 0 Thermo Fisher BpiI + T4 DNA Ligase",
)

LIGATION_TABLES: Dict[str, Dict[str, Any]] = {
    LEVEL_MINUS1_LIGATION: {
        **deepcopy(_BSAI_HFV2),
        "proxy_for": "GRASP Level −1 BsaI-HFv2 + T4 DNA Ligase",
        "cloning_level": "level_minus1",
    },
    LEVEL0_LIGATION: {
        **deepcopy(_BBSI_HF),
        "cloning_level": "level0",
    },
    LEVEL1_LIGATION: {
        **deepcopy(_BSAI_HFV2),
        "proxy_for": "GRASP Level 1 BsaI-HFv2 + T4 DNA Ligase",
        "cloning_level": "level1",
    },
    "T4 ligase only · 18 h · 25 °C (Potapov 2018; validated cycling proxy)": (
        _potapov_static(25, 18, None)
    ),
    "T4 ligase only · 1 h · 25 °C (Potapov 2018)": _potapov_static(
        25, 1, None
    ),
    "T4 ligase only · 1 h · 37 °C (Potapov 2018)": _potapov_static(
        37, 1, "FileS_T4_01h_37C.csv"
    ),
    "T4 ligase only · 18 h · 37 °C (Potapov 2018)": _potapov_static(
        37, 18, None
    ),
}

# Enzyme-matched defaults for each physical GRASP cloning stage.
GRASP_LIGATION_BY_LEVEL: Dict[str, str] = {
    "level_minus1": LEVEL_MINUS1_LIGATION,
    "level0": LEVEL0_LIGATION,
    "level1": LEVEL1_LIGATION,
}


def vendor_names() -> List[str]:
    return list(SYNTHESIS_VENDORS.keys())


def enzyme_names() -> List[str]:
    return list(ASSEMBLY_ENZYMES.keys())


def ligation_table_names() -> List[str]:
    return list(LIGATION_TABLES.keys())


def redesign_ligation_table_names() -> List[str]:
    """Matrices offered for Level 0 overhang-redesign scoring.

    Levels −1 and 1 always use the BsaI-HFv2 cycling proxy; they are not
    selectable as the Level 0 redesign objective.
    """
    return [
        name
        for name, meta in LIGATION_TABLES.items()
        if meta.get("cloning_level") in {None, "level0"}
    ]


def _protocol_payload(table_name: str) -> Dict[str, Any]:
    meta = LIGATION_TABLES[table_name]
    return {
        "table_name": table_name,
        "temperature": meta["temperature"],
        "hours": meta["hours"],
        "ligation_table": meta["table"],
        "protocol_metadata": {
            key: deepcopy(value)
            for key, value in meta.items()
            if key not in {"temperature", "hours", "table"}
        },
    }


def grasp_ligation_by_level(
    *,
    level0_override: str | None = None,
) -> Dict[str, Dict[str, Any]]:
    """Return enzyme-matched Pryor proxies for each GRASP cloning stage."""
    mapping = dict(GRASP_LIGATION_BY_LEVEL)
    if level0_override is not None:
        if level0_override not in LIGATION_TABLES:
            raise ValueError(f"Unknown ligation table: {level0_override!r}")
        mapping["level0"] = level0_override
    return {
        level: _protocol_payload(table_name)
        for level, table_name in mapping.items()
    }


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
        "machine_hard_constraints": profile["machine_hard_constraints"],
        "manual_review_rules": profile["manual_review_rules"],
    }
    if profile["synthesis"] is not None:
        # Keep unknown keys; overwrite known synthesis constraints
        synth = dict(updated.get("synthesis", {}))
        synth.update(profile["synthesis"])
        updated["synthesis"] = synth
    return updated


def apply_enzyme_to_config(config: Dict[str, Any], enzyme_name: str) -> Dict[str, Any]:
    """Select internal-site domestication rules without changing cloning flanks.

    Existing GRASP modules and vectors have fixed Type IIS architecture. This
    selector therefore controls only which internal recognition sites are
    forbidden during sequence optimization; it must not imply that orderable
    flanks or destination vectors have been converted to a different enzyme.
    """
    updated = deepcopy(config)
    if enzyme_name not in ASSEMBLY_ENZYMES:
        raise ValueError(f"Unknown domestication enzyme filter: {enzyme_name!r}")
    updated["assembly_enzyme"] = enzyme_name
    updated["domestication_enzyme_filter"] = enzyme_name
    updated["assembly_enzyme_semantics"] = "internal_site_filter_only"
    updated["assembly_flanks_modified"] = False
    updated["forbidden_sites"] = dict(ASSEMBLY_ENZYMES[enzyme_name])
    return updated


def apply_ligation_table_to_config(
    config: Dict[str, Any], table_name: str
) -> Dict[str, Any]:
    """Apply a Level 0 redesign matrix and freeze enzyme-matched stage proxies.

    Selecting a Potapov/T4-only table overrides only the Level 0 redesign
    objective. Level −1 and Level 1 remain BsaI-HFv2 37↔16 °C cycling proxies;
    Level 0 remains BbsI-HF unless the selected table itself is that proxy or
    another Level 0 override.
    """
    if table_name not in LIGATION_TABLES:
        raise ValueError(f"Unknown ligation table: {table_name!r}")
    updated = deepcopy(config)
    meta = LIGATION_TABLES[table_name]
    level0_name = (
        table_name
        if meta.get("cloning_level") in {None, "level0"}
        else LEVEL0_LIGATION
    )
    # Selecting the Level −1 / Level 1 BsaI proxy must not silently retarget
    # Level 0 redesign scoring away from BbsI-HF / BpiI.
    if meta.get("cloning_level") in {"level_minus1", "level1"}:
        level0_name = LEVEL0_LIGATION
        primary_name = LEVEL0_LIGATION
    else:
        primary_name = table_name

    by_level = grasp_ligation_by_level(level0_override=level0_name)
    primary = by_level["level0"] if primary_name == level0_name else _protocol_payload(
        primary_name
    )
    # Primary ligation block stays the Level 0 redesign matrix for backwards
    # compatibility with workflows that read config["ligation"].
    lig = dict(updated.get("ligation", {}))
    lig.update(primary)
    lig["by_level"] = by_level
    lig["redesign_level"] = "level0"
    updated["ligation"] = lig
    return updated


def ligation_protocol_for_level(
    config: Dict[str, Any], level: str = "level0"
) -> Dict[str, Any]:
    """Return the protocol block for one cloning stage."""
    level = str(level).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "minus1": "level_minus1",
        "level_1": "level1",
        "l0": "level0",
        "l1": "level1",
        "lm1": "level_minus1",
        "entry": "level_minus1",
    }
    level = aliases.get(level, level)
    lig = dict(config.get("ligation", {}))
    by_level = lig.get("by_level")
    if not isinstance(by_level, Mapping) or level not in by_level:
        by_level = grasp_ligation_by_level(
            level0_override=lig.get("table_name")
            if lig.get("table_name") in redesign_ligation_table_names()
            else None
        )
    if level not in by_level:
        raise ValueError(
            f"Unknown cloning level {level!r}; expected one of "
            f"{sorted(GRASP_LIGATION_BY_LEVEL)}"
        )
    return deepcopy(by_level[level])


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
