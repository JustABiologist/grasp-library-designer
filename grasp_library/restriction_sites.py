"""Curated restriction-site blacklist used during sequence domestication.

The optimizer already treats ``config["forbidden_sites"]`` as a hard filter.
This module supplies the extra user-facing blacklist: ~100 common cloning
enzymes, with SapI, BsaI, and BpiI selected by default. The assembly-enzyme
filter is merged on top so GRASP Type IIS sites stay forbidden unless the
user turns that filter off.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple, Union

# Recognition sites from REBASE / Biopython Restriction. Four-cutters that
# appear constantly in CDS (AluI, HaeIII, MboI, …) are omitted on purpose:
# blacklisting them makes synonymous design fail more often than it helps.
COMMON_RESTRICTION_SITES: Dict[str, str] = {
    "AarI": "CACCTGC",
    "AatII": "GACGTC",
    "Acc65I": "GGTACC",
    "AccI": "GTMKAC",
    "AflII": "CTTAAG",
    "AgeI": "ACCGGT",
    "AleI": "CACNNNNGTG",
    "ApaI": "GGGCCC",
    "ApaLI": "GTGCAC",
    "AscI": "GGCGCGCC",
    "AseI": "ATTAAT",
    "AsiSI": "GCGATCGC",
    "AvrII": "CCTAGG",
    "BamHI": "GGATCC",
    "BanI": "GGYRCC",
    "BbsI": "GAAGAC",
    "BclI": "TGATCA",
    "BcoDI": "GTCTC",
    "BfuAI": "ACCTGC",
    "BglI": "GCCNNNNNGGC",
    "BglII": "AGATCT",
    "BlpI": "GCTNAGC",
    "BpiI": "GAAGAC",
    "BsaI": "GGTCTC",
    "BsiWI": "CGTACG",
    "BsmAI": "GTCTC",
    "BsmBI": "CGTCTC",
    "BsmFI": "GGGAC",
    "BspEI": "TCCGGA",
    "BspHI": "TCATGA",
    "BspMI": "ACCTGC",
    "BspQI": "GCTCTTC",
    "BsrGI": "TGTACA",
    "BssHII": "GCGCGC",
    "BstBI": "TTCGAA",
    "BstEII": "GGTNACC",
    "BstXI": "CCANNNNNNTGG",
    "BstZ17I": "GTATAC",
    "BtgZI": "GCGATG",
    "ClaI": "ATCGAT",
    "DraI": "TTTAAA",
    "EagI": "CGGCCG",
    "EarI": "CTCTTC",
    "Eco53kI": "GAGCTC",
    "EcoNI": "CCTNNNNNAGG",
    "EcoRI": "GAATTC",
    "EcoRV": "GATATC",
    "Esp3I": "CGTCTC",
    "FokI": "GGATG",
    "FseI": "GGCCGGCC",
    "FspI": "TGCGCA",
    "HincII": "GTYRAC",
    "HindIII": "AAGCTT",
    "HpaI": "GTTAAC",
    "KasI": "GGCGCC",
    "KpnI": "GGTACC",
    "LguI": "GCTCTTC",
    "MfeI": "CAATTG",
    "MluI": "ACGCGT",
    "NaeI": "GCCGGC",
    "NarI": "GGCGCC",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NgoMIV": "GCCGGC",
    "NheI": "GCTAGC",
    "NotI": "GCGGCCGC",
    "NsiI": "ATGCAT",
    "PacI": "TTAATTAA",
    "PaqCI": "CACCTGC",
    "PflMI": "CCANNNNNTGG",
    "PmeI": "GTTTAAAC",
    "PsiI": "TTATAA",
    "PspOMI": "GGGCCC",
    "PstI": "CTGCAG",
    "PvuI": "CGATCG",
    "PvuII": "CAGCTG",
    "RsrII": "CGGWCCG",
    "SacI": "GAGCTC",
    "SacII": "CCGCGG",
    "SalI": "GTCGAC",
    "SapI": "GCTCTTC",
    "SbfI": "CCTGCAGG",
    "ScaI": "AGTACT",
    "SexAI": "ACCWGGT",
    "SfiI": "GGCCNNNNNGGCC",
    "SfoI": "GGCGCC",
    "SgrAI": "CRCCGGYG",
    "SmaI": "CCCGGG",
    "SnaBI": "TACGTA",
    "SpeI": "ACTAGT",
    "SphI": "GCATGC",
    "SrfI": "GCCCGGGC",
    "StuI": "AGGCCT",
    "SwaI": "ATTTAAAT",
    "Tth111I": "GACNNNGTC",
    "XbaI": "TCTAGA",
    "XhoI": "CTCGAG",
    "XmaI": "CCCGGG",
    "XmnI": "GAANNNNTTC",
    "ZraI": "GACGTC",
}

DEFAULT_SITE_BLACKLIST: Tuple[str, ...] = ("SapI", "BsaI", "BpiI")

_NAME_LOOKUP = {name.lower(): name for name in COMMON_RESTRICTION_SITES}

SiteBlacklistInput = Union[str, Sequence[str], None]


def restriction_site_names() -> List[str]:
    return list(COMMON_RESTRICTION_SITES.keys())


def restriction_site_options() -> List[Tuple[str, str]]:
    """Return ``(label, enzyme_name)`` pairs for a multi-select widget."""
    return [
        (f"{name} ({site})", name)
        for name, site in COMMON_RESTRICTION_SITES.items()
    ]


def parse_site_blacklist(value: SiteBlacklistInput) -> List[str]:
    """Normalize a form string or widget tuple into canonical enzyme names."""
    if value is None:
        return []
    if isinstance(value, str):
        parts = [item.strip() for item in value.replace(";", ",").split(",")]
        tokens = [item for item in parts if item]
    else:
        tokens = [str(item).strip() for item in value if str(item).strip()]

    resolved: List[str] = []
    seen = set()
    unknown: List[str] = []
    for token in tokens:
        name = _canonical_enzyme_name(token)
        if name is None:
            unknown.append(token)
            continue
        if name not in seen:
            seen.add(name)
            resolved.append(name)
    if unknown:
        sample = ", ".join(restriction_site_names()[:8])
        raise ValueError(
            "Unknown restriction enzyme(s) in site blacklist: "
            f"{', '.join(unknown)}. Use names such as {sample}, …"
        )
    return resolved


def resolve_site_blacklist(names: Iterable[str]) -> Dict[str, str]:
    return {
        name: COMMON_RESTRICTION_SITES[name]
        for name in parse_site_blacklist(list(names))
    }


def apply_site_blacklist_to_config(
    config: Mapping[str, Any],
    names: SiteBlacklistInput,
) -> Dict[str, Any]:
    """Merge the extra cut-site blacklist into ``forbidden_sites``.

    Assembly-enzyme sites already written by ``apply_enzyme_to_config`` are
    kept. Deselected blacklist enzymes are dropped so the optimizer only
    forbids the current union.
    """
    from .synthesis_vendors import ASSEMBLY_ENZYMES

    updated = deepcopy(dict(config))
    selected = parse_site_blacklist(names)
    updated["site_blacklist"] = selected
    assembly_key = updated.get("assembly_enzyme")
    assembly_sites = dict(ASSEMBLY_ENZYMES.get(assembly_key, {}))
    assembly_sites.update(resolve_site_blacklist(selected))
    updated["forbidden_sites"] = assembly_sites
    return updated


def _canonical_enzyme_name(token: str) -> str | None:
    cleaned = str(token).strip()
    if not cleaned:
        return None
    # Allow "SapI (GCTCTTC)" or "SapI / BspQI (GCTCTTC)" widget/form labels.
    head = cleaned.split("(", 1)[0]
    for part in head.split("/"):
        key = part.strip().lower()
        if key in _NAME_LOOKUP:
            return _NAME_LOOKUP[key]
    return _NAME_LOOKUP.get(cleaned.lower())
