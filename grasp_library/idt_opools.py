"""IDT oPools DNA list prices (EUR).

Tiers match the IDT oPools quote table (one tube, cumulative bases).
This is a list-price estimate: no VAT, shipping, or currency conversion.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Sequence

POOL_FLOOR_EUR = 109.00
POOL_FLOOR_BASES = 3300
TIER2_END = 50_000
TIER3_END = 100_000
TIER2_EUR_PER_BASE = 0.038
TIER3_EUR_PER_BASE = 0.025
TIER4_EUR_PER_BASE = 0.013
PHOSPHO_EUR_PER_OLIGO = 1.63
MIN_OLIGO_NT = 40
MAX_OLIGO_NT = 350

SCALES: List[Dict[str, Any]] = [
    {"name": "1 pmol", "min_oligos": 100, "max_oligos": 20_000},
    {"name": "10 pmol", "min_oligos": 10, "max_oligos": 2_000},
    {"name": "50 pmol", "min_oligos": 2, "max_oligos": 384},
]


def _round_eur(value: float) -> float:
    return round(float(value) + 1e-12, 2)


def _dna_eur(total_bases: int) -> tuple[float, Dict[str, int]]:
    remaining = max(int(total_bases), 0)
    b1 = min(remaining, POOL_FLOOR_BASES)
    remaining -= b1
    b2 = min(remaining, TIER2_END - POOL_FLOOR_BASES)
    remaining -= b2
    b3 = min(remaining, TIER3_END - TIER2_END)
    remaining -= b3
    b4 = max(0, remaining)
    dna = 0.0
    if total_bases > 0:
        dna += POOL_FLOOR_EUR
        dna += b2 * TIER2_EUR_PER_BASE
        dna += b3 * TIER3_EUR_PER_BASE
        dna += b4 * TIER4_EUR_PER_BASE
    return _round_eur(dna), {
        "tier1_bases": b1,
        "tier2_bases": b2,
        "tier3_bases": b3,
        "tier4_bases": b4,
    }


def eligible_opool_scales(n_oligos: int) -> List[str]:
    n = int(n_oligos)
    return [
        scale["name"]
        for scale in SCALES
        if scale["min_oligos"] <= n <= scale["max_oligos"]
    ]


def price_idt_opool(
    oligo_lengths: Sequence[int],
    *,
    phosphorylate_5prime: bool = False,
) -> Dict[str, Any]:
    """Price one IDT oPools tube from wrapped oligo lengths (nt)."""
    lengths = [int(n) for n in oligo_lengths]
    n = len(lengths)
    total_bases = int(sum(lengths))
    dna_eur, tiers = _dna_eur(total_bases)
    phospho_eur = _round_eur(PHOSPHO_EUR_PER_OLIGO * n) if n else 0.0
    warnings: List[str] = []
    if n and (min(lengths) < MIN_OLIGO_NT or max(lengths) > MAX_OLIGO_NT):
        warnings.append(
            f"Oligo length must be {MIN_OLIGO_NT}–{MAX_OLIGO_NT} nt for oPools"
        )
    scales = eligible_opool_scales(n)
    if n and not scales:
        warnings.append(
            f"{n} oligos is outside every listed oPools scale "
            "(1 pmol needs ≥100; 10 pmol needs 10–2000; 50 pmol needs 2–384)"
        )
    return {
        "vendor": "IDT",
        "product": "oPools DNA",
        "currency": "EUR",
        "n_oligos": n,
        "total_bases": total_bases,
        "min_oligo_nt": min(lengths) if lengths else 0,
        "max_oligo_nt": max(lengths) if lengths else 0,
        "dna_eur": dna_eur,
        "phospho_eur": phospho_eur,
        "total_eur": _round_eur(dna_eur + (phospho_eur if phosphorylate_5prime else 0.0)),
        "phosphorylate_5prime": bool(phosphorylate_5prime),
        "eligible_scales": scales,
        "warnings": warnings,
        "list_price_not_a_quote": True,
        **tiers,
    }


def quote_from_oligo_table(
    oligos: Mapping[str, Any],
    *,
    phosphorylate_5prime: bool = False,
) -> Dict[str, Any]:
    lengths = list(oligos["oligo_length"])
    return price_idt_opool(lengths, phosphorylate_5prime=phosphorylate_5prime)


def format_idt_opool_status(quote: Mapping[str, Any]) -> str:
    scales = ", ".join(quote.get("eligible_scales") or ()) or "no listed scale"
    return (
        f"IDT oPools (one tube, list EUR): <b>€{quote['dna_eur']:.0f}</b> "
        f"(~€110 at the pool floor) · "
        f"{quote['n_oligos']} oligos · {quote['total_bases']:,} bases · "
        f"5′ P would add €{quote['phospho_eur']:.2f} · scale {scales}"
    )
