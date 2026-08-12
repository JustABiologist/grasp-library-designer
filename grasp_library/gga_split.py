"""Split an optimized CDS into Golden Gate fragments at chosen cut sites."""

from __future__ import annotations

import math
import random
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from .dna import clean_dna, is_self_reverse_complement, reverse_complement
from .ligation_fidelity import LigationFidelityCalculator


def suggest_fragment_count(
    cds_length: int,
    *,
    max_fragment_cds: int = 240,
    min_fragments: int = 2,
    max_fragments: int = 8,
) -> int:
    n = max(min_fragments, math.ceil(cds_length / max_fragment_cds))
    return int(min(max_fragments, n))


def _default_cut_aa_positions(n_aa: int, n_fragments: int) -> List[int]:
    """Exclusive AA ends of fragments 0..n-2 (codon-aligned cuts)."""
    if n_fragments < 2:
        return []
    if n_aa < n_fragments * 4:
        raise ValueError(f"Protein too short ({n_aa} aa) for {n_fragments} fragments")
    raw = [round(n_aa * (i + 1) / n_fragments) for i in range(n_fragments - 1)]
    cuts = []
    prev = 3
    for i, c in enumerate(raw):
        c = max(prev, min(int(c), n_aa - 3 * (n_fragments - 1 - i)))
        # keep overhang fully inside upstream codon block: cut after ≥2 aa from ends
        cuts.append(c)
        prev = c + 3
    return cuts


def overhangs_at_cuts(cds: str, cut_aa: Sequence[int]) -> List[str]:
    """4-nt overhangs = last 4 nt of each upstream fragment (codon-aligned cuts)."""
    cds = clean_dna(cds)
    ohs = []
    for aa_end in cut_aa:
        dna_end = 3 * int(aa_end)
        if dna_end < 4 or dna_end > len(cds):
            raise ValueError(f"Cut aa={aa_end} out of range for CDS length {len(cds)}")
        ohs.append(cds[dna_end - 4 : dna_end])
    return ohs


def _valid_overhang_set(overhangs: Sequence[str]) -> bool:
    used = set()
    for oh in overhangs:
        oh = clean_dna(oh)
        if len(oh) != 4 or is_self_reverse_complement(oh):
            return False
        rc = reverse_complement(oh)
        if oh in used or rc in used:
            return False
        used.add(oh)
        used.add(rc)
    return True


def search_cut_sites(
    cds: str,
    n_fragments: int,
    fidelity: LigationFidelityCalculator,
    *,
    seed: int = 42,
    window: int = 6,
    attempts: int = 800,
) -> Tuple[List[int], List[str], float]:
    """
    Search codon-aligned cut positions so the DNA overhangs at those cuts
    form a high-fidelity Golden Gate set.
    """
    cds = clean_dna(cds)
    n_aa = len(cds) // 3
    base = _default_cut_aa_positions(n_aa, n_fragments)
    if not base:
        return [], [], 1.0

    rng = random.Random(seed)
    best_cuts = list(base)
    best_ohs = overhangs_at_cuts(cds, best_cuts)
    try:
        best_f = (
            fidelity.set_fidelity(best_ohs)
            if _valid_overhang_set(best_ohs)
            else -1.0
        )
    except (KeyError, ValueError):
        best_f = -1.0

    n_junc = len(base)
    for _ in range(attempts):
        cuts = []
        prev = 3
        for i in range(n_junc):
            center = base[i]
            lo = max(prev, center - window)
            hi = min(n_aa - 3 * (n_junc - i), center + window)
            if lo > hi:
                lo, hi = prev, max(prev, n_aa - 3 * (n_junc - i))
            c = rng.randint(lo, hi) if hi >= lo else prev
            cuts.append(c)
            prev = c + 3
        ohs = overhangs_at_cuts(cds, cuts)
        if not _valid_overhang_set(ohs):
            continue
        try:
            f = fidelity.set_fidelity(ohs)
        except (KeyError, ValueError):
            continue
        if f > best_f:
            best_f = f
            best_cuts = cuts
            best_ohs = ohs
            if best_f >= 0.999:
                break

    if best_f < 0:
        raise RuntimeError(
            "Could not find cut sites yielding a valid overhang set. "
            "Try fewer fragments or a longer anneal."
        )
    return best_cuts, best_ohs, float(best_f)


def plan_fragments_from_cds(
    cds: str,
    aa_sequence: str,
    cut_aa: Sequence[int],
    overhangs: Sequence[str],
    *,
    ligation_fidelity: float,
) -> pd.DataFrame:
    """Build fragment rows (AA + coding_mask + DNA slice) from a full CDS."""
    cds = clean_dna(cds)
    aa = str(aa_sequence).upper().replace(" ", "")
    if 3 * len(aa) != len(cds):
        raise ValueError("CDS/AA length mismatch")

    bounds = [0] + list(cut_aa) + [len(aa)]
    n = len(bounds) - 1
    rows = []
    for i in range(n):
        aa0, aa1 = bounds[i], bounds[i + 1]
        dna0, dna1 = 3 * aa0, 3 * aa1
        # Each fragment includes its flanking junction overhangs
        if i > 0:
            dna0 = dna0 - 4  # start at shared 5′ overhang
        frag_dna = cds[dna0:dna1]
        frag_aa = aa[aa0:aa1]
        # Mask locks overhang bases only
        mask = ["N"] * len(frag_dna)
        oh5 = overhangs[i - 1] if i > 0 else ""
        oh3 = overhangs[i] if i < len(overhangs) else ""
        if oh5:
            mask[:4] = list(oh5)
        if oh3:
            mask[-4:] = list(oh3)

        # AA for the fragment DNA: if we prepended 4 nt overhang from previous
        # codon, those 4 nt are the end of the previous AA's codon — they should
        # already encode the start of this fragment's first AA when cut is
        # codon-aligned. When i>0, dna0 = 3*aa0 - 4, so frag_dna length =
        # 3*(aa1-aa0) + 4. That means DNA is NOT a clean ORF for frag_aa alone.
        #
        # Cleaner model: fragments are ORF slices without double-counting:
        #   frag i DNA = cds[3*aa0 : 3*aa1]  (codon complete)
        #   overhang for junction i is cds[3*cut-4:3*cut] which is the last 4 of upstream
        # For GGA synthesis, upstream oligo ends with overhang, downstream begins with it.
        # So downstream synthesizes from (cut-4) or we use Type IIS to expose overhangs
        # from longer flanks. Match GRASP deposited style: coding window includes overhangs.
        #
        # Use codon-complete fragments for optimization; attach overhangs only as
        # sticky-end design metadata, and build oligo CDS as:
        #   F1: cds[0:cut1]
        #   Fi: cds[cut_{i-1}-4 : cut_i]  so 5' overhang included
        #   Fn: cds[cut_{n-1}-4 : end]
        rows.append(
            {
                "fragment_id": f"F{i+1}",
                "assembly_order": i + 1,
                "aa_start_0based": aa0,
                "aa_end_0based": aa1,
                "aa_sequence": frag_aa,
                "cds_slice_start": dna0,
                "cds_slice_end": dna1,
                "fragment_cds": frag_dna,
                "coding_mask": "".join(mask),
                "oh5": oh5 or None,
                "oh3": oh3 or None,
                "ligation_fidelity_set": ligation_fidelity,
            }
        )
    return pd.DataFrame(rows)


def plan_gga_from_optimized_cds(
    cds: str,
    aa_sequence: str,
    *,
    n_fragments: Optional[int] = None,
    max_fragment_cds: int = 240,
    fidelity: Optional[LigationFidelityCalculator] = None,
    seed: int = 42,
) -> pd.DataFrame:
    """Find high-fidelity cuts in an optimized CDS and return fragment plan + DNA."""
    cds = clean_dna(cds)
    if n_fragments is None:
        n_fragments = suggest_fragment_count(len(cds), max_fragment_cds=max_fragment_cds)
    calc = fidelity or LigationFidelityCalculator(25, 18)
    cuts, ohs, fid = search_cut_sites(cds, int(n_fragments), calc, seed=seed)
    return plan_fragments_from_cds(cds, aa_sequence, cuts, ohs, ligation_fidelity=fid)
