"""Joint Golden Gate gene design: cuts, overhangs, and coding sequence together.

One-shot assembly no longer borrows GRASP library modules or ARELF junctions.
Given a protein, this module searches codon-aligned cut sites and the 4-nt
overhangs those cuts can encode, then recodes the rest of the CDS for codon
optimality, cut-site depletion, and synthesis — in one pass.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .dna import (
    clean_dna,
    is_self_reverse_complement,
    reverse_complement,
)
from .ligation_fidelity import LigationFidelityCalculator
from .optimizer import codon_score, optimize_coding_sequence


TYPE_IIS_GEOMETRY: Dict[str, Dict[str, object]] = {
    "BsaI": {
        "enzyme": "BsaI",
        "recognition_site": "GGTCTC",
        "spacer_len": 1,
        "overhang_len": 4,
        "clamp_5p": "TTT",
        "clamp_3p": "AAA",
        "fidelity_level": "level_minus1",
    },
    "BbsI": {
        "enzyme": "BbsI",
        "recognition_site": "GAAGAC",
        "spacer_len": 2,
        "overhang_len": 4,
        "clamp_5p": "TTT",
        "clamp_3p": "AAA",
        "fidelity_level": "level0",
    },
    "BpiI": {
        "enzyme": "BpiI",
        "recognition_site": "GAAGAC",
        "spacer_len": 2,
        "overhang_len": 4,
        "clamp_5p": "TTT",
        "clamp_3p": "AAA",
        "fidelity_level": "level0",
    },
    "BsmBI": {
        "enzyme": "BsmBI",
        "recognition_site": "CGTCTC",
        "spacer_len": 1,
        "overhang_len": 4,
        "clamp_5p": "TTT",
        "clamp_3p": "AAA",
        "fidelity_level": "level_minus1",
    },
}

_WRAP_ALIASES = {
    "bsai": "BsaI",
    "bsai (ggtctc)": "BsaI",
    "bbsi": "BbsI",
    "bpii": "BpiI",
    "bpii / bbsi": "BpiI",
    "bpii / bbsi (gaagac)": "BpiI",
    "bbsi-hf": "BbsI",
    "bsmbi": "BsmBI",
    "esp3i": "BsmBI",
    "bsmbi / esp3i": "BsmBI",
    "bsmbi / esp3i (cgtctc)": "BsmBI",
}


def resolve_wrap_enzyme(name: str | None) -> str:
    if not name:
        return "BsaI"
    key = str(name).strip()
    if key in TYPE_IIS_GEOMETRY:
        return key
    alias = _WRAP_ALIASES.get(key.lower())
    if alias is None:
        raise ValueError(
            f"Unsupported Golden Gate wrap enzyme {name!r}; "
            f"expected one of {sorted(TYPE_IIS_GEOMETRY)}"
        )
    return alias


def wrap_geometry(name: str | None) -> Dict[str, object]:
    return dict(TYPE_IIS_GEOMETRY[resolve_wrap_enzyme(name)])


def oligo_flank_overhead(geometry: Mapping[str, object]) -> int:
    """Bases added around the payload (clamps, sites, spacers — not overhangs)."""
    site = len(str(geometry["recognition_site"]))
    spacer = int(geometry["spacer_len"])
    return (
        len(str(geometry["clamp_5p"]))
        + site
        + spacer
        + spacer
        + site
        + len(str(geometry["clamp_3p"]))
    )


def wrap_payload(payload: str, geometry: Mapping[str, object]) -> str:
    """Inward-facing Type IIS arms around a payload that already starts/ends with overhangs."""
    payload = clean_dna(payload)
    site = clean_dna(str(geometry["recognition_site"]))
    spacer = int(geometry["spacer_len"])
    prefix = (
        clean_dna(str(geometry["clamp_5p"])) + site + ("A" * spacer)
    )
    suffix = (
        ("T" * spacer) + reverse_complement(site) + clean_dna(str(geometry["clamp_3p"]))
    )
    return prefix + payload + suffix


def suggest_fragment_count(
    cds_length: int,
    *,
    max_fragment_cds: int = 240,
    min_fragments: int = 1,
    max_fragments: int = 8,
) -> int:
    n = max(min_fragments, math.ceil(max(cds_length, 1) / max(max_fragment_cds, 1)))
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

    Sequential fallback used when the CDS is already fixed. Prefer
    ``co_design_gene_assembly`` when synonymous recoding is allowed.
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
    """Build fragment rows (AA + coding_mask + DNA slice) from a full CDS.

    ``oh5``/``oh3`` are four-base sites as written on the assembled coding
    strand.  They are not independently oriented physical terminal labels;
    callers that model a particular vector/enzyme geometry must derive and
    report those terminals separately.
    """
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
        if i > 0:
            dna0 = dna0 - 4
        frag_dna = cds[dna0:dna1]
        frag_aa = aa[aa0:aa1]
        mask = ["N"] * len(frag_dna)
        oh5 = overhangs[i - 1] if i > 0 else ""
        oh3 = overhangs[i] if i < len(overhangs) else ""
        if oh5:
            mask[:4] = list(oh5)
        if oh3:
            mask[-4:] = list(oh3)

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
                "oh5_coding_site_5to3": oh5 or None,
                "oh3_coding_site_5to3": oh3 or None,
                "overhang_notation": "assembled_coding_strand_site",
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


# ---------------------------------------------------------------------------
# Joint co-design
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JunctionOption:
    overhang: str
    codon_score: float
    min_adaptiveness: float
    gc: int


@dataclass
class AssemblyLayout:
    cuts: List[int]
    overhangs: List[str]
    junction_codon_score: float
    ligation_fidelity: float = 0.0


@dataclass
class GeneAssembly:
    cds: str
    aa_sequence: str
    cuts: List[int]
    overhangs: List[str]
    destination_5prime: str
    destination_3prime: str
    destination_3prime_coding: str
    ligation_fidelity: float
    codon_score: float
    sequence_score: float
    coding_mask: str
    wrap_enzyme: str
    geometry: Dict[str, object]
    payloads: List[str]
    oligos: List[str]
    rescue_notes: List[dict] = field(default_factory=list)


def _clean_overhang(value: str, label: str) -> str:
    cleaned = clean_dna(value)
    if len(cleaned) != 4:
        raise ValueError(f"{label} must be exactly four DNA bases (got {value!r})")
    return cleaned


def _overhang_gc(overhang: str) -> int:
    return overhang.count("G") + overhang.count("C")


def junction_options_at(
    aa_sequence: str,
    aa_end: int,
    codon_data: Mapping[str, Sequence[dict]],
    *,
    min_adaptiveness: float = 0.0,
) -> List[JunctionOption]:
    """Synonymous 4-nt overhangs encodable at a codon-aligned cut after ``aa_end``."""
    if aa_end < 2 or aa_end > len(aa_sequence):
        return []
    wobble_aa = aa_sequence[aa_end - 2]
    full_aa = aa_sequence[aa_end - 1]
    if wobble_aa not in codon_data or full_aa not in codon_data:
        raise ValueError(
            f"Codon table missing amino acids {wobble_aa}/{full_aa} at cut {aa_end}"
        )
    best: Dict[str, JunctionOption] = {}
    for left in codon_data[wobble_aa]:
        for right in codon_data[full_aa]:
            overhang = str(left["codon"])[2] + str(right["codon"])
            left_w = float(left["relative_adaptiveness"])
            right_w = float(right["relative_adaptiveness"])
            if min(left_w, right_w) < min_adaptiveness:
                continue
            score = math.log(max(left_w, 1e-12)) + math.log(max(right_w, 1e-12))
            current = best.get(overhang)
            if current is None or score > current.codon_score:
                best[overhang] = JunctionOption(
                    overhang=overhang,
                    codon_score=score,
                    min_adaptiveness=min(left_w, right_w),
                    gc=_overhang_gc(overhang),
                )
    return list(best.values())


def _rank_junction_options(
    options: Sequence[JunctionOption],
    *,
    used: set[str],
    prefer_adaptiveness: float,
) -> List[JunctionOption]:
    ranked = []
    for option in options:
        oh = option.overhang
        if len(oh) != 4 or is_self_reverse_complement(oh):
            continue
        rc = reverse_complement(oh)
        if oh in used or rc in used:
            continue
        ranked.append(option)
    ranked.sort(
        key=lambda opt: (
            opt.min_adaptiveness >= prefer_adaptiveness,
            1 <= opt.gc <= 3,
            opt.codon_score,
        ),
        reverse=True,
    )
    return ranked


def _position_candidates(lo: int, hi: int, ideal: int, limit: int = 14) -> List[int]:
    if hi < lo:
        return []
    if hi - lo + 1 <= limit:
        return list(range(lo, hi + 1))
    pts = {ideal} if lo <= ideal <= hi else set()
    # Denser near the even split, sparse toward the edges.
    for fraction in (0.0, 0.15, 0.3, 0.45, 0.5, 0.55, 0.7, 0.85, 1.0):
        pts.add(int(round(lo + fraction * (hi - lo))))
    extra = max(0, limit - len(pts))
    if extra:
        step = max(1, (hi - lo) // (extra + 1))
        pts.update(range(lo, hi + 1, step))
    return sorted(p for p in pts if lo <= p <= hi)[:limit]


def _length_penalty(cuts: Sequence[int], n_aa: int) -> float:
    bounds = [0, *cuts, n_aa]
    lengths = [b - a for a, b in zip(bounds, bounds[1:])]
    if not lengths:
        return 0.0
    mean = sum(lengths) / len(lengths)
    return float(sum((length - mean) ** 2 for length in lengths) / len(lengths))


def _reaction_overhangs(
    destination_5prime: str,
    destination_3prime_coding: str,
    junctions: Sequence[str],
) -> List[str]:
    return [destination_5prime, *junctions, destination_3prime_coding]


def _score_fidelity(
    fidelity: LigationFidelityCalculator, overhangs: Sequence[str]
) -> float:
    if not _valid_overhang_set(overhangs):
        return -1.0
    try:
        return float(fidelity.set_fidelity(overhangs))
    except (KeyError, ValueError):
        return -1.0


def _pair_ok(
    overhang: str,
    chosen: Sequence[str],
    fidelity: LigationFidelityCalculator,
    cache: Dict[Tuple[str, str], float],
    *,
    min_pair: float = 0.05,
) -> bool:
    for other in chosen:
        key = (overhang, other) if overhang <= other else (other, overhang)
        score = cache.get(key)
        if score is None:
            try:
                score = float(fidelity.pair_fidelity(overhang, other))
            except (KeyError, ValueError):
                score = 0.0
            cache[key] = score
        if score < min_pair:
            return False
    return True


def _aa_limits(
    n_aa: int,
    n_fragments: int,
    min_aa: int,
    max_aa: int,
) -> Tuple[int, int]:
    if n_aa < n_fragments * min_aa:
        raise ValueError(
            f"Protein too short ({n_aa} aa) for {n_fragments} fragments "
            f"of at least {min_aa} aa"
        )
    needed = max(1, math.ceil(n_aa / max(max_aa, 1)))
    if n_fragments < needed:
        raise ValueError(
            f"{n_fragments} fragments cannot keep every oligo within the length "
            f"cap ({n_aa} aa, max {max_aa} aa / fragment). Use at least "
            f"{needed} fragments or raise max_oligo_length."
        )
    return min_aa, max_aa


def _feasible_cut_window(
    *,
    index: int,
    n_junc: int,
    last_cut: int,
    n_aa: int,
    min_aa: int,
    max_aa: int,
) -> Tuple[int, int]:
    remaining_after = n_junc - index
    lo = max(last_cut + min_aa, n_aa - max_aa * remaining_after)
    hi = min(last_cut + max_aa, n_aa - min_aa * remaining_after)
    return lo, hi


def _build_junction_catalog(
    aa_sequence: str,
    codon_data: Mapping[str, Sequence[dict]],
    min_adaptiveness: float,
) -> Dict[int, List[JunctionOption]]:
    catalog: Dict[int, List[JunctionOption]] = {}
    # Cuts need two upstream residues to form a 4-nt overhang.
    for aa_end in range(2, len(aa_sequence)):
        preferred = junction_options_at(
            aa_sequence, aa_end, codon_data, min_adaptiveness=min_adaptiveness
        )
        fallback = (
            junction_options_at(aa_sequence, aa_end, codon_data, min_adaptiveness=0.0)
            if not preferred
            else preferred
        )
        if fallback:
            catalog[aa_end] = fallback
    return catalog


def _beam_search_layouts(
    *,
    n_aa: int,
    n_fragments: int,
    min_aa: int,
    max_aa: int,
    catalog: Mapping[int, Sequence[JunctionOption]],
    dest5: str,
    dest3_coding: str,
    fidelity: LigationFidelityCalculator,
    prefer_adaptiveness: float,
    beam_width: int,
    options_per_cut: int,
) -> List[AssemblyLayout]:
    n_junc = n_fragments - 1
    if n_junc <= 0:
        reaction = _reaction_overhangs(dest5, dest3_coding, [])
        return [
            AssemblyLayout(
                cuts=[],
                overhangs=[],
                junction_codon_score=0.0,
                ligation_fidelity=_score_fidelity(fidelity, reaction),
            )
        ]

    used0 = {dest5, reverse_complement(dest5), dest3_coding, reverse_complement(dest3_coding)}
    pair_cache: Dict[Tuple[str, str], float] = {}
    chosen0 = [dest5, dest3_coding]
    beam: List[Tuple[float, List[int], List[str], List[float], set[str]]] = [
        (0.0, [], [], [], set(used0))
    ]
    ideals = _default_cut_aa_positions(n_aa, n_fragments)

    for index in range(n_junc):
        nxt: List[Tuple[float, List[int], List[str], List[float], set[str]]] = []
        for _, cuts, ohs, scores, used in beam:
            last = cuts[-1] if cuts else 0
            lo, hi = _feasible_cut_window(
                index=index,
                n_junc=n_junc,
                last_cut=last,
                n_aa=n_aa,
                min_aa=min_aa,
                max_aa=max_aa,
            )
            ideal = ideals[index] if index < len(ideals) else (lo + hi) // 2
            for pos in _position_candidates(lo, hi, ideal):
                ranked = _rank_junction_options(
                    catalog.get(pos, ()),
                    used=used,
                    prefer_adaptiveness=prefer_adaptiveness,
                )[:options_per_cut]
                for option in ranked:
                    if not _pair_ok(option.overhang, chosen0 + ohs, fidelity, pair_cache):
                        continue
                    new_ohs = ohs + [option.overhang]
                    new_cuts = cuts + [pos]
                    new_scores = scores + [option.codon_score]
                    new_used = set(used)
                    new_used.add(option.overhang)
                    new_used.add(reverse_complement(option.overhang))
                    heuristic = (
                        sum(new_scores)
                        - 0.02 * _length_penalty(new_cuts, n_aa)
                        + (0.15 if 1 <= option.gc <= 3 else 0.0)
                    )
                    nxt.append((heuristic, new_cuts, new_ohs, new_scores, new_used))
        if not nxt:
            break
        nxt.sort(key=lambda item: item[0], reverse=True)
        beam = nxt[:beam_width]

    layouts: List[AssemblyLayout] = []
    for _, cuts, ohs, scores, _used in beam:
        if len(cuts) != n_junc:
            continue
        reaction = _reaction_overhangs(dest5, dest3_coding, ohs)
        fid = _score_fidelity(fidelity, reaction)
        if fid < 0:
            continue
        layouts.append(
            AssemblyLayout(
                cuts=list(cuts),
                overhangs=list(ohs),
                junction_codon_score=float(sum(scores) / max(len(scores), 1)),
                ligation_fidelity=fid,
            )
        )
    layouts.sort(
        key=lambda layout: (
            layout.ligation_fidelity,
            layout.junction_codon_score,
            -_length_penalty(layout.cuts, n_aa),
        ),
        reverse=True,
    )
    return layouts


def coding_mask_for_cuts(
    n_aa: int, cuts: Sequence[int], overhangs: Sequence[str]
) -> str:
    mask = ["N"] * (3 * n_aa)
    for aa_end, overhang in zip(cuts, overhangs):
        start = 3 * int(aa_end) - 4
        oh = clean_dna(overhang)
        if start < 0 or start + 4 > len(mask):
            raise ValueError(f"Cut aa={aa_end} does not fit a 4-nt overhang")
        mask[start : start + 4] = list(oh)
    return "".join(mask)


def fragment_payloads(
    cds: str,
    cuts: Sequence[int],
    *,
    destination_5prime: str,
    destination_3prime_coding: str,
) -> List[str]:
    cds = clean_dna(cds)
    dest5 = clean_dna(destination_5prime)
    dest3 = clean_dna(destination_3prime_coding)
    if not cuts:
        return [dest5 + cds + dest3]
    payloads = [dest5 + cds[: 3 * cuts[0]]]
    for left, right in zip(cuts, cuts[1:]):
        payloads.append(cds[3 * left - 4 : 3 * right])
    payloads.append(cds[3 * cuts[-1] - 4 :] + dest3)
    return payloads


def assemble_payloads_to_cds(
    payloads: Sequence[str],
    *,
    destination_5prime: str,
    destination_3prime_coding: str,
) -> str:
    if not payloads:
        raise ValueError("No fragments to assemble")
    dest5 = clean_dna(destination_5prime)
    dest3 = clean_dna(destination_3prime_coding)
    if not payloads[0].startswith(dest5):
        raise AssertionError("first fragment missing the destination 5′ overhang")
    if not payloads[-1].endswith(dest3):
        raise AssertionError("last fragment missing the destination 3′ coding site")
    gene = payloads[0][len(dest5) :]
    for payload in payloads[1:]:
        if payload[:4] != gene[-4:]:
            raise AssertionError(
                f"Golden Gate junction mismatch {gene[-4:]} != {payload[:4]}"
            )
        gene += payload[4:]
    if not gene.endswith(dest3):
        raise AssertionError("assembled gene missing the destination 3′ coding site")
    return gene[: -len(dest3)]


def _max_aa_per_fragment(
    *,
    max_oligo: int,
    geometry: Mapping[str, object],
    n_fragments: int,
) -> int:
    overhead = oligo_flank_overhead(geometry)
    # Payload includes two 4-nt overhangs (destination and/or junction).
    max_payload = max_oligo - overhead
    # Internal fragment payload = 3*aa + 4 (5′ overlap). First adds dest5 (+4)
    # already counted; last adds dest3 (+4).
    return max(8, (max_payload - 8) // 3)


def _min_aa_per_fragment(*, min_oligo: int, geometry: Mapping[str, object]) -> int:
    overhead = oligo_flank_overhead(geometry)
    min_payload = max(0, min_oligo - overhead)
    return max(8, (min_payload - 8 + 2) // 3)


def _choose_fragment_count(
    n_aa: int,
    *,
    n_fragments: Optional[int],
    min_aa: int,
    max_aa: int,
    max_fragments: int = 8,
) -> int:
    if n_fragments is not None:
        n = int(n_fragments)
        if n < 1:
            raise ValueError("n_fragments must be ≥ 1")
        return n
    n = max(1, math.ceil(n_aa / max(max_aa, 1)))
    n = min(max_fragments, n)
    if n_aa > n * max_aa:
        n = min(max_fragments, math.ceil(n_aa / max_aa))
    if n > 1 and n_aa < n * min_aa:
        n = max(1, n_aa // min_aa)
    return int(n)


def _refine_layout(
    layout: AssemblyLayout,
    *,
    n_aa: int,
    min_aa: int,
    max_aa: int,
    catalog: Mapping[int, Sequence[JunctionOption]],
    dest5: str,
    dest3_coding: str,
    fidelity: LigationFidelityCalculator,
    prefer_adaptiveness: float,
    used_fixed: set[str],
) -> AssemblyLayout:
    """Move each cut ±3 aa and keep the first fidelity/codon improvement."""
    best = layout
    n_junc = len(layout.cuts)
    for index in range(n_junc):
        left = 0 if index == 0 else best.cuts[index - 1]
        right = n_aa if index == n_junc - 1 else best.cuts[index + 1]
        lo = max(left + min_aa, best.cuts[index] - 3)
        hi = min(right - min_aa, best.cuts[index] + 3, left + max_aa)
        others = [oh for j, oh in enumerate(best.overhangs) if j != index]
        used = set(used_fixed)
        for oh in others:
            used.add(oh)
            used.add(reverse_complement(oh))
        for pos in range(lo, hi + 1):
            ranked = _rank_junction_options(
                catalog.get(pos, ()),
                used=used,
                prefer_adaptiveness=prefer_adaptiveness,
            )[:4]
            for option in ranked:
                cuts = list(best.cuts)
                ohs = list(best.overhangs)
                cuts[index] = pos
                ohs[index] = option.overhang
                if any(
                    nxt - prev < min_aa or nxt - prev > max_aa
                    for prev, nxt in zip([0, *cuts], [*cuts, n_aa])
                ):
                    continue
                reaction = _reaction_overhangs(dest5, dest3_coding, ohs)
                fid = _score_fidelity(fidelity, reaction)
                if fid < best.ligation_fidelity - 1e-9 and fid < 0.99:
                    continue
                codon = float(
                    sum(
                        next(
                            (
                                opt.codon_score
                                for opt in catalog.get(cut, ())
                                if opt.overhang == oh
                            ),
                            0.0,
                        )
                        for cut, oh in zip(cuts, ohs)
                    )
                    / max(len(ohs), 1)
                )
                candidate = AssemblyLayout(
                    cuts=cuts,
                    overhangs=ohs,
                    junction_codon_score=codon,
                    ligation_fidelity=fid,
                )
                better = (
                    fid > best.ligation_fidelity + 1e-6
                    or (
                        abs(fid - best.ligation_fidelity) <= 1e-6
                        and codon > best.junction_codon_score
                    )
                )
                if better:
                    best = candidate
    return best


def co_design_gene_assembly(
    aa_sequence: str,
    codon_data: Mapping[str, Sequence[dict]],
    config: Mapping,
    *,
    n_fragments: Optional[int] = None,
    destination_5prime: str = "CTCA",
    destination_3prime: str = "CTCG",
    wrap_enzyme: str = "BsaI",
    fidelity: Optional[LigationFidelityCalculator] = None,
    max_oligo_length: Optional[int] = None,
    min_oligo_length: Optional[int] = None,
    beam_width: int = 24,
    layouts_to_encode: int = 8,
    iterations: Optional[int] = None,
) -> GeneAssembly:
    """Co-design cut sites, overhangs, and a synonymous CDS for one-pot assembly.

    Destination overhangs are cloning adapters (not part of the protein). Internal
    junctions are 4-nt coding-strand sites produced by synonymous recoding at
    codon-aligned cuts. The body of the gene is then codon-optimized under that
    mask, depleting blacklisted cut sites and applying synthesis heuristics.
    """
    aa = str(aa_sequence).upper().replace(" ", "").replace("*", "")
    if not aa:
        raise ValueError("Empty amino-acid sequence")
    n_aa = len(aa)
    geometry = wrap_geometry(wrap_enzyme)
    dest5 = _clean_overhang(destination_5prime, "destination 5′ overhang")
    dest3_physical = _clean_overhang(destination_3prime, "destination 3′ overhang")
    dest3_coding = reverse_complement(dest3_physical)
    if is_self_reverse_complement(dest5) or is_self_reverse_complement(dest3_coding):
        raise ValueError("Destination overhangs must not be self-reverse-complementary")
    if dest5 == dest3_coding or dest5 == dest3_physical:
        raise ValueError("Destination 5′ and 3′ overhangs must be distinct")

    synthesis = config.get("synthesis") or {}
    max_oligo = int(max_oligo_length or synthesis.get("max_oligo_length") or 300)
    min_oligo = int(min_oligo_length or synthesis.get("min_oligo_length") or 40)
    max_aa = _max_aa_per_fragment(
        max_oligo=max_oligo, geometry=geometry, n_fragments=n_fragments or 1
    )
    min_aa = _min_aa_per_fragment(min_oligo=min_oligo, geometry=geometry)
    n_frag = _choose_fragment_count(
        n_aa, n_fragments=n_fragments, min_aa=min_aa, max_aa=max_aa
    )
    min_aa, max_aa = _aa_limits(n_aa, n_frag, min_aa, max_aa)

    calc = fidelity or LigationFidelityCalculator(25, 18)
    prefer = float(
        (config.get("codon_optimization") or {}).get("minimum_relative_adaptiveness", 0.20)
    )
    catalog = _build_junction_catalog(aa, codon_data, prefer)
    layouts = _beam_search_layouts(
        n_aa=n_aa,
        n_fragments=n_frag,
        min_aa=min_aa,
        max_aa=max_aa,
        catalog=catalog,
        dest5=dest5,
        dest3_coding=dest3_coding,
        fidelity=calc,
        prefer_adaptiveness=prefer,
        beam_width=int(beam_width),
        options_per_cut=4,
    )
    if not layouts:
        raise RuntimeError(
            "Could not co-design a valid Golden Gate overhang set for this "
            "protein. Try fewer fragments, a longer max oligo, or different "
            "destination overhangs."
        )

    used_fixed = {
        dest5,
        reverse_complement(dest5),
        dest3_coding,
        reverse_complement(dest3_coding),
    }
    layouts = [
        _refine_layout(
            layout,
            n_aa=n_aa,
            min_aa=min_aa,
            max_aa=max_aa,
            catalog=catalog,
            dest5=dest5,
            dest3_coding=dest3_coding,
            fidelity=calc,
            prefer_adaptiveness=prefer,
            used_fixed=used_fixed,
        )
        for layout in layouts[: max(layouts_to_encode, 1)]
    ]
    layouts.sort(
        key=lambda layout: (layout.ligation_fidelity, layout.junction_codon_score),
        reverse=True,
    )

    best: Optional[GeneAssembly] = None
    best_total = -1e18
    search_iters = 0
    final_iters = config["optimizer"]["iterations_per_part"] if iterations is None else iterations

    for layout in layouts[: max(layouts_to_encode, 1)]:
        mask = coding_mask_for_cuts(n_aa, layout.cuts, layout.overhangs)
        rescue_notes: List[dict] = []
        try:
            cds, seq_score = optimize_coding_sequence(
                aa,
                mask,
                codon_data,
                config,
                iterations=search_iters,
                rescue_log=rescue_notes,
            )
        except (RuntimeError, ValueError, AssertionError):
            continue
        total = (
            8.0 * layout.ligation_fidelity
            + float(seq_score)
            + 0.5 * layout.junction_codon_score
        )
        if total > best_total:
            best_total = total
            best = GeneAssembly(
                cds=cds,
                aa_sequence=aa,
                cuts=list(layout.cuts),
                overhangs=list(layout.overhangs),
                destination_5prime=dest5,
                destination_3prime=dest3_physical,
                destination_3prime_coding=dest3_coding,
                ligation_fidelity=float(layout.ligation_fidelity),
                codon_score=codon_score(cds, aa, codon_data),
                sequence_score=float(seq_score),
                coding_mask=mask,
                wrap_enzyme=str(geometry["enzyme"]),
                geometry=dict(geometry),
                payloads=[],
                oligos=[],
                rescue_notes=list(rescue_notes),
            )

    if best is None:
        raise RuntimeError(
            "Every high-fidelity layout produced a CDS that still contains a "
            "blacklisted cut site. Relax the site blacklist or move the "
            "destination overhangs."
        )

    if final_iters and final_iters > 0:
        rescue_notes = []
        cds, seq_score = optimize_coding_sequence(
            aa,
            best.coding_mask,
            codon_data,
            config,
            iterations=int(final_iters),
            rescue_log=rescue_notes,
        )
        best.cds = cds
        best.sequence_score = float(seq_score)
        best.codon_score = codon_score(cds, aa, codon_data)
        best.rescue_notes = list(rescue_notes)

    for aa_end, overhang in zip(best.cuts, best.overhangs):
        observed = best.cds[3 * aa_end - 4 : 3 * aa_end]
        if observed != overhang:
            raise AssertionError(
                f"Optimized CDS lost junction overhang {overhang} at aa {aa_end}"
            )

    best.payloads = fragment_payloads(
        best.cds,
        best.cuts,
        destination_5prime=dest5,
        destination_3prime_coding=dest3_coding,
    )
    assembled = assemble_payloads_to_cds(
        best.payloads, destination_5prime=dest5, destination_3prime_coding=dest3_coding
    )
    if assembled != best.cds:
        raise AssertionError("Fragment payloads do not reconstruct the designed CDS")
    best.oligos = [wrap_payload(payload, geometry) for payload in best.payloads]
    for oligo in best.oligos:
        if len(oligo) > max_oligo:
            raise RuntimeError(
                f"Designed oligo is {len(oligo)} bp, above max_oligo_length={max_oligo}. "
                "Increase n_fragments or the oligo length cap."
            )
    return best


def assembly_plan_frame(design: GeneAssembly) -> pd.DataFrame:
    aa = design.aa_sequence
    bounds = [0, *design.cuts, len(aa)]
    n = len(bounds) - 1
    reaction = _reaction_overhangs(
        design.destination_5prime,
        design.destination_3prime_coding,
        design.overhangs,
    )
    rows = []
    for i in range(n):
        aa0, aa1 = bounds[i], bounds[i + 1]
        oh5 = design.destination_5prime if i == 0 else design.overhangs[i - 1]
        oh3 = (
            design.destination_3prime_coding
            if i == n - 1
            else design.overhangs[i]
        )
        rows.append(
            {
                "fragment_id": f"F{i + 1}",
                "assembly_order": i + 1,
                "aa_start_0based": aa0,
                "aa_end_0based": aa1,
                "aa_length": aa1 - aa0,
                "aa_sequence": aa[aa0:aa1],
                "cds_start": 3 * aa0,
                "cds_end": 3 * aa1,
                "payload_5to3": design.payloads[i],
                "oh5_coding_site_5to3": oh5,
                "oh3_coding_site_5to3": oh3,
                "five_prime_end_overhang": oh5,
                "three_prime_end_overhang": reverse_complement(oh3),
                "wrap_enzyme": design.wrap_enzyme,
                "ligation_fidelity_set": design.ligation_fidelity,
                "reaction_overhangs": ";".join(reaction),
            }
        )
    return pd.DataFrame(rows)
