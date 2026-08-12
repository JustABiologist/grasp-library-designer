# GRASP profile — Farley et al., NAR 2025

Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212  
Data: https://github.com/farleykvdg/GRASP

## Contents

- `genbank/GRASP_-1.gb` — 42 level −1 modules (multi-record GenBank)
- `genbank/pPR-1_*.gb` — same modules as individual plasmids
- `GRASP_AP.py` — upstream Assembly Planner (GAP) part-picking logic

## Regenerating notebook inputs

```bash
python -m grasp_library.import_grasp
```

Writes into `grasp_library_project/input/`:

| File | Role |
|---|---|
| `parts.csv` | AA, coding masks, pAGM1311 BsaI order arms, and overhang coordinates |
| `parts_full.csv` | Sidecar with native CDS + overhang coordinates |
| `junction_map.csv` | Native cut coordinates from the deposited modules |
| `overhang_candidates.csv` | Native fixed-cut candidates retained for provenance / compatibility |
| `target_map.csv` | Module catalog (plans via `compile_target_gap`) |

## Native 9S overhang set

`AGGT – ACTC – AAGA – GCAC – TGAA – CTTC – ACTC – AAGA – GCAC – TGAA – TTCG`

B/C/D are shared across every five-part Level 0 block. The runtime redesign is
limited to cuts fully contained within invariant `ARELF`, but it is not limited
to the paper's native cut indices: all motif-relative starts 0–11 are explored.
Candidate identity includes both the four-base overhang and its ARELF offset.
Prefer `1A_*_AGGT` for the deposited 5′ fusion; `AATG` variants remain
available as alternate start-compatible parts.

Deposited subsequent MoClo block chains are:

- 9S: `AGGT–CDS1–CTTC–CDS2–TTCG`
- 14S: `AGGT–CDS1–GTGA–CDS14–CTTC–CDS2–TTCG`
- 19S: `AGGT–CDS1–GTGA–CDS14–CACG–CDS19–CTTC–CDS2–TTCG`

## Physical cloning path

Order each part as a double-stranded synthesis fragment with inward-facing BsaI recognition sites. The following top strand is derived from the paper's pAGM1311 primer geometry and verified against the deposited vector maps:

`TTTGGTCTCAACAT{pAGM1311 insert}TTGTTGAGACCAAA`

The dashboard asks only for physical 5′ and 3′ overhangs at Level −1, Level 0,
and Level 1. Every value is written 5′→3′. Deposited defaults are `ACAT/ACAA`,
`CTCA/CTCG`, and `GGAG/AGCG`; retained coding-strand sites are tracked
separately. Without a custom acceptor sequence, the software can check insert
geometry and interface requirements but cannot claim that the custom vector
backbone was sequence-simulated.
