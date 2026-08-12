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

B/C/D are shared across every five-part Level 0 block. Runtime redesign is
restricted to cuts fully contained within invariant `ARELF`, with every
motif-relative start 0–11 eligible; it is not tied to the paper's native cut
indices. Candidate identity is the overhang plus its ARELF offset.

Deposited subsequent MoClo block chains are:

- 9S: `AGGT–CDS1–CTTC–CDS2–TTCG`
- 14S: `AGGT–CDS1–GTGA–CDS14–CTTC–CDS2–TTCG`
- 19S: `AGGT–CDS1–GTGA–CDS14–CACG–CDS19–CTTC–CDS2–TTCG`

Order each part as a double-stranded synthesis fragment with inward-facing BsaI recognition sites. The following top strand is derived from the paper's pAGM1311 primer geometry and verified against the deposited vector maps:

`TTTGGTCTCAACAT{pAGM1311 insert}TTGTTGAGACCAAA`

BsaI clones into pAGM1311 in the deposited preset. The retained BpiI sites then
release the module for five-part assembly into pAGM9121 with external overhangs
`CTCA` and `CGAG`. Dashboard-defined custom interfaces are also supported and
are labelled as requirements-checked, not vector-sequence-verified, unless the
matching acceptor sequence is supplied.
