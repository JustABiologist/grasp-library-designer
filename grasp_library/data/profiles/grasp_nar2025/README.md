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
| `parts.csv` | AA, coding masks, BpiI oligo flanks |
| `parts_full.csv` | Sidecar with native CDS + overhang coordinates |
| `junction_map.csv` | Fixed `mask_start_0based` for 7 unique 9S overhangs |
| `overhang_candidates.csv` | Native + synonym-compatible 4-mers |
| `target_map.csv` | Module catalog (plans via `compile_target_gap`) |

## 9S overhang set (do not move cut indices)

`AGGT – ACTC – AAGA – GCAC – TGAA – CTTC – ACTC – AAGA – GCAC – TGAA – TTCG`

B/C/D are shared across CDS1/CDS2, so redesign uses 7 unique junction variables (`J_Nterm` … `J_Cterm`). Prefer `1A_*_AGGT` for MoClo N-terminal fusion; `AATG` variants are kept as alternate parts.
