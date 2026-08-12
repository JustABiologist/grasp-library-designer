# Third-party components

## Fleishman-Lab / GGAssembler (dawdlib golden-gate fidelity)

- Path: `third_party/dawdlib_golden_gate/`
- Upstream: https://github.com/Fleishman-Lab/GGAssembler
- Related package: https://github.com/Fleishman-Lab/dawdlib
- License: AGPL-3.0
- What was vendored: `gate_data.py` and Potapov/NEB ligation frequency CSVs used by `GGData.reaction_fidelity`

Four-base T4-ligase-only frequency data originate from Potapov et al.,
*ACS Synthetic Biology* (2018), DOI
[`10.1021/acssynbio.8b00333`](https://doi.org/10.1021/acssynbio.8b00333),
including the bundled 25 °C and 37 °C, 1 h and 18 h datasets.

The bundled BsaI-HFv2 and BbsI-HF whole-Golden-Gate 37↔16 °C cycling matrices
originate from Pryor et al., *PLOS ONE* (2020), DOI
[`10.1371/journal.pone.0238592`](https://doi.org/10.1371/journal.pone.0238592),
supplementary datasets S1 (BsaI-HFv2) and S4 (BbsI-HF). GRASP uses them as
enzyme-matched stage proxies:

- Level −1 entry cloning → BsaI-HFv2
- Level 0 five-part assemblies → BbsI-HF (published isoschizomer of BpiI)
- Level 1 block joining → BsaI-HFv2

They are labelled proxies, not exact measurements of the GRASP reaction
formulation (GRASP used 26× 3 min/4 min cycles rather than Pryor's 30× 5 min).
