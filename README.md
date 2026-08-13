# GRASP Library Designer

Codon-optimize [GRASP](https://academic.oup.com/nar/article/53/20/gkaf1169/8321212) (Farley et al., *NAR* 2025) binder DNA for Golden Gate assembly.

**PyPI:** [`grasp-library-designer`](https://pypi.org/project/grasp-library-designer/) · **Import:** `grasp_library`

---

## Open in Google Colab

Click a badge → run **0 · Install** (PyPI) → fill the forms top to bottom. No GitHub token needed.

| Notebook | Open |
|---|---|
| **One-shot** (one RNA → configured Level −1 order fragments) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_oneshot_designer.ipynb) |
| **Library** (42-module redesign → GAP compile) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_library_designer.ipynb) |

Direct links:

- One-shot: https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_oneshot_designer.ipynb
- Library: https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_library_designer.ipynb

Each notebook installs with:

```python
%pip install -q -U "grasp-library-designer>=0.1.12"
```

Bundled GenBank modules, Potapov ligase-only matrices, and Pryor Golden Gate
cycling matrices ship inside the package (`materialize_project()`).

---

## Install locally

```bash
pip install grasp-library-designer
# optional notebook extras
pip install "grasp-library-designer[notebook]"
```

```python
from grasp_library import materialize_project, build_default_config, LigationFidelityCalculator

project = materialize_project()  # ./grasp_library_project + GenBank
config = build_default_config(project / "input")
print(LigationFidelityCalculator(25, 18).set_fidelity(["AATG", "GATA"]))
```

From a blank Colab / Jupyter, you can also drop the Forms notebooks onto disk:

```python
%pip install -q -U grasp-library-designer
from grasp_library import write_notebook
write_notebook("oneshot")   # or "library"
# then open the written .ipynb from the file browser
```

---

## What each notebook does

| Notebook | Purpose |
|---|---|
| [`grasp_oneshot_designer.ipynb`](grasp_oneshot_designer.ipynb) | One target RNA → target-specific GRASP modules → BsaI order fragments for the configured Level −1 entry vector → BpiI Level 0 blocks |
| [`grasp_library_designer.ipynb`](grasp_library_designer.ipynb) | Redesign / anneal the 42-module combinatorial library, then GAP-compile a target |

Hard constraints (library path): the protein sequence is fixed and every movable
four-base cut is restricted to the invariant `ARELF` motif. The search explores
all motif-relative offsets 0–11 rather than only the four cut positions chosen
in the paper. A candidate is therefore an `(overhang, ARELF offset)` pair, and
each part is rematerialized before codon optimization. Objectives are ligation
fidelity, codon optimality, and synthesis fitness.

Ligation fidelity is reported per physical six-overhang Level 0 reaction (and
optionally as an explicitly labelled product across independently transformed
blocks). The scalar is the orientation-invariant geometric mean of the two
directional products. Stage-matched Pryor et al. 37↔16 °C Golden Gate cycling
matrices are used by default: **BsaI-HFv2** for Levels −1 and 1, and **BbsI-HF**
(BpiI isoschizomer) for Level 0 redesign scoring. Potapov’s ligase-only data
remain available as optional Level 0 surrogates; they contain no measured 16 °C
matrix, so the program does not interpolate or blend static temperature
matrices. These scores are optimization surrogates, not cloning guarantees.
Synthesis QC distinguishes `PASS`, `WARNING`, and `FAIL`; vendor profiles remain
transparent heuristics with `vendor_acceptance_confirmed=False`.

The order file contains double-stranded synthesis fragments with paired,
inward-facing BsaI sites. The dashboard exposes exactly one physical 5′/3′
overhang pair for each cloning level. Every overhang is written 5′→3′. The
deposited GRASP toolbox defaults are:

- Level −1: `ACAT / ACAA`.
- Level 0: `CTCA / CTCG`.
- Level 1: `GGAG / AGCG`.

The 3′ sticky ends are reverse-complemented internally when constructing the
coding-oriented sequence. Thus the retained 3′ coding sites are `TTGT`,
`CGAG`, and `CGCT`, respectively. Internal five-part and ARELF junctions are
derived from the GRASP architecture rather than presented as extra dashboard
overhang fields. When no custom acceptor sequence is provided, the exporter
validates interface requirements but does not claim backbone simulation.

The exported GRASP tract is a PPR block set, not a standalone expression
plasmid. The PPR block-chain check does not validate an entire Level 1
expression construct; promoter, upstream domain, effector, terminator, and
acceptor context must be supplied separately.

For 14S and 19S, intermediate junctions are generated inside the invariant
`ARELF` motif. They are architecture-derived and can be explored by the
overhang redesign search without adding more dashboard fields.

---

## Develop from source

```bash
git clone https://github.com/JustABiologist/grasp-library-designer.git
cd grasp-library-designer
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[notebook,dev]"
```

---

## Package layout

```
grasp_library/                 # installable Python package
  data/profiles/.../genbank/   # bundled GRASP GenBank modules
  notebooks/                   # Colab Forms notebooks (also at repo root)
  paths.py                     # materialize_project()
  ...
third_party/dawdlib_golden_gate/   # Potapov ligation fidelity (AGPL)
```

---

## License

**AGPL-3.0** (required by the vendored GGAssembler / dawdlib ligation engine). See [`LICENSE`](LICENSE) and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

GRASP sequences: Farley et al., *Nucleic Acids Res.* 2025.
