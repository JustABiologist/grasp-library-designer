# GRASP Library Designer

Codon-optimize [GRASP](https://academic.oup.com/nar/article/53/20/gkaf1169/8321212) (Farley et al., *NAR* 2025) binder DNA for Golden Gate assembly.

**PyPI:** [`grasp-library-designer`](https://pypi.org/project/grasp-library-designer/) · **Import:** `grasp_library`

---

## Open in Google Colab

Click a badge → run **0 · Install** (PyPI) → fill the forms top to bottom. No GitHub token needed.

| Notebook | Open |
|---|---|
| **One-shot** (one RNA → free GGA oligos) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_oneshot_designer.ipynb) |
| **Library** (42-module redesign → GAP compile) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_library_designer.ipynb) |

Direct links:

- One-shot: https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_oneshot_designer.ipynb
- Library: https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/main/grasp_library_designer.ipynb

Each notebook installs with:

```python
%pip install -q -U "grasp-library-designer>=0.1.5"
```

Bundled GenBank modules and Potapov ligation tables ship inside the PyPI package (`materialize_project()`).

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
| [`grasp_oneshot_designer.ipynb`](grasp_oneshot_designer.ipynb) | One target RNA → binder protein → free GGA cut sites → oligos |
| [`grasp_library_designer.ipynb`](grasp_library_designer.ipynb) | Redesign / anneal the 42-module combinatorial library, then GAP-compile a target |

Hard constraints (library path): protein fixed (synonymous codons only); coding Golden Gate overhang bases stay locked in `coding_mask`. Objectives: ligation fidelity (Potapov / GGAssembler), codon optimality, synthesis fitness.

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
