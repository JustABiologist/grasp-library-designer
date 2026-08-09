# GRASP Library Designer

Codon-optimize [GRASP](https://academic.oup.com/nar/article/53/20/gkaf1169/8321212) (Farley et al., *NAR* 2025) binder DNA for Golden Gate assembly.

**PyPI package:** `grasp-library-designer`  
**Import name:** `grasp_library`

Two Colab Forms notebooks:

| Notebook | Purpose |
|---|---|
| [`grasp_oneshot_designer.ipynb`](grasp_oneshot_designer.ipynb) | **One target RNA** → binder protein → free GGA cut sites → oligos |
| [`grasp_library_designer.ipynb`](grasp_library_designer.ipynb) | Redesign / anneal the **42-module combinatorial library**, then GAP-compile a target |

Hard constraints (library path): protein sequence fixed (synonymous codons only); coding Golden Gate overhang bases stay locked in `coding_mask`. Objectives: ligation fidelity (Potapov / GGAssembler), codon optimality, synthesis fitness.

> **License:** AGPL-3.0 (required by the vendored GGAssembler / dawdlib ligation engine). See [`LICENSE`](LICENSE) and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

---

## Install (PyPI)

```bash
pip install grasp-library-designer
# optional notebook extras
pip install "grasp-library-designer[notebook]"
```

Minimal API check:

```python
from grasp_library import materialize_project, build_default_config, LigationFidelityCalculator

project = materialize_project()          # creates ./grasp_library_project + GenBank
config = build_default_config(project / "input")
print(LigationFidelityCalculator(25, 18).set_fidelity(["AATG", "GATA"]))
```

Until the package is published on PyPI, install from GitHub (private repo needs a PAT):

```bash
pip install "git+https://<TOKEN>@github.com/JustABiologist/grasp-library-designer.git@main"
```

Or clone and install editable:

```bash
git clone https://github.com/JustABiologist/grasp-library-designer.git
cd grasp-library-designer
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[notebook,dev]"
```

---

## Run in Google Colab

Both notebooks use **Colab Forms** (`#@title` / `#@param`, `{display-mode: "form"}`).

### 1. Open a notebook

From GitHub → **Open in Colab**, or upload the `.ipynb`.

### 2. Install

In **0 · Install**, choose:

| Mode | When |
|---|---|
| **PyPI** | After the package is on PyPI (`pip install grasp-library-designer`) |
| **Private GitHub** | Clone this private repo with a `repo`-scoped PAT |
| **Local editable** | Notebook already lives in a checkout |

### 3a. One-shot

Settings → Preview binder → Design oligos → Export Excel  
Outputs: `grasp_library_project/output/oneshot/{RNA}/`

### 3b. Combinatorial library

Settings → Import → Redesign overhangs → Anneal → Pareto plot → Export → Compile target  
Outputs: `grasp_library_project/output/`

Bundled Farley et al. GenBank modules are copied into the project folder on first run via `materialize_project()`.

---

## Run locally (Cursor / Jupyter / VS Code)

```bash
pip install -e ".[notebook]"
python -m ipykernel install --user --name grasp-library-designer --display-name "grasp-library-designer"
```

Select that kernel, open either notebook, run top-to-bottom.

---

## Package layout

```
grasp_library/                 # installable Python package
  data/profiles/.../genbank/   # bundled GRASP GenBank modules
  paths.py                     # materialize_project()
  ...
third_party/dawdlib_golden_gate/   # Potapov ligation fidelity (AGPL; also installed)
grasp_*_designer.ipynb             # Colab Forms UIs (also in sdist)
```

---

## Build / publish (maintainers)

```bash
pip install -e ".[dev]"
python -m build
twine check dist/*
# Test PyPI first (recommended):
twine upload --repository testpypi dist/*
# Production:
twine upload dist/*
```

Requires a PyPI API token (`TWINE_USERNAME=__token__`, `TWINE_PASSWORD=pypi-...`).

---

## Hiding code (Colab / Jupyter / VS Code)

| Frontend | How |
|---|---|
| **Google Colab** | Forms: `#@title … {display-mode: "form"}` + `#@param` |
| **Cursor / VS Code** | **Notebook: Collapse All Cell Inputs** |

**Hide ≠ protect.** Source remains in the `.ipynb`.

---

## License notes

- Distributed package license: **AGPL-3.0** (see [`LICENSE`](LICENSE)).
- Vendored ligation engine under `third_party/dawdlib_golden_gate/` is AGPL-3.0 (Fleishman-Lab / GGAssembler).
- GRASP sequences: Farley et al., *Nucleic Acids Res.* 2025.

---

## Quick smoke test

```bash
python - <<'PY'
from pathlib import Path
from grasp_library import (
    materialize_project,
    build_default_config,
    run_oneshot_design,
    LigationFidelityCalculator,
)
from grasp_library.codon_tables import apply_organism_codon_table, load_codon_usage

project = materialize_project()
input_dir = project / "input"
cfg = build_default_config(input_dir)
cfg["optimizer"]["iterations_per_part"] = 200
apply_organism_codon_table("Escherichia coli (Kazusa)", input_dir / "codon_usage.csv")
_, codon_data = load_codon_usage(input_dir / "codon_usage.csv", genetic_code=1)
run_oneshot_design(
    target_rna="UUACACGUG",
    codon_data=codon_data,
    config=cfg,
    output_dir=project / "output" / "oneshot" / "UUACACGUG",
    n_fragments=4,
    fidelity=LigationFidelityCalculator(25, 18),
)
print("ok")
PY
```
