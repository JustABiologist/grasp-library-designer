# GRASP Library Designer

Codon-optimize [GRASP](https://academic.oup.com/nar/article/53/20/gkaf1169/8321212) (Farley et al., *NAR* 2025) binder DNA for Golden Gate assembly.

Two notebooks:

| Notebook | Purpose |
|---|---|
| [`grasp_oneshot_designer.ipynb`](grasp_oneshot_designer.ipynb) | **One target RNA** → binder protein → free GGA cut sites → oligos (**Colab Forms UI**) |
| [`grasp_library_designer.ipynb`](grasp_library_designer.ipynb) | Redesign / anneal the **42-module combinatorial library**, then compile targets |

Hard constraints: protein sequence fixed (synonymous codons only); coding Golden Gate overhang bases stay locked in `coding_mask`. Objectives: ligation fidelity (Potapov / GGAssembler), codon optimality, synthesis fitness.

---

## Run in Google Colab (recommended for biologists)

These notebooks are set up for **Colab Forms** (`#@title` / `#@param`). Code stays hidden by default via `{display-mode: "form"}`. You can also use the cell ⋮ menu → **Form → Hide code**.

### 1. Open the notebook in Colab

From GitHub (private repo): open the `.ipynb` on GitHub → **Open in Colab**, or upload the notebook file to Colab.

Primary entry point: **`grasp_oneshot_designer.ipynb`**.

### 2. Install the private package

In the **0 · Install** cell:

1. Set `REPO_SLUG` to `JustABiologist/grasp-library-designer` (or whatever you named the repo).
2. Paste a GitHub **personal access token** with `repo` scope into `GITHUB_TOKEN` (Colab needs this to clone a private repo).
3. Set `BRANCH` to `main` if that is your default branch.
4. Run the cell.

The cell clones into `/content/grasp-library-designer` and runs `pip install -e .`.

> Never commit a token into the notebook. Treat Colab “hide code” as presentation only — anyone with the `.ipynb` can reveal source.

### 3. Settings → preview → design → export

1. **1 · Settings** — organism, vendor, ligation table, target RNA, anneal depth, fragment count (`0` = auto).
2. **2 · Preview binder protein** — shows PPR code and AA length.
3. **3 · Design oligos** — full CDS anneal, high-fidelity cut search, oligo table.
4. **4 · Export Excel** — writes under `grasp_library_project/output/oneshot/{RNA}/` and triggers a Colab download.

### Combinatorial library in Colab

Open `grasp_library_designer.ipynb`, run **Colab · Install**, then **Colab · Settings (Forms)**. Continue with the import / redesign / anneal cells further down (those still use the shared `CONFIG` / `CODON_DATA`).

---

## Run locally (Cursor / Jupyter / VS Code)

```bash
git clone https://github.com/JustABiologist/grasp-library-designer.git
cd grasp-library-designer
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[notebook]"
python -m ipykernel install --user --name grasp-library-designer --display-name "grasp-library-designer"
```

Select the `grasp-library-designer` kernel, open either notebook, run top-to-bottom.

- **One-shot:** Colab Forms cells also work locally (they are plain Python assignments).
- **Library:** uses an `ipywidgets` control panel; collapse long helper cells with **Notebook: Collapse All Cell Inputs** if you want a cleaner UI.

---

## Package layout

```
grasp_library/           # Python package
  binder.py              # RNA → PPR → binder AA (no library parts)
  oneshot.py             # protein-first one-shot pipeline
  gga_split.py           # free GGA cut / overhang search
  optimizer.py           # masked codon + synthesis SA
  import_grasp.py        # GenBank → parts / GAP compile (library path)
  ...
grasp_library_project/
  input/                 # parts.csv, junction_map, codon tables, …
  profiles/grasp_nar2025/genbank/   # deposited GRASP modules
  output/                # gitignored runtime results
third_party/dawdlib_golden_gate/    # Potapov ligation fidelity (AGPL)
```

---

## Hiding code (Colab / Jupyter / VS Code)

| Frontend | How |
|---|---|
| **Google Colab** | Forms: `#@title … {display-mode: "form"}` + `#@param`. Cell ⋮ → **Form → Hide code**. |
| **JupyterLab** | Collapse inputs; optional metadata `"jupyter": {"source_hidden": true}`. |
| **Cursor / VS Code** | Command Palette → **Notebook: Collapse All Cell Inputs**. |
| **Classic Notebook 6** | Often needs extensions (`hide_input`, etc.). |

**Hide ≠ protect.** Source remains in the `.ipynb`. For a true app UI with no notebook chrome, use Voilà, Streamlit, Panel, or Gradio.

---

## License notes

- Project code: use as you see fit in this private repo unless you add a license file.
- Vendored ligation engine under `third_party/dawdlib_golden_gate/` is **AGPL-3.0** (Fleishman-Lab / GGAssembler). See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
- GRASP sequences: Farley et al., *Nucleic Acids Res.* 2025; GenBank deposit linked from the paper.

---

## Quick smoke test (CLI)

```bash
source .venv/bin/activate
python - <<'PY'
from pathlib import Path
from grasp_library import build_default_config, run_oneshot_design, LigationFidelityCalculator
from grasp_library.codon_tables import apply_organism_codon_table, load_codon_usage

input_dir = Path("grasp_library_project/input")
cfg = build_default_config(input_dir)
cfg["optimizer"]["iterations_per_part"] = 200
apply_organism_codon_table("Escherichia coli (Kazusa)", input_dir / "codon_usage.csv")
_, codon_data = load_codon_usage(input_dir / "codon_usage.csv", genetic_code=1)
run_oneshot_design(
    target_rna="UUACACGUG",
    codon_data=codon_data,
    config=cfg,
    output_dir=Path("grasp_library_project/output/oneshot/UUACACGUG"),
    n_fragments=4,
    fidelity=LigationFidelityCalculator(25, 18),
)
print("ok")
PY
```
