# Developer documentation

This file is for people changing the code. User-facing install and Colab
links stay in [`README.md`](README.md). Scientific caveats (ligation scores are
surrogates, vendor QC is not an acceptance letter) also stay in the README.

PyPI name: `grasp-library-designer`  
Import name: `grasp_library`  
Current version: `0.1.11` (see `pyproject.toml`)  
License: AGPL-3.0 (required by the vendored GGAssembler ligation engine)

---

## 1. What this repo does

Two design paths share one library:

| Path | Entry point | What it produces |
|---|---|---|
| **One-shot** | `run_oneshot_design()` / `grasp_oneshot_designer.ipynb` | One target RNA → GRASP modules → BsaI Level −1 order fragments → BpiI Level 0 blocks |
| **Library** | `run_library_redesign_and_anneal()` / `grasp_library_designer.ipynb` | Redesign / anneal the 42-module combinatorial library, then GAP-compile a target |

Hard constraints on the library path: protein sequence is fixed; movable
four-base cuts stay inside the invariant `ARELF` motif (offsets 0–11).
Objectives are ligation fidelity, codon optimality, and synthesis fitness.

Default physical overhangs (all written 5′→3′):

- Level −1: `ACAT / ACAA`
- Level 0: `CTCA / CTCG`
- Level 1: `GGAG / AGCG`

The 3′ sticky ends are reverse-complemented internally when building the
coding-oriented sequence. Retained 3′ coding sites are therefore `TTGT`,
`CGAG`, and `CGCT`.

---

## 2. Repository layout

```
grasp-library-designer/
├── grasp_library/                 # installable package (this is the product)
│   ├── data/profiles/grasp_nar2025/   # bundled GenBank (shipped on PyPI)
│   ├── notebooks/                 # same Colab notebooks, packaged for write_notebook()
│   └── *.py                       # library code
├── grasp_library_project/         # writable example project tree (also used by tests)
│   ├── input/                     # CSVs + config.yaml produced by import
│   ├── output/                    # gitignored runtime output
│   └── profiles/grasp_nar2025/    # duplicate of bundled GenBank (see §7)
├── tests/                         # pytest
├── third_party/dawdlib_golden_gate/   # vendored GGAssembler fidelity engine (AGPL)
├── grasp_oneshot_designer.ipynb   # Colab badge target (repo root)
├── grasp_library_designer.ipynb
├── pyproject.toml
├── README.md
└── DEV_DOC.md                     # this file
```

Two copies exist on purpose today, not because they should:

1. **Notebooks** live at repo root (Colab `blob/main/...` badges) **and** in
   `grasp_library/notebooks/` (`write_notebook()` after `pip install`).
   They are currently byte-identical. Edit one, copy to the other, or you
   will ship a stale Forms notebook.
2. **GenBank** lives in `grasp_library/data/...` (wheel) **and**
   `grasp_library_project/profiles/...` (checked-in project tree).
   `materialize_project()` copies bundled → project if the dest is missing.

Canonical source of sequence data is the bundled package copy. The project
tree is a working directory.

---

## 3. Runtime data flow

```
bundled GenBank (GRASP_-1.gb, 42 modules)
        │
        ▼
import_grasp.import_grasp_profile()     → parts.csv, parts_full.csv,
                                          junction_map.csv, target_map.csv, …
        │
        ▼
build_default_config() + Colab forms    → nested dict CONFIG
        │
        ├─ one-shot ─► run_oneshot_design()
        │                 compile_target_gap → pick modules
        │                 optimize_library (masked codon SA)
        │                 build_order_fragment + in-silico checks
        │
        └─ library ──► run_library_redesign_and_anneal()
                          run_overhang_redesign (Pareto, optional)
                          run_library_optimize (anneal every module)
                          compile_and_assemble_target (GAP)
```

Config is a nested `dict`, not a typed model. Keys that matter most:

- `assembly_interfaces` — Level −1 / 0 / 1 overhangs and preset
- `overhang_redesign.cut_mode` — `movable_arelf` vs `native_fixed`
- `ligation` / per-level protocols in `synthesis_vendors`
- `optimizer.iterations_per_part` — `0` = greedy only (used during Pareto)
- `synthesis` + `synthesis_vendor_meta` — QC heuristics
- `forbidden_sites` — BsaI / BpiI / BsmBI by default
- `genetic_code`, `target_rna`, `architecture`

`build_default_config(input_dir)` is the constructor. Colab Forms go through
`apply_form_settings()`. The ipywidgets dashboard mutates the same dict.

---

## 4. File guide — `grasp_library/`

Grouped by job, not alphabetically. Starred files are the ones most PRs
touch.

### 4.1 Public surface

| File | Role |
|---|---|
| `__init__.py` | Re-exports ~90 names. Prefer importing from submodules in new code; do not grow `__all__` without a reason. |
| `paths.py` | `bundled_profile_genbank()`, `materialize_project()`, `project_paths()`. Creates `input/`, `output/`, copies GenBank. |
| `colab.py` | `write_notebook("oneshot" \| "library")` — copies packaged `.ipynb` onto disk. |

### 4.2 Domain: RNA → protein → modules

| File | Role |
|---|---|
| `binder.py` | Target RNA (length 9 / 14 / 19) → PPR recognition code → binder amino acids. No library parts required. |
| `import_grasp.py` * | Parse deposited GenBank, build parts tables, GAP `compile_target_gap` / `pick_parts_for_target`, pAGM1311 order-fragment arms. Runnable as `python -m grasp_library.import_grasp`. ~955 lines. |
| `arelf.py` * | Invariant `ARELF` motif, offsets 0–11, candidate `(overhang, offset)` pairs, `materialize_arelf_parts()`. |
| `gga_split.py` | Legacy / generic Golden Gate split of an optimized CDS. One-shot no longer uses configurable fragment counts. |

### 4.3 Domain: cloning geometry

| File | Role |
|---|---|
| `assembly_interfaces.py` * | Presets (`deposited_grasp`, `custom`), physical 5′/3′ overhangs vs assembled coding sites, `build_order_fragment` / `extract_order_payload`. Has its own `reverse_complement` (duplicate of `dna.py`). |
| `oneshot.py` * | `run_oneshot_design()`, in-silico order-fragment + PPR block-chain checks. Does not claim whole-vector or wet-lab validation. |
| `dna.py` | Clean DNA/mask, RC, GC, homopolymer, k-mer penalties, forbidden sites, mask overlay. |

### 4.4 Domain: scoring and optimization

| File | Role |
|---|---|
| `ligation_fidelity.py` * | Wrapper around vendored `GGData.reaction_fidelity`. Stage-matched calculators. `set_fidelity` = geometric mean of both directional products (NEB-viewer style). |
| `synthesis_vendors.py` | Named vendor heuristics, enzyme lists, Pryor vs Potapov ligation tables, `apply_*_to_config`. |
| `objectives.py` | `ObjectiveScores`, codon / synthesis / junction window scores, `evaluate_design()`. |
| `pareto.py` | Overhang-set search, dominance, knee point. Greedy CDS (`iterations=0`) during search. |
| `optimizer.py` * | Masked codon simulated annealing, `optimize_library`, `synthesis_qc` (`PASS` / `WARNING` / `FAIL`), `simulate_assembled_cds`. ~750 lines. |
| `workflows.py` * | Headless library pipeline (no widgets): redesign → anneal → export → compile. ~914 lines. Safe to call from notebooks or a future CLI. |

### 4.5 Codon tables

| File | Role |
|---|---|
| `codon_tables.py` | Load CSV, apply organism table to parts, validate. |
| `codon_validation.py` | Translation vs organism table, cut-site AA risk, issue formatting. |
| `sample_codon_tables.py` | Built-in Kazusa-derived tables shipped in-package. |
| `kazusa.py` | Fetch CUTG HTML by species ID. |
| `codon_upload.py` | Parse uploaded tables; Colab file prompt. |

### 4.6 Notebook UI (not required for library logic)

| File | Role |
|---|---|
| `control_panel.py` * | `GraspControlPanel` ipywidgets dashboard + `build_default_config()`. ~946 lines. |
| `colab_forms.py` | Maps Colab Form fields → CONFIG. Shared by both notebooks. |
| `notebook_ui.py` | Inline-CSS HTML for VS Code / Colab (no `<style>` tags). |
| `plotting.py` | Pareto scatter. Pulls in matplotlib. |

### 4.7 Packaged data

| Path | Role |
|---|---|
| `data/profiles/grasp_nar2025/genbank/GRASP_-1.gb` | 42 Level −1 modules (multi-record). Primary import source. |
| `data/profiles/grasp_nar2025/genbank/pPR-1_*.gb` | Same modules as individual plasmids. |
| `data/profiles/grasp_nar2025/README.md` | Profile notes + how to regenerate CSVs. |
| `notebooks/*.ipynb` | Forms notebooks copied by `write_notebook()`. |
| `py.typed` | Marker only; there is no mypy job yet. |

---

## 5. File guide — everything else

### `tests/`

Domain tests, not coverage theater. Run with `pytest` from the repo root
after `pip install -e ".[dev]"`.

| File | What it locks in |
|---|---|
| `test_oneshot.py` | Deposited parts reconstruct valid pAGM1311 fragments; 9S / 14S / 19S export requirements |
| `test_strand_semantics.py` | Physical sticky end vs coding site; directional terminal pairs |
| `test_assembly_interfaces.py` | Preset geometry, validation, order-fragment arms |
| `test_dashboard_interfaces.py` | Forms + widgets write the same overhang fields |
| `test_arelf_movable.py` | Motif-bounded cuts, rematerialized parts |
| `test_architectures.py` | GAP compile for 9S / 14S / 19S five-part blocks |
| `test_ligation_protocols.py` | Pryor/Potapov table selection; numeric fidelity vs NEB viewer |
| `test_qc_fidelity.py` | `PASS` / `WARNING` / `FAIL`, vendor hard rules, sequence_kind length bounds |

If you change overhang orientation, ligation matrices, ARELF offsets, or
order-fragment arms, these are the tests that should fail first.

### `grasp_library_project/`

Checked-in working tree. `input/config.yaml` is a snapshot of dashboard
defaults, not the runtime source of truth (`build_default_config` is).
Regenerate CSVs with:

```bash
python -m grasp_library.import_grasp
```

Do not commit `output/`.

### `third_party/dawdlib_golden_gate/`

Vendored slice of [GGAssembler](https://github.com/Fleishman-Lab/GGAssembler):
`gate_data.py` plus Potapov / Pryor CSV matrices. Do not “clean up” this
tree without re-checking AGPL and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
`ligation_fidelity.py` imports `dawdlib_golden_gate` from the installed
package, falling back to `third_party/` in a raw checkout.

### Root notebooks

`grasp_oneshot_designer.ipynb` and `grasp_library_designer.ipynb` are what
the README Colab badges open. Cell **0 · Install** currently pins PyPI
(`grasp-library-designer>=0.1.11`) and purges stale `grasp_library` modules
from `sys.modules`. Keep that purge if you change the install cell.

---

## 6. Local development

```bash
git clone https://github.com/JustABiologist/grasp-library-designer.git
cd grasp-library-designer
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[notebook,dev]"
pytest
```

Python ≥ 3.10. There is no CLI yet; call `run_oneshot_design` /
`run_library_redesign_and_anneal` from a notebook or a short script.

Smoke a one-shot without Colab:

```python
from pathlib import Path
from grasp_library import (
    materialize_project,
    build_default_config,
    run_oneshot_design,
)
from grasp_library.codon_tables import load_codon_usage

project = materialize_project()
config = build_default_config(project / "input")
config["optimizer"]["iterations_per_part"] = 0  # greedy; fast
_, codon_data = load_codon_usage(project / "input" / "codon_usage.csv")
run_oneshot_design(
    target_rna="UUACACGUG",
    codon_data=codon_data,
    config=config,
    output_dir=Path("grasp_library_project/output/dev_oneshot"),
)
```

---

## 7. Suggested refactor — small PRs, one concern each

Do not boil the ocean. The domain code is the asset; the structure around
it is what hurts. Each item below is a **separate branch and PR**. Merge
order is the numbered list.

### Ground rules for every change

- One concern per PR. If a review needs a novel, split it.
- Do not bump `pyproject.toml` version on refactor / CI / docs branches.
  Version bumps are release PRs only (see §9).
- Do not mix notebook Form copy with ligation math.
- If you touch overhangs, masks, or ligation tables, add or extend a test
  in the matching `tests/test_*.py` file. Do not “fix” a failing strand
  test by loosening the assertion.
- Keep AGPL attribution when moving vendored files.

### Branch naming

```
feat/<short-topic>        new user-visible behavior
fix/<short-topic>         correctness bug
refactor/<short-topic>    structure, no intended behavior change
ci/<short-topic>          GitHub Actions, hooks, linters
docs/<short-topic>        README / DEV_DOC / comments
release/v0.1.x            version bump + notebook pin only
```

Examples: `ci/pytest-gha`, `refactor/notebook-extra`, `fix/level0-rc-overhang`.

Work from `main`, keep branches short-lived, open a PR even if you are the
only reviewer. That is how CI (once it exists) actually runs.

### PR 1 — `ci/pytest-gha` (do this first)

Add `.github/workflows/test.yml`:

- Trigger: pull requests and pushes to `main`
- Matrix: Python 3.10, 3.11, 3.12, 3.13
- `pip install -e ".[dev]"` then `pytest -q`

Until this lands, every later refactor is unguarded. Do not add ruff/mypy
in the same PR.

### PR 2 — `ci/ruff` (optional, immediately after)

Add ruff with a loose first config (line length 100, no isort war). Fix
only what the first run flags in files you already touch, or land a
dedicated format PR. Do not reformat `third_party/`.

### PR 3 — `refactor/notebook-extra`

Move `ipywidgets`, `IPython`, and `matplotlib` out of core
`[project].dependencies` into `[project.optional-dependencies] notebook`.
Keep `pandas`, `numpy`, `biopython`, `pyyaml`, `openpyxl` in core.

Guard widget imports so `from grasp_library import run_oneshot_design`
does not require a Jupyter kernel. Update Colab install cells to:

```python
%pip install -q -U --force-reinstall "grasp-library-designer[notebook] @ git+..."
```

(or the PyPI equivalent after the next release).

This is the highest-leverage packaging change. It will break Colab if the
notebook extra is forgotten — test with §8 before merging.

### PR 4 — `refactor/single-genbank-source`

Stop tracking `grasp_library_project/profiles/grasp_nar2025/genbank/*.gb`
once `materialize_project()` can always copy from package data. Keep the
directory in `.gitignore` or a one-line README pointer. Same for root vs
`grasp_library/notebooks/` duplication: pick a generation step (e.g. a
small `scripts/sync_notebooks.py` run in CI) so humans edit one place.

### PR 5 — `feat/cli`

`workflows.py` already claims “notebook or CLI”. Add

```toml
[project.scripts]
grasp-library = "grasp_library.cli:main"
```

Thin wrapper: `oneshot` and `library` subcommands, argparse, write to
`--output`. No widgets. Tests: invoke help + one greedy oneshot on a tmp
path.

### PR 6 — `refactor/typed-config`

Introduce a dataclass (or pydantic, if you accept the extra dep) for the
config that `build_default_config` already builds. Keep a
`to_dict()` / `from_dict()` so notebooks do not break. Convert one
consumer per follow-up PR (`oneshot`, then `optimizer`, then widgets).
Do not convert everything in one diff.

### PR 7 — split the 750–950 line modules

Only after CI exists. Suggested seams:

| Module | Split toward |
|---|---|
| `import_grasp.py` | GenBank parse / tables / GAP compile / order-fragment builder |
| `control_panel.py` | widget construction vs `build_default_config` (config can move next to the typed model) |
| `workflows.py` | redesign vs anneal vs export vs compile (already almost that) |
| `optimizer.py` | SA loop vs `synthesis_qc` vs `simulate_assembled_cds` |

Each split is its own PR. Public imports in `__init__.py` must keep working
(`from grasp_library import optimize_library` etc.).

### Later, not now

- Shrink `__all__` (breaking)
- Replace simulated annealing with an external codon optimizer
- Sphinx / mkdocs API docs
- mypy in CI (wait until config is typed, or you will drown in `dict[str, Any]`)

---

## 8. Test a new version on Colab **without** a PyPI bump

PyPI is only for tagged releases. Feature branches, ligation-matrix
tweaks, and notebook copy should be tried on Colab from git (or a wheel)
first.

Colab badges on the README always open **`main`**. To test a branch you
must both (a) open the notebook from that branch if Forms cells changed,
and (b) install the package from that branch.

### 8.1 Fastest: install the branch into a notebook already on Colab

1. Open either notebook (main or branch URL).
2. **Do not run** the stock **0 · Install** cell.
3. Insert a cell **above** it (or replace it) with one of the blocks below.
4. Runtime → Restart session after a successful install if imports look stale.
5. Run Settings and the rest of the notebook as usual.

**Public branch** (no token):

```python
#@title 0 · Install (git branch, not PyPI)
BRANCH = "feat/my-topic"  # or a commit SHA

%pip install -q -U --force-reinstall \
  "grasp-library-designer @ git+https://github.com/JustABiologist/grasp-library-designer.git@{BRANCH}"

import sys
for _m in [m for m in list(sys.modules)
           if m == "grasp_library" or m.startswith("grasp_library.")]:
    del sys.modules[_m]

from importlib.metadata import version
import grasp_library
print("grasp-library-designer", version("grasp-library-designer"))
print("file", grasp_library.__file__)
```

Replace the GitHub org/repo if you are on a fork:

```text
git+https://github.com/<user>/grasp-library-designer.git@<branch>
```

**Pinned commit** (reproducible review):

```python
%pip install -q -U --force-reinstall \
  "grasp-library-designer @ git+https://github.com/JustABiologist/grasp-library-designer.git@abc1234"
```

Confirm you are not on the PyPI wheel:

```python
import grasp_library, pathlib
print(pathlib.Path(grasp_library.__file__).resolve())
# expect something under /usr/local/lib/python3.*/dist-packages/grasp_library/
# and version() may still print 0.1.11 — that is the pyproject version on the
# branch, not proof of PyPI. The path + git URL in the pip log are the proof.
```

### 8.2 Open the **branch** notebook (needed when Forms cells changed)

```
https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/<BRANCH>/grasp_oneshot_designer.ipynb
https://colab.research.google.com/github/JustABiologist/grasp-library-designer/blob/<BRANCH>/grasp_library_designer.ipynb
```

Then still point **0 · Install** at `git+...@<BRANCH>`. Otherwise Colab
runs new Forms against the old PyPI package.

### 8.3 Blank Colab: install git, then drop the packaged notebook

Useful when you only changed library code, not Forms:

```python
%pip install -q -U --force-reinstall \
  "grasp-library-designer[notebook] @ git+https://github.com/JustABiologist/grasp-library-designer.git@BRANCH"

from grasp_library import write_notebook
write_notebook("oneshot")   # or "library"
```

Open the written `.ipynb` from the Files panel. **Skip its Install cell**
(it still says PyPI) or you will overwrite the git install.

Until PR 3 lands, drop `[notebook]` from the pip spec if the extra does
not exist yet.

### 8.4 Wheel upload (no GitHub from Colab)

On your laptop:

```bash
pip install -e ".[dev]"
python -m build
# dist/grasp_library_designer-0.1.11-py3-none-any.whl
```

In Colab: upload the wheel →

```python
%pip install -q -U --force-reinstall /content/grasp_library_designer-0.1.11-py3-none-any.whl
```

Bump the *filename* if you need to distinguish builds, but you do **not**
need to bump the project version for this. Two wheels with the same
version are fine as long as you `--force-reinstall` the file you uploaded.

### 8.5 Private fork

Colab cannot clone a private repo without a token. Prefer §8.4 (wheel
upload) or a short-lived GitHub PAT in the pip URL. Do not commit tokens
into notebooks.

### 8.6 Stale-import checklist

Colab runtimes cache `sys.modules`. After any reinstall:

1. The `del sys.modules[...]` loop from the install cell (keep it).
2. Runtime → Restart session if `hasattr` checks fail (the notebooks
   already warn when `kazusa_codon_reminder` is missing).
3. Re-run Install, then Settings, then everything below. Do not resume
   mid-notebook across a package swap.

### 8.7 What not to do

- Do not `pip install -U grasp-library-designer` while testing a branch.
  That pulls PyPI and hides your git install.
- Do not bump the version in `pyproject.toml` “so Colab notices”.
  `--force-reinstall` from git/wheel is enough.
- Do not publish to PyPI from a feature branch.

---

## 9. Releases (when you *do* bump PyPI)

Separate `release/v0.1.x` PR, after the feature is on `main` and Colab-tested
via §8:

1. Set `version` in `pyproject.toml` and the fallback in
   `grasp_library/__init__.py`.
2. Update the notebook Install pin (`>=0.1.x`) in **both** root notebooks
   and `grasp_library/notebooks/`.
3. Update the README pip snippet.
4. Merge, tag `v0.1.x`, build, `twine upload`.

If a Colab user still sees old behavior: they need `--force-reinstall` and
a session restart. Version pins like `>=0.1.11` will not reinstall 0.1.11
if that version is already in the runtime.

---

## 10. Conventions when editing domain code

**Overhangs.** Always write physical sticky ends 5′→3′. Compatible ends are
reverse complements. If the bases retained on the coding strand differ from
the sticky-end label, use `*_assembled_coding_site` — do not overload one
field.

**Ligation.** Default scoring is stage-matched Pryor cycling: BsaI-HFv2 for
Levels −1 and 1, BbsI-HF for Level 0. Potapov ligase-only tables are Level 0
overrides. Do not interpolate a 16 °C matrix that was not measured.

**QC.** `synthesis_qc` is a heuristic. Keep
`vendor_acceptance_confirmed=False`. Warnings are not failures;
`hard_constraints_passed` can be true when `passed` is false.

**Masks.** Coding-mask compatibility beats codon frequency
(`build_allowed_codons`). Optimization must still
`cds_matches_organism` and `mask_matches`.

**ARELF.** Library redesign explores motif-relative offsets 0–11, not only
the paper’s native cut indices. A candidate is `(overhang, offset)`.

---

## 11. “I want to change X — which files?”

| Change | Start here | Tests |
|---|---|---|
| Colab Form field | `colab_forms.py`, both notebooks (root + packaged) | `test_dashboard_interfaces.py` |
| Widget dashboard | `control_panel.py`, `notebook_ui.py` | `test_dashboard_interfaces.py` |
| Default overhangs / presets | `assembly_interfaces.py`, `build_default_config` | `test_assembly_interfaces.py`, `test_strand_semantics.py` |
| Order-fragment arms | `assembly_interfaces.py`, `import_grasp.py` | `test_oneshot.py` |
| One-shot pipeline | `oneshot.py` | `test_oneshot.py`, `test_architectures.py` |
| Library pipeline | `workflows.py`, `pareto.py`, `arelf.py` | `test_arelf_movable.py` |
| Codon SA / QC | `optimizer.py`, `objectives.py` | `test_qc_fidelity.py` |
| Ligation matrix / enzyme | `synthesis_vendors.py`, `ligation_fidelity.py` | `test_ligation_protocols.py` |
| GAP part picking | `import_grasp.py` | `test_architectures.py` |
| PPR code / RNA length | `binder.py` | oneshot tests that pass 9 / 14 / 19 nt |
| Packaged GenBank | `grasp_library/data/...` then `python -m grasp_library.import_grasp` | oneshot + strand tests |
| PyPI / extras | `pyproject.toml` | install + Colab §8 |
