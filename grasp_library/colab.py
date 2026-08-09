"""Helpers to drop Colab notebooks into the current working directory."""

from __future__ import annotations

import shutil
from importlib import resources
from pathlib import Path
from typing import Literal

NotebookName = Literal["oneshot", "library"]

_NOTEBOOK_FILES = {
    "oneshot": "grasp_oneshot_designer.ipynb",
    "library": "grasp_library_designer.ipynb",
}


def write_notebook(
    which: NotebookName = "oneshot",
    dest_dir: str | Path = ".",
    *,
    overwrite: bool = True,
) -> Path:
    """
    Copy a bundled Colab notebook next to the current working directory.

    Useful in a blank Colab after ``pip install grasp-library-designer``:
    then open the written ``.ipynb`` from the Files panel.
    """
    filename = _NOTEBOOK_FILES[which]
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / filename

    # Prefer package data (installed wheel); fall back to repo checkout.
    try:
        pkg = resources.files("grasp_library") / "notebooks" / filename
        with resources.as_file(pkg) as src:
            if not overwrite and dest.exists():
                return dest.resolve()
            shutil.copy2(src, dest)
            return dest.resolve()
    except (FileNotFoundError, ModuleNotFoundError, TypeError):
        pass

    repo_candidate = Path(__file__).resolve().parents[1] / filename
    if not repo_candidate.exists():
        raise FileNotFoundError(
            f"Notebook {filename} not found in package data or repo root."
        )
    if not overwrite and dest.exists():
        return dest.resolve()
    shutil.copy2(repo_candidate, dest)
    return dest.resolve()
