"""Locate bundled data and materialize a writable project folder."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]

PACKAGE_DIR = Path(__file__).resolve().parent
DATA_DIR = PACKAGE_DIR / "data"
BUNDLED_PROFILE_GENBANK = DATA_DIR / "profiles" / "grasp_nar2025" / "genbank"
DEFAULT_PROJECT_NAME = "grasp_library_project"


def bundled_profile_genbank() -> Path:
    """GenBank modules shipped inside the installed package."""
    return BUNDLED_PROFILE_GENBANK


def materialize_project(
    root: Optional[PathLike] = None,
    *,
    force_profile: bool = False,
) -> Path:
    """
    Create a writable project tree and copy bundled GRASP GenBank if needed.

    Layout:
      {root}/input/
      {root}/output/
      {root}/profiles/grasp_nar2025/genbank/
    """
    project = Path(root) if root is not None else Path.cwd() / DEFAULT_PROJECT_NAME
    input_dir = project / "input"
    output_dir = project / "output"
    profile_gb = project / "profiles" / "grasp_nar2025" / "genbank"
    for path in (input_dir, output_dir, profile_gb):
        path.mkdir(parents=True, exist_ok=True)

    src = bundled_profile_genbank()
    if src.is_dir():
        for gb in sorted(src.glob("*.gb")):
            dest = profile_gb / gb.name
            if force_profile or not dest.exists():
                shutil.copy2(gb, dest)

    return project.resolve()


def project_paths(root: Optional[PathLike] = None) -> dict[str, Path]:
    """Return standard project subpaths (materializes directories + GenBank)."""
    project = materialize_project(root)
    return {
        "project": project,
        "input": project / "input",
        "output": project / "output",
        "profile_genbank": project / "profiles" / "grasp_nar2025" / "genbank",
    }
