"""Parse / save user-supplied codon tables (CSV or Kazusa text)."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import pandas as pd
from Bio.Seq import Seq

from .sample_codon_tables import parse_frequency_block

PathLike = Union[str, Path]
BytesLike = Union[bytes, bytearray, memoryview]


def _decode_upload(data: Union[str, BytesLike]) -> str:
    if isinstance(data, str):
        return data
    raw = bytes(data)
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def frequencies_to_dataframe(
    frequencies: Dict[str, float],
    *,
    genetic_code: int = 1,
) -> pd.DataFrame:
    rows = []
    for codon, frequency in sorted(frequencies.items()):
        aa = str(Seq(codon).translate(table=genetic_code))
        rows.append({"codon": codon, "aa": aa, "frequency": float(frequency)})
    return pd.DataFrame(rows)


def parse_codon_table_text(
    text: str,
    *,
    genetic_code: int = 1,
) -> Tuple[pd.DataFrame, dict]:
    """Parse CSV (`codon,frequency`) or a Kazusa frequency block / HTML dump."""
    content = text.strip()
    if not content:
        raise ValueError("Uploaded codon table is empty.")

    # HTML dump from Kazusa
    if "<PRE" in content.upper() or "showcodon" in content.lower():
        from .kazusa import parse_kazusa_html

        frequencies, meta = parse_kazusa_html(content)
        table = frequencies_to_dataframe(frequencies, genetic_code=genetic_code)
        return table, {
            **meta,
            "genetic_code": int(meta.get("genetic_code", genetic_code)),
            "source": meta.get("source", "Uploaded Kazusa HTML"),
        }

    # CSV / TSV with header
    first = content.splitlines()[0].lower()
    if "codon" in first and ("freq" in first or "frequency" in first or "count" in first):
        table = pd.read_csv(io.StringIO(content), sep=None, engine="python")
        table.columns = [str(c).strip().lower() for c in table.columns]
        # Accept common aliases
        rename = {}
        for col in table.columns:
            if col in {"freq", "frequency_per_thousand", "per_thousand", "count"}:
                rename[col] = "frequency"
            elif col in {"amino_acid", "aminoacid", "aa_code"}:
                rename[col] = "aa"
        if rename:
            table = table.rename(columns=rename)
        if "codon" not in table.columns or "frequency" not in table.columns:
            raise ValueError(
                "CSV needs columns `codon` and `frequency` "
                f"(found {list(table.columns)})."
            )
        table = table[["codon", "frequency"] + (["aa"] if "aa" in table.columns else [])].copy()
        table["codon"] = (
            table["codon"]
            .astype(str)
            .str.upper()
            .str.replace("U", "T", regex=False)
            .str.replace(r"\s+", "", regex=True)
        )
        table["frequency"] = pd.to_numeric(table["frequency"], errors="raise")
        if "aa" not in table.columns:
            table["aa"] = [
                str(Seq(c).translate(table=genetic_code)) for c in table["codon"]
            ]
        return table, {
            "organism": "uploaded CSV",
            "clade": "custom",
            "genetic_code": int(genetic_code),
            "source": "User-uploaded CSV codon table",
            "url": "",
        }

    # Bare Kazusa frequency block
    try:
        frequencies = parse_frequency_block(content)
    except ValueError as exc:
        raise ValueError(
            "Could not parse upload. Provide a CSV with columns "
            "`codon,frequency` or a Kazusa-style frequency block "
            "(64 triplets like `UUU 24.4 …`)."
        ) from exc

    table = frequencies_to_dataframe(frequencies, genetic_code=genetic_code)
    return table, {
        "organism": "uploaded Kazusa text",
        "clade": "custom",
        "genetic_code": int(genetic_code),
        "source": "User-uploaded Kazusa frequency block",
        "url": "",
    }


def save_codon_table(
    table: pd.DataFrame,
    path: PathLike,
) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in ("codon", "aa", "frequency") if c in table.columns]
    table[cols].to_csv(out, index=False)
    return out


def write_uploaded_codon_table(
    data: Union[str, BytesLike],
    path: PathLike,
    *,
    genetic_code: int = 1,
    filename: Optional[str] = None,
) -> Tuple[pd.DataFrame, dict]:
    """Parse upload bytes/text and write canonical codon_usage.csv."""
    table, meta = parse_codon_table_text(
        _decode_upload(data), genetic_code=genetic_code
    )
    save_codon_table(table, path)
    if filename:
        meta["filename"] = filename
    meta["path"] = str(path)
    return table, meta


def prompt_colab_codon_upload(
    path: PathLike,
    *,
    genetic_code: int = 1,
) -> Tuple[pd.DataFrame, dict]:
    """Open Colab's file picker and save the chosen codon table."""
    try:
        from google.colab import files  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Colab file upload is only available in Google Colab. "
            "Use the control-panel Upload widget, or place a file at "
            f"{path}."
        ) from exc

    print("Upload a codon table (CSV with codon,frequency — or Kazusa text/HTML).")
    uploaded = files.upload()
    if not uploaded:
        raise ValueError("No file was uploaded.")
    filename, raw = next(iter(uploaded.items()))
    return write_uploaded_codon_table(
        raw, path, genetic_code=genetic_code, filename=str(filename)
    )
