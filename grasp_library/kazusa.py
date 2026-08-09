"""Fetch codon-usage tables from the Kazusa CUTG webserver."""

from __future__ import annotations

import re
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, Optional, Tuple

from .sample_codon_tables import parse_frequency_block

KAZUSA_HOME = "https://www.kazusa.or.jp/codon/"
KAZUSA_SHOW = "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi"

# CUTG species accessions: digits, optionally with organelle suffix
# e.g. 37762, 3055.chloroplast, 4932
_SPECIES_RE = re.compile(r"^\d+(?:\.[A-Za-z0-9_-]+)?$")

_USER_AGENT = "grasp-library-designer/1.0 (+https://github.com/grasp-library-designer)"


def normalize_kazusa_species_id(value: str) -> str:
    """Normalize a Kazusa species accession from a typed id or pasted URL."""
    text = str(value or "").strip()
    if not text:
        raise ValueError("Kazusa species ID is empty.")

    if "showcodon.cgi" in text or "species=" in text.lower():
        match = re.search(r"[?&]species=([^&\s#]+)", text, flags=re.I)
        if not match:
            raise ValueError(
                "Could not find species=… in the Kazusa URL. "
                f"Paste a link like {KAZUSA_SHOW}?species=37762"
            )
        text = urllib.parse.unquote(match.group(1))

    text = text.strip().strip("/")
    if not _SPECIES_RE.match(text):
        raise ValueError(
            f"Invalid Kazusa species ID {text!r}. "
            "Use the numeric CUTG accession from Kazusa "
            f"({KAZUSA_HOME}), e.g. 37762 or 3055.chloroplast."
        )
    return text


def kazusa_table_url(species_id: str, *, with_amino_acids: bool = False) -> str:
    sid = normalize_kazusa_species_id(species_id)
    url = (
        f"{KAZUSA_SHOW}?species={urllib.parse.quote(sid, safe='.')}&style=N"
    )
    if with_amino_acids:
        url += "&aa=1"
    return url


def _http_get(url: str, *, timeout: float = 30.0) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw = response.read()
    except urllib.error.HTTPError as exc:
        raise ValueError(
            f"Kazusa returned HTTP {exc.code} for {url}. "
            "Check the species ID on the Codon Usage Database."
        ) from exc
    except urllib.error.URLError as exc:
        raise ValueError(
            f"Could not reach Kazusa ({exc.reason}). "
            "Check your network connection and try again."
        ) from exc

    for encoding in ("utf-8", "latin-1", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def _strip_tags(html: str) -> str:
    return re.sub(r"<[^>]+>", "", html)


def _extract_organism_name(html: str) -> str:
    strong = re.search(
        r"<STRONG>\s*(?:<i>)?(.*?)(?:</i>)?\s*\[",
        html,
        flags=re.I | re.S,
    )
    if strong:
        name = " ".join(_strip_tags(strong.group(1)).split())
        if name:
            return name
    return "Kazusa organism"


def _extract_genetic_code(html: str) -> Optional[int]:
    match = re.search(r"Genetic code\s+(\d+)\s*:", html, flags=re.I)
    if match:
        return int(match.group(1))
    return None


def _extract_pre_block(html: str) -> str:
    match = re.search(r"<PRE[^>]*>(.*?)</PRE>", html, flags=re.I | re.S)
    if not match:
        raise ValueError(
            "Kazusa response did not contain a codon table (<PRE> block). "
            "The species ID may be wrong or the page format changed."
        )
    return match.group(1)


def _pre_to_frequency_block(pre: str) -> str:
    """Normalize Kazusa <PRE> text to `CODON freq` pairs for parsing."""
    text = pre.upper()
    # aa=1 style: UUU F 0.64 24.4 ( 56791)  → keep per-thousand frequency
    text = re.sub(
        r"([ACGTU]{3})\s+[A-Z*]\s+[0-9.]+\s+([0-9]+(?:\.[0-9]+)?)\s*\(",
        r"\1 \2 (",
        text,
    )
    text = re.sub(r"\([^)]*\)", " ", text)
    return text


def parse_kazusa_html(html: str) -> Tuple[Dict[str, float], dict]:
    """Parse a Kazusa showcodon.cgi HTML page into frequencies + metadata."""
    if re.search(r"not found|no data", html, flags=re.I) and "<PRE" not in html.upper():
        raise ValueError("Kazusa reported no codon-usage data for that species ID.")

    pre = _extract_pre_block(html)
    frequencies = parse_frequency_block(_pre_to_frequency_block(pre))

    meta = {
        "organism": _extract_organism_name(html),
        "clade": "Kazusa CUTG",
        "source": "Kazusa Codon Usage Database (CUTG)",
        "url": "",
        "frequencies": frequencies,
    }
    code = _extract_genetic_code(html)
    if code is not None:
        meta["genetic_code"] = code
    return frequencies, meta


def fetch_kazusa_codon_table(
    species_id: str,
    *,
    timeout: float = 30.0,
) -> Tuple[Dict[str, float], dict]:
    """Download a codon table by Kazusa species accession.

    Returns (frequencies, meta). Frequencies are per-thousand CUTG values
    keyed by DNA triplets (T, not U).
    """
    sid = normalize_kazusa_species_id(species_id)
    # Prefer aa=1 so NCBI genetic-code hints are present; parser handles both.
    url = kazusa_table_url(sid, with_amino_acids=True)
    html = _http_get(url, timeout=timeout)
    try:
        frequencies, meta = parse_kazusa_html(html)
    except ValueError:
        url = kazusa_table_url(sid, with_amino_acids=False)
        html = _http_get(url, timeout=timeout)
        frequencies, meta = parse_kazusa_html(html)

    meta["url"] = kazusa_table_url(sid, with_amino_acids=False)
    meta["species_id"] = sid
    meta["source"] = f"Kazusa CUTG (species {sid})"
    meta["genetic_code"] = int(meta.get("genetic_code", 1))
    return frequencies, meta
