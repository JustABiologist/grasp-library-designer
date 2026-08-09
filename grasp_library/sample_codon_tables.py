"""Built-in codon-usage samples (Kazusa CUTG, frequency per thousand)."""

from __future__ import annotations

import re
from typing import Dict

import pandas as pd

# Frequency per thousand; U/T interchangeable. Source: Kazusa Codon Usage Database.
_ECOLI_KAZUSA = """
UUU 24.4  UCU 13.1  UAU 21.6  UGU  5.9
UUC 13.9  UCC  9.7  UAC 11.7  UGC  5.5
UUA 17.4  UCA 13.1  UAA  2.0  UGA  1.1
UUG 12.9  UCG  8.2  UAG  0.3  UGG 13.4
CUU 14.5  CCU  9.5  CAU 12.4  CGU 15.9
CUC  9.5  CCC  6.2  CAC  7.3  CGC 14.0
CUA  5.6  CCA  9.1  CAA 14.4  CGA  4.8
CUG 37.4  CCG 14.5  CAG 26.7  CGG  7.9
AUU 29.6  ACU 13.1  AAU 29.3  AGU 13.2
AUC 19.4  ACC 18.9  AAC 20.3  AGC 14.3
AUA 13.3  ACA 15.1  AAA 37.2  AGA  7.1
AUG 23.7  ACG 13.6  AAG 15.3  AGG  4.0
GUU 21.6  GCU 18.9  GAU 33.7  GGU 23.7
GUC 13.1  GCC 21.6  GAC 17.9  GGC 20.6
GUA 13.1  GCA 23.0  GAA 35.1  GGA 13.6
GUG 19.9  GCG 21.1  GAG 19.4  GGG 12.3
"""

_SCEREVISIAE_KAZUSA = """
UUU 26.1  UCU 23.5  UAU 18.8  UGU  8.1
UUC 18.4  UCC 14.2  UAC 14.8  UGC  4.8
UUA 26.2  UCA 18.7  UAA  1.1  UGA  0.7
UUG 27.2  UCG  8.6  UAG  0.5  UGG 10.4
CUU 12.3  CCU 13.5  CAU 13.6  CGU  6.4
CUC  5.4  CCC  6.8  CAC  7.8  CGC  2.6
CUA 13.4  CCA 18.3  CAA 27.3  CGA  3.0
CUG 10.5  CCG  5.3  CAG 12.1  CGG  1.7
AUU 30.1  ACU 20.3  AAU 35.7  AGU 14.2
AUC 17.2  ACC 12.7  AAC 24.8  AGC  9.8
AUA 17.8  ACA 17.8  AAA 41.9  AGA 21.3
AUG 20.9  ACG  8.0  AAG 30.8  AGG  9.2
GUU 22.1  GCU 21.2  GAU 37.6  GGU 23.9
GUC 11.8  GCC 12.6  GAC 20.2  GGC  9.8
GUA 11.8  GCA 16.2  GAA 45.6  GGA 10.9
GUG 10.8  GCG  6.2  GAG 19.2  GGG  6.0
"""

_HSAPIENS_KAZUSA = """
UUU 17.6  UCU 15.2  UAU 12.2  UGU 10.6
UUC 20.3  UCC 17.7  UAC 15.3  UGC 12.6
UUA  7.7  UCA 12.2  UAA  1.0  UGA  1.6
UUG 12.9  UCG  4.4  UAG  0.8  UGG 13.2
CUU 13.2  CCU 17.5  CAU 10.9  CGU  4.5
CUC 19.6  CCC 19.8  CAC 15.1  CGC 10.4
CUA  7.2  CCA 16.9  CAA 12.3  CGA  6.2
CUG 39.6  CCG  6.9  CAG 34.2  CGG 11.4
AUU 16.0  ACU 13.1  AAU 17.0  AGU 12.1
AUC 20.8  ACC 18.9  AAC 19.1  AGC 19.5
AUA  7.5  ACA 15.1  AAA 24.4  AGA 12.2
AUG 22.0  ACG  6.1  AAG 31.9  AGG 12.0
GUU 11.0  GCU 18.4  GAU 21.8  GGU 10.8
GUC 14.5  GCC 27.7  GAC 25.1  GGC 22.2
GUA  7.1  GCA 15.8  GAA 29.0  GGA 16.5
GUG 28.1  GCG  7.4  GAG 39.6  GGG 16.5
"""

# Kazusa CUTG nuclear CDS sums (frequency per thousand)
_EUGLENA_GRACILIS_KAZUSA = """
UUU  8.9  UCU 12.2  UAU  9.0  UGU  3.6
UUC 28.1  UCC 20.7  UAC 22.1  UGC 14.3
UUA  1.0  UCA  5.9  UAA  0.7  UGA  0.7
UUG 19.1  UCG  6.2  UAG  0.3  UGG 11.4
CUU  9.3  CCU 11.6  CAU  6.2  CGU  9.1
CUC 13.7  CCC 20.8  CAC 13.9  CGC 18.3
CUA  0.8  CCA 10.9  CAA  8.9  CGA  5.2
CUG 37.6  CCG  9.8  CAG 32.3  CGG 10.9
AUU 17.5  ACU 12.7  AAU  8.8  AGU  4.9
AUC 31.7  ACC 25.2  AAC 28.0  AGC 13.8
AUA  0.7  ACA 10.6  AAA  7.5  AGA  2.5
AUG 27.6  ACG 10.3  AAG 45.9  AGG  4.3
GUU 19.1  GCU 22.9  GAU 21.4  GGU 18.3
GUC 22.3  GCC 37.0  GAC 33.4  GGC 29.6
GUA  1.8  GCA 16.7  GAA 13.4  GGA 12.8
GUG 38.5  GCG 14.5  GAG 46.5  GGG 16.1
"""

_CHLAMYDOMONAS_REINHARDTII_KAZUSA = """
UUU  5.0  UCU  4.7  UAU  2.6  UGU  1.4
UUC 27.1  UCC 16.1  UAC 22.8  UGC 13.1
UUA  0.6  UCA  3.2  UAA  1.0  UGA  0.5
UUG  4.0  UCG 16.1  UAG  0.4  UGG 13.2
CUU  4.4  CCU  8.1  CAU  2.2  CGU  4.9
CUC 13.0  CCC 29.5  CAC 17.2  CGC 34.9
CUA  2.6  CCA  5.1  CAA  4.2  CGA  2.0
CUG 65.2  CCG 20.7  CAG 36.3  CGG 11.2
AUU  8.0  ACU  5.2  AAU  2.8  AGU  2.6
AUC 26.6  ACC 27.7  AAC 28.5  AGC 22.8
AUA  1.1  ACA  4.1  AAA  2.4  AGA  0.7
AUG 25.7  ACG 15.9  AAG 43.3  AGG  2.7
GUU  5.1  GCU 16.7  GAU  6.7  GGU  9.5
GUC 15.4  GCC 54.6  GAC 41.7  GGC 62.0
GUA  2.0  GCA 10.6  GAA  2.8  GGA  5.0
GUG 46.5  GCG 44.4  GAG 53.5  GGG  9.7
"""

# Chloroplast of C. reinhardtii (bacterial-like; NCBI genetic code 11)
_CHLAMY_CHLOROPLAST_KAZUSA = """
UUU 33.4  UCU 17.0  UAU 24.6  UGU  7.6
UUC 17.1  UCC  2.8  UAC 10.0  UGC  1.5
UUA 77.7  UCA 22.0  UAA  2.9  UGA  0.1
UUG  4.3  UCG  4.0  UAG  0.4  UGG 13.5
CUU 14.3  CCU 15.5  CAU 10.1  CGU 32.4
CUC  1.0  CCC  3.4  CAC  8.8  CGC  4.1
CUA  6.4  CCA 23.6  CAA 38.4  CGA  3.4
CUG  3.7  CCG  2.4  CAG  4.1  CGG  0.5
AUU 51.4  ACU 24.4  AAU 42.1  AGU 16.0
AUC  8.2  ACC  5.1  AAC 17.7  AGC  5.4
AUA  6.9  ACA 32.4  AAA 69.1  AGA  5.3
AUG 22.3  ACG  3.9  AAG  6.2  AGG  0.9
GUU 29.3  GCU 34.0  GAU 25.3  GGU 44.0
GUC  2.5  GCC  5.9  GAC  9.8  GGC  6.4
GUA 26.0  GCA 20.7  GAA 41.1  GGA  8.6
GUG  5.6  GCG  3.3  GAG  5.7  GGG  3.7
"""


def parse_frequency_block(block: str) -> Dict[str, float]:
    """Parse a Kazusa-style `UUU 24.4 …` frequency-per-thousand block."""
    pairs = re.findall(r"([ACGTU]{3})\s+([0-9]+(?:\.[0-9]+)?)", block.upper())
    if len(pairs) != 64:
        raise ValueError(f"Expected 64 codons, found {len(pairs)}")
    return {codon.replace("U", "T"): float(freq) for codon, freq in pairs}


# Back-compat alias
_parse_frequency_block = parse_frequency_block


SAMPLE_CODON_TABLES: Dict[str, dict] = {
    "Escherichia coli (Kazusa)": {
        "organism": "Escherichia coli",
        "clade": "Bacteria",
        "genetic_code": 1,
        "source": "Kazusa CUTG / GenBank bacterial CDS (species 37762)",
        "url": "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37762",
        "frequencies": _parse_frequency_block(_ECOLI_KAZUSA),
    },
    "Saccharomyces cerevisiae (Kazusa)": {
        "organism": "Saccharomyces cerevisiae",
        "clade": "Fungi",
        "genetic_code": 1,
        "source": "Kazusa CUTG / GenBank (species 4932)",
        "url": "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932",
        "frequencies": _parse_frequency_block(_SCEREVISIAE_KAZUSA),
    },
    "Homo sapiens (Kazusa)": {
        "organism": "Homo sapiens",
        "clade": "Metazoa",
        "genetic_code": 1,
        "source": "Kazusa CUTG approximate human nuclear codon usage",
        "url": "https://www.kazusa.or.jp/codon/",
        "frequencies": _parse_frequency_block(_HSAPIENS_KAZUSA),
    },
    "Euglena gracilis nuclear (Kazusa)": {
        "organism": "Euglena gracilis",
        "clade": "Euglenozoa / Excavata",
        "genetic_code": 1,
        "source": "Kazusa CUTG nuclear CDS sum (species 3039; 74 CDS / 43653 codons)",
        "url": "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3039",
        "frequencies": _parse_frequency_block(_EUGLENA_GRACILIS_KAZUSA),
    },
    "Chlamydomonas reinhardtii nuclear (Kazusa)": {
        "organism": "Chlamydomonas reinhardtii",
        "clade": "Chlorophyta / Viridiplantae",
        "genetic_code": 1,
        "source": "Kazusa CUTG nuclear CDS sum (species 3055; 846 CDS / 420455 codons)",
        "url": "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3055",
        "frequencies": _parse_frequency_block(_CHLAMYDOMONAS_REINHARDTII_KAZUSA),
    },
    "Chlamydomonas reinhardtii chloroplast (Kazusa)": {
        "organism": "Chlamydomonas reinhardtii chloroplast",
        "clade": "Chlorophyta chloroplast",
        "genetic_code": 11,
        "source": "Kazusa CUTG chloroplast CDS sum (species 3055.chloroplast; 93 CDS)",
        "url": "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3055.chloroplast",
        "frequencies": _parse_frequency_block(_CHLAMY_CHLOROPLAST_KAZUSA),
        "notes": (
            "Use only for proteins expressed in the chloroplast genome. "
            "Nuclear-encoded chloroplast-targeted proteins should use the nuclear table. "
            "NCBI genetic code 11 (bacterial)."
        ),
    },
    "Custom file (codon_usage.csv)": {
        "organism": "custom",
        "clade": "custom",
        "genetic_code": 1,
        "source": "User-supplied grasp_library_project/input/codon_usage.csv",
        "url": "",
        "frequencies": None,
    },
    "Fetch from Kazusa by species ID": {
        "organism": "kazusa_fetch",
        "clade": "Kazusa CUTG",
        "genetic_code": 1,
        "source": "Fetched live from https://www.kazusa.or.jp/codon/",
        "url": "https://www.kazusa.or.jp/codon/",
        "frequencies": None,
    },
    "Upload your own codon table": {
        "organism": "upload",
        "clade": "custom",
        "genetic_code": 1,
        "source": "User-uploaded codon table (CSV or Kazusa text)",
        "url": "",
        "frequencies": None,
    },
}

# Special organism keys used by forms / control panel
FETCH_FROM_KAZUSA = "Fetch from Kazusa by species ID"
UPLOAD_OWN_TABLE = "Upload your own codon table"
CUSTOM_FILE = "Custom file (codon_usage.csv)"

KAZUSA_REMINDER_HTML = (
    "Browse codon tables at "
    '<a href="https://www.kazusa.or.jp/codon/" target="_blank" rel="noopener">'
    "Kazusa Codon Usage Database</a> "
    "(search organism → note the <code>species=</code> accession in the URL). "
    "If it is not in the built-in list, choose <b>Fetch from Kazusa by species ID</b> "
    "or <b>Upload your own codon table</b> (CSV with "
    "<code>codon,frequency</code> or a Kazusa frequency block)."
)


def sample_names() -> list[str]:
    return list(SAMPLE_CODON_TABLES.keys())


def builtin_sample_names() -> list[str]:
    """Built-in tables only (excludes custom / fetch / upload sentinels)."""
    return [
        name
        for name, meta in SAMPLE_CODON_TABLES.items()
        if meta.get("frequencies") is not None
    ]


def codon_table_dataframe(name: str) -> pd.DataFrame:
    from Bio.Seq import Seq

    entry = SAMPLE_CODON_TABLES[name]
    frequencies = entry["frequencies"]
    if frequencies is None:
        raise ValueError(
            "Custom file selected — load codon_usage.csv instead of a built-in sample."
        )
    code = int(entry.get("genetic_code", 1))
    rows = []
    for codon, frequency in sorted(frequencies.items()):
        aa = str(Seq(codon).translate(table=code))
        rows.append({"codon": codon, "aa": aa, "frequency": frequency})
    return pd.DataFrame(rows)


def write_sample_codon_table(name: str, path) -> pd.DataFrame:
    table = codon_table_dataframe(name)
    table.to_csv(path, index=False)
    return table
