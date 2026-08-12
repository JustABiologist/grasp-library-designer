"""Biologist-facing ipywidgets control panel for the GRASP notebook."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import ipywidgets as widgets
import yaml
from IPython.display import clear_output, display

from . import notebook_ui as ui
from .codon_tables import apply_organism_codon_table, validate_parts_for_organism
from .codon_validation import format_issues
from .sample_codon_tables import (
    FETCH_FROM_KAZUSA,
    KAZUSA_REMINDER_HTML,
    SAMPLE_CODON_TABLES,
    UPLOAD_OWN_TABLE,
    sample_names,
)
from .synthesis_vendors import (
    apply_enzyme_to_config,
    apply_ligation_table_to_config,
    apply_vendor_to_config,
    enzyme_names,
    ligation_table_names,
    vendor_names,
)

GENETIC_CODE_OPTIONS = [
    (1, "1 · Standard"),
    (11, "11 · Bacterial / chloroplast"),
    (4, "4 · Mold / protozoan / Coelenterate mito."),
    (6, "6 · Ciliate macronuclear"),
]

# Widget chrome that renders in VS Code (Layout/Style APIs, not <style> tags)
_DD = dict(
    style={"description_width": "140px"},
    layout=widgets.Layout(width="460px"),
)
_DD_WIDE = dict(
    style={"description_width": "140px"},
    layout=widgets.Layout(width="520px"),
)


def _label(text: str) -> widgets.HTML:
    return widgets.HTML(
        value=(
            f"<div style='font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;"
            f"font-size:11px;letter-spacing:0.12em;text-transform:uppercase;"
            f"color:#0f6b4c;margin:10px 0 4px 0'>{text}</div>"
        )
    )


def build_default_config(input_dir: Path) -> Dict[str, Any]:
    return {
        "project_name": "GRASP_9S_organism_optimized",
        "genetic_code": 1,
        "architecture": "9S",
        "nterm_overhang": "AGGT",
        "assembly_interfaces": {
            "overhang_notation": "directional_terminal_5p",
            "level_minus1_entry": {
                "profile": "custom",
                "vector_name": "Custom Level -1 entry vector",
                "enzyme": "BsaI",
                "n_terminal_overhang": "AACA",
                "c_terminal_overhang": "GGAG",
            },
            "level0": {
                "acceptor_name": "Custom Level 0 acceptor",
                "release_enzyme": "BpiI / BbsI",
                "acceptor_n_terminal_overhang": "CTCA",
                "acceptor_c_terminal_overhang": "CGAG",
                "block_junctions": {
                    "cds1_to_cds14": {
                        "upstream_c": "GTGA",
                        "downstream_n": "TCAC",
                        "arelf_offset_nt": 4,
                    },
                    "cds14_to_cds19": {
                        "upstream_c": "CACG",
                        "downstream_n": "CGTG",
                        "arelf_offset_nt": 1,
                    },
                    "cds1_to_cds2": {
                        "upstream_c": "CTTC",
                        "downstream_n": "GAAG",
                        "arelf_offset_nt": 11,
                    },
                },
            },
            "level1": {
                "acceptor_name": "Custom Level 1 acceptor",
                "n_terminal_overhang": "GCCC",
                "c_terminal_overhang": "GCGA",
            },
        },
        "assembly_enzyme": "GRASP default · BsaI + BpiI + BsmBI",
        "synthesis_vendor": "Twist · Standard gene guidelines",
        "codon_usage_file": input_dir / "codon_usage.csv",
        "parts_file": input_dir / "parts.csv",
        "target_map_file": input_dir / "target_map.csv",
        "forbidden_sites": {
            "BsaI": "GGTCTC",
            "BpiI": "GAAGAC",
            "BsmBI": "CGTCTC",
        },
        "synthesis": {
            "global_gc_min": 0.25,
            "global_gc_max": 0.65,
            "window_size": 50,
            "window_gc_min": 0.15,
            "window_gc_max": 0.85,
            "max_homopolymer": 3,
            "repeat_k": 16,
            "max_repeat_count": 1,
            "min_oligo_length": 20,
            "max_oligo_length": 300,
            "min_gene_length": 300,
            "max_gene_length": 5000,
        },
        "codon_optimization": {"minimum_relative_adaptiveness": 0.20},
        "weights": {
            "codon": 4.0,
            "global_gc": 3.0,
            "local_gc": 4.0,
            "homopolymer": 8.0,
            "internal_repeat": 5.0,
            "library_similarity": 3.0,
        },
        "optimizer": {
            "iterations_per_part": 2_000,  # Quick default; raise via control panel
            "initial_temperature": 2.0,
            "final_temperature": 0.02,
            "orthogonal_versions_per_part": 1,
        },
        "ligation": {
            # A cycling matrix is not a constant-temperature/time assay.
            "temperature": None,
            "hours": None,
            "min_efficiency": 0.25,
            "min_fidelity": 0.9,
            "table_name": (
                "GRASP Level 0 proxy · BbsI-HF + T4 · "
                "37↔16 °C cycling (Pryor 2020)"
            ),
            "ligation_table": "BbsI-HF.csv",
        },
        "overhang_redesign": {
            "enabled": False,
            "selection": "knee",
            "cut_mode": "movable_arelf",
            "allowed_arelf_offsets_nt": list(range(12)),
        },
        # Pareto CDS scoring uses greedy (0 iters). Keep evaluations modest.
        "pareto": {
            "max_evaluations": 40,
            "beam_width": 40,
            "iterations_per_part": 0,
            "junction_flank": 9,
        },
        "target_rna": "UUACACGUG",
    }


class GraspControlPanel:
    """Interactive dropdown panel that mutates a shared CONFIG dict."""

    def __init__(
        self,
        config: Dict[str, Any],
        *,
        input_dir: Path,
        load_codon_usage: Optional[Callable] = None,
        on_apply: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> None:
        self.config = config
        self.input_dir = Path(input_dir)
        self.load_codon_usage = load_codon_usage
        self.on_apply = on_apply
        self.CODON_TABLE = None
        self.CODON_DATA: Dict = {}
        self.SELECTED_ORGANISM = "Escherichia coli (Kazusa)"
        self._build_widgets()

    def _build_widgets(self) -> None:
        assembly = self.config.get("assembly_interfaces", {})
        entry = assembly.get("level_minus1_entry", {})
        level0 = assembly.get("level0", {})
        level1 = assembly.get("level1", {})
        block_junctions = level0.get("block_junctions", {})

        def _junction(name: str, key: str, default):
            return block_junctions.get(name, {}).get(key, default)

        self.organism = widgets.Dropdown(
            options=sample_names(),
            value="Escherichia coli (Kazusa)",
            description="Codon table",
            **_DD_WIDE,
        )
        self.kazusa_id = widgets.Text(
            value=str(self.config.get("kazusa_species_id", "")),
            description="Kazusa species",
            placeholder="e.g. 37762 or paste showcodon.cgi URL",
            **_DD_WIDE,
        )
        self.upload = widgets.FileUpload(
            accept=".csv,.txt,.html,.htm",
            multiple=False,
            description="Upload table",
            layout=widgets.Layout(width="520px"),
        )
        self.kazusa_help = widgets.HTML(
            value=(
                "<div style='font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;"
                "font-size:12px;color:#3d5248;line-height:1.45;margin:2px 0 8px 0;"
                "padding:8px 10px;background:#fff6e8;border-left:4px solid #e0b872'>"
                f"{KAZUSA_REMINDER_HTML}</div>"
            )
        )
        self.genetic_code = widgets.Dropdown(
            options=[(label, code) for code, label in GENETIC_CODE_OPTIONS],
            value=int(self.config.get("genetic_code", 1)),
            description="Genetic code",
            **_DD_WIDE,
        )
        self._sync_codon_source_visibility()
        self.vendor = widgets.Dropdown(
            options=vendor_names(),
            value=self.config.get(
                "synthesis_vendor", "Twist · Standard gene guidelines"
            ),
            description="Synthesis vendor",
            **_DD_WIDE,
        )
        self.enzyme = widgets.Dropdown(
            options=enzyme_names(),
            value=self.config.get(
                "assembly_enzyme", "GRASP default · BsaI + BpiI + BsmBI"
            ),
            description="Internal-site filter",
            **_DD_WIDE,
        )
        self.ligation = widgets.Dropdown(
            options=ligation_table_names(),
            value=self.config.get("ligation", {}).get(
                "table_name",
                "GRASP Level 0 proxy · BbsI-HF + T4 · "
                "37↔16 °C cycling (Pryor 2020)",
            ),
            description="Ligation table",
            **_DD_WIDE,
        )
        self.architecture = widgets.Dropdown(
            options=["9S", "14S", "19S"],
            value=self.config.get("architecture", "9S"),
            description="Architecture",
            **_DD,
        )
        self.nterm = widgets.Dropdown(
            options=["AGGT (MoClo N-fusion)", "AATG (Met start)"],
            value=(
                "AGGT (MoClo N-fusion)"
                if self.config.get("nterm_overhang", "AGGT") == "AGGT"
                else "AATG (Met start)"
            ),
            description="N-term overhang",
            **_DD,
        )
        self.redesign = widgets.Dropdown(
            options=[
                ("On — explore synonymous cuts within ARELF", True),
                ("Off — keep native GRASP overhangs", False),
            ],
            value=bool(self.config.get("overhang_redesign", {}).get("enabled", False)),
            description="Overhang redesign",
            **_DD_WIDE,
        )
        self.selection = widgets.Dropdown(
            options=["knee", "max_fidelity"],
            value=self.config.get("overhang_redesign", {}).get("selection", "knee"),
            description="Pick from Pareto",
            **_DD,
        )
        self.depth = widgets.Dropdown(
            options=[
                ("Quick preview (2k iters)", 2_000),
                ("Standard (25k iters)", 25_000),
                ("Thorough (80k iters)", 80_000),
            ],
            value=int(self.config.get("optimizer", {}).get("iterations_per_part", 2_000)),
            description="Optimize depth",
            **_DD,
        )
        self.target_rna = widgets.Text(
            value=str(self.config.get("target_rna", "UUACACGUG")),
            description="Target RNA",
            placeholder="e.g. UUACACGUG",
            **_DD,
        )
        profile_value = str(entry.get("profile", "custom"))
        if profile_value not in {"custom", "deposited_grasp"}:
            profile_value = "custom"
        self.interface_preset = widgets.Dropdown(
            options=[
                ("Custom editable interfaces (default)", "custom"),
                ("Deposited GRASP · pAGM1311 → pAGM9121", "deposited_grasp"),
            ],
            value=profile_value,
            description="Vector preset",
            **_DD_WIDE,
        )
        self.entry_vector_name = widgets.Text(
            value=str(entry.get("vector_name", "Custom Level -1 entry vector")),
            description="Level -1 vector",
            **_DD_WIDE,
        )
        self.entry_n = widgets.Text(
            value=str(entry.get("n_terminal_overhang", "AACA")),
            description="Entry N overhang",
            **_DD,
        )
        self.entry_c = widgets.Text(
            value=str(entry.get("c_terminal_overhang", "GGAG")),
            description="Entry C overhang",
            **_DD,
        )
        self.level0_acceptor_name = widgets.Text(
            value=str(level0.get("acceptor_name", "Custom Level 0 acceptor")),
            description="Level 0 acceptor",
            **_DD_WIDE,
        )
        self.level0_acceptor_n = widgets.Text(
            value=str(level0.get("acceptor_n_terminal_overhang", "CTCA")),
            description="L0 acceptor N",
            **_DD,
        )
        self.level0_acceptor_c = widgets.Text(
            value=str(level0.get("acceptor_c_terminal_overhang", "CGAG")),
            description="L0 acceptor C",
            **_DD,
        )
        self.cds1_c = widgets.Text(
            value=str(_junction("cds1_to_cds2", "upstream_c", "CTTC")),
            description="CDS1 C overhang",
            **_DD,
        )
        self.cds2_n = widgets.Text(
            value=str(_junction("cds1_to_cds2", "downstream_n", "GAAG")),
            description="CDS2 N overhang",
            **_DD,
        )
        self.cds1_c_offset = widgets.IntSlider(
            value=int(_junction("cds1_to_cds2", "arelf_offset_nt", 11)),
            min=0,
            max=11,
            step=1,
            description="CDS1/2 ARELF cut",
            continuous_update=False,
            **_DD_WIDE,
        )
        self.cds1_14_c = widgets.Text(
            value=str(_junction("cds1_to_cds14", "upstream_c", "GTGA")),
            description="CDS1→14 C",
            **_DD,
        )
        self.cds14_n = widgets.Text(
            value=str(_junction("cds1_to_cds14", "downstream_n", "TCAC")),
            description="CDS14 N",
            **_DD,
        )
        self.cds1_14_offset = widgets.IntSlider(
            value=int(_junction("cds1_to_cds14", "arelf_offset_nt", 4)),
            min=0,
            max=11,
            step=1,
            description="CDS1/14 ARELF cut",
            continuous_update=False,
            **_DD_WIDE,
        )
        self.cds14_19_c = widgets.Text(
            value=str(_junction("cds14_to_cds19", "upstream_c", "CACG")),
            description="CDS14→19 C",
            **_DD,
        )
        self.cds19_n = widgets.Text(
            value=str(_junction("cds14_to_cds19", "downstream_n", "CGTG")),
            description="CDS19 N",
            **_DD,
        )
        self.cds14_19_offset = widgets.IntSlider(
            value=int(_junction("cds14_to_cds19", "arelf_offset_nt", 1)),
            min=0,
            max=11,
            step=1,
            description="CDS14/19 ARELF cut",
            continuous_update=False,
            **_DD_WIDE,
        )
        self.level1_acceptor_name = widgets.Text(
            value=str(level1.get("acceptor_name", "Custom Level 1 acceptor")),
            description="Level 1 acceptor",
            **_DD_WIDE,
        )
        self.level1_n = widgets.Text(
            value=str(level1.get("n_terminal_overhang", "GCCC")),
            description="Cassette N overhang",
            **_DD,
        )
        self.level1_c = widgets.Text(
            value=str(level1.get("c_terminal_overhang", "GCGA")),
            description="Cassette C overhang",
            **_DD,
        )
        self.apply_btn = widgets.Button(
            description="Apply settings",
            button_style="success",
            icon="check",
            layout=widgets.Layout(width="180px", height="36px"),
        )
        self.reload_btn = widgets.Button(
            description="Reload GRASP GenBank",
            button_style="info",
            icon="refresh",
            layout=widgets.Layout(width="200px", height="36px"),
        )
        # Persistent HTML status (Output widgets often fail to refresh in Cursor/VS Code)
        self.status_html = widgets.HTML(
            value="<div style='font-family:-apple-system,sans-serif;font-size:13px;"
            "color:#3d5248;padding:8px;background:#eef2ef;border:1px solid #c9d0c8'>"
            "Select an organism — status updates when the dropdown changes.</div>"
        )
        self.out = widgets.Output()  # kept for reload logs only
        self._applying = False
        self._wire_events()

    def _wire_events(self) -> None:
        self.apply_btn.on_click(lambda _: self.apply())
        self.organism.observe(self._on_organism_change, names="value")
        self.kazusa_id.on_submit(lambda _: self.apply())
        self.upload.observe(self._on_upload_change, names="value")
        self.interface_preset.observe(
            self._on_interface_preset_change, names="value"
        )
        # Auto-apply on change — button clicks are unreliable in Cursor
        for w in (
            self.genetic_code,
            self.architecture,
            self.nterm,
            self.vendor,
            self.enzyme,
            self.ligation,
            self.redesign,
            self.selection,
            self.depth,
            self.cds1_c_offset,
            self.cds1_14_offset,
            self.cds14_19_offset,
        ):
            w.observe(self._on_setting_change, names="value")
        self.target_rna.on_submit(lambda _: self.apply())

    def _sync_codon_source_visibility(self) -> None:
        choice = self.organism.value
        self.kazusa_id.layout.display = (
            "" if choice == FETCH_FROM_KAZUSA else "none"
        )
        self.upload.layout.display = (
            "" if choice == UPLOAD_OWN_TABLE else "none"
        )

    def _on_upload_change(self, change=None) -> None:
        if not change or change.get("name") != "value":
            return
        if not change.get("new"):
            return
        if self.organism.value != UPLOAD_OWN_TABLE:
            self.organism.value = UPLOAD_OWN_TABLE
            return
        self.apply()

    def _on_interface_preset_change(self, change=None) -> None:
        if self._applying or not change or change.get("name") != "value":
            return
        if change.get("old") == change.get("new"):
            return
        from .assembly_interfaces import resolve_assembly_interfaces

        profile = resolve_assembly_interfaces(preset=str(change["new"]))
        entry = profile["level_minus1_entry"]
        level0 = profile["level0"]
        final = profile["final_cassette"]
        outer = level0.get("acceptor_outer") or {
            "n_overhang_5p": "CTCA",
            "c_overhang_5p": "CGAG",
        }
        self._applying = True
        try:
            self.entry_vector_name.value = entry["vector_id"]
            self.entry_n.value = entry["n_overhang_5p"]
            self.entry_c.value = entry["c_overhang_5p"]
            self.level0_acceptor_name.value = level0["acceptor_id"]
            self.level0_acceptor_n.value = outer["n_overhang_5p"]
            self.level0_acceptor_c.value = outer["c_overhang_5p"]
            self.cds1_c.value = profile["junctions"]["terminal_to_cds2"][
                "upstream_c_5p"
            ]
            self.cds2_n.value = profile["junctions"]["terminal_to_cds2"][
                "downstream_n_5p"
            ]
            self.cds1_14_c.value = profile["junctions"]["cds1_to_cds14"][
                "upstream_c_5p"
            ]
            self.cds14_n.value = profile["junctions"]["cds1_to_cds14"][
                "downstream_n_5p"
            ]
            self.cds14_19_c.value = profile["junctions"]["cds14_to_cds19"][
                "upstream_c_5p"
            ]
            self.cds19_n.value = profile["junctions"]["cds14_to_cds19"][
                "downstream_n_5p"
            ]
            self.level1_acceptor_name.value = final["vector_id"]
            self.level1_n.value = final["n_overhang_5p"]
            self.level1_c.value = final["c_overhang_5p"]
        finally:
            self._applying = False
        self.apply()

    def widget(self) -> widgets.Widget:
        """Return the full control panel as a VBox."""
        header = widgets.HTML(
            value=(
                "<div style='font-family:Georgia,Times New Roman,serif;font-size:20px;"
                "color:#1c2b24;font-weight:700;margin:0 0 6px 0'>Run setup</div>"
                "<div style='font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;"
                "font-size:13px;color:#3d5248;margin-bottom:8px'>"
                "Changing <b>Codon table</b> updates status immediately. "
                "Amino acids are re-derived from the genetic code. "
                "Re-run optimize cells after you change host."
                "</div>"
            )
        )
        box = widgets.VBox(
            [
                header,
                _label("Expression host"),
                self.kazusa_help,
                self.organism,
                self.kazusa_id,
                self.upload,
                self.genetic_code,
                self.architecture,
                self.nterm,
                self.target_rna,
                _label("Synthesis & assembly"),
                self.vendor,
                self.enzyme,
                self.ligation,
                _label("Oligo → Level -1 entry (directional 5′ overhangs)"),
                self.interface_preset,
                self.entry_vector_name,
                widgets.HBox([self.entry_n, self.entry_c]),
                _label("Level -1 → Level 0 (BpiI release / acceptor)"),
                self.level0_acceptor_name,
                widgets.HBox([self.level0_acceptor_n, self.level0_acceptor_c]),
                _label("Level 0 block junctions (C/N reverse-complement pairs)"),
                widgets.HBox([self.cds1_c, self.cds2_n]),
                self.cds1_c_offset,
                widgets.HBox([self.cds1_14_c, self.cds14_n]),
                self.cds1_14_offset,
                widgets.HBox([self.cds14_19_c, self.cds19_n]),
                self.cds14_19_offset,
                _label("Resulting Level 1 cassette"),
                self.level1_acceptor_name,
                widgets.HBox([self.level1_n, self.level1_c]),
                _label("Optimizer"),
                self.redesign,
                self.selection,
                self.depth,
                widgets.HBox(
                    [self.apply_btn, self.reload_btn],
                    layout=widgets.Layout(margin="12px 0 0 0", gap="8px"),
                ),
                self.status_html,
                self.out,
            ],
            layout=widgets.Layout(
                border="1px solid #c9d0c8",
                padding="14px 16px",
                margin="4px 0 12px 0",
            ),
        )
        return box

    def display(self) -> None:
        ui.section(
            "Control panel",
            "Status updates when you change the codon-table dropdown (no button needed).",
        )
        display(self.widget())
        self.apply()

    def _set_status(self, body_html: str, *, warn_html: str = "", err_html: str = "") -> None:
        warn_block = (
            f"<div style='margin-top:8px;padding:8px;background:#fff6e8;"
            f"border-left:4px solid #e0b872;font-size:12px'>{warn_html}</div>"
            if warn_html
            else ""
        )
        err_block = (
            f"<div style='margin-top:8px;padding:8px;background:#fdecea;"
            f"border-left:4px solid #c0392b;font-size:12px'>{err_html}</div>"
            if err_html
            else ""
        )
        self.status_html.value = (
            "<div style='font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;"
            "font-size:13px;color:#1c2b24;padding:10px 12px;background:#d8ebe3;"
            f"border:1px solid #c9d0c8;line-height:1.45'>{body_html}{warn_block}{err_block}</div>"
        )

    def _on_organism_change(self, change=None) -> None:
        if not change or change.get("name") != "value":
            return
        if change.get("old") == change.get("new"):
            return
        self._sync_codon_source_visibility()
        meta = SAMPLE_CODON_TABLES.get(self.organism.value, {})
        code = int(meta.get("genetic_code", 1))
        # Update genetic code without double-firing apply from that observe
        self._applying = True
        try:
            if (
                self.organism.value not in {FETCH_FROM_KAZUSA, UPLOAD_OWN_TABLE}
                and code in {c for c, _ in GENETIC_CODE_OPTIONS}
            ):
                self.genetic_code.value = code
        finally:
            self._applying = False
        # Wait for species id / file when those modes are selected
        if self.organism.value == FETCH_FROM_KAZUSA and not str(
            self.kazusa_id.value or ""
        ).strip():
            self._set_status(
                "Enter a Kazusa <b>species</b> accession (or paste the "
                "showcodon.cgi URL), then press Enter / Apply."
            )
            return
        if self.organism.value == UPLOAD_OWN_TABLE and not self.upload.value:
            self._set_status(
                "Drag a codon table onto <b>Upload table</b> "
                "(CSV with <code>codon,frequency</code> or Kazusa text/HTML)."
            )
            return
        self.apply()

    def _on_setting_change(self, change=None) -> None:
        if self._applying:
            return
        if not change or change.get("name") != "value":
            return
        if change.get("old") == change.get("new"):
            return
        self.apply()

    def apply(self) -> Dict[str, Any]:
        """Push widget values into CONFIG, load codon table, validate cut-site AAs."""
        if self._applying:
            return self.config
        self._applying = True
        try:
            return self._apply_impl()
        finally:
            self._applying = False

    def _apply_impl(self) -> Dict[str, Any]:
        cfg = deepcopy(self.config)

        def _overhang(widget: widgets.Text, label: str) -> str:
            value = str(widget.value).strip().upper().replace("U", "T")
            if len(value) != 4 or set(value) - set("ACGT"):
                raise ValueError(f"{label} must be exactly four DNA bases (ACGT)")
            return value

        def _paired(upstream: widgets.Text, downstream: widgets.Text, label: str):
            from .dna import reverse_complement

            up = _overhang(upstream, f"{label} upstream C overhang")
            down = _overhang(downstream, f"{label} downstream N overhang")
            if reverse_complement(up) != down:
                raise ValueError(
                    f"{label}: directional terminal overhangs must be reverse "
                    f"complements ({up} pairs with {reverse_complement(up)}, not {down})"
                )
            return up, down

        cfg = apply_vendor_to_config(cfg, self.vendor.value)
        cfg = apply_enzyme_to_config(cfg, self.enzyme.value)
        cfg = apply_ligation_table_to_config(cfg, self.ligation.value)

        cfg["genetic_code"] = int(self.genetic_code.value)
        cfg["architecture"] = self.architecture.value
        cfg["nterm_overhang"] = "AGGT" if self.nterm.value.startswith("AGGT") else "AATG"
        cfg["target_rna"] = self.target_rna.value.strip().upper().replace("T", "U")
        cds1_c, cds2_n = _paired(self.cds1_c, self.cds2_n, "CDS1→CDS2")
        cds1_14_c, cds14_n = _paired(
            self.cds1_14_c, self.cds14_n, "CDS1→CDS14"
        )
        cds14_19_c, cds19_n = _paired(
            self.cds14_19_c, self.cds19_n, "CDS14→CDS19"
        )
        selected_profile = str(self.interface_preset.value)
        if selected_profile == "deposited_grasp":
            deposited_values = (
                str(self.entry_vector_name.value).strip() == "pAGM1311"
                and _overhang(self.entry_n, "Entry N overhang") == "ACAT"
                and _overhang(self.entry_c, "Entry C overhang") == "ACAA"
                and str(self.level0_acceptor_name.value).strip() == "pAGM9121"
                and _overhang(self.level0_acceptor_n, "Level 0 N") == "CTCA"
                and _overhang(self.level0_acceptor_c, "Level 0 C") == "CGAG"
                and cds1_c == "CTTC"
                and cds2_n == "GAAG"
                and cds1_14_c == "GTGA"
                and cds14_n == "TCAC"
                and cds14_19_c == "CACG"
                and cds19_n == "CGTG"
                and str(self.level1_acceptor_name.value).strip()
                == "modified_1-1R_pICH47802_lc_p15A_ori_"
                and _overhang(self.level1_n, "Level 1 N") == "GGAG"
                and _overhang(self.level1_c, "Level 1 C") == "CGCT"
            )
            if not deposited_values:
                # Never retain deposited-vector completion contexts after the
                # user changes an interface. The design becomes a custom vector.
                selected_profile = "custom"
                self.interface_preset.value = "custom"
                print("▶ Edited deposited interfaces; switched vector preset to Custom")
        cfg["assembly_interfaces"] = {
            "overhang_notation": "directional_terminal_5p",
            "level_minus1_entry": {
                "profile": selected_profile,
                "vector_name": str(self.entry_vector_name.value).strip(),
                "enzyme": "BsaI",
                "n_terminal_overhang": _overhang(self.entry_n, "Entry N overhang"),
                "c_terminal_overhang": _overhang(self.entry_c, "Entry C overhang"),
            },
            "level0": {
                "acceptor_name": str(self.level0_acceptor_name.value).strip(),
                "release_enzyme": "BpiI / BbsI",
                "acceptor_n_terminal_overhang": _overhang(
                    self.level0_acceptor_n, "Level 0 acceptor N overhang"
                ),
                "acceptor_c_terminal_overhang": _overhang(
                    self.level0_acceptor_c, "Level 0 acceptor C overhang"
                ),
                "block_junctions": {
                    "cds1_to_cds14": {
                        "upstream_c": cds1_14_c,
                        "downstream_n": cds14_n,
                        "arelf_offset_nt": int(self.cds1_14_offset.value),
                    },
                    "cds14_to_cds19": {
                        "upstream_c": cds14_19_c,
                        "downstream_n": cds19_n,
                        "arelf_offset_nt": int(self.cds14_19_offset.value),
                    },
                    "cds1_to_cds2": {
                        "upstream_c": cds1_c,
                        "downstream_n": cds2_n,
                        "arelf_offset_nt": int(self.cds1_c_offset.value),
                    },
                },
            },
            "level1": {
                "acceptor_name": str(self.level1_acceptor_name.value).strip(),
                "n_terminal_overhang": _overhang(
                    self.level1_n, "Level 1 cassette N overhang"
                ),
                "c_terminal_overhang": _overhang(
                    self.level1_c, "Level 1 cassette C overhang"
                ),
            },
        }
        cfg["overhang_redesign"] = {
            "enabled": bool(self.redesign.value),
            "selection": self.selection.value,
            "cut_mode": "movable_arelf",
            "allowed_arelf_offsets_nt": list(range(12)),
        }
        cfg["optimizer"] = dict(cfg.get("optimizer", {}))
        cfg["optimizer"]["iterations_per_part"] = int(self.depth.value)
        cfg["codon_usage_file"] = self.input_dir / "codon_usage.csv"
        cfg["parts_file"] = self.input_dir / "parts.csv"
        cfg["target_map_file"] = self.input_dir / "target_map.csv"
        cfg["selected_organism"] = self.organism.value
        cfg["kazusa_species_id"] = str(self.kazusa_id.value or "").strip()

        self.config.clear()
        self.config.update(cfg)
        self.SELECTED_ORGANISM = self.organism.value

        upload_bytes = None
        upload_filename = None
        if self.organism.value == UPLOAD_OWN_TABLE and self.upload.value:
            entry = self._first_upload_entry(self.upload.value)
            if entry is not None:
                upload_bytes = entry["content"]
                upload_filename = entry.get("name")

        try:
            table, codon_data, meta, issues = apply_organism_codon_table(
                self.SELECTED_ORGANISM,
                Path(self.config["codon_usage_file"]),
                genetic_code=self.config["genetic_code"],
                kazusa_species_id=self.config.get("kazusa_species_id"),
                upload_bytes=upload_bytes,
                upload_filename=upload_filename,
            )
            self.CODON_TABLE = table
            self.CODON_DATA = codon_data
            if meta.get("species_id"):
                self.config["kazusa_species_id"] = meta["species_id"]
                if not str(self.kazusa_id.value or "").strip():
                    self.kazusa_id.value = str(meta["species_id"])
            if (
                meta.get("genetic_code")
                and self.organism.value == FETCH_FROM_KAZUSA
                and int(meta["genetic_code"]) in {c for c, _ in GENETIC_CODE_OPTIONS}
            ):
                # _applying is already True inside apply(); observers no-op
                self.genetic_code.value = int(meta["genetic_code"])
                self.config["genetic_code"] = int(meta["genetic_code"])

            parts_path = self.input_dir / "parts_full.csv"
            if not parts_path.exists():
                parts_path = Path(self.config["parts_file"])
            if parts_path.exists():
                import pandas as pd

                parts_df = pd.read_csv(parts_path)
                issues.extend(
                    validate_parts_for_organism(
                        parts_df,
                        codon_data,
                        genetic_code=self.config["genetic_code"],
                        keep_cut_sites=not self.config["overhang_redesign"]["enabled"],
                        minimum_relative_adaptiveness=self.config[
                            "codon_optimization"
                        ]["minimum_relative_adaptiveness"],
                    )
                )

            errors = [i for i in issues if i.level == "error"]
            warnings = [i for i in issues if i.level == "warning"]

            clade = meta.get("clade", "")
            link = (
                f' · <a href="{meta["url"]}" target="_blank">Kazusa</a>'
                if meta.get("url")
                else ""
            )
            synth = self.config["synthesis"]
            ligation = self.config.get("ligation", {})
            protocol = ligation.get("protocol_metadata", {})
            grasp_protocol = protocol.get("grasp_reference_protocol", {})
            if grasp_protocol:
                cycle_steps = " / ".join(
                    f"{step['temperature_c']} °C {step['minutes']} min"
                    for step in grasp_protocol.get("steps", [])
                )
                protocol_html = (
                    f"<br/>Ligation data: <b>{ligation.get('table_name')}</b><br/>"
                    f"GRASP reference cycling: <b>{grasp_protocol.get('cycles')}×</b> "
                    f"({cycle_steps}); fidelity matrix is a labelled Pryor proxy"
                )
            else:
                protocol_html = (
                    f"<br/>Ligation data: <b>{ligation.get('table_name')}</b>"
                )
            body = (
                f"<div style='font-size:10px;letter-spacing:0.12em;text-transform:uppercase;"
                f"color:#0f6b4c;margin-bottom:4px'>Active host</div>"
                f"<b>{meta.get('organism')}</b> ({clade}) · "
                f"genetic code <b>{self.config['genetic_code']}</b>{link}<br/>"
                f"Codons: <b>{len(table)}</b> sense · AA labels forced from genetic code<br/>"
                f"Architecture: <b>{self.config['architecture']}</b> · "
                f"N-term: <b>{self.config['nterm_overhang']}</b> · "
                f"Target: <code>{self.config['target_rna']}</code><br/>"
                f"Synthesis: <b>{self.config['synthesis_vendor']}</b> · "
                f"GC {synth['global_gc_min']*100:.0f}–{synth['global_gc_max']*100:.0f}% · "
                f"max HP {synth['max_homopolymer']}<br/>"
                f"Cut sites: <b>"
                f"{'KEEP native (AA risk check on)' if not self.config['overhang_redesign']['enabled'] else 'redesign ON (synonymous)'}"
                f"</b> · depth <b>{self.config['optimizer']['iterations_per_part']:,}</b>"
                f"{protocol_html}"
            )
            self._set_status(
                body,
                warn_html="<br/>".join(w.message for w in warnings),
                err_html="<br/>".join(e.message for e in errors),
            )
            print(
                f"▶ Host set to {meta.get('organism')} "
                f"(code {self.config['genetic_code']}, {len(table)} sense codons)"
            )
            if warnings or errors:
                print(format_issues(warnings + errors))

            with open(self.input_dir / "config.yaml", "w") as handle:
                yaml.safe_dump(
                    __import__("json").loads(
                        __import__("json").dumps(self.config, default=str)
                    ),
                    handle,
                    sort_keys=False,
                )
        except Exception as exc:
            self._set_status(
                f"<b>Could not apply settings</b><br/>{exc}",
                err_html=str(exc),
            )
            print(f"✗ Apply failed: {exc}")
            raise

        if self.on_apply:
            self.on_apply(self.config)
        return self.config


    @staticmethod
    def _first_upload_entry(value) -> Optional[dict]:
        """Normalize ipywidgets FileUpload payload across v7/v8 shapes."""
        if not value:
            return None
        if isinstance(value, dict):
            # v7: {filename: {"metadata":…, "content": b"…"}}
            name, payload = next(iter(value.items()))
            if isinstance(payload, dict) and "content" in payload:
                return {
                    "name": payload.get("metadata", {}).get("name", name),
                    "content": payload["content"],
                }
            return None
        # v8: tuple/list of {name, type, size, content}
        first = value[0]
        if isinstance(first, dict) and "content" in first:
            return {"name": first.get("name"), "content": first["content"]}
        return None


def wire_reload_button(panel: GraspControlPanel, import_fn: Callable) -> None:
    """Attach GenBank reload to the panel reload button."""

    def _reload(_):
        try:
            result = import_fn()
            n = len(result["parts"]) if isinstance(result, dict) else "?"
            panel._set_status(f"Reloaded GRASP GenBank → <b>{n}</b> modules in parts.csv")
            print(f"▶ Reloaded GRASP GenBank ({n} modules)")
        except Exception as exc:
            panel._set_status(f"Reload failed: {exc}", err_html=str(exc))
            print(f"✗ Reload failed: {exc}")

    panel.reload_btn.on_click(_reload)
