"""Biologist-facing ipywidgets control panel for the GRASP notebook."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import ipywidgets as widgets
import yaml
from IPython.display import clear_output, display

from . import notebook_ui as ui
from .assembly_interfaces import (
    CANONICAL_NOTATION,
    FIVE_PRIME_CODING_SITE,
    FIVE_PRIME_END,
    THREE_PRIME_CODING_SITE,
    THREE_PRIME_END,
    deposited_grasp_interface_preset,
    reverse_complement as interface_reverse_complement,
)
from .codon_tables import apply_organism_codon_table, validate_parts_for_organism
from .codon_validation import format_issues
from .sample_codon_tables import (
    FETCH_FROM_KAZUSA,
    KAZUSA_REMINDER_HTML,
    SAMPLE_CODON_TABLES,
    UPLOAD_OWN_TABLE,
    sample_names,
)
from .restriction_sites import (
    DEFAULT_SITE_BLACKLIST,
    apply_site_blacklist_to_config,
    restriction_site_options,
)
from .synthesis_vendors import (
    GRASP_STAGE_MATCHED_LIGATION,
    apply_enzyme_to_config,
    apply_ligation_table_to_config,
    apply_vendor_to_config,
    enzyme_names,
    redesign_ligation_table_names,
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
    toolbox = deposited_grasp_interface_preset()
    entry_defaults = toolbox["level_minus1_entry"]
    level0_defaults = toolbox["level0"]["acceptor_outer"]
    level1_defaults = toolbox["final_cassette"]
    junction_offsets = {
        "terminal_to_cds2": 11,
        "cds1_to_cds14": 4,
        "cds14_to_cds19": 1,
    }
    cfg: Dict[str, Any] = {
        "project_name": "GRASP_9S_organism_optimized",
        "genetic_code": 1,
        "architecture": "9S",
        "ppr_5prime_fusion_site": toolbox["level0"]["ppr_outer"][
            FIVE_PRIME_CODING_SITE
        ],
        "assembly_interfaces": {
            # Every overhang is named only by cloning level and sequence end.
            # All sequence values are written 5'->3'.
            "terminal_side_convention": "construct_ends_5prime_and_3prime",
            "overhang_sequence_notation": "5prime_to_3prime",
            "notation": CANONICAL_NOTATION,
            "coding_strand_direction": "5prime_to_3prime",
            "level_minus1_entry": {
                "profile": "deposited_grasp",
                "vector_name": entry_defaults["vector_id"],
                "enzyme": "BsaI",
                FIVE_PRIME_END: entry_defaults[FIVE_PRIME_END],
                THREE_PRIME_END: entry_defaults[THREE_PRIME_END],
                FIVE_PRIME_CODING_SITE: entry_defaults[FIVE_PRIME_CODING_SITE],
                THREE_PRIME_CODING_SITE: entry_defaults[THREE_PRIME_CODING_SITE],
            },
            "level0": {
                "acceptor_name": toolbox["level0"]["acceptor_id"],
                "release_enzyme": "BpiI / BbsI",
                "acceptor_outer": {
                    FIVE_PRIME_END: level0_defaults[FIVE_PRIME_END],
                    THREE_PRIME_END: level0_defaults[THREE_PRIME_END],
                    FIVE_PRIME_CODING_SITE: level0_defaults[FIVE_PRIME_CODING_SITE],
                    THREE_PRIME_CODING_SITE: level0_defaults[THREE_PRIME_CODING_SITE],
                },
            },
            "junctions": {
                name: {**values, "arelf_offset_nt": junction_offsets[name]}
                for name, values in toolbox["junctions"].items()
            },
            "level1": {
                "acceptor_name": level1_defaults["vector_id"],
                FIVE_PRIME_END: level1_defaults[FIVE_PRIME_END],
                THREE_PRIME_END: level1_defaults[THREE_PRIME_END],
                FIVE_PRIME_CODING_SITE: level1_defaults[FIVE_PRIME_CODING_SITE],
                THREE_PRIME_CODING_SITE: level1_defaults[THREE_PRIME_CODING_SITE],
            },
        },
        "assembly_enzyme": "GRASP default · BsaI + BpiI + BsmBI",
        "synthesis_vendor": "Twist · Standard gene guidelines",
        "codon_usage_file": input_dir / "codon_usage.csv",
        "parts_file": input_dir / "parts.csv",
        "target_map_file": input_dir / "target_map.csv",
        "site_blacklist": list(DEFAULT_SITE_BLACKLIST),
        "forbidden_sites": {
            "BsaI": "GGTCTC",
            "BpiI": "GAAGAC",
            "BsmBI": "CGTCTC",
            "SapI": "GCTCTTC",
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
        "overhang_redesign": {
            "plasmid_overhangs": False,
            "level0_junctions": False,
            # Backward-compatible alias for Level 0 junction redesign.
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
    cfg = apply_ligation_table_to_config(cfg, GRASP_STAGE_MATCHED_LIGATION)
    cfg["ligation"]["min_efficiency"] = 0.25
    cfg["ligation"]["min_fidelity"] = 0.9
    return cfg


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
        toolbox = deposited_grasp_interface_preset()
        entry_defaults = toolbox["level_minus1_entry"]
        level0_defaults = toolbox["level0"]["acceptor_outer"]
        level1_defaults = toolbox["final_cassette"]
        entry = assembly.get("level_minus1_entry", {})
        level0 = assembly.get("level0", {})
        level1 = assembly.get("level1", {})
        acceptor_outer = level0.get("acceptor_outer") or {}

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
        known_sites = {name for _, name in restriction_site_options()}
        selected_blacklist = [
            name
            for name in self.config.get("site_blacklist", DEFAULT_SITE_BLACKLIST)
            if name in known_sites
        ]
        self.site_blacklist = widgets.SelectMultiple(
            options=restriction_site_options(),
            value=selected_blacklist or list(DEFAULT_SITE_BLACKLIST),
            description="Cut-site blacklist",
            rows=10,
            style={"description_width": "140px"},
            layout=widgets.Layout(width="520px"),
        )
        self.ligation = widgets.Dropdown(
            options=redesign_ligation_table_names(),
            value=self.config.get("ligation", {}).get(
                "table_name",
                GRASP_STAGE_MATCHED_LIGATION,
            ),
            description="Ligation scoring",
            **_DD_WIDE,
        )
        self.architecture = widgets.Dropdown(
            options=["9S", "14S", "19S"],
            value=self.config.get("architecture", "9S"),
            description="Architecture",
            **_DD,
        )
        redesign_cfg = self.config.get("overhang_redesign", {})
        self.redesign_plasmid = widgets.Checkbox(
            value=bool(redesign_cfg.get("plasmid_overhangs", False)),
            description="Redesign plasmid overhangs",
            indent=False,
            layout=widgets.Layout(width="520px"),
        )
        self.redesign_level0 = widgets.Checkbox(
            value=bool(
                redesign_cfg.get(
                    "level0_junctions",
                    redesign_cfg.get("enabled", False),
                )
            ),
            description="Redesign Level 0 overhangs (assembled parts)",
            indent=False,
            layout=widgets.Layout(width="520px"),
        )
        self.selection = widgets.Dropdown(
            options=["knee", "max_fidelity"],
            value=redesign_cfg.get("selection", "knee"),
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
        self.level_minus1_5prime = widgets.Text(
            value=str(entry.get(FIVE_PRIME_END, entry_defaults[FIVE_PRIME_END])),
            description="5′ overhang",
            **_DD,
        )
        self.level_minus1_3prime = widgets.Text(
            value=str(entry.get(THREE_PRIME_END, entry_defaults[THREE_PRIME_END])),
            description="3′ overhang",
            **_DD,
        )
        self.level0_5prime = widgets.Text(
            value=str(
                acceptor_outer.get(
                    FIVE_PRIME_END,
                    level0_defaults[FIVE_PRIME_END],
                )
            ),
            description="5′ overhang",
            **_DD,
        )
        self.level0_3prime = widgets.Text(
            value=str(
                acceptor_outer.get(
                    THREE_PRIME_END,
                    level0_defaults[THREE_PRIME_END],
                )
            ),
            description="3′ overhang",
            **_DD,
        )
        self.level1_5prime = widgets.Text(
            value=str(level1.get(FIVE_PRIME_END, level1_defaults[FIVE_PRIME_END])),
            description="5′ overhang",
            **_DD,
        )
        self.level1_3prime = widgets.Text(
            value=str(level1.get(THREE_PRIME_END, level1_defaults[THREE_PRIME_END])),
            description="3′ overhang",
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
        # Auto-apply on change — button clicks are unreliable in Cursor
        for w in (
            self.genetic_code,
            self.architecture,
            self.vendor,
            self.enzyme,
            self.site_blacklist,
            self.ligation,
            self.redesign_plasmid,
            self.redesign_level0,
            self.selection,
            self.depth,
        ):
            w.observe(self._on_setting_change, names="value")
        self.target_rna.on_submit(lambda _: self.apply())
        self._sync_redesign_visibility()
        self.redesign_plasmid.observe(self._on_redesign_toggle, names="value")
        self.redesign_level0.observe(self._on_redesign_toggle, names="value")

    def _on_redesign_toggle(self, change=None) -> None:
        if not change or change.get("name") != "value":
            return
        self._sync_redesign_visibility()
        self._on_setting_change(change)

    def _sync_redesign_visibility(self) -> None:
        """Show Pareto picker only when Level 0 junction redesign is on."""
        self.selection.layout.display = "" if self.redesign_level0.value else "none"

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
                self.target_rna,
                _label("Synthesis & assembly"),
                self.vendor,
                self.enzyme,
                self.site_blacklist,
                widgets.HTML(
                    value=(
                        "<div style='font-family:-apple-system,sans-serif;font-size:12px;"
                        "color:#3d5248;margin:2px 0 8px 0'>"
                        "⌘/Ctrl-click to select several. Default: <b>SapI, BsaI, BpiI</b>. "
                        "These sites are removed from designed CDS, in addition to the "
                        "assembly-enzyme filter above."
                        "</div>"
                    )
                ),
                self.ligation,
                widgets.HTML(
                    value=(
                        "<div style='font-family:-apple-system,sans-serif;font-size:12px;"
                        "color:#3d5248;margin:2px 0 8px 0;padding:8px 10px;"
                        "background:#eef2ef;border-left:4px solid #0f6b4c'>"
                        "<b>One ligation menu, two enzyme scores.</b><br/>"
                        "• Levels −1 &amp; 1 are scored with <b>BsaI-HFv2</b> "
                        "37↔16 °C (Pryor) and reported on the Pareto table.<br/>"
                        "• Level 0 five-part redesign is scored with <b>BbsI-HF</b> "
                        "(BpiI isoschizomer) — that Level 0 score is the ligation "
                        "objective used for search, unless you pick a Potapov "
                        "Level 0 override."
                        "</div>"
                    )
                ),
                _label("Plasmid / acceptor overhangs (typed · not A–E junctions)"),
                widgets.HTML(
                    value=(
                        "<div style='font-family:-apple-system,sans-serif;font-size:12px;"
                        "color:#3d5248;margin:2px 0 8px 0'>"
                        "These are the <b>vector-facing</b> sticky ends for each cloning "
                        "level (defaults = deposited GRASP toolbox). They are "
                        "<b>not</b> the internal A–E Level 0 junctions.</div>"
                    )
                ),
                _label("Level −1 entry vector"),
                widgets.HBox(
                    [self.level_minus1_5prime, self.level_minus1_3prime]
                ),
                _label("Level 0 acceptor outer"),
                widgets.HBox([self.level0_5prime, self.level0_3prime]),
                _label("Level 1 acceptor outer"),
                widgets.HBox([self.level1_5prime, self.level1_3prime]),
                _label("Overhang redesign"),
                widgets.HTML(
                    value=(
                        "<div style='font-family:-apple-system,sans-serif;font-size:12px;"
                        "color:#3d5248;margin:2px 0 8px 0;padding:8px 10px;"
                        "background:#fff6e8;border-left:4px solid #e0b872'>"
                        "<b>Redesign plasmid overhangs</b> (default off): keep the "
                        "typed Level −1 / 0 / 1 fields above as fixed backbone ends. "
                        "Turn on only when you intentionally change those fields from "
                        "the deposited toolbox.<br/><br/>"
                        "<b>Redesign Level 0 overhangs</b> (default off): Pareto search "
                        "over the <b>assembled five-part junctions</b> inside ARELF:<br/>"
                        "A 3′ (ACTC) · B both · C both · D both · E 5′ (TGAA).<br/>"
                        "Does <b>not</b> move A 5′ / E 3′ acceptor outers or MoClo "
                        "sites (AGGT / CTTC / TTCG).</div>"
                    )
                ),
                self.redesign_plasmid,
                self.redesign_level0,
                self.selection,
                _label("Optimizer"),
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

        cfg = apply_vendor_to_config(cfg, self.vendor.value)
        cfg = apply_enzyme_to_config(cfg, self.enzyme.value)
        cfg = apply_site_blacklist_to_config(cfg, self.site_blacklist.value)
        cfg = apply_ligation_table_to_config(cfg, self.ligation.value)

        cfg["genetic_code"] = int(self.genetic_code.value)
        cfg["architecture"] = self.architecture.value
        toolbox_defaults = deposited_grasp_interface_preset()
        ppr_5prime = toolbox_defaults["level0"]["ppr_outer"][
            FIVE_PRIME_CODING_SITE
        ]
        cfg["ppr_5prime_fusion_site"] = ppr_5prime
        cfg["target_rna"] = self.target_rna.value.strip().upper().replace("T", "U")
        level_minus1_5prime = _overhang(
            self.level_minus1_5prime, "Level −1 5′ overhang"
        )
        level_minus1_3prime = _overhang(
            self.level_minus1_3prime, "Level −1 3′ overhang"
        )
        level0_5prime = _overhang(self.level0_5prime, "Level 0 5′ overhang")
        level0_3prime = _overhang(self.level0_3prime, "Level 0 3′ overhang")
        level1_5prime = _overhang(self.level1_5prime, "Level 1 5′ overhang")
        level1_3prime = _overhang(self.level1_3prime, "Level 1 3′ overhang")

        entry_defaults = toolbox_defaults["level_minus1_entry"]
        level0_defaults = toolbox_defaults["level0"]["acceptor_outer"]
        level1_defaults = toolbox_defaults["final_cassette"]
        standard_values = (
            level_minus1_5prime == entry_defaults[FIVE_PRIME_END]
            and level_minus1_3prime == entry_defaults[THREE_PRIME_END]
            and level0_5prime == level0_defaults[FIVE_PRIME_END]
            and level0_3prime == level0_defaults[THREE_PRIME_END]
            and level1_5prime == level1_defaults[FIVE_PRIME_END]
            and level1_3prime == level1_defaults[THREE_PRIME_END]
        )
        selected_profile = "deposited_grasp" if standard_values else "custom"
        entry_vector = entry_defaults["vector_id"] if (
            level_minus1_5prime,
            level_minus1_3prime,
        ) == (
            entry_defaults[FIVE_PRIME_END],
            entry_defaults[THREE_PRIME_END],
        ) else "custom_level_minus1_entry"
        level0_vector = toolbox_defaults["level0"]["acceptor_id"] if (
            level0_5prime,
            level0_3prime,
        ) == (
            level0_defaults[FIVE_PRIME_END],
            level0_defaults[THREE_PRIME_END],
        ) else "custom_level0_acceptor"
        level1_vector = (
            level1_defaults["vector_id"]
            if (level1_5prime, level1_3prime)
            == (
                level1_defaults[FIVE_PRIME_END],
                level1_defaults[THREE_PRIME_END],
            )
            else "custom_level1_acceptor"
        )
        junctions = {
            name: {
                "upstream_three_prime_end_overhang": values[
                    "upstream_three_prime_end_overhang"
                ],
                "downstream_five_prime_end_overhang": values[
                    "downstream_five_prime_end_overhang"
                ],
                "assembled_coding_site": values["assembled_coding_site"],
                "arelf_offset_nt": values.get("arelf_offset_nt", offset),
            }
            for (name, values), offset in zip(
                toolbox_defaults["junctions"].items(), (11, 4, 1)
            )
        }
        cfg["assembly_interfaces"] = {
            "preset": selected_profile,
            "terminal_side_convention": "construct_ends_5prime_and_3prime",
            "overhang_sequence_notation": "5prime_to_3prime",
            "notation": CANONICAL_NOTATION,
            "coding_strand_direction": "5prime_to_3prime",
            "level_minus1_entry": {
                "profile": selected_profile,
                "vector_name": entry_vector,
                "enzyme": "BsaI",
                FIVE_PRIME_END: level_minus1_5prime,
                THREE_PRIME_END: level_minus1_3prime,
                FIVE_PRIME_CODING_SITE: level_minus1_5prime,
                THREE_PRIME_CODING_SITE: interface_reverse_complement(
                    level_minus1_3prime
                ),
            },
            "level0": {
                "acceptor_name": level0_vector,
                "release_enzyme": "BpiI / BbsI",
                "acceptor_outer": {
                    FIVE_PRIME_END: level0_5prime,
                    THREE_PRIME_END: level0_3prime,
                    FIVE_PRIME_CODING_SITE: level0_5prime,
                    THREE_PRIME_CODING_SITE: interface_reverse_complement(
                        level0_3prime
                    ),
                },
            },
            "junctions": junctions,
            "level1": {
                "acceptor_name": level1_vector,
                FIVE_PRIME_END: level1_5prime,
                THREE_PRIME_END: level1_3prime,
                FIVE_PRIME_CODING_SITE: level1_5prime,
                THREE_PRIME_CODING_SITE: interface_reverse_complement(
                    level1_3prime
                ),
            },
        }
        cfg["overhang_redesign"] = {
            "plasmid_overhangs": bool(self.redesign_plasmid.value),
            "level0_junctions": bool(self.redesign_level0.value),
            "enabled": bool(self.redesign_level0.value),
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
                        keep_cut_sites=not bool(
                            self.config.get("overhang_redesign", {}).get(
                                "level0_junctions",
                                self.config.get("overhang_redesign", {}).get(
                                    "enabled", False
                                ),
                            )
                        ),
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
            by_level = ligation.get("by_level", {})
            protocol = ligation.get("protocol_metadata", {})
            grasp_protocol = protocol.get("grasp_reference_protocol", {})
            level_lines = []
            for key, label in (
                ("level_minus1", "Level −1 · BsaI-HFv2"),
                ("level0", "Level 0 · BbsI-HF / BpiI"),
                ("level1", "Level 1 · BsaI-HFv2"),
            ):
                block = by_level.get(key, {})
                enzyme = (block.get("protocol_metadata") or {}).get(
                    "restriction_enzyme", "?"
                )
                table = block.get("ligation_table", "?")
                level_lines.append(f"{label}: <code>{table}</code> ({enzyme})")
            levels_html = "<br/>".join(level_lines)
            if grasp_protocol:
                cycle_steps = " / ".join(
                    f"{step['temperature_c']} °C {step['minutes']} min"
                    for step in grasp_protocol.get("steps", [])
                )
                protocol_html = (
                    f"<br/>Level 0 redesign matrix: <b>{ligation.get('table_name')}</b>"
                    f"<br/>{levels_html}"
                    f"<br/>GRASP reference cycling: <b>{grasp_protocol.get('cycles')}×</b> "
                    f"({cycle_steps}); Pryor matrices are labelled proxies"
                )
            else:
                protocol_html = (
                    f"<br/>Level 0 redesign matrix: <b>{ligation.get('table_name')}</b>"
                    f"<br/>{levels_html}"
                )
            body = (
                f"<div style='font-size:10px;letter-spacing:0.12em;text-transform:uppercase;"
                f"color:#0f6b4c;margin-bottom:4px'>Active host</div>"
                f"<b>{meta.get('organism')}</b> ({clade}) · "
                f"genetic code <b>{self.config['genetic_code']}</b>{link}<br/>"
                f"Codons: <b>{len(table)}</b> sense · AA labels forced from genetic code<br/>"
                f"Architecture: <b>{self.config['architecture']}</b> · "
                f"PPR 5′ fusion site: <b>{self.config['ppr_5prime_fusion_site']}</b> · "
                f"Target: <code>{self.config['target_rna']}</code><br/>"
                f"Synthesis: <b>{self.config['synthesis_vendor']}</b> · "
                f"GC {synth['global_gc_min']*100:.0f}–{synth['global_gc_max']*100:.0f}% · "
                f"max HP {synth['max_homopolymer']}<br/>"
                f"Cut-site blacklist: <b>"
                f"{', '.join(self.config.get('site_blacklist') or []) or 'none'}</b> · "
                f"forbidden: <code>"
                f"{', '.join(self.config.get('forbidden_sites') or []) or 'none'}</code><br/>"
                f"Plasmid overhangs: <b>"
                f"{'CUSTOM fields' if self.config['overhang_redesign'].get('plasmid_overhangs') else 'fields as fixed'}"
                f"</b> · Level 0 junctions: <b>"
                f"{'REDESIGN A3′/B/C/D/E5′' if self.config['overhang_redesign'].get('level0_junctions', self.config['overhang_redesign'].get('enabled')) else 'KEEP native ACTC/AAGA/GCAC/TGAA'}"
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
