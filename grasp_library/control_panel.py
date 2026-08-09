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
from .sample_codon_tables import SAMPLE_CODON_TABLES, sample_names
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
            "max_homopolymer": 13,
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
            "temperature": 25,
            "hours": 18,
            "min_efficiency": 0.25,
            "min_fidelity": 0.9,
            "table_name": "T4 · 18 h · 25 °C (Potapov)",
            "ligation_table": None,
        },
        "overhang_redesign": {"enabled": True, "selection": "knee"},
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
        self.organism = widgets.Dropdown(
            options=sample_names(),
            value="Escherichia coli (Kazusa)",
            description="Codon table",
            **_DD_WIDE,
        )
        self.genetic_code = widgets.Dropdown(
            options=[(label, code) for code, label in GENETIC_CODE_OPTIONS],
            value=int(self.config.get("genetic_code", 1)),
            description="Genetic code",
            **_DD_WIDE,
        )
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
            description="Assembly enzyme",
            **_DD_WIDE,
        )
        self.ligation = widgets.Dropdown(
            options=ligation_table_names(),
            value=self.config.get("ligation", {}).get(
                "table_name", "T4 · 18 h · 25 °C (Potapov)"
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
                ("On — Pareto redesign at fixed cuts", True),
                ("Off — keep native GRASP overhangs", False),
            ],
            value=bool(self.config.get("overhang_redesign", {}).get("enabled", True)),
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
        ):
            w.observe(self._on_setting_change, names="value")
        self.target_rna.on_submit(lambda _: self.apply())

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
                self.organism,
                self.genetic_code,
                self.architecture,
                self.nterm,
                self.target_rna,
                _label("Synthesis & assembly"),
                self.vendor,
                self.enzyme,
                self.ligation,
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
        meta = SAMPLE_CODON_TABLES.get(self.organism.value, {})
        code = int(meta.get("genetic_code", 1))
        # Update genetic code without double-firing apply from that observe
        self._applying = True
        try:
            if code in {c for c, _ in GENETIC_CODE_OPTIONS}:
                self.genetic_code.value = code
        finally:
            self._applying = False
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

        cfg = apply_vendor_to_config(cfg, self.vendor.value)
        cfg = apply_enzyme_to_config(cfg, self.enzyme.value)
        cfg = apply_ligation_table_to_config(cfg, self.ligation.value)

        cfg["genetic_code"] = int(self.genetic_code.value)
        cfg["architecture"] = self.architecture.value
        cfg["nterm_overhang"] = "AGGT" if self.nterm.value.startswith("AGGT") else "AATG"
        cfg["target_rna"] = self.target_rna.value.strip().upper().replace("T", "U")
        cfg["overhang_redesign"] = {
            "enabled": bool(self.redesign.value),
            "selection": self.selection.value,
        }
        cfg["optimizer"] = dict(cfg.get("optimizer", {}))
        cfg["optimizer"]["iterations_per_part"] = int(self.depth.value)
        cfg["codon_usage_file"] = self.input_dir / "codon_usage.csv"
        cfg["parts_file"] = self.input_dir / "parts.csv"
        cfg["target_map_file"] = self.input_dir / "target_map.csv"
        cfg["selected_organism"] = self.organism.value

        self.config.clear()
        self.config.update(cfg)
        self.SELECTED_ORGANISM = self.organism.value

        try:
            table, codon_data, meta, issues = apply_organism_codon_table(
                self.SELECTED_ORGANISM,
                Path(self.config["codon_usage_file"]),
                genetic_code=self.config["genetic_code"],
            )
            self.CODON_TABLE = table
            self.CODON_DATA = codon_data

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
