"""Notebook presentation that actually renders in VS Code / Cursor / Jupyter.

VS Code notebook webviews strip bare <style> tags. Inline styles on elements
do render — so every visual here uses inline CSS only.
"""

from __future__ import annotations

from IPython.display import HTML, Markdown, display

# Lab-notebook palette (inline only — no purple AI defaults)
INK = "#1c2b24"
INK_SOFT = "#3d5248"
PAPER = "#f4f1ea"
PANEL = "#ffffff"
RULE = "#c9d0c8"
ACCENT = "#0f6b4c"
ACCENT_SOFT = "#d8ebe3"
WARN_BG = "#fff6e8"
WARN_BORDER = "#e0b872"
CODE_BG = "#eef2ef"


def _html(fragment: str) -> None:
    display(HTML(fragment))


def inject_styles() -> None:
    """Kept for API compatibility; styling is inline."""
    return None


def hero() -> None:
    _html(
        f"""
<div style="
  font-family: Georgia, 'Times New Roman', serif;
  background: {PAPER};
  border: 1px solid {RULE};
  border-left: 6px solid {ACCENT};
  padding: 22px 26px 18px 26px;
  margin: 0 0 14px 0;
">
  <div style="
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    font-size: 11px;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    color: {ACCENT};
    margin-bottom: 8px;
  ">Golden Gate · codon · synthesis</div>
  <div style="
    font-size: 34px;
    line-height: 1.15;
    color: {INK};
    margin: 0 0 10px 0;
    font-weight: 700;
  ">GRASP Library Designer</div>
  <div style="
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    font-size: 14px;
    line-height: 1.5;
    color: {INK_SOFT};
    max-width: 52rem;
  ">
    Redesign GRASP module DNA for your organism without moving cut sites or
    changing proteins. Pick options below — then run the cells beneath.
  </div>
  <div style="
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    margin-top: 14px;
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
  ">
    {_chip("fixed cut indices")}
    {_chip("overhang sequences")}
    {_chip("Pareto front")}
    {_chip("organism codon tables")}
  </div>
</div>
"""
    )


def _chip(text: str) -> str:
    return f"""
<span style="
  font-size: 12px;
  color: {INK};
  background: {ACCENT_SOFT};
  border: 1px solid {RULE};
  padding: 4px 10px;
">{text}</span>
"""


def section(title: str, subtitle: str = "") -> None:
    sub = (
        f"""<div style="
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          font-size: 13px;
          color: {INK_SOFT};
          margin-top: 4px;
        ">{subtitle}</div>"""
        if subtitle
        else ""
    )
    _html(
        f"""
<div style="
  margin: 22px 0 10px 0;
  padding: 0 0 8px 0;
  border-bottom: 2px solid {ACCENT};
">
  <div style="
    font-family: Georgia, 'Times New Roman', serif;
    font-size: 22px;
    color: {INK};
    font-weight: 700;
  ">{title}</div>
  {sub}
</div>
"""
    )


def note(text: str) -> None:
    _html(
        f"""
<div style="
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  background: {WARN_BG};
  border: 1px solid {WARN_BORDER};
  border-left: 4px solid {WARN_BORDER};
  padding: 10px 14px;
  margin: 8px 0 12px 0;
  color: {INK};
  font-size: 13px;
  line-height: 1.45;
"><b style="color:{ACCENT}">Note.</b> {text}</div>
"""
    )


def status(text: str) -> None:
    cleaned = (
        text.replace("<br/>", "<br>")
        .replace("<b>", f"<b style='color:{ACCENT}'>")
        .replace("<code>", f"<code style='background:{CODE_BG};padding:1px 5px;border-radius:3px;font-size:12px'>")
    )
    _html(
        f"""
<div style="
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  background: {ACCENT_SOFT};
  border: 1px solid {RULE};
  padding: 10px 14px;
  margin: 8px 0 12px 0;
  color: {INK};
  font-size: 13px;
  line-height: 1.5;
"><span style="
  display:inline-block;
  font-size:10px;
  letter-spacing:0.12em;
  text-transform:uppercase;
  color:{ACCENT};
  margin-right:8px;
">Status</span>{cleaned}</div>
"""
    )


def panel(title: str, body_html: str) -> None:
    _html(
        f"""
<div style="
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  background: {PANEL};
  border: 1px solid {RULE};
  padding: 16px 18px;
  margin: 8px 0 14px 0;
">
  <div style="
    font-family: Georgia, 'Times New Roman', serif;
    font-size: 16px;
    color: {INK};
    margin-bottom: 10px;
    font-weight: 700;
  ">{title}</div>
  {body_html}
</div>
"""
    )


def markdown(text: str) -> None:
    display(Markdown(text))
