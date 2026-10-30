#!/usr/bin/env python
"""Render the three reports to A4 PDFs (pandoc + tectonic, Libertinus fonts from the TeX bundle).

Usage:  PATH=/tmp/render/bin:$PATH python reports/render/render_pdfs.py
(scratch env: mamba create -p /tmp/render -c conda-forge pandoc tectonic poppler)
The manuscript keeps its figure legends at the end; the PDF places each figure image above its legend.
"""
import pathlib
import re
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = pathlib.Path(__file__).resolve().parent
FIGS = {str(i): n for i, n in enumerate(["", "fig1_system_overview", "fig2_speed", "fig3_sensitivity",
                                         "fig4_error_control", "fig5_ptm_quant", "fig6_gate_history"]) if n}
DOCS = [("briefing/BRIEFING.md", "11pt"), ("technical/TECHNICAL_REPORT.md", "10pt"), ("manuscript/MANUSCRIPT.md", "10pt")]


def embed_manuscript_figures(md: str) -> str:
    def rep(m):
        k = m.group(1)
        return f"![](../figures/{FIGS[k]}.png){{width=100%}}\n\n{m.group(0)}"
    return re.sub(r"^\*\*Figure (\d) \|", rep, md, flags=re.M)


def render(rel: str, size: str):
    src = ROOT / "reports" / rel
    md = src.read_text()
    if "MANUSCRIPT" in rel:
        md = embed_manuscript_figures(md)
    out = src.with_suffix(".pdf")
    with tempfile.NamedTemporaryFile("w", suffix=".md", dir=src.parent, delete=False) as tmp:
        tmp.write(md)
    try:
        cmd = ["pandoc", tmp.name, "-f", "markdown-implicit_figures+pipe_tables+lists_without_preceding_blankline", "-o", str(out),
               "--pdf-engine=tectonic", f"--resource-path={src.parent}",
               "--lua-filter", str(HERE / "colwidths.lua"), "-H", str(HERE / "header.tex"),
               "-V", "papersize=a4", "-V", "geometry:margin=2cm", "-V", f"fontsize={size}",
               "-V", "mainfont=LibertinusSerif-Regular.otf",
               "-V", "mainfontoptions=BoldFont=LibertinusSerif-Bold.otf,ItalicFont=LibertinusSerif-Italic.otf,BoldItalicFont=LibertinusSerif-BoldItalic.otf",
               "-V", "sansfont=LibertinusSans-Regular.otf",
               "-V", "sansfontoptions=BoldFont=LibertinusSans-Bold.otf,ItalicFont=LibertinusSans-Italic.otf",
               "-V", "monofont=LibertinusMono-Regular.otf", "-V", "monofontoptions=Scale=0.9",
               "-V", "colorlinks=true", "-V", "linkcolor=blue!50!black", "-V", "urlcolor=blue!50!black",
               "-V", "linestretch=1.08"]
        subprocess.run(cmd, check=True, cwd=src.parent)
    finally:
        pathlib.Path(tmp.name).unlink()
    print("wrote", out.relative_to(ROOT))


if __name__ == "__main__":
    for rel, size in DOCS:
        render(rel, size)
