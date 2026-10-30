#!/usr/bin/env python
"""Publication figures for the Vulcan reports (briefing, technical report, manuscript).

Reproduce:  source env/activate.sh && python reports/figures/make_figures.py
Writes fig1..fig6 as PDF + PNG (300 dpi) into reports/figures/ and regenerates
reports/figures/SOURCES.md (every plotted number -> file + JSON key + value).

Data policy: only committed dev results are read (results/2026-10-2*/ on master,
E12 from commit 7f64a72 on branch worktree-agent-a0df6277cef40ec85, and REPORT.md
for the 1 TB estimate). Nothing under data/sealed is touched. There is no held-out
test-set result; no figure shows one.
"""
from __future__ import annotations

import json
import subprocess
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
R = "results"
E12_COMMIT = "7f64a72"  # branch worktree-agent-a0df6277cef40ec85 (E12, unmerged)
E12_PATH = "results/2026-10-28_E12_speed_abundance_ptm/metrics.json"

# ---------------------------------------------------------------- palette ---
# dataviz reference palette, light mode, validated (slots 1-3 all-pairs) on #ffffff.
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
VULCAN = "#2a78d6"  # slot 1 blue   - Vulcan (frozen v1 / default scorer)
COMP = "#eb6834"  # slot 2 orange   - comparator (Sage / xiSEARCH / pyAscore)
ALT = "#1baf7a"  # slot 3 aqua     - Vulcan candidate v2_sep / secondary series
BLUE_TINT = "#cde2fb"  # blue ramp step 100 (GPU blocks in schematic)
NEUTRAL_TINT = "#f0efec"  # diverging-neutral grey (CPU blocks)
ST_GOOD, ST_WARN, ST_CRIT = "#0ca30c", "#fab219", "#d03b3b"  # status (fixed)

MM = 1 / 25.4
DOUBLE = 183 * MM  # Nature double-column width

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.5,
    "grid.linestyle": "-",
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "lines.solid_capstyle": "round",
})

# ------------------------------------------------------- source registry ---
SOURCES: "OrderedDict[str, list]" = OrderedDict()
_cache: dict = {}


def _load(path: str):
    if path not in _cache:
        if path.startswith("git:"):
            _, commit, p = path.split(":", 2)
            try:
                txt = subprocess.check_output(["git", "-C", str(ROOT), "show", f"{commit}:{p}"],
                                              stderr=subprocess.DEVNULL)
            except (subprocess.CalledProcessError, FileNotFoundError):
                # public export: the development history is not included; the E12 result files from that
                # commit are shipped in the working tree under the same path
                txt = (ROOT / p).read_text()
            _cache[path] = json.loads(txt)
        else:
            _cache[path] = json.loads((ROOT / path).read_text())
    return _cache[path]


def get(panel: str, label: str, path: str, *keys, scale: float = 1.0):
    """Read one value from a committed JSON file and record its provenance."""
    obj = _load(path)
    for k in keys:
        obj = obj[k]
    val = obj
    shown = [v * scale for v in val] if isinstance(val, list) else (val * scale if scale != 1.0 else val)
    SOURCES.setdefault(panel, []).append((label, path, "/" + "/".join(str(k) for k in keys), val))
    return shown


def text_value(panel: str, label: str, path: str, needle: str, value):
    """A number quoted from committed prose (REPORT.md); asserts the quote exists."""
    txt = (ROOT / path).read_text()
    assert needle in txt, f"{needle!r} not found in {path}"
    SOURCES.setdefault(panel, []).append((label, path, f'text: "{needle}"', value))
    return value


def note(panel: str, label: str, value: str):
    SOURCES.setdefault(panel, []).append((label, "(derived)", value, ""))


# ------------------------------------------------------------- helpers ---
def style_axes(ax, grid_axis=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)


def panel_label(ax, letter, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="right", color=INK)


def dot(ax, x, y, color, hollow=False, marker="o", size=4.2, z=4, **kw):
    ax.plot([x], [y], marker=marker, ms=size, mfc=SURFACE if hollow else color, mec=color if hollow else SURFACE,
            mew=1.1 if hollow else 0.9, color=color, linestyle="none", zorder=z, **kw)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


DEV_TAG = "dev data only · no held-out test-set result"


# ===================================================================== FIG 1
def fig1():
    P = "Fig 1 (schematic)"
    fig = plt.figure(figsize=(DOUBLE, 3.55))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    def box(x, y, w, h, fc, text=None, sub=None, tag=None, fs=6.6, bold=True, ec="none", lw=0, r=0.8):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                    fc=fc, ec=ec, lw=lw, zorder=2))
        if text:
            ty = y + h / 2 + (0.95 if sub else 0)
            ax.text(x + 1.0, ty, text, fontsize=fs, fontweight="bold" if bold else "normal",
                    va="center", ha="left", color=INK, zorder=3)
        if sub:
            ax.text(x + 1.0, y + h / 2 - 1.15, sub, fontsize=5.6, va="center", ha="left", color=INK2, zorder=3)
        if tag:
            ax.text(x + w - 0.9, y + h / 2, tag, fontsize=5.6, va="center", ha="right", color=INK2,
                    fontweight="bold", zorder=3)

    def arrow(x0, y0, x1, y1, color=MUTED, lw=0.9):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=7,
                                    shrinkA=0, shrinkB=0), zorder=1)

    # ---- left: layer stack (DESIGN §3-4 / REPORT §1) ----
    ax.text(1, 48.3, "a  Layers and device placement", fontsize=7.5, fontweight="bold", va="center")
    layers = [
        ("L0  Ingest", "RawFileReader / Bruker readers → Arrow store; persistent decode workers", "CPU", NEUTRAL_TINT),
        ("L0.5  Per-run pilot", "calibration, analyzer/activation, reporter checks", "CPU", NEUTRAL_TINT),
        ("L1  Proteome index", "content-hashed; paired decoys + entrapment; tiled inverted index", "CPU build · GPU resident", BLUE_TINT),
        ("L2  Candidate generation", "K1 CUDA kernel, pair-closed top-K = 500", "GPU", BLUE_TINT),
        ("L3  Exact rescoring", "K2 kernel, label-blind + predicted features, scorer apply, localization", "GPU", BLUE_TINT),
        ("L4  Statistics", "paired double competition, per-stratum TDC, picked protein, FLR", "CPU", NEUTRAL_TINT),
        ("L5  Cross-run", "MS1 features (GPU opt-in), alignment, MBR, LFQ / TMT", "CPU (+GPU opt-in)", NEUTRAL_TINT),
        ("L6  Repository", "best evidence per peptidoform → global paired TDC → picked protein", "CPU", NEUTRAL_TINT),
    ]
    x0, w, h, gap = 1, 52, 4.55, 0.85
    ytop = 45.6
    for i, (name, sub, tag, fc) in enumerate(layers):
        y = ytop - (i + 1) * h - i * gap
        box(x0, y, w, h, fc, name, sub, tag)
        if i < len(layers) - 1:
            arrow(x0 + 5, y, x0 + 5, y - gap)
    leg_y = 1.0
    ax.add_patch(FancyBboxPatch((x0, leg_y), 2.2, 1.6, boxstyle="round,pad=0,rounding_size=0.3", fc=BLUE_TINT, ec="none"))
    ax.text(x0 + 3, leg_y + 0.8, "GPU (A100 80 GB)", fontsize=5.8, va="center", color=INK2)
    ax.add_patch(FancyBboxPatch((x0 + 17, leg_y), 2.2, 1.6, boxstyle="round,pad=0,rounding_size=0.3", fc=NEUTRAL_TINT, ec="none"))
    ax.text(x0 + 20, leg_y + 0.8, "CPU (32 cores, EPYC 9354)", fontsize=5.8, va="center", color=INK2)

    # ---- right: cohort process ----
    cx = 58
    ax.text(cx, 48.3, "b  One long-lived cohort process", fontsize=7.5, fontweight="bold", va="center")
    # RAW input
    box(cx, 39.6, 41, 5.4, NEUTRAL_TINT, "RAW files", "staged to local NVMe (Lustre mmap stalled) or streamed from PRIDE")
    # container
    ax.add_patch(FancyBboxPatch((cx, 11.3), 41, 26.6, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=SURFACE, ec=AXIS, lw=0.8, zorder=1))
    ax.text(cx + 9.0, 36.6, "vulcan.cohort: index + predicted library resident on GPU",
            fontsize=5.8, color=INK2, va="center")
    box(cx + 1.5, 29.0, 38, 5.6, NEUTRAL_TINT, "Decode workers ×24 (CPU)", "L0 decode, L0.5 pilot, MS1 features / hypotheses")
    box(cx + 1.5, 21.0, 38, 5.6, BLUE_TINT, "GPU-stage thread(s)", "L2 search, L3 features + scorer apply (chimeric: secondary search)")
    box(cx + 1.5, 13.0, 38, 5.6, NEUTRAL_TINT, "Post workers ×3 (CPU)", "L4 competition (q), picked protein groups, quant")
    arrow(cx + 6, 39.6, cx + 6, 34.6)
    arrow(cx + 6, 29.0, cx + 6, 26.6)
    arrow(cx + 6, 21.0, cx + 6, 18.6)
    arrow(cx + 6, 11.3, cx + 6, 8.9)
    box(cx, 3.0, 41, 5.9, NEUTRAL_TINT, "Per-run outputs → L6 repository", "global FDR across runs; per-run tables kept for audit")
    fig.text(0.995, -0.02, "Architecture as described in REPORT.md §1 / DESIGN §3-4. Thread counts from E7b/E12 cohort configs.",
             fontsize=5.2, color=MUTED, ha="right", va="bottom")
    note(P, "content", "Qualitative schematic; layer list and device placement from REPORT.md §1 'Approach' and "
         "'CPU/GPU co-design'; worker counts (24 decode, 3 post) from results/2026-10-28_E7b_fair_speed/NOTES.md "
         "§2 cohort setup and E12 NOTES §1 (same design). No data values plotted.")
    save(fig, "fig1_system_overview")


# ===================================================================== FIG 2
def fig2():
    fig = plt.figure(figsize=(DOUBLE, 3.9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1.35, 1], wspace=0.55, hspace=0.75)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    # ---- a: CPU core-s per unit of work, log scale ----
    P = "Fig 2a (core-s per run)"
    E12 = f"git:{E12_COMMIT}:{E12_PATH}"
    H = lambda i: get(P, f"Vulcan E12 headline[{i}] after core-s/run", E12, "headline", i, "after", "core_s_per_run")  # noqa: E731
    rows = []  # (label, unit, vulcan, comparator, comparator_with_conversion, comp_name)
    sage_closed = get(P, "Sage closed D0, 32 thr, core-s/run", E12, "comparators", "sage_closed_D0_per_run", "t32", "core_s")
    trfp = get(P, "TRFP RAW→mzML conversion core-s/run", E12, "comparators", "sage_closed_D0_per_run", "trfp_mzml_conversion_core_s")
    rows.append(("Closed DDA search", "RAW → protein groups, per run", H(0), sage_closed, sage_closed + trfp, "Sage"))
    sage_lfq = get(P, "Sage search+LFQ D9, 32 thr, core-s/run (mzML given)", E12, "comparators", "sage_lfq_D9_per_run", "t32", "core_s")
    rows.append(("Label-free quant (LFQ)", "RAW → protein abundance, per run", H(1), sage_lfq, sage_lfq + trfp, "Sage"))
    sp = get(P, "Sage phospho D0 (2 runs) core-s incl. DB build, no localization", E12, "comparators", "sage_phospho_D0",
             "E12_sage_phospho_D0_t32", 0, "core_s")
    rows.append(("Phospho + localization + FLR", "per run (Sage: no localization)", H(4), sp / 2, None, "Sage"))
    o1 = get(P, "Vulcan open search D0 run 1 core-s", f"{R}/2026-10-27_M4_ptm/metrics.json", "speed_ours", "d0_open", "warm_runs", 0, "core_s")
    o2 = get(P, "Vulcan open search D0 run 2 core-s", f"{R}/2026-10-27_M4_ptm/metrics.json", "speed_ours", "d0_open", "warm_runs", 1, "core_s")
    so = get(P, "Sage open D0 (2 runs), 32 thr core-s", f"{R}/2026-10-27_M4_ptm/metrics.json", "speed_sage", 0, "core_s")
    rows.append(("Open search [−150, +500] Da", "search stage, per run", (o1 + o2) / 2, so / 2, None, "Sage"))
    lv = get(P, "Vulcan semi-specific (LiP) warm per-run core-s median", f"{R}/2026-10-27_M5_lip_cf/metrics.json", "lip", "search", "warm_per_run_core_s_median")
    ls = get(P, "Sage semi LiP arm total core-s (6 runs)", f"{R}/2026-10-27_M5_lip_cf/metrics.json", "lip", "sage_semi", "LiP", "timing", "core_s")
    note(P, "Sage semi per-run divisor", "6 LiP runs (D7 LiP arm; M5_lip_cf NOTES: 446 core-s per run)")
    rows.append(("Semi-specific (LiP)", "search stage, per run", lv, ls / 6, None, "Sage"))
    iv = get(P, "Vulcan IT-CID D4 per-fraction core-s (warm, locked)", f"{R}/2026-10-27_3c_itcid/metrics.json", "runtime", "bench_locked_warm_cache", "D4", "median", "total_core_s")
    isg = get(P, "Sage D4 prod 16 thr core-s (12 fractions, unlocked)", f"{R}/2026-10-27_3c_itcid/metrics.json", "runtime", "comparators_unlocked_16thr", "3c_sage:D4:prod@t16", "core_s")
    note(P, "Sage IT-CID per-fraction divisor", "12 D4 fractions (REPORT §2.1: 'Sage 68 s / 643 core-s for 12 files')")
    rows.append(("Ion-trap CID / TMT", "search + scoring, per fraction", iv, isg / 12, None, "Sage"))
    xv = get(P, "Vulcan XL DSS D3 search (warm) core-s", f"{R}/2026-10-27_M5_xl/metrics.json", "datasets", "D3", "runtime", "search_warm", "core_s")
    xs = get(P, "xiSEARCH D3 core-s (32 thr)", f"{R}/2026-10-27_M5_xl/comparator_runtime.json", "runs", "xiSEARCH_D3", "timing", "core_s")
    rows.append(("Cross-linking (DSS)", "search stage, one run", xv, xs, None, "xiSEARCH"))

    n = len(rows)
    for i, (lab, unit, v, c, cc, cname) in enumerate(rows):
        y = n - 1 - i
        hi = cc if cc else c
        axA.plot([v, hi], [y, y], color=GRID, lw=1.6, zorder=1, solid_capstyle="round")
        dot(axA, v, y, VULCAN, size=5)
        dot(axA, c, y, COMP, size=5)
        if cc:
            dot(axA, cc, y, COMP, hollow=True, size=4.6)
        fold = c / v
        axA.text(hi * 1.6, y, f"{fold:,.0f}×" if fold >= 10 else f"{fold:.0f}×", va="center", ha="left",
                 fontsize=6.3, color=INK, fontweight="bold")
        axA.text(0.02, y + 0.03, lab, va="bottom", ha="right", fontsize=6.6, color=INK,
                 transform=axA.get_yaxis_transform())
        axA.text(0.02, y - 0.05, unit, va="top", ha="right", fontsize=5.5, color=MUTED,
                 transform=axA.get_yaxis_transform())
        if cname != "Sage":
            axA.text(c, y + 0.3, cname, ha="center", va="bottom", fontsize=5.5, color=INK2)
    axA.set_xscale("log")
    axA.set_xlim(0.2, 2e6)
    axA.set_ylim(-0.7, n - 0.3)
    axA.set_yticks([])
    axA.spines["left"].set_visible(False)
    style_axes(axA, "x")
    axA.set_xlabel("CPU core-seconds (user + sys; Vulcan GPU time not included), log scale")
    axA.text(-0.58, 1.075, "CPU cost per unit of work (fold = comparator / Vulcan)", transform=axA.transAxes,
             fontsize=7.5, va="bottom", ha="left")
    handles = [Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=5, label="Vulcan (32 cores + A100)"),
               Line2D([], [], marker="o", ls="none", mfc=COMP, mec=SURFACE, ms=5, label="Comparator (32 or 16 threads)"),
               Line2D([], [], marker="o", ls="none", mfc=SURFACE, mec=COMP, mew=1.1, ms=4.6, label="+ RAW→mzML conversion")]
    axA.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.62, 1.0), ncol=3, fontsize=5.8,
               handletextpad=0.2, columnspacing=0.8, borderaxespad=0.2)
    panel_label(axA, "a", x=-0.62, y=1.07)

    # ---- b: runs/h before -> after (E12) ----
    P = "Fig 2b (runs/h, E12)"
    names = ["Closed", "LFQ, MBR off", "LFQ, MBR on", "Chimeric v2_sep", "Phospho + FLR"]
    for i, nm in enumerate(names):
        b = get(P, f"{nm}: before runs/h", E12, "headline", i, "before", "runs_per_h")
        a = get(P, f"{nm}: after runs/h", E12, "headline", i, "after", "runs_per_h")
        y = len(names) - 1 - i
        axB.plot([b, a], [y, y], color=VULCAN, lw=1.4, alpha=0.35, zorder=1)
        dot(axB, b, y, MUTED, size=4.2)
        dot(axB, a, y, VULCAN, size=5)
        axB.text(a + 180, y, f"{a:,.0f}", va="center", ha="left", fontsize=6.1, color=INK)
    axB.set_yticks(range(len(names)))
    axB.set_yticklabels(names[::-1])
    axB.set_xlim(0, 7000)
    axB.set_xticks([0, 2000, 4000, 6000])
    axB.set_xticklabels(["0", "2,000", "4,000", "6,000"])
    axB.tick_params(axis="y", length=0)
    axB.spines["left"].set_visible(False)
    style_axes(axB, "x")
    axB.set_xlabel("RAW → output, runs per hour (96-run cohorts)")
    axB.set_title("Throughput before → after exact speed-ups (E12)", loc="left")
    axB.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=MUTED, mec=SURFACE, ms=4.2, label="before"),
                        Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=5, label="after (opt-in, unmerged)")],
               loc="lower left", bbox_to_anchor=(-0.02, 0.98), ncol=2, fontsize=5.8, handletextpad=0.2, borderaxespad=0.1)
    axB.set_title("Throughput before → after exact speed-ups (E12)", loc="left", pad=14)
    panel_label(axB, "b", x=-0.36, y=1.1)

    # ---- c: 1 TB LFQ estimate ----
    P = "Fig 2c (1 TB LFQ estimate)"
    v_lo = text_value(P, "Vulcan wall time low (h)", "REPORT.md", "**about 0.8–1 h**", 0.8)
    v_hi = text_value(P, "Vulcan wall time high (h)", "REPORT.md", "**about 0.8–1 h**", 1.0)
    s_lo = text_value(P, "Sage + conversion wall time low (h)", "REPORT.md", "**about 4.5–7 h**", 4.5)
    s_hi = text_value(P, "Sage + conversion wall time high (h)", "REPORT.md", "**about 4.5–7 h**", 7.0)
    text_value(P, "assumption: ~1 GB per run, ~1,000 runs", "REPORT.md", "Assumes about 1 GB per run, so about 1,000 runs.", "")
    for y, lo, hi, col, lab in [(1, v_lo, v_hi, VULCAN, "Vulcan"), (0, s_lo, s_hi, COMP, "Sage + mzML conversion")]:
        axC.plot([lo, hi], [y, y], color=col, lw=6, solid_capstyle="round", zorder=2)
        axC.text(hi + 0.25, y, f"{lo:g}–{hi:g} h", va="center", ha="left", fontsize=6.3, color=INK)
    axC.set_yticks([0, 1])
    axC.set_yticklabels(["Sage + RAW→mzML\n(32 cores, CPU only)", "Vulcan\n(32 cores + A100)"])
    axC.tick_params(axis="y", length=0)
    axC.spines["left"].set_visible(False)
    axC.set_ylim(-0.6, 1.6)
    axC.set_xlim(0, 8.5)
    style_axes(axC, "x")
    axC.set_xlabel("Estimated wall time for ~1,000 RAW files (h)")
    axC.set_title("ESTIMATE: 1 TB LFQ on one machine", loc="left")
    axC.text(1.0, -0.42, "From per-run costs, not a measured end-to-end run", transform=axC.transAxes,
             ha="right", va="top", fontsize=5.5, color=MUTED)
    panel_label(axC, "c", x=-0.36)

    fig.text(0.01, -0.02, f"{DEV_TAG} · shared node (load 14–110), runs/h vary ±25% between rounds",
             fontsize=5.5, color=MUTED, ha="left")
    save(fig, "fig2_speed")


# ===================================================================== FIG 3
def fig3():
    P = "Fig 3 (identification ratios)"
    T1 = f"{R}/2026-10-27_T1driver_dryrun"
    E7B = f"{R}/2026-10-28_E7b_fair_speed/t1"
    IT = f"{R}/2026-10-27_3c_itcid/followup/ratios.json"
    M7 = f"{R}/2026-10-28_M7_encoders/metrics.json"

    def nonin(path, comp, level, lab):
        r = get(P, f"{lab} ratio", f"{path}/metrics.json", "noninferiority", comp, level, "matched_fdp", "0.01", "ratio")
        ci = get(P, f"{lab} 95% CI", f"{path}/metrics.json", "noninferiority", comp, level, "matched_fdp", "0.01", "ci95")
        return r, ci

    groups = []
    g = ("Standard DDA · frozen v1", "vs Sage + predicted features + mokapot", [])
    g[2].append(("D0", "peptidoforms", *nonin(f"{T1}/D0", "sage+pred+mokapot", "peptidoform", "v1 D0 peptidoform"), "v1"))
    g[2].append(("D1", "peptidoforms", *nonin(f"{T1}/D1", "sage+pred+mokapot", "peptidoform", "v1 D1 peptidoform"), "v1"))
    groups.append(g)
    g = ("Chimeric · v2_sep (candidate)", "vs Sage chimera + predicted features + mokapot", [])
    g[2].append(("D0 (R = 5)", "peptidoforms", *nonin(f"{E7B}/D0", "sage_chim5+pred+mokapot", "peptidoform", "v2_sep D0 peptidoform vs R=5"), "v2"))
    g[2].append(("D1 (R = 2)", "peptidoforms", *nonin(f"{E7B}/D1", "sage_chim2+pred+mokapot", "peptidoform", "v2_sep D1 peptidoform vs R=2"), "v2"))
    g[2].append(("D8h dev holdout (R = 2)", "peptidoforms", *nonin(f"{E7B}/D8h", "sage_chim2+pred+mokapot", "peptidoform", "v2_sep D8h peptidoform"), "v2"))
    g[2].append(("D8h dev holdout (R = 2)", "protein groups", *nonin(f"{E7B}/D8h", "sage_chim2+pred+mokapot", "protein", "v2_sep D8h protein"), "v2p"))
    groups.append(g)
    g = ("Ion-trap CID / TMT · itcid_linear_pred", "vs Sage + Prosit features + mokapot", [])
    for ds in ("D4", "D10cf"):
        r = get(P, f"IT-CID {ds} sequence ratio", IT, "datasets", ds, "ratios", "sage+pred+mokapot", "sequence", "matched_fdp", "0.01", "ratio")
        ci = get(P, f"IT-CID {ds} sequence 95% CI", IT, "datasets", ds, "ratios", "sage+pred+mokapot", "sequence", "matched_fdp", "0.01", "ci95")
        g[2].append((ds, "sequences", r, ci, "v1"))
    groups.append(g)
    g = ("Cross-encoder ce_s (research arm) · NO-GO", "vs frozen v1 (static_gbm_pred), not vs Sage", [])
    r = get(P, "ce_s / v1 D8 ratio", M7, "E6", "D8", "ratio_ce_vs", "static_gbm_pred", "P_0.01", "ratio")
    ci = get(P, "ce_s / v1 D8 95% CI", M7, "E6", "D8", "ratio_ce_vs", "static_gbm_pred", "P_0.01", "ci95")
    g[2].append(("D8h dev holdout", "peptidoforms", r, ci, "nogo"))
    groups.append(g)

    nrows = sum(len(g[2]) for g in groups) + 2 * len(groups)
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.62, 3.9))
    y = nrows
    ylabels = []
    for title, sub, rows in groups:
        y -= 0.35
        ax.text(-0.02, y, title, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.7, fontweight="bold", color=INK)
        y -= 0.62
        ax.text(-0.02, y, sub, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.7, color=MUTED)
        y -= 0.03
        for ds, unit, r, ci, kind in rows:
            y -= 1
            col = {"v1": VULCAN, "v2": ALT, "v2p": ALT, "nogo": MUTED}[kind]
            ax.plot(ci, [y, y], color=col, lw=1.6, zorder=2, solid_capstyle="round")
            dot(ax, r, y, col, hollow=(kind == "v2p"), size=5, z=4)
            ax.text(max(ci[1], r) + 0.008, y, f"{r:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]", va="center", ha="left",
                    fontsize=5.8, color=INK2)
            ax.text(-0.02, y, f"{ds} · {unit}", transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=6.3, color=INK)
            ylabels.append(y)
        y -= 0.4
    ax.axvline(1.0, color=INK2, lw=0.8, zorder=1)
    ax.axvline(0.97, color=MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=1)
    ax.text(1.0, nrows - 0.02, "parity", ha="center", va="bottom", fontsize=5.8, color=INK2)
    ax.text(0.97, nrows - 0.02, "non-inferiority\nmargin 0.97", ha="right", va="bottom", fontsize=5.3, color=MUTED)
    ax.set_xlim(0.93, 1.31)
    ax.set_ylim(y + 0.3, nrows + 0.9)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    style_axes(ax, "x")
    ax.set_xlabel("Ratio of identifications, Vulcan / comparator,\nat matched entrapment-estimated FDP 1% (arm P), 95% CI")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=5, label="default scorer (frozen / per class)"),
                       Line2D([], [], marker="o", ls="none", mfc=ALT, mec=SURFACE, ms=5, label="chimeric candidate v2_sep"),
                       Line2D([], [], marker="o", ls="none", mfc=SURFACE, mec=ALT, mew=1.1, ms=5, label="protein groups (hollow)"),
                       Line2D([], [], marker="o", ls="none", mfc=MUTED, mec=SURFACE, ms=5, label="NO-GO research arm")],
              loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=5.6, handletextpad=0.2, borderaxespad=0.1)
    fig.text(0.0, -0.17, f"{DEV_TAG} · D0/D1 are selection sets · FDR claims not gate-validated",
             fontsize=5.5, color=MUTED, ha="left", transform=fig.transFigure)
    save(fig, "fig3_sensitivity")


# ===================================================================== FIG 4
def fig4():
    fig = plt.figure(figsize=(DOUBLE, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.95, 1.55], wspace=0.38)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---- a: per-run FDP-hat at nominal 1% ----
    P = "Fig 4a (per-run entrapment FDP at nominal 1%)"
    T1 = f"{R}/2026-10-27_T1driver_dryrun"
    E7 = f"{R}/2026-10-28_E7_chimeric_v2/t1"
    sets = [("D0", {"v1": (f"{T1}/D0", "vulcan"), "v2": (f"{E7}/D0", "vulcan_v2_sep"), "sage": (f"{T1}/D0", "sage+pred+mokapot")}),
            ("D1", {"v1": (f"{T1}/D1", "vulcan"), "v2": (f"{E7}/D1", "vulcan_v2_sep"), "sage": (f"{T1}/D1", "sage+pred+mokapot")}),
            ("D8h", {"v1": (f"{E7}/D8h", "vulcan"), "v2": (f"{E7}/D8h", "vulcan_v2_sep"), "sage": (f"{E7}/D8h", "sage+pred+mokapot")})]
    cols = {"v1": VULCAN, "v2": ALT, "sage": COMP}
    offs = {"v1": -0.22, "v2": 0.0, "sage": 0.22}
    for xi, (ds, m) in enumerate(sets):
        for k, (path, meth) in m.items():
            base = ("methods", meth, "verdicts", "peptidoform", "0.01")
            fp = get(P, f"{ds} {meth} arm P mean FDP", f"{path}/metrics.json", *base, "P", "a", "mean_fdp", scale=100)
            ci = get(P, f"{ds} {meth} arm P 95% CI", f"{path}/metrics.json", *base, "P", "a", "ci95", scale=100)
            fc = get(P, f"{ds} {meth} arm C mean FDP", f"{path}/metrics.json", *base, "C", "a", "mean_fdp", scale=100)
            x = xi + offs[k]
            axA.plot([x, x], ci, color=cols[k], lw=1.5, zorder=2)
            dot(axA, x, fp, cols[k], size=4.8)
            dot(axA, x, fc, cols[k], hollow=True, marker="D", size=3.4)
    axA.axhline(1.0, color=INK2, lw=0.8, zorder=1)
    axA.text(1.48, 1.005, "nominal", ha="center", va="bottom", fontsize=5.6, color=INK2)
    axA.set_xticks(range(3))
    axA.set_xticklabels(["D0\n(2 runs)", "D1\n(1 run)", "D8h dev holdout\n(10 runs)"])
    axA.set_xlim(-0.5, 2.5)
    axA.set_ylim(0.6, 1.3)
    axA.set_ylabel("Entrapment-estimated FDP (%)")
    style_axes(axA, "y")
    axA.set_title("Per-run lists at nominal 1%", loc="left")
    axA.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=4.8, label="Vulcan v1"),
                        Line2D([], [], marker="o", ls="none", mfc=ALT, mec=SURFACE, ms=4.8, label="v2_sep"),
                        Line2D([], [], marker="o", ls="none", mfc=COMP, mec=SURFACE, ms=4.8, label="Sage+pred+mokapot"),
                        Line2D([], [], marker="D", ls="none", mfc=SURFACE, mec=INK2, mew=1.0, ms=3.4, label="arm C (natural)")],
               loc="upper left", fontsize=5.4, ncol=2, handletextpad=0.15, columnspacing=0.6, borderaxespad=0.1)
    panel_label(axA, "a", x=-0.2)

    # ---- b: E9 global vs naive union ----
    P = "Fig 4b (L6 global FDR vs naive union, E9)"
    tab = _load(f"{R}/2026-10-27_E9_scale/metrics.json")["e9_table"]
    series = {}
    for idx, row in enumerate(tab):
        if row["score"] == "score_raw" and row["alpha"] == 0.01 and row["arm"] == "P":
            for kind in ("global", "naive"):
                key = (row["level"], kind)
                v = get(P, f"N={row['N']} {row['level']} {kind} arm P mean FDP (score_raw = static_gbm_pred)",
                        f"{R}/2026-10-27_E9_scale/metrics.json", "e9_table", idx, kind, "mean_fdp", scale=100)
                series.setdefault(key, []).append((row["N"], v))
    style = {("peptidoform", "global"): (VULCAN, "o", "L6 global · peptide"),
             ("protein_group", "global"): (VULCAN, "s", "L6 global · protein"),
             ("peptidoform", "naive"): (COMP, "o", "naive union · peptide"),
             ("protein_group", "naive"): (COMP, "s", "naive union · protein")}
    for key, pts in series.items():
        pts.sort()
        col, mk, lab = style[key]
        xs, ys = zip(*pts)
        axB.plot(xs, ys, color=col, lw=1.6, zorder=2, label=lab, marker=mk, ms=4.2, mec=SURFACE, mew=0.8,
                 mfc=col if key[0] == "peptidoform" else SURFACE)
        if key[0] == "protein_group":
            axB.lines[-1].set_markeredgecolor(col)
            axB.lines[-1].set_markeredgewidth(1.0)
    axB.axhline(1.0, color=INK2, lw=0.8, zorder=1)
    axB.set_xscale("log")
    axB.set_xticks([1, 10, 100, 1000])
    axB.set_xticklabels(["1", "10", "100", "1,000"])
    axB.set_xlabel("Runs combined, N (D8, dev)")
    axB.set_ylabel("Entrapment-estimated FDP (%), arm P")
    axB.set_ylim(0, 21)
    style_axes(axB, "y")
    axB.set_title("Repository scale (L6), 1%", loc="left")
    pep = dict(series[("peptidoform", "naive")])
    prot = dict(series[("protein_group", "naive")])
    axB.text(1000, prot[1000] + 0.9, f"{prot[1000]:.1f}%", ha="center", fontsize=5.8, color=INK)
    axB.text(1000, pep[1000] + 0.9, f"{pep[1000]:.1f}%", ha="center", fontsize=5.8, color=INK)
    axB.text(1150, 1.0 + 0.35, "≈1.0–1.2%", ha="right", va="bottom", fontsize=5.6, color=INK)
    axB.legend(loc="upper left", fontsize=5.4, handletextpad=0.3, borderaxespad=0.1)
    panel_label(axB, "b", x=-0.2)

    # ---- c: chimeric per-stratum ----
    P = "Fig 4c (chimeric per-stratum FDP at 1%)"
    E7S = f"{R}/2026-10-28_E7_chimeric_v2"
    E7B = f"{R}/2026-10-28_E7b_fair_speed"
    blocks = [("Group-walk\nD0", f"{E7S}/strata_D0_vulcan_chim_v1.json", MUTED),
              ("v2_sep\nD0", f"{E7S}/strata_D0_vulcan_v2_sep.json", ALT),
              ("v2_sep\nD8h", f"{E7S}/strata_D8h_vulcan_v2_sep.json", ALT)]
    x = 0
    ticks, ticklabs, groupcent = [], [], []
    for name, path, col in blocks:
        xs0 = x
        for s, sl in (("1", "1"), ("2", "2"), (">=3", "≥3")):
            v = get(P, f"{name.replace(chr(10), ' ')} stratum {sl} arm P FDP", path, "by_stratum", s, "verdicts", "0.01", "P", "a", "mean_fdp", scale=100)
            ci = get(P, f"{name.replace(chr(10), ' ')} stratum {sl} arm P 95% CI", path, "by_stratum", s, "verdicts", "0.01", "P", "a", "ci95", scale=100)
            ver = get(P, f"{name.replace(chr(10), ' ')} stratum {sl} verdict (a)", path, "by_stratum", s, "verdicts", "0.01", "overall_a")
            axC.plot([x, x], ci, color=col, lw=1.5, zorder=2)
            dot(axC, x, v, col, size=4.8, marker="X" if ver == "violated" else "o")
            if ver == "violated":
                axC.text(x, ci[1] + 0.5, "violated", ha="center", fontsize=5.6, color=INK, fontweight="bold")
            ticks.append(x)
            ticklabs.append(sl)
            x += 1
        groupcent.append(((xs0 + x - 1) / 2, name))
        x += 0.9
    for ds, lab in (("D0", "Sage\nD0"), ("D8h", "Sage\nD8h")):
        meth = "sage_chim5+pred+mokapot" if ds == "D0" else "sage_chim2+pred+mokapot"
        path = f"{E7B}/sage_rank_strata_{ds}.json"
        xs0 = x
        for s, sl in (("rank_1", "1"), ("rank_ge2", "≥2")):
            v = get(P, f"{lab.replace(chr(10), ' ')} {sl} arm P FDP", path, "methods", meth, s, "0.01", "mean_fdp", scale=100)
            ci = get(P, f"{lab.replace(chr(10), ' ')} {sl} 95% CI", path, "methods", meth, s, "0.01", "ci95", scale=100)
            ver = get(P, f"{lab.replace(chr(10), ' ')} {sl} verdict", path, "methods", meth, s, "0.01", "verdict")
            axC.plot([x, x], ci, color=COMP, lw=1.5, zorder=2)
            dot(axC, x, v, COMP, size=4.8, marker="X" if ver == "violated" else "o")
            if ver == "violated":
                axC.text(x, ci[1] + 0.5, "violated", ha="center", fontsize=5.6, color=INK, fontweight="bold")
            elif v > 3:
                axC.text(x, ci[1] + 0.5, "not\nevaluable", ha="center", fontsize=5.0, color=INK2, linespacing=0.9)
            ticks.append(x)
            ticklabs.append(sl)
            x += 1
        groupcent.append(((xs0 + x - 1) / 2, lab))
        x += 0.9
    axC.axhline(1.0, color=INK2, lw=0.8, zorder=1)
    axC.set_xticks(ticks)
    axC.set_xticklabels(ticklabs)
    for cx, name in groupcent:
        axC.text(cx, -0.13, name, transform=axC.get_xaxis_transform(), ha="center", va="top", fontsize=5.2,
                 color=INK2, linespacing=0.95)
    axC.text(-0.9, -0.075, "rank", transform=axC.get_xaxis_transform(), ha="right", va="top", fontsize=5.5, color=MUTED)
    axC.set_xlim(-0.8, x - 0.6)
    axC.set_ylim(0, 12.5)
    axC.set_ylabel("Arm-P FDP within stratum (%)")
    style_axes(axC, "y")
    axC.set_title("Chimeric strata, nominal 1%", loc="left")
    axC.set_xlabel(" ")
    panel_label(axC, "c", x=-0.12)

    fig.text(0.01, -0.13, f"{DEV_TAG} · entrapment-estimated (not verified) FDP · no null-spectrum gate passed",
             fontsize=5.5, color=MUTED, ha="left")
    save(fig, "fig4_error_control")


# ===================================================================== FIG 5
def fig5():
    fig = plt.figure(figsize=(DOUBLE, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.45)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
    T2 = f"{R}/2026-10-28_T2_driver_dryrun/D2/metrics.json"

    # ---- a: true FLR vs nominal with one-sided 95% UCB ----
    P = "Fig 5a (Ala-decoy FLR, true vs nominal)"
    cells = [("PXD000138", "0.01"), ("PXD009449", "0.01"), ("PXD000138", "0.05"), ("PXD009449", "0.05")]
    for i, (ds, a) in enumerate(cells):
        nominal = float(a) * 100
        for k, off, col, hollow in (("ours_A", -0.13, VULCAN, False), ("ours_L", 0.13, VULCAN, True)):
            tf = get(P, f"{ds} {k} @{a} true FLR", T2, "datasets", ds, "t2", k, a, "true_flr", scale=100)
            ucb = get(P, f"{ds} {k} @{a} one-sided 95% UCB", T2, "datasets", ds, "t2", k, a, "ucb", scale=100)
            ver = get(P, f"{ds} {k} @{a} verdict", T2, "datasets", ds, "t2", k, a, "verdict")
            x = i + off
            axA.plot([x, x], [tf, ucb], color=col, lw=1.4, alpha=0.55 if hollow else 1, zorder=2)
            axA.plot([x - 0.06, x + 0.06], [ucb, ucb], color=col, lw=1.2, alpha=0.55 if hollow else 1, zorder=2)
            dot(axA, x, tf / nominal * nominal, col, hollow=hollow, size=4.8)
            if k == "ours_A":
                axA.text(x - 0.08, ucb, {"PASS": "pass", "inconclusive": "inconcl."}.get(ver, ver),
                         ha="right", va="center", fontsize=5.2, color=INK2)
        axA.plot([i - 0.38, i + 0.38], [nominal, nominal], color=INK2, lw=0.9, zorder=1)
    axA.set_yscale("log")
    axA.set_ylim(0.15, 12)
    axA.set_yticks([0.2, 0.5, 1, 2, 5, 10])
    axA.set_yticklabels(["0.2", "0.5", "1", "2", "5", "10"])
    axA.set_xticks(range(4))
    axA.set_xticklabels(["PXD\n000138\n1%", "PXD\n009449\n1%", "PXD\n000138\n5%", "PXD\n009449\n5%"])
    axA.set_xlim(-0.55, 3.55)
    axA.set_ylabel("True FLR of accepted sites (%)")
    style_axes(axA, "y")
    axA.set_title("True vs nominal site FLR", loc="left")
    axA.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=4.8, label="Ala decoy (primary)"),
                        Line2D([], [], marker="o", ls="none", mfc=SURFACE, mec=VULCAN, mew=1.1, ms=4.8, label="Leu decoy (secondary)"),
                        Line2D([], [], color=INK2, lw=0.9, label="nominal"),
                        Line2D([], [], color=VULCAN, lw=1.2, label="to one-sided 95% UCB")],
               loc="lower right", fontsize=5.2, handletextpad=0.3, borderaxespad=0.1)
    panel_label(axA, "a", x=-0.2)

    # ---- b: recall difference vs pyAscore ----
    P = "Fig 5b (site recall vs pyAscore, paired bootstrap)"
    cells_b = [("PXD000138", "true_0.01", "PXD\n000138\n1%"), ("PXD009449", "true_0.01", "PXD\n009449\n1%"),
               ("PXD000138", "true_0.05", "PXD\n000138\n5%"), ("PXD009449", "true_0.05", "PXD\n009449\n5%")]
    for i, (ds, key, lab) in enumerate(cells_b):
        d = get(P, f"{ds} {key} recall diff (ours − pyAscore), Ala", T2, "datasets", ds, "t2", "diff_A", key, "diff", scale=100)
        ci = get(P, f"{ds} {key} 95% CI", T2, "datasets", ds, "t2", "diff_A", key, "ci95", scale=100)
        ver = get(P, f"{ds} {key} verdict", T2, "datasets", ds, "t2", "diff_A", key, "verdict")
        get(P, f"{ds} {key} ours recall", T2, "datasets", ds, "t2", "diff_A", key, "ours")
        get(P, f"{ds} {key} pyAscore recall", T2, "datasets", ds, "t2", "diff_A", key, "pyascore")
        axB.plot([i, i], ci, color=VULCAN, lw=1.5, zorder=2)
        dot(axB, i, d, VULCAN, size=4.8)
        axB.text(i + 0.12, d, {"no difference shown": "n.d."}.get(ver, ver), va="center", ha="left", fontsize=5.5, color=INK2)
    axB.axhline(0, color=INK2, lw=0.8, zorder=1)
    axB.set_xticks(range(4))
    axB.set_xticklabels([c[2] for c in cells_b], fontsize=5.8)
    axB.set_xlim(-0.5, 3.8)
    axB.set_ylabel("Site recall difference, Vulcan − pyAscore\n(percentage points, 95% CI)")
    style_axes(axB, "y")
    axB.set_title("Recall at matched true FLR", loc="left")
    panel_label(axB, "b", x=-0.24)

    # ---- c: LFQ accuracy vs Sage ----
    P = "Fig 5c (LFQ accuracy vs Sage, D9 HYE, MBR off)"
    M6D = f"{R}/2026-10-27_M6d_quant_accuracy/metrics.json"
    base = ("before_after", "d9_all_m6d_off", "vs_sage", "matched_this_version")
    labs = [("all", "All"), ("HUMAN", "Human\n(1:1)"), ("YEAST", "Yeast\n(2:1)"), ("ECOLI", "E. coli\n(1:4)")]
    for i, (sp, lab) in enumerate(labs):
        if sp == "all":
            o = get(P, "Vulcan median |eps| matched set", M6D, *base, "ours_eps")
            s = get(P, "Sage median |eps| matched set", M6D, *base, "sage_eps")
        else:
            o = get(P, f"Vulcan median |eps| {sp}", M6D, *base, "ours_eps_species", sp)
            s = get(P, f"Sage median |eps| {sp}", M6D, *base, "sage_eps_species", sp)
        axC.plot([i - 0.12, i + 0.12], [o, s], color=GRID, lw=1.2, zorder=1)
        dot(axC, i - 0.12, o, VULCAN, size=5)
        dot(axC, i + 0.12, s, COMP, size=5)
    get(P, "matched peptidoforms n", M6D, *base, "n")
    axC.set_xticks(range(4))
    axC.set_xticklabels([l for _, l in labs])
    axC.set_xlim(-0.5, 3.5)
    axC.set_ylim(0.09, 0.175)
    axC.set_ylabel("Median |log2 ratio error|")
    style_axes(axC, "y")
    axC.set_title("LFQ accuracy, HYE (D9)", loc="left")
    axC.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=VULCAN, mec=SURFACE, ms=5, label="Vulcan (MBR off)"),
                        Line2D([], [], marker="o", ls="none", mfc=COMP, mec=SURFACE, ms=5, label="Sage LFQ")],
               loc="upper left", fontsize=5.6, handletextpad=0.2, borderaxespad=0.1)
    axC.text(1.0, 0.02, "lower is better", transform=axC.transAxes, ha="right", va="bottom", fontsize=5.4, color=MUTED)
    panel_label(axC, "c", x=-0.24)

    fig.text(0.01, -0.15, f"{DEV_TAG} · PXD000138/PXD009449 settings tuned on dev · T2 sealed, not run",
             fontsize=5.5, color=MUTED, ha="left")
    save(fig, "fig5_ptm_quant")


# ===================================================================== FIG 6
def fig6():
    P = "Fig 6 (statistical-gate history)"
    # Events: (decision id, date) in decision-log order; x is the ordinal position.
    events = ["D3\n10-26", "D9\n10-27", "D12\n10-27", "D15\n10-27", "D16\n10-27", "D20\n10-28", "D22\n10-28",
              "D24\n10-28", "D28\n10-28", "D29\n10-28", "D30\n10-28", "D32\n10-28", "D34\n10-28", "D36\n10-28",
              "D38\n10-30"]
    X = {e.split("\n")[0]: i for i, e in enumerate(events)}
    PROC, GOOD, WARN, FAIL = "proc", "good", "warn", "fail"
    rows = [
        ("Pre-registration v0.4 (DESIGN §6)", [("D3", PROC, "frozen")], "dev/sealed split fixed", PROC),
        ("A1: tail-G2 / G1 / G3 on v1", [("D9", PROC, "adopted"), ("D12", FAIL, "")], "FAILED", FAIL),
        ("A2 / DN2 fresh null spectra", [("D15", PROC, "adopted"), ("D16", PROC, "bug, re-walk"), ("D24", WARN, "")], "not evaluable (6 of ≥10 projects)", WARN),
        ("DNIT: ion-trap CID nulls", [("D20", WARN, "")], "futile (581 vs 120k spectra)", WARN),
        ("A3: v2_sep claimable", [("D22", PROC, "adopted"), ("D24", WARN, "")], "DN2 tiers expected not evaluable", WARN),
        ("A4 / DN3 new null set", [("D28", WARN, "")], "not run (P(v1 min) 0.08–0.43)", WARN),
        ("C9 pilot (simulation, not a gate)", [("D30", GOOD, "")], "GO (P = 0.97–1.00)", GOOD),
        ("A5: relaxed claim conditions", [("D29", PROC, "adopted")], "single T1+T5 look allowed", PROC),
        ("A6 / HN1 human removed-peptide nulls", [("D32", PROC, "adopted"), ("D34", WARN, "")], "not tagged: infeasible (P = 0.63 < 0.8)", WARN),
        ("Sealed look T5 + T1", [("D36", PROC, "opened 14:23 CET"), ("D38", FAIL, "")], "failed (unterminated), unread", FAIL),
    ]
    for r in rows:
        note(P, r[0], f"events {[e[0] for e in r[1]]}; final status '{r[2]}' — RESEARCH_TREE.md Decisions and REPORT.md §4")
    text_value(P, "A1 FAILED", "REPORT.md", "v1 then **failed** tail-G2, G1 and G3", "FAILED")
    text_value(P, "A2 not evaluable", "REPORT.md", "so A2 is **not evaluable**", "not evaluable")
    text_value(P, "A4 not run", "REPORT.md", "**A4 not run**", "not run")
    text_value(P, "A6 not tagged", "REPORT.md", "**Not tagged: infeasible**", "not tagged: infeasible")
    text_value(P, "look failed", "REPORT.md", "'failed (unterminated)'", "failed (unterminated)")
    text_value(P, "C9 GO", "REPORT.md", "human removed-peptide nulls GO, P(v1 minimum) 0.97-1.00", "GO")

    mk = {PROC: ("s", MUTED, 3.6), GOOD: ("o", ST_GOOD, 6.5), WARN: ("D", ST_WARN, 5.6), FAIL: ("X", ST_CRIT, 7.5)}
    fig, ax = plt.subplots(figsize=(DOUBLE, 3.1))
    n = len(rows)
    for i, (name, evs, status, final) in enumerate(rows):
        y = n - 1 - i
        xs = [X[e[0]] for e in evs]
        if len(xs) > 1:
            ax.plot([min(xs), max(xs)], [y, y], color=AXIS, lw=1.2, zorder=1)
        for e, kind, lab in evs:
            m, c, s = mk[kind]
            ax.plot(X[e], y, marker=m, ms=s, color=c, mec=SURFACE, mew=0.8, ls="none", zorder=3)
            if lab:
                ax.text(X[e], y + (-0.3 if lab.startswith("bug") else 0.28), lab, va="top" if lab.startswith("bug") else "bottom", ha="center", fontsize=5.0, color=MUTED)
        ax.text(-0.01, y, name, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.4, color=INK)
        icon = {FAIL: "✕ ", WARN: "◆ ", GOOD: "● ", PROC: "■ "}[final]
        ax.text(1.01, y, icon, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.4,
                color={FAIL: ST_CRIT, WARN: ST_WARN, GOOD: ST_GOOD, PROC: MUTED}[final])
        ax.text(1.03, y, status, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.3,
                color=INK, fontweight="bold" if final in (FAIL, WARN) else "normal")
    ax.set_xticks(range(len(events)))
    ax.set_xticklabels(events, fontsize=5.6)
    ax.set_xlim(-0.6, len(events) - 0.4)
    ax.set_ylim(-0.7, n - 0.2)
    ax.set_yticks([])
    for s in ("left", "top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Decision-log entry (RESEARCH_TREE.md) and date, 2026 — ordinal, not to time scale")
    ax.text(1.01, n - 0.15, "Outcome (literal)", transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=6.6, fontweight="bold", color=INK)
    ax.legend(handles=[Line2D([], [], marker="X", ls="none", color=ST_CRIT, mec=SURFACE, ms=7, label="failed"),
                       Line2D([], [], marker="D", ls="none", color=ST_WARN, mec=SURFACE, ms=5.6, label="not evaluable / not run / infeasible"),
                       Line2D([], [], marker="o", ls="none", color=ST_GOOD, mec=SURFACE, ms=6, label="GO (simulation only)"),
                       Line2D([], [], marker="s", ls="none", color=MUTED, mec=SURFACE, ms=3.6, label="process step")],
              loc="upper center", bbox_to_anchor=(0.45, -0.25), ncol=4, fontsize=5.8, handletextpad=0.2, columnspacing=1.2)
    fig.text(0.0, -0.2, "No null-spectrum gate was passed. No held-out test-set (T1–T5) number has been seen. "
             "FDR claims for v1, v2_sep and IT-CID are not gate-supported.", fontsize=5.6, color=INK2, ha="left")
    save(fig, "fig6_gate_history")


# ============================================================= SOURCES.md
def write_sources():
    lines = ["# Figure data sources", "",
             "Generated by `reports/figures/make_figures.py` — do not edit by hand. Every plotted number is read at "
             "run time from the file and JSON key listed here (keys as JSON-pointer-like paths; list indices are "
             "integers). `git:7f64a72:<path>` means the file as committed in E12 commit 7f64a72 on branch "
             "`worktree-agent-a0df6277cef40ec85` (unmerged). Values are raw (fractions, not %) unless the label says "
             "otherwise. All data are **dev**; no sealed file was read and no held-out test-set result exists.", ""]
    for panel, items in SOURCES.items():
        lines += [f"## {panel}", "", "| quantity | file | key | value |", "|---|---|---|---|"]
        for lab, path, key, val in items:
            if isinstance(val, float):
                v = f"{val:.6g}"
            elif isinstance(val, list):
                v = "[" + ", ".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in val) + "]"
            else:
                v = str(val)
            lines.append(f"| {lab} | `{path}` | `{key}` | {v} |")
        lines.append("")
    (OUT / "SOURCES.md").write_text("\n".join(lines))


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    write_sources()
    print("wrote", sorted(p.name for p in OUT.glob("fig*")))
