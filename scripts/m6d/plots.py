#!/usr/bin/env python3
"""M6d diagnostic plots from OUT/diag_d9.json and OUT/diag_tables (PNG, light theme).
    scripts/m2b/vpy.sh scripts/m6d/plots.py OUT"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

OUT = sys.argv[1]
P = os.path.join(OUT, "plots")
os.makedirs(P, exist_ok=True)
OURS, SAGE, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#a3a29c", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#ecebe7", "grid.linewidth": 0.6, "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb"})
d = json.load(open(os.path.join(OUT, "diag_d9.json")))
dec = d["decomposition"]

# 1. decomposition of the matched gap
items = [("Sage LFQ", dec["sage"]["eps"], SAGE), ("Sage, masked to our missing values", dec["sage_masked_to_our_missingness"]["eps"], SAGE),
         ("Vulcan M6", dec["ours_M6"]["eps"], OURS), ("Vulcan M6, dominant charge only", dec["ours_charge_rule_dominant"]["eps"], OURS)]
fig, ax = plt.subplots(figsize=(6.4, 2.4))
for i, (lab, v, c) in enumerate(items[::-1]):
    ax.barh(i, v, color=c, height=0.6)
    ax.text(v + 0.002, i, f"{v:.3f}", va="center", color=INK)
ax.set_yticks(range(len(items)), [x[0] for x in items[::-1]])
ax.set_xlim(0.09, 0.145)
ax.set_xlabel("median |epsilon| on 14,154 matched D9 peptidoforms")
ax.set_title("Where the matched-peptide gap comes from (D9, MBR off)", loc="left", color=INK)
fig.tight_layout(); fig.savefig(os.path.join(P, "gap_decomposition.png"), dpi=150); plt.close(fig)

# 2. intensity dependence
t = pd.read_csv(os.path.join(OUT, "diag_tables", "species_by_intensity.tsv"), sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(8.4, 2.6), sharey=True)
for ax, sp in zip(axs, ("HUMAN", "YEAST", "ECOLI")):
    x = t[t.species == sp]
    ax.plot(x.sage_intensity_quintile, x.ours_median_abs_eps, "-o", color=OURS, lw=2, ms=5, label="Vulcan M6")
    ax.plot(x.sage_intensity_quintile, x.sage_median_abs_eps, "-o", color=SAGE, lw=2, ms=5, label="Sage LFQ")
    ax.set_title(sp.capitalize(), loc="left", color=INK)
    ax.set_xlabel("Sage intensity quintile")
axs[0].set_ylabel("median |epsilon|")
axs[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(P, "eps_by_intensity.png"), dpi=150); plt.close(fig)

# 3. integration variants (precursor level, all quantified precursors)
v = pd.read_csv(os.path.join(OUT, "diag_tables", "integration_variants.tsv"), sep="\t")
v = v[v.variant != "feature hill area (variable bounds)"]
fig, ax = plt.subplots(figsize=(7.2, 5.0))
y = range(len(v))
ax.scatter(v.prec_eps_H, y, color="#1baf7a", s=22, label="Human", zorder=3)
ax.scatter(v.prec_eps_Y, y, color="#eda100", s=22, label="Yeast", zorder=3)
ax.scatter(v.prec_eps_E, y, color="#4a3aa7", s=22, label="E. coli", zorder=3)
ax.scatter(v.prec_eps, y, color=INK, marker="|", s=120, label="all", zorder=4)
ax.set_yticks(list(y), v.variant)
ax.invert_yaxis()
ax.set_xlabel("median |epsilon| (precursor level, D9 6 runs, MBR off)")
ax.set_title("XIC re-integration variants (D9)", loc="left", color=INK)
ax.legend(frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0, 1.04))
fig.tight_layout(); fig.savefig(os.path.join(P, "integration_variants.png"), dpi=150); plt.close(fig)

# 4. paired residual vs apex deviation (feature-link errors)
s = pd.DataFrame(d["paired_residuals"]["strata"])
fig, axs = plt.subplots(1, 2, figsize=(8.4, 2.6))
for ax, f, lab in ((axs[0], "d_rt_cons_sd", "feature apex - consensus RT (alignment SDs)"), (axs[1], "n_scans", "MS1 scans with mono signal")):
    x = s[s.factor == f]
    ax.bar(range(len(x)), x["frac_gt_0.5"], color=OURS, width=0.6)
    ax.set_xticks(range(len(x)), x["bin"], rotation=30, fontsize=7)
    ax.set_xlabel(lab)
axs[0].set_ylabel("share of values with |paired residual| > 0.5")
fig.suptitle("Value-level disagreement with Sage (same peptide, same run)", x=0.01, ha="left", color=INK)
fig.tight_layout(); fig.savefig(os.path.join(P, "paired_residual_strata.png"), dpi=150); plt.close(fig)
print("ok")
