#!/usr/bin/env python3
"""3g POST HOC diagnostic (after the step-3 B result showed true presence of foreign material): is the E9 arm-C
Arabidopsis excess spread evenly over projects / proteins (a per-competition tilt) or concentrated (true presence of
plant-like material)?  Label-symmetric descriptors: per project E vs dec(E) counts and a chi-square heterogeneity
test; leave-one-project-out coin range; protein clustering (number of base proteins with >= 2 E winners vs >= 2
dec(E) winners).  Changes no verdict."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "entrapment"))
import numpy as np, pandas as pd  # noqa
from scipy.stats import chi2_contingency  # noqa
import common as C  # noqa
from peptides import read_fasta, digest, il  # noqa
prot = {}
for h, s in read_fasta(os.path.join(C.ROOT, "data", "fasta", "arath_UP000006548.fasta")):
    g = h.split("GN=")[1].split()[0] if "GN=" in h else h.split()[0]
    for p in digest(s):
        prot.setdefault(il(p), set()).add(g)
out = {}
for col in ("score", "hyperscore"):
    t = pd.read_parquet(os.path.join(C.SCR, f"e9tail_C_{col}.parquet"))
    t["prot"] = t.base_seq.map(lambda b: sorted(prot.get(il(b), {"?"}))[0])
    r = {}
    for thr, m in (("0.05", np.ones(len(t), bool)), ("0.01", t.above1.values)):
        s = t[m]
        pj = s.groupby("project").E_won.agg(["sum", "count"])
        pj["exc"] = pj["sum"] - pj["count"] / 2
        ct = np.column_stack([pj["sum"], pj["count"] - pj["sum"]])
        ok = pj["count"] >= 10
        loo = [C.coin(s.E_won.values[s.project.values != p])["coin"] for p in pj.index]
        top = pj.sort_values("exc", ascending=False)
        cum = top.exc.cumsum() / pj.exc.sum()
        pe = s[s.E_won].prot.value_counts()
        pd_ = s[~s.E_won].prot.value_counts()
        r[thr] = {"n_projects": int(len(pj)), "heterogeneity_p": float(chi2_contingency(ct[ok.values])[1]),
                  "loo_coin_range": [float(min(loo)), float(max(loo))],
                  "n_projects_for_half_excess": int((cum < 0.5).sum() + 1),
                  "top_projects": top.head(8)[["sum", "count", "exc"]].to_dict("index"),
                  "coin_without_top5_excess_projects": C.coin(s.E_won.values[~s.project.isin(top.index[:5]).values]),
                  "proteins_ge2_E": int((pe >= 2).sum()), "proteins_ge2_dec": int((pd_ >= 2).sum()),
                  "winners_in_ge2_proteins_E": int(pe[pe >= 2].sum()), "winners_in_ge2_proteins_dec": int(pd_[pd_ >= 2].sum()),
                  "top_E_proteins": pe.head(12).to_dict(), "top_dec_proteins": pd_.head(8).to_dict()}
        print(col, thr, json.dumps(r[thr], default=float)[:1500], flush=True)
    out[col] = r
json.dump(out, open(os.path.join(C.SCR, "diag_a.json"), "w"), indent=1, default=float)
