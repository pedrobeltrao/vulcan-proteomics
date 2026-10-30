#!/usr/bin/env python3
"""3g POST HOC: permutation null for protein clustering of arm-C tail winners.  Under a per-pair tilt that is spread
evenly over pairs, E-labels are exchangeable across pairs given their total; we permute E_won across tail pairs
(1,000x, seed 11) and recompute D = (#E winners in base proteins holding >= 2 E winners) - (#dec winners in proteins
holding >= 2 dec winners).  Observed D far above the permutation distribution = excess concentrated in proteins
(coherent real-evidence matches: homologous or truly present plant material), not a uniform naturalness tilt.
Also: the coin when the 2 B-contaminated projects of step 3 (PXD069854, PXD024123: B true presence) are excluded,
and the coin outside plastid/photosynthesis + highly conserved proteins named in diag_a (descriptive only)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "entrapment"))
import numpy as np, pandas as pd  # noqa
import common as C  # noqa
from peptides import read_fasta, digest, il  # noqa
prot = {}
for h, s in read_fasta(os.path.join(C.ROOT, "data", "fasta", "arath_UP000006548.fasta")):
    g = h.split("GN=")[1].split()[0] if "GN=" in h else h.split()[0]
    for p in digest(s):
        prot.setdefault(il(p), set()).add(g)
rng = np.random.default_rng(C.SEED)
def D(w, pr):
    e = pd.Series(pr[w]).value_counts(); d = pd.Series(pr[~w]).value_counts()
    return int(e[e >= 2].sum() - d[d >= 2].sum())
out = {}
for col in ("score", "hyperscore"):
    t = pd.read_parquet(os.path.join(C.SCR, f"e9tail_C_{col}.parquet"))
    t["prot"] = t.base_seq.map(lambda b: sorted(prot.get(il(b), {"?"}))[0])
    for thr, m in (("0.05", np.ones(len(t), bool)), ("0.01", t.above1.values)):
        s = t[m]
        w, pr = s.E_won.values, s.prot.values
        obs = D(w, pr)
        perm = np.array([D(rng.permutation(w), pr) for _ in range(1000)])
        excess = int(w.sum() - (~w).sum())
        keep = ~s.project.isin(["PXD069854", "PXD024123"]).values
        out[f"{col}_{thr}"] = {"D_obs": obs, "D_perm_mean": float(perm.mean()), "D_perm_q975": float(np.quantile(perm, .975)),
                               "p_perm": float((perm >= obs).mean()), "total_E_minus_dec": excess,
                               "coin_excl_2_Bcontaminated_projects": C.coin(w[keep])}
        print(col, thr, out[f"{col}_{thr}"], flush=True)
json.dump(out, open(os.path.join(C.SCR, "diag_cluster_perm.json"), "w"), indent=1, default=float)
