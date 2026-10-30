#!/usr/bin/env python3
"""Follow-up to g2_null_tail.py: break the null winners above the q<=1% threshold down by provenance and by
near-homology to the foreign proteome (Hamming distance <= 2 after I->L, or a shared 6-mer), to tell an imperfect
null (near-true matches that a better scorer legitimately ranks up) from a label artifact.  Dev only."""
import os, sys, json
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import m3_common as C, m3_data as MD
import gates as GT
from vulcan.rescore.static import StaticModel
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import double_competition
from scipy.stats import binomtest

STATIC = ["static_linear_base", "static_gbm_base", "static_linear_pred", "static_gbm_pred"]
ms = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in STATIC}
lab = MD.labels(C.index_dir("human_contam", "prod"))

def kmers(F, k=6):
    S = set()
    for f in F:
        f = f.replace("I", "L")
        for i in range(len(f) - k + 1):
            S.add(f[i:i + k])
    return S
def ham2keys(F):
    S = set()
    for f in F:
        f = f.replace("I", "L"); n = len(f)
        for i in range(n):
            for j in range(i + 1, n):
                S.add((f[:i], f[i + 1:j], f[j + 1:], n))
    return S
FK, FH = {}, {}
def near(seq, src):
    if src not in FK:
        F = GT.foreign_sets(src)[0]
        FK[src] = kmers(F); FH[src] = ham2keys(F)
    s = seq.replace("I", "L"); n = len(s)
    if any(s[i:i + 6] in FK[src] for i in range(n - 5)): return "shared6mer"
    if any((s[:i], s[i + 1:j], s[j + 1:], n) in FH[src] for i in range(n) for j in range(i + 1, n)): return "ham2"
    return "far"

res = {}
rows_out = []
for run, psm, feat in MD.load_runs(C.psm_dir("G2_mix")):
    prov = np.load(os.path.join(C.WORK, "store_null", run + ".provenance.npy"))
    rows, spec = psm["pf_row"], psm["spec_idx"]
    dec = lab["is_decoy"][rows]; tseq = lab["target_seq"][rows]; pv = prov[spec]
    bad = np.zeros(prov.size, bool)
    for src in (2, 3):
        m = pv == src
        b = GT.homolog_excluded_seqs(set(tseq[m].tolist()), src)
        bad[spec[m & np.isin(tseq, list(b))]] = True
    tb = hash_tiebreak(psm["pf_id"], 11)
    key = np.array([f"{run}:{x}" for x in spec])
    pair = lab["tdc_pair"][rows]
    scores = {"hyperscore": psm["hyperscore"].astype(float)}
    scores.update({s: m.apply(psm, feat) for s, m in ms.items()})
    for s, sc in scores.items():
        r = double_competition(key, psm["pf_id"], pair, sc, dec, seed=11, tiebreak=tb)
        ok = r.q <= 0.01
        tmin = r.score[ok & ~r.is_decoy].min()
        above = (r.score >= tmin) & (pv[r.row] > 0) & ~bad[spec[r.row]]
        for i in np.flatnonzero(above):
            row = r.row[i]; p = int(pv[row])
            # for decoys, homology is assessed on the pair's TARGET sequence (label-symmetric)
            nh = near(tseq[row], p) if p in (2, 3) else "offset"
            d = res.setdefault(s, {}).setdefault(f"prov{p}:{nh}", {"T": 0, "D": 0})
            d["D" if r.is_decoy[i] else "T"] += 1
for s in res:
    tot = {"T": 0, "D": 0}
    for k in sorted(res[s]):
        tot["T"] += res[s][k]["T"]; tot["D"] += res[s][k]["D"]
    far = {"T": sum(v["T"] for k, v in res[s].items() if k.endswith("far")), "D": sum(v["D"] for k, v in res[s].items() if k.endswith("far"))}
    n = far["T"] + far["D"]
    res[s]["_far_only"] = dict(far, binom_p=float(binomtest(far["T"], n, 0.5).pvalue) if n else None)
    res[s]["_total"] = tot
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(os.path.dirname(__file__), "g2_null_tail_breakdown.json"), "w"), indent=1)
