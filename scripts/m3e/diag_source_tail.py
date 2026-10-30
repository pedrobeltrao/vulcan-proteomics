#!/usr/bin/env python3
"""Diagnostic (not a gate; post hoc): which human proteins do the tail winners of each null source map to?
Recomputes, per mixture, the frozen scorer's and hyperscore's PSM winners and own 1% thresholds (as gates_a1.py), keeps
A1.3-clean species-null tail winners (target or decoy) and reports, per source, the genes of target winners
(human_sp_canonical GN=) -- e.g. plasma / serum / keratin proteins would indicate real human material in that sample.
    python scripts/m3e/diag_source_tail.py"""
import collections, json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3")); sys.path.insert(0, HERE)
import m3_common as C, m3_data as MD  # noqa: E402
import numpy as np  # noqa: E402
import evaluate as EV  # noqa: E402
from gates_a1 import _load, homolog_filters, PROV_NAME  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.index.build import L1Index  # noqa: E402
from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition, psm_competition  # noqa: E402
W = "/tmp/vulcan_m3e"
idx = C.index_dir("human_contam", "prod"); lab = MD.labels(idx); l1 = L1Index(idx)
model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
hfs = homolog_filters()
gn = {}
for line in open(os.path.join(C.ROOT, "data", "fasta", "human_sp_canonical.fasta")):
    if line.startswith(">"):
        acc = line.split("|")[1]
        m = re.search(r"GN=(\S+)", line)
        gn[acc] = m.group(1) if m else acc
pend = []
names = sorted(os.path.basename(f)[:-len(".feat.parquet")] for f in __import__("glob").glob(os.path.join(W, "psms", "G2mix", "*.feat.parquet")))
res = {p: collections.defaultdict(collections.Counter) for p in ("hyperscore", "gbm")}
cnt = {p: collections.defaultdict(lambda: [0, 0]) for p in ("hyperscore", "gbm")}
for name in names:
    run, psm, feat = _load(os.path.join(W, "psms", "G2mix"), [name])[0]
    prov = np.load(os.path.join(W, "null_mix", name + ".prov.npy"))
    rows, spec = psm["pf_row"], psm["spec_idx"]; dec = lab["is_decoy"][rows]; tb = hash_tiebreak(psm["pf_id"], 11)
    for p, sc in (("hyperscore", psm["hyperscore"].astype(float)), ("gbm", model.apply(psm, feat))):
        w = psm_competition(spec, sc, tiebreak=tb)
        r = double_competition(spec.astype(np.int64), psm["pf_id"], lab["tdc_pair"][rows], sc, dec, seed=11, tiebreak=tb)
        thr = GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01)
        sel = w[(prov[spec[w]] >= 2) & (sc[w] >= thr)]
        for i in sel:
            row = int(rows[i])
            pend.append((p, int(prov[spec[i]]), row, int(l1.partner[row]), bool(dec[i]), str(lab["target_seq"][row])))
    print(name, flush=True)
near_rows = {}
for src, hf in hfs.items():
    rr = sorted({x for q in pend if q[1] == src for x in q[2:4] if x >= 0})
    near_rows[src] = set(np.array(rr)[hf.near([l1.sequence(r) for r in rr])].tolist()) if rr else set()
for p, src, row, part, isdec, ts in pend:
    if row in near_rows[src] or part in near_rows[src]:
        continue
    cnt[p][PROV_NAME[src]][int(isdec)] += 1
    if not isdec:
        for a in EV.protein_maps(idx, {ts})[ts][:1]:
            res[p][PROV_NAME[src]][gn.get(a.split("|")[1] if "|" in a else a, a)] += 1
out = {p: {s: {"targets_decoys": cnt[p][s], "top_genes": res[p][s].most_common(25)} for s in cnt[p]} for p in cnt}
json.dump(out, open(os.path.join(C.ROOT, "results", "2026-10-27_M3e_hardening", "diag_source_tail.json"), "w"), indent=1)
for p in out:
    for s, v in out[p].items():
        print(p, s, v["targets_decoys"], v["top_genes"][:12])
