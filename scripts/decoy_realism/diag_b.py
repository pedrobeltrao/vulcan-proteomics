#!/usr/bin/env python3
"""3g POST HOC diagnostic (written after seeing c_B ~ 0.92 at 1%): which B organisms / proteins / runs carry the B
target excess in the 3g search.  Changes no verdict.  Reports E-winner vs dec-winner counts by organism, the top
proteins and the project spread at the global 1%/5% cuts (frozen score)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "entrapment"))
import numpy as np, pandas as pd  # noqa
import common as C  # noqa
from peptides import read_fasta, digest, il, accession  # noqa
from build_3g_db import B_FASTAS  # noqa
W = "/tmp/vulcan_3g"
lab = C.labels(os.path.join(W, "index", "human_contam_3g"))
fam = dict(pd.read_csv(os.path.join(W, "family_3g.tsv"), sep="\t")[["seq", "family"]].values)
src, prot = {}, {}
for f in B_FASTAS:
    for h, s in read_fasta(os.path.join(C.ROOT, "data", "fasta", f)):
        g = h.split("GN=")[1].split()[0] if "GN=" in h else accession(h)
        for p in digest(s):
            q = il(p)
            src.setdefault(q, set()).add(f.split("_")[0] + ("" if not f.startswith("dn") else f[3:16]))
            prot.setdefault(q, set()).add(g)
man = C.manifest()[:100]
runs = [r["run"] for r in man]
proj = np.array([r["accession"] for r in man])
d = C.load_winners(W, "P", runs, "score")
g, cut = C.compete(d, lab)
rows = d["pf_row"][g.row]
f = lab["family"][rows]
out = {}
for al in (0.01, 0.05):
    m = (g.score >= cut[al]) & np.isin(f, (1, 3))
    base = lab["target_seq"][rows[m]]
    isE = f[m] == 1
    fb = np.array([fam[b] for b in base])
    sel = fb == "B"
    t = pd.DataFrame({"E": isE[sel], "base": base[sel], "proj": proj[d["run"][g.row[m]][sel]]})
    t["org"] = t.base.map(lambda b: "|".join(sorted(src.get(b, {"?"}))))
    t["prot"] = t.base.map(lambda b: "|".join(sorted(prot.get(b, {"?"})))[:40])
    r = {"by_org": t.groupby("org").E.agg(["sum", "count"]).sort_values("count", ascending=False).head(12).to_dict("index"),
         "by_project": t.groupby("proj").E.agg(["sum", "count"]).sort_values("count", ascending=False).head(12).to_dict("index"),
         "top_E_proteins": t[t.E].prot.value_counts().head(25).to_dict(),
         "n_proteins_E": int(t[t.E].prot.nunique()), "n_proteins_dec": int(t[~t.E].prot.nunique()),
         "n_E": int(t.E.sum()), "n_dec": int((~t.E).sum())}
    out[str(al)] = r
    print(al, json.dumps(r, indent=0, default=int)[:3000])
json.dump(out, open(os.path.join(C.SCR, "diag_b.json"), "w"), indent=1, default=int)
