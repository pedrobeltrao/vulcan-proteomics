#!/usr/bin/env python3
"""A2.9.3 tail proxy threshold: the hyperscore 1% threshold of each D0 host run (production search on the prod gate
index, paired double competition, hash tie-break seed 11; threshold = lowest accepted target hyperscore at q <= 1%,
vulcan.gates.stats.tail_threshold), and tail_proxy_threshold = the MAX over the two hosts (stricter: fewer A spectra
count towards the proxy, so the walk can only run longer).

This is the ONLY script under scripts/a2/ that reads gate-index labels before step 5, and only for the two human HOST
runs (D0), never for a null spectrum; it stores two scalars and no count (whitelisted in
tests/a2/test_no_label_counts.py).  Run once at step 3 before the first admit.

    python scripts/a2/host_threshold.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")


def main():  # pragma: no cover (GPU)
    import m3_common as C
    import m3_data as MD
    import numpy as np

    from vulcan.gates import stats as GS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    idx = C.index_dir("human_contam", "prod")
    lab = MD.labels(idx)
    l1 = L1Index(idx)
    p = E.SearchParams(report_psms=C.REPORT)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    thr = {}
    for raw in C.DEV["D0"]["raw"]:
        sp = E.load_store(C.E2.store_path(raw))
        tbl, _, _ = E.search(sp, gi, l1, p)
        spec = tbl.column("spec_idx").to_numpy()
        pf = tbl.column("pf_id").to_numpy()
        rows = tbl.column("pf_row").to_numpy()
        hs = tbl.column("hyperscore").to_numpy().astype(float)
        tb = hash_tiebreak(pf, 11)
        r = double_competition(spec.astype(np.int64), pf, lab["tdc_pair"][rows], hs, lab["is_decoy"][rows], seed=11,
                               tiebreak=tb)
        thr[sp.run] = GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01)
    os.makedirs(W, exist_ok=True)
    out = {"host_hyperscore_1pct": thr, "tail_proxy_threshold": max(thr.values()), "rule": "max over D0 hosts"}
    json.dump(out, open(os.path.join(W, "host_threshold.json"), "w"), indent=1)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
