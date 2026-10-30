#!/usr/bin/env python3
"""Reference: Sage 0.14.7's own LDA discriminant (E2 Sage runs on the SAME index/spectra, /tmp/vulcan_e2/sage_runs/
<ds>_prod_oracle, top-1 PSM) put through OUR estimator (paired double competition, peptidoform units, hash tie-break
11), next to Sage's own peptide_q.  No new Sage run (the baselines agent owns Sage+LDA / Sage+mokapot with entrapment).

    python scripts/m3/sage_lda_ours.py   -> results/2026-10-27_M3_scoring/sage_lda_ours.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.index.build import L1Index  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402


def main():
    out = {}
    for ds, fam in (("D0", "human_contam"), ("D1", "pfu_contam")):
        d = os.path.join(C.E2.WORK, "sage_runs", f"{ds}_prod_oracle")
        idx = C.index_dir(fam, "prod")
        l1 = L1Index(idx)
        lab = MD.labels(idx)
        s = SB.read_results(d).to_pandas()
        s = s[s["rank"] == 1]
        pm = [SB.parse_peptide(x, C.VAR) for x in s.peptide]
        need = set(pm)
        key = {}
        for i in range(len(l1)):
            k = (l1.sequence(i), int(l1.mask[i]))
            if k in need:
                key[k] = i
        rows = np.array([key[k] for k in pm])
        dec = lab["is_decoy"][rows]
        assert np.array_equal(dec, s.label.values == -1)
        spec = (s.filename + ":" + s.scannr.astype(str)).values
        unit = l1.pf_id[rows]
        r = double_competition(spec, unit, lab["tdc_pair"][rows], s.sage_discriminant_score.values, dec, seed=11,
                               tiebreak=hash_tiebreak(unit, 11))
        tseq = lab["target_seq"][rows]
        res = {"n_psms": int(len(rows))}
        for a in (0.001, 0.01, 0.05):
            acc = (r.q <= a) & ~r.is_decoy
            tg = s[(s.label == 1) & (s.peptide_q <= a)]
            res[str(a)] = {"sage_lda_vulcan_paired_tdc_peptidoforms": int(acc.sum()),
                           "sage_lda_vulcan_unique_sequences": int(np.unique(tseq[r.row[acc]]).size),
                           "sage_own_peptide_q_peptidoforms": int(tg.peptide.nunique())}
        out[ds] = res
        print(ds, json.dumps(res))
    json.dump(out, open(os.path.join(C.RES, "sage_lda_ours.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
