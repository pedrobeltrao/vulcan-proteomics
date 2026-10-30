#!/usr/bin/env python3
"""A3.2.2 / A3.7: the frozen dev table config/a3_secondary_offsets.tsv (label-free; built before the tag).

Source: every SECONDARY precursor group of the D0 (2 runs) and D1 (1 run) vulcan.ms1 hypothesis tables, as the E7
chimeric search uses them (vulcan.chimeric.search.precursor_groups: feature groups that are not the instrument's
precursor; table defaults charge 2-6, MS1 chimeric_rank <= 5, isolation fraction >= 1%; the isotope-error-0 row of
each group).  Per group:
    delta_mz = secondary mono m/z - instrument precursor m/z of the MS2 spectrum (store precursor_mz), Th
    charge   = the secondary group's charge
No label, no search result and no DN2 file is read.  The MS1 tables are regenerated with the committed vulcan.ms1
(bit-identical to E7's, sha256 recorded in the .json sidecar).
Row order (= the empirical CDF used by the A-syn draw, vulcan.gates.a3_asyn): sorted by (delta_mz, charge, dataset,
run, spectrum index, group); delta_mz written with 6 decimals (%.6f) and re-read from the file, so the draw uses
exactly the committed values.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/a3/build_offsets.py --out config/a3_secondary_offsets.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

T1W = os.environ.get("VULCAN_T1_WORK", f"/tmp/vulcan_t1_{os.getuid()}")
A3W = os.environ.get("VULCAN_A3_WORK", "/tmp/vulcan_a3")
SOURCES = [("D0", "b1906_293T_proteinID_01A_QE3_122212"), ("D0", "b1937_293T_proteinID_01B_QE3_122212"),
           ("D1", "Velos005137")]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "config", "a3_secondary_offsets.tsv"))
    a = ap.parse_args()
    from vulcan.chimeric import search as CS
    from vulcan.gpu.engine import load_store
    rows, side = [], {"sources": []}
    for ds, run in SOURCES:
        store = os.path.join(T1W, ds, "stage", run + ".vspec.arrow")
        hyp = os.path.join(A3W, "ms1", ds, run + ".hypotheses.arrow")
        sp = load_store(store)
        g = CS.precursor_groups(CS.read_hypotheses(hyp), sp)
        sec = g["sec"]
        d = sec["mz"] - np.asarray(sp.precursor_mz, np.float64)[sec["spec"]]
        for i in range(len(d)):
            rows.append((round(float(d[i]), 6), int(sec["charge"][i]), ds, run, int(sec["spec"][i]), int(sec["group"][i])))
        side["sources"].append({"dataset": ds, "run": run, "hypotheses_sha256": sha(hyp), "store_sha256": sha(store),
                                "n_spectra": len(sp), "n_secondary_groups": int(len(d))})
    rows.sort()
    with open(a.out, "w") as fh:
        fh.write("delta_mz\tcharge\tdataset\trun\tspec_idx\tgroup\n")
        for r in rows:
            fh.write(f"{r[0]:.6f}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n")
    side["n_rows"] = len(rows)
    side["sha256"] = sha(a.out)
    dm = np.array([r[0] for r in rows])
    z = np.array([r[1] for r in rows])
    side["summary"] = {"delta_mz_quantiles_1_5_25_50_75_95_99": np.quantile(dm, [.01, .05, .25, .5, .75, .95, .99]).round(4).tolist(),
                       "charge_counts": {int(k): int(v) for k, v in zip(*np.unique(z, return_counts=True))}}
    json.dump(side, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(side, indent=1))


if __name__ == "__main__":
    main()
