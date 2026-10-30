#!/usr/bin/env python3
"""G4 (DESIGN §6.4, A1.4): calibration across strata, decided on TRAINING decoys only.
Partition (frozen before estimation): instrument class OT-HCD x activation HCD x precursor charge {2, 3, 4+}
(the only class/activation the frozen scorer serves).  Decoys = spectra whose PSM-level winner (top-5 candidates x
hypotheses, hash tie-break seed 11) under the scorer is a decoy, on the 12 OT-HCD training runs (train.tsv; the IT-CID
run MN-WT1-Frac8 excluded as in training).  Statistic and tolerance: vulcan.gates.stats.g4_gate (tau = 95th pct of
within-stratum between-run KS; pass iff every stratum's KS vs the rest <= tau).  Reported for static_gbm_pred, its
RT-neutralised ablation, and hyperscore (reference).
    conda_envs/pred/bin/python scripts/m3e/g4.py --out results/.../g4.json"""
import argparse, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3")); sys.path.insert(0, HERE)
import m3_common as C, m3_data as MD  # noqa: E402
import numpy as np  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import psm_competition  # noqa: E402
from gates_a1 import neutral_rt  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    acc = {p: {"score": [], "stratum": [], "run": []} for p in ("hyperscore", "gbm", "gbm_noRT")}
    for sp in ("human", "mouse", "ecoli"):
        lab = MD.labels(C.index_dir(f"{sp}_contam", "prod"))
        for run, psm, feat in MD.load_runs(os.path.join(C.WORK, "psms", f"train_{sp}")):
            if run == "MN-WT1-Frac8":
                continue
            dec = lab["is_decoy"][psm["pf_row"]]
            tb = hash_tiebreak(psm["pf_id"], 11)
            scores = {"hyperscore": psm["hyperscore"].astype(float), "gbm": model.apply(psm, feat),
                      "gbm_noRT": model.apply(psm, neutral_rt(feat))}
            for p, s in scores.items():
                w = psm_competition(psm["spec_idx"], s, tiebreak=tb)
                w = w[dec[w]]
                z = psm["charge"][w].astype(int)
                acc[p]["score"].append(s[w]); acc[p]["stratum"].append(np.minimum(z, 4)); acc[p]["run"].append(np.full(w.size, run))
            print(sp, run, int(w.size), flush=True)
    out = {"partition": "OT-HCD x HCD x charge {2,3,4+} (pre-specified) + charge 1 found in the data, kept as its own stratum (stricter)", "runs": "12 OT-HCD training runs (train.tsv minus TR10 IT-CID)"}
    for p, d in acc.items():
        out[p] = GS.g4_gate(np.concatenate(d["score"]), np.concatenate(d["stratum"]), np.concatenate(d["run"]))
        print(p, out[p]["pass"], round(out[p]["tau"], 4), {k: round(v["ks_vs_rest"], 4) for k, v in out[p]["strata"].items()})
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
