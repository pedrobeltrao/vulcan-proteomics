#!/usr/bin/env python3
"""A3.1.8(i): G4 (DESIGN §6.4, A1.4) for v2_sep with the partition extended to chimeric stratum {1, 2, >=3}.  Dev only.

PARTITION AND DECOY DEFINITION (frozen here and in results/2026-10-28_A3_preconditions/g4_partition.json, both committed
BEFORE this script is run for the first time; the run refuses if the partition file is not committed unmodified):
  * data: the 40 frozen v2 training runs (results/2026-10-28_E7_chimeric_v2/train_v2_folds.tsv minus the 3 runs E7's
    prep excluded label-independently), i.e. "training decoys" (§6.4: G4 is decided on training or dev decoys only).
    Their chimeric PSM tables and v2 features are E7's ($VULCAN_E7_WORK/train/<run>.chim.{psms,feat}.parquet, produced
    by scripts/e7/prep_train.py at the E7 commits); labels from each run's species prod index (scripts/e7/train_v2.INDEX).
  * score: the frozen v2 model (static_v2.pkl, sha256 3f2c69cc...) applied as bench/t1/ours_chimeric (StaticModel.apply).
  * decoys: STAGE-1 UNIT WINNERS that are decoys -- per query unit (run, spectrum, precursor group) the best PSM by the
    score with the production hash tie-break (pf_id, seed 11), as vulcan.chimeric.compete stage 1 and as A3.2.3 uses
    stage-1 winners for G1/G2/G5.  (scripts/m3e/g4.py used the PSM-level winner per spectrum for v1; the unit is the
    chimeric analogue.)
  * strata: instrument class OT-HCD x activation HCD (every training run; the class the scorer serves, as the A1 G4 of
    v1) x chimeric stratum of the unit {1, 2, >=3} (chim_rank pooled at 3, vulcan.stats.groupwalk.pool_chimeric_rank)
    x precursor charge of the PSM {2, 3, 4+}; charge 1 (if present) is kept as its own charge class (stricter, as
    m3e).  Label: "r{1|2|3}_z{1|2|3|4}".
  * statistic and tolerance: vulcan.gates.stats.g4_gate unchanged (tau = 95th percentile of within-stratum
    between-run KS over runs with >= 200 decoys in the stratum; pass iff every stratum's KS(stratum vs all other
    strata) <= tau).
  * verdict: G4(v2_sep) = g4_gate on the partition above for the v2 score.  Reported only (never decide): hyperscore
    reference; the A1 charge-only partition for v2; the pair-winner decoys after the full v2_sep competition (stratum =
    pair winner's); the partition split by instrument model instead of class.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/a3/dev_g4.py --out results/2026-10-28_A3_preconditions/g4.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "e7"))

import numpy as np  # noqa: E402

RES = os.path.join(ROOT, "results", "2026-10-28_A3_preconditions")
PARTITION = os.path.join(RES, "g4_partition.json")
V2 = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl")
V2_SHA = "3f2c69cc6cad002e42c51475dc1ad52b9bfbbc84f3dd51997b22cb24a70f0186"


def stratum_labels(chim_rank, charge):
    r = np.minimum(np.asarray(chim_rank, np.int64), 3)
    z = np.minimum(np.maximum(np.asarray(charge, np.int64), 1), 4)
    return np.char.add(np.char.add("r", r.astype(str)), np.char.add("_z", z.astype(str)))


def committed_unmodified(path):
    rel = os.path.relpath(path, ROOT)
    ls = subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", rel], capture_output=True)
    st = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return ls.returncode == 0 and not st.stdout.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if not committed_unmodified(PARTITION) or not committed_unmodified(os.path.abspath(__file__)):
        raise SystemExit("dev_g4: the partition file and this script must be committed (frozen) before estimation")
    import train_v2 as T
    from vulcan.chimeric.compete import chimeric_competition
    from vulcan.gates import stats as GS
    from vulcan.rescore.features_v2 import model_psm
    from vulcan.rescore.static import StaticModel
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import compete
    assert hashlib.sha256(open(V2, "rb").read()).hexdigest() == V2_SHA
    v2 = StaticModel.load(V2)
    runs, excluded = T.load()
    inst = {r["filename"][:-4]: r["instrument"] for r in
            csv.DictReader(open(os.path.join(T.RES, "train_v2_folds.tsv")), delimiter="\t")}
    acc = {k: {"score": [], "stratum": [], "charge_only": [], "run": [], "inst": []} for k in ("v2", "hyperscore")}
    pw = {"score": [], "stratum": [], "run": []}
    for r in runs:
        psm = r["psm"]
        s2 = v2.apply(model_psm(psm), r["feat"]).astype(np.float64)
        tb = hash_tiebreak(psm["pf_id"], 11)
        unit = psm["spec_idx"].astype(np.int64) * 64 + psm["group"].astype(np.int64)
        for name, s in (("v2", s2), ("hyperscore", psm["hyperscore"].astype(np.float64))):
            w = compete(unit, s, tiebreak=tb)
            w = w[psm["is_decoy"][w]]
            acc[name]["score"].append(s[w])
            acc[name]["stratum"].append(stratum_labels(psm["chim_rank"][w], psm["charge"][w]))
            acc[name]["charge_only"].append(np.minimum(psm["charge"][w].astype(int), 4))
            acc[name]["run"].append(np.full(w.size, r["run"]))
            acc[name]["inst"].append(np.full(w.size, inst[r["run"]]))
        cc = chimeric_competition(psm["spec_idx"].astype(np.int64), psm["group"], psm["chim_rank"], psm["pf_id"],
                                  psm["pair"], s2, psm["is_decoy"], psm["pf_row"], r["ions"], seed=11)
        res, win = cc["res"], cc["win"]
        d = res.is_decoy
        pw["score"].append(res.score[d])
        pw["stratum"].append(stratum_labels(psm["chim_rank"][win][d], psm["charge"][win][d]))
        pw["run"].append(np.full(int(d.sum()), r["run"]))
        print(r["run"], int(acc["v2"]["score"][-1].size), flush=True)
    cat = {k: {kk: np.concatenate(vv) for kk, vv in v.items()} for k, v in acc.items()}
    pwc = {k: np.concatenate(v) for k, v in pw.items()}
    part = json.load(open(PARTITION))
    out = {"partition": part, "n_runs": len(runs), "excluded_runs": excluded, "v2_sha256": V2_SHA,
           "decision": {}, "reported": {}}
    c = cat["v2"]
    out["decision"]["v2_sep_G4"] = GS.g4_gate(c["score"], c["stratum"], c["run"])
    out["reported"]["hyperscore_reference"] = GS.g4_gate(cat["hyperscore"]["score"], cat["hyperscore"]["stratum"],
                                                         cat["hyperscore"]["run"])
    out["reported"]["v2_A1_charge_only_partition"] = GS.g4_gate(c["score"], c["charge_only"], c["run"])
    out["reported"]["v2_pair_winner_decoys"] = GS.g4_gate(pwc["score"], pwc["stratum"], pwc["run"])
    out["reported"]["v2_by_instrument_model"] = GS.g4_gate(c["score"], np.char.add(np.char.add(c["inst"], "|"),
                                                                                    c["stratum"]), c["run"])
    out["pass"] = bool(out["decision"]["v2_sep_G4"]["pass"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    g = out["decision"]["v2_sep_G4"]
    print("G4 v2_sep", g["pass"], "tau", round(g["tau"], 4),
          {k: (v["n"], round(v["ks_vs_rest"], 4), v["pass_"]) for k, v in g["strata"].items()})
    for k, v in out["reported"].items():
        print(k, v["pass"], round(v["tau"], 4))


if __name__ == "__main__":
    main()
