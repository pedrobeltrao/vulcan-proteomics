#!/usr/bin/env python3
"""A1.1 confirmation on the D8 holdout: the first 10 human runs of the E9 manifest (data/manifests/e9_d8.tsv, manifest
order; fetched by scripts/m3e/fetch_d8.py, sha256 == E9 record).  Not used for any selection.  Dev role.
Pipeline exactly as M3/D0 (scripts/m3): production search (prod index, report_psms 5) + arm P (E2001-E2003) + arm C
(Arabidopsis) indexes, L3 features (scripts/m3/features.py), frozen static_gbm_pred vs hyperscore.
Each run is its own experiment: paired double competition per run (peptidoform units), IDs at q <= 0.1/1/5%,
picked protein groups; entrapment counts per run; pooled verdict = union of the per-run lists (keys prefixed by run).
    --search (GPU) | --features (pred env) | --eval
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
KINDS = ("prod", "E2001", "E2002", "E2003", "C")
ALPHAS = (0.001, 0.01, 0.05)
RES = os.path.join(C.ROOT, "results", "2026-10-27_M3e_hardening")


def runs():
    return sorted(os.path.basename(p)[:-len(".vspec.arrow")] for p in glob.glob(os.path.join(W, "store_d8", "*.vspec.arrow")))


def pdir(kind):
    return os.path.join(W, "psms", f"D8_{kind}")


def search():
    import pyarrow.parquet as pq
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=C.REPORT)
    os.makedirs(os.path.join(W, "spectra_d8"), exist_ok=True)
    stores = []
    for r in runs():
        sp = E.load_store(os.path.join(W, "store_d8", r + ".vspec.arrow"))
        pk = os.path.join(W, "spectra_d8", r + ".spectra.pkl")
        if not os.path.exists(pk):
            pickle.dump(sp, open(pk, "wb"), protocol=5)
        stores.append(sp)
    for kind in KINDS:
        os.makedirs(pdir(kind), exist_ok=True)
        l1 = L1Index(C.index_dir("human_contam", kind))
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        for sp in stores:
            out = os.path.join(pdir(kind), sp.run + ".psms.parquet")
            if not os.path.exists(out):
                tbl, st, _ = E.search(sp, gi, l1, p)
                pq.write_table(tbl, out)
                print(kind, sp.run, len(sp), tbl.num_rows, flush=True)
        del gi
        torch.cuda.empty_cache()


def features():
    py = os.path.join(C.ROOT, "conda_envs", "pred", "bin", "python")
    subprocess.run([py, os.path.join(C.ROOT, "scripts", "m3", "features.py"), "--psms", *[pdir(k) for k in KINDS],
                    "--index", *[C.index_dir("human_contam", k) for k in KINDS], "--dataset", "D8",
                    "--spectra", os.path.join(W, "spectra_d8")], check=True)


class RunSearch:
    """evaluate.Search equivalent for one D8 run and one index kind."""

    def __init__(self, run, kind, lab):
        import m3_data as MD
        from vulcan.stats._util import hash_tiebreak
        self.ds, self.kind = "D8", kind
        self.index = C.index_dir("human_contam", kind)
        self.runs = [x for x in MD.load_runs(pdir(kind)) if x[0] == run]
        self.psm = self.runs[0][1]
        self.feat = self.runs[0][2]
        self.lab = lab
        r = self.psm["pf_row"]
        self.is_decoy = lab["is_decoy"][r]
        self.is_ent = lab["is_entrapment"][r]
        self.pair = lab["tdc_pair"][r]
        self.unit = self.psm["pf_id"]
        self.spec = self.psm["spec_idx"].astype(np.int64)
        self.tb = hash_tiebreak(self.unit, 11)

    def static_score(self, model):
        return np.concatenate([model.apply(p, f) for _, p, f in self.runs])

    def compete(self, score):
        from vulcan.stats.tdc import double_competition
        return double_competition(self.spec, self.unit, self.pair, score, self.is_decoy, seed=11, tiebreak=self.tb)


def evaluate():
    import evaluate as EV
    import m3_data as MD
    from vulcan.rescore.static import StaticModel
    from vulcan.stats import entrapment as ENT
    models = {"static_gbm_pred": StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))}
    labs = {k: MD.labels(C.index_dir("human_contam", k)) for k in KINDS}
    r_c = C.c_arm_r("human_contam")
    out = {"runs": runs(), "per_run": {}, "pooled": {}, "r_C": r_c}
    pooled = {sc: {a: {"P": {s: [] for s in (2001, 2002, 2003)}, "C": []} for a in ALPHAS} for sc in ("hyperscore", "static_gbm_pred")}
    for run in runs():
        searches = {k: RunSearch(run, k, labs[k]) for k in KINDS}
        rr = out["per_run"].setdefault(run, {})
        for sc in ("hyperscore", "static_gbm_pred"):
            ws = {k: EV.scored(s, sc, models) for k, s in searches.items()}
            w = ws["prod"]
            tseq = searches["prod"].lab["target_seq"][searches["prod"].psm["pf_row"][w["row"]]]
            rec = {"ids": {}, "ids_unique_sequences": {}, "proteins": {}, "entrapment": {}}
            for a in ALPHAS:
                acc = (w["q"] <= a) & ~w["is_decoy"]
                rec["ids"][str(a)] = int(acc.sum())
                rec["ids_unique_sequences"][str(a)] = int(np.unique(tseq[acc]).size)
            pl = EV.protein_level(searches["prod"], w, ALPHAS)
            rec["proteins"] = {str(a): v["n_groups"] for a, v in pl.items()}
            for a in ALPHAS:
                pcs = {}
                for k in ("E2001", "E2002", "E2003"):
                    pc = EV.peptide_entrapment(searches[k], ws[k], a, "paired")
                    pcs[k] = pc
                    pc2 = ENT.EntrapmentCounts(key=np.array([f"{run}|{x}" for x in pc.key]), t=pc.t, e=pc.e,
                                               e_only=pc.e_only, e_gt_t=pc.e_gt_t)
                    pooled[sc][a]["P"][int(k[1:])].append(pc2)
                cc = EV.peptide_entrapment(searches["C"], ws["C"], a, "combined")
                pooled[sc][a]["C"].append(ENT.EntrapmentCounts(key=np.array([f"{run}|{x}" for x in cc.key]), t=cc.t,
                                                               e=cc.e, e_only=cc.e_only, e_gt_t=cc.e_gt_t))
                armP = {cv: ENT.evaluate_arm([ENT.Replicate(pcs[k], int(k[1:]), 1001) for k in pcs], a, "paired",
                                             n_boot=1000, ci_convention=cv) for cv in ("a", "b")}
                armC = ENT.evaluate_arm([ENT.Replicate(cc, -1, 1001)], a, "combined", r=r_c, n_boot=1000)
                rec["entrapment"][str(a)] = {"P": {cv: {k: armP[cv][k] for k in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}
                                                   for cv in armP},
                                             "C": {k: armC[k] for k in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}}
            rr[sc] = rec
            print(run, sc, rec["ids"], rec["proteins"], {a: (round(v["P"]["a"]["mean_fdp"], 4), v["P"]["a"]["verdict"],
                                                            round(v["C"]["mean_fdp"], 4), v["C"]["verdict"])
                                                        for a, v in rec["entrapment"].items()}, flush=True)
    cat = lambda L: ENT.EntrapmentCounts(key=np.concatenate([x.key for x in L]), t=np.concatenate([x.t for x in L]),  # noqa: E731
                                         e=np.concatenate([x.e for x in L]), e_only=np.concatenate([x.e_only for x in L]),
                                         e_gt_t=np.concatenate([x.e_gt_t for x in L]))
    for sc, d in pooled.items():
        for a, arms in d.items():
            reps = [ENT.Replicate(cat(v), s, 1001) for s, v in arms["P"].items()]
            P = {cv: ENT.evaluate_arm(reps, a, "paired", ci_convention=cv) for cv in ("a", "b")}
            Cc = ENT.evaluate_arm([ENT.Replicate(cat(arms["C"]), -1, 1001)], a, "combined", r=r_c)
            out["pooled"].setdefault(sc, {})[str(a)] = {"P": P, "C": Cc}
            print("POOLED", sc, a, "P", round(P["a"]["mean_fdp"], 5), P["a"]["verdict"], P["b"]["verdict"],
                  "C", round(Cc["mean_fdp"], 5), Cc["verdict"], flush=True)
    json.dump(out, open(os.path.join(RES, "d8_confirm.json"), "w"), indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--features", action="store_true")
    ap.add_argument("--eval", action="store_true")
    a = ap.parse_args()
    if a.search:
        search()
    if a.features:
        features()
    if a.eval:
        evaluate()


if __name__ == "__main__":
    main()
