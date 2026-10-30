#!/usr/bin/env python3
"""M3 / E4 evaluation on dev data (D0, D1): IDs at q <= {0.1, 1, 5}% (peptidoform units, paired double competition)
and picked protein groups, with §6.3 dev entrapment verdicts (arm P: seeds 2001-2003, pseudo-reverse decoys; arm C:
natural foreign proteome) for every scorer.

Scorers
  hyperscore              fixed hand-crafted score (required arms: P or C)
  static_{linear,gbm}_{base,pred}   FROZEN static-first models trained on project-disjoint data (required: P and C);
                          plain paired TDC + picked protein (DESIGN §6.4: not fitted on evaluated data -> no RESET)
  reset_lda_{base,pred}   RESET (vulcan.stats.reset, single-decoy, LDA, seed 11) on the hyperscore double-competition
                          winners; pooled peptide-level FDR ONLY (§6.4); fallback rate reported
    conda_envs/pred/bin/python scripts/m3/evaluate.py [--datasets D0 D1] [--scorers ...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402

from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, StaticModel, feature_matrix  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.protein import picked_group_fdr  # noqa: E402
from vulcan.stats.reset import reset  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402

ALPHAS = (0.001, 0.01, 0.05)
SEED = 11
FAM = {"D0": "human_contam", "D1": "pfu_contam"}
STATIC = ["static_linear_base", "static_gbm_base", "static_linear_pred", "static_gbm_pred"]
RESETS = {"reset_lda_base": BASE_FEATURES, "reset_lda_pred": BASE_FEATURES + PRED_FEATURES}


class Search:
    """One search (dataset x index kind): PSMs of all runs (all ranks), features, labels (L4)."""

    def __init__(self, ds, kind):
        self.ds, self.kind = ds, kind
        self.index = C.index_dir(FAM[ds], kind)
        self.runs = MD.load_runs(C.psm_dir(f"{ds}_{kind}"))
        self.psm = MD.concat([p for _, p, _ in self.runs])
        self.feat = MD.concat([f for _, _, f in self.runs])
        self.lab = MD.labels(self.index)
        r = self.psm["pf_row"]
        self.is_decoy = self.lab["is_decoy"][r]
        self.is_ent = self.lab["is_entrapment"][r]
        self.pair = self.lab["tdc_pair"][r]
        self.unit = self.psm["pf_id"]
        self.spec = np.array([f"{a}:{b}" for a, b in zip(self.psm["run"], self.psm["spec_idx"])])
        self.tb = hash_tiebreak(self.unit, SEED)

    def static_score(self, model):
        return np.concatenate([model.apply(p, f) for _, p, f in self.runs])

    def matrix(self, names):
        return np.vstack([feature_matrix(p, f, names) for _, p, f in self.runs])

    def compete(self, score):
        return double_competition(self.spec, self.unit, self.pair, score, self.is_decoy, seed=SEED, tiebreak=self.tb)


def scored(search: Search, scorer, models):
    """-> dict(row, score, is_decoy, q, extra) over pair-competition winners."""
    if scorer == "hyperscore" or scorer in STATIC:
        s = search.psm["hyperscore"].astype(float) if scorer == "hyperscore" else search.static_score(models[scorer])
        r = search.compete(s)
        return {"row": r.row, "score": r.score, "is_decoy": r.is_decoy, "q": r.q, "extra": {}}
    r = search.compete(search.psm["hyperscore"].astype(float))
    X = search.matrix(RESETS[scorer])[r.row]
    res = reset(X, r.is_decoy, init_score=r.score, mode="single", seed=SEED, alpha_train=0.01)
    fallback, retry = False, False
    if not res.n_positives or res.n_positives[0] < 10:
        retry = True
        res = reset(X, r.is_decoy, init_score=r.score, mode="single", seed=SEED, alpha_train=0.05)
        fallback = not res.n_positives or res.n_positives[0] < 10
    keep = ~res.is_training_decoy
    return {"row": r.row[keep], "score": res.score[keep], "is_decoy": r.is_decoy[keep], "q": res.q[keep],
            "extra": {"n_positives_per_iter": res.n_positives, "retry_alpha_train_5pct": retry,
                      "fallback_fixed_score": fallback, "n_training_decoys": int(res.is_training_decoy.sum())}}


# ------------------------------------------------------------------------------------------- proteins
_prot_cache = {}


def protein_maps(index_dir, seqs_needed):
    key = index_dir
    if key not in _prot_cache:
        with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
            t = rd.read_all()
        pep = t.column("peptide").to_pylist()
        off = t.column("prot_off").to_numpy()
        ln = t.column("prot_len").to_numpy()
        pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
        acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
        _prot_cache[key] = (dict(zip(pep, range(len(pep)))), off, ln, pidx, acc)
    pos, off, ln, pidx, acc = _prot_cache[key]
    out = {}
    for s in seqs_needed:
        i = pos.get(s)
        out[s] = [] if i is None else [acc[j] for j in pidx[off[i]:off[i] + ln[i]]]
    return out


def protein_level(search: Search, w, alphas, entrapment=False):
    """Picked protein groups from the pair-competition winners (best peptide score per group)."""
    rows = search.psm["pf_row"][w["row"]]
    fam = search.lab["family"][rows]
    tseq = search.lab["target_seq"][rows]
    eseed = search.lab["entrap_seed"][rows]
    keep = eseed >= 0          # arm-C foreign peptides have no protein (protein entrapment is arm P only, §6.3)
    m = protein_maps(search.index, set(tseq[keep].tolist()))
    pre = {0: "", 1: "E_", 2: "rev_", 3: "rev_E_"}
    pep2prot, pep_score = {}, {}
    for i in np.flatnonzero(keep):
        prots = m[tseq[i]]
        if not prots:
            continue
        u = int(search.unit[w["row"][i]])
        pep2prot[u] = [pre[int(fam[i])] + p for p in prots]
        pep_score[u] = max(pep_score.get(u, -np.inf), float(w["score"][i]))
    allp = {p for v in pep2prot.values() for p in v}
    base = lambda p: p[4:] if p.startswith("rev_") else p  # noqa: E731
    pair_of = {p: base(p) for p in allp}
    is_dec = {p: p.startswith("rev_") for p in allp}
    res = picked_group_fdr(pep2prot, pep_score, pair_of, is_dec, seed=SEED)
    out = {}
    for a in alphas:
        acc = (res.q <= a) & ~res.is_decoy & ~res.is_obsolete
        if not entrapment:
            out[a] = {"n_groups": int(acc.sum())}
        else:
            is_e = {p: p.startswith("E_") for p in allp}
            partner = {p: p[2:] for p in allp if p.startswith("E_")}
            tside = ~res.is_decoy & ~res.is_obsolete
            groups = [g for g, t in zip(res.groups, tside) if t]
            out[a] = ENT.protein_paired_counts(groups, res.score[tside], acc[tside], is_e, partner)
    return out


# ------------------------------------------------------------------------------------------- entrapment
def peptide_entrapment(search: Search, w, alpha, estimator):
    tside = ~w["is_decoy"]
    rows = w["row"][tside]
    acc = w["q"][tside] <= alpha
    ent = search.is_ent[rows]
    if estimator == "paired":
        return ENT.paired_counts(search.lab["entrap_pair"][search.psm["pf_row"][rows]], ent, w["score"][tside], acc)
    return ENT.combined_counts(search.unit[rows], ent, acc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1"])
    ap.add_argument("--scorers", nargs="+", default=["hyperscore"] + STATIC + list(RESETS))
    ap.add_argument("--models", default=os.path.join(C.RES, "models"))
    ap.add_argument("--out", default=os.path.join(C.RES, "eval.json"))
    a = ap.parse_args()
    models = {s: StaticModel.load(os.path.join(a.models, s + ".pkl")) for s in a.scorers if s in STATIC}
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for ds in a.datasets:
        t0 = time.perf_counter()
        searches = {k: Search(ds, k) for k in ("prod", "E2001", "E2002", "E2003", "C")}
        r_c = C.c_arm_r(FAM[ds])
        print(ds, "loaded", round(time.perf_counter() - t0, 1), "s", flush=True)
        res_ds = out.setdefault(ds, {})
        for sc in a.scorers:
            t1 = time.perf_counter()
            mclass = "fixed" if sc == "hyperscore" else "static_low_capacity"
            is_reset = sc in RESETS
            ws = {k: scored(s, sc, models) for k, s in searches.items()}
            w = ws["prod"]
            rec = {"model_class": mclass, "ids": {}, "proteins": {}, "entrapment": {}, "protein_entrapment": {}}
            rec["ids_unique_sequences"] = {}
            tseq_w = searches["prod"].lab["target_seq"][searches["prod"].psm["pf_row"][w["row"]]]
            for al in ALPHAS:
                acc = (w["q"] <= al) & ~w["is_decoy"]
                rec["ids"][str(al)] = int(acc.sum())
                rec["ids_unique_sequences"][str(al)] = int(np.unique(tseq_w[acc]).size)
            if not is_reset:
                pl = protein_level(searches["prod"], w, ALPHAS)
                rec["proteins"] = {str(al): v["n_groups"] for al, v in pl.items()}
                pe = {k: protein_level(searches[k], ws[k], ALPHAS, entrapment=True) for k in ("E2001", "E2002", "E2003")}
            else:
                rec["proteins"] = "not reported: RESET mode gives pooled peptide-level FDR only (DESIGN §6.4)"
                rec["reset"] = {k: v["extra"] for k, v in ws.items()}
            for al in ALPHAS:
                reps = [ENT.Replicate(peptide_entrapment(searches[k], ws[k], al, "paired"), int(k[1:]), 0)
                        for k in ("E2001", "E2002", "E2003")]
                armP = ENT.evaluate_arm(reps, al, "paired")
                armC = ENT.evaluate_arm([ENT.Replicate(peptide_entrapment(searches["C"], ws["C"], al, "combined"), -1, 0)],
                                        al, "combined", r=r_c)
                vs = {"P": armP["verdict"], "C": armC["verdict"]}
                # §6.3 weaker claim when inconclusive: "consistent with control" if the lower CI of paired FDP <= alpha
                armP["consistent_with_control"] = bool(armP["ci95"][0] <= al)
                rec["entrapment"][str(al)] = {"P": armP, "C": armC, "overall": ENT.combine_arms(vs, mclass)}
                if not is_reset:
                    preps = [ENT.Replicate(pe[k][al], int(k[1:]), 0) for k in ("E2001", "E2002", "E2003")]
                    rec["protein_entrapment"][str(al)] = ENT.evaluate_arm(preps, al, "paired")
            rec["eval_wall_s"] = time.perf_counter() - t1
            res_ds[sc] = rec
            e1 = rec["entrapment"]["0.01"]
            print(ds, sc, rec["ids"], rec["proteins"] if not is_reset else "-", "P", round(e1["P"]["mean_fdp"], 4),
                  e1["P"]["verdict"], "C", round(e1["C"]["mean_fdp"], 4), e1["C"]["verdict"], "->", e1["overall"], flush=True)
            json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
