#!/usr/bin/env python3
"""E9 analysis (DESIGN §7 E9): global FDP-hat vs N, naive merge vs the L6 global procedure (vulcan.stats.global_fdr).

    conda_envs/pe/bin/python scripts/e9/analyze.py --ns 1 10 100 [1000] [--work /tmp/vulcan_e9] --out results/.../e9.json

Inputs: work/out/{P,C}/<run>.scored.parquet (frozen static_gbm_pred score; scripts/e9/score_worker.py) for the first N
runs of data/manifests/e9_d8.tsv (nested prefixes).  Labels are read here (L4 only).

Per N and arm:
  naive   per-run paired double competition (run q), accepted units at q_run <= alpha, UNION over runs
  global  best evidence per unit across runs -> global paired competition -> global q (L6)
Peptidoform level: arm P (paired, E_2001, r = 1) and arm C (combined, Arabidopsis, r from the M1e C database).
Protein level (arm P only, §6.3): picked protein groups (vulcan.stats.protein.picked_group_fdr, best peptide score,
  rescued subsets) on (a) every run's own pair winners -> run lists at alpha -> union of accepted groups (naive), and
  (b) the global pair winners (global).  Paired protein entrapment (FDRBench convention, entrapment.protein_paired_counts).
Verdict machinery = vulcan.stats.entrapment.evaluate_arm (1 replicate: 1 entrapment x 1 decoy seed, as pre-registered
for E9 at N = 1000; CI = pair bootstrap only).  Also: G4-style decoy score comparability across runs, L6 outputs
(master list / per-run lists / expected false counts, written to work/l6/), and runtimes.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from vulcan.stats import entrapment as ENT  # noqa: E402
from vulcan.stats import global_fdr as G  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.protein import picked_group_fdr  # noqa: E402

SEED = 11
ALPHAS = (0.001, 0.01, 0.05)
IDX = {"P": "/tmp/vulcan_e2/index/human_contam_E2001", "C": "/tmp/vulcan_m3/index/human_contam_C"}
C_JSON = "/tmp/vulcan_m1e/entrapment/db/human_contam.C.json"
MAN = os.path.join(ROOT, "data", "manifests", "e9_d8.tsv")
VARIANTS = ("score_raw", "score_cal", "score_bulk", "hyperscore")


def labels(index_dir):
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    cols = ("family", "entrap_seed", "is_decoy", "is_entrapment", "tdc_pair", "entrap_pair", "target_seq")
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in cols}


def load_arm(work, arm, runs, col="score"):
    """PSM winners (step 1, per run, by score column `col`) of the given runs -> pooled arrays."""
    parts = []
    t0 = time.perf_counter()
    for i, r in enumerate(runs):
        f = os.path.join(work, "out", arm, r + ".scored.parquet")
        if not os.path.exists(f):
            continue
        t = pq.read_table(f, columns=["spec_idx", "pf_row", "pf_id", col, "charge"])
        spec = t.column("spec_idx").to_numpy()
        pf_id = t.column("pf_id").to_numpy()
        sc = t.column(col).to_numpy().astype(np.float64)
        w = G.psm_winners(spec, sc, hash_tiebreak(pf_id, SEED))
        loc, scale = G.runner_up_stats(spec, sc)          # label-blind null proxy (runner-up candidates)
        parts.append({"run": np.full(len(w), i, np.int32), "pf_row": t.column("pf_row").to_numpy()[w],
                      "pf_id": pf_id[w], "score_raw": sc[w], "charge": t.column("charge").to_numpy()[w],
                      "score_cal": (sc[w] - loc) / scale,
                      # label-blind bulk calibration: low quantiles of the run's PSM winners (null-dominated region)
                      "score_bulk": (sc[w] - np.quantile(sc[w], 0.25)) / max(np.quantile(sc[w], 0.25)
                                                                             - np.quantile(sc[w], 0.05), 1e-6)})
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]} if parts else None
    return out, time.perf_counter() - t0, len(parts)


# ------------------------------------------------------------------------------------------- proteins
_PM = {}


def protein_map(index_dir):
    if index_dir not in _PM:
        with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
            t = rd.read_all()
        pep = t.column("peptide").to_pylist()
        off = t.column("prot_off").to_numpy()
        ln = t.column("prot_len").to_numpy()
        pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
        acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
        _PM[index_dir] = (dict(zip(pep, range(len(pep)))), off, ln, pidx, acc)
    return _PM[index_dir]


PRE = {0: "", 1: "E_", 2: "rev_", 3: "rev_E_"}


def protein_groups(index_dir, units, scores, fam, tseq, alphas):
    """Picked protein groups from pair winners (target and decoy units); returns {alpha: (groups, scores)} of
    accepted target-side groups plus the full result for entrapment counting."""
    pos, off, ln, pidx, acc = protein_map(index_dir)
    pep2prot, pep_score = {}, {}
    for u, s, f, ts in zip(units.tolist(), scores.tolist(), fam.tolist(), tseq):
        i = pos.get(ts)
        if i is None:
            continue
        pep2prot[u] = [PRE[int(f)] + acc[j] for j in pidx[off[i]:off[i] + ln[i]]]
        if s > pep_score.get(u, -np.inf):
            pep_score[u] = s
    allp = {p for v in pep2prot.values() for p in v}
    pair_of = {p: (p[4:] if p.startswith("rev_") else p) for p in allp}
    is_dec = {p: p.startswith("rev_") for p in allp}
    res = picked_group_fdr(pep2prot, pep_score, pair_of, is_dec, seed=SEED)
    out = {}
    tside = ~res.is_decoy & ~res.is_obsolete
    for a in alphas:
        accm = (res.q <= a) & tside
        out[a] = ([res.groups[i] for i in np.flatnonzero(accm)], res.score[accm])
    return out


def _run_proteins(args):
    index_dir, units, scores, fam, tseq, alphas = args
    return protein_groups(index_dir, units, scores, fam, tseq, alphas)


def protein_entrapment(groups, scores, alpha):
    allp = {p for g in groups for p in g}
    is_e = {p: p.startswith("E_") for p in allp}
    partner = {p: p[2:] for p in allp if p.startswith("E_")}
    c = ENT.protein_paired_counts(groups, scores, np.ones(len(groups), bool), is_e, partner)
    return ENT.evaluate_arm([ENT.Replicate(c, 2001, 1001)], alpha, "paired", n_boot=1000)


def union_groups(lists):
    best = {}
    for groups, scores in lists:
        for g, s in zip(groups, scores):
            if s > best.get(g, -np.inf):
                best[g] = s
    gs = list(best)
    return gs, np.array([best[g] for g in gs])


def slim(e):
    return {k: e[k] for k in ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "mean_n_t", "verdict")}


# ------------------------------------------------------------------------------------------- main
def analyze_arm(arm, d, lab, nruns, r_c, do_protein, pool):
    rows = d["pf_row"]
    is_dec = lab["is_decoy"][rows]
    is_ent = lab["is_entrapment"][rows]
    pair = lab["tdc_pair"][rows]
    unit = d["pf_id"]
    tb = hash_tiebreak(unit, SEED)
    t0 = time.perf_counter()
    glob = G.global_competition(d["run"], unit, pair, d["score"], is_dec, tb)
    t_glob = time.perf_counter() - t0
    t0 = time.perf_counter()
    per = G.run_competitions(d["run"], unit, pair, d["score"], is_dec, tb)
    t_per = time.perf_counter() - t0
    res = {"timing_s": {"global_competition": t_glob, "per_run_competitions": t_per}, "peptide": {}}
    ent_pair = lab["entrap_pair"]
    for a in ALPHAS:
        rec = {}
        # global
        gt = ~glob.is_decoy
        grow = rows[glob.row[gt]]
        acc = glob.q[gt] <= a
        if arm == "P":
            cg = ENT.paired_counts(ent_pair[grow], is_ent[glob.row[gt]], glob.score[gt], acc)
        else:
            cg = ENT.combined_counts(glob.unit[gt], is_ent[glob.row[gt]], acc)
        # naive union
        u, s, nr = G.naive_union(per, a)
        # entrapment flags per union unit: look up a representative row
        pt = per.accepted(a)
        prow = per.row[pt]
        o = np.argsort(per.unit[pt], kind="stable")
        uu = per.unit[pt][o]
        rep = prow[o][np.searchsorted(uu, u)] if len(u) else np.zeros(0, np.int64)
        if arm == "P":
            cn = ENT.paired_counts(ent_pair[rows[rep]], is_ent[rep], s, np.ones(len(u), bool))
        else:
            cn = ENT.combined_counts(u, is_ent[rep], np.ones(len(u), bool))
        est, r = ("paired", 1.0) if arm == "P" else ("combined", r_c)
        eg = ENT.evaluate_arm([ENT.Replicate(cg, 2001, 1001)], a, est, r=r, n_boot=1000)
        en = ENT.evaluate_arm([ENT.Replicate(cn, 2001, 1001)], a, est, r=r, n_boot=1000)
        n_t_glob = int((acc & ~is_ent[glob.row[gt]]).sum())
        n_t_naive = int((~is_ent[rep]).sum())
        per_run_n = np.bincount(per.run[pt].astype(np.int64), minlength=nruns)
        per_run_n_t = np.bincount(per.run[pt][~is_ent[prow]].astype(np.int64), minlength=nruns)
        rec["global"] = dict(slim(eg), n_targets=n_t_glob, n_entrapment=int((acc & is_ent[glob.row[gt]]).sum()))
        rec["naive"] = dict(slim(en), n_targets=n_t_naive, n_entrapment=int(is_ent[rep].sum()),
                            mean_runs_per_unit=float(nr.mean()) if len(nr) else 0.0)
        rec["per_run_list_targets_mean"] = float(per_run_n_t.mean())
        rec["per_run_list_size_median"] = float(np.median(per_run_n))
        res["peptide"][str(a)] = rec
    # Fair-coin diagnostic on the entrapment pairs (L4 labels): above the global / per-run threshold, entrapment
    # winners N_E vs decoy-of-entrapment winners D_E (exchangeable by construction -> N_E/(N_E+D_E) ~ 1/2), and the
    # decoys of original targets D_T.  If N_E ~ D_E while D_E/(D_E+D_T) grows with N, the rising entrapment lower
    # bound reflects E-type nulls taking a growing share of the nulls (E is mass-matched to present targets), not a
    # broken coin.
    fam = lab["family"]
    res["family_counts"] = {}
    for a in ALPHAS:
        acc = glob.accepted(a)
        s_cut = glob.score[acc].min() if acc.any() else np.inf
        above = glob.score >= s_cut
        f = fam[rows[glob.row[above]]]
        cnt = {int(k): int(v) for k, v in zip(*np.unique(f, return_counts=True))}
        n_e, d_e, d_t = cnt.get(1, 0), cnt.get(3, 0), cnt.get(2, 0)
        pa = per.accepted(a)
        pf = fam[rows[per.row]]
        # per-run lists: decoys at/above each run's own threshold
        cut_r = np.full(nruns, np.inf)
        np.minimum.at(cut_r, per.run[pa].astype(np.int64), per.score[pa])
        abv = per.score >= cut_r[per.run.astype(np.int64)]
        pc = {int(k): int(v) for k, v in zip(*np.unique(pf[abv], return_counts=True))}
        res["family_counts"][str(a)] = {
            "global": {"T": cnt.get(0, 0), "E": n_e, "dec_T": d_t, "dec_E": d_e,
                       "coin_E_vs_decE": n_e / max(n_e + d_e, 1), "share_decE_of_decoys": d_e / max(d_e + d_t, 1)},
            "per_run_sum": {"T": pc.get(0, 0), "E": pc.get(1, 0), "dec_T": pc.get(2, 0), "dec_E": pc.get(3, 0),
                            "coin_E_vs_decE": pc.get(1, 0) / max(pc.get(1, 0) + pc.get(3, 0), 1),
                            "share_decE_of_decoys": pc.get(3, 0) / max(pc.get(3, 0) + pc.get(2, 0), 1)}}
    master, perl, exp = G.global_outputs(glob, per, ALPHAS)
    res["expected_false"] = exp
    res["n_pair_winners_global"] = int(len(glob.row))
    res["n_psm_winners"] = int(len(rows))
    # per-run calibration diagnostics: run-level score cutoff at q_run <= 1%
    cut = []
    for rr in range(nruns):
        m = (per.run == rr) & per.accepted(0.01)
        if m.sum() >= 20:
            cut.append(float(per.score[m].min()))
    res["run_cutoff_1pct"] = {"n_runs": len(cut), "median": float(np.median(cut)) if cut else None,
                              "iqr": [float(np.quantile(cut, 0.25)), float(np.quantile(cut, 0.75))] if cut else None,
                              "min": float(min(cut)) if cut else None, "max": float(max(cut)) if cut else None}
    gcut = glob.score[glob.accepted(0.01)]
    res["global_cutoff_1pct"] = float(gcut.min()) if len(gcut) else None
    if do_protein and arm == "P":
        t0 = time.perf_counter()
        fam = lab["family"]
        tseq = lab["target_seq"]
        gp = protein_groups(IDX["P"], glob.unit, glob.score, fam[rows[glob.row]], tseq[rows[glob.row]], ALPHAS)
        t_gp = time.perf_counter() - t0
        t0 = time.perf_counter()
        jobs = []
        for rr in range(nruns):
            m = per.run == rr
            if not m.any():
                continue
            rw = rows[per.row[m]]
            jobs.append((IDX["P"], per.unit[m], per.score[m], fam[rw], tseq[rw], ALPHAS))
        runp = list(pool.map(_run_proteins, jobs, chunksize=1)) if pool else [_run_proteins(j) for j in jobs]
        t_rp = time.perf_counter() - t0
        res["protein"] = {}
        for a in ALPHAS:
            g_groups, g_scores = gp[a]
            ng, ns = union_groups([rp[a] for rp in runp])
            eg = protein_entrapment(g_groups, g_scores, a) if g_groups else None
            en = protein_entrapment(ng, ns, a) if ng else None
            is_t = lambda gs: int(sum(any(not p.startswith("E_") for p in g) for g in gs))  # noqa: E731
            res["protein"][str(a)] = {
                "global": dict(slim(eg) if eg else {}, n_groups=len(g_groups), n_target_groups=is_t(g_groups)),
                "naive": dict(slim(en) if en else {}, n_groups=len(ng), n_target_groups=is_t(ng)),
                "per_run_groups_mean": float(np.mean([len(rp[a][0]) for rp in runp])) if runp else 0.0}
        res["timing_s"]["protein_global"] = t_gp
        res["timing_s"]["protein_per_run_all"] = t_rp
    return res, glob, per, master, perl


def confirm_first_runs(work, want, labs, r_c, n=10):
    """DESIGN §6.0 A1.1: dev numbers of the frozen default on runs not used for selection (first 10 D8 runs).
    Per run: peptidoform IDs at 1% and per-run entrapment (arm P paired / arm C combined, 1 replicate), for
    static_gbm_pred and the hyperscore reference, from the same augmented searches."""
    out = {}
    for arm in ("P", "C"):
        lab = labs[arm]
        for scorer, col in (("static_gbm_pred", "score"), ("hyperscore", "hyperscore")):
            rows_out = []
            for r in want[:n]:
                f = os.path.join(work, "out", arm, r + ".scored.parquet")
                if not os.path.exists(f):
                    continue
                t = pq.read_table(f)
                spec, pf_id = t.column("spec_idx").to_numpy(), t.column("pf_id").to_numpy()
                sc = t.column(col).to_numpy().astype(np.float64)
                tb = hash_tiebreak(pf_id, SEED)
                w = G.psm_winners(spec, sc, tb)
                rows = t.column("pf_row").to_numpy()[w]
                dec = lab["is_decoy"][rows]
                g = G.global_competition(np.zeros(len(w), np.int32), pf_id[w], lab["tdc_pair"][rows], sc[w], dec, tb[w])
                tt = ~g.is_decoy
                acc = g.q[tt] <= 0.01
                ent = lab["is_entrapment"][rows[g.row[tt]]]
                if arm == "P":
                    c = ENT.paired_counts(lab["entrap_pair"][rows[g.row[tt]]], ent, g.score[tt], acc)
                    e = ENT.evaluate_arm([ENT.Replicate(c, 2001, 1001)], 0.01, "paired", n_boot=500)
                else:
                    c = ENT.combined_counts(g.unit[tt], ent, acc)
                    e = ENT.evaluate_arm([ENT.Replicate(c, -1, 1001)], 0.01, "combined", r=r_c, n_boot=500)
                rows_out.append({"run": r, "ids_1pct": int((acc & ~ent).sum()), **slim(e)})
            out[f"{arm}_{scorer}"] = {"runs": rows_out,
                                      "mean_fdp": float(np.mean([x["mean_fdp"] for x in rows_out])) if rows_out else None,
                                      "sum_ids_1pct": int(sum(x["ids_1pct"] for x in rows_out))}
    return out


def comparability(d, lab, run_names, instr):
    rows = d["pf_row"]
    dec = lab["is_decoy"][rows]
    out = {"per_run_charge": G.score_comparability(d["run"], d["score"], dec, strata=np.minimum(d["charge"], 4)),
           "per_run": G.score_comparability(d["run"], d["score"], dec)}
    inst = np.array([instr[run_names[i]] for i in d["run"]])
    out["per_instrument"] = G.score_comparability(inst, d["score"], dec, min_n=500)
    # same check on targets of the low-score bulk is NOT done (labels: decoys only, §6.4 G4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=int, nargs="*", default=[1, 10, 100])
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--arms", nargs="+", default=["P", "C"])
    ap.add_argument("--no-protein-above", type=int, default=100000)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--confirm", action="store_true", help="A1.1 per-run confirmation on the first 10 runs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    man = list(csv.DictReader(open(MAN), delimiter="\t"))
    import re
    name = {r["id"]: f"{r['id']}__{re.sub(r'[^A-Za-z0-9._-]', '_', r['filename'])[:-4]}" for r in man}
    instr = {name[r["id"]]: r["instrument"] for r in man}
    proj = {name[r["id"]]: r["accession"] for r in man}
    r_c = float(json.load(open(C_JSON))["r"])
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    out.setdefault("meta", {})
    out["meta"].update({"seed": SEED, "alphas": ALPHAS, "index": IDX, "r_C": r_c, "unit": "peptidoform (pf_id)",
                        "score": "static_gbm_pred (frozen, M3 D7)", "replicates": "1 entrapment (2001) x 1 decoy (1001)"})
    labs = {arm: labels(IDX[arm]) for arm in a.arms}
    pool = ProcessPoolExecutor(a.procs)
    if a.confirm:
        want = [name[r["id"]] for r in man[:10]]
        out["A1_1_confirmation_first10"] = confirm_first_runs(a.work, want, labs, r_c)
        json.dump(out, open(a.out, "w"), indent=1, default=float)
    for n in a.ns:
        want = [name[r["id"]] for r in man[:n]]
        for arm in a.arms:
            d, t_load, n_present = load_arm(a.work, arm, want)
            if d is None:
                continue
            present = [r for r in want if os.path.exists(os.path.join(a.work, "out", arm, r + ".scored.parquet"))]
            # run index in d = position in `want`
            # G4-style comparability (decoys, dev data) of the raw frozen score and of two label-blind per-run
            # calibrations.  Primary = the raw frozen score (D7); calibrations are diagnostics (see NOTES: the KS-min
            # variant 'score_bulk' loses IDs at equal validity, and no label-blind variant passes G4).
            comp = {}
            dh = None
            if "hyperscore" in a.variants:       # fixed hand-crafted reference score (own spectrum-level winners)
                dh, _, _ = load_arm(a.work, arm, want, col="hyperscore")
                dh["hyperscore"] = dh["score_raw"]
            for v in a.variants:
                dd = dh if v == "hyperscore" else d
                dd["score"] = dd[v]
                comp[v] = comparability(dd, labs[arm], want, instr)
            ks = {v: comp[v]["per_run"].get("summary", {}).get("excess_median", 0.0) for v in a.variants}
            primary = "score_raw"
            out.setdefault(f"N{n}", {})[arm] = {"comparability": comp, "primary": primary, "ks_excess_median": ks,
                                                "ks_min_variant": min(a.variants, key=lambda v: ks[v])}
            for v in a.variants:
                dd = dh if v == "hyperscore" else d
                dd["score"] = dd[v]
                t0 = time.perf_counter()
                res, glob, per, master, perl = analyze_arm(arm, dd, labs[arm], len(want), r_c,
                                                           do_protein=n <= a.no_protein_above, pool=pool)
                res["timing_s"]["load"] = t_load
                res["timing_s"]["analysis_total"] = time.perf_counter() - t0
                out[f"N{n}"][arm][v] = res
                p1 = res["peptide"]["0.01"]
                print(f"N={n} arm={arm} {v} runs={n_present} pep@1%: global FDP={p1['global']['mean_fdp']:.4f} "
                      f"(T={p1['global']['n_targets']}) naive FDP={p1['naive']['mean_fdp']:.4f} (T={p1['naive']['n_targets']})"
                      + (f" | prot@1%: global {res['protein']['0.01']['global'].get('mean_fdp')} naive "
                         f"{res['protein']['0.01']['naive'].get('mean_fdp')}" if "protein" in res else ""), flush=True)
                if arm == "P" and v == primary:
                    l6 = os.path.join(a.work, "l6", f"N{n}")
                    os.makedirs(l6, exist_ok=True)
                    pq.write_table(pa.table({"pf_id": master["unit"], "best_run": np.asarray(want)[master["best_run"]],
                                             "score": master["score"], "q_global": master["q_global"]}),
                                   os.path.join(l6, "master.parquet"))
                    pq.write_table(pa.table({"run": np.asarray(want)[perl["run"]], "pf_id": perl["unit"],
                                             "q_run": perl["q_run"], "q_global": perl["q_global"]}),
                                   os.path.join(l6, "per_run.parquet"))
            res = out[f"N{n}"][arm]
            res["n_runs_requested"] = n
            res["n_runs_present"] = n_present
            res["n_projects"] = len({proj[r] for r in present})
            json.dump(out, open(a.out, "w"), indent=1, default=float)
    pool.shutdown()


if __name__ == "__main__":
    main()
