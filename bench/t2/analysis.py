"""E8 analysis for the T2 driver: the M4 E8 computation (scripts/m4/e8.py functions, imported unchanged) plus the
pre-registered additions that M4 did not compute (DESIGN §6.5 close-acceptor strata, D8 site-level FLR̂) and the
uncertainty used by the claim rules in config/t2.yaml (Clopper-Pearson and sequence-cluster bootstrap UCB of true FLR,
paired cluster bootstrap of the recall difference vs pyAscore, recall at estimated FLR).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
M4 = os.path.join(ROOT, "scripts", "m4")
if M4 not in sys.path:
    sys.path.insert(0, M4)

import common as C  # noqa: E402  (scripts/m4)
import e8 as M4E8  # noqa: E402  (scripts/m4/e8.py, unchanged)

from vulcan.stats import flr as FLR  # noqa: E402

ALPHAS = (0.01, 0.05)


def set_pyascore_path(pylib: str):
    """scripts/m4/e8.run_pyascore imports pyascore from <C.WORK>/pylib; point it at the pinned persistent copy."""
    C.WORK = os.path.dirname(os.path.abspath(pylib))
    assert os.path.basename(pylib.rstrip("/")) == "pylib"


# ------------------------------------------------------------------------------------------------ PSM list (M4 E8 step 1-2)
def psm_list(acc_runs: list[str], psm_dir: str, index: str, l1, q: float = 0.01, seed: int = 11):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition

    tbl = pa.concat_tables([pq.read_table(os.path.join(psm_dir, f"{r}.psms.parquet")) for r in acc_runs])
    tbl = tbl.filter(pc.equal(tbl.column("rank"), 1))
    rows = tbl.column("pf_row").to_numpy()
    lab = M4E8.read_labels(index, rows)
    seqs, phos, oxs, nts = M4E8.decode(l1, rows)
    comp = [f"p{len(ph)}o{len(ox)}n{int(nt)}" for ph, ox, nt in zip(phos, oxs, nts)]
    unit = np.array([f"{s}|{c}" for s, c in zip(seqs, comp)])
    pair = np.array([f"{f}|{c}" for f, c in zip(lab["fwd_seq"], comp)])
    spec = np.array([f"{r}:{i}" for r, i in zip(tbl.column("run").to_pylist(), tbl.column("spec_idx").to_pylist())])
    r = double_competition(spec, unit, pair, tbl.column("hyperscore").to_numpy(), lab["is_decoy"], seed=seed,
                           tiebreak=hash_tiebreak(unit, seed))
    acc_units = set(unit[r.accepted(q)].tolist())
    keep = np.array([(u in acc_units) and (not d) and len(ph) > 0 for u, d, ph in zip(unit, lab["is_decoy"], phos)])
    idx = np.flatnonzero(keep)
    run_col = np.array(tbl.column("run").to_pylist())
    psm = {"seq": [seqs[i] for i in idx], "phos": [phos[i] for i in idx], "ox": [oxs[i] for i in idx],
           "nterm": [nts[i] for i in idx], "charge": tbl.column("charge").to_numpy()[idx],
           "run": run_col[idx], "spec_idx": tbl.column("spec_idx").to_numpy()[idx]}
    n_target_units = int(len(acc_units))
    per_run = {}
    for rn in acc_runs:
        m = psm["run"] == rn
        per_run[rn] = {"n_phospho_psms_localized": int(m.sum())}
    counts = {"n_rank1_psms": int(len(rows)), "n_peptidoforms_q01": int(len(r.accepted(q))),
              "n_target_units_q01": n_target_units, "n_phospho_psms_localized": int(len(idx)),
              "n_mono": int(sum(len(x) == 1 for x in psm["phos"])), "n_multi": int(sum(len(x) > 1 for x in psm["phos"])),
              "per_run": per_run}
    return psm, counts


# ------------------------------------------------------------------------------------------------ statistics helpers
def cp_upper(k: int, n: int, conf: float = 0.95) -> float | None:
    from scipy.stats import beta
    if n == 0:
        return None
    if k >= n:
        return 1.0
    return float(beta.ppf(conf, k + 1, n - k))


def _cluster_weights(rng, n_clusters: int, n_boot: int) -> np.ndarray:
    """(n_boot, n_clusters) multiplicities of a cluster bootstrap (resample clusters with replacement)."""
    draws = rng.integers(0, n_clusters, size=(n_boot, n_clusters))
    w = np.zeros((n_boot, n_clusters), np.int32)
    for b in range(n_boot):
        w[b] = np.bincount(draws[b], minlength=n_clusters)
    return w


def _recall_matched(score_sorted, sty_sorted, correct_sorted, last_tie, w_sites, n_true_w, alpha):
    """Weighted version of M4 evaluate's recall at matched true FLR (weights = cluster multiplicities)."""
    wrong = np.cumsum(w_sites * (sty_sorted & ~correct_sorted), axis=-1)
    nsty = np.cumsum(w_sites * sty_sorted, axis=-1)
    good = np.cumsum(w_sites * (sty_sorted & correct_sorted), axis=-1)
    tf = np.where(nsty > 0, wrong / np.maximum(nsty, 1), 0.0)
    ok = last_tie & (tf <= alpha)
    g = np.where(ok, good, 0).max(axis=-1)
    return np.where(n_true_w > 0, g / np.maximum(n_true_w, 1), np.nan)


class ToolSites:
    """Mono-phospho sites of one tool and decoy residue, with everything needed for the bootstrap."""

    def __init__(self, sites, psm, truth_of, decoy, cluster_of_seq):
        self.ok = bool(sites)
        if not self.ok:
            return
        pi = np.array([s[0] for s in sites])
        pos = np.array([s[1] for s in sites])
        sc = np.array([s[2] for s in sites])
        seqs = psm["seq"]
        res = np.array([seqs[i][q] for i, q in zip(pi, pos)])
        nph = np.array([len(psm["phos"][i]) for i in pi])
        msk = nph == 1
        self.pi, self.pos, self.sc, self.res = pi[msk], pos[msk], sc[msk], res[msk]
        self.f = FLR.flr_curve(self.sc, self.res, self.pi, seqs, decoy_residue=decoy)   # = M4 evaluate's f
        rr = np.where(np.isin(self.res, ["I"]) & (decoy == "L"), "L", self.res)
        self.is_sty = np.isin(rr, list("STY"))
        tp = [truth_of[i] for i in self.pi]
        self.elig = np.array([t is not None for t in tp])
        self.correct = np.array([t is not None and t == q for t, q in zip(tp, self.pos)])
        self.cl = np.array([cluster_of_seq[seqs[i]] for i in self.pi])
        ep = sorted({int(i) for i in self.pi[self.elig]})
        self.elig_psm_cl = np.array([cluster_of_seq[seqs[i]] for i in ep], dtype=np.int64)
        # matched-true-FLR arrays (eligible sites, sorted by score desc, stable = M4)
        e = self.elig
        o = np.argsort(-self.sc[e], kind="stable")
        self.m_sc = self.sc[e][o]
        self.m_sty = self.is_sty[e][o]
        self.m_cor = self.correct[e][o]
        self.m_cl = self.cl[e][o]
        self.m_last = np.r_[self.m_sc[1:] != self.m_sc[:-1], True]
        self.decoy = decoy
        self.seqs = seqs
        self.nph_all = nph

    def recall_matched(self, alpha, W=None):
        if W is None:
            ws = np.ones(len(self.m_sc))
            nt = float(len(self.elig_psm_cl))
            return float(_recall_matched(self.m_sc, self.m_sty, self.m_cor, self.m_last, ws, nt, alpha))
        ws = W[:, self.m_cl]
        nt = W[:, self.elig_psm_cl].sum(axis=1)
        return _recall_matched(self.m_sc, self.m_sty, self.m_cor, self.m_last, ws, nt, alpha)

    def at_estimated(self, alpha, W=None):
        """(true FLR, recall) at FLR̂ <= alpha (threshold from all mono sites), on eligible sites; weighted if W."""
        acc = (self.f <= alpha) & self.is_sty & self.elig
        if W is None:
            n = int(acc.sum())
            k = int((acc & ~self.correct).sum())
            g = int((acc & self.correct).sum())
            nt = len(self.elig_psm_cl)
            return {"n": n, "k_wrong": k, "true_flr": (k / n) if n else None, "correct": g,
                    "recall": (g / nt) if nt else None}
        ws = W[:, self.cl]
        n = (ws * acc).sum(axis=1)
        k = (ws * (acc & ~self.correct)).sum(axis=1)
        g = (ws * (acc & self.correct)).sum(axis=1)
        nt = W[:, self.elig_psm_cl].sum(axis=1)
        return {"true_flr": np.where(n > 0, k / np.maximum(n, 1), np.nan), "recall": np.where(nt > 0, g / np.maximum(nt, 1), np.nan)}

    def strata(self, alpha):
        """§6.5: sites with / without another candidate acceptor within 2 residues (global FLR̂ threshold)."""
        st = FLR.flr_strata(np.ones(len(self.pos), int), self.pos, self.seqs, self.pi, decoy_residue=self.decoy)
        out = {}
        for name, m in (("close", st["close_site"]), ("not_close", ~st["close_site"])):
            acc = (self.f <= alpha) & self.is_sty & m
            acc_e = acc & self.elig
            n = int(acc_e.sum())
            own = FLR.flr_curve(self.sc[m], self.res[m], self.pi[m], self.seqs, decoy_residue=self.decoy) if m.any() else np.array([])
            acc_own = (own <= alpha) & self.is_sty[m] & self.elig[m]
            n2 = int(acc_own.sum())
            out[name] = {"n_sites": int(m.sum()), "n_sty_at_global_thr": int(acc.sum()), "n_eligible_at_global_thr": n,
                         "true_flr_at_global_thr": float((acc_e & ~self.correct).sum() / n) if n else None,
                         "n_eligible_at_own_stratum_thr": n2,
                         "true_flr_at_own_stratum_thr": float((acc_own & ~self.correct[m]).sum() / n2) if n2 else None}
        return out

    def site_level(self, alpha):
        """D8 secondary: unit = distinct (sequence, position); score = max over PSMs; acceptor counts per distinct
        sequence; truth = the sequence's truth position (if any PSM of it is eligible)."""
        key = {}
        for j, (i, q) in enumerate(zip(self.pi, self.pos)):
            k = (self.seqs[i], int(q))
            if k not in key or self.sc[j] > self.sc[key[k]]:
                key[k] = j
        js = np.array([key[k] for k in sorted(key)])
        seq_ids = {}
        pidx = np.array([seq_ids.setdefault(self.seqs[self.pi[j]], len(seq_ids)) for j in js])
        uniq = sorted(seq_ids, key=seq_ids.get)
        f = FLR.flr_curve(self.sc[js], self.res[js], pidx, uniq, decoy_residue=self.decoy)
        sty, el, cor = self.is_sty[js], self.elig[js], self.correct[js]
        out = {"n_sites": int(len(js))}
        for a in ALPHAS:
            acc = (f <= a) & sty
            ae = acc & el
            n = int(ae.sum())
            out[f"{a}"] = {"n_sty_sites": int(acc.sum()), "n_eligible": n,
                           "true_flr": float((ae & ~cor).sum() / n) if n else None}
        return out


def cluster_index(psm) -> dict:
    return {s: i for i, s in enumerate(sorted(set(psm["seq"])))}


def run_ours_deterministic(psm, spectra, acc_set, loc_params):
    """NOT USED by the driver since PI decision D35 (the look runs the frozen, non-deterministic path); kept for the
    dev diagnosis in results/2026-10-28_T2_driver_dryrun.  scripts/m4/e8.run_ours with torch deterministic algorithms on.  The CUDA index_add_ in
    vulcan.ptm.localize.score_configs is otherwise atomic and order-dependent, so tied configuration scores (e.g. S vs
    A with identical ion sets) resolve differently from run to run (measured on dev: n accepted +-3, true FLR +-0.2 pp).
    Deterministic mode changes only the summation order (no parameter), and is switched off again afterwards so the
    search engine is untouched (IMPLEMENTATION_NOTES §3)."""
    import torch
    prev = torch.are_deterministic_algorithms_enabled()
    torch.use_deterministic_algorithms(True)
    try:
        return M4E8.run_ours(psm, spectra, acc_set, loc_params)
    finally:
        torch.use_deterministic_algorithms(prev)


def e8_component(acc, psm, truth_of, spectra, loc_params, run_pyascore: bool, boot_n: int, boot_seed: int,
                 n_min: int, log=print, repeats: int = 0):
    """M4 E8 per accession (unchanged eval) + T2 additions.  Returns the metrics dict for this component."""
    ds = {"n_truth_eligible_psms": int(sum(t is not None for t in truth_of))}
    own_correct = [t is not None and len(ph) == 1 and ph[0] == t for t, ph in zip(truth_of, psm["phos"])]
    ds["search_rank1_site_correct_frac_on_eligible"] = (
        float(np.mean([c for c, t in zip(own_correct, truth_of) if t is not None]))
        if any(t is not None for t in truth_of) else None)
    cl = cluster_index(psm)
    tools = {}
    for decoy in ("A", "L"):
        acc_set = "STY" + decoy
        sites, dt, nskip = M4E8.run_ours(psm, spectra, acc_set, loc_params)     # frozen path, unchanged (PI D35)
        ds[f"ours_{decoy}"] = {"eval": M4E8.evaluate(sites, psm, truth_of, decoy), "time_s": dt, "n_skipped": nskip,
                               "psms_per_s": len(psm["seq"]) / dt if dt else None}
        tools[("ours", decoy)] = ToolSites(sites, psm, truth_of, decoy, cl)
        log(f"{acc} ours_{decoy} done")
        if run_pyascore:
            try:
                sites, dt, nf = M4E8.run_pyascore(psm, spectra, acc_set, loc_params)
                ds[f"pyascore_{decoy}"] = {"eval": M4E8.evaluate(sites, psm, truth_of, decoy), "time_s": dt,
                                           "n_failed": nf, "psms_per_s": len(psm["seq"]) / dt if dt else None}
                tools[("pyascore", decoy)] = ToolSites(sites, psm, truth_of, decoy, cl)
            except Exception as e:  # completeness: recorded, comparator counts as not evaluable
                ds[f"pyascore_{decoy}"] = {"error": type(e).__name__}
            log(f"{acc} pyascore_{decoy} done")
    ds["t2"] = t2_additions(tools, cl, boot_n, boot_seed, n_min, ds["n_truth_eligible_psms"])
    # run-to-run variability of the frozen localizer (CUDA atomic sums break exact ties at random; PI D35): repeat our
    # localization `repeats` more times on the same PSMs and spectra (pyAscore is deterministic and is reused); every
    # claim number of run 0 is printed with its [min, max] over all runs.  Run 0 is the result of record.
    reps = [ds["t2"]]
    for k in range(repeats):
        tk = {key: v for key, v in tools.items() if key[0] == "pyascore"}
        for decoy in ("A", "L"):
            sites, _, _ = M4E8.run_ours(psm, spectra, "STY" + decoy, loc_params)
            tk[("ours", decoy)] = ToolSites(sites, psm, truth_of, decoy, cl)
        reps.append(t2_additions(tk, cl, boot_n, boot_seed, n_min, ds["n_truth_eligible_psms"]))
        log(f"{acc} repeat {k + 1}/{repeats} done")
    ds["t2_band"] = band(reps)
    return ds


def _leaves(d, pre=()):
    if isinstance(d, dict):
        for k, v in d.items():
            yield from _leaves(v, pre + (k,))
    else:
        yield pre, d


def band(reps):
    """{path: {'min','max','n'}} for numeric leaves, {path: {'counts': {value: n}}} for verdict strings."""
    acc = {}
    for r in reps:
        for path, v in _leaves(r):
            if path[-1] in ("ci95",) or path[0] == "bootstrap":
                continue
            acc.setdefault("/".join(path), []).append(v)
    out = {}
    for k, vs in acc.items():
        if len(vs) != len(reps):
            out[k] = {"incomplete": len(vs)}
        elif all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in vs):
            out[k] = {"min": min(vs), "max": max(vs), "n": len(vs)}
        elif all(isinstance(v, str) for v in vs) and k.endswith("verdict"):
            c = {}
            for v in vs:
                c[v] = c.get(v, 0) + 1
            out[k] = {"counts": c, "n": len(vs)}
    return out


def t2_additions(tools, cl, boot_n, boot_seed, n_min, n_elig_psms):
    out = {"bootstrap": {"n": boot_n, "seed": boot_seed, "n_clusters": len(cl)}}
    rng = np.random.default_rng(boot_seed)
    W = _cluster_weights(rng, max(len(cl), 1), boot_n) if len(cl) else None

    def chunked(fn, *args):     # bounded memory: evaluate the bootstrap in blocks of 200 replicates
        parts = [fn(*args, W[b:b + 200]) for b in range(0, W.shape[0], 200)]
        if isinstance(parts[0], dict):
            return {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
        return np.concatenate(parts)
    for (tool, decoy), T in sorted(tools.items()):
        key = f"{tool}_{decoy}"
        if not T.ok:
            out[key] = {"no_sites": True}
            continue
        d = {}
        for a in ALPHAS:
            pt = T.at_estimated(a)
            cell = {"n_accepted_eligible": pt["n"], "k_wrong": pt["k_wrong"], "true_flr": pt["true_flr"],
                    "recall_at_estimated_flr": pt["recall"], "recall_at_true_flr": T.recall_matched(a) if n_elig_psms else None}
            if n_elig_psms and W is not None and pt["n"]:
                bs = chunked(T.at_estimated, a)
                cpu = cp_upper(pt["k_wrong"], pt["n"])
                bu = float(np.nanpercentile(bs["true_flr"], 95))
                cell.update(ucb_cp=cpu, ucb_boot=bu, ucb=max(cpu, bu),
                            n_boot_nan=int(np.isnan(bs["true_flr"]).sum()))
            cell["verdict"] = e8_verdict(cell, a, n_min, n_elig_psms)
            d[f"{a}"] = cell
            d[f"strata_{a}"] = T.strata(a)
        d["site_level"] = T.site_level(0.01)
        d["n_mono_sites"] = int(len(T.pos))
        out[key] = d
    # paired recall differences ours - pyAscore (same eligible PSMs, same cluster draws)
    for decoy in ("A", "L"):
        o, p = tools.get(("ours", decoy)), tools.get(("pyascore", decoy))
        key = f"diff_{decoy}"
        if o is None or p is None or not (o.ok and p.ok) or not n_elig_psms or W is None:
            out[key] = {"not_evaluable": "comparator output missing or no truth"}
            continue
        d = {}
        for a in ALPHAS:
            for kind, fo, fp in (("true", o.recall_matched, p.recall_matched),
                                 ("estimated", lambda a_, W_=None: o.at_estimated(a_, W_)["recall"],
                                  lambda a_, W_=None: p.at_estimated(a_, W_)["recall"])):
                pt = fo(a) - fp(a)
                bo, bp = chunked(fo, a), chunked(fp, a)
                diff = bo - bp
                nn = int(np.isnan(diff).sum())
                lo, hi = (float(np.nanpercentile(diff, 2.5)), float(np.nanpercentile(diff, 97.5))) if nn < len(diff) else (None, None)
                ver = ("not evaluable" if nn > 0.05 * len(diff) or lo is None else
                       "higher" if lo > 0 else "lower" if hi < 0 else "no difference shown")
                d[f"{kind}_{a}"] = {"ours": fo(a), "pyascore": fp(a), "diff": pt, "ci95": [lo, hi],
                                    "n_boot_nan": nn, "verdict": ver}
        out[key] = d
    return out


def e8_verdict(cell, alpha, n_min, n_elig_psms):
    if not n_elig_psms:
        return "not evaluable (no truth)"
    if cell["n_accepted_eligible"] < n_min:
        return f"not evaluable (n {cell['n_accepted_eligible']} < {n_min})"
    if cell["true_flr"] > alpha:
        return "FAILED"
    if cell.get("ucb") is None or cell["ucb"] > alpha:
        return "inconclusive"
    return "PASS"


def load_spectra(stores, keys):
    t = time.perf_counter()
    sp = M4E8.load_spectra(stores, keys)
    return sp, time.perf_counter() - t
