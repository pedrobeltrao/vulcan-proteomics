"""Monte-Carlo simulators for estimator verification (DESIGN.md §6.7).

Every simulator returns per-replicate FDP (false discoveries / max(discoveries, 1),
0 for an empty list) and discovery counts for each alpha, using the library code
paths (``tdc_qvalues`` with replicate blocks, ``double_competition``-equivalent
staged competitions, ``group_walk_batch``, ``reset``, picked protein code).
Replicates are split evenly over the pre-registered general seeds
``SEEDS = (11, 23, 47, 83, 131)``; per-seed means are reported.
"""
from __future__ import annotations

import time

import numpy as np

from ._util import SEEDS, make_rng
from .groupwalk import group_walk, group_walk_batch
from .protein import picked_group_fdr, picked_protein_fdr_simple
from .reset import reset
from .tdc import compete, multi_decoy_competition, multi_decoy_qvalues, tdc_qvalues

ALPHAS = (0.001, 0.005, 0.01, 0.05)


def _summ(fdp: np.ndarray, ndisc: np.ndarray, seed_of: np.ndarray, alpha: float) -> dict:
    R = fdp.size
    m = float(fdp.mean())
    se = float(fdp.std(ddof=1) / np.sqrt(R)) if R > 1 else float("nan")
    per_seed = {int(s): float(fdp[seed_of == s].mean()) for s in np.unique(seed_of)}
    return {"alpha": alpha, "n_rep": int(R), "mean_fdp": m, "se": se, "mean_fdp_over_alpha": m / alpha,
            "upper_99": m + 2.576 * se, "frac_rep_fdp_gt_alpha": float((fdp > alpha).mean()),
            "mean_discoveries": float(ndisc.mean()), "per_seed_mean_fdp": per_seed,
            "controls": bool(m <= alpha), "controls_within_mc_error": bool(m - 2.576 * se <= alpha)}


def _split(R: int, seeds=SEEDS):
    base, extra = divmod(R, len(seeds))
    return [(s, base + (1 if i < extra else 0)) for i, s in enumerate(seeds)]


def _fdp_from_q(q, is_decoy, is_false, block, R, alpha):
    acc = (q <= alpha) & ~is_decoy
    nd = np.bincount(block[acc], minlength=R)
    nf = np.bincount(block[acc & is_false], minlength=R)
    return nf / np.maximum(nd, 1), nd


# ----------------------------------------------------------------------------- TDC

def sim_tdc(R: int, m: int = 2000, pi0: float = 0.5, mu: float = 3.0, p_null_target: float = 0.5,
            alphas=ALPHAS, chunk: int = 500, sim_id: int = 1, true_always_win: bool = False) -> dict:
    """Independent competitions (peptide-level TDC). Null: target and decoy scores
    iid N(0,1) (winner = max, label target w.p. ``p_null_target`` independent of
    the winning score). True targets: N(mu,1) vs decoy N(0,1).

    ``true_always_win``: correct targets always beat their decoy (decoy wins come
    only from null competitions). TDC is then nearly tight (FDP ~ alpha), the
    regime in which a biased null coin (p > 1/2) inflates the FDP by ~p/(1-p)."""
    res = {a: ([], [], []) for a in alphas}
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            n0 = int(round(m * pi0))
            null = np.zeros((r, m), bool)
            null[:, :n0] = True
            a = rng.standard_normal((r, m))
            b = rng.standard_normal((r, m))
            t = np.where(null, a, a + mu)
            w = np.maximum(t, b)
            dec = b > t
            if true_always_win:
                w = np.where(null, w, t)
                dec = np.where(null, dec, False)
            if p_null_target != 0.5:
                w = np.where(null, np.maximum(a, b), w)
                dec = np.where(null, rng.random((r, m)) >= p_null_target, dec)
            blk = np.repeat(np.arange(r), m)
            q = tdc_qvalues(w.ravel(), dec.ravel(), block=blk)
            false = null.ravel() & ~dec.ravel()
            for al in alphas:
                f, nd = _fdp_from_q(q, dec.ravel(), false, blk, r, al)
                res[al][0].append(f); res[al][1].append(nd); res[al][2].append(np.full(r, seed))
    return {str(a): _summ(np.concatenate(v[0]), np.concatenate(v[1]), np.concatenate(v[2]), a)
            for a, v in res.items()}


# --------------------------------------------- spectrum-level + double competition

def _gen_spectra(rng, n_pairs, n_present, mu, n_foreign, mult_present, mult_foreign, n_cand,
                 sigma_analyte, sigma_spec):
    """One replicate of candidate rows. Returns dict of arrays (rows) and spectrum truth."""
    # present peptides: spectra counts
    kp = 1 + rng.poisson(mult_present, n_present)
    kf = 1 + rng.poisson(mult_foreign, n_foreign)
    sp_pep = np.repeat(np.arange(n_present), kp)          # true pair of present spectra
    sp_for = np.repeat(np.arange(n_foreign), kf)          # foreign analyte of foreign spectra
    ns_p, ns_f = sp_pep.size, sp_for.size
    # candidate pairs
    cand_p = rng.integers(0, n_pairs, size=(ns_p, n_cand))
    cand_p[:, 0] = sp_pep
    f_cands = rng.integers(0, n_pairs, size=(n_foreign, n_cand))
    cand_f = f_cands[sp_for]
    # scores (per pair side): present spectra: independent per-spectrum effects
    sp_t = rng.standard_normal((ns_p, n_cand)) * np.hypot(sigma_analyte, sigma_spec)
    sp_d = rng.standard_normal((ns_p, n_cand)) * np.hypot(sigma_analyte, sigma_spec)
    sp_t[:, 0] += mu
    # foreign spectra: analyte-level shared effect + spectrum noise (multiplicity)
    a_t = rng.standard_normal((n_foreign, n_cand)) * sigma_analyte
    a_d = rng.standard_normal((n_foreign, n_cand)) * sigma_analyte
    sf_t = a_t[sp_for] + rng.standard_normal((ns_f, n_cand)) * sigma_spec
    sf_d = a_d[sp_for] + rng.standard_normal((ns_f, n_cand)) * sigma_spec
    cand = np.vstack([cand_p, cand_f])
    st = np.vstack([sp_t, sf_t])
    sd = np.vstack([sp_d, sf_d])
    nspec = ns_p + ns_f
    true_pair = np.r_[sp_pep, np.full(ns_f, -1)]
    spec = np.repeat(np.arange(nspec), 2 * n_cand)
    pair = np.repeat(cand, 2, axis=1).ravel()
    side = np.tile(np.array([0, 1]), (nspec, n_cand)).ravel()  # 1 = decoy
    score = np.stack([st, sd], axis=2).reshape(nspec, -1).ravel()
    return spec, pair, side.astype(bool), score, true_pair


def sim_double_competition(R: int, n_pairs: int = 1500, n_present: int = 450, mu: float = 3.5,
                           n_foreign: int = 300, mult_present: float = 1.0, mult_foreign: float = 2.0,
                           n_cand: int = 2, sigma_analyte: float = 0.85, sigma_spec: float = 0.5,
                           alphas=ALPHAS, chunk: int = 200, sim_id: int = 2) -> dict:
    """Spectra -> PSM-level competition (all candidates of a spectrum) -> best PSM
    per peptide -> paired target/decoy peptide competition -> (D+1)/T.

    Foreign analytes (absent from the DB) produce several spectra that share an
    analyte-level random effect on the same candidates (multiplicity, Lin et al.
    2022). Peptide-level FDP is evaluated for double competition; PSM-level FDP
    (a PSM is false if its target is not the spectrum's true peptide) is reported
    as "not FDR-controlled" context.
    """
    out_pep = {a: ([], [], []) for a in alphas}
    out_psm = {a: ([], [], []) for a in alphas}
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            parts = [_gen_spectra(rng, n_pairs, n_present, mu, n_foreign, mult_present, mult_foreign, n_cand,
                                  sigma_analyte, sigma_spec) for _ in range(r)]
            offs_s = np.cumsum([0] + [p[4].size for p in parts])
            spec = np.concatenate([p[0] + offs_s[i] for i, p in enumerate(parts)])
            rep = np.concatenate([np.full(p[0].size, i) for i, p in enumerate(parts)])
            pair = np.concatenate([p[1] for p in parts]) + rep * n_pairs
            dec = np.concatenate([p[2] for p in parts])
            score = np.concatenate([p[3] for p in parts])
            true_pair = np.concatenate([np.where(p[4] >= 0, p[4] + i * n_pairs, -1) for i, p in enumerate(parts)])
            tb = rng.random(score.size)
            # stage 1: PSM competition (all candidates/hypotheses of a spectrum)
            w1 = compete(spec, score, tiebreak=tb)
            # PSM-level TDC (context only)
            rep1 = rep[w1]
            q1 = tdc_qvalues(score[w1], dec[w1], block=rep1)
            psm_false = pair[w1] != true_pair[spec[w1]]
            # stage 2: best PSM per peptide (pair, side)
            unit = pair[w1] * 2 + dec[w1]
            w2 = w1[compete(unit, score[w1], tiebreak=tb[w1])]
            # stage 3: paired competition
            w3 = w2[compete(pair[w2], score[w2], tiebreak=tb[w2])]
            rep3 = rep[w3]
            q3 = tdc_qvalues(score[w3], dec[w3], block=rep3)
            pep_false = (pair[w3] - rep3 * n_pairs) >= n_present
            for al in alphas:
                f, nd = _fdp_from_q(q3, dec[w3], pep_false, rep3, r, al)
                out_pep[al][0].append(f); out_pep[al][1].append(nd); out_pep[al][2].append(np.full(r, seed))
                f, nd = _fdp_from_q(q1, dec[w1], psm_false, rep1, r, al)
                out_psm[al][0].append(f); out_psm[al][1].append(nd); out_psm[al][2].append(np.full(r, seed))
    fin = lambda o: {str(a): _summ(np.concatenate(v[0]), np.concatenate(v[1]), np.concatenate(v[2]), a)
                     for a, v in o.items()}
    return {"peptide_double_competition": fin(out_pep), "psm_level_not_fdr_controlled": fin(out_psm)}


# ------------------------------------------------------------------------ Group-walk

def sim_groupwalk(R: int, sizes=(1000, 1000, 1000), frac_true=(0.2, 0.4, 0.6), mus=(2.5, 3.0, 3.5),
                  alphas=ALPHAS, K: int = 40, chunk: int = 2000, sim_id: int = 3, factor: float = 1.0) -> dict:
    """Groups with different true fractions and alternative distributions
    (like the groupwalk R-package example). Also reports pooled TDC on the same
    data for the power comparison."""
    G = len(sizes)
    L = max(sizes)
    gw = {a: ([], [], []) for a in alphas}
    pooled = {a: ([], [], []) for a in alphas}
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            labs = np.zeros((r, G, L), bool)
            nul = np.zeros((r, G, L), bool)
            flat_s, flat_d, flat_n, flat_b = [], [], [], []
            for g, (m, ft, mu) in enumerate(zip(sizes, frac_true, mus)):
                nt = int(round(m * ft))
                t = rng.standard_normal((r, m))
                t[:, :nt] += mu
                d = rng.standard_normal((r, m))
                w = np.maximum(t, d)
                isd = d > t
                isnull = np.zeros((r, m), bool)
                isnull[:, nt:] = True
                o = np.argsort(w, axis=1)
                labs[:, g, :m] = np.take_along_axis(isd, o, 1)
                nul[:, g, :m] = np.take_along_axis(isnull, o, 1)
                flat_s.append(w); flat_d.append(isd); flat_n.append(isnull & ~isd)
            out = group_walk_batch(labs, nul, np.asarray(sizes), alphas, K=K, factor=factor)
            W = np.concatenate(flat_s, 1); Dm = np.concatenate(flat_d, 1); Nm = np.concatenate(flat_n, 1)
            blk = np.repeat(np.arange(r), W.shape[1])
            q = tdc_qvalues(W.ravel(), Dm.ravel(), block=blk, factor=factor)
            for a in alphas:
                T_a, F_a = out[a]
                gw[a][0].append(F_a / np.maximum(T_a, 1)); gw[a][1].append(T_a); gw[a][2].append(np.full(r, seed))
                f, nd = _fdp_from_q(q, Dm.ravel(), Nm.ravel(), blk, r, a)
                pooled[a][0].append(f); pooled[a][1].append(nd); pooled[a][2].append(np.full(r, seed))
    fin = lambda o: {str(a): _summ(np.concatenate(v[0]), np.concatenate(v[1]), np.concatenate(v[2]), a)
                     for a, v in o.items()}
    return {"group_walk": fin(gw), "pooled_tdc": fin(pooled)}


# ----------------------------------------------------------------------------- RESET

def _reset_features(rng, n, p_inf, mu_vec, null):
    X = rng.standard_normal((n, p_inf))
    X[~null] += mu_vec
    return X


def sim_reset(R: int, mode: str = "single", m: int = 5000, pi0: float = 0.5,
              mu=(3.0, 2.0, 1.5, 0.0, 0.0), alphas=ALPHAS, sim_id: int = 4, learner_factory=None,
              n_iter: int = 5, then_groupwalk: bool = False, n_groups: int = 2, naive: bool = False,
              n_noise: int = 0) -> dict:
    """Peptide pairs with p features; feature 0 is the search-engine score used
    for competition (and RESET's first iteration). Null targets and decoys share
    the feature distribution (exchangeable). ``mode='two'``: two decoys per
    target from independent draws. ``then_groupwalk``: exploratory RESET +
    Group-walk with decoy counts x c_final (unproven adaptation, DESIGN §6.4).
    ``naive``: negative control -- the same decoys train the learner and estimate
    the FDR (factor 1, no split), i.e. the Percolator-without-CV leak RESET fixes.
    ``n_noise``: extra pure-noise features (invites overfitting)."""
    mu = np.r_[np.asarray(mu, float), np.zeros(n_noise)]
    p = mu.size
    res = {a: ([], [], []) for a in alphas}
    not_started = 0
    for seed, Rs in _split(R):
        for i in range(Rs):
            rng = make_rng(seed, sim_id, i)
            n0 = int(round(m * pi0))
            null = np.zeros(m, bool)
            null[:n0] = True
            XT = _reset_features(rng, m, p, mu, null)
            XD = rng.standard_normal((m, p))
            if mode == "single":
                tb = rng.random(m) < 0.5
                dwin = (XD[:, 0] > XT[:, 0]) | ((XD[:, 0] == XT[:, 0]) & tb)
                X = np.where(dwin[:, None], XD, XT)
                isd = dwin
            else:
                XD2 = rng.standard_normal((m, p))
                S = np.column_stack([XT[:, 0], XD[:, 0], XD2[:, 0]])
                which = S.argmax(axis=1)
                X = np.where((which == 0)[:, None], XT, np.where((which == 1)[:, None], XD, XD2))
                isd = which > 0
            learner = learner_factory() if learner_factory is not None else None
            if naive:
                q = _naive_selftrained(X, isd, learner, n_iter)
                for a in alphas:
                    acc = (q <= a) & ~isd
                    nd = int(acc.sum()); nf = int((acc & null).sum())
                    res[a][0].append(nf / max(nd, 1)); res[a][1].append(nd); res[a][2].append(seed)
                continue
            rr = reset(X, isd, X[:, 0], mode=mode, seed=make_rng(seed, sim_id, i, 7), learner=learner,
                       n_iter=n_iter)
            not_started += int(not rr.n_positives or rr.n_positives[0] < 10)
            if then_groupwalk:
                keep = ~rr.is_training_decoy
                grp = (np.arange(m) % n_groups)[keep]   # label-blind arbitrary partition
                qg = np.full(m, np.nan)
                qg[keep] = group_walk(rr.score[keep], isd[keep], grp, factor=rr.c_final,
                                      seed=make_rng(seed, sim_id, i, 9), tie_break="random")
                q = qg
            else:
                q = rr.q
            for a in alphas:
                acc = (q <= a) & ~isd & ~rr.is_training_decoy
                nd = int(acc.sum())
                nf = int((acc & null).sum())
                res[a][0].append(nf / max(nd, 1)); res[a][1].append(nd); res[a][2].append(seed)
    out = {str(a): _summ(np.asarray(v[0]), np.asarray(v[1]), np.asarray(v[2]), a) for a, v in res.items()}
    if not naive:
        out["frac_training_not_started"] = not_started / max(R, 1)
    return out


def _naive_selftrained(X, isd, learner, n_iter, alpha_train=0.01):
    """Negative control: iterative training with ALL decoys as negatives and the
    same decoys used for estimation (factor 1). Not a valid procedure."""
    from .reset import LDALearner, _standardize

    Xs = _standardize(X)
    learner = LDALearner() if learner is None else learner
    score = X[:, 0].copy()
    for _ in range(n_iter):
        qt = tdc_qvalues(score, isd)
        pos = ~isd & (qt <= alpha_train)
        if pos.sum() < 10:
            break
        learner.fit(Xs[pos], Xs[isd])
        score = learner.decision_function(Xs)
    return tdc_qvalues(score, isd)


# ------------------------------------------------------------- picked protein FDR

def sim_picked_protein(R: int, n_prot: int = 2000, frac_present: float = 0.25, mean_pep: float = 4.0,
                       detect: float = 0.6, mu: float = 3.0, null_psm_rate: float = 0.4, alphas=ALPHAS,
                       chunk: int = 250, sim_id: int = 5) -> dict:
    """One-protein-per-group picked protein FDR (Savitski 2015) vs classic.

    Each target protein has 1+Poisson(mean_pep) peptides; its decoy mirrors the
    count. Peptides of present proteins are detected with probability ``detect``
    (score N(mu,1)); every other peptide (absent proteins, undetected peptides of
    present proteins, decoys) has a null PSM with probability ``null_psm_rate``
    (score N(0,1)), otherwise no score. Protein score = best peptide score.
    A target protein is false if absent.
    """
    out = {"picked": {a: ([], [], []) for a in alphas}, "classic": {a: ([], [], []) for a in alphas}}
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            npep = 1 + rng.poisson(mean_pep, size=(r, n_prot))
            tot = int(npep.sum())
            prot_of = np.repeat(np.arange(r * n_prot), npep.ravel())
            present = np.zeros((r, n_prot), bool)
            present[:, : int(round(n_prot * frac_present))] = True
            pres_pep = present.ravel()[prot_of]
            det = pres_pep & (rng.random(tot) < detect)
            tscore = np.where(det, rng.standard_normal(tot) + mu,
                              np.where(rng.random(tot) < null_psm_rate, rng.standard_normal(tot), -np.inf))
            dscore = np.where(rng.random(tot) < null_psm_rate, rng.standard_normal(tot), -np.inf)
            tmax = np.full(r * n_prot, -np.inf)
            dmax = np.full(r * n_prot, -np.inf)
            np.maximum.at(tmax, prot_of, tscore)
            np.maximum.at(dmax, prot_of, dscore)
            blk = np.repeat(np.arange(r), n_prot)
            false_t = ~present.ravel()
            for name, picked in (("picked", True), ("classic", False)):
                s, isd, q, pi, bb = picked_protein_fdr_simple(tmax, dmax, block=blk, picked=picked,
                                                              seed=make_rng(seed, sim_id, c0, 1))
                for a in alphas:
                    f, nd = _fdp_from_q(q, isd, false_t[pi], bb, r, a)
                    out[name][a][0].append(f); out[name][a][1].append(nd); out[name][a][2].append(np.full(r, seed))
    return {k: {str(a): _summ(np.concatenate(v[0]), np.concatenate(v[1]), np.concatenate(v[2]), a)
                for a, v in o.items()} for k, o in out.items()}


def _gen_proteome_with_isoforms(rng, n_genes, frac_present, mean_pep, detect, mu, null_psm_rate,
                                iso_frac=0.3):
    """Genes with optional isoform pairs sharing all but one peptide each."""
    pep2prot, pep_score, pair_of, is_dec, present = {}, {}, {}, {}, {}
    pid = 0
    for gi in range(n_genes):
        gp = rng.random() < frac_present
        n = 1 + rng.poisson(mean_pep)
        prots = [f"G{gi}"]
        if rng.random() < iso_frac:
            prots.append(f"G{gi}-2")
        shared = [f"p{pid + k}" for k in range(n)]
        pid += n
        uniq = {}
        if len(prots) == 2:
            for pr in prots:
                uniq[pr] = f"p{pid}"
                pid += 1
        # which isoform is expressed (if gene present): first one, the other only shares
        for side in (False, True):
            pre = "rev_" if side else ""
            for pr in prots:
                pair_of[pre + pr] = pr
                is_dec[pre + pr] = side
                present[pre + pr] = (not side) and gp and pr == prots[0]
            for pep in shared:
                pep2prot[pre + pep] = [pre + pr for pr in prots]
            for pr, up in uniq.items():
                pep2prot[pre + up] = [pre + pr]
        for side in (False, True):
            pre = "rev_" if side else ""
            for pep in shared + [uniq[pr] for pr in prots if pr in uniq]:
                owner_present = (not side) and gp and (pep in shared or pep == uniq.get(prots[0]))
                if owner_present and rng.random() < detect:
                    pep_score[pre + pep] = rng.standard_normal() + mu
                elif rng.random() < null_psm_rate:
                    pep_score[pre + pep] = rng.standard_normal()
    return pep2prot, pep_score, pair_of, is_dec, present


def sim_picked_group_general(R: int, n_genes: int = 600, frac_present: float = 0.3, mean_pep: float = 3.0,
                             detect: float = 0.6, mu: float = 3.0, null_psm_rate: float = 0.4,
                             alphas=ALPHAS, sim_id: int = 6, rescue: bool = True) -> dict:
    """Full picked-protein-group code path with isoform groups and shared peptides.
    A reported target group is false if none of its proteins is present."""
    out = {"picked_group": {a: ([], [], []) for a in alphas}, "classic": {a: ([], [], []) for a in alphas}}
    for seed, Rs in _split(R):
        for i in range(Rs):
            rng = make_rng(seed, sim_id, i)
            p2p, ps, po, isd, pres = _gen_proteome_with_isoforms(rng, n_genes, frac_present, mean_pep, detect, mu,
                                                                null_psm_rate)
            for name in ("picked_group", "classic"):
                if name == "picked_group":
                    res = picked_group_fdr(p2p, ps, po, isd, rescue=rescue, seed=make_rng(seed, sim_id, i, 1))
                else:
                    res = picked_group_fdr(p2p, ps, po, isd, rescue=False, picked=False,
                                           seed=make_rng(seed, sim_id, i, 1))
                for a in alphas:
                    acc = res.accepted(a)
                    nd = len(acc)
                    nf = sum(1 for g in acc if not any(pres[p] for p in g))
                    out[name][a][0].append(nf / max(nd, 1)); out[name][a][1].append(nd); out[name][a][2].append(seed)
    return {k: {str(a): _summ(np.asarray(v[0]), np.asarray(v[1]), np.asarray(v[2]), a) for a, v in o.items()}
            for k, o in out.items()}


# ----------------------------------------------------------------- multi-decoy TDC

def sim_multi_decoy(R: int, d: int = 3, m: int = 2000, pi0: float = 0.5, mu: float = 3.0, alphas=ALPHAS,
                    chunk: int = 500, sim_id: int = 7) -> dict:
    res = {a: ([], [], []) for a in alphas}
    kt = kd = (d + 1) // 2
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            n0 = int(round(m * pi0))
            t = rng.standard_normal((r, m))
            t[:, n0:] += mu
            D = rng.standard_normal((r * m, d))
            w, lab = multi_decoy_competition(t.ravel(), D, kt, kd, seed=make_rng(seed, sim_id, c0, 1))
            blk = np.repeat(np.arange(r), m)
            null = np.zeros((r, m), bool)
            null[:, :n0] = True
            keep = lab != 0
            q = np.ones(r * m)
            q[keep] = tdc_qvalues(w[keep], lab[keep] < 0, factor=kt / kd, block=blk[keep])
            for a in alphas:
                f, nd = _fdp_from_q(q, lab < 0, null.ravel(), blk, r, a)
                res[a][0].append(f); res[a][1].append(nd); res[a][2].append(np.full(r, seed))
    return {str(a): _summ(np.concatenate(v[0]), np.concatenate(v[1]), np.concatenate(v[2]), a)
            for a, v in res.items()}


# ------------------------------------------------------------ entrapment estimators

def sim_entrapment_estimators(R: int, m: int = 3000, frac_present: float = 0.4, mu: float = 3.0,
                              alphas=ALPHAS, chunk: int = 500, sim_id: int = 8) -> dict:
    """Peptide-level TDC on T u E u dec(T) u dec(E) with r = 1 paired shuffled
    entrapment. Compares the mean of each entrapment estimator with the mean true
    FDP of the reported list (false = absent target or entrapment)."""
    from .entrapment import combined, lower_bound, paired

    acc_res = {a: {"true": [], "lb": [], "combined": [], "paired": [], "nT": []} for a in alphas}
    for seed, Rs in _split(R):
        for c0 in range(0, Rs, chunk):
            r = min(chunk, Rs - c0)
            rng = make_rng(seed, sim_id, c0)
            npres = int(round(m * frac_present))
            t = rng.standard_normal((r, m))
            t[:, :npres] += mu
            e = rng.standard_normal((r, m))
            dt = rng.standard_normal((r, m))
            de = rng.standard_normal((r, m))
            # target-side vs own decoy (pair competition), for T and E separately
            wt = np.maximum(t, dt); t_dec = dt > t
            we = np.maximum(e, de); e_dec = de > e
            W = np.concatenate([wt, we], 1)
            Dm = np.concatenate([t_dec, e_dec], 1)
            isE = np.zeros((r, 2 * m), bool); isE[:, m:] = True
            absent = np.ones((r, m), bool); absent[:, :npres] = False
            false = np.concatenate([absent, np.ones((r, m), bool)], 1) & ~Dm
            blk = np.repeat(np.arange(r), 2 * m)
            q = tdc_qvalues(W.ravel(), Dm.ravel(), block=blk).reshape(r, 2 * m)
            for a in alphas:
                acc = (q <= a) & ~Dm
                nd = acc.sum(1)
                nf = (acc & false).sum(1)
                accT, accE = acc[:, :m], acc[:, m:]
                nT = accT.sum(1); nE = accE.sum(1)
                e_only = (accE & ~accT).sum(1)
                e_gt = (accE & accT & (we >= wt)).sum(1)
                acc_res[a]["true"].append(nf / np.maximum(nd, 1))
                acc_res[a]["lb"].append(lower_bound(nT, nE))
                acc_res[a]["combined"].append(combined(nT, nE, 1.0))
                acc_res[a]["paired"].append(paired(nT, nE, e_only, e_gt))
                acc_res[a]["nT"].append(nT)
    out = {}
    for a, d in acc_res.items():
        out[str(a)] = {k: float(np.concatenate(v).mean()) for k, v in d.items()}
        out[str(a)]["n_rep"] = int(sum(len(x) for x in d["true"]))
    return out


def timed(fn, *args, **kw):
    t0 = time.perf_counter()
    r = fn(*args, **kw)
    return r, time.perf_counter() - t0


# ------------------------------------------------- entrapment decision-rule demo

def sim_entrapment_verdict(alpha: float = 0.01, m: int = 30000, frac_present: float = 0.4, mu: float = 3.0,
                           p_target_vs_decoy: float = 0.5, p_entrap_vs_decoy: float | None = None,
                           entrap_seeds=None, decoy_seed: int = 1001, n_boot: int = 1000,
                           q_scale: float = 1.0) -> dict:
    """End-to-end use of the §6.3 decision rule on simulated augmented searches.

    One dataset of ``m`` original targets (``frac_present`` correct), paired
    shuffled entrapment (r = 1) and decoys for both. Null competitions are won
    by the target side with probability ``p_target_vs_decoy`` (absent targets)
    and ``p_entrap_vs_decoy`` (entrapments; default = same). Replicates =
    entrapment seeds (ENTRAP_SEEDS_DEV) x one decoy seed. Returns the paired
    arm's evaluation (mean FDP, CI, verdict) and the true mean FDP.
    ``q_scale`` < 1 mimics a tool with conservative q-values (the list at nominal
    alpha is the TDC list at q <= q_scale * alpha).
    """
    from ._util import ENTRAP_SEEDS_DEV
    from .entrapment import Replicate, evaluate_arm, paired_counts

    entrap_seeds = ENTRAP_SEEDS_DEV if entrap_seeds is None else entrap_seeds
    pe = p_target_vs_decoy if p_entrap_vs_decoy is None else p_entrap_vs_decoy
    base = make_rng(decoy_seed, 99)
    npres = int(round(m * frac_present))
    t_true = base.standard_normal(m) + mu
    reps, true_fdp = [], []
    for s in entrap_seeds:
        rng = make_rng(s, decoy_seed)
        # null competitions: winner score ~ max of two N(0,1); label by biased coin
        wt = np.where(np.arange(m) < npres, t_true, np.maximum(rng.standard_normal(m), rng.standard_normal(m)))
        t_dec = np.where(np.arange(m) < npres, rng.standard_normal(m) > t_true, rng.random(m) >= p_target_vs_decoy)
        we = np.maximum(rng.standard_normal(m), rng.standard_normal(m))
        e_dec = rng.random(m) >= pe
        W = np.r_[wt, we]
        D = np.r_[t_dec, e_dec]
        q = tdc_qvalues(W, D)
        acc = (q <= alpha * q_scale) & ~D
        pair = np.r_[np.arange(m), np.arange(m)]
        ent = np.r_[np.zeros(m, bool), np.ones(m, bool)]
        tside = ~D
        reps.append(Replicate(paired_counts(pair[tside], ent[tside], W[tside], acc[tside]), entrap_seed=s,
                              decoy_seed=decoy_seed))
        false = np.r_[np.arange(m) >= npres, np.ones(m, bool)]
        true_fdp.append((acc & false).sum() / max(acc.sum(), 1))
    ev = evaluate_arm(reps, alpha, "paired", n_boot=n_boot, seed=11)
    ev["true_mean_fdp"] = float(np.mean(true_fdp))
    ev["p_target_vs_decoy"] = p_target_vs_decoy
    ev["p_entrap_vs_decoy"] = pe
    ev["q_scale"] = q_scale
    return ev
