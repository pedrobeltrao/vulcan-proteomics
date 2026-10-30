"""T1 driver, Amendment A5 (DESIGN §6.0d, reviews/amendment_A5_draft.md A5.11.2.5): additive read-outs that the claim
script (bench/t1/claims_a5.py) needs.  Called after `a3.secondary`, from the same outputs; nothing is re-searched,
re-scored or re-weighted, and no number computed before A5 is changed (the module only reads the driver's outputs,
`metrics` and `evaluate.LAST_CURVES`; `evaluate.noninferiority` is called on a private copy so `metrics` keeps its
own `matched_fdp_curves`).

metrics["a5"] = {
  "version": "a5-1",
  "stratum1":   A5.5.5 (descriptive; NOT family S).  The v2_sep accepted units of MERGED stratum 1 (the label-blind
                merged block that contains chimeric rank 1, per replicate, as bench/t1/a3.py) at the common q_sep
                cutoff.  The replicate tables are restricted to stratum-1 pair winners (targets, entrapments and their
                decoys), so matched arm-P FDP uses that list's OWN paired counts (evaluate.paired_grid, first-crossing
                rule of evaluate.matched_ids).  Ratios + bootstrap = evaluate.noninferiority, LOO = a3.loo_ratios,
                against every other method of the run (every A3.3.3 comparator, plus v1).
                {"noninferiority": {comp: {level: {...}}}, "loo": {comp: {level: {...}}}, "ids": {level: {alpha: {
                "per_decoy_seed": [...], "mean": x}}}, "merged_block": {rep: label}}
  "armC_matched": A5.7.1 (a)-(f), for EVERY method (numerator or comparator; printed for both sides, the claim script
                decides which side may void).  For level in {peptidoform, sequence} and f in {1%, 5%}:
                (b) threshold = `ours_q_at_match`: first crossing of the 15-replicate mean arm-P FDP-hat above f on the
                    grid, interpolated (identical to evaluate.noninferiority's point estimate; cross-checked);
                    applied to each of the 5 arm-C replicates' OWN q;
                (c) LB = N_E / (N_T + N_E) per arm-C replicate at that threshold (combined counts, evaluate.
                    combined_official), mean and 95% CI under convention (a) (deciding) and (b) (reported);
                    r = metrics.methods[m].arm_C_r[level] is used only for the reported FDP-hat_C (`mean_fdp_c`).
                (d) the claim script voids a ratio iff lb_ci95[0] > f.
                {method: {level: {"0.01"|"0.05": {"ours_q_at_match", "fdp_reached", "lb_reps", "mean_lb", "lb_ci95",
                 "lb_ci95_b", "mean_fdp_c", "r"}}}}
  "rank_change": A1.8, label-free, on the 5 production replicates (mean over replicates; per replicate listed):
                v1 (`vulcan`): fraction of spectra whose best re-scored PSM has engine rank > 1, and fraction of
                accepted units (target, own q <= 1%, paired double competition exactly as evaluate.ours_units) whose
                winning PSM has engine rank > 1;  v2_sep: fraction of accepted units (target, q_sep <= 1%) whose
                winning PSM had engine rank > 1 within its group.
  "screened_armC": D18 secondary arm C, REPORTED ONLY (never changes a verdict).  Foreign (arm-C) peptides are
                excluded if (i) LADDER: they share a >= 4-ion terminal b or y ladder with any TARGET peptide of the
                search database whose precursor lies in the window -- implemented as an identical I/L-collapsed
                4-residue N-terminal prefix (b1-b4) or C-terminal suffix (y1-y4), with |M_E - M_T - k*1.003355| <=
                10 ppm for k in {-1, 0, 1} on unmodified + Met-ox (0..3) neutral masses (CAM); "any target" rather than
                "any human target" (on T1, HYE, this is stricter: yeast and E. coli targets screen too); or (ii) TAXON:
                they are a tryptic peptide (scripts/entrapment/peptides.digest, I/L-collapsed) of a foreign protein
                whose FASTA header matches the fixed plastid / photosynthesis / food-allergen list TAXON_PATTERNS
                below (3g remedy 1; written before any T1 output).  Arm C is then recomputed on the 5 arm-C
                replicates (combined estimator, r recomputed on the screened foreign set), at alpha in {0.1, 1, 5}%,
                peptidoform and sequence, conventions (a) and (b).
                {"n_excluded": int, "n_excluded_ladder": int, "n_excluded_taxon": int, "n_foreign": int,
                 "definition": str, "methods": {method: {level: {alpha: {"a": {...}, "b": {...}}}}}}
  "complete": bool, "missing": [...]}
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

import numpy as np

from . import a3 as A3
from . import evaluate as EV

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VERSION = "a5-1"
S1 = A3.SECONDARY + "_stratum1"
ARMC_LEVELS = ("peptidoform", "sequence")
ARMC_F = (0.01, 0.05)
SCREEN_ALPHAS = (0.001, 0.01, 0.05)
RANK_ALPHA = 0.01
LADDER_K = 4
LADDER_PPM = 10.0
ISO = 1.0033548378
ISO_OFFSETS = (-1, 0, 1)
OX = 15.9949146
MAX_OX = 3
# D18 / 3g remedy 1 taxon screen (fixed before any T1 output): plastid, photosynthesis and food-allergen proteins.
TAXON_PATTERNS = [r"chloroplast", r"plastid", r"photosystem", r"photosynth", r"ribulose", r"rubisco",
                  r"light-harvesting", r"chlorophyll", r"thylakoid", r"allergen",
                  r"GN=(rbcL|RBCS\S*|LHC[AB]\S*|TUFA|atp[AB]|psa\S+|psb\S+)\b"]
_TAXON = re.compile("|".join(TAXON_PATTERNS), re.IGNORECASE)
_ARM_KEYS = ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")


def _peptides_mod():
    p = os.path.join(ROOT, "scripts", "entrapment")
    if p not in sys.path:
        sys.path.insert(0, p)
    import peptides as PEP          # scripts/entrapment/peptides.py, unchanged
    return PEP


# ------------------------------------------------------------------------------------------------ helpers
def q_at_match(curve: dict, level: str, f: float):
    """Point-estimate matched threshold (A5.7.1 b): exactly evaluate.noninferiority's `matched` on all runs."""
    c = curve["cnt"][level]
    x = curve["ids"][level]
    fdp = EV.fdp_paired(c.sum(0)).mean(axis=(0, 1))
    _, q, reached = EV.matched_ids(fdp, x.sum(0).mean(0), f)
    return q, reached


def armc_lb(tabs_c: list, level: str, q: float, r: float, f: float) -> dict:
    """A5.7.1 (c): per arm-C replicate LB = N_E/(N_T+N_E) at threshold q (each replicate's own q); mean, CI (a)/(b)."""
    from vulcan.stats import entrapment as ENT
    reps = [ENT.Replicate(EV.combined_official(t, level, q), -1, j) for j, t in tabs_c]
    e = ENT.evaluate_arm(reps, f, "combined", r=r, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention="a")
    return {"lb_reps": [float(x) for x in e["lb_per_replicate"]], "mean_lb": float(e["mean_lb"]),
            "lb_ci95": [float(x) for x in e["lb_ci95_a"]], "lb_ci95_b": [float(x) for x in e["lb_ci95_b"]],
            "mean_fdp_c": float(e["mean_fdp"]), "r": float(r)}


def _level_tab(pf, level):
    return pf if level == "peptidoform" else EV.seq_units(pf)


# ------------------------------------------------------------------------------------------------ stratum 1
def _s1_view(tab, cc):
    lab = "_".join(map(str, cc["merge"][1]))
    merged = np.asarray(cc["merged"])
    assert len(merged) == len(tab)
    return tab[merged == lab].reset_index(drop=True), lab


def stratum1(work, runs, decoy_seeds, entrap_seeds, metrics, curves, methods, log=print):
    from . import ours_chimeric as OC
    from .dbs import parse_rep, replicate_names
    reps = replicate_names(decoy_seeds, entrap_seeds)
    R, J, S = len(runs), len(decoy_seeds), len(entrap_seeds)
    pm = EV.ProteinMap(os.path.join(work, "index", "full"))
    ids_rj = {lv: np.zeros((R, J, len(EV.GRID))) for lv in EV.LEVELS}
    cnt = {lv: np.zeros((R, J, S, len(EV.GRID), 4)) for lv in EV.LEVELS}
    ds_ids = {lv: {str(a): [] for a in EV.ALPHAS} for lv in EV.LEVELS}
    blocks = {}
    rank = {"per_replicate": {}}
    for rep in reps:
        kind, j, s = parse_rep(rep)
        if kind == "C":
            continue
        ji = list(decoy_seeds).index(j)
        rl = EV.RepLabels(os.path.join(work, "index", rep))
        tab, cc = OC.chim_units(os.path.join(work, "ours_chim", rep), rl, runs, "score_v2", return_detail=True,
                                qcol="q_sep")
        if kind == "prod":
            rank["per_replicate"][rep] = _rank_v2(os.path.join(work, "ours_chim", rep), runs, tab, cc)
        s1, lab = _s1_view(tab, cc)
        blocks[rep] = lab
        tabs = {"peptidoform": s1, "sequence": EV.seq_units(s1)}
        tabs["protein"] = EV.protein_table(tabs["sequence"], pm)
        if kind == "prod":
            for lv in EV.LEVELS:
                ids_rj[lv][:, ji] = EV.ids_per_run(tabs[lv], R, lv)
                for a in EV.ALPHAS:
                    ds_ids[lv][str(a)].append(EV.ids_dataset(tabs[lv], lv, a))
        else:
            si = list(entrap_seeds).index(s)
            for lv in EV.LEVELS:
                cnt[lv][:, ji, si] = EV.paired_grid(tabs[lv], lv, R)
        log(f"a5 stratum1 {rep}")
    c1 = {S1: {"ids": ids_rj, "cnt": cnt}}
    comps = [m for m in methods if m != A3.SECONDARY and m in curves]
    allc = {**{m: curves[m] for m in comps}, **c1}
    priv = {"methods": metrics["methods"]}              # private: noninferiority writes matched_fdp_curves into it
    ni = EV.noninferiority(priv, allc, [S1] + comps, S1)
    loo = A3.loo_ratios(allc, S1, comps)
    A3.check_loo_full_equals_noninferiority(loo, ni)
    ids = {lv: {a: {"per_decoy_seed": v, "mean": float(np.mean(v)) if v else None} for a, v in d.items()}
           for lv, d in ds_ids.items()}
    v = list(rank["per_replicate"].values())
    rank["frac_accepted_units_rank_gt1"] = float(np.mean([x["frac_accepted_units_rank_gt1"] for x in v])) if v else None
    return ({"noninferiority": ni, "loo": loo, "ids": ids, "merged_block": blocks,
             "matched_fdp_curve": priv["matched_fdp_curves"].get(S1)}, rank)


# ------------------------------------------------------------------------------------------------ A1.8 rank change
def _rank_v2(score_dir, runs, tab, cc):
    import pandas as pd
    import pyarrow.parquet as pq
    parts = [pq.read_table(os.path.join(score_dir, f"{r}.chim.scored.parquet"), columns=["rank"]).to_pandas()
             for r in runs]
    rk = pd.concat(parts, ignore_index=True)["rank"].values[np.asarray(cc["win"])]
    acc = (~tab.is_decoy.values) & (tab.q.values <= RANK_ALPHA)
    return {"n_accepted": int(acc.sum()), "frac_accepted_units_rank_gt1": float((rk[acc] > 1).mean()) if acc.any() else None}


def rank_v1(score_dir, rl, runs) -> dict:
    """A1.8 for v1 on one production replicate (label-free apart from the production target/decoy competition)."""
    import pandas as pd
    import pyarrow.parquet as pq

    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    parts = [pq.read_table(os.path.join(score_dir, f"{r}.scored.parquet")).to_pandas() for r in runs]
    psm = pd.concat(parts, ignore_index=True)
    # best re-scored PSM per spectrum (stable order on ties: first row, i.e. the engine order)
    o = np.lexsort((np.arange(len(psm)), -psm.score.values, psm.spec_idx.values, psm.run.astype(str).values))
    p = psm.iloc[o]
    first = ~p.duplicated(["run", "spec_idx"]).values
    best_rank = p["rank"].values[first]
    # accepted units: exactly evaluate.ours_units' competition
    row = psm.pf_row.values
    dec = rl.is_decoy_row[row]
    unit = psm.pf_id.values
    spec = (psm.run.astype(str) + ":" + psm.spec_idx.astype(str)).values
    res = double_competition(spec, unit, rl.tdc_pair[row], psm.score.values, dec, seed=EV.SEED,
                             tiebreak=hash_tiebreak(unit, EV.SEED))
    acc = (~np.asarray(res.is_decoy)) & (np.asarray(res.q) <= RANK_ALPHA)
    wr = psm["rank"].values[np.asarray(res.row)]
    return {"n_spectra": int(first.sum()), "frac_spectra_rank_gt1": float((best_rank > 1).mean()),
            "n_accepted": int(acc.sum()),
            "frac_accepted_units_rank_gt1": float((wr[acc] > 1).mean()) if acc.any() else None}


# ------------------------------------------------------------------------------------------------ D18 screen
def _variant_masses(seqs, PEP):
    ms, idx = [], []
    for i, s in enumerate(seqs):
        m0 = PEP.mass(s, cam=True)
        for k in range(min(s.count("M"), MAX_OX) + 1):
            ms.append(m0 + k * OX)
            idx.append(i)
    return np.asarray(ms), np.asarray(idx)


def ladder_screen(t_seqs, f_seqs, PEP=None) -> np.ndarray:
    """bool[len(f_seqs)]: F shares an I/L-collapsed 4-residue terminal prefix or suffix with a target whose precursor
    mass is within LADDER_PPM (isotope offsets ISO_OFFSETS, Met-ox 0..MAX_OX on both sides)."""
    import pandas as pd
    PEP = PEP or _peptides_mod()
    t_seqs = [PEP.il(s) for s in t_seqs]
    f_seqs = [PEP.il(s) for s in f_seqs]
    tm, ti = _variant_masses(t_seqs, PEP)
    fm, fi = _variant_masses(f_seqs, PEP)
    hit = np.zeros(len(f_seqs), bool)
    for end in ("pre", "suf"):
        cut = (lambda s: s[:LADDER_K]) if end == "pre" else (lambda s: s[-LADDER_K:])
        tk = np.array([cut(s) for s in t_seqs], dtype=object)[ti]
        fk = np.array([cut(s) for s in f_seqs], dtype=object)[fi]
        codes, uniq = pd.factorize(np.concatenate([tk, fk]))
        tc, fc = codes[:len(tk)].astype(np.float64), codes[len(tk):].astype(np.float64)
        span = 1.0e4                                         # masses < 1e4 Da: key*1e4 + mass is unique per key
        tkey = np.sort(tc * span + tm)
        for k in ISO_OFFSETS:
            mf = fm - k * ISO
            tol = mf * LADDER_PPM * 1e-6
            lo = np.searchsorted(tkey, fc * span + mf - tol, side="left")
            hi = np.searchsorted(tkey, fc * span + mf + tol, side="right")
            h = hi > lo
            hit[fi[h]] = True
    return hit


def taxon_peptides(work, PEP=None) -> set:
    """I/L-collapsed tryptic peptides of foreign proteins whose header matches TAXON_PATTERNS."""
    PEP = PEP or _peptides_mod()
    meta = json.load(open(os.path.join(work, "armC", "ds.C.json")))
    out = set()
    for fa in meta.get("foreign_fastas", []):
        p = fa if os.path.isabs(fa) else os.path.join(ROOT, fa)
        for h, s in PEP.read_fasta(p):
            if _TAXON.search(h):
                out.update(PEP.il(x) for x in PEP.digest(s))
    return out


def screen_foreign(work, rl, log=print) -> tuple[set, dict]:
    PEP = _peptides_mod()
    t = [s for s, (fam, dec, _) in rl.seq.items() if fam == 0 and not dec]
    f = [s for s, (fam, dec, _) in rl.seq.items() if fam == -1 and not dec]
    t0 = time.perf_counter()
    lad = ladder_screen(t, f, PEP)
    tax_set = taxon_peptides(work, PEP)
    tax = np.array([PEP.il(s) in tax_set for s in f], bool)
    excl = {s for s, a, b in zip(f, lad, tax) if a or b}
    info = {"n_foreign": len(f), "n_excluded": len(excl), "n_excluded_ladder": int(lad.sum()),
            "n_excluded_taxon": int(tax.sum()), "screen_s": time.perf_counter() - t0}
    log(f"a5 D18 screen: {info}")
    return excl, info


def _screened_tab(pf, excl):
    drop = (pf.fam.values == -1) & pf.seq.isin(excl).values
    return pf[~drop].reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ main
def additions(work, runs, decoy_seeds, entrap_seeds, metrics, curves, methods, meta, log=print, loader=None) -> dict:
    """See module docstring.  `loader(method, rep, RepLabels) -> peptidoform DataFrame | None` is the driver's
    evaluate loader (the same function evaluate used); it is needed for the arm-C replicates."""
    from vulcan.stats import entrapment as ENT

    from .dbs import replicate_names
    t0 = time.perf_counter()
    out = {"version": VERSION, "missing": []}
    if A3.SECONDARY in methods and A3.SECONDARY in curves:
        out["stratum1"], rank2 = stratum1(work, runs, decoy_seeds, entrap_seeds, metrics, curves, methods, log)
    else:
        out["stratum1"], rank2 = None, None
        out["missing"].append("stratum1: v2_sep curves absent")
    # A1.8 rank change, v1
    rank = {}
    if "vulcan" in methods:
        per = {}
        for rep in [r for r in replicate_names(decoy_seeds, entrap_seeds) if r.startswith("prod_")]:
            try:
                per[rep] = rank_v1(os.path.join(work, "ours", rep),
                                   EV.RepLabels(os.path.join(work, "index", rep)), runs)
            except (FileNotFoundError, KeyError) as e:
                out["missing"].append(f"rank_change vulcan {rep}: {type(e).__name__}")
        if per:
            rank["vulcan"] = {"frac_spectra_rank_gt1": float(np.mean([v["frac_spectra_rank_gt1"] for v in per.values()])),
                              "frac_accepted_units_rank_gt1":
                                  float(np.mean([v["frac_accepted_units_rank_gt1"] for v in per.values()
                                                 if v["frac_accepted_units_rank_gt1"] is not None] or [np.nan])),
                              "alpha": RANK_ALPHA, "per_replicate": per}
    if rank2 is not None:
        rank[A3.SECONDARY] = {"frac_accepted_units_rank_gt1": rank2["frac_accepted_units_rank_gt1"],
                              "alpha": RANK_ALPHA, "per_replicate": rank2["per_replicate"]}
    out["rank_change"] = rank
    # arm-C replicates: matched-threshold LB (A5.7.1) and the D18 screened arm C, one load per method x replicate
    creps = [r for r in replicate_names(decoy_seeds, entrap_seeds) if r.startswith("C_")]
    rls = {r: EV.RepLabels(os.path.join(work, "index", r)) for r in creps}
    excl, sinfo = screen_foreign(work, rls[creps[0]], log)
    rl0 = rls[creps[0]]
    n_t = {"peptidoform": rl0.n_fwd_pf[0], "sequence": rl0.n_fwd_seq[0]}
    f_kept_pf = int(sum(1 for s, (fam, dec, _) in rl0.seq.items() if fam == -1 and not dec and s not in excl))
    armc, scr = {}, {}
    if loader is None:
        out["missing"].append("armC: no loader")
    for m in (methods if loader is not None else []):
        if m not in curves or "arm_C_r" not in metrics["methods"].get(m, {}):
            out["missing"].append(f"armC {m}: method not evaluated")
            continue
        tabs = []
        for rep in creps:
            pf = loader(m, rep, rls[rep])
            if pf is None:
                out["missing"].append(f"armC {m} {rep}")
                break
            tabs.append((int(rep.split("_")[1]), pf))
        if len(tabs) != len(creps):
            continue
        rC = metrics["methods"][m]["arm_C_r"]
        armc[m] = {}
        for lv in ARMC_LEVELS:
            lt = [(j, _level_tab(pf, lv)) for j, pf in tabs]
            armc[m][lv] = {}
            for f in ARMC_F:
                q, reached = q_at_match(curves[m], lv, f)
                armc[m][lv][str(f)] = {"ours_q_at_match": float(q), "fdp_reached": bool(reached),
                                       **armc_lb(lt, lv, q, rC[lv], f)}
        # D18 screened arm C (reported only)
        scr[m] = {}
        st = [(j, _screened_tab(pf, excl)) for j, pf in tabs]
        # r on the screened foreign set: sequence-level count from the same rule as RepLabels.n_fwd_seq
        f_kept = {"peptidoform": _n_fwd_pf_kept(rl0, excl), "sequence": f_kept_pf}
        for lv in ARMC_LEVELS:
            lt = [(j, _level_tab(pf, lv)) for j, pf in st]
            r_s = f_kept[lv] / n_t[lv]
            scr[m][lv] = {}
            for al in SCREEN_ALPHAS:
                reps_ = [ENT.Replicate(EV.combined_official(t, lv, al), -1, j) for j, t in lt]
                d = {}
                for cv in ENT.CI_CONVENTIONS:
                    e = ENT.evaluate_arm(reps_, al, "combined", r=r_s, n_boot=EV.N_BOOT, seed=EV.SEED, ci_convention=cv)
                    d[cv] = {k: e[k] for k in _ARM_KEYS}
                d["r_screened"] = r_s
                scr[m][lv][str(al)] = d
        log(f"a5 armC {m}")
    out["armC_matched"] = armc
    out["screened_armC"] = {**sinfo, "definition": __doc__.split('"screened_armC":')[1].split('"complete"')[0].strip(),
                            "reported_only": True, "methods": scr}
    # cross-check (b): v1's matched threshold equals evaluate.noninferiority's ours_q_at_match
    ni = metrics.get("noninferiority", {})
    for comp, byl in ni.items():
        for lv in ARMC_LEVELS:
            for f in ARMC_F:
                ref = byl.get(lv, {}).get("matched_fdp", {}).get(str(f), {}).get("ours_q_at_match")
                got = armc.get("vulcan", {}).get(lv, {}).get(str(f), {}).get("ours_q_at_match")
                if ref is not None and got is not None and not np.isclose(ref, got, rtol=1e-12, atol=0):
                    raise AssertionError(f"a5 matched threshold != evaluate ({lv}, {f}): {got} vs {ref}")
        break
    want = [m for m in A3.A3_METHODS if m in methods]
    out["complete"] = bool(not out["missing"] and out["stratum1"] is not None and set(want) <= set(armc)
                           and set(want) <= set(scr) and "vulcan" in rank and A3.SECONDARY in rank)
    out["wall_s"] = time.perf_counter() - t0
    return out


def _n_fwd_pf_kept(rl, excl) -> int:
    """Forward foreign peptidoform rows (as RepLabels.n_fwd_pf[-1]) whose sequence survives the screen."""
    fwd = (~rl.is_decoy_row) & (rl.fwdfam_row == -1)
    rows = np.flatnonzero(fwd)
    if not excl:
        return int(len(rows))
    seqs = rl.row_seq(rows)
    return int(sum(1 for s in seqs if s not in excl))
