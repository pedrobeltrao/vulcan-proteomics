"""T1 driver, stage `evaluate`: the pre-registered T1 comparison (DESIGN §1 "no loss", §6.3, §6.0 A1.6/A1.9).

UNITS (identical for every method; every method's output is mapped onto the replicate DB's sequences):
  peptidoform  (I/L-collapsed sequence, Met-ox position set).  The method's own peptide-level q (Sage peptide_q; mokapot
               peptide q-value; Vulcan paired double competition, peptidoform units, hash tie-break 11).
  sequence     I/L-collapsed stripped sequence; q = min over its peptidoforms, score = max (same rule for all).
  protein      picked protein groups (vulcan.stats.protein.picked_group_fdr, subset grouping, rescue, seed 11) over the
               method's sequence-level best scores, proteins from the SAME FASTA mapping for every method; target
               groups at picked q <= alpha.  (Our estimator for every method: Sage's own protein_q is not used.)
  A unit is "identified in run r" if it is accepted (dataset-level q) and has >= 1 top-1 PSM in run r.

ENTRAPMENT (augmented searches; sensitivity from production searches, §6.3 "Reporting"):
  arm P  paired, replicates = 5 decoy seeds x 3 entrapment seeds; pair key = target sequence (sequence level),
         (target sequence, number of Met-ox) (peptidoform level), target protein (protein level, FDRBench rule).
  arm C  combined (natural foreign proteome), 5 replicates (decoy seeds); r = |F| / |T| forward units on the
         effective search space (index), at the level's unit.
  Verdicts: vulcan.stats.entrapment.evaluate_arm / combine_arms (static_low_capacity: P and C), under BOTH CI
  conventions of A1.6 ((a) bootstrap + between-seed, (b) bootstrap only); 3e's coverage simulation chose (a), which
  decides; the more conservative of the two is reported alongside.  Arm M (A1.5) is not part of T1 (T5 / dev D7).

NON-INFERIORITY (§1): ratio R = ours / comparator of run-level IDs summed over runs and averaged over decoy shuffles;
  95% CI = percentile bootstrap (2000, seed 11) resampling runs x decoy shuffles (both with replacement, the same draw
  for both methods).  Pass iff lower CI >= 0.97.
  (i) at each method's own q <= alpha  -- counts only if the comparator's verdict at alpha is "supported";
  (ii) at matched entrapment FDP f in {1%, 5%}: each method's threshold (its own q, on a fixed grid) is moved until
      its arm-P FDP estimate (mean over the 15 replicates of run-level-summed paired counts) reaches f (first
      crossing, linear interpolation between grid points); IDs at that threshold from the production searches.
      The FDP curve is re-estimated inside every bootstrap draw.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import pyarrow as pa

from vulcan.stats import entrapment as ENT
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.protein import picked_group_fdr
from vulcan.stats.tdc import double_competition

from . import comparators as CMP
from .dbs import parse_rep, replicate_names, seq_bytes

ALPHAS = (0.001, 0.005, 0.01, 0.05)
FMATCH = (0.01, 0.05)
if os.environ.get("VULCAN_T1_FMATCH_0001") == "1":     # E7b (dev): also report matched arm-P FDP 0.1% (opt-in)
    FMATCH = FMATCH + (0.001,)
NI_MARGIN = 0.97
N_BOOT = 2000
SEED = 11
GRID = np.unique(np.concatenate([np.geomspace(1e-4, 0.25, 160), ALPHAS]))
LEVELS = ("peptidoform", "sequence", "protein")
CONSERVATIVE = ["violated", "not-evaluable", "inconclusive", "supported-noisy", "supported"]   # worst first
# A1.6: the CI convention is chosen by Monte-Carlo coverage (scripts/m3e/ci_convention_sim.py,
# results/2026-10-27_M3e_hardening/ci_convention_sim.json: chosen = "a", coverage 0.94-0.95 vs 0.83 for "b").
# Both conventions are still reported; "deciding" = convention (a); "deciding_conservative" = worse of (a), (b).
CI_CONVENTION_DECIDING = "a"
# A3 (bench/t1/a3.py): the per-method curves of the last evaluate() call (run x decoy-seed IDs and arm-P counts on the
# grid), exposed read-only for the A3.3.4 leave-one-run-out backstop.  Nothing else reads or changes them.
LAST_CURVES: dict = {}


def _read(path):
    with pa.OSFile(path, "rb") as src:
        return pa.ipc.open_file(src).read_all()


# ------------------------------------------------------------------------------------------- replicate labels
class RepLabels:
    """sequence -> (forward family, is_decoy, target_seq, n of forward rows) for one replicate index (L4 only)."""

    def __init__(self, index_dir):
        pf = _read(os.path.join(index_dir, "peptidoforms.arrow"))
        lab = _read(os.path.join(index_dir, "labels.arrow"))
        res = np.load(os.path.join(index_dir, "residues.npy"), mmap_mode="r")
        off = pf.column("seq_off").to_numpy()
        ln = pf.column("length").to_numpy()
        self.mask = pf.column("mask").to_numpy()
        self.off = off
        self.ln = ln
        self.res = res
        fam = lab.column("family").to_numpy()
        es = lab.column("entrap_seed").to_numpy()
        self.is_decoy_row = lab.column("is_decoy").to_numpy(zero_copy_only=False)
        self.fwdfam_row = np.where(fam % 2 == 0, 0, np.where(es == -1, -1, 1))       # 0 T, 1 E, -1 F
        self.tdc_pair = lab.column("tdc_pair").to_numpy()
        tseq = lab.column("target_seq")
        self.sb = seq_bytes(res, off, ln)
        u, first = np.unique(self.sb, return_index=True)
        ts = tseq.take(pa.array(first)).to_pylist()
        self.seq = {b.decode(): (int(self.fwdfam_row[i]), bool(self.is_decoy_row[i]), t)
                    for b, i, t in zip(u, first, ts)}
        fwd = ~self.is_decoy_row
        self.n_fwd_pf = {k: int(((self.fwdfam_row == k) & fwd).sum()) for k in (0, 1, -1)}
        ff = np.array([v[0] for v in self.seq.values()])
        dd = np.array([v[1] for v in self.seq.values()])
        self.n_fwd_seq = {k: int(((ff == k) & ~dd).sum()) for k in (0, 1, -1)}

    def row_seq(self, rows):
        return [b.decode() for b in self.sb[np.asarray(rows)]]


# ------------------------------------------------------------------------------------------- unit tables
ORIENTATION_FLIPS = []


def _orient(df: pd.DataFrame) -> pd.DataFrame:
    """The estimators compare scores (paired E>T, picked protein); every method's q must be non-increasing in its
    score. mokapot 0.10.0 sometimes writes the 'mokapot score' column with the opposite sign while its q-values are
    correct (seen on one D1 decoy seed, and as the M1e "collapsed protein count"). Orientation is decided from the
    method's own q-values only (label-blind): if score and q are positively rank-correlated, the score is negated."""
    if len(df) > 10:
        rq = df.q.rank().values
        rs = df.score.rank().values
        if np.corrcoef(rq, rs)[0, 1] > 0:
            df["score"] = -df["score"]
            ORIENTATION_FLIPS.append(int(len(df)))
    return df


def _finish(df: pd.DataFrame, rl: RepLabels, runs: list[str]) -> pd.DataFrame:
    """df: seq, mask, is_decoy, q, score, runbits (one row per peptidoform) -> + fam, target_seq, nvar."""
    df = _orient(df)
    info = [rl.seq.get(s) for s in df.seq]
    miss = sum(x is None for x in info)
    if miss:
        raise RuntimeError(f"{miss} peptides not in the replicate DB")
    df["fam"] = np.array([x[0] for x in info], np.int8)
    dec_db = np.array([x[1] for x in info])
    if not np.array_equal(dec_db, df.is_decoy.values):
        raise RuntimeError(f"decoy label mismatch on {(dec_db != df.is_decoy.values).sum()} peptides")
    df["target_seq"] = [x[2] for x in info]
    df["nvar"] = [bin(int(m)).count("1") for m in df["mask"]]
    return df


def ours_units(score_dir, rl: RepLabels, runs) -> pd.DataFrame:
    import pyarrow.parquet as pq
    parts = [pq.read_table(os.path.join(score_dir, f"{r}.scored.parquet")).to_pandas() for r in runs]
    psm = pd.concat(parts, ignore_index=True)
    row = psm.pf_row.values
    dec = rl.is_decoy_row[row]
    unit = psm.pf_id.values
    spec = (psm.run.astype(str) + ":" + psm.spec_idx.astype(str)).values
    res = double_competition(spec, unit, rl.tdc_pair[row], psm.score.values, dec, seed=SEED,
                             tiebreak=hash_tiebreak(unit, SEED))
    rix = {r: i for i, r in enumerate(runs)}
    w = res.psm_winners
    wb = pd.DataFrame({"unit": unit[w], "bit": [1 << rix[r] for r in psm.run.values[w]]})
    bits = wb.groupby("unit").bit.agg(lambda x: int(np.bitwise_or.reduce(x.values)))
    rows = row[res.row]
    df = pd.DataFrame({"seq": rl.row_seq(rows), "mask": rl.mask[rows].astype(np.uint64), "is_decoy": res.is_decoy,
                       "q": res.q, "score": res.score})
    df["runbits"] = bits.reindex(unit[res.row]).fillna(0).astype(np.int64).values
    return _finish(df, rl, runs)


def _agg(rows, runs, run_of):
    """rows: iterable of (peptide string, is_decoy, q, score, run) -> peptidoform DataFrame."""
    rix = {r: i for i, r in enumerate(runs)}
    d = {}
    for pep, dec, q, s, run in rows:
        k = CMP.parse_peptide(pep) + (bool(dec),)
        b = 1 << rix[run_of(run)] if run is not None else 0
        cur = d.get(k)
        if cur is None:
            d[k] = [q, s, b]
        else:
            cur[0] = min(cur[0], q)
            cur[1] = max(cur[1], s)
            cur[2] |= b
    df = pd.DataFrame([(k[0], k[1], k[2], v[0], v[1], v[2]) for k, v in d.items()],
                      columns=["seq", "mask", "is_decoy", "q", "score", "runbits"])
    df["mask"] = df["mask"].astype(np.uint64)
    return df


def sage_units(sage_dir, rl, runs, all_ranks=False):
    """all_ranks (E7b, Sage chimera mode): keep every rank (Sage's chimeric PSMs); default rank 1 only (T1)."""
    t = pd.read_csv(os.path.join(sage_dir, "results.sage.tsv"), sep="\t",
                    usecols=["peptide", "label", "peptide_q", "sage_discriminant_score", "filename", "rank"])
    if not all_ranks:
        t = t[t["rank"] == 1]
    rows = zip(t.peptide.values, (t.label.values == -1), t.peptide_q.values, t.sage_discriminant_score.values,
               t.filename.values)
    return _finish(_agg(rows, runs, lambda f: f[:-5] if f.endswith(".mzML") else f), rl, runs)


def mokapot_units(mk_dir, rl, runs, tool):
    """Peptide-level q/score from mokapot.{,decoy.}peptides.txt; run membership from mokapot.{,decoy.}psms.txt."""
    if tool == "comet":
        run_of = lambda spec: os.path.basename(spec).rsplit("_", 3)[0]          # noqa: E731  <path>/<run>_<scan>_<z>_<rank>
    else:
        run_of = lambda f: f[:-5] if f.endswith(".mzML") else f                  # noqa: E731
    pep = {}
    for fn, dec in (("mokapot.peptides.txt", False), ("mokapot.decoy.peptides.txt", True)):
        p = os.path.join(mk_dir, fn)
        if not os.path.exists(p):
            continue
        for x in CMP.read_ragged(p):
            pep[(x["Peptide"], dec)] = (float(x["mokapot q-value"]), float(x["mokapot score"]))
    rows = []
    for fn, dec in (("mokapot.psms.txt", False), ("mokapot.decoy.psms.txt", True)):
        p = os.path.join(mk_dir, fn)
        if not os.path.exists(p):
            continue
        for x in CMP.read_ragged(p):
            key = (x["Peptide"], dec)
            if key not in pep:
                continue
            run = x["SpecId"] if tool == "comet" else x["FileName"]
            rows.append((x["Peptide"], dec, pep[key][0], pep[key][1], run))
    return _finish(_agg(rows, runs, run_of), rl, runs)


def seq_units(pf: pd.DataFrame) -> pd.DataFrame:
    g = pf.groupby(["seq", "is_decoy"], sort=False)
    out = g.agg(q=("q", "min"), score=("score", "max"), runbits=("runbits", lambda x: int(np.bitwise_or.reduce(x.values))),
                fam=("fam", "first"), target_seq=("target_seq", "first")).reset_index()
    return out


# ------------------------------------------------------------------------------------------- proteins
class ProteinMap:
    def __init__(self, index_dir):
        t = _read(os.path.join(index_dir, "proteins.arrow"))
        pep = t.column("peptide").to_pylist()
        self.off = t.column("prot_off").to_numpy()
        self.ln = t.column("prot_len").to_numpy()
        self.pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
        self.acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
        self.pos = dict(zip(pep, range(len(pep))))

    def prots(self, tseq):
        i = self.pos.get(tseq)
        return [] if i is None else [self.acc[j] for j in self.pidx[self.off[i]:self.off[i] + self.ln[i]]]


def protein_table(sq: pd.DataFrame, pm: ProteinMap):
    """Picked groups over sequence-level scores. -> DataFrame(group, score, q, is_decoy, is_obsolete, is_target_grp,
    rep_target (representative original target protein), e_partner (paired target protein of an entrapment group),
    runbits)."""
    pre = {(0, False): "", (1, False): "E_", (0, True): "rev_", (1, True): "rev_E_"}
    pep2prot, pep_score, bits = {}, {}, {}
    for s, dec, sc, fam, ts, rb in zip(sq.seq.values, sq.is_decoy.values, sq.score.values, sq.fam.values,
                                       sq.target_seq.values, sq.runbits.values):
        if fam == -1:            # natural foreign peptides: no protein (protein entrapment is arm P only, §6.3)
            continue
        pr = pm.prots(ts)
        if not pr:
            continue
        k = ("D:" if dec else "T:") + s
        pp = [pre[(int(fam), bool(dec))] + p for p in pr]
        pep2prot[k] = pp
        pep_score[k] = float(sc)
        for p in pp:
            bits[p] = bits.get(p, 0) | int(rb)
    allp = {p for v in pep2prot.values() for p in v}
    base = lambda p: p[4:] if p.startswith("rev_") else p            # noqa: E731
    res = picked_group_fdr(pep2prot, pep_score, {p: base(p) for p in allp}, {p: p.startswith("rev_") for p in allp},
                           seed=SEED)
    is_e = lambda p: p.startswith("E_")                               # noqa: E731
    rows = []
    for g, s, q, d, o in zip(res.groups, res.score, res.q, res.is_decoy, res.is_obsolete):
        tgt = any(not is_e(p) for p in g)
        rep = next((p for p in g if not is_e(p)), None)
        rb = 0
        for p in g:
            rb |= bits.get(p, 0)
        rows.append((tuple(g), float(s), float(q), bool(d), bool(o), tgt, rep, (g[0][2:] if not tgt else None), rb))
    return pd.DataFrame(rows, columns=["group", "score", "q", "is_decoy", "is_obsolete", "is_target_grp",
                                       "rep_target", "e_partner", "runbits"])


# ------------------------------------------------------------------------------------------- counts
def ids_per_run(tab: pd.DataFrame, nruns: int, level: str) -> np.ndarray:
    """[run, grid] run-level IDs (targets of the original target family only)."""
    if level == "protein":
        t = tab[~tab.is_decoy & ~tab.is_obsolete & tab.is_target_grp]
    else:
        t = tab[~tab.is_decoy & (tab.fam == 0)]
    q = t.q.values
    rb = t.runbits.values.astype(np.int64)
    out = np.zeros((nruns, len(GRID)))
    for r in range(nruns):
        qr = np.sort(q[(rb >> r) & 1 == 1])
        out[r] = np.searchsorted(qr, GRID, side="right")
    return out


def ids_dataset(tab, level, alpha) -> int:
    if level == "protein":
        return int(((tab.q <= alpha) & ~tab.is_decoy & ~tab.is_obsolete & tab.is_target_grp).sum())
    return int(((tab.q <= alpha) & ~tab.is_decoy & (tab.fam == 0)).sum())


def _pair_key(tab, level):
    if level == "peptidoform":
        return tab.target_seq.astype(str) + "|" + tab.nvar.astype(str)
    return tab.target_seq.astype(str)


def paired_official(tab, level, alpha):
    """vulcan.stats.entrapment counts at the dataset level (the estimator of record)."""
    if level == "protein":
        t = tab[~tab.is_decoy & ~tab.is_obsolete]
        is_e = {p: p.startswith("E_") for g in t.group for p in g}
        partner = {p: p[2:] for p in is_e if is_e[p]}
        return ENT.protein_paired_counts(list(t.group), t.score.values, (t.q <= alpha).values, is_e, partner)
    t = tab[~tab.is_decoy & (tab.fam != -1)]
    return ENT.paired_counts(_pair_key(t, level).values, (t.fam == 1).values, t.score.values, (t.q <= alpha).values)


def combined_official(tab, level, alpha):
    t = tab[~tab.is_decoy & (tab.fam != 1)]
    key = (t.seq.astype(str) + "|" + t["mask"].astype(str)) if level == "peptidoform" else t.seq
    return ENT.combined_counts(key.values, (t.fam == -1).values, (t.q <= alpha).values)


def paired_grid(tab, level, nruns, check=True) -> np.ndarray:
    """[run, grid, 4] totals (n_t, n_e, n_e_only, n_e_gt_t) of the paired estimator on run-level lists; vectorized
    (q is non-increasing in score for every method, so a side's best accepted score is its overall best); checked
    against the official counts on the pooled list at alpha = 1% and 5%."""
    if level == "protein":
        t = tab[~tab.is_decoy & ~tab.is_obsolete].copy()
        t["side_e"] = ~t.is_target_grp
        t["key"] = np.where(t.side_e, t.e_partner.astype(object), t.rep_target.astype(object))
    else:
        t = tab[~tab.is_decoy & (tab.fam != -1)].copy()
        t["side_e"] = t.fam == 1
        t["key"] = _pair_key(t, level)
    out = np.zeros((nruns + 1, len(GRID), 4))
    allbits = (1 << nruns) - 1
    for r in range(nruns + 1):
        m = ((t.runbits.values.astype(np.int64) >> r) & 1 == 1) if r < nruns else (t.runbits.values >= 0)
        s = t[m]
        if level == "protein":
            # targets: distinct representative proteins; entrapment groups compared with the best target group
            # containing their partner (any member)
            tg = s[~s.side_e]
            rep_q = tg.groupby("rep_target").q.min()
            memb = [(p, q, sc) for g, q, sc in zip(tg.group, tg.q, tg.score) for p in g if not p.startswith("E_")]
            mt = pd.DataFrame(memb, columns=["p", "q", "s"]) if memb else pd.DataFrame(columns=["p", "q", "s"])
            tbest = mt.groupby("p").agg(q=("q", "min"), s=("s", "max"))
            eg = s[s.side_e]
            n_t = (rep_q.values[:, None] <= GRID[None]).sum(0)
            eq = eg.q.values[:, None] <= GRID[None]
            tq = tbest.q.reindex(eg.e_partner.values).fillna(np.inf).values
            ts = tbest.s.reindex(eg.e_partner.values).fillna(-np.inf).values
            tacc = tq[:, None] <= GRID[None]
            n_e = eq.sum(0)
            n_eo = (eq & ~tacc).sum(0)
            n_egt = (eq & tacc & (eg.score.values >= ts)[:, None]).sum(0)
        else:
            g = s.groupby(["key", "side_e"]).agg(q=("q", "min"), s=("score", "max")).reset_index()
            T = g[~g.side_e].set_index("key")
            E = g[g.side_e].set_index("key")
            tq = T.q.reindex(E.index).fillna(np.inf).values
            ts = T.s.reindex(E.index).fillna(-np.inf).values
            eq = E.q.values[:, None] <= GRID[None]
            tacc = tq[:, None] <= GRID[None]
            n_t = (T.q.values[:, None] <= GRID[None]).sum(0)
            n_e = eq.sum(0)
            n_eo = (eq & ~tacc).sum(0)
            n_egt = (eq & tacc & (E.s.values >= ts)[:, None]).sum(0)
        out[r] = np.stack([n_t, n_e, n_eo, n_egt], axis=1)
    if check:
        for a in (0.01, 0.05):
            k = int(np.flatnonzero(np.isclose(GRID, a))[0])
            off = np.array(paired_official(tab, level, a).totals(), float)
            if not np.allclose(out[nruns, k], off):
                raise AssertionError(f"vectorized paired counts differ from vulcan.stats ({level}, {a}): "
                                     f"{out[nruns, k]} vs {off}")
    _ = allbits
    return out[:nruns]


def fdp_paired(c):
    n_t, n_e, n_eo, n_egt = c[..., 0], c[..., 1], c[..., 2], c[..., 3]
    den = n_t + n_e
    return np.where(den > 0, (n_e + n_eo + 2 * n_egt) / np.maximum(den, 1), 0.0)


def matched_ids(fdp, ids, f):
    """First crossing of fdp above f along the grid, linear interpolation of ids. -> (ids, q_threshold, reached)."""
    over = np.flatnonzero(fdp > f)
    if len(over) == 0:
        return float(ids[-1]), float(GRID[-1]), False
    k = int(over[0])
    if k == 0:
        return 0.0, 0.0, True
    w = (f - fdp[k - 1]) / (fdp[k] - fdp[k - 1])
    return float(ids[k - 1] + w * (ids[k] - ids[k - 1])), float(GRID[k - 1] + w * (GRID[k] - GRID[k - 1])), True


# ------------------------------------------------------------------------------------------- verdicts
def arm_eval(reps, alpha, estimator, r=1.0):
    out = {}
    for c in ENT.CI_CONVENTIONS if hasattr(ENT, "CI_CONVENTIONS") else ("a",):
        kw = {"ci_convention": c} if hasattr(ENT, "CI_CONVENTIONS") else {}
        out[c] = ENT.evaluate_arm(reps, alpha, estimator, r=r, n_boot=N_BOOT, seed=SEED, **kw)
    return out


def worst(vs):
    return min(vs, key=CONSERVATIVE.index)


# ------------------------------------------------------------------------------------------- main
def evaluate(work, runs, decoy_seeds, entrap_seeds, methods, unit_loader, log=print, ours: str = "vulcan") -> dict:
    """unit_loader(method, rep) -> peptidoform DataFrame (or None if missing)."""
    reps = replicate_names(decoy_seeds, entrap_seeds)
    R, J, S = len(runs), len(decoy_seeds), len(entrap_seeds)
    pm = ProteinMap(os.path.join(work, "index", "full"))
    rls = {}
    out = {"units": __doc__.split("ENTRAPMENT")[0].strip(), "runs": runs, "decoy_seeds": list(decoy_seeds),
           "entrap_seeds": list(entrap_seeds), "grid_q": GRID.tolist(), "methods": {}}
    curves = {}
    for m in methods:
        log(f"evaluate {m}")
        ids_rj = {lv: np.zeros((R, J, len(GRID))) for lv in LEVELS}
        cnt = {lv: np.zeros((R, J, S, len(GRID), 4)) for lv in LEVELS}
        ds_ids = {lv: {str(a): [] for a in ALPHAS} for lv in LEVELS}
        pairs = {lv: {a: {} for a in ALPHAS} for lv in LEVELS}
        comb = {lv: {a: [] for a in ALPHAS} for lv in ("peptidoform", "sequence")}
        rC = {}
        missing = []
        for rep in reps:
            kind, j, s = parse_rep(rep)
            ji = list(decoy_seeds).index(j)
            if rep not in rls:
                if os.environ.get("VULCAN_T1_NO_LABEL_CACHE") == "1":   # memory-bounded mode (large dev sets)
                    rls.clear()
                rls[rep] = RepLabels(os.path.join(work, "index", rep))
            rl = rls[rep]
            pf = unit_loader(m, rep, rl)
            if pf is None:
                missing.append(rep)
                continue
            tabs = {"peptidoform": pf, "sequence": seq_units(pf)}
            if kind != "C":
                tabs["protein"] = protein_table(tabs["sequence"], pm)
            if kind == "prod":
                for lv in LEVELS:
                    ids_rj[lv][:, ji] = ids_per_run(tabs[lv], R, lv)
                    for a in ALPHAS:
                        ds_ids[lv][str(a)].append(ids_dataset(tabs[lv], lv, a))
            elif kind == "P":
                si = list(entrap_seeds).index(s)
                for lv in LEVELS:
                    cnt[lv][:, ji, si] = paired_grid(tabs[lv], lv, R)
                    for a in ALPHAS:
                        pairs[lv][a][(j, s)] = paired_official(tabs[lv], lv, a)
            else:
                rC = {"peptidoform": rl.n_fwd_pf[-1] / rl.n_fwd_pf[0], "sequence": rl.n_fwd_seq[-1] / rl.n_fwd_seq[0]}
                for lv in ("peptidoform", "sequence"):
                    for a in ALPHAS:
                        comb[lv][a].append(ENT.Replicate(combined_official(tabs[lv], lv, a), -1, j))
        n_flips = len(ORIENTATION_FLIPS)
        ORIENTATION_FLIPS.clear()
        if missing:
            log(f"  {m}: missing replicates {missing}")
            out["methods"][m] = {"missing_replicates": missing}
            continue
        rec = {"ids_dataset_level": {}, "verdicts": {}, "arm_C_r": rC, "n_replicates_score_sign_flipped": n_flips}
        for lv in LEVELS:
            rec["ids_dataset_level"][lv] = {a: {"per_decoy_seed": v, "mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1))}
                                            for a, v in ds_ids[lv].items()}
            rec["verdicts"][lv] = {}
            for a in ALPHAS:
                reps_p = [ENT.Replicate(c, s, j) for (j, s), c in pairs[lv][a].items()]
                P = arm_eval(reps_p, a, "paired")
                v = {"P": {c: {k: P[c][k] for k in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}
                           for c in P}}
                if lv != "protein":
                    C = arm_eval(comb[lv][a], a, "combined", r=rC[lv])
                    v["C"] = {c: {k: C[c][k] for k in ("mean_fdp", "ci95", "verdict", "mean_lb", "lb_ci95", "mean_n_t")}
                              for c in C}
                    v["overall"] = {c: ENT.combine_arms({"P": P[c]["verdict"], "C": C[c]["verdict"]},
                                                        "static_low_capacity") for c in P}
                else:
                    v["overall"] = {c: P[c]["verdict"] for c in P}
                v["deciding"] = v["overall"][CI_CONVENTION_DECIDING]
                v["deciding_conservative"] = worst(list(v["overall"].values()))
                v["consistent_with_control_P"] = {c: bool(P[c]["ci95"][0] <= a) for c in P}
                rec["verdicts"][lv][str(a)] = v
        out["methods"][m] = rec
        curves[m] = {"ids": ids_rj, "cnt": cnt}
    numerators = ours.split(",")          # E7/v2 (dev): several numerators; the first is "noninferiority"
    out["noninferiority"] = noninferiority(out, curves, methods, numerators[0])
    out["noninferiority_numerator"] = numerators[0]
    if len(numerators) > 1:
        out["noninferiority_by_numerator"] = {o: noninferiority(out, curves, methods, o) for o in numerators}
    LAST_CURVES.clear()
    LAST_CURVES.update(curves)
    return out


def noninferiority(out, curves, methods, ours: str = "vulcan"):
    if ours not in curves:
        return {}
    rng = np.random.default_rng(SEED)
    R, J = curves[ours]["ids"]["sequence"].shape[:2]
    draws = [(rng.integers(0, R, R), rng.integers(0, J, J)) for _ in range(N_BOOT)]
    kq = {a: int(np.flatnonzero(np.isclose(GRID, a))[0]) for a in ALPHAS}

    def own_q(m, lv, a, ri=None, ji=None):
        x = curves[m]["ids"][lv][:, :, kq[a]]
        if ri is not None:
            x = x[ri][:, ji]
        return x.sum(0).mean()

    def matched(m, lv, f, ri=None, ji=None):
        c = curves[m]["cnt"][lv]
        x = curves[m]["ids"][lv]
        if ri is not None:
            c, x = c[ri][:, ji], x[ri][:, ji]
        fdp = fdp_paired(c.sum(0)).mean(axis=(0, 1))          # sum over runs, mean over (decoy, entrapment) replicates
        return matched_ids(fdp, x.sum(0).mean(0), f)

    res = {}
    for m in methods:
        if m == ours or m not in curves:
            continue
        res[m] = {}
        for lv in LEVELS:
            d = {"own_q": {}, "matched_fdp": {}}
            for a in (0.01, 0.05, 0.001):
                o, c = own_q(ours, lv, a), own_q(m, lv, a)
                bs = np.array([own_q(ours, lv, a, ri, ji) / max(own_q(m, lv, a, ri, ji), 1e-9) for ri, ji in draws])
                lo, hi = np.percentile(bs, [2.5, 97.5])
                vv = out["methods"][m]["verdicts"][lv].get(str(a), {})
                verdict = vv.get("deciding")
                d["own_q"][str(a)] = {"ours_ids": o, "comparator_ids": c, "ratio": o / c if c else None,
                                      "ci95": [float(lo), float(hi)], "lower_ci_ge_0.97": bool(lo >= NI_MARGIN),
                                      "comparator_verdict_deciding": verdict,
                                      "counts_for_noninferiority": verdict == "supported",
                                      "comparator_verdict_conservative": vv.get("deciding_conservative")}
            for f in FMATCH:
                (o, qo, ro), (c, qc, rc) = matched(ours, lv, f), matched(m, lv, f)
                bs = []
                for ri, ji in draws:
                    ob, cb = matched(ours, lv, f, ri, ji)[0], matched(m, lv, f, ri, ji)[0]
                    bs.append(ob / cb if cb > 0 else np.nan)
                bs = np.array(bs)
                lo, hi = np.nanpercentile(bs, [2.5, 97.5])
                d["matched_fdp"][str(f)] = {"ours_ids": o, "comparator_ids": c, "ratio": o / c if c else None,
                                            "ci95": [float(lo), float(hi)], "lower_ci_ge_0.97": bool(lo >= NI_MARGIN),
                                            "ours_q_at_match": qo, "comparator_q_at_match": qc,
                                            "fdp_reached": [ro, rc], "n_boot_nan": int(np.isnan(bs).sum())}
            res[m][lv] = d
    # matched-FDP curves (point estimates) for plotting / inspection
    out["matched_fdp_curves"] = {}
    for m in curves:
        out["matched_fdp_curves"][m] = {}
        for lv in LEVELS:
            c = curves[m]["cnt"][lv]
            out["matched_fdp_curves"][m][lv] = {
                "fdp_P_mean": fdp_paired(c.sum(0)).mean(axis=(0, 1)).round(6).tolist(),
                "ids_run_summed_mean": curves[m]["ids"][lv].sum(0).mean(0).round(2).tolist()}
    return res
