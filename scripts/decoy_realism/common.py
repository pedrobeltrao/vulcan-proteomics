"""3g shared helpers (analysis only; L4-side: labels are read here, never by scoring code)."""
from __future__ import annotations

import csv
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from vulcan.stats import global_fdr as G  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402

SEED = 11
MAN = os.path.join(ROOT, "data", "manifests", "e9_d8.tsv")
IDX = {"P": "/tmp/vulcan_e2/index/human_contam_E2001", "C": "/tmp/vulcan_m3/index/human_contam_C"}
OUT = os.path.join(ROOT, "results", "2026-10-27_3g_decoy_realism")
SCR = "/tmp/vulcan_3g"


def manifest():
    man = list(csv.DictReader(open(MAN), delimiter="\t"))
    for r in man:
        r["run"] = f"{r['id']}__{re.sub(r'[^A-Za-z0-9._-]', '_', r['filename'])[:-4]}"
    return man


def labels(index_dir):
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    cols = ("pf_id", "family", "is_decoy", "is_entrapment", "tdc_pair", "target_seq")
    lab = {k: t.column(k).to_numpy(zero_copy_only=False) for k in cols}
    with ipc.open_file(os.path.join(index_dir, "peptidoforms.arrow")) as r:
        p = r.read_all()
    lab["seq_off"] = p.column("seq_off").to_numpy()
    lab["length"] = p.column("length").to_numpy()
    lab["missed"] = p.column("missed").to_numpy()
    lab["n_var"] = p.column("n_var").to_numpy()
    lab["mass"] = p.column("mass32").to_numpy()
    lab["residues"] = np.load(os.path.join(index_dir, "residues.npy"), mmap_mode="r")
    return lab


def row_seq(lab, row):
    o, n = int(lab["seq_off"][row]), int(lab["length"][row])
    return bytes(lab["residues"][o:o + n]).decode()


def _load_one(args):
    f, i, col, extra = args
    if not os.path.exists(f):
        return None
    t = pq.read_table(f, columns=["spec_idx", "pf_row", "pf_id", col, "charge"] + extra)
    spec = t.column("spec_idx").to_numpy()
    pf_id = t.column("pf_id").to_numpy()
    sc = t.column(col).to_numpy().astype(np.float64)
    w = G.psm_winners(spec, sc, hash_tiebreak(pf_id, SEED))
    out = {"run": np.full(len(w), i, np.int32), "spec_idx": spec[w].astype(np.int64),
           "pf_row": t.column("pf_row").to_numpy()[w].astype(np.int64), "pf_id": pf_id[w], "score": sc[w],
           "charge": t.column("charge").to_numpy()[w]}
    for e in extra:
        out[e] = t.column(e).to_numpy()[w]
    return out


def load_winners(work, arm, runs, col="score", procs=16, extra=()):
    jobs = [(os.path.join(work, "out", arm, r + ".scored.parquet"), i, col, list(extra)) for i, r in enumerate(runs)]
    with ProcessPoolExecutor(procs) as ex:
        parts = [p for p in ex.map(_load_one, jobs, chunksize=4) if p is not None]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def compete(d, lab, alphas=(0.01, 0.05)):
    """L6 global competition (exactly scripts/e9/analyze.py).  Returns glob, the threshold per alpha, and the
    index (into glob arrays) of pair winners at/above the loosest threshold."""
    rows = d["pf_row"]
    is_dec = lab["is_decoy"][rows]
    pair = lab["tdc_pair"][rows]
    tb = hash_tiebreak(d["pf_id"], SEED)
    g = G.global_competition(d["run"], d["pf_id"], pair, d["score"], is_dec, tb)
    cut = {}
    for a in alphas:
        acc = g.accepted(a)
        cut[a] = float(g.score[acc].min()) if acc.any() else math.inf
    return g, cut


def wilson(k, n, z=1.959964):
    if n == 0:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [c - h, c + h]


def coin(win_e, z=1.959964):
    """win_e: bool array, True = E (forward) side won the pair."""
    n = int(len(win_e))
    k = int(np.sum(win_e))
    return {"n": n, "n_E": k, "n_decE": n - k, "coin": k / n if n else None, "ci95": wilson(k, n, z),
            "z": (k - n / 2) / math.sqrt(n / 4) if n else None}


def boot_ratio(win_ref, win_sub, n_boot=2000, seed=SEED):
    """Retained fraction of the excess, (c_sub - .5)/(c_ref - .5); pair bootstrap (sub is resampled within itself and
    ref within itself, independently: conservative when sub is a subset of ref)."""
    rng = np.random.default_rng(seed)
    a, b = np.asarray(win_ref, float), np.asarray(win_sub, float)
    if len(b) == 0:
        return None
    est = (b.mean() - .5) / (a.mean() - .5)
    rs = []
    for _ in range(n_boot):
        ra = a[rng.integers(0, len(a), len(a))].mean() - .5
        rb = b[rng.integers(0, len(b), len(b))].mean() - .5
        rs.append(rb / ra if ra != 0 else np.nan)
    rs = np.array(rs)
    return {"retained": float(est), "ci95": [float(np.nanquantile(rs, .025)), float(np.nanquantile(rs, .975))]}


def boot_diff(w1, w2, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(w1, float), np.asarray(w2, float)
    rs = [a[rng.integers(0, len(a), len(a))].mean() - b[rng.integers(0, len(b), len(b))].mean() for _ in range(n_boot)]
    return {"diff": float(a.mean() - b.mean()), "ci95": [float(np.quantile(rs, .025)), float(np.quantile(rs, .975))]}


def cluster_boot(win_e, cluster, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    u, inv = np.unique(cluster, return_inverse=True)
    k = np.bincount(inv, weights=win_e.astype(float), minlength=len(u))
    n = np.bincount(inv, minlength=len(u)).astype(float)
    rs = []
    for _ in range(n_boot):
        s = rng.integers(0, len(u), len(u))
        rs.append(k[s].sum() / n[s].sum())
    return [float(np.quantile(rs, .025)), float(np.quantile(rs, .975))]


def std_coin(win_e, strata, ref_strata, n_boot=1000, seed=SEED):
    """Length-adjusted coin: direct standardisation of the per-stratum coins of a subset to the stratum mix of the
    reference set (strata absent from the subset are dropped and the weights renormalised).  Bootstrap CI."""
    win_e = np.asarray(win_e, float)
    us, cnt = np.unique(ref_strata, return_counts=True)
    wref = dict(zip(us.tolist(), (cnt / cnt.sum()).tolist()))

    def est(w, s):
        num = den = 0.0
        for k in np.unique(s):
            m = s == k
            num += wref.get(int(k), 0.0) * w[m].mean()
            den += wref.get(int(k), 0.0)
        return num / den if den else np.nan
    e = est(win_e, strata)
    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(win_e), len(win_e))
        rs.append(est(win_e[i], strata[i]))
    return {"coin_std": float(e), "ci95": [float(np.nanquantile(rs, .025)), float(np.nanquantile(rs, .975))]}


def len_stratum(L):
    L = np.asarray(L)
    return np.select([L <= 8, L <= 10, L <= 13, L <= 17], [0, 1, 2, 3], 4)


def auc(pos, neg):
    """Mann-Whitney AUC P(pos > neg) + 0.5 P(tie)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return None
    allv = np.concatenate([pos, neg])
    from scipy.stats import rankdata
    r = rankdata(allv)
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))
