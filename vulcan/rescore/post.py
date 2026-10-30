"""Post-scoring stage of production scoring, torch-free (runs in the GPU process or in CPU worker processes):
paired double competition + q-values (vulcan.stats.tdc.double_competition) and picked protein groups
(vulcan.stats.protein.picked_group_fdr), with EXACTLY the inputs of the frozen M3/M3e path
(scripts/m3e/bench_prod.py, scripts/m3/evaluate.py:protein_level): seed 11, hash tie-break on pf_id, spectrum =
spec_idx of the run, unit = pf_id, protein prefixes {0: '', 1: 'E_', 2: 'rev_', 3: 'rev_E_'}, arm-C foreign
peptides (entrap_seed < 0) without proteins.  Estimators are called unchanged; only the glue is faster
(per-index protein map loaded once, per-sequence protein lists memoized across runs).
"""
from __future__ import annotations

import os

import numpy as np

SEED = 11
_PRE = {0: "", 1: "E_", 2: "rev_", 3: "rev_E_"}


def load_labels(index_dir):
    import pyarrow.ipc as ipc
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    cols = ("pf_id", "family", "entrap_seed", "is_decoy", "is_entrapment", "tdc_pair", "entrap_pair", "target_seq")
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in cols}


class ProteinMap:
    """target peptide sequence -> protein accessions of one index (as evaluate.protein_maps), memoized."""

    def __init__(self, index_dir):
        import pyarrow.ipc as ipc
        with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
            t = rd.read_all()
        pep = t.column("peptide").to_pylist()
        self.pos = dict(zip(pep, range(len(pep))))
        self.off = t.column("prot_off").to_numpy()
        self.ln = t.column("prot_len").to_numpy()
        self.pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
        self.acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
        self.memo = {}

    def prefixed(self, seq, fam):
        k = (seq, fam)
        v = self.memo.get(k)
        if v is None:
            i = self.pos.get(seq)
            prots = [] if i is None else [self.acc[j] for j in self.pidx[self.off[i]:self.off[i] + self.ln[i]]]
            v = tuple(_PRE[fam] + p for p in prots)
            self.memo[k] = v
        return v


class IndexCtx:
    """Per-index state of the post stage (labels L4-only, protein map)."""

    def __init__(self, index_dir):
        self.index_dir = index_dir
        self.lab = load_labels(index_dir)
        self.pmap = ProteinMap(index_dir)


def compete(ctx: IndexCtx, pf_row, pf_id, spec_idx, score):
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    lab = ctx.lab
    return double_competition(np.asarray(spec_idx).astype(np.int64), pf_id, lab["tdc_pair"][pf_row], score,
                              lab["is_decoy"][pf_row], seed=SEED, tiebreak=hash_tiebreak(pf_id, SEED))


def protein_groups(ctx: IndexCtx, pf_row, unit, r, alphas=(0.01,)):
    """Same pep2prot / pep_score / pair_of / is_dec (incl. insertion order) as evaluate.protein_level."""
    from vulcan.stats.protein import picked_group_fdr
    lab = ctx.lab
    rows = pf_row[r.row]
    fam = lab["family"][rows]
    tseq = lab["target_seq"][rows]
    eseed = lab["entrap_seed"][rows]
    u_all = unit[r.row]
    sc = r.score
    pep2prot, pep_score = {}, {}
    pf = ctx.pmap.prefixed
    for i in np.flatnonzero(eseed >= 0).tolist():
        prots = pf(tseq[i], int(fam[i]))
        if not prots:
            continue
        u = int(u_all[i])
        pep2prot[u] = list(prots)
        s = float(sc[i])
        if u in pep_score:
            pep_score[u] = max(pep_score[u], s)
        else:
            pep_score[u] = max(-np.inf, s)
    allp = {p for v in pep2prot.values() for p in v}
    pair_of = {p: (p[4:] if p.startswith("rev_") else p) for p in allp}
    is_dec = {p: p.startswith("rev_") for p in allp}
    res = picked_group_fdr(pep2prot, pep_score, pair_of, is_dec, seed=SEED)
    out = {}
    for a in alphas:
        acc = (res.q <= a) & ~res.is_decoy & ~res.is_obsolete
        out[a] = {"n_groups": int(acc.sum())}
    return res, out


def post_stage(ctx: IndexCtx, pf_row, pf_id, spec_idx, score, alphas=(0.01,)):
    """-> (competition result, protein result, summary)."""
    r = compete(ctx, pf_row, pf_id, spec_idx, score)
    res, prot = protein_groups(ctx, pf_row, pf_id, r, alphas)
    summ = {f"ids_{a}": int(((r.q <= a) & ~r.is_decoy).sum()) for a in alphas}
    summ.update({f"proteins_{a}": prot[a]["n_groups"] for a in alphas})
    return r, res, summ
