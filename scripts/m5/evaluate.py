#!/usr/bin/env python3
"""M5 L4: XL FDR (vulcan.stats.xl, TT/TD/DD, merge-then-estimate, self vs heteromeric) at CSM -> peptide pair ->
residue pair -> PPI, and the true FDP of the accepted target-target entities from the library design (truth.py).

    python scripts/m5/evaluate.py --which D3b --csm X.csm.parquet --index DIR --fasta F --out metrics_part.json

Estimator (DESIGN L4, frozen): FDR = (TD + 1 - DD)/TT clipped at 0 (xl_qvalues plus_one=True, `q_tddd`), always
reported with the conservative (TD + 1)/TT (`q_td`). A level/class is *reportable* only if it accepts TT >= 50 at
the nominal level; other cells are kept in the JSON with reportable=false and are not headline numbers.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from truth import Truth  # noqa: E402

from vulcan.stats.xl import LEVELS, xl_level_fdr  # noqa: E402

ALPHAS = (0.01, 0.05)
MIN_TT = 50


def read_fasta(path):
    out, acc, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if acc:
                out[acc] = "".join(buf)
            acc, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if acc:
        out[acc] = "".join(buf)
    return out


class ProteinMap:
    """target peptide (index key) -> accessions, and positions (I = L)."""

    def __init__(self, index_dir, fasta):
        from vulcan.index.build import _read
        pr = _read(os.path.join(index_dir, "proteins.arrow"))
        pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
        accs = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
        off = pr.column("prot_off").to_numpy()
        ln = pr.column("prot_len").to_numpy()
        self.prots = {}
        for p, o, n in zip(pr.column("peptide").to_pylist(), off, ln):
            self.prots[p.replace("I", "L")] = sorted({accs[i] for i in pidx[o:o + n]})
        self.seqs = {a: s.replace("I", "L") for a, s in read_fasta(fasta).items()}
        self.raw_seqs = read_fasta(fasta)

    def proteins(self, pep):
        p = pep.replace("I", "L")
        if p in self.prots:
            return self.prots[p]
        return sorted(a for a, s in self.seqs.items() if p in s)

    def position(self, acc, pep):
        return self.seqs[acc].find(pep.replace("I", "L"))


def assign_proteins(pm: ProteinMap, pep_a, pep_b):
    """Shared protein if the two peptides share one (-> self link), else the first accession of each."""
    pa, pb, oa, ob = [], [], [], []
    for x, y in zip(pep_a, pep_b):
        A, B = pm.proteins(x), pm.proteins(y)
        common = sorted(set(A) & set(B))
        a = common[0] if common else (A[0] if A else "?")
        b = common[0] if common else (B[0] if B else "?")
        pa.append(a); pb.append(b)
        oa.append(pm.position(a, x) if a != "?" else -1)
        ob.append(pm.position(b, y) if b != "?" else -1)
    return np.array(pa), np.array(pb), np.array(oa), np.array(ob)


def residue_groups(truth: Truth, pm: ProteinMap):
    """(protein, 0-based residue) -> set of library groups whose peptide covers it."""
    rg = {}
    for p, g in truth.lib:
        for a, s in pm.seqs.items():
            i = s.find(p)
            while i >= 0:
                for k in range(i, i + len(p)):
                    rg.setdefault((a, k), set()).add(g)
                i = s.find(p, i + 1)
    return rg


def evaluate(which, score, pep_a, pep_b, dec_a, dec_b, site_a, site_b, pm: ProteinMap, csm_id=None, label=""):
    """All inputs per CSM (one row per spectrum, already the winning explanation). pep_* are TARGET-side sequences
    (the decoy's own target for decoy peptides); dec_* booleans; site_* 0-based within peptide."""
    truth = Truth(which, pm.raw_seqs)
    rg = residue_groups(truth, pm)
    pa, pb, oa, ob = assign_proteins(pm, pep_a, pep_b)
    abs_a, abs_b = oa + np.asarray(site_a), ob + np.asarray(site_b)
    prot_a = np.where(dec_a, np.char.add("DECOY_", pa.astype(str)), pa)
    prot_b = np.where(dec_b, np.char.add("DECOY_", pb.astype(str)), pb)
    base_of = {f"DECOY_{x}": x for x in set(pa.tolist()) | set(pb.tolist())}
    # truth per CSM (TT only meaningful)
    cls = np.array([truth.link_class(x, y) for x, y in zip(pep_a, pep_b)])
    res_true = np.array([bool(rg.get((a, i), set()) & rg.get((b, j), set()))
                         for a, i, b, j in zip(pa, abs_a, pb, abs_b)])
    ppi_true = np.array([truth.ppi_true(a, b) for a, b in zip(pa, pb)])
    out = {"label": label, "n_csm_input": int(len(score)),
           "n_tt": int((~dec_a & ~dec_b).sum()), "n_td": int((dec_a ^ dec_b).sum()), "n_dd": int((dec_a & dec_b).sum())}
    pep_key_a = np.char.add(np.asarray(pep_a).astype(str), "")
    lv = {}
    for level in LEVELS:
        r = xl_level_fdr(level, score, prot_a, prot_b, dec_a, dec_b, base_of=base_of, plus_one=True,
                         pep_a=pep_key_a, pep_b=np.asarray(pep_b).astype(str),
                         site_a=abs_a, site_b=abs_b, csm_id=csm_id)
        # truth per merged entity: take it from the best-scoring member (merge_level keeps the best member)
        key_of = _entity_keys(level, prot_a, prot_b, dec_a, dec_b, pep_key_a, np.asarray(pep_b).astype(str),
                              abs_a, abs_b, csm_id, len(score))
        best_member = {}
        o = np.argsort(-np.asarray(score), kind="stable")
        for i in o:
            best_member.setdefault(key_of[i], i)
        lv[level] = {}
        for cl in ("self", "hetero"):
            e = r[cl]
            tt = e["n_decoy"] == 0
            if level in ("csm", "peptide_pair"):
                true = np.array([cls[best_member[k]] == "true" for k in e["key"]])
            elif level == "residue_pair":
                true = np.array([res_true[best_member[k]] for k in e["key"]])
            else:
                true = np.array([ppi_true[best_member[k]] for k in e["key"]])
            cell = {"n_entities": int(len(e["key"])), "n_tt_total": int(tt.sum())}
            for a in ALPHAS:
                for est, q in (("tddd_plus1", e["q_tddd"]), ("td_plus1", e["q_td"])):
                    acc = tt & (q <= a)
                    n = int(acc.sum())
                    nf = int((acc & ~true).sum())
                    c = {"accepted_tt": n, "false": nf, "true_fdp": (nf / n) if n else None,
                         "accepted_td": int(((e["n_decoy"] == 1) & (q <= a)).sum()),
                         "accepted_dd": int(((e["n_decoy"] == 2) & (q <= a)).sum()),
                         "reportable": n >= MIN_TT}
                    if level in ("csm", "peptide_pair"):
                        kinds = np.array([cls[best_member[k]] for k in e["key"]])
                        c["false_breakdown"] = {k: int((acc & (kinds == k)).sum())
                                                for k in ("between_group", "non_library")}
                    cell[f"{est}@{a}"] = c
            lv[level][cl] = cell
    out["levels"] = lv
    return out


def _entity_keys(level, pa, pb, da, db, pep_a, pep_b, sa, sb, csm_id, n):
    """Re-derive merge_level's key per CSM (same construction as vulcan.stats.xl.merge_level)."""
    from vulcan.stats.xl import _unordered
    pa, pb = np.asarray(pa).astype(str), np.asarray(pb).astype(str)
    da, db = np.asarray(da, bool), np.asarray(db, bool)
    if level == "csm":
        return np.asarray(csm_id if csm_id is not None else np.arange(n)).astype(str)
    if level == "peptide_pair":
        x, y = _unordered(np.char.add(np.char.add(pep_a, "|"), da.astype(str)),
                          np.char.add(np.char.add(pep_b, "|"), db.astype(str)))
    elif level == "residue_pair":
        x, y = _unordered(np.char.add(np.char.add(pa, ":"), np.asarray(sa).astype(str)),
                          np.char.add(np.char.add(pb, ":"), np.asarray(sb).astype(str)))
    else:
        x, y = _unordered(pa, pb)
    return np.char.add(np.char.add(x, "~"), y)


def ours_inputs(csm_path, index_dir, homo=False):
    """Our CSM table -> evaluate() inputs (rank-1 XL explanations that beat the linear/mono/loop competitor).
    homo=False: hetero-family pairs only (the TT/TD/DD-estimable set); True: homo-family only; None: all."""
    import pyarrow.parquet as pq

    from vulcan.index.build import _read
    t = pq.read_table(csm_path).to_pandas()
    t = t[(t["rank"] == 1) & t["xl_wins"] & (t["score"] > 0)]
    lab = _read(os.path.join(index_dir, "labels.arrow"))
    idx = _read(os.path.join(index_dir, "peptidoforms.arrow"))
    assert np.array_equal(lab.column("pf_id").to_numpy(), idx.column("pf_id").to_numpy())
    dec = lab.column("is_decoy").to_numpy(zero_copy_only=False)
    tseq = np.array(lab.column("target_seq").to_pylist(), dtype=object)
    ra, rb = t["pf_row_a"].values, t["pf_row_b"].values
    if homo is not None:
        m = t["homo_family"].values == homo
        t, ra, rb = t[m], ra[m], rb[m]
    return dict(score=t["score"].values, pep_a=tseq[ra].astype(str), pep_b=tseq[rb].astype(str),
                dec_a=dec[ra], dec_b=dec[rb], site_a=t["site_a"].values, site_b=t["site_b"].values,
                csm_id=(t["run"] + ":" + t["scan"].astype(str)).values), t


def homo_summary(which, pm, inp, t):
    """Homo-family CSMs (same peptide on both sides up to target<->decoy pairing): not estimable by TT/TD/DD.
    Counts, true fraction, and an EXPLORATORY (non-preregistered, not used for any claim) DD/TT-threshold curve."""
    truth = Truth(which, pm.raw_seqs)
    da, db = np.asarray(inp["dec_a"]), np.asarray(inp["dec_b"])
    tt = ~da & ~db
    cls = np.array([truth.link_class(x, y) for x, y in zip(inp["pep_a"], inp["pep_b"])])
    out = {"n": int(len(da)), "n_tt": int(tt.sum()), "n_td": int((da ^ db).sum()), "n_dd": int((da & db).sum()),
           "tt_true": int((tt & (cls == "true")).sum()), "tt_false": int((tt & (cls != "true")).sum())}
    s = np.asarray(inp["score"])
    o = np.argsort(-s, kind="stable")
    ctt = np.cumsum(tt[o]); cdd = np.cumsum((da & db)[o]); cf = np.cumsum((tt & (cls != "true"))[o])
    q = np.minimum.accumulate((cdd / np.maximum(ctt, 1))[::-1])[::-1]
    for a in ALPHAS:
        k = np.flatnonzero(q <= a)
        n = int(ctt[k[-1]]) if len(k) else 0
        out[f"exploratory_dd_over_tt@{a}"] = {"accepted_tt": n, "true_fdp": (float(cf[k[-1]]) / n) if n else None}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", required=True)
    ap.add_argument("--csm", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pm = ProteinMap(a.index, a.fasta)
    inp, t = ours_inputs(a.csm, a.index, homo=False)
    res = evaluate(a.which, pm=pm, label="vulcan", **inp)
    res["homo_family"] = homo_summary(a.which, pm, *ours_inputs(a.csm, a.index, homo=True))
    json.dump(res, open(a.out, "w"), indent=1)
    for lvl, d in res["levels"].items():
        for cl, c in d.items():
            x1, x5 = c["tddd_plus1@0.01"], c["tddd_plus1@0.05"]
            print(f"{lvl:13s} {cl:6s} 1%: TT={x1['accepted_tt']:5d} FDP={x1['true_fdp']}  "
                  f"5%: TT={x5['accepted_tt']:5d} FDP={x5['true_fdp']}")


if __name__ == "__main__":
    main()
