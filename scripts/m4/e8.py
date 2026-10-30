#!/usr/bin/env python3
"""E8 on dev (D2): FLR̂ (decoy-amino-acid, DESIGN §6.5, vulcan.stats.flr) vs true FLR, site recall at matched true
FLR, our L3 operator vs pyAscore (same PSMs, same decoy-AA candidates, same neutral-loss treatment).

    python scripts/m4/e8.py --psms /tmp/vulcan_m4/psm_d2_phospho --index /tmp/vulcan_m4/index/d2_phospho_prod --out RES

Pipeline per accession (PXD000138 4 pools, PXD009449 Ymod, PXD007058 4 pools; all runs of an accession pooled):
 1. Rank-1 PSMs of the phospho closed search (STY phospho, M-ox, pep-N-term acetyl; T + pseudo-reverse decoys).
 2. Localization-agnostic peptidoform FDR: paired double competition, unit = (sequence, composition), pair = (forward
    sequence, composition), score = hyperscore, hash tie-break seed 11; keep target PSMs of units with q <= 1% that
    carry >= 1 phospho.
 3. Localization (L3) with candidates S/T/Y + Ala (primary) and, separately, S/T/Y + Leu/Ile (secondary).
 4. FLR̂ curve (flr.flr_curve, monotone, PSM-level acceptor counts), true FLR from scripts/m4/truth.py.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
import truth as TR  # noqa: E402

from vulcan.index.build import L1Index  # noqa: E402
from vulcan.ptm import localize as LZ  # noqa: E402
from vulcan.stats import flr as FLR  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402

NTERM_BIT = np.uint64(1) << np.uint64(63)
ACETYL = 42.010565
OX = 15.994915
ALPHAS = (0.01, 0.05)


def read_labels(index, rows):
    rows = np.asarray(rows, np.int64)
    out = {"is_decoy": np.zeros(len(rows), bool), "fwd_seq": np.empty(len(rows), object)}
    with pa.memory_map(os.path.join(index, "labels.arrow")) as src:
        f = pa.ipc.open_file(src)
        base = 0
        for b in range(f.num_record_batches):
            rb = f.get_batch(b)
            n = rb.num_rows
            m = (rows >= base) & (rows < base + n)
            if m.any():
                loc = pa.array(rows[m] - base)
                out["is_decoy"][m] = rb.column(rb.schema.get_field_index("is_decoy")).take(loc).to_numpy(
                    zero_copy_only=False)
                out["fwd_seq"][m] = rb.column(rb.schema.get_field_index("fwd_seq")).take(loc).to_pylist()
            base += n
    return out


def decode(l1, rows):
    seqs, phos, oxs, nts = [], [], [], []
    for r in rows:
        s = l1.sequence(int(r))
        m = int(l1.mask[r])
        nt = bool(m >> 63 & 1)
        bits = [k for k in range(len(s)) if m >> k & 1]
        seqs.append(s)
        phos.append([k for k in bits if s[k] in "STY"])
        oxs.append([k for k in bits if s[k] == "M"])
        nts.append(nt)
    return seqs, phos, oxs, nts


def load_spectra(stores, keys):
    """keys: list of (run, spec_idx) -> list of (mz, intensity) arrays (raw centroids from the store)."""
    from vulcan.gpu.engine import load_store
    out = {}
    need = {}
    for run, si in keys:
        need.setdefault(run, set()).add(int(si))
    for s in stores:
        sp = load_store(s)
        if sp.run not in need:
            continue
        for si in need[sp.run]:
            a, b = sp.offsets[si], sp.offsets[si + 1]
            out[(sp.run, si)] = (sp.mz[a:b].astype(np.float64), sp.intensity[a:b].astype(np.float64))
    return [out[k] for k in keys]


def run_ours(psm, spectra, acceptors, p: LZ.LocParams):
    rm = [LZ.residue_masses(s, C.STATIC, {k: OX for k in ox}) for s, ox in zip(psm["seq"], psm["ox"])]
    cp = [LZ.candidate_positions(s, acceptors) for s in psm["seq"]]
    nt = np.array([ACETYL if x else 0.0 for x in psm["nterm"]])
    k = np.array([len(x) for x in psm["phos"]])
    peaks = LZ.prepare_peaks([x[0] for x in spectra], [x[1] for x in spectra], p)
    t = time.perf_counter()
    sc = LZ.score_configs(rm, nt, cp, k, psm["charge"], np.arange(len(rm)), peaks, p)
    sm = LZ.summarize(sc, cp)
    dt = time.perf_counter() - t
    sites = []
    for i, s in enumerate(sm):
        if s is None:
            continue
        for pos, d in zip(s["best_cfg"], s["site_delta"]):
            sites.append((i, int(pos), float(d)))
    return sites, dt, sum(1 for s in sm if s is None)


def run_pyascore(psm, spectra, acceptors, p: LZ.LocParams):
    sys.path.insert(0, os.path.join(C.WORK, "pylib"))
    from pyascore import PyAscore
    grp = acceptors.replace("L", "LI") if "L" in acceptors else acceptors
    sc = PyAscore(bin_size=100., n_top=10, mod_group=grp, mod_mass=LZ.PHOSPHO, mz_error=0.02)
    sc.add_neutral_loss(grp.lower(), LZ.H3PO4)      # same residue-blind neutral loss as our operator
    sites = []
    t = time.perf_counter()
    n_fail = 0
    for i, (seq, ox, nt, ph, z) in enumerate(zip(psm["seq"], psm["ox"], psm["nterm"], psm["phos"], psm["charge"])):
        pos = [j + 1 for j, c in enumerate(seq) if c == "C"] + [j + 1 for j in ox] + ([0] if nt else [])
        mass = [57.021464] * seq.count("C") + [OX] * len(ox) + ([ACETYL] if nt else [])
        o = np.argsort(pos, kind="stable")
        mz, it = spectra[i]
        try:
            sc.score(mz_arr=np.ascontiguousarray(mz), int_arr=np.ascontiguousarray(it), peptide=seq,
                     n_of_mod=len(ph), max_fragment_charge=max(1, min(int(z) - 1, p.max_frag_charge)),
                     aux_mod_pos=np.array(pos, np.uint32)[o], aux_mod_mass=np.array(mass, np.float32)[o])
        except Exception:
            n_fail += 1
            continue
        best = sc.best_sequence
        # positions of [80] marks in the best sequence
        ppos, j, k = [], 0, 0
        while j < len(best):
            c = best[j]
            if c == "[":
                e = best.index("]", j)
                if best[j + 1:e] == "80":
                    ppos.append(k - 1)
                j = e + 1
                continue
            if c.isalpha():
                k += 1
            j += 1
        asc = list(np.asarray(sc.ascores, float))
        for q, a in zip(sorted(ppos), asc):
            sites.append((i, int(q), float(a)))
    return sites, time.perf_counter() - t, n_fail


def evaluate(sites, psm, truth_of, decoy):
    """FLR̂ curve over all sites; true FLR / recall on truth-eligible mono-phospho PSMs."""
    if not sites:
        return {}
    pi = np.array([s[0] for s in sites])
    pos = np.array([s[1] for s in sites])
    sc = np.array([s[2] for s in sites])
    seqs = psm["seq"]
    res = np.array([seqs[i][q] for i, q in zip(pi, pos)])
    nph = np.array([len(psm["phos"][i]) for i in pi])
    out = {}
    for cls, msk in (("mono", nph == 1), ("multi", nph > 1)):
        if not msk.any():
            continue
        f = FLR.flr_curve(sc[msk], res[msk], pi[msk], seqs, decoy_residue=decoy)
        r = res[msk]
        rr = np.where(np.isin(r, ["I"]) & (decoy == "L"), "L", r)
        is_sty = np.isin(rr, list("STY"))
        is_dec = rr == decoy
        d = {"n_sites": int(msk.sum()), "n_sty": int(is_sty.sum()), "n_decoy": int(is_dec.sum())}
        for a in ALPHAS:
            acc = (f <= a) & is_sty
            d[f"n_sty_at_flrhat_{a}"] = int(acc.sum())
        if cls == "mono":
            tp = np.array([truth_of[i] for i in pi[msk]], dtype=object)       # truth position or None
            elig = np.array([t is not None for t in tp])
            correct = np.array([t is not None and t == q for t, q in zip(tp, pos[msk])])
            d["n_eligible_psms"] = int(len({int(i) for i in pi[msk][elig]}))
            # (a) FLR̂ threshold from all mono sites, truth on eligible
            for a in ALPHAS:
                acc = (f <= a) & is_sty & elig
                n = int(acc.sum())
                d[f"all_thr_{a}"] = {"n_sty_eligible": n, "true_flr": float((acc & ~correct).sum() / n) if n else None}
            # (b) FLR̂ computed on the eligible subset only
            fe = FLR.flr_curve(sc[msk][elig], r[elig], pi[msk][elig], seqs, decoy_residue=decoy)
            ce, se = correct[elig], is_sty[elig]
            for a in ALPHAS:
                acc = (fe <= a) & se
                n = int(acc.sum())
                d[f"elig_thr_{a}"] = {"n_sty": n, "true_flr": float((acc & ~ce).sum() / n) if n else None,
                                      "flrhat_at_thr": float(fe[acc].max()) if n else None}
            # overall (no threshold) on eligible: true FLR of all STY sites, FLR̂ of the whole list
            fl_all = FLR.decoy_aa_flr(r[elig], pi[msk][elig], seqs, decoy_residue=decoy)
            d["elig_whole_list"] = {"true_flr": float((se & ~ce).sum() / max(se.sum(), 1)), "flrhat": fl_all["flr"],
                                    "n_sty": int(se.sum()), "n_decoy": int((~se).sum())}
            # recall at matched true FLR (largest score threshold set with true FLR <= a among STY sites)
            o = np.argsort(-sc[msk][elig], kind="stable")
            s_sorted = sc[msk][elig][o]
            wrong = np.cumsum(se[o] & ~ce[o])
            nsty = np.cumsum(se[o])
            good = np.cumsum(se[o] & ce[o])
            last_of_tie = np.r_[s_sorted[1:] != s_sorted[:-1], True]
            tf = np.where(nsty > 0, wrong / np.maximum(nsty, 1), 0.0)
            n_true = d["n_eligible_psms"]
            for a in ALPHAS:
                okk = last_of_tie & (tf <= a)
                g = int(good[okk].max()) if okk.any() else 0
                d[f"recall_at_true_flr_{a}"] = g / n_true if n_true else None
                d[f"correct_sites_at_true_flr_{a}"] = g
            d["recall_any_threshold"] = float(good[-1] / n_true) if n_true else None
            d["best_site_correct_frac"] = float(ce.mean()) if len(ce) else None
        out[cls] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", default=os.path.join(C.WORK, "psm_d2_phospho"))
    ap.add_argument("--index", default=C.index_dir("d2_phospho_prod"))
    ap.add_argument("--out", default=C.RES)
    ap.add_argument("--skip-pyascore", action="store_true")
    ap.add_argument("--loc-params", default="{}", help="JSON overrides of LocParams (dev grid)")
    ap.add_argument("--datasets", nargs="+", default=["PXD000138", "PXD009449", "PXD007058"])
    a = ap.parse_args()
    t0 = time.perf_counter()
    l1 = L1Index(a.index)
    t_l1 = time.perf_counter() - t0
    truths = TR.all_truth()
    p = LZ.LocParams(**json.loads(a.loc_params))
    report = {"truth_reports": {k: v["report"] for k, v in truths.items()}, "loc_params": p.__dict__,
              "datasets": {}, "timing_s": {"l1_load": t_l1}}
    stores_of = C.d2_stores()
    for acc in a.datasets:
        t_acc = time.perf_counter()
        runs = [C.run_name(r) for r in C.D2_RAW[acc]]
        tbl = pa.concat_tables([pq.read_table(os.path.join(a.psms, f"{r}.psms.parquet")) for r in runs])
        tbl = tbl.filter(pc.equal(tbl.column("rank"), 1))
        rows = tbl.column("pf_row").to_numpy()
        lab = read_labels(a.index, rows)
        seqs, phos, oxs, nts = decode(l1, rows)
        comp = [f"p{len(ph)}o{len(ox)}n{int(nt)}" for ph, ox, nt in zip(phos, oxs, nts)]
        unit = np.array([f"{s}|{c}" for s, c in zip(seqs, comp)])
        pair = np.array([f"{f}|{c}" for f, c in zip(lab["fwd_seq"], comp)])
        spec = np.array([f"{r}:{i}" for r, i in zip(tbl.column("run").to_pylist(), tbl.column("spec_idx").to_pylist())])
        r = double_competition(spec, unit, pair, tbl.column("hyperscore").to_numpy(), lab["is_decoy"], seed=11,
                               tiebreak=hash_tiebreak(unit, 11))
        acc_units = set(unit[r.accepted(0.01)].tolist())
        keep = np.array([(u in acc_units) and (not d) and len(ph) > 0 for u, d, ph in zip(unit, lab["is_decoy"], phos)])
        idx = np.flatnonzero(keep)
        run_col = np.array(tbl.column("run").to_pylist())
        psm = {"seq": [seqs[i] for i in idx], "phos": [phos[i] for i in idx], "ox": [oxs[i] for i in idx],
               "nterm": [nts[i] for i in idx], "charge": tbl.column("charge").to_numpy()[idx],
               "run": run_col[idx], "spec_idx": tbl.column("spec_idx").to_numpy()[idx]}
        # truth position per PSM (mono-phospho, by pool/run)
        T = truths[acc]["truth"]
        truth_of = []
        for rn, s, ph in zip(psm["run"], psm["seq"], psm["phos"]):
            key = (rn, s) if acc != "PXD000138" else (rn, s)
            truth_of.append(T.get(key) if len(ph) == 1 else None)
        t_sp = time.perf_counter()
        spectra = load_spectra(stores_of[acc], list(zip(psm["run"], psm["spec_idx"])))
        t_sp = time.perf_counter() - t_sp
        ds = {"n_rank1_psms": int(len(rows)), "n_peptidoforms_q01": int(len(r.accepted(0.01))),
              "n_phospho_psms_localized": int(len(idx)),
              "n_mono": int(sum(len(x) == 1 for x in psm["phos"])), "n_multi": int(sum(len(x) > 1 for x in psm["phos"])),
              "n_truth_eligible_psms": int(sum(t is not None for t in truth_of)),
              "search_mode_counts": {"nterm_acetyl": int(sum(psm["nterm"]))},
              "timing_s": {"load_spectra": t_sp}}
        # baseline: the search's own site placement (the identification's best peptidoform), STY only
        own_correct = [t is not None and len(ph) == 1 and ph[0] == t for t, ph in zip(truth_of, psm["phos"])]
        ds["search_rank1_site_correct_frac_on_eligible"] = (float(np.mean([c for c, t in zip(own_correct, truth_of)
                                                                           if t is not None]))
                                                            if any(t is not None for t in truth_of) else None)
        for decoy in ("A", "L"):
            acc_set = "STY" + decoy
            sites, dt, nskip = run_ours(psm, spectra, acc_set, p)
            ds[f"ours_{decoy}"] = {"eval": evaluate(sites, psm, truth_of, decoy), "time_s": dt, "n_skipped": nskip,
                                   "psms_per_s": len(idx) / dt if dt else None}
            if not a.skip_pyascore:
                sites, dt, nf = run_pyascore(psm, spectra, acc_set, p)
                ds[f"pyascore_{decoy}"] = {"eval": evaluate(sites, psm, truth_of, decoy), "time_s": dt,
                                           "n_failed": nf, "psms_per_s": len(idx) / dt if dt else None}
        ds["timing_s"]["total"] = time.perf_counter() - t_acc
        report["datasets"][acc] = ds
        print(acc, json.dumps({k: v for k, v in ds.items() if not k.startswith(("ours", "pyascore"))}), flush=True)
    from vulcan.hw import fingerprint
    report["hw"] = fingerprint()
    C.dump(report, os.path.join(a.out, "e8.json"))
    print(json.dumps(report["datasets"], indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
