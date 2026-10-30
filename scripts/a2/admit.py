#!/usr/bin/env python3
"""Amendment A2.2 spectrum admission (stratum A / B candidates) for one DN2 project, label-free.

Steps per certified spectrum (codes and order: vulcan.gates.a2_admission):
  (a) certification            scripts/a2/own_search.py (own-species search; never the gate DB labels)
  (b) identity screen          max over c in (T∪D)∩W of f(c) >= 0.20, or isobaric Hamming <= 2 to a T∪D peptide
                               (vulcan.gates.a2_identity; T∪D = GateDB, label-free)
  (d) family screen            exhaustive hyperscore search of the spectrum with its certified precursor over W:
                               max over ALL candidates (no label read) >= hyperscore(p) -> removed
  (c-pre) A1.3                 exhaustive top candidate and its pair partner (l1.partner, label-free) A1.3-near the
                               source proteome ∪ cRAP ∪ Frankenfield ∪ Panel H ∪ Panel B -> removed
  (e) single precursor         recorded here (certified mono m/z, z) and applied by the mixture builder
Label-free counts per project (stopping rule): n_A, n_A_with_candidate (A spectra with >= 1 exhaustive candidate),
tail_proxy (A spectra whose max T∪D hyperscore >= the host hyperscore 1% threshold, max over the two D0 hosts;
scripts/a2/host_threshold.py).  Output per project: <work>/admission/<acc>/<run>.npz (codes; A spectra with certified
precursor) and counts.json; the sha256 of every codes file and the swap-test result are committed (A2.9.4 step 3).

P1 swap test (--swap-test PROJECT): the gate index is rebuilt with vulcan.index.build.build(..., swap_labels=True)
(T and D files exchanged), its label-free side is checked byte-identical, and the whole admission of PROJECT is re-run
on it; PASS iff the membership sha256 is identical.  tests/a2/test_admit_swap.py runs the same check on synthetic T/D
lists through admit_project().

No per-label (target / decoy) quantity of the gate search is read, printed or stored here (A2.9.4).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from vulcan.gates import a2_admission as AD  # noqa: E402
from vulcan.gates import a2_config  # noqa: E402
from vulcan.gates.a2_identity import GateDB, PeptideHamming, hamming_isobaric_near, max_f_in_window  # noqa: E402

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")


# ------------------------------------------------------------------------------------------------ core (testable)
def admit_project(cert_by_run: dict, db: GateDB, ph: PeptideHamming, exhaustive, near_a13, thr_proxy: float,
                  cfg: dict) -> dict:
    """Label-free admission of one project.

    cert_by_run: {run: {"spec", "p_seq", "p_mask", "z", "mz_mono", "hyperscore_p", "n_spectra"}} (own_search output)
    exhaustive(run, cert) -> {"max_hs": float array (-inf if no candidate), "top_seq": list (None if none),
                              "partner_seq": list} -- max over ALL T∪D candidates in W; label-free.
    near_a13(list[str]) -> bool array (A1.3 reference of this project).
    Returns {"codes": {run: uint8}, "A": {run: {...}}, "counts": label-free counts}."""
    ad = cfg["admission"]
    codes, A_tab, has_c, proxy = {}, {}, {}, {}
    fcache: dict = {}
    for run in sorted(cert_by_run):
        c = cert_by_run[run]
        n = len(c["spec"])
        mf = np.array([max_f_in_window(db, c["p_seq"][i], int(c["p_mask"][i]), float(c["mz_mono"][i]),
                                       int(c["z"][i]), fcache, ad["identity_frag_ppm"], ad["identity_prec_ppm"],
                                       tuple(ad["identity_isotope_errors"])) for i in range(n)])
        rem_f = mf >= ad["identity_f_threshold"]
        rem_h = np.array([hamming_isobaric_near(ph, s, ad["hamming_variant_cap"]) for s in c["p_seq"]], bool)
        ex = exhaustive(run, c)
        mx = np.asarray(ex["max_hs"], float)
        rem_d = mx >= np.asarray(c["hyperscore_p"], float)
        seqs = [s for pair in zip(ex["top_seq"], ex["partner_seq"]) for s in pair if s is not None]
        near = dict(zip(seqs, np.asarray(near_a13(seqs), bool))) if seqs else {}
        rem_c = np.array([bool(near.get(t, False) or near.get(p, False)) if t is not None else False
                          for t, p in zip(ex["top_seq"], ex["partner_seq"])], bool)
        cd = AD.combine(int(c["n_spectra"]), c["spec"], rem_f, rem_h, rem_d, rem_c)
        codes[run] = cd
        a_local = cd[np.asarray(c["spec"], np.int64)] == AD.A
        A_tab[run] = {"spec": np.asarray(c["spec"], np.int64)[a_local], "mz_mono": np.asarray(c["mz_mono"])[a_local],
                      "z": np.asarray(c["z"])[a_local]}
        hc = np.zeros(cd.size, bool)
        hc[np.asarray(c["spec"], np.int64)] = np.isfinite(mx)
        tp = np.zeros(cd.size, bool)
        tp[np.asarray(c["spec"], np.int64)] = mx >= thr_proxy
        has_c[run], proxy[run] = hc, tp
    allc = np.concatenate([codes[r] for r in sorted(codes)]) if codes else np.zeros(0, np.uint8)
    counts = AD.label_free_counts(allc, np.concatenate([has_c[r] for r in sorted(codes)]) if codes else None,
                                  np.concatenate([proxy[r] for r in sorted(codes)]) if codes else None)
    return {"codes": codes, "A": A_tab, "counts": counts, "membership_sha256": AD.membership_hash(codes)}


# ------------------------------------------------------------------------------------------------ production
_GATE = {}


def gate_objects(idx_dir: str):  # pragma: no cover (heavy)
    if idx_dir not in _GATE:
        db = GateDB.from_l1(idx_dir)
        _GATE[idx_dir] = (db, PeptideHamming(db.peptides()))
    return _GATE[idx_dir]


def gpu_exhaustive(idx_dir: str, stores: dict):  # pragma: no cover (GPU)
    """exhaustive(run, cert): K2 over every T∪D candidate in W of the certified precursor (production isotope set),
    rank-1 by hyperscore (engine rank; labels never joined)."""
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
    import m3_common as C
    p = E.SearchParams()
    l1 = L1Index(idx_dir)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)

    def run_one(run, cert):
        sp = stores[run]
        sub = certified_precursor(sp, cert["spec"], cert["mz_mono"], cert["z"])
        tbl, _, _ = E.search_exhaustive(sub, gi, l1, p, report=1)
        s = tbl.column("spec_idx").to_numpy()
        hs = tbl.column("hyperscore").to_numpy().astype(float)
        rows = tbl.column("pf_row").to_numpy()
        n = len(sub)
        mx = np.full(n, -np.inf)
        top = [None] * n
        par = [None] * n
        o = np.lexsort((-hs, s))
        first = np.r_[True, s[o][1:] != s[o][:-1]]
        for i in o[first]:
            mx[s[i]] = hs[i]
            top[s[i]] = l1.sequence(int(rows[i]))
            par[s[i]] = l1.sequence(int(l1.partner[rows[i]]))
        return {"max_hs": mx, "top_seq": top, "partner_seq": par}
    return run_one


def run_project(acc: str, files: list, prot: dict, cfg: dict) -> dict:  # pragma: no cover (GPU, data)
    """Own search + marker screen; if the project passes, admission.  Returns {"marker", "counts"} (label-free)."""
    import own_search as OS
    from vulcan.gates.homolog import HomologFilter
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))
    import ingest as ING
    stores = {}
    for f in files:
        out, st = ING.conv((f["path"], os.path.join(W, "store_dn2")))
        if not out or "FAILED" in st:
            raise RuntimeError(f"ingest failed {f['path']}: {st}")
        stores[os.path.basename(out)[:-len(".vspec.arrow")]] = out
    idx, src = OS.build_own_index(prot["upid"], cfg)
    runs = {r: OS.search_run(path, idx) for r, path in stores.items()}
    cm = OS.certify_and_markers(runs, idx, src, cfg)
    out = {"marker": cm["marker"], "counts": None}
    d = os.path.join(W, "admission", acc)
    os.makedirs(d, exist_ok=True)
    json.dump(cm["marker"], open(os.path.join(d, "marker.json"), "w"), indent=1)
    if cm["marker"]["exclude"]:
        return out
    gate_idx = cfg_gate_index()
    db, ph = gate_objects(gate_idx)
    sp_by_run = {r: runs[r][0] for r in runs}
    near = HomologFilter(src["own"] + src["contam"]).near
    thr = json.load(open(os.path.join(W, "host_threshold.json")))["tail_proxy_threshold"]
    res = admit_project(cm["certified"], db, ph, gpu_exhaustive(gate_idx, sp_by_run), near, thr, cfg)
    for run, cd in res["codes"].items():
        a = res["A"][run]
        np.savez(os.path.join(d, run + ".npz"), codes=cd, A_spec=a["spec"], A_mz_mono=a["mz_mono"], A_z=a["z"],
                 store=np.array(stores[run]))
    json.dump({"counts": res["counts"], "membership_sha256": res["membership_sha256"], "upid": prot["upid"],
               "runs": sorted(res["codes"])}, open(os.path.join(d, "counts.json"), "w"), indent=1)
    out["counts"] = res["counts"]
    return out


def cfg_gate_index() -> str:  # pragma: no cover
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
    import m3_common as C
    return C.index_dir("human_contam", "prod")


def swap_test(acc: str, cfg: dict) -> dict:  # pragma: no cover (heavy; run at step 3)
    """P1 on real data: rebuild the gate index with T and D swapped (vulcan.index.build.build(swap_labels=True), the G0
    knob) and re-run the whole admission of one already-admitted project on it."""
    import tempfile

    import own_search as OS
    from vulcan.gates.homolog import HomologFilter
    from vulcan.index.build import build
    gate = cfg_gate_index()
    meta = json.load(open(os.path.join(gate, "meta.json")))
    d = tempfile.mkdtemp(prefix="a2_swap_", dir=W)
    fasta = meta["fasta"] if isinstance(meta["fasta"], list) else [meta["fasta"]]
    build(fasta, d, params=meta["params"], swap_labels=True)
    same = {f: a2_config.sha256_file(os.path.join(gate, f)) == a2_config.sha256_file(os.path.join(d, f))
            for f in ("peptidoforms.arrow", "residues.npy")}
    ref = json.load(open(os.path.join(W, "admission", acc, "counts.json")))
    runs = sorted(ref["runs"])
    stores = {r: str(np.load(os.path.join(W, "admission", acc, r + ".npz"))["store"]) for r in runs}
    idx, src = OS.build_own_index(ref["upid"], cfg)
    rr = {r: OS.search_run(stores[r], idx) for r in runs}
    cm = OS.certify_and_markers(rr, idx, src, cfg)
    db, ph = gate_objects(d)
    thr = json.load(open(os.path.join(W, "host_threshold.json")))["tail_proxy_threshold"]
    res = admit_project(cm["certified"], db, ph, gpu_exhaustive(d, {r: rr[r][0] for r in runs}),
                        HomologFilter(src["own"] + src["contam"]).near, thr, cfg)
    return {"project": acc, "label_free_files_identical": same, "membership_sha256_ref": ref["membership_sha256"],
            "membership_sha256_swapped": res["membership_sha256"],
            "pass": bool(all(same.values()) and res["membership_sha256"] == ref["membership_sha256"])}


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--swap-test", metavar="ACCESSION")
    ap.add_argument("--out")
    a = ap.parse_args()
    cfg = a2_config.load()
    if a.swap_test:
        r = swap_test(a.swap_test, cfg)
        json.dump(r, open(a.out or os.path.join(W, "swap_test.json"), "w"), indent=1)
        print("P1 swap test", "PASS" if r["pass"] else "FAIL")


if __name__ == "__main__":
    main()
