#!/usr/bin/env python3
"""A6 drafting check (label-free on null candidates; DEV data only): yield of the FINAL A6 admission rules on the
C9 pilot's 12 DEV runs (D0 2 runs, D8h 10 runs).

Differences from scripts/c9/pilot_yield.py (which is imported and reused unchanged):
  * ladder screen f_max = 0.05 (A6.3.2(ii)), the only primary value;
  * contaminant pairs dropped from the gate DB (A6.3.2(iii)): a T∪D row is removed iff its target peptide (or the
    pair partner's) occurs in any protein of the contaminant FASTA (label-free: sequence -> source, pair-closed);
  * certification variants (A6.3.1), counted on the SAME null stage (the null stage does not depend on which
    certified spectra are kept, except through I, which is unchanged):
      cert_hs     hyperscore, per-run paired double competition, PSM q <= 1e-3, peptidoform q <= 1e-2 (C9 pilot rule);
      cert_agree  cert_hs AND a second engine (Sage, frozen comparator binary) reports the same I/L-collapsed stripped
                  peptide at the same charge as its rank-1 target PSM with spectrum_q <= 1e-3.  The Sage searches
                  used here are the existing DEV arm-C runs (human + contaminant + Arabidopsis entrapment, decoy seed
                  1001; $VULCAN_T1_WORK/{D0,D8h}/cmp/sage/C_1001), a conservative proxy (larger DB) for the
                  A6 rule's human + contaminant DB;
      cert_pool4  project-pooled hyperscore TDC at PSM q <= 1e-4 (runs of a project concatenated; only D0 has 2 runs).
Labels of the full index are read ONLY for certification, the host threshold and I (exactly as the pilot).  No
target/decoy label of any null candidate is read; every null quantity is label-collapsed (counts, max over T∪D).

    conda_envs/pe/bin/python reviews/checks/a6_cert_yield.py   # -> reviews/checks/a6_cert_yield.json
"""
from __future__ import annotations

import csv
import json
import os
T1_WORK = os.environ.get("VULCAN_T1_WORK", f"/tmp/vulcan_t1_{os.getuid()}")  # bench/t1_driver.py default work root
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "c9"))
sys.path.insert(0, ROOT)
import pilot_yield as PY  # noqa: E402

OUT = os.path.join(HERE, "a6_cert_yield.json")
SAGE = {"D0": T1_WORK + "/D0/cmp/sage/C_1001/results.sage.tsv",
        "D8h": T1_WORK + "/D8h/cmp/sage/C_1001/results.sage.tsv"}
F_MAX = 0.05
SAGE_Q = 1e-3


def sage_rank1(path):
    """(mzML stem, scan) -> (I/L stripped peptide, charge) for rank-1 TARGET PSMs with spectrum_q <= SAGE_Q."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["rank"] != "1" or r["label"] != "1" or float(r["spectrum_q"]) > SAGE_Q:
                continue
            m = re.search(r"scan=(\d+)", r["scannr"])
            pep = re.sub(r"[^A-Z]", "", r["peptide"]).replace("I", "L")
            out[(r["filename"][:-len(".mzML")], int(m.group(1)))] = (pep, int(r["charge"]))
    return out


class NearPlusContam:
    """HomologFilter(I).near OR 'sequence occurs in a contaminant-FASTA protein' (pair closure is done by the
    pilot's null stage, which flags a row if the row OR its partner is flagged)."""

    def __init__(self, hf, contam_seqs):
        self.hf, self.cs = hf, contam_seqs

    def near(self, seqs):
        a = np.asarray(self.hf.near(seqs), bool)
        b = np.array([s.replace("I", "L") in self.cs for s in seqs], bool)
        return a | b


def main():
    import m3_common as C
    from vulcan.gates import a2_admission as AD
    from vulcan.gates.a2_identity import isobaric_variants
    from vulcan.gates.homolog import HomologFilter, read_fasta_seqs
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.io import store
    from vulcan.stats._util import hash_tiebreak
    PY.F_THR = F_MAX
    t0 = time.perf_counter()
    l1 = L1Index(PY.IDX)
    p = E.SearchParams(report_psms=C.REPORT)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    lab, src = PY.full_labels(PY.IDX)
    near_contam = HomologFilter(read_fasta_seqs([PY.CONTAM_FASTA])).near
    contam_seqs = {k.replace("I", "L") for k, v in src.items() if any(a.startswith("CON|") for a in v)}
    sage = {k: sage_rank1(v) for k, v in SAGE.items()}
    e9 = {r["id"]: r["accession"] for r in csv.DictReader(open(os.path.join(ROOT, "data/manifests/e9_d8.tsv")),
                                                           delimiter="\t")}
    projects = [("PXD001468", "D0", [(os.path.basename(x)[:-len(".vspec.arrow")], x) for x in PY.D0_STORES])]
    for x in PY.D8H_STORES:
        run = os.path.basename(x)[:-len(".vspec.arrow")]
        projects.append((e9[run.split("__")[0]], "D8h", [(run, x)]))
    res = {"rules": {"f_max": F_MAX, "contaminant_pairs_dropped": True, "sage_spectrum_q": SAGE_Q,
                     "n_contaminant_target_peptides": len(contam_seqs)}, "per_run": {}, "per_project_pool4": {}}
    for name, ds, runs in projects:
        per, ident_proj, pooled = {}, set(), []
        for run, path in runs:
            clk = PY.Clock()
            sp = E.load_store(path)
            cert, ident, thr, info = PY.certify(run, sp, gi, l1, p, lab, src, near_contam, clk)
            ident_proj |= ident
            t = store.read_table(path)
            lvl = t.column("ms_level").to_numpy()
            scans = t.column("scan").to_numpy()[lvl == 2].astype(np.int64)
            agree = np.zeros(len(cert["spec"]), bool)
            for i, (s, r, z) in enumerate(zip(cert["spec"], cert["row"], cert["z"])):
                h = sage[ds].get((run, int(scans[int(s)])))
                agree[i] = h is not None and h[0] == l1.sequence(int(r)).replace("I", "L") and h[1] == int(z)
            per[run] = (sp, cert, thr, info, agree)
            # project-pooled PSM q <= 1e-4 (hyperscore TDC over the concatenated runs of the project)
            tbl, _, _ = E.search(sp, gi, l1, p)
            pooled.append((tbl, len(sp)))
        off, cols = 0, {"spec": [], "rows": [], "pf": [], "hs": []}
        for tbl, n in pooled:
            cols["spec"].append(tbl.column("spec_idx").to_numpy().astype(np.int64) + off)
            cols["rows"].append(tbl.column("pf_row").to_numpy().astype(np.int64))
            cols["pf"].append(tbl.column("pf_id").to_numpy())
            cols["hs"].append(tbl.column("hyperscore").to_numpy().astype(float))
            off += n
        spec, rows, pf, hs = (np.concatenate(cols[k]) for k in ("spec", "rows", "pf", "hs"))
        useq = {int(r): l1.sequence(int(r)) for r in np.unique(rows)}
        il = np.array([useq[int(r)].replace("I", "L") for r in rows], dtype=object)
        r4 = AD.certify_run(spec, pf, lab["tdc_pair"][rows], hs, lab["is_decoy"][rows], il, src, near_contam,
                            hash_tiebreak(pf, PY.TB_SEED), 1e-4, PY.Q_PEP)
        res["per_project_pool4"][name] = {"n_runs": len(runs), "n_certified_pool_q1e-4": r4["n_certified"]}
        extra = []
        for q in ident_proj:
            v = isobaric_variants(q, 64)
            if v:
                extra += [x for x in v if x != q]
        hf = NearPlusContam(HomologFilter(sorted(ident_proj) + extra), contam_seqs)
        for run, (sp, cert, thr, info, agree) in per.items():
            clk = PY.Clock()
            mx, _, _, rinfo = PY.null_stage(sp, gi, l1, p, cert, hf, clk)
            cand = np.isfinite(mx)
            tail = mx >= thr
            rec = {"project": name, "run": run, "n_ms2": info["n_ms2"],
                   "cert_hs": {"n_certified": int(len(cert["spec"])), "A_cand": int(cand.sum()),
                               "tail_proxy_own_thr": int(tail.sum())},
                   "cert_agree": {"n_certified": int(agree.sum()), "A_cand": int((cand & agree).sum()),
                                  "tail_proxy_own_thr": int((tail & agree).sum())},
                   "removed_rows_homolog_or_contam": rinfo["removed_rows_homolog"],
                   "removed_rows_ladder": rinfo["removed_rows_ladder"]}
            res["per_run"][run] = rec
            print(json.dumps(rec), flush=True)
    tot = {v: {k: int(sum(r[v][k] for r in res["per_run"].values())) for k in ("n_certified", "A_cand",
                                                                               "tail_proxy_own_thr")}
           for v in ("cert_hs", "cert_agree")}
    res["totals"] = tot
    res["retention_agree_vs_hs"] = {k: tot["cert_agree"][k] / max(tot["cert_hs"][k], 1) for k in tot["cert_hs"]}
    res["wall_s"] = time.perf_counter() - t0
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    print(json.dumps({"totals": tot, "retention": res["retention_agree_vs_hs"],
                      "pool4": res["per_project_pool4"]}), flush=True)


if __name__ == "__main__":
    main()
