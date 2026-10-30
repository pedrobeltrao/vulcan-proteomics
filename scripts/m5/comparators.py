#!/usr/bin/env python3
"""M5 comparators: true FDP (library design, scripts/m5/truth.py) of third-party XL results on the same dev files.

* submitter xiSEARCH 1.7.6.5 + xiFDR 2.1.5.5 result for D3b (PRIDE PXD029252 `..._KK_1FDR_linkandCSM_Xi1.7.6.5.mzid`,
  data/truth/PXD029252; 1% link + CSM FDR): passing cross-link CSMs -> CSM / peptide pair / residue pair / PPI sets.
* our own xiSEARCH/xiFDR runs (normalized TSVs from /tmp/vulcan_m5/comparators), if present.

    python scripts/m5/comparators.py --out RES/comparators.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from truth import Truth  # noqa: E402

NS = "{http://psidev.info/psi/pi/mzIdentML/1.2}"
MZID_D3B = os.path.join(ROOT, "data/truth/PXD029252/"
                        "20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2_KK_1FDR_linkandCSM_Xi1.7.6.5.mzid")


def parse_xi_mzid(path):
    """Passing cross-link CSMs: list of dict(scan, pep_a, pep_b, site_a, site_b (0-based), prot_a, prot_b, pos_a,
    pos_b (0-based protein offsets of the peptides), decoy_a, decoy_b, score)."""
    root = ET.parse(path).getroot()
    dbseq = {}
    for d in root.iter(NS + "DBSequence"):
        s = d.find(NS + "Seq")
        dbseq[d.get("id")] = (d.get("accession"), s.text if s is not None else "")
    peps = {}
    for p in root.iter(NS + "Peptide"):
        seq = p.find(NS + "PeptideSequence").text
        loc = None
        for m in p.findall(NS + "Modification"):
            for cv in m.findall(NS + "cvParam"):
                if cv.get("accession") in ("MS:1002509", "MS:1002510"):
                    loc = int(m.get("location"))
        peps[p.get("id")] = (seq, loc)
    pev = {}
    for e in root.iter(NS + "PeptideEvidence"):
        acc, _ = dbseq.get(e.get("dBSequence_ref"), ("?", ""))
        pev[e.get("id")] = (acc, e.get("isDecoy") == "true", int(e.get("start") or 0) - 1)
    out = []
    for sir in root.iter(NS + "SpectrumIdentificationResult"):
        scan = None
        for cv in sir.findall(NS + "cvParam"):
            if cv.get("accession") == "MS:1000797":
                scan = int(cv.get("value"))
        groups = {}
        for sii in sir.findall(NS + "SpectrumIdentificationItem"):
            if sii.get("passThreshold") != "true":
                continue
            xl = score = None
            for cv in sii.findall(NS + "cvParam"):
                if cv.get("accession") == "MS:1002511":
                    xl = cv.get("value")
                if cv.get("accession") == "MS:1002545":
                    score = float(cv.get("value"))
            if xl is None:
                continue
            ref = sii.find(NS + "PeptideEvidenceRef").get("peptideEvidence_ref")
            groups.setdefault(xl, []).append((peps[sii.get("peptide_ref")], pev[ref], score))
        for g in groups.values():
            if len(g) != 2:
                continue
            (sa, la), (pa, da, oa), sc = g[0]
            (sb, lb), (pb, db, ob), _ = g[1]
            out.append({"scan": scan, "pep_a": sa, "pep_b": sb, "site_a": (la or 1) - 1, "site_b": (lb or 1) - 1,
                        "prot_a": pa, "prot_b": pb, "pos_a": oa, "pos_b": ob, "decoy_a": da, "decoy_b": db,
                        "score": sc})
    return out


def _acc(x: str) -> str:
    """'sg|P0A7N1|RL31B_ECOLI' / 'P0A7N1' -> 'P0A7N1'."""
    return x.split("|")[1] if "|" in x else x


def true_fdp_of_sets(which, rows, prot_seqs):
    """rows: accepted target-target CSMs (one per spectrum). Returns per level: n, false, true_fdp, and self/hetero
    split (self = same protein)."""
    truth = Truth(which, prot_seqs)
    seqs = {_acc(a): s.replace("I", "L") for a, s in prot_seqs.items()}
    rg = {}
    for p, g in truth.lib:
        for a, s in seqs.items():
            i = s.find(p)
            while i >= 0:
                for k in range(i, i + len(p)):
                    rg.setdefault((a, k), set()).add(g)
                i = s.find(p, i + 1)
    pg = {}
    for p, g in truth.lib:
        for a, s in seqs.items():
            if p in s:
                pg.setdefault(a, set()).add(g)
    res = {}
    lv = {"csm": {}, "peptide_pair": {}, "residue_pair": {}, "ppi": {}}
    for r in rows:
        pa, pb = _acc(r["prot_a"]), _acc(r["prot_b"])
        selfl = pa == pb
        ra = (pa, r["pos_a"] + r["site_a"])
        rb = (pb, r["pos_b"] + r["site_b"])
        keys = {"csm": r["scan"], "peptide_pair": tuple(sorted((r["pep_a"], r["pep_b"]))),
                "residue_pair": tuple(sorted((ra, rb))), "ppi": tuple(sorted((pa, pb)))}
        tru = {"csm": truth.link_class(r["pep_a"], r["pep_b"]) == "true",
               "peptide_pair": truth.link_class(r["pep_a"], r["pep_b"]) == "true",
               "residue_pair": bool(rg.get(ra, set()) & rg.get(rb, set())),
               "ppi": bool(pg.get(pa, set()) & pg.get(pb, set()))}
        for k in lv:
            lv[k].setdefault(keys[k], (tru[k], selfl))
    for k, d in lv.items():
        vals = list(d.values())
        cell = {}
        for name, sel in (("all", lambda s: True), ("self", lambda s: s), ("hetero", lambda s: not s)):
            v = [t for t, s in vals if sel(s)]
            n, nf = len(v), sum(1 for t in v if not t)
            cell[name] = {"n": n, "false": nf, "true_fdp": (nf / n) if n else None}
        res[k] = cell
    return res


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--comp-dir", default="/tmp/vulcan_m5/comparators")
    a = ap.parse_args()
    res = {}
    fa = {"D3b": "/tmp/vulcan_m5/fasta/D3b_ribo_human671.fasta", "D3": "/tmp/vulcan_m5/fasta/D3_cas9_contam.fasta"}
    rows = parse_xi_mzid(MZID_D3B)
    tt = [r for r in rows if not r["decoy_a"] and not r["decoy_b"]]
    res["submitter_xisearch_1.7.6.5_xifdr_D3b"] = {
        "source": os.path.relpath(MZID_D3B, ROOT), "fdr_setting": "xiFDR 1% link + CSM (file name; PRIDE)",
        "n_passing_xl_csms": len(rows), "n_tt": len(tt), "n_with_decoy": len(rows) - len(tt),
        "true_fdp": true_fdp_of_sets("D3b", tt, read_fasta(fa["D3b"])),
        "true_fdp_non_homo": true_fdp_of_sets("D3b", [r for r in tt if r["pep_a"].replace("I", "L") !=
                                                      r["pep_b"].replace("I", "L")], read_fasta(fa["D3b"])),
        "note": "submitter's search DB (proteins in the mzid) may differ from our 671-sequence DB; truth is "
                "sequence-based (library groups), so it applies regardless."}
    for f in sorted(glob.glob(os.path.join(a.comp_dir, "*_csm.tsv"))):
        import pandas as pd
        name = os.path.basename(f)[:-8]
        which = "D3b" if name.endswith("D3b") else "D3"
        t = pd.read_csv(f, sep="\t")
        if "kind" in t.columns:
            t = t[t["kind"] == "crosslink"]
        t = t[t.site_a.notna() & t.site_b.notna()]
        for lvl_col in [c for c in t.columns if c.startswith("passes_")]:
            sel = t[(t[lvl_col].astype(str).str.lower() == "true") & ~t.decoy_a.astype(bool) & ~t.decoy_b.astype(bool)]
            pos = {}
            seqs = read_fasta(fa[which])
            rr = []
            for _, r in sel.iterrows():
                pa_, pb_ = str(r.prot_a).split(";")[0], str(r.prot_b).split(";")[0]
                sa = next((s for k, s in seqs.items() if _acc(k) == _acc(pa_)), "")
                sb = next((s for k, s in seqs.items() if _acc(k) == _acc(pb_)), "")
                rr.append({"scan": r.scan, "pep_a": r.pep_a, "pep_b": r.pep_b, "site_a": int(r.site_a) - 1,
                           "site_b": int(r.site_b) - 1, "prot_a": pa_, "prot_b": pb_,
                           "pos_a": sa.replace("I", "L").find(str(r.pep_a).replace("I", "L")),
                           "pos_b": sb.replace("I", "L").find(str(r.pep_b).replace("I", "L"))})
            nh = [x for x in rr if x["pep_a"].replace("I", "L") != x["pep_b"].replace("I", "L")]
            res[f"{name}:{lvl_col}"] = {"n_tt": len(rr), "true_fdp": true_fdp_of_sets(which, rr, seqs),
                                        "n_tt_non_homo": len(nh),
                                        "true_fdp_non_homo": true_fdp_of_sets(which, nh, seqs)}
            del pos
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    for k, v in res.items():
        print(k, json.dumps({kk: vv for kk, vv in v.items() if kk in ("n_passing_xl_csms", "n_tt")}),
              json.dumps(v["true_fdp"]))


if __name__ == "__main__":
    main()
