#!/usr/bin/env python3
"""Amendment A2.2(a) own-species search + certification, and A2.3 step 3 marker screen (per DN2 project).
Runs inside the seeded walk (scripts/a2/select_dn2.py -> admit.run_project) after the tag prereg-v0.4.2.

Own database: the project's UniProt reference proteome (upid recorded at selection; streamed from UniProt REST
proteome:<upid>) ∪ GPM cRAP ∪ Frankenfield 2022 contaminants ∪ Panel H (incl. keratins) ∪ Panel B, with pseudo-reverse
decoys of all entries (vulcan.index.build.build, default_params()).  First header tokens are rewritten to OWN|<acc>,
CRAP|<tok>, CON|<tok>, and PNL|<panel>|<symbol>|<taxid>|<acc> so that the source of every peptide is known from
proteins.arrow.  meta.json gets "a2_role": "own" (vulcan.gates.a2_admission.own_labels refuses other indexes).
Search: production closed search (E2 SearchParams defaults, report_psms 5), hyperscore; estimator per run = paired
double competition (unit peptidoform, hash tie-break seed 11); PSM-level q = TDC over the run's rank-1 PSMs.
The human gate database is NOT touched here.
"""
from __future__ import annotations

import os
import sys
import urllib.parse
import urllib.request

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import a2_admission as AD  # noqa: E402
from vulcan.gates import a2_config  # noqa: E402

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")


def _records(path):
    hdr, cur = None, []
    for line in open(path):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(cur)
            hdr, cur = line[1:].strip(), []
        else:
            cur.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(cur)


def fetch_proteome(upid: str, out: str) -> str:  # pragma: no cover (network)
    if not os.path.exists(out):
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
            {"query": f"proteome:{upid}", "format": "fasta"})
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "vulcan-a2"}), timeout=600) as r:
            open(out + ".part", "wb").write(r.read())
        os.replace(out + ".part", out)
    return out


def write_own_fasta(own_fasta: str, cfg: dict, out: str) -> dict:
    """Combined, prefixed FASTA.  Returns the sequences per source (for the A1.3-near filters)."""
    p = cfg["panels"]
    src = {"own": [], "contam": []}
    with open(out, "w") as fh:
        for hdr, seq in _records(own_fasta):
            tok = hdr.split()[0]
            acc = tok.split("|")[1] if tok.count("|") >= 2 else tok
            fh.write(f">OWN|{acc} {hdr}\n{seq}\n")
            src["own"].append(seq)
        for path, pre in ((p["crap_fasta"], ""), (p["frankenfield_fasta"], "CON|"), (p["fasta"], "")):
            for hdr, seq in _records(os.path.join(ROOT, path)):
                fh.write(f">{pre}{hdr}\n{seq}\n")
                src["contam"].append(seq)
    return src


def build_own_index(upid: str, cfg: dict) -> tuple[str, dict]:  # pragma: no cover (heavy)
    import json

    from vulcan.index.build import build, default_params
    d = os.path.join(W, "own", upid)
    os.makedirs(d, exist_ok=True)
    own = fetch_proteome(upid, os.path.join(d, f"{upid}.fasta"))
    comb = os.path.join(d, "own_db.fasta")
    src = write_own_fasta(own, cfg, comb)
    idx = os.path.join(d, "index")
    if not os.path.exists(os.path.join(idx, "meta.json")):
        build([comb], idx, params=default_params())
        m = json.load(open(os.path.join(idx, "meta.json")))
        m.update(a2_role="own", upid=upid, own_fasta_sha256=a2_config.sha256_file(own),
                 own_db_sha256=a2_config.sha256_file(comb))
        json.dump(m, open(os.path.join(idx, "meta.json"), "w"), indent=1)
    return idx, src


def search_run(store_path: str, idx: str):  # pragma: no cover (GPU)
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
    import m3_common as C
    p = E.SearchParams(report_psms=C.REPORT)
    l1 = L1Index(idx)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    sp = E.load_store(store_path)
    tbl, _, _ = E.search(sp, gi, l1, p)
    return sp, l1, {c: tbl.column(c).to_numpy() for c in ("spec_idx", "pf_row", "pf_id", "charge", "isotope_error",
                                                         "hyperscore")}


def certify_and_markers(runs: dict, idx: str, src: dict, cfg: dict) -> dict:  # pragma: no cover (needs searches)
    """runs: {run: (Spectra, L1Index, psm dict)} of one project.  -> per-run certified tables + marker screen."""
    from vulcan.gates.homolog import HomologFilter
    from vulcan.gpu.host import NEUTRON
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    lab = AD.own_labels(idx)
    sources = AD.peptide_sources(lab)
    hf_contam = HomologFilter(src["contam"])
    hf_own = HomologFilter(src["own"])
    oc = cfg["own_search"]["certification"]
    seed = cfg["seeds"]["tiebreak"]
    cert, accepted = {}, {}
    for run, (sp, l1, psm) in runs.items():
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        seq = {int(r): l1.sequence(int(r)) for r in np.unique(rows)}
        il = {r: s.replace("I", "L") for r, s in seq.items()}
        il_by_row = np.array([il[int(r)] for r in rows], dtype=object)
        tb = hash_tiebreak(psm["pf_id"], seed)
        res = AD.certify_run(psm["spec_idx"], psm["pf_id"], lab["tdc_pair"][rows], psm["hyperscore"], dec,
                             il_by_row, sources, hf_contam.near, tb, oc["q_psm"], oc["q_peptide"])
        w = res["rows"]
        z = psm["charge"][w].astype(np.int64)
        iso = psm["isotope_error"][w].astype(np.int64)
        s_idx = psm["spec_idx"][w].astype(np.int64)
        mz_mono = sp.precursor_mz[s_idx] - iso * float(NEUTRON) / z
        r_ = rows[w]
        cert[run] = {"spec": s_idx, "p_seq": np.array([seq[int(r)] for r in r_], dtype=object),
                     "p_mask": l1.mask[r_].astype(np.uint64), "z": z, "iso": iso, "mz_mono": mz_mono,
                     "hyperscore_p": psm["hyperscore"][w].astype(np.float64), "n_spectra": len(sp)}
        dc = double_competition(psm["spec_idx"], psm["pf_id"], lab["tdc_pair"][rows], psm["hyperscore"], dec,
                                tiebreak=tb)
        ok = (dc.q <= cfg["marker_screen"]["q_peptide"]) & ~dc.is_decoy
        accepted[run] = {il[int(rows[i])] for i in dc.row[ok]}
    mk = AD.marker_screen(accepted, sources, hf_own.near, tuple(cfg["marker_screen"]["trigger_panels"]),
                          cfg["marker_screen"]["min_distinct_peptides"],
                          cfg["marker_screen"]["exclude_if_trigger_members_at_least"])
    return {"certified": cert, "marker": mk}
