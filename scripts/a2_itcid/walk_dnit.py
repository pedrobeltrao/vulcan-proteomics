#!/usr/bin/env python3
"""DNIT (IT-CID fresh nulls) steps 2-3: seeded walk, PRIDE-SHA-1-verified download, IT-CID class check, own-species
search + certification + marker screen, label-free admission (A2.2), stopping rule.  Rules (normative):
results/2026-10-27_3c_itcid/followup/NULL_RULES.md, config/a2_itcid.yaml (committed 2814f06 before any listing).

Tagged A2 code is imported, never edited:  scripts/a2/select_dn2.py (frame filter, protocol / SDRF rule, file pick,
reference proteome, PRIDE I/O), scripts/a2/own_search.py (own DB, certification, marker screen), scripts/a2/admit.py
(admit_project), vulcan.gates.a2_* .  IT-CID adaptations live here only:
  * search parameters = 3c production IT-CID (fragment +-0.5 Da) for the own search and the (d) exhaustive search;
  * (b) f(c): ions coincide within +-0.5 Da (max_f_in_window_da below; replaces admit.max_f_in_window);
  * gate index = human_contam prod; tail proxy threshold = max over the 29 D10cf hosts;
  * the IT-CID class check after ingest; downloads verified by PRIDE SHA-1 (project skipped if unverifiable);
  * a run with zero certified spectra returns empty exhaustive arrays (A2 DN2 deviation item 17, same guard).
No per-label quantity of any null spectrum is computed, printed or stored here (labels are read only for the 29 host
runs in host_threshold(), as scripts/a2/host_threshold.py does for D0).

    python scripts/a2_itcid/walk_dnit.py host-threshold
    python scripts/a2_itcid/walk_dnit.py walk
    python scripts/a2_itcid/walk_dnit.py swap-test ACC
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

W = os.environ.setdefault("VULCAN_A2_WORK", "/tmp/vulcan_dnit")       # own indexes etc. (never DN2's /tmp/vulcan_a2)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (ROOT, os.path.join(ROOT, "scripts", "a2"), os.path.join(ROOT, "scripts", "m3c"),
          os.path.join(ROOT, "scripts", "m3"), os.path.join(ROOT, "scripts", "m3e")):
    sys.path.insert(0, p)

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import select_dn2 as S  # noqa: E402
from vulcan.gates import a2_config  # noqa: E402
from vulcan.gates import a2_identity as AI  # noqa: E402
from vulcan.gates.a2_admission import stopping_state  # noqa: E402

MAN = os.path.join(ROOT, "data", "manifests")
RAW = os.path.join(W, "raw")
STORE = os.path.join(W, "store")
EVENTS = os.path.join(W, "walk_events.jsonl")
CFG_IT = os.path.join(ROOT, "config", "a2_itcid.yaml")
GATE_INDEX = "/tmp/vulcan_e2/index/human_contam_prod"
SEARCH_IT = {"report_psms": 5, "fragment_unit": "da", "fragment_ppm": (-0.5, 0.5), "precursor_ppm": (-10.0, 10.0),
             "bin_ppm": 200.0}
FRAG_DA = 0.5


def event(**kw):
    kw = {"t": time.strftime("%Y-%m-%dT%H:%M:%S"), **kw}
    os.makedirs(W, exist_ok=True)
    with open(EVENTS, "a") as fh:
        fh.write(json.dumps(kw, default=str) + "\n")
    print(json.dumps(kw, default=str), flush=True)


# ------------------------------------------------------------------------------------------------ config
def load_cfg():
    base = a2_config.load()
    it = yaml.safe_load(open(CFG_IT))
    cfg = copy.deepcopy(base)
    cfg["seeds"].update({k: it["seeds"][k] for k in ("walk", "file_pick_base", "mixture_permutation",
                                                      "offsets_a1_recipe")})
    cfg["selection"]["allow"] = it["selection"]["allow_it"]
    cfg["selection"]["frame_filter"]["instruments_regex"] = it["selection"]["frame_filter"]["instruments_regex"]
    cfg["selection"]["raw_size_bytes"] = it["selection"]["raw_size_bytes"]
    sp = it["selection"]["exclusions"]["species"]
    cfg["selection"]["exclusions"]["training_species"] = sorted(set(cfg["selection"]["exclusions"]["training_species"])
                                                                | set(sp))
    return base, it, cfg


def exclusions(base, frame):
    """A2.9.1 set (manifests + adoption list) ∪ DN2-seen ∪ DN2-reachable (frame rows passing DN2's frame filter)."""
    excl = S.load_exclusions(base)
    seen = set()
    for f in [os.path.join(MAN, x) for x in os.listdir(MAN) if x.startswith("dn2")] + \
            [os.path.join(ROOT, "results", "2026-10-27_A2_dn2", f) for f in ("selection_log.tsv", "walk_events.jsonl")]:
        if os.path.exists(f):
            seen |= set(re.findall(r"PXD\d{6}", open(f).read()))
    reach = {r["accession"] for r in S.frame_candidates(frame, base, set())[0]}
    return excl | seen | reach, {"a2_manifest_union": len(excl), "dn2_seen": len(seen), "dn2_reachable": len(reach)}


# ------------------------------------------------------------------------------------------------ IT-CID admission
def f_overlap_da(p_ions, c_ions, tol=FRAG_DA):
    if p_ions.size == 0:
        return 0.0
    s = np.sort(c_ions)
    lo = np.searchsorted(s, p_ions - tol, side="left")
    hi = np.searchsorted(s, p_ions + tol, side="right")
    return float(np.count_nonzero(hi > lo)) / p_ions.size


def max_f_in_window_da(db, p_seq, p_mask, mz_mono, z, cache=None, ppm_frag=None, ppm_prec=AI.PREC_PPM,
                       isotopes=AI.ISOTOPE_ERRORS, static=AI.STATIC, var=AI.VAR):
    """vulcan.gates.a2_identity.max_f_in_window with the IT-CID fragment tolerance (+-0.5 Da; ppm_frag ignored)."""
    p_ions = AI.ion_ladder(p_seq, p_mask, first=2, static=static, var=var)
    best = 0.0
    for r in db.window_rows(mz_mono, z, ppm_prec, isotopes):
        key = (db.seq(r), int(db.mask[r]))
        if cache is not None and key in cache:
            c_ions = cache[key]
        else:
            c_ions = AI.ion_ladder(key[0], key[1], first=1, static=static, var=var, nterm=db.nterm)
            if cache is not None:
                cache[key] = c_ions
        f = f_overlap_da(p_ions, c_ions)
        if f > best:
            best = f
            if best >= 1.0:
                break
    return best


def search_run_it(store_path, idx):
    import common as C3
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(**SEARCH_IT)
    l1 = L1Index(idx)
    st, var = C3.mods("human")
    gi = E.GpuIndex(l1, p, st, var)
    sp = E.load_store(store_path)
    tbl, _, _ = E.search(sp, gi, l1, p)
    del gi
    _free()
    return sp, l1, {c: tbl.column(c).to_numpy() for c in ("spec_idx", "pf_row", "pf_id", "charge", "isotope_error",
                                                         "hyperscore")}


def _free():
    import torch
    torch.cuda.empty_cache()


def gpu_exhaustive_it(idx_dir, stores):
    """admit.gpu_exhaustive with the IT-CID search parameters (+ the zero-certified-run guard)."""
    import common as C3
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(**SEARCH_IT)
    l1 = L1Index(idx_dir)
    st, var = C3.mods("human")
    gi = E.GpuIndex(l1, p, st, var)

    def run_one(run, cert):
        if len(cert["spec"]) == 0:
            return {"max_hs": np.zeros(0), "top_seq": [], "partner_seq": []}
        sp = stores[run]
        sub = certified_precursor(sp, cert["spec"], cert["mz_mono"], cert["z"])
        tbl, _, _ = E.search_exhaustive(sub, gi, l1, p, report=1)
        s = tbl.column("spec_idx").to_numpy()
        hs = tbl.column("hyperscore").to_numpy().astype(float)
        rows = tbl.column("pf_row").to_numpy()
        n = len(sub)
        mx = np.full(n, -np.inf)
        top, par = [None] * n, [None] * n
        o = np.lexsort((-hs, s))
        first = np.r_[True, s[o][1:] != s[o][:-1]] if len(o) else np.zeros(0, bool)
        for i in o[first]:
            mx[s[i]] = hs[i]
            top[s[i]] = l1.sequence(int(rows[i]))
            par[s[i]] = l1.sequence(int(l1.partner[rows[i]]))
        return {"max_hs": mx, "top_seq": top, "partner_seq": par}
    return run_one, gi


def class_check(store):
    import fetch_train_itcid as F
    return F.validate({"stratum": "LF"}, store)


def ingest(path):
    import ingest as ING
    out, st = ING.conv((path, STORE))
    if not out or "FAILED" in st:
        raise RuntimeError(f"ingest failed {path}: {st}")
    return out


_GATE = {}


def gate_objects(idx_dir):
    if idx_dir not in _GATE:
        db = AI.GateDB.from_l1(idx_dir)
        _GATE[idx_dir] = (db, AI.PeptideHamming(db.peptides()))
    return _GATE[idx_dir]


def run_project(acc, files, prot, cfg, gate_index=GATE_INDEX, swap=False):
    """Ingest + class check + own search + marker screen + admission of one project.  -> dict (label-free)."""
    import admit as ADM
    import own_search as OS
    from vulcan.gates.homolog import HomologFilter
    stores = {}
    for f in files:
        s = ingest(f["path"])
        stores[os.path.basename(s)[:-len(".vspec.arrow")]] = s
    cc = {r: class_check(s) for r, s in stores.items()}
    if not all(v["accept"] for v in cc.values()):
        return {"class_check": cc, "marker": None, "counts": None}
    idx, src = OS.build_own_index(prot["upid"], cfg)
    runs = {r: search_run_it(path, idx) for r, path in stores.items()}
    cm = OS.certify_and_markers(runs, idx, src, cfg)
    out = {"class_check": cc, "marker": cm["marker"], "counts": None}
    d = os.path.join(W, "admission" + ("_swap" if swap else ""), acc)
    os.makedirs(d, exist_ok=True)
    json.dump(cm["marker"], open(os.path.join(d, "marker.json"), "w"), indent=1)
    if cm["marker"]["exclude"]:
        return out
    db, ph = gate_objects(gate_index)
    near = HomologFilter(src["own"] + src["contam"]).near
    thr = json.load(open(os.path.join(W, "host_threshold.json")))["tail_proxy_threshold"]
    ADM.max_f_in_window = max_f_in_window_da                       # IT-CID (b) tolerance
    exh, gi = gpu_exhaustive_it(gate_index, {r: runs[r][0] for r in runs})
    res = ADM.admit_project(cm["certified"], db, ph, exh, near, thr, cfg)
    del gi
    _free()
    for run, cd in res["codes"].items():
        a = res["A"][run]
        np.savez(os.path.join(d, run + ".npz"), codes=cd, A_spec=a["spec"], A_mz_mono=a["mz_mono"], A_z=a["z"],
                 store=np.array(stores[run]))
    json.dump({"counts": res["counts"], "membership_sha256": res["membership_sha256"], "upid": prot["upid"],
               "runs": sorted(res["codes"])}, open(os.path.join(d, "counts.json"), "w"), indent=1)
    out["counts"] = res["counts"]
    out["membership_sha256"] = res["membership_sha256"]
    return out


# ------------------------------------------------------------------------------------------------ PRIDE I/O
def sha_both(p):
    h256, h1 = hashlib.sha256(), hashlib.sha1()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h256.update(b)
            h1.update(b)
    return h256.hexdigest(), h1.hexdigest()


class IO(S.PrideIO):
    _cks: dict = {}

    def _checksum_txt(self, acc, url):
        """PRIDE's per-project checksum.txt (SHA-1 written by PRIDE), as scripts/a2_ops/run_dn2_walk.py (NULL_RULES
        dated note 2)."""
        if acc not in self._cks:
            m = {}
            try:
                txt = self._get(url.rsplit("/", 1)[0] + "/checksum.txt").decode(errors="replace")
                for line in txt.splitlines():
                    parts = line.split()
                    if len(parts) >= 2 and len(parts[-1]) == 40:
                        m[parts[0]] = parts[-1].lower()
                    elif len(parts) >= 2 and len(parts[0]) == 40:
                        m[parts[-1]] = parts[0].lower()
            except Exception:  # noqa: BLE001
                pass
            self._cks[acc] = m
        return self._cks[acc]

    def list_files(self, acc):
        out = [{"fileName": f["fileName"], "fileSizeBytes": f["fileSizeBytes"], "url": self._url(f),
                "pride_sha1": (f.get("checksum") or "").strip().lower()}
               for f in self._files(acc) if f.get("fileCategory", {}).get("value") == "RAW"]
        event(ev="listed", accession=acc, n_raw_records=len(out))
        return out

    def download(self, acc, files):
        """-> rows, or None if a file cannot be verified by PRIDE SHA-1 (project skipped, NULL_RULES §5)."""
        d = os.path.join(RAW, acc)
        os.makedirs(d, exist_ok=True)

        def one(f):
            p = os.path.join(d, f["fileName"])
            exp1 = f.get("pride_sha1", "") or self._checksum_txt(acc, f["url"]).get(f["fileName"], "")  # note 2
            if not exp1:
                return None, "no PRIDE checksum (API or checksum.txt)"
            for attempt in range(3):
                if not (os.path.exists(p) and os.path.getsize(p) == f["fileSizeBytes"]):
                    if attempt > 0 and os.path.exists(p):
                        os.remove(p)
                    subprocess.run(["curl", "-sSf", "--retry", "5", "--speed-limit", "20000", "--speed-time", "120",
                                    "-C", "-", "-o", p, f["url"]], check=False)
                if not os.path.exists(p) or os.path.getsize(p) != f["fileSizeBytes"]:
                    continue
                s256, s1 = sha_both(p)
                if s1 == exp1:
                    return {"filename": f["fileName"], "url": f["url"], "expected_size": f["fileSizeBytes"],
                            "sha256": s256, "sha1": s1, "pride_sha1": exp1, "path": p}, "ok"
                event(ev="sha1_mismatch", accession=acc, file=f["fileName"], attempt=attempt)
                os.remove(p)
            return None, "unverifiable after 3 attempts"
        with ThreadPoolExecutor(3) as ex:
            res = list(ex.map(one, files))
        for (row, why), f in zip(res, files):
            event(ev="downloaded" if row else "download_failed", accession=acc, file=f["fileName"], status=why,
                  **({k: row[k] for k in ("sha256", "sha1")} if row else {}))
        if any(r is None for r, _ in res):
            return None
        return [r for r, _ in res]


def walk(cands, cfg, io, max_projects=None):
    taken, admitted, rows, log = set(), [], [], []
    state = stopping_state([], cfg["selection"]["stopping"])

    def L(**kw):
        log.append(kw)
        event(ev="decision", **kw)
    for r in cands:
        acc, org = r["accession"], r["organisms"].strip()
        slot = S.genus_slot(org, cfg)
        if slot in taken:
            L(accession=acc, organism=org, result="genus slot taken", slot=slot)
            continue
        why = S.protocol_reject(io.project_record(acc), io.sdrf(acc), cfg)
        if why:
            L(accession=acc, organism=org, result="protocol", reason=why)
            continue
        raws = S.usable_raws(io.list_files(acc), cfg)
        if len(raws) < cfg["selection"]["runs_per_project"]:
            L(accession=acc, organism=org, result=f"{len(raws)} usable raw")
            continue
        pick = S.pick_files(raws, acc, cfg)
        prot = io.proteome(r)
        if prot is None:
            L(accession=acc, organism=org, result="no reference proteome")
            continue
        got = io.download(acc, pick)
        if got is None:
            L(accession=acc, organism=org, result="download not verifiable by PRIDE SHA-1")
            continue
        t0 = time.time()
        res = run_project(acc, got, prot, cfg)
        if res["marker"] is None:
            L(accession=acc, organism=org, result="IT-CID class check failed",
              class_check={k: {kk: v[kk] for kk in ("n_ms2", "n_ms3", "it_cid_fraction", "reason")}
                           for k, v in res["class_check"].items()})
            continue
        mk = res["marker"]
        if mk["exclude"]:
            L(accession=acc, organism=org, result="marker screen exclude", upid=prot["upid"],
              trigger_members=mk["trigger_members"], n_trigger_members_detected=mk["n_trigger_members_detected"],
              n_panel_members_detected=mk["n_panel_members_detected"])
            continue
        counts = res["counts"]
        taken.add(slot)
        admitted.append(counts)
        for f in got:
            rows.append({"accession": acc, "organism": org, "upid": prot["upid"], **f})
        state = stopping_state(admitted, cfg["selection"]["stopping"])
        L(accession=acc, organism=org, result="admitted", slot=slot, upid=prot["upid"], counts=counts,
          membership_sha256=res["membership_sha256"], seconds=round(time.time() - t0),
          class_check={k: v["it_cid_fraction"] for k, v in res["class_check"].items()},
          marker={k: mk[k] for k in ("n_panel_members_detected", "n_trigger_members_detected", "detected_members")},
          stopping=state)
        if state["stop"]:
            break
    return {"admitted": [x["accession"] for x in log if x["result"] == "admitted"], "rows": rows,
            "final_state": state, "log": log, "exhausted": not state["stop"]}


# ------------------------------------------------------------------------------------------------ commands
def host_threshold():
    """Max over the 29 D10cf hosts of the host hyperscore 1% threshold (3c production IT-CID search of the host alone,
    $VULCAN_3C_WORK/psms/D10cf_prod; paired double competition, tie-break 11; vulcan.gates.stats.tail_threshold)."""
    import common as C3
    import m3_data as MD
    from vulcan.gates import stats as GS
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    lab = MD.labels(C3.index_dir("human", "prod"))
    assert os.path.realpath(C3.index_dir("human", "prod")) == os.path.realpath(GATE_INDEX)
    thr = {}
    for run, psm, _ in MD.load_runs(C3.psm_dir("D10cf_prod"), with_feat=False):
        rows = psm["pf_row"]
        tb = hash_tiebreak(psm["pf_id"], 11)
        r = double_competition(psm["spec_idx"].astype(np.int64), psm["pf_id"], lab["tdc_pair"][rows],
                               psm["hyperscore"].astype(float), lab["is_decoy"][rows], seed=11, tiebreak=tb)
        thr[run] = float(GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01))
    fin = [v for v in thr.values() if np.isfinite(v)]          # NULL_RULES deviation note 1 (hosts without a threshold)
    out = {"host_hyperscore_1pct": {k: (v if np.isfinite(v) else None) for k, v in thr.items()},
           "tail_proxy_threshold": max(fin), "n_hosts_defined": len(fin),
           "rule": "max over the D10cf hosts whose 1% threshold is defined (>= 1 target accepted at 1%)"}
    os.makedirs(W, exist_ok=True)
    json.dump(out, open(os.path.join(W, "host_threshold.json"), "w"), indent=1)
    print(out["tail_proxy_threshold"])
    return out


def main_walk():
    base, it, cfg = load_cfg()
    a2_config.check_panel_files(base)
    frame = [r for r in csv.DictReader(open(os.path.join(ROOT, base["selection"]["frame"]), newline=""),
                                       delimiter="\t") if r.get("eligible") == it["selection"]["frame_filter"]["eligible"]]
    excl, ex_info = exclusions(base, list(csv.DictReader(open(os.path.join(ROOT, base["selection"]["frame"]),
                                                              newline=""), delimiter="\t")))
    cands, drops = S.frame_candidates(frame, cfg, excl)
    import select_train_itcid as ST
    n_dn2 = ST.dn2_guard([r["accession"] for r in cands], [r["organisms"] for r in cands])
    order = S.shuffled(cands, cfg)
    event(ev="start", walk_run=os.environ.get("DNIT_WALK_RUN", "1"), n_candidates=len(cands), frame_drops=drops, exclusions=ex_info, n_excl=len(excl),
          dn2_guard_checked=n_dn2, order=[r["accession"] for r in order],
          git_head=subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip())
    w = walk(order, cfg, IO(cfg, RAW))
    cols = ["id", "accession", "filename", "url", "expected_size", "purpose", "role", "set_id", "member"]
    with open(os.path.join(MAN, "dnit.tsv"), "w") as fh, open(os.path.join(MAN, "dnit.sha256"), "w") as sh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(w["rows"]):
            fh.write("\t".join([f"DNIT_{i:02d}", r["accession"], r["filename"], r["url"], str(r["expected_size"]),
                                f"IT-CID fresh null (DNIT): {r['organism']}; own proteome {r['upid']}; PRIDE SHA-1 "
                                f"{r['sha1']}", "dev", "DNIT", ""]) + "\n")
            sh.write(f"{r['sha256']}  {r['accession']}/{r['filename']}\n")
    json.dump({"seed": cfg["seeds"]["walk"], "n_frame_candidates": len(cands), "frame_drops": drops,
               "exclusions": ex_info, **{k: w[k] for k in ("admitted", "final_state", "log", "exhausted")}},
              open(os.path.join(MAN, "dnit_select.log.json"), "w"), indent=1, default=str)
    event(ev="end", admitted=w["admitted"], final_state=w["final_state"], exhausted=w["exhausted"])


def swap_test(acc):
    """P1 on the real gate index: rebuild with swap_labels=True, re-run the admission of ACC, compare membership."""
    import tempfile
    from vulcan.index.build import build
    base, it, cfg = load_cfg()
    meta = json.load(open(os.path.join(GATE_INDEX, "meta.json")))
    d = tempfile.mkdtemp(prefix="dnit_swap_", dir=W)
    fasta = meta["fasta"] if isinstance(meta["fasta"], list) else [meta["fasta"]]
    build(fasta, d, params=meta["params"], swap_labels=True)
    same = {f: a2_config.sha256_file(os.path.join(GATE_INDEX, f)) == a2_config.sha256_file(os.path.join(d, f))
            for f in ("peptidoforms.arrow", "residues.npy")}
    ref = json.load(open(os.path.join(W, "admission", acc, "counts.json")))
    files = [{"path": str(np.load(os.path.join(W, "admission", acc, r + ".npz"))["store"])} for r in ref["runs"]]
    # stores already exist; ingest() returns them unchanged
    rows = [{"path": os.path.join(RAW, acc, os.path.basename(f["path"])[:-len(".vspec.arrow")] + ".raw")}
            for f in files]
    for r, f in zip(rows, files):
        if not os.path.exists(r["path"]):
            r["path"] = f["path"].replace(".vspec.arrow", ".raw")
    res = run_project(acc, rows, {"upid": ref["upid"]}, cfg, gate_index=d, swap=True)
    out = {"project": acc, "label_free_files_identical": same, "membership_sha256_ref": ref["membership_sha256"],
           "membership_sha256_swapped": res.get("membership_sha256"),
           "pass": bool(all(same.values()) and res.get("membership_sha256") == ref["membership_sha256"])}
    json.dump(out, open(os.path.join(W, "swap_test.json"), "w"), indent=1)
    event(ev="swap_test", **out)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "host-threshold":
        host_threshold()
    elif cmd == "walk":
        main_walk()
    elif cmd == "plan":
        base, it, cfg = load_cfg()
        fr = list(csv.DictReader(open(os.path.join(ROOT, base["selection"]["frame"]), newline=""), delimiter="\t"))
        excl, info = exclusions(base, fr)
        cands, drops = S.frame_candidates([r for r in fr if r.get("eligible") == "1"], cfg, excl)
        print(info, drops, len(cands))
        for r in S.shuffled(cands, cfg):
            print(r["accession"], r["organisms"], r["instruments"])
    elif cmd == "swap-test":
        swap_test(sys.argv[2])
