#!/usr/bin/env python3
"""A6.8.2 step 2 dry run on DEV data only (D0, D8h): end-to-end admission + throughput r.

DEV projects only (never an HN1 candidate; no PRIDE listing of any HN1 candidate).  Per run, pipelined exactly as the
walk would run it (agent S): PRIDE listing of the (dev) project -> download (size + PRIDE SHA-1) -> decode to a
vspec store + run checks (A6.3.4(d)) -> TRFP mzML -> Sage per project -> [GPU, /tmp/vulcan_a6.lock] certification,
I, R_P, freeze, reduced index, null search -> swap test (every project).  Per-stage wall and core seconds and an
event timeline are written to results/2026-10-28_A6_hn1/dry_run/.

Label discipline: identical to admit_hn1 (labels only inside certification, at host rank-1 rows).  The outputs
are label-free counts and memberships.

    scripts/m2b/vpy.sh scripts/a6/dry_run.py --set D8h --workers 6
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import glob
import hashlib
import json
import os
T1_WORK = os.environ.get("VULCAN_T1_WORK", f"/tmp/vulcan_t1_{os.getuid()}")  # bench/t1_driver.py default work root
import queue
import resource
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "e9"))
OUT = os.path.join(ROOT, "results", "2026-10-28_A6_hn1", "dry_run")
WORK = "/tmp/vulcan_a6/dry"
LOCK = "/tmp/vulcan_a6.lock"
TRFP = os.path.join(ROOT, "conda_envs", "trfp", "bin", "ThermoRawFileParser")
SAGE_FASTA = T1_WORK + "/D8h/fasta/prod_1001.fasta"   # T1 prod DB (prod index T + dec_1001), rev_ decoys
VPY = os.path.join(ROOT, "scripts", "m2b", "vpy.sh")
TMT = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825, 130.141145,
       131.138180]
ITRAQ = [113.107873, 114.111228, 115.108263, 116.111618, 117.114973, 118.112008, 119.115363, 121.122072]
T0 = time.time()
EV = []
_evl = threading.Lock()


def ev(**k):
    k["t"] = round(time.time() - T0, 2)
    with _evl:
        EV.append(k)
    print(json.dumps(k), flush=True)


def child_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def timed(cmd, log):
    t = log + ".time"
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.call(["/usr/bin/time", "-v", "-o", t] + cmd, stdout=fh, stderr=subprocess.STDOUT)
    core = None
    try:
        import re
        txt = open(t).read()
        core = float(re.search(r"User time \(seconds\): ([\d.]+)", txt).group(1)) + \
            float(re.search(r"System time \(seconds\): ([\d.]+)", txt).group(1))
    except Exception:  # noqa: BLE001
        pass
    return rc, time.time() - t0, core


# ------------------------------------------------------------------------------------------------ stages (per run)
def listing(acc):
    import select_e9 as SE
    t0 = time.time()
    try:
        fl = SE.pride_files(acc)
        n = len(fl)
    except Exception as e:  # noqa: BLE001
        n = f"error {e}"
    return {"wall_s": time.time() - t0, "n_files": n}


def download(url, size, sha1, dst):
    t0 = time.time()
    if os.path.exists(dst) and os.path.getsize(dst) == size:
        return {"wall_s": 0.0, "reused": True}
    part = dst + ".part"
    subprocess.run(["curl", "-s", "-S", "--retry", "3", "--connect-timeout", "60", "-o", part, url])
    ok = os.path.exists(part) and os.path.getsize(part) == size
    h = hashlib.sha1()
    if ok and sha1:
        with open(part, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 24), b""):
                h.update(b)
        ok = h.hexdigest() == sha1
    if ok:
        os.replace(part, dst)
    return {"wall_s": time.time() - t0, "ok": ok, "bytes": size, "sha1_checked": bool(sha1)}


def decode(raw, store_path, run):
    if os.path.exists(store_path):
        return {"wall_s": 0.0, "reused": True}
    code = (f"from vulcan.io.thermo_direct import convert; convert({raw!r}, {store_path!r}, run_name={run!r})")
    rc, w, core = timed([VPY, "-c", code], store_path + ".decode.log")
    return {"wall_s": w, "core_s": core, "rc": rc}


def run_checks(store_path):
    """A6.3.4(d): >= 99% MS2 FTMS HCD; no TMT/iTRAQ reporters; MS1 present; decodable."""
    from vulcan.io import store
    t = store.read_table(store_path)
    lvl = t.column("ms_level").to_numpy()
    act = t.column("activation").to_numpy()[lvl == 2]
    ana = t.column("analyzer").to_numpy()[lvl == 2]
    frac = float(((act == 2) & (ana == 1)).mean()) if act.size else 0.0
    t2 = t.filter(__import__("pyarrow").array(lvl == 2))
    pk = store.peaks_numpy(t2)
    n = len(pk.offsets) - 1
    rng = np.random.default_rng(0)
    sel = rng.choice(n, min(n, 3000), replace=False) if n else np.zeros(0, int)

    def frac_hit(targets):
        hits = []
        for mz0 in targets:
            c = 0
            for i in sel:
                m = pk.mz[pk.offsets[i]:pk.offsets[i + 1]]
                j = np.searchsorted(m, mz0 * (1 - 1e-5))
                c += int(j < m.size and m[j] <= mz0 * (1 + 1e-5))
            hits.append(c / max(sel.size, 1))
        return float(np.median(hits))
    rep = {}
    for name, tg in (("TMT", TMT), ("iTRAQ", ITRAQ)):
        f, ctl = frac_hit(tg), frac_hit([x + 0.3 for x in tg])
        rep[name] = {"median_frac": f, "control_frac": ctl, "detected": bool(f >= 0.3 and f >= 5 * max(ctl, 1e-3))}
    ok = frac >= 0.99 and int((lvl == 1).sum()) > 0 and not any(v["detected"] for v in rep.values())
    return {"frac_ftms_hcd": frac, "n_ms1": int((lvl == 1).sum()), "reporters": rep, "pass": bool(ok)}


def trfp(raw, mzml_dir, run):
    mz = os.path.join(mzml_dir, run + ".mzML")
    if os.path.exists(mz):
        return {"wall_s": 0.0, "reused": True}
    rc, w, core = timed([TRFP, "-i", raw, "-o", mzml_dir, "-f", "2", "-l", "3"], os.path.join(mzml_dir, run + ".trfp.log"))
    return {"wall_s": w, "core_s": core, "rc": rc}


# ------------------------------------------------------------------------------------------------ main
def projects_for(setname):
    rows = list(csv.DictReader(open(os.path.join(ROOT, "data/manifests/e9_d8.tsv")), delimiter="\t"))
    if setname == "D0":
        return [{"acc": "PXD001468", "runs": [
            {"run": os.path.basename(r)[:-4], "local_raw": r} for r in
            sorted(glob.glob("/tmp/vulcan_m1e/raw/b19*_293T_proteinID_*.raw"))]}]
    if setname == "DEV12":   # yields on the pilot's 12 DEV runs: existing stores + existing mzML (no download)
        d0 = [{"run": os.path.basename(x)[:-len(".vspec.arrow")], "local_store": x,
               "local_mzml": T1_WORK + "/D0/mzml/" + os.path.basename(x)[:-len(".vspec.arrow")] + ".mzML"}
              for x in sorted(glob.glob("/tmp/vulcan_e2/store/b19*_293T_proteinID_*.vspec.arrow"))]
        out = [{"acc": "PXD001468", "runs": d0}]
        acc = {r["id"]: r["accession"] for r in rows}
        for x in sorted(glob.glob(T1_WORK + "/D8h/stage/E9_*.vspec.arrow")):
            run = os.path.basename(x)[:-len(".vspec.arrow")]
            out.append({"acc": acc[run.split("__")[0]], "runs": [{
                "run": run, "local_store": x, "local_mzml": f"{T1_WORK}/D8h/mzml/{run}.mzML"}]})
        return out
    ids = {"D8h2": ["E9_0001", "E9_0005"], "D8h": [f"E9_{i:04d}" for i in range(10)]}[setname]
    out = []
    for r in rows:
        if r["id"] in ids:
            out.append({"acc": r["accession"], "runs": [{"run": f"{r['id']}__{r['filename'][:-4]}", "url": r["url"],
                                                          "size": int(r["expected_size"]), "sha1": r["pride_sha1"]}]})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="D8h2", choices=["D0", "D8h2", "D8h", "DEV12"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--sage-threads", type=int, default=8)
    ap.add_argument("--no-swap", action="store_true")
    ap.add_argument("--exclude", nargs="*", default=[], help="accessions to skip (already processed in a part-1 run)")
    ap.add_argument("--out-tag", default=None, help="output sub-directory name (default: the set name)")
    ap.add_argument("--lock", default=LOCK, help="A6 GPU lock file (never the bench lock)")
    a = ap.parse_args()
    out = os.path.join(OUT, a.out_tag or a.set)
    work = os.path.join(WORK, a.set)
    for d in (out, work, os.path.join(work, "raw"), os.path.join(work, "stage"), os.path.join(work, "mzml")):
        os.makedirs(d, exist_ok=True)
    projs = [pr for pr in projects_for(a.set) if pr["acc"] not in set(a.exclude)]
    ev(ev="start", set=a.set, n_projects=len(projs), n_runs=sum(len(p["runs"]) for p in projs),
       loadavg=os.getloadavg())
    ready = queue.Queue()

    def prep_project(pr):
        acc = pr["acc"]
        if any("url" in r for r in pr["runs"]):
            li = listing(acc)
            ev(ev="listing", acc=acc, **li)
        res = []
        for r in pr["runs"]:
            run = r["run"]
            if "local_store" in r:
                st = os.path.join(work, "stage", run + ".vspec.arrow")
                mz = os.path.join(work, "mzml", run + ".mzML")
                for src, dst in ((r["local_store"], st), (r["local_mzml"], mz)):
                    if not os.path.lexists(dst):
                        os.symlink(src, dst)
                ck = run_checks(st)
                ev(ev="run_checks", run=run, **ck)
                res.append((run, st, mz, ck["pass"]))
                continue
            raw = r.get("local_raw") or os.path.join(work, "raw", run + ".raw")
            if "url" in r:
                d = download(r["url"], r["size"], r["sha1"], raw)
                ev(ev="download", run=run, **d)
            st = os.path.join(work, "stage", run + ".vspec.arrow")
            de = decode(raw, st, run)
            ev(ev="decode", run=run, **de)
            ck = run_checks(st)
            ev(ev="run_checks", run=run, **ck)
            tr = trfp(raw, os.path.join(work, "mzml"), run)
            ev(ev="trfp", run=run, **tr)
            res.append((run, st, os.path.join(work, "mzml", run + ".mzML"), ck["pass"]))
        runs_ok = [x for x in res if x[3]]
        if runs_ok:
            sys.path.insert(0, ROOT)
            import bench.t1.comparators as CMP
            t0 = time.time()
            rec = CMP.run_sage(acc, SAGE_FASTA, [x[2] for x in runs_ok], os.path.join(work, "sage"),
                               {"prec_ppm": 10.0, "frag_ppm": 20.0}, a.sage_threads, tag="sage")
            ev(ev="sage", acc=acc, n_runs=len(runs_ok), wall_s=time.time() - t0, core_s=rec.get("core_s"))
        ready.put((acc, runs_ok))

    ex = ThreadPoolExecutor(a.workers)
    for pr in projs:
        ex.submit(prep_project, pr)
    # GPU stage in this process, under the A6 lock (not the bench lock)
    import admit_hn1 as AH
    lockf = open(a.lock, "w")
    fcntl.flock(lockf, fcntl.LOCK_EX)
    ev(ev="a6_lock_acquired")
    t0 = time.time()
    c0 = resource.getrusage(resource.RUSAGE_SELF)
    ctx = AH.Ctx()
    ev(ev="ctx_setup", wall_s=time.time() - t0)
    ctx_sw = None
    if not a.no_swap:
        t0 = time.time()
        ctx_sw = AH.Ctx(AH.IDX_SWAPPED)
        ev(ev="ctx_swapped_setup", wall_s=time.time() - t0)
    summary = {"projects": {}}
    for _ in projs:
        acc, runs_ok = ready.get()
        if not runs_ok:
            ev(ev="project_not_admitted", acc=acc)
            continue
        s1 = AH.sage_rank1(os.path.join(work, "sage", "sage", acc, "results.sage.tsv"))
        t0 = time.time()
        r0 = resource.getrusage(resource.RUSAGE_SELF)
        res = AH.admit_project(ctx, acc, [(r, st) for r, st, _, _ in runs_ok], s1,
                               os.path.join(out, "freeze", acc))
        r1 = resource.getrusage(resource.RUSAGE_SELF)
        ev(ev="admit", acc=acc, wall_s=time.time() - t0,
           cpu_s=(r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime), timing=res["timing_s"],
           removal=res["removal"], counts={r: v["counts"] for r, v in res["per_run"].items()})
        rec = {"counts": {r: v["counts"] for r, v in res["per_run"].items()}, "removal": res["removal"],
               "timing_s": res["timing_s"], "n_I_pairs": int(res["I_keys"].size), "thr": res["thr"]}
        if ctx_sw is not None:
            t0 = time.time()
            sw = AH.swap_test(ctx_sw, acc, [(r, st) for r, st, _, _ in runs_ok], s1, os.path.join(out, "freeze", acc),
                              res)
            ev(ev="swap_test", acc=acc, wall_total_s=time.time() - t0, **{k: v for k, v in sw.items() if k != "project"})
            rec["swap_test"] = sw
        n_lab_rows = int(sum(x.size for x in ctx.labels.rows_read))
        rec["label_rows_read_cumulative"] = n_lab_rows
        summary["projects"][acc] = rec
        ev(ev="admitted", acc=acc, n_runs=len(runs_ok))
    ex.shutdown()
    fcntl.flock(lockf, fcntl.LOCK_UN)
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    summary["wall_s_total"] = time.time() - T0
    summary["self_cpu_s"] = (c1.ru_utime + c1.ru_stime) - (c0.ru_utime + c0.ru_stime)
    summary["children_cpu_s"] = child_cpu()
    summary["loadavg_end"] = os.getloadavg()
    summary["events"] = EV
    json.dump(summary, open(os.path.join(out, "dry_run.json"), "w"), indent=1, default=float)
    ev(ev="done", wall_s=summary["wall_s_total"])


if __name__ == "__main__":
    main()
