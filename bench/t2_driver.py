#!/usr/bin/env python3
"""Pre-registered T2 driver: PTM localization / FLR (DESIGN §6.5, §6.6 T2; RESEARCH_TREE D4, D8; M4 E8 pipeline).

    scripts/m2b/vpy.sh bench/t2_driver.py --dataset D2  --look-id dev-YYYYmmdd-HHMM   # reproduce M4 E8 (dev)
    scripts/m2b/vpy.sh bench/t2_driver.py --dataset D2T --look-id dev-...             # T2-shaped dev analogue
    python3 bench/sealed_run.py --dataset T2 --verify --purpose "T2 look (FLR, DESIGN 6.5)" -- \
        scripts/m2b/vpy.sh bench/t2_driver.py --dataset T2 --look-id t2-look-1         # the one sealed look

Stages (one process, fresh work dir, nothing reused from earlier runs):
  preflight   config/t2.yaml, git sha, clean tree + single T2 ledger line (sealed), out dir must not exist, sha256 of
              every library/truth input against the committed manifests (sealed: mismatch aborts)
  convert     RAW -> node-local scratch -> vulcan.io.convert store (4 threads per run)
  index       FASTA = human SwissProt + contaminants + the dataset's libraries (M4 d2_fasta rule); L1 phospho index
  search      GPU closed search with STY phospho / M-ox / pep-N-term acetyl, report 1 (M4 d2_phospho search)
  e8          per accession: hyperscore paired double competition q <= 1% -> L3 localization (A primary, L secondary)
              and pyAscore on the same PSMs -> M4 evaluate + T2 additions (bench/t2/analysis.py)
  claims      bench/t2/claims.py -> claims.json / claims.md (deterministic, labels from config/t2.yaml)
  finish      predictor guard, completeness rule, manifest.sha256, leak test, deletion of all sealed-derived work
Outputs: <root>/<look-id>/out/{metrics.json, claims.json, claims.md, config.yaml, git_sha, timing.json,
manifest.sha256}; logs in <root>/<look-id>/logs/.  The terminal only shows lines matching config console_whitelist.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")   # required by torch deterministic mode (localization)
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")   # as scripts/m4/bench.sh
sys.dont_write_bytecode = True          # nothing may be written into the repository during a look (leak test)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench"))

CFG_PATH = os.path.join(ROOT, "config", "t2.yaml")


class Console:
    """Redirect fds 1/2 to a log file; only whitelisted lines reach the original terminal."""

    def __init__(self, log_path: str, whitelist: list[str]):
        self.pat = [re.compile(p) for p in whitelist]
        sys.stdout.flush(); sys.stderr.flush()
        self.term = os.fdopen(os.dup(1), "w", buffering=1)
        fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        os.dup2(fd, 1); os.dup2(fd, 2); os.close(fd)

    def say(self, line: str):
        print(line, flush=True)                                   # always to the log
        if any(p.match(line) for p in self.pat):
            self.term.write(line + "\n")


def git(*a):
    return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True).stdout.strip()


def core_s():
    import resource
    r = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime + c.ru_utime + c.ru_stime


def preflight(a, cfg):
    sealed = a.dataset == "T2"
    sha = git("rev-parse", "HEAD")
    dirty = bool(git("status", "--porcelain", "--untracked-files=no"))
    if sealed:
        if os.environ.get("SEALED_DATASET") != "T2" or not os.environ.get("SEALED_FILES"):
            raise SystemExit("t2: REFUSED not run through sealed_run")
        if git("status", "--porcelain"):
            raise SystemExit("t2: REFUSED tree not clean")
        led = os.path.join(ROOT, "results", "SEALED_LEDGER.md")
        n_t2 = sum(1 for l in open(led) if re.match(r"^\| [0-9-]+ [0-9:]+Z \| \w+ \| [^|]+ \| T2 \|", l))
        hist = git("log", "--format=%s", "--", "results/SEALED_LEDGER.md")
        n_hist = sum(1 for l in hist.splitlines() if l.startswith("sealed-ledger: T2 at"))
        if n_t2 != 1 or n_hist != 1:
            raise SystemExit("t2: REFUSED one look only")
    root = cfg["look_root"] if sealed else (a.dev_root or cfg["dev_root"])
    base = os.path.join(root, a.look_id)
    if os.path.exists(base):
        raise SystemExit("t2: REFUSED out dir exists")
    return sealed, sha, dirty, base


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, choices=["D2", "D2T", "T2"])
    ap.add_argument("--look-id", required=True)
    ap.add_argument("--keep-work", action="store_true", help="dev only: keep the scratch work dir")
    ap.add_argument("--skip-pyascore", action="store_true", help="dev only (completeness test): drop the comparator")
    ap.add_argument("--reuse-work", default=None, help="dev only: reuse stores + index from an earlier dev scratch dir")
    ap.add_argument("--dev-root", default=None, help="dev only: output root instead of config dev_root")
    ap.add_argument("--check", action="store_true",
                    help="operator pre-check BEFORE sealed_run (reads no sealed RAW): look root writable, config, "
                         "comparator hash, sealed truth/library hashes, GPU memory free")
    a = ap.parse_args()
    import yaml
    cfg = yaml.safe_load(open(CFG_PATH))
    for _k in ("look_root", "dev_root", "scratch_prefix"):   # paths may hold $VULCAN_LOOK_ROOT etc.
        if isinstance(cfg.get(_k), str):
            cfg[_k] = os.path.expandvars(cfg[_k])
    if a.check:
        sys.exit(check(a, cfg))
    if a.dataset == "T2" and (a.dev_root or a.reuse_work):
        raise SystemExit("t2: REFUSED dev root in sealed mode")
    sealed, sha, dirty, base = preflight(a, cfg)
    if sealed and (a.keep_work or a.skip_pyascore):
        raise SystemExit("t2: REFUSED dev-only option in sealed mode")
    out, logs = os.path.join(base, "out"), os.path.join(base, "logs")
    os.makedirs(base, mode=0o700)
    os.chmod(base, 0o700)
    os.makedirs(out); os.makedirs(logs)
    scratch = cfg["scratch_prefix"] + a.look_id
    if os.path.exists(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch, mode=0o700)
    stamp = os.path.join(base, ".start_stamp")
    open(stamp, "w").close()
    time.sleep(1.1)
    con = Console(os.path.join(logs, "driver.log"), cfg["console_whitelist"])
    start_guard(cfg, scratch)
    con.say(f"t2: start {a.dataset}")
    rc = 1
    try:
        run(a, cfg, sealed, sha, dirty, base, out, logs, scratch, con)
        rc = 0
    except SystemExit as e:
        traceback.print_exc()
        con.say("t2: FAILED " + re.sub(r"[^a-z_ ]", "", str(e).lower())[:60].strip() if str(e) else "t2: FAILED exit")
    except BaseException:
        traceback.print_exc()
        con.say("t2: FAILED memory guard" if GUARD["tripped"] else "t2: FAILED exception")
    finally:
        if sealed or not a.keep_work:
            shutil.rmtree(scratch, ignore_errors=True)
            if sealed:
                # leak test after deletion: nothing newer than the start outside the look root and scratch
                leaks = leak_test(stamp, scratch)
                json.dump({"leaks": leaks}, open(os.path.join(out, "leak_test.json"), "w"), indent=1)
                if leaks:
                    rc = rc or 3
                    con.say("t2: FAILED leak test")
        con.say(f"t2: exit {rc}")
    sys.exit(rc)


def check(a, cfg) -> int:
    """Pre-look check (no ledger line, no sealed RAW read)."""
    from t2 import datasets as DS
    ok = True
    root = cfg["look_root"] if a.dataset == "T2" else (a.dev_root or cfg["dev_root"])
    try:
        os.makedirs(root, exist_ok=True)
        t = os.path.join(root, ".write_test")
        open(t, "w").close(); os.remove(t)
        print(f"t2 check: look root writable: {root}")
    except OSError as e:
        ok = False
        print(f"t2 check: FAIL look root not writable: {root} ({type(e).__name__})")
    if os.path.exists(os.path.join(root, a.look_id)):
        ok = False
        print("t2 check: FAIL out dir for this look id exists")
    tree = pyascore_tree_sha(os.path.join(ROOT, cfg["comparator"]["pylib"]))
    ok &= tree == cfg["comparator"]["tree_sha256"]
    print(f"t2 check: pyAscore tree sha256 {'ok' if tree == cfg['comparator']['tree_sha256'] else 'MISMATCH'}")
    if a.dataset == "T2":
        fake = [os.path.join(ROOT, "data/sealed", acc, f) for acc, f in DS.sealed_raw_names()]
        ds = DS.dataset("T2", fake)
        try:
            DS.verify_inputs(DS.library_inputs(ds), strict=True)
            print("t2 check: sealed truth/library inputs match manifests (hashes only; contents not printed)")
        except SystemExit:
            ok = False
            print("t2 check: FAIL sealed truth/library hash mismatch")
        print("t2 check: sealed RAW files present: " + str(all(os.path.exists(f) for f in fake)))
    try:
        free = subprocess.run(["nvidia-smi", "--query-gpu=memory.used,memory.total", "--format=csv,noheader,nounits"],
                              capture_output=True, text=True).stdout.split(",")
        used = int(free[0])
        print(f"t2 check: GPU memory used {used} MiB (search peaks near 58 GB; needs an idle GPU)")
        ok &= used < 8000
    except Exception:
        ok = False
        print("t2 check: FAIL no GPU")
    st = os.statvfs(os.path.dirname(cfg["scratch_prefix"]))
    free_gb = st.f_bavail * st.f_frsize / 1e9
    need = cfg["memory_guard"]["scratch_min_free_gb"]
    print(f"t2 check: node-local scratch free {free_gb:.0f} GB (index ~30 GB + stores ~5 GB; need >= {need})")
    ok &= free_gb >= need
    h = cgroup_headroom_gb()
    need = cfg["memory_guard"]["headroom_before_index_gb"]
    print(f"t2 check: job cgroup memory headroom {h} GB (index build waits for >= {need}; D2T peak RSS 115 GB)")
    ok &= h is None or h >= need
    print(f"t2 check: load average {os.getloadavg()}")
    print("t2 check: " + ("OK" if ok else "FAIL"))
    return 0 if ok else 1


def leak_test(stamp: str, scratch: str) -> list[str]:
    dirs = [os.path.join(ROOT, d) for d in ("results", "data/cache", "data/truth", "bench", "vulcan", "scripts", "config")]
    dirs += [d for d in ("/tmp/vulcan_m4",) if os.path.isdir(d)]
    t0 = os.path.getmtime(stamp)
    leaks = []
    for d in dirs:
        for dp, dn, fn in os.walk(d):
            if "torch_extensions" in dp:             # JIT build cache: lock files only, no data
                continue
            for f in fn:
                p = os.path.join(dp, f)
                try:
                    if os.path.getmtime(p) > t0 and not p.endswith("results/SEALED_LEDGER.md"):
                        leaks.append(p)
                except OSError:
                    pass
    if os.path.exists(scratch):
        leaks.append(scratch)
    return leaks


def convert_all(raws, scratch, threads, con, logs):
    """Stage each RAW to node-local scratch, convert to a store, delete the staged copy (LESSONS 4)."""
    import concurrent.futures as cf
    from vulcan.hw import affinity_cpus
    store_dir = os.path.join(scratch, "store")
    stage = os.path.join(scratch, "rawstage")
    os.makedirs(store_dir); os.makedirs(stage)
    vpy = os.path.join(ROOT, "scripts", "m2b", "vpy.sh")
    workers = max(1, min(len(raws), len(affinity_cpus()) // threads))
    res = {}

    def one(raw):
        t0 = time.time()
        loc = os.path.join(stage, os.path.basename(raw))
        shutil.copyfile(raw, loc)
        t_copy = time.time() - t0
        name = os.path.basename(raw)[:-4]
        with open(os.path.join(logs, f"convert_{name}.json"), "w") as fh:
            p = subprocess.run([vpy, "-m", "vulcan.io.convert", loc, store_dir, "--threads", str(threads)],
                               stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
        os.remove(loc)
        if p.returncode:
            raise RuntimeError(f"convert failed for run {name}")
        meta = {}
        try:
            txt = open(os.path.join(logs, f"convert_{name}.json")).read()
            meta = json.loads(txt[txt.index("{"):txt.rindex("}") + 1])
        except Exception:
            pass
        return name, {"stage_copy_s": t_copy, "wall_s": time.time() - t0, "n_ms2": meta.get("n_ms2"),
                      "convert_timing_s": meta.get("timing_s"), "convert_cpu_s": meta.get("cpu_s"),
                      "store": os.path.join(store_dir, name + ".vspec.arrow")}

    with cf.ThreadPoolExecutor(workers) as ex:
        for k, (name, r) in enumerate(ex.map(one, raws)):
            res[name] = r
            con.say(f"t2: run {k + 1}/{len(raws)} convert done")
    return res, workers


def run(a, cfg, sealed, sha, dirty, base, out, logs, scratch, con):
    import numpy as np
    from t2 import datasets as DS
    from t2 import analysis as AN
    from t2 import claims as CL

    t_all, c_all = time.time(), core_s()
    timing = {"stages": {}}
    from vulcan.hw import fingerprint
    meta = {"dataset": a.dataset, "role": "sealed" if sealed else "dev", "look_id": a.look_id, "git_sha": sha,
            "tree_dirty_tracked": dirty, "config_version": cfg["version"],
            "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "hw": fingerprint(), "loadavg_start": os.getloadavg()}
    open(os.path.join(out, "git_sha"), "w").write(sha + "\n")
    shutil.copyfile(CFG_PATH, os.path.join(out, "config.yaml"))

    # ---------------------------------------------------------------- preflight: datasets, inputs, comparator, params
    con.say("t2: stage preflight start")
    sealed_files = os.environ.get("SEALED_FILES", "").split() if sealed else None
    ds = DS.dataset(a.dataset, sealed_files)
    lib_files = DS.library_inputs(ds)
    meta["inputs"] = DS.verify_inputs(lib_files, strict=sealed)
    if not sealed:
        raws = [r for d in ds["accessions"].values() for r in d["raws"]]
        meta["raw_inputs"] = DS.verify_inputs(raws, strict=False)
    pylib = os.path.join(ROOT, cfg["comparator"]["pylib"])
    tree = pyascore_tree_sha(pylib)
    meta["comparator"] = {"name": "pyAscore", "tree_sha256": tree, "expected": cfg["comparator"]["tree_sha256"]}
    if tree != cfg["comparator"]["tree_sha256"]:
        raise SystemExit("t2: REFUSED comparator hash mismatch")
    AN.set_pyascore_path(pylib)
    from vulcan.ptm import localize as LZ
    lp = LZ.LocParams(**{k: v for k, v in cfg["localization"].items()})
    if {k: v for k, v in lp.__dict__.items()} != cfg["localization"]:
        raise SystemExit("t2: REFUSED localization params differ from config")
    meta["loc_params"] = dict(lp.__dict__)
    guard_check()
    con.say("t2: stage preflight done")

    # ---------------------------------------------------------------- convert
    con.say("t2: stage convert start")
    t0 = time.time()
    raws = [r for d in ds["accessions"].values() for r in d["raws"]]
    if a.reuse_work:
        conv, workers = reuse_stores(raws, a.reuse_work, scratch), 0
    else:
        conv, workers = convert_all(raws, scratch, cfg["convert_threads"], con, logs)
    timing["stages"]["convert"] = {"wall_s": time.time() - t0, "parallel_workers": workers,
                                   "threads_per_run": cfg["convert_threads"]}
    guard_check()
    con.say("t2: stage convert done")

    # ---------------------------------------------------------------- index
    con.say("t2: stage index start")
    t0, c0 = time.time(), core_s()
    from vulcan.index.build import build, L1Index
    fasta = DS.write_fasta(ds, os.path.join(scratch, "search.fasta"))
    meta["search_fasta_sha256"] = DS.sha256(fasta)
    idx = os.path.join(scratch, "index")
    s = cfg["search"]
    if a.reuse_work:
        prev = os.path.join(a.reuse_work, "index")
        if DS.sha256(os.path.join(a.reuse_work, "search.fasta")) != meta["search_fasta_sha256"]:
            raise SystemExit("t2: REFUSED reused index from another fasta")
        os.symlink(prev, idx)
        meta["reused_work"] = a.reuse_work
    else:
        wait_mem(con, cfg["memory_guard"]["headroom_before_index_gb"])
        build([fasta], idx, {"variable_mods": s["variable"], "max_variable_mods": s["max_variable_mods"],
                             "static_mods": s["static"], "entrap_seeds": []})
    l1 = L1Index(idx)
    meta["index_meta"] = {k: v for k, v in l1.meta.items() if k not in ("fasta",)}
    timing["stages"]["index_build"] = {"wall_s": time.time() - t0, "core_s": core_s() - c0}
    guard_check()
    con.say("t2: stage index done")

    # ---------------------------------------------------------------- search
    con.say("t2: stage search start")
    t0 = time.time()
    psm_dir = os.path.join(scratch, "psms")
    os.makedirs(psm_dir)
    per_run = search_all(conv, idx, l1, s, psm_dir, con)
    timing["stages"]["search"] = {"wall_s": time.time() - t0}
    guard_check()
    con.say("t2: stage search done")

    # ---------------------------------------------------------------- E8
    con.say("t2: stage e8 start")
    t_e8 = time.time()
    truths = DS.truth_for(ds)
    st = cfg["stats"]
    metrics = {"meta": meta, "truth_reports": {k: v["report"] for k, v in truths.items()}, "datasets": {}}
    for acc, d in ds["accessions"].items():
        t_acc = time.time()
        runs = [DS.run_name(r) for r in d["raws"]]
        psm, counts = AN.psm_list(runs, psm_dir, idx, l1, q=cfg["fdr"]["q"], seed=cfg["fdr"]["tiebreak_seed"])
        T = truths[acc]["truth"]
        truth_of = [T.get((rn, sq)) if len(ph) == 1 else None for rn, sq, ph in zip(psm["run"], psm["seq"], psm["phos"])]
        spectra, t_sp = AN.load_spectra([conv[r]["store"] for r in runs], list(zip(psm["run"], psm["spec_idx"])))
        comp = AN.e8_component(acc, psm, truth_of, spectra, lp, run_pyascore=not a.skip_pyascore,
                               boot_n=st["bootstrap"]["n"], boot_seed=st["bootstrap"]["seed"],
                               n_min=st["n_min_accepted_eligible_sites"], log=print,
                               repeats=st["localization_repeats"])
        comp.update(counts)
        comp["search_mode_counts"] = {"nterm_acetyl": int(sum(psm["nterm"]))}
        comp["timing_s"] = {"load_spectra": t_sp, "total": time.time() - t_acc}
        comp["runs"] = runs
        metrics["datasets"][acc] = comp
    timing["stages"]["e8"] = {"wall_s": time.time() - t_e8}
    guard_check()
    con.say("t2: stage e8 done")

    # ---------------------------------------------------------------- runtime per run
    metrics["runtime_per_run"] = {r: {"convert_wall_s": conv[r]["wall_s"], "stage_copy_s": conv[r]["stage_copy_s"],
                                      "convert_timing_s": conv[r]["convert_timing_s"], "n_ms2": conv[r]["n_ms2"],
                                      **per_run[r]} for r in conv}
    timing["total_wall_s"] = time.time() - t_all
    timing["total_core_s"] = core_s() - c_all
    import resource
    timing["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    meta["loadavg_end"] = os.getloadavg()
    metrics["timing"] = timing

    # ---------------------------------------------------------------- guards + outputs
    bad = sorted(m for m in sys.modules if m.split(".")[0] in set(cfg["predictors"]["forbidden_modules"]))
    metrics["predictor_guard"] = {"forbidden_imported": bad, "used": cfg["predictors"]["used"]}
    if bad:
        raise SystemExit("t2: FAILED forbidden predictor imported")
    with open(os.path.join(out, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1, default=_jd, sort_keys=True)
    json.dump(timing, open(os.path.join(out, "timing.json"), "w"), indent=1, default=_jd)
    con.say("t2: stage claims start")
    rep = CL.write_claims(os.path.join(out, "metrics.json"), cfg, out)
    con.say("t2: stage claims done")
    write_manifest(out)
    guard_check()
    if rep["incomplete"]:
        raise SystemExit("t2: FAILED look incomplete")
    con.say(f"t2: done {a.dataset}")


def _jd(x):
    import numpy as np
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, tuple):
        return list(x)
    return str(x)


def search_all(conv, idx, l1, s, psm_dir, con):
    import pyarrow.parquet as pq
    import torch
    from vulcan.gpu import engine as E
    wait_gpu(con)
    p = E.SearchParams(**s["params"])
    p_report = s["report"]
    t0 = time.time()
    torch.zeros(1, device="cuda")
    gi = E.GpuIndex(l1, p, s["static"], s["variable"])
    torch.cuda.synchronize()
    setup = time.time() - t0
    out = {}
    names = sorted(conv)
    for k, name in enumerate(names):
        t0, c0 = time.time(), core_s()
        sp = E.load_store(conv[name]["store"])
        tbl, st, _ = E.search(sp, gi, l1, p, report=p_report)
        torch.cuda.synchronize()
        pq.write_table(tbl, os.path.join(psm_dir, f"{sp.run}.psms.parquet"))
        out[name] = {"search_wall_s": time.time() - t0, "search_core_s": core_s() - c0,
                     "search_engine_wall_s": st.wall_s, "n_psms": tbl.num_rows}
        con.say(f"t2: run {k + 1}/{len(names)} search done")
    out["_setup"] = {"gpu_index_build_s": setup, "gpu_max_mem_bytes": torch.cuda.max_memory_allocated()}
    del gi
    torch.cuda.empty_cache()
    setup_info = out.pop("_setup")
    for v in out.values():
        v["gpu_index_setup_s_shared"] = setup_info["gpu_index_build_s"]
    return out


def wait_gpu(con, need_gb: float = 64.0, timeout_s: float = 5400.0):
    """The phospho index needs ~58 GB peak on the GPU; wait (logged) until that much is free instead of failing."""
    import torch
    t0 = time.time()
    while True:
        free, total = torch.cuda.mem_get_info()
        if free / 1e9 >= need_gb:
            print(f"gpu free {free / 1e9:.1f} GB after {time.time() - t0:.0f} s wait", flush=True)
            return time.time() - t0
        if time.time() - t0 > timeout_s:
            raise SystemExit("t2: FAILED gpu busy")
        time.sleep(15)


def cgroup_headroom_gb():
    """Free memory (GB) under the job's cgroup limit (batch scheduler), or None outside a limited cgroup."""
    try:
        rel = open("/proc/self/cgroup").read().strip().split("::", 1)[1]
        best = None
        parts = rel.strip("/").split("/")
        for k in range(len(parts), 0, -1):
            d = "/sys/fs/cgroup/" + "/".join(parts[:k])
            mx = open(os.path.join(d, "memory.max")).read().strip()
            if mx != "max":
                h = (int(mx) - int(open(os.path.join(d, "memory.current")).read())) / 1e9
                best = h if best is None else min(best, h)
        return best
    except Exception:
        return None


def start_guard(cfg, scratch):
    """Watchdog (LESSONS 7): stop this driver if the job cgroup headroom or scratch free space falls below the floor,
    so a shared-node squeeze fails this run cleanly instead of OOM-killing another agent's process."""
    import _thread
    import threading
    g = cfg["memory_guard"]

    def loop():
        while True:
            h = cgroup_headroom_gb()
            st = os.statvfs(scratch)
            free = st.f_bavail * st.f_frsize / 1e9
            if (h is not None and h < g["abort_floor_gb"]) or free < 2:
                print(f"guard: abort (cgroup headroom {h} GB, scratch free {free:.1f} GB)", flush=True)
                GUARD["tripped"] = True
                _thread.interrupt_main()
                return
            time.sleep(2)
    threading.Thread(target=loop, daemon=True).start()


GUARD = {"tripped": False}


def guard_check():
    """interrupt_main() cannot interrupt a long native call (the index build), and a KeyboardInterrupt raised in a
    library may be swallowed (dev D2T run dev-d2t-3: the guard tripped and the run still exited 0).  So every stage
    boundary and the end of the run check the flag; a tripped guard always fails the run."""
    if GUARD["tripped"]:
        raise SystemExit("t2: FAILED memory guard")


def wait_mem(con, need_gb: float = 170.0, timeout_s: float = 5400.0):
    """The phospho index build needs tens of GB of RAM; a shared job memory cgroup OOM-killed the first D2T dry run.
    Wait (logged) until the job cgroup has need_gb of headroom."""
    t0 = time.time()
    while True:
        h = cgroup_headroom_gb()
        if h is None or h >= need_gb:
            print(f"cgroup memory headroom {h} GB after {time.time() - t0:.0f} s wait", flush=True)
            return
        if time.time() - t0 > timeout_s:
            raise SystemExit("t2: FAILED memory busy")
        time.sleep(15)


def reuse_stores(raws, prev, scratch):
    """dev only: link stores of an earlier dev run (convert timings are then not measured)."""
    os.makedirs(os.path.join(scratch, "store"))
    out = {}
    for r in raws:
        name = os.path.basename(r)[:-4]
        src = os.path.join(prev, "store", name + ".vspec.arrow")
        dst = os.path.join(scratch, "store", name + ".vspec.arrow")
        os.symlink(src, dst)
        out[name] = {"stage_copy_s": None, "wall_s": None, "n_ms2": None, "convert_timing_s": None,
                     "convert_cpu_s": None, "store": dst, "reused": True}
    return out


def pyascore_tree_sha(pylib: str) -> str:
    import hashlib
    lines = []
    for sub in ("pyascore", "pyascore-1.0.1.dist-info"):
        for dp, dn, fn in os.walk(os.path.join(pylib, sub)):
            dn[:] = [d for d in dn if d != "__pycache__"]
            for f in fn:
                if f.endswith(".pyc"):
                    continue
                p = os.path.join(dp, f)
                h = hashlib.sha256(open(p, "rb").read()).hexdigest()
                lines.append((os.path.relpath(p, pylib), h))
    txt = "".join(f"{h}  {r}\n" for r, h in sorted(lines))
    return hashlib.sha256(txt.encode()).hexdigest()


def write_manifest(out):
    import hashlib
    lines = []
    for f in sorted(os.listdir(out)):
        if f == "manifest.sha256":
            continue
        lines.append(f"{hashlib.sha256(open(os.path.join(out, f), 'rb').read()).hexdigest()}  {f}\n")
    open(os.path.join(out, "manifest.sha256"), "w").write("".join(lines))


if __name__ == "__main__":
    main()
