"""Cohort driver (M2b, DESIGN §4 L0-L2 micro-cohort batching, §10 E3 groundwork).

One long-lived process turns a list of Thermo RAW files into per-run PSM tables, fully pipelined:

    N decode workers (processes, persistent: .NET runtime + JIT + imports paid once per worker)
        RAW (Lustre, sequential prefetch) -> Arrow store staged on local NVMe
    loader threads      store (mmap, local) -> MS2 Spectra (numpy)
    GPU stage (main)    L1 index built ONCE and kept resident; K0..K2 + features per run
    writer threads      PSM parquet (+ optional family top-2 table), staged store deleted

Decode of run i+1.. overlaps GPU search of run i. Jobs are submitted largest-file-first (LPT) so the last wave
does not leave cores idle behind one big file.

    python -m vulcan.cohort --index DIR --raws a.raw b.raw ... --out OUT [--workers 24] [--stage /tmp/...]
           [--raw-list FILE] [--repeat-to 50] [--keep-stores] [--family-top2] [--torch-threads 4]

Writes OUT/<run>.psms.parquet (the vulcan.gpu.engine.psm_table schema, unchanged; see vulcan/gpu/PSM_SCHEMA.md) and
OUT/cohort_stats.json: per-run stage walls / core-s / GPU ms, timeline, GPU utilization samples, bottleneck summary.
Nothing here reads labels (family top-2 uses only the opaque per-row group id, see vulcan/gpu/family.py).

Optional production scoring (`--score-lib DIR`, vulcan/rescore/cohort_scoring.py): after each search the GPU-stage
thread computes the L3 predicted features (per-run label-blind NCE choice, resident predicted library) and the
frozen static_gbm_pred score on the GPU (`score` column added to the PSM table); paired double competition + picked
protein groups (these read the index's decoy/pair labels, as every FDR step must) run in torch-free post worker
processes and write OUT/<run>.peptides.parquet and OUT/<run>.proteins.parquet.
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import queue
import resource
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor



def process_age_s() -> float:
    """Seconds since this process was exec'd (interpreter start + imports so far), from /proc."""
    try:
        with open("/proc/self/stat") as fh:
            start_ticks = int(fh.read().rsplit(")", 1)[1].split()[19])
        with open("/proc/uptime") as fh:
            up = float(fh.read().split()[0])
        return up - start_ticks / os.sysconf("SC_CLK_TCK")
    except Exception:
        return float("nan")


def _rus(who=resource.RUSAGE_SELF):
    r = resource.getrusage(who)
    return r.ru_utime + r.ru_stime


# ------------------------------------------------------------------------------------------- decode workers
_PREP = None
_MS1_BACKEND = "cpu"


def _worker_init(prep: dict | None = None, ms1_cfg: dict | None = None):
    """Runs once per decode worker: 1-thread pools, .NET runtime + RawFileReader assemblies loaded once; with
    `prep` ({index, params}) also the index masses, so the worker can emit search-ready MS2 arrays + hypotheses.
    ms1_cfg {"backend": "gpu", "gate": Semaphore} (E11, opt-in): MS1 hills + isotopes on the GPU through the torch-free
    vulcan.ms1.gpu library (own CUDA context per worker, created here); the semaphore bounds concurrent GPU MS1 runs
    and each run's cached device memory is returned to the driver afterwards (VULCAN_MS1_GPU_TRIM=1)."""
    global _PREP, _MS1_BACKEND
    os.environ["VULCAN_THREADS"] = "1"
    from vulcan.hw import pin_thread_pools
    pin_thread_pools(1)
    if ms1_cfg and ms1_cfg.get("backend", "cpu") == "gpu":
        os.environ["VULCAN_MS1_GPU_TRIM"] = "1"
        from vulcan.ms1 import gpu as G
        G.set_gate(ms1_cfg.get("gate"))
        G.session()
        _MS1_BACKEND = "gpu"
    from vulcan.io.store import configure_threads
    configure_threads(1)
    from vulcan.io import thermo_direct as T
    T.load_runtime()
    if prep:
        import types

        import pyarrow as pa
        with pa.memory_map(os.path.join(prep["index"], "peptidoforms.arrow")) as src:
            m32 = pa.ipc.open_file(src).read_all().column("mass32").to_numpy().copy()
        _PREP = {"gi": types.SimpleNamespace(mass32_np=m32), "params": prep["params"]}


PREP_ARRAYS = ("scan", "index", "rt", "precursor_mz", "charge", "mz", "intensity", "offsets")
HYP_ARRAYS = ("spec", "z", "iso", "pm", "M", "lo", "hi", "maxfz")


def _prep_search(store: str, run_name: str) -> tuple[str, float, float]:
    """In the decode worker: MS2 Spectra + precursor hypotheses (engine.pre_hypotheses) written as .npy files next to
    the staged store, so the GPU process only memory-maps them (no Arrow flatten / numpy in its serial path)."""
    import numpy as np

    from vulcan.gpu.host import SearchParams, pre_hypotheses
    t0, c0 = time.time(), _rus()
    sp = load_ms2(store)
    valid, hyp = pre_hypotheses(sp, _PREP["gi"], SearchParams(**_PREP["params"]))
    d = store[: -len(".vspec.arrow")] + ".prep"
    os.makedirs(d, exist_ok=True)
    for k in PREP_ARRAYS:
        np.save(os.path.join(d, k + ".npy"), getattr(sp, k))
    np.save(os.path.join(d, "valid.npy"), valid)
    for k in HYP_ARRAYS:
        np.save(os.path.join(d, "hyp_" + k + ".npy"), hyp[k])
    with open(os.path.join(d, "run.txt"), "w") as fh:
        fh.write(sp.run)
    return d, time.time() - t0, _rus() - c0


def load_prep(d: str):
    """GPU process side of _prep_search: memory-mapped arrays (local NVMe / page cache)."""
    import numpy as np

    from vulcan.gpu.host import Spectra
    a = {k: np.load(os.path.join(d, k + ".npy"), mmap_mode="r") for k in PREP_ARRAYS}
    sp = Spectra(open(os.path.join(d, "run.txt")).read(), **a)
    hyp = {k: np.load(os.path.join(d, "hyp_" + k + ".npy")) for k in HYP_ARRAYS}
    return sp, (np.load(os.path.join(d, "valid.npy")), hyp)


def _ms1_stage(store: str, run_name: str, out_dir: str) -> dict:
    """Optional MS1 stage (vulcan.ms1, CPU, 1 numba thread in this worker): features + precursor-hypothesis table,
    written next to the PSM tables. Not yet consumed by the GPU search (engine.make_hypotheses still uses the
    instrument charge); measured here for its share of the cohort budget."""
    from vulcan.ms1.features import FeatureParams, find_features, write_features
    from vulcan.ms1.hypotheses import build_hypotheses, write_table
    t0, c0 = time.time(), _rus()
    ref = os.environ.get("VULCAN_MS1_REFERENCE") == "1"      # E7b comparison only: the original (slower) kernels
    fr = find_features(store, FeatureParams(), device="cpu", exact_reference=ref, backend=_MS1_BACKEND,
                       hill_peaks=False)     # GPU only: per-peak hill arrays stay on the device (unused here)
    write_features(fr.features, os.path.join(out_dir, f"{run_name}.features.arrow"))
    tbl, _ = build_hypotheses(store, fr)
    write_table(tbl, os.path.join(out_dir, f"{run_name}.hypotheses.arrow"))
    rec = {"ms1_wall_s": time.time() - t0, "ms1_core_s": _rus() - c0, "n_features": fr.features.num_rows,
           "n_hypotheses_ms1": tbl.num_rows}
    if fr.timings.get("backend"):            # gpu, or cpu_fallback (device memory short; bitwise-equal tables)
        rec["ms1_backend"] = fr.timings["backend"]
        rec["ms1_gpu_s"] = sum(fr.timings.get("gpu_ms", {}).values()) / 1e3
        rec["ms1_gpu_gate_wait_s"] = fr.timings.get("gpu_gate_wait", 0.0)
    return rec


def _ingest(raw: str, store: str, run_name: str, ms1_out: str | None = None, stage_raw: bool = True) -> dict:
    """One RAW -> staged store (+ optional MS1 stage), in a persistent worker. Returns timings and core-s.

    stage_raw: copy the RAW to local NVMe first (one sequential Lustre read, in-kernel copy) and decode the local copy.
    RawFileReader memory-maps the RAW and reads it out of order; on Lustre, page-cache reclaim under the job's memory
    limit turned that into random Lustre reads (decode 5 s -> 77 s stalls). DESIGN §4 L0: never mmap on Lustre."""
    from vulcan.io import thermo_direct as T
    t0, c0 = time.time(), _rus()
    src, t_copy = raw, 0.0
    if stage_raw:
        src = store[: -len(".vspec.arrow")] + ".raw"
        shutil.copyfile(raw, src)
        t_copy = time.time() - t0
    try:
        meta = T.convert(src, store, run_name=run_name, do_prefetch=not stage_raw,
                         source_path=os.path.abspath(raw))
    finally:
        if stage_raw and os.path.exists(src):
            os.remove(src)
    t1, c1 = time.time(), _rus()
    meta["timing_s"]["raw_stage_copy"] = t_copy
    prep_dir = None
    if _PREP is not None:
        prep_dir, tp, cp = _prep_search(store, run_name)
        t1, c1 = time.time(), _rus()
        meta["timing_s"]["search_prep"] = tp
        meta["cpu_s"]["search_prep"] = cp
    rec = {"run": run_name, "raw": raw, "store": store, "pid": os.getpid(), "t_start": t0, "t_end": t1,
           "wall_s": t1 - t0, "core_s": c1 - c0, "n_spectra": meta["n_spectra"], "n_ms2": meta["n_ms2"],
           "n_peaks": meta["n_peaks"], "store_bytes": os.path.getsize(store), "raw_bytes": meta["source_size"],
           "timing_s": meta["timing_s"], "cpu_s": meta["cpu_s"], "instrument": meta.get("instrument_model"),
           "prep_dir": prep_dir}
    if ms1_out:
        rec.update(_ms1_stage(store, run_name, ms1_out))
        rec["t_end_ms1"] = time.time()
    return rec


def _worker_warm() -> dict:
    """No-op job: forces worker start (+ runtime load) and reports its cost."""
    return {"pid": os.getpid(), "core_s": _rus()}


def _worker_main(wid: int, conn, prep, ms1_cfg=None):
    """Decode worker process: init once, then run jobs sent over its own pipe (no shared queue, so killing one
    hung worker cannot corrupt another's channel)."""
    import traceback
    try:
        _worker_init(prep, ms1_cfg)
        conn.send(("ready", wid, os.getpid(), _rus()))
    except BaseException:
        conn.send(("init_error", wid, traceback.format_exc(), None))
        return
    while True:
        msg = conn.recv()
        if msg is None:
            conn.send(("bye", wid, _rus(), None))
            return
        jid, args = msg
        conn.send(("start", wid, jid, time.time()))
        try:
            rec = _ingest(*args)
            conn.send(("done", wid, jid, rec))
        except BaseException:
            conn.send(("error", wid, jid, traceback.format_exc()))


class DecodePool:
    """Supervised decode workers (spawn): per-job timeout -> SIGKILL the worker, respawn it, retry the job (at most
    `retries` times); a worker that dies is replaced the same way. Jobs are dispatched in the given (LPT) order.
    PRIDE-scale input contains corrupt/odd files, and M2b hit sporadic RawFileReader hangs, so one stuck file must
    cost one retry, not the cohort."""

    @staticmethod
    def _proc_cpu(pid):
        try:
            with open(f"/proc/{pid}/stat") as fh:
                f = fh.read().rsplit(")", 1)[1].split()
            return int(f[11]) + int(f[12])
        except Exception:
            return None

    def __init__(self, n, prep=None, timeout_fn=None, retries=2, log=print, stall_s=30.0, ms1_backend="cpu",
                 ms1_gpu_slots=2):
        import multiprocessing as mp
        self.ctx = mp.get_context("spawn")
        self.prep, self.retries, self.log, self.stall_s = prep, retries, log, stall_s
        self.ms1_cfg = None
        if ms1_backend == "gpu":      # E11 opt-in; the default (cpu) leaves the workers exactly as before
            self.ms1_cfg = {"backend": "gpu", "gate": self.ctx.Semaphore(max(1, int(ms1_gpu_slots)))}
        self.timeout_fn = timeout_fn or (lambda job: 120.0)
        self.w: dict[int, dict] = {}
        self.events = {"timeouts": 0, "stalls": 0, "worker_deaths": 0, "respawns": 0, "job_errors": 0}
        self.init_core_s = []
        for i in range(n):
            self._spawn(i)

    def _spawn(self, wid):
        a, b = self.ctx.Pipe(duplex=True)
        args = (wid, b, self.prep) if self.ms1_cfg is None else (wid, b, self.prep, self.ms1_cfg)
        pr = self.ctx.Process(target=_worker_main, args=args, daemon=True)
        pr.start()
        b.close()
        self.w[wid] = {"proc": pr, "conn": a, "job": None, "t0": None, "ready": False}

    def _kill(self, wid, why):
        st = self.w[wid]
        try:
            st["proc"].kill()
            st["proc"].join(5)
        except Exception:
            pass
        try:
            st["conn"].close()
        except Exception:
            pass
        self.log(f"[cohort] decode worker {wid} {why}; respawning", file=sys.stderr)
        self.events["respawns"] += 1
        job = st["job"]
        self._spawn(wid)
        return job

    def run(self, jobs: list[tuple[int, tuple, dict]]):
        """jobs: (jid, args, meta). Yields (meta, rec | Exception) as jobs complete."""
        from collections import deque
        from multiprocessing.connection import wait
        pending = deque((jid, args, meta, 0) for jid, args, meta in jobs)
        info = {jid: (args, meta) for jid, args, meta in jobs}
        tries = {jid: 0 for jid, _, _ in jobs}
        left = len(jobs)
        while left:
            for wid, st in self.w.items():
                if st["ready"] and st["job"] is None and pending:
                    jid, args, meta, _ = pending.popleft()
                    st["job"], st["t0"] = jid, None
                    st["conn"].send((jid, args))
            conns = {st["conn"]: wid for wid, st in self.w.items()}
            sents = {st["proc"].sentinel: wid for wid, st in self.w.items()}
            for obj in wait(list(conns) + list(sents), timeout=1.0):
                if obj in sents:
                    wid = sents[obj]
                    if self.w[wid]["proc"].is_alive() or self.w[wid]["conn"].poll():
                        continue
                    self.events["worker_deaths"] += 1
                    jid = self._kill(wid, "died")
                    if jid is not None:
                        r = self._retry(jid, tries, pending, info, "worker died")
                        if r is not None:
                            left -= 1
                            yield r
                    continue
                wid = conns[obj]
                try:
                    msg = obj.recv()
                except (EOFError, OSError):
                    continue
                kind = msg[0]
                if kind == "ready":
                    self.w[wid]["ready"] = True
                    self.init_core_s.append(msg[3])
                elif kind == "init_error":
                    raise RuntimeError(f"decode worker init failed:\n{msg[2]}")
                elif kind == "start":
                    self.w[wid]["t0"] = time.time()
                    self.w[wid]["cpu"] = (self._proc_cpu(self.w[wid]["proc"].pid), time.time())
                elif kind in ("done", "error"):
                    jid = msg[2]
                    self.w[wid]["job"] = None
                    if kind == "done":
                        msg[3]["attempts"] = tries[jid] + 1
                        left -= 1
                        yield info[jid][1], msg[3]
                    else:
                        self.events["job_errors"] += 1
                        r = self._retry(jid, tries, pending, info, msg[3])
                        if r is not None:
                            left -= 1
                            yield r
            now = time.time()
            for wid, st in list(self.w.items()):
                if st["job"] is not None and st["t0"] is not None:
                    jid = st["job"]
                    # stall: no CPU progress for stall_s (the sporadic RawFileReader hang is a futex wait at 0% CPU)
                    c = self._proc_cpu(st["proc"].pid)
                    c0, tc = st.get("cpu", (None, now))
                    if c is not None and c != c0:
                        st["cpu"] = (c, now)
                        stalled = False
                    else:
                        stalled = now - tc > self.stall_s
                    if stalled or now - st["t0"] > self.timeout_fn(info[jid][1]):
                        self.events["stalls" if stalled else "timeouts"] += 1
                        self._kill(wid, f"{'stalled (no CPU progress)' if stalled else 'timed out'} on "
                                        f"{info[jid][1].get('run')} after {now - st['t0']:.0f}s")
                        r = self._retry(jid, tries, pending, info, "timeout")
                        if r is not None:
                            left -= 1
                            yield r

    def _retry(self, jid, tries, pending, info, why):
        tries[jid] += 1
        if tries[jid] > self.retries:
            return info[jid][1], RuntimeError(f"ingest failed for {info[jid][1].get('run')} after "
                                              f"{tries[jid]} attempts: {why}")
        pending.appendleft((jid, info[jid][0], info[jid][1], tries[jid]))
        return None

    def shutdown(self):
        for st in self.w.values():
            try:
                st["conn"].send(None)
            except Exception:
                pass
        self.final_core_s = []
        for st in self.w.values():
            try:
                while st["conn"].poll(5):
                    m = st["conn"].recv()
                    if m[0] == "bye":
                        self.final_core_s.append(m[2])
                        break
            except Exception:
                pass
        for st in self.w.values():
            st["proc"].join(10)
            if st["proc"].is_alive():
                st["proc"].kill()
                st["proc"].join(5)


# ------------------------------------------------------------------------------------------- spectra loading
SEARCH_COLUMNS = ["scan", "index", "ms_level", "rt_min", "precursor_mz", "precursor_charge", "mz", "intensity"]


def load_ms2(path: str):
    """MS2 Spectra from a local store; same arrays as vulcan.gpu.engine.load_store but reads only the search
    columns and gathers MS2 peaks per record batch in one vectorized pass (no Table.filter + re-flatten)."""
    import numpy as np
    import pyarrow as pa

    from vulcan.gpu.host import Spectra
    from vulcan.io import store
    r = store.open_store(path)
    scans, idxs, rts, pmzs, zs, mzs, its, lens = [], [], [], [], [], [], [], []
    for i in range(r.num_record_batches):
        b = r.get_batch(i)
        lvl = b.column("ms_level").to_numpy()
        sel = np.nonzero(lvl == 2)[0]
        if len(sel) == 0:
            continue
        scans.append(b.column("scan").to_numpy()[sel])
        idxs.append(b.column("index").to_numpy()[sel])
        rts.append(np.asarray(b.column("rt_min").to_numpy(), np.float64)[sel])
        pm = b.column("precursor_mz")
        pmzs.append(np.asarray(pm.to_numpy(zero_copy_only=False), np.float64)[sel])
        z = b.column("precursor_charge")
        zs.append(np.nan_to_num(np.asarray(z.to_numpy(zero_copy_only=False), np.float64), nan=0)[sel].astype(np.int8))
        mzc, itc = b.column("mz"), b.column("intensity")
        off = mzc.offsets.to_numpy().astype(np.int64)
        st, ln = off[sel], off[sel + 1] - off[sel]
        tot = int(ln.sum())
        cl = np.cumsum(ln) - ln
        g = np.repeat(st - cl, ln) + np.arange(tot, dtype=np.int64)
        mzs.append(mzc.values.to_numpy()[g])
        its.append(itc.values.to_numpy()[g])
        lens.append(ln)
    ln = np.concatenate(lens) if lens else np.zeros(0, np.int64)
    off = np.zeros(len(ln) + 1, np.int64)
    np.cumsum(ln, out=off[1:])
    cat = lambda xs, dt: np.concatenate(xs) if xs else np.zeros(0, dt)  # noqa: E731
    name = store.run_meta(path).get("run_name", path)
    del pa
    return Spectra(name, cat(scans, np.uint32), cat(idxs, np.uint32), cat(rts, np.float64), cat(pmzs, np.float64),
                   cat(zs, np.int8), cat(mzs, np.float32), cat(its, np.float32), off)


# ------------------------------------------------------------------------------------------- GPU utilization
class GpuSampler(threading.Thread):
    """NVML (ctypes, no pynvml dependency) utilization.gpu / memory.used sampler."""

    def __init__(self, period_s=0.1, device=0):
        super().__init__(daemon=True)
        self.period, self.samples, self._stop = period_s, [], threading.Event()
        self.ok = False
        try:
            self.nv = ctypes.CDLL("libnvidia-ml.so.1")
            assert self.nv.nvmlInit_v2() == 0
            self.h = ctypes.c_void_p()
            assert self.nv.nvmlDeviceGetHandleByIndex_v2(device, ctypes.byref(self.h)) == 0
            self.ok = True
        except Exception:
            pass

    def run(self):
        if not self.ok:
            return

        class U(ctypes.Structure):
            _fields_ = [("gpu", ctypes.c_uint), ("memory", ctypes.c_uint)]

        class M(ctypes.Structure):
            _fields_ = [("total", ctypes.c_ulonglong), ("free", ctypes.c_ulonglong), ("used", ctypes.c_ulonglong)]
        u, m = U(), M()
        while not self._stop.is_set():
            if self.nv.nvmlDeviceGetUtilizationRates(self.h, ctypes.byref(u)) == 0:
                self.nv.nvmlDeviceGetMemoryInfo(self.h, ctypes.byref(m))
                self.samples.append((time.time(), u.gpu, m.used / 2 ** 30))
            time.sleep(self.period)

    def stop(self):
        self._stop.set()


class CpuSampler(threading.Thread):
    """Host-wide busy cores from /proc/stat (the node is shared: other agents' load is included)."""

    def __init__(self, period_s=0.5):
        super().__init__(daemon=True)
        self.period, self.samples, self._stop = period_s, [], threading.Event()

    @staticmethod
    def _read():
        with open("/proc/stat") as fh:
            v = [int(x) for x in fh.readline().split()[1:]]
        idle = v[3] + v[4]
        return sum(v), idle

    def run(self):
        hz = os.sysconf("SC_CLK_TCK")
        a = self._read()
        while not self._stop.is_set():
            time.sleep(self.period)
            b = self._read()
            busy = ((b[0] - a[0]) - (b[1] - a[1])) / hz / self.period
            self.samples.append((time.time(), busy))
            a = b

    def stop(self):
        self._stop.set()


# ------------------------------------------------------------------------------------------- driver
def plan_jobs(raws: list[str], repeat_to: int | None, stage: str) -> list[dict]:
    """Unique run names (repeats get _rNN), LPT order (largest RAW first)."""
    base = list(raws)
    jobs = []
    names: dict[str, int] = {}
    n = max(len(base), repeat_to or 0)
    for i in range(n):
        raw = base[i % len(base)]
        nm = os.path.basename(raw).rsplit(".", 1)[0]
        k = names.get(nm, 0)
        names[nm] = k + 1
        run = nm if k == 0 else f"{nm}_r{k:02d}"
        jobs.append({"raw": raw, "run": run, "repeat": k, "size": os.path.getsize(raw),
                     "store": os.path.join(stage, run + ".vspec.arrow")})
    jobs.sort(key=lambda j: (-j["size"], j["run"]))
    return jobs


def run_cohort(jobs, index_dir, out_dir, workers, stage, keep_stores=False, torch_threads=4, loader_threads=2,
               params=None, family_top2=False, family_groups=None, write_psms=True, ms1=False, log=print,
               watchdog_s=900, stage_raw=True, search_threads=1, hyp_in_loader=True, worker_prep=True,
               skip_failed=False, scoring=None, chimeric=False, ms1_backend="cpu", ms1_gpu_slots=2):
    """skip_failed: a run whose ingest fails after all retries, or whose search raises, is recorded in
    out["failed"] and skipped instead of aborting the cohort (repository streams contain broken/odd RAWs, E9).
    scoring: None, or dict(lib_root=, model_path=, post_workers=) -> production scoring per run
    (vulcan.rescore.cohort_scoring.CohortScorer); not with family_top2.
    chimeric (E7b): E7 chimeric search + scorer v2_sep per run.  Needs ms1 (the decode workers write the MS1
    hypothesis table) and scoring with the v2 model; precursor groups are built in the loader threads, the secondary
    search + v2 features/model run in the GPU-stage thread, the chimeric competition (per-stratum TDC) + picked
    proteins in the post workers.  Writes OUT/<run>.chim.psms.parquet (+ score) instead of <run>.psms.parquet.
    ms1_backend (E11, opt-in): "cpu" (default, the reference) or "gpu" (hills + isotope grouping on the GPU in the
    decode workers, bitwise-equal tables; at most ms1_gpu_slots workers hold GPU MS1 buffers at once)."""
    if ms1_backend not in ("cpu", "gpu"):
        raise ValueError(f"ms1_backend must be cpu or gpu, got {ms1_backend!r}")
    if chimeric and (not ms1 or not scoring or family_top2):
        raise ValueError("chimeric mode needs ms1=True, scoring (v2 model) and no family_top2")
    t0 = time.time()
    failed: list = []
    c0_self, c0_child = _rus(), _rus(resource.RUSAGE_CHILDREN)
    os.makedirs(stage, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    cpu_s, gpu_s = CpuSampler(), GpuSampler()
    cpu_s.start()
    gpu_s.start()
    import multiprocessing as mp
    p0 = dict(params or {})
    prep = {"index": index_dir, "params": p0} if (worker_prep and not family_top2 and
                                                   p0.get("peak_mode", "global") == "global") else None
    tmul = (6.0 if ms1 else 1.0) * float(os.environ.get("VULCAN_COHORT_TIMEOUT_SCALE", "1"))

    def _timeout(job):   # ~5-10x the measured decode time (D0 3.4 GB: ~8 s -> 86 s); hangs are caught by the stall check
        return tmul * (45.0 + 12.0 * job["size"] / 1e9)
    pool = DecodePool(workers, prep=prep, timeout_fn=_timeout, log=log,
                      ms1_backend=ms1_backend if ms1 else "cpu", ms1_gpu_slots=ms1_gpu_slots)
    djobs = [(i, (j["raw"], j["store"], j["run"], out_dir if ms1 else None, stage_raw), j) for i, j in enumerate(jobs)]

    # ---- GPU setup overlaps the first decodes
    ts = time.perf_counter()
    os.environ["VULCAN_THREADS"] = str(torch_threads)
    import numpy as np
    import torch

    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    torch.set_num_threads(torch_threads)
    t_imp = time.perf_counter() - ts
    p = E.SearchParams(**(params or {}))
    torch.zeros(1, device="cuda")
    t_cuda = time.perf_counter() - ts - t_imp
    l1 = L1Index(index_dir)
    t_l1 = time.perf_counter() - ts - t_imp - t_cuda
    gi = E.GpuIndex(l1, p, {"C": 57.021464}, {"M": 15.9949})
    fam = None
    if family_top2:
        from vulcan.gpu import family as F
        fam = F.FamilyIndex.from_l1(l1, gi, family_groups)
    torch.cuda.synchronize()
    setup = {"torch_import_s": t_imp, "cuda_init_s": t_cuda, "l1_load_s": t_l1, "gpu_index_build_s": gi.build_s,
             "setup_total_s": time.perf_counter() - ts, "gpu_index_bytes": gi.device_bytes,
             "n_postings": gi.n_postings, "n_rows": len(l1), "setup_done_at_s": time.time() - t0}
    log(f"[cohort] {len(jobs)} runs, {workers} decode workers; GPU setup {setup['setup_total_s']:.2f}s (import {t_imp:.2f}, cuda {t_cuda:.2f}, "
        f"l1 {t_l1:.2f}, index {gi.build_s:.2f})")
    if scoring:
        from vulcan.rescore.cohort_scoring import REPORT_PSMS
        if fam is not None or p.report_psms != REPORT_PSMS:
            pool.shutdown()
            raise ValueError(f"production scoring needs report_psms={REPORT_PSMS} (the frozen scorer's PSM table) "
                             f"and no family_top2; got report_psms={p.report_psms}: pass "
                             f"--params '{{\"report_psms\": {REPORT_PSMS}}}'")

    import pyarrow.parquet as pq
    loaders = ThreadPoolExecutor(loader_threads)
    writers = ThreadPoolExecutor(2)
    loaded: queue.Queue = queue.Queue()
    runs = {}

    def _load(rec):
        try:
            a = time.time()
            if rec.get("prep_dir"):                 # worker already prepared arrays + hypotheses: mmap only
                sp, rec["_pre_hyp"] = load_prep(rec["prep_dir"])
            else:
                sp = load_ms2(rec["store"])
                if fam is None and p.peak_mode == "global" and hyp_in_loader:
                    rec["_pre_hyp"] = E.pre_hypotheses(sp, gi, p)     # host-side, off the serial GPU path
            if chimeric:                            # E7 precursor groups from the worker's MS1 hypothesis table
                from vulcan.chimeric import search as CS
                b0 = time.time()
                rec["_groups"] = CS.precursor_groups(
                    CS.read_hypotheses(os.path.join(out_dir, f"{rec['run']}.hypotheses.arrow")), sp)
                rec["chim_groups_s"] = time.time() - b0
            rec["load_t"] = (a - t0, time.time() - t0)
            loaded.put((rec, sp))
        except BaseException as e:  # surface in the GPU loop instead of hanging it
            loaded.put((None, e))

    n_jobs = len(jobs)

    def _feed():
        try:
            it = pool.run(djobs)
            for job, rec in it:
                if isinstance(rec, BaseException):
                    if skip_failed:
                        failed.append({"stage": "ingest", "run": job.get("run"), "error": repr(rec)[-1500:]})
                        loaded.put(("skip", rec))
                        continue
                    loaded.put((None, rec))
                    return
                _handle(rec)
        except BaseException as e:  # noqa: BLE001
            loaded.put((None, e))

    def _handle(rec):
        rec["t_start"] -= t0
        rec["t_end"] -= t0
        if "t_end_ms1" in rec:
            rec["t_end_ms1"] -= t0
        runs[rec["run"]] = rec
        loaders.submit(_load, rec)
    feeder = threading.Thread(target=_feed, daemon=True)
    feeder.start()
    scorer = None
    if scoring:      # after the feeder: the first decodes overlap the scorer setup
        try:
            from vulcan.rescore.cohort_scoring import CohortScorer
            scorer = CohortScorer(l1, index_dir, out_dir, log=log, chimeric=chimeric, **scoring)
        except BaseException:
            pool.shutdown()
            raise
        setup["scorer"] = dict(scorer.setup, ready_at_s=time.time() - t0)

    writes = []

    def _write(tbl, path, rec):
        a = time.time()
        if write_psms:
            pq.write_table(tbl, path)
        if not keep_stores:
            for fp in (rec["store"], rec["store"] + ".json"):
                try:
                    os.remove(fp)
                except OSError:
                    pass
        if rec.get("prep_dir"):
            shutil.rmtree(rec["prep_dir"], ignore_errors=True)
        rec["write_t"] = (a - t0, time.time() - t0)

    lock = threading.Lock()
    state = {"taken": 0, "done": 0, "err": None}

    def _gpu_worker(stream):
        """One GPU-stage thread on its own CUDA stream; >1 threads overlap one run's host-side work (numpy, Python)
        with another run's kernels."""
        with torch.cuda.stream(stream):
            while True:
                with lock:
                    if state["taken"] >= n_jobs or state["err"] is not None:
                        return
                    state["taken"] += 1
                a = time.time()
                try:
                    rec, sp = loaded.get(timeout=watchdog_s)
                except queue.Empty:
                    alive = [st_["proc"].is_alive() for st_ in pool.w.values()]
                    state["err"] = RuntimeError(f"cohort watchdog: no run ready for {watchdog_s} s; "
                                                f"decode workers alive: {sum(alive)}/{len(alive)}")
                    return
                if rec is None:
                    state["err"] = sp
                    return
                if isinstance(rec, str) and rec == "skip":
                    continue
                try:
                    rec["gpu_wait_s"] = time.time() - a
                    b = time.time()
                    if fam is not None:      # replicate-reconstructable output instead of the plain PSM table
                        tbl, st = F.search_family(sp, gi, l1, fam, p)
                        out_name = f"{rec['run']}.family_top2.parquet"
                    else:
                        tbl, st, _ = E.search(sp, gi, l1, p, pre_hyp=rec.pop("_pre_hyp", None))
                        out_name = f"{rec['run']}.psms.parquet"
                    if chimeric:
                        from vulcan.chimeric import search as CS
                        tbl, cst = CS.search_chimeric(sp, gi, l1, p, rec.pop("_groups"), primary=tbl)
                        if "secondary" in cst:
                            st.gpu_ms["chim_secondary"] = cst["secondary"]["gpu_s"] * 1e3
                        rec["chim"] = {"n_secondary_groups": cst["n_secondary_groups"], "n_psms": tbl.num_rows,
                                       "secondary_wall_s": cst.get("secondary", {}).get("wall_s")}
                        out_name = f"{rec['run']}.chim.psms.parquet"
                    stream.synchronize()
                    e = time.time()
                    if scorer is not None:   # features + frozen GBM here; competition/proteins in post workers
                        tbl, _ = scorer(rec, tbl, sp)
                        stream.synchronize()
                        rec["score_t"] = (e - t0, time.time() - t0)
                except BaseException as ex:  # noqa: BLE001
                    if skip_failed and isinstance(ex, Exception):
                        failed.append({"stage": "search", "run": rec.get("run"), "error": repr(ex)[-1500:]})
                        runs.pop(rec.get("run"), None)
                        for fp in (rec["store"], rec["store"] + ".json"):
                            if not keep_stores and os.path.exists(fp):
                                os.remove(fp)
                        if rec.get("prep_dir"):
                            shutil.rmtree(rec["prep_dir"], ignore_errors=True)
                        continue
                    state["err"] = ex
                    return
                rec["search_t"] = (b - t0, e - t0)
                rec["search_wall_s"] = e - b
                rec["gpu_ms"] = st.gpu_ms
                rec["h2d_bytes"] = st.h2d_bytes
                rec["n_psms"] = tbl.num_rows
                rec["n_searched"] = st.counts["n_searched"]
                rec["pre_hyp_used"] = st.counts.get("pre_hyp_used")
                rec["gpu_thread"] = threading.get_ident()
                fut = writers.submit(_write, tbl, os.path.join(out_dir, out_name), rec)
                del sp
                with lock:
                    writes.append(fut)
                    state["done"] += 1
                    k = state["done"]
                if k % 25 == 0 or k == n_jobs:
                    log(f"[cohort] {k}/{n_jobs} runs searched at {time.time() - t0:.1f}s", file=sys.stderr)

    gthreads = [threading.Thread(target=_gpu_worker, args=(torch.cuda.Stream(),), daemon=True)
                for _ in range(max(1, search_threads))]
    for th in gthreads:
        th.start()
    for th in gthreads:
        th.join()
    if state["err"] is not None:
        pool.shutdown()
        raise state["err"]
    gpu_busy_ms = sum(sum(r["gpu_ms"].values()) for r in runs.values())
    for w in writes:
        w.result()
    scoring_out = scorer.collect(runs) if scorer is not None else None   # waits for the post workers
    t_end = time.time() - t0
    if scorer is not None:
        scorer.shutdown()
    pool.shutdown()
    loaders.shutdown()
    writers.shutdown()
    cpu_s.stop()
    gpu_s.stop()
    core_self, core_child = _rus() - c0_self, _rus(resource.RUSAGE_CHILDREN) - c0_child
    out = summarize(list(runs.values()), setup, t_end, workers, core_self, core_child, gpu_busy_ms,
                    gpu_s.samples, cpu_s.samples, t0, search_threads)
    out["params"] = p.to_dict()
    out["decode_pool_events"] = pool.events
    out["decode_worker_init_core_s"] = pool.init_core_s
    out["decode_worker_final_core_s"] = getattr(pool, "final_core_s", [])
    out["family_top2"] = bool(family_top2)
    out["ms1_stage"] = bool(ms1)
    out["chimeric"] = bool(chimeric)
    if ms1 and ms1_backend == "gpu":
        rs = list(runs.values())
        out["ms1_backend"] = "gpu"
        out["ms1_gpu_slots"] = int(ms1_gpu_slots)
        out["ms1_gpu"] = {"gpu_s_sum": float(sum(r.get("ms1_gpu_s", 0.0) for r in rs)),
                          "gpu_s_per_run": float(sum(r.get("ms1_gpu_s", 0.0) for r in rs) / max(1, len(rs))),
                          "gate_wait_s_sum": float(sum(r.get("ms1_gpu_gate_wait_s", 0.0) for r in rs)),
                          "n_runs_gpu": sum(r.get("ms1_backend") == "gpu" for r in rs),
                          "n_runs_cpu_fallback": sum(r.get("ms1_backend") == "cpu_fallback" for r in rs)}
    out["stage_raw"] = bool(stage_raw)
    out["search_threads"] = int(search_threads)
    out["hyp_in_loader"] = bool(hyp_in_loader)
    out["worker_prep"] = prep is not None
    out["gpu_max_mem_bytes"] = torch.cuda.max_memory_allocated()
    out["failed"] = failed
    if scoring_out is not None:
        out["scoring"] = scoring_out
    return out


def _pct(x, q):
    import numpy as np
    return float(np.percentile(x, q)) if len(x) else None


def summarize(runs, setup, makespan, workers, core_self, core_child, gpu_busy_ms, gsamp, csamp, t0,
              search_threads=1) -> dict:
    import numpy as np
    runs = sorted(runs, key=lambda r: r["search_t"][1])
    n = len(runs)
    ing_wall = np.array([r["wall_s"] for r in runs])
    ing_core = np.array([r["core_s"] for r in runs])
    ms1_wall = np.array([r.get("ms1_wall_s", 0.0) for r in runs])
    ms1_core = np.array([r.get("ms1_core_s", 0.0) for r in runs])
    srch = np.array([r["search_wall_s"] for r in runs])
    load = np.array([r["load_t"][1] - r["load_t"][0] for r in runs])
    wr = np.array([r["write_t"][1] - r["write_t"][0] for r in runs])
    gms = np.array([sum(r["gpu_ms"].values()) for r in runs])
    done = np.array([r["search_t"][1] for r in runs])
    ing_end = np.sort([r["t_end"] for r in runs])
    n_ms2 = int(sum(r["n_ms2"] for r in runs))
    # steady state: completions between the 10th and 90th percentile of ingest completion
    i10, i90 = int(0.1 * n), max(int(0.9 * n) - 1, int(0.1 * n) + 1)
    steady = (i90 - i10) / max(1e-9, ing_end[min(i90, n - 1)] - ing_end[i10]) * 3600 if n >= 5 else None
    gu = [s[1] for s in gsamp if s[0] - t0 >= setup["setup_done_at_s"]]
    first_ing = min(r["t_start"] for r in runs)
    stage_busy = {
        # fraction of the makespan each stage (with its parallelism) was busy
        "decode_pool": float((ing_wall.sum() + ms1_wall.sum()) / (workers * makespan)),
        "loader": float(load.sum() / (2 * makespan)),
        "gpu_stage_wall": float(srch.sum() / (search_threads * makespan)),
        "gpu_kernels_events": float(gms.sum() / 1e3 / makespan),
        "writer": float(wr.sum() / (2 * makespan)),
    }
    return {
        "n_runs": n, "workers": workers, "makespan_s": makespan,
        "runs_per_hour": n / makespan * 3600, "runs_per_hour_steady": steady,
        "ms2_per_s": n_ms2 / makespan, "n_ms2": n_ms2,
        "time_to_first_result_s": float(done.min()), "first_ingest_start_s": first_ing,
        "setup": setup,
        "core_s": {"main_process": core_self, "decode_workers": core_child, "total": core_self + core_child,
                   "per_run": (core_self + core_child) / n, "per_1e5_ms2": (core_self + core_child) * 1e5 / n_ms2,
                   "decode_jobs_sum": float(ing_core.sum()), "ms1_jobs_sum": float(ms1_core.sum()),
                   "decode_worker_fixed": float(core_child - ing_core.sum() - ms1_core.sum())},
        "ms1_share_of_worker_busy": float(ms1_wall.sum() / max(1e-9, ing_wall.sum() + ms1_wall.sum())),
        "gpu_s": {"events_total": float(gms.sum() / 1e3), "per_run": float(gms.mean() / 1e3),
                  "util_mean_pct_after_setup": float(np.mean(gu)) if gu else None,
                  "util_samples": len(gsamp)},
        "host_busy_cores_mean": float(np.mean([s[1] for s in csamp])) if csamp else None,
        "stage_busy_fraction": stage_busy,
        "bottleneck": max(("decode_pool", "loader", "gpu_stage_wall", "writer"), key=lambda k: stage_busy[k]),
        "per_run_median": {"ingest_wall_s": _pct(ing_wall, 50), "ingest_core_s": _pct(ing_core, 50),
                           "ms1_wall_s": _pct(ms1_wall, 50), "ms1_core_s": _pct(ms1_core, 50),
                           "load_s": _pct(load, 50), "search_wall_s": _pct(srch, 50), "gpu_ms": _pct(gms, 50),
                           "write_s": _pct(wr, 50)},
        "per_run_sum": {"ingest_wall_s": float(ing_wall.sum()), "ingest_core_s": float(ing_core.sum()),
                        "ms1_wall_s": float(ms1_wall.sum()), "ms1_core_s": float(ms1_core.sum()),
                        "load_s": float(load.sum()), "search_wall_s": float(srch.sum()),
                        "gpu_s": float(gms.sum() / 1e3), "write_s": float(wr.sum())},
        "runs": runs,
        "gpu_util_samples": [(round(s[0] - t0, 2), s[1], round(s[2], 2)) for s in gsamp],
        "host_busy_samples": [(round(s[0] - t0, 2), round(s[1], 1)) for s in csamp],
    }


def main(argv=None):
    import signal
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))   # run finally/atexit: pool workers get reaped
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True)
    ap.add_argument("--raws", nargs="*", default=[])
    ap.add_argument("--raw-list", help="file with one RAW path per line")
    ap.add_argument("--repeat-to", type=int, default=None, help="cycle the list to this many runs (unique names)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--stage", default=f"/tmp/vulcan_cohort_{os.getuid()}/stage")
    ap.add_argument("--workers", type=int, default=None, help="decode workers (default: affinity - 8)")
    ap.add_argument("--torch-threads", type=int, default=4)
    ap.add_argument("--loader-threads", type=int, default=2)
    ap.add_argument("--search-threads", type=int, default=1, help="GPU-stage threads (own CUDA stream each)")
    ap.add_argument("--no-worker-prep", action="store_true",
                    help="flatten MS2 + build hypotheses in the GPU process loaders instead of the decode workers")
    ap.add_argument("--hyp-in-gpu-thread", action="store_true",
                    help="build precursor hypotheses in the GPU-stage thread (E2 behaviour) instead of the loader")
    ap.add_argument("--keep-stores", action="store_true")
    ap.add_argument("--no-psms", action="store_true", help="do not write PSM tables (pure throughput)")
    ap.add_argument("--family-top2", action="store_true")
    ap.add_argument("--no-stage-raw", action="store_true", help="decode RAW in place on Lustre (mmap; not advised)")
    ap.add_argument("--ms1", action="store_true", help="optional MS1 stage (features + hypotheses) in the workers")
    ap.add_argument("--params", default="{}")
    ap.add_argument("--stats", default=None, help="stats JSON path (default OUT/cohort_stats.json)")
    ap.add_argument("--skip-failed", action="store_true",
                    help="record and skip runs whose ingest/search fails instead of aborting the cohort")
    ap.add_argument("--score-lib", default=None,
                    help="production scoring: consolidated predicted-library root (scripts/scoring_speed/build_library.py)")
    ap.add_argument("--score-model", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                                                          "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl"))
    ap.add_argument("--score-post-workers", type=int, default=2)
    ap.add_argument("--score-impl", choices=["fast", "frozen"], default="fast",
                    help="frozen = the pre-optimization scoring components (comparison only)")
    ap.add_argument("--score-placement", choices=["gpu", "host", "mmap"], default="gpu",
                    help="predicted library in GPU memory (default), process memory, or page cache only")
    ap.add_argument("--chimeric", action="store_true",
                    help="E7 chimeric search + scorer v2_sep (implies --ms1; needs --score-lib; default model: v2)")
    ap.add_argument("--ms1-backend", choices=["cpu", "gpu"], default="cpu",
                    help="E11: MS1 hills + isotope grouping on the CPU (default, reference) or the GPU (opt-in; "
                         "bitwise-equal feature/hypothesis tables)")
    ap.add_argument("--ms1-gpu-slots", type=int, default=2,
                    help="max decode workers running GPU MS1 at once (bounds GPU memory, ~4.5 GiB per slot for an 80M-peak "
                         "run + 0.4 GiB CUDA context per worker; a run that finds too little free memory uses the "
                         "bitwise-equal CPU path; --ms1-backend gpu)")
    a = ap.parse_args(argv)
    if a.chimeric:
        a.ms1 = True
        if not a.score_lib:
            ap.error("--chimeric needs --score-lib")
        if "--score-model" not in " ".join(sys.argv if argv is None else argv):
            a.score_model = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                                         "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl")
    age_at_main = process_age_s()
    raws = list(a.raws)
    if a.raw_list:
        raws += [x.strip() for x in open(a.raw_list) if x.strip() and not x.startswith("#")]
    from vulcan.hw import fingerprint, n_threads
    workers = a.workers or max(1, n_threads() - 8)
    jobs = plan_jobs(raws, a.repeat_to, a.stage)
    out = run_cohort(jobs, a.index, a.out, workers, a.stage, keep_stores=a.keep_stores,
                     torch_threads=a.torch_threads, loader_threads=a.loader_threads, params=json.loads(a.params),
                     family_top2=a.family_top2, write_psms=not a.no_psms, ms1=a.ms1,
                     stage_raw=not a.no_stage_raw, search_threads=a.search_threads,
                     hyp_in_loader=not a.hyp_in_gpu_thread, worker_prep=not a.no_worker_prep,
                     skip_failed=a.skip_failed, chimeric=a.chimeric, ms1_backend=a.ms1_backend,
                     ms1_gpu_slots=a.ms1_gpu_slots,
                     scoring=None if not a.score_lib else {"lib_root": a.score_lib, "model_path": a.score_model,
                                                           "post_workers": a.score_post_workers,
                                                           "impl": a.score_impl, "placement": a.score_placement})
    out["process"] = {"age_at_main_s": age_at_main, "argv": sys.argv, "python": sys.executable,
                      "age_at_end_s": process_age_s()}
    out["hw"] = fingerprint()
    path = a.stats or os.path.join(a.out, "cohort_stats.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("n_runs", "workers", "makespan_s", "runs_per_hour", "runs_per_hour_steady",
                                          "bottleneck", "stage_busy_fraction", "core_s", "gpu_s")}, default=float))
    try:
        shutil.rmtree(a.stage) if not a.keep_stores and not os.listdir(a.stage) else None
    except OSError:
        pass


if __name__ == "__main__":
    main()
