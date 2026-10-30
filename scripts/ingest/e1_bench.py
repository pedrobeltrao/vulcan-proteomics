"""E1 ingest benchmark driver (DESIGN.md section 10, E1).

MUST run under the node-wide benchmark lock (bench/README.md):
    flock /tmp/vulcan_bench.lock python scripts/ingest/e1_bench.py --out results/2026-10-26_E1_ingest

Every measured step is a separate process wrapped in /usr/bin/time -v (core-s = user+sys of
that process tree). Raw inputs are evicted from the page cache before "cold" steps
(posix_fadvise DONTNEED), so raw->X steps include reading the RAW from Lustre.
Thread pools are pinned (OMP/MKL/NUMBA/RAYON=1, Arrow=1, .NET processor count 4).
Writes <out>/runs.jsonl incrementally and <out>/metrics_raw.json at the end.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PY = sys.executable
TRFP = os.path.join(ROOT, "conda_envs", "trfp", "bin", "ThermoRawFileParser")
STAGE = "/tmp/vulcan_stage/e1"
STORE = "/tmp/vulcan_store/e1"
LUSTRE_OUT = os.path.join(ROOT, ".e1_lustre_scratch")

F = {
    "01A": os.path.join(ROOT, "data/raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw"),
    "01B": os.path.join(ROOT, "data/raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw"),
    "S4": os.path.join(ROOT, "data/raw/PXD029252/20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw"),
    "XL": os.path.join(ROOT, "data/raw/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"),
}
BRUKER_ZIP = os.path.join(ROOT, "data/raw/ProteoBench-DIA-plasma/A9_G_DIA_nLC_tTOF_R1.d.zip")

ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMBA_NUM_THREADS="1", RAYON_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", VULCAN_THREADS="1", DOTNET_PROCESSOR_COUNT=os.environ.get("DOTNET_PROCESSOR_COUNT", "4"), DOTNET_gcServer="0",
           PYTHONPATH=ROOT)

RUNS = []
OUT = None
RUNS_NAME = "runs.jsonl"


def evict(*paths):
    for p in paths:
        if os.path.isdir(p):
            for dp, _, fs in os.walk(p):
                evict(*[os.path.join(dp, f) for f in fs])
            continue
        if not os.path.exists(p):
            continue
        fd = os.open(p, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def cpu_pressure_us():
    """cgroup 'some' CPU-pressure stall total (us): time any task in our cgroup waited for CPU."""
    try:
        with open("/sys/fs/cgroup/cpu.pressure") as fh:
            return int(fh.readline().split("total=")[1])
    except (OSError, IndexError, ValueError):
        return -1


def _parse_time(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            k, _, v = line.strip().partition(": ")
            d[k] = v
    def hms(s):
        parts = [float(x) for x in s.split(":")]
        return sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
    return dict(user_s=float(d["User time (seconds)"]), sys_s=float(d["System time (seconds)"]),
                wall_s=hms(d["Elapsed (wall clock) time (h:mm:ss or m:ss)"]), max_rss_kb=int(d["Maximum resident set size (kbytes)"]),
                exit=int(d.get("Exit status", 0)))


_HW = None


def _record(rec):
    global _HW
    if _HW is None:
        sys.path.insert(0, ROOT)
        from vulcan import hw
        _HW = hw.fingerprint()
    rec["hw"] = _HW
    rec["core_s"] = rec.get("user_s", 0) + rec.get("sys_s", 0)
    RUNS.append(rec)
    with open(os.path.join(OUT, RUNS_NAME), "a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")
    print(f"[E1] {rec['name']}: wall {rec.get('wall_s', 0):.1f}s core {rec['core_s']:.1f}s", flush=True)


def _launch(cmd, tf):
    return subprocess.Popen(["/usr/bin/time", "-v", "-o", tf, *cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            env=ENV, text=True)


def _finish(p, tf, name, meta):
    so, se = p.communicate()
    rec = dict(name=name, cmd=" ".join(p.args[3:]), **meta, **_parse_time(tf))
    if rec["exit"] != 0:
        rec["stderr"] = se[-2000:]
    try:
        rec["result"] = json.loads(so.strip().splitlines()[-1]) if so.strip() else None
    except (json.JSONDecodeError, IndexError):
        rec["stdout_tail"] = so[-500:]
    return rec


def timed(name, cmd, **meta):
    tf = os.path.join(STAGE, f"time_{name}.txt")
    p0 = cpu_pressure_us()
    rec = _finish(_launch(cmd, tf), tf, name, meta)
    rec["cgroup_cpu_stall_s"] = (cpu_pressure_us() - p0) / 1e6
    _record(rec)
    return rec


def concurrent(name, cmds, **meta):
    t0 = time.perf_counter()
    p0 = cpu_pressure_us()
    ps = [(_launch(c, os.path.join(STAGE, f"time_{name}_{i}.txt")), os.path.join(STAGE, f"time_{name}_{i}.txt"))
          for i, c in enumerate(cmds)]
    recs = [_finish(p, tf, f"{name}[{i}]", meta) for i, (p, tf) in enumerate(ps)]
    wall = time.perf_counter() - t0
    agg = dict(name=name, n_files=len(cmds), batch_wall_s=wall, wall_s=wall,
               user_s=sum(r["user_s"] for r in recs), sys_s=sum(r["sys_s"] for r in recs),
               members=recs, cgroup_cpu_stall_s=(cpu_pressure_us() - p0) / 1e6, **meta)
    _record(agg)
    return agg


def convert_cmd(inp, out, *extra):
    return [PY, "-m", "vulcan.io.convert", inp, out, *extra]


def rm(p):
    if os.path.isdir(p):
        shutil.rmtree(p, ignore_errors=True)
    elif os.path.exists(p):
        os.remove(p)
    for side in (p + ".json",):
        if os.path.exists(side):
            os.remove(side)


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip", default="", help="comma list of sections to skip: direct,conc,trfp_mzml,trfp_parquet,"
                    "trfp_none,trfp_conc,bruker,readback")
    ap.add_argument("--runs-name", default="runs.jsonl", help="append-only run log name (follow-up runs use another name)")
    ap.add_argument("--warm-trfp", action="store_true",
                    help="prefetch the RAW (sequential cat, recorded as its own step) before TRFP: RawFileReader mmap "
                         "faults on cold Lustre inflate TRFP wall time ~2-10x without changing its core-s")
    a = ap.parse_args()
    OUT = os.path.abspath(a.out)
    skip = set(filter(None, a.skip.split(",")))
    os.makedirs(OUT, exist_ok=True)
    for d in (STAGE, STORE, LUSTRE_OUT):
        os.makedirs(d, exist_ok=True)
    global RUNS_NAME
    RUNS_NAME = a.runs_name
    open(os.path.join(OUT, RUNS_NAME), "w").close()
    files = {k: v for k, v in F.items() if os.path.exists(v)}
    print("[E1] files:", files, flush=True)
    t_start = time.time()

    # ---------------- Thermo direct (in-process RawFileReader)
    if "direct" not in skip:
        # warm-up of the .NET/python import path (not recorded; caches site-packages on the Lustre client)
        subprocess.run(convert_cmd(files["S4"] if "S4" in files else files["01A"], os.path.join(STORE, "warmup.vspec.arrow")),
                       env=ENV, capture_output=True)
        rm(os.path.join(STORE, "warmup.vspec.arrow"))
        for k, f in files.items():
            evict(f)
            timed(f"direct_{k}", convert_cmd(f, os.path.join(STORE, f"{k}_direct.vspec.arrow")), file=k, path="direct",
                  raw="lustre_cold_prefetch", out="tmp")
        f = files["01A"]
        evict(f)
        timed("direct_01A_noprefetch_cold", convert_cmd(f, os.path.join(STORE, "01A_direct_nopf.vspec.arrow"), "--no-prefetch"),
              file="01A", path="direct", raw="lustre_cold_noprefetch", out="tmp")
        rm(os.path.join(STORE, "01A_direct_nopf.vspec.arrow"))
        evict(f)
        staged = os.path.join(STAGE, os.path.basename(f))
        timed("stage_cp_01A", ["cp", f, staged], file="01A", path="stage_copy", raw="lustre_cold", out="tmp")
        timed("direct_01A_from_tmp", convert_cmd(staged, os.path.join(STORE, "01A_direct_tmpraw.vspec.arrow"), "--no-prefetch"),
              file="01A", path="direct", raw="tmp_staged", out="tmp")
        rm(staged); rm(os.path.join(STORE, "01A_direct_tmpraw.vspec.arrow"))
        evict(f)
        timed("direct_01A_lustre_out", convert_cmd(f, os.path.join(LUSTRE_OUT, "01A_direct.vspec.arrow")),
              file="01A", path="direct", raw="lustre_cold_prefetch", out="lustre")
        evict(f)
        timed("direct_01A_lz4", convert_cmd(f, os.path.join(STORE, "01A_direct_lz4.vspec.arrow"), "--compression", "lz4"),
              file="01A", path="direct", raw="lustre_cold_prefetch", out="tmp", compression="lz4")
        evict(f)
        timed("direct_01A_nofilter", convert_cmd(f, os.path.join(STORE, "01A_direct_nofilt.vspec.arrow"), "--no-filter-string"),
              file="01A", path="direct", raw="lustre_cold_prefetch", out="tmp", filter_string=False)
        rm(os.path.join(STORE, "01A_direct_nofilt.vspec.arrow"))

    # ---------------- concurrency (direct)
    if "conc" not in skip:
        third = "S4" if "S4" in files else "XL"
        if "direct" in skip and third in files:  # single-file baseline for the third file
            evict(files[third])
            timed(f"direct_{third}", convert_cmd(files[third], os.path.join(STORE, f"{third}_direct.vspec.arrow")), file=third,
                  path="direct", raw="lustre_cold_prefetch", out="tmp")
        for combo in (["01A", "01B"], ["01A", "01B", third]):
            combo = [k for k in combo if k in files]
            if len(combo) < len(set(combo)) or len(combo) < 2 or (len(combo) == 2 and any(r["name"] == "direct_conc2" for r in RUNS)):
                continue
            evict(*[files[k] for k in combo])
            concurrent(f"direct_conc{len(combo)}", [convert_cmd(files[k], os.path.join(STORE, f"{k}_conc.vspec.arrow")) for k in combo],
                       files=combo, path="direct", raw="lustre_cold_prefetch", out="tmp")
            for k in combo:
                rm(os.path.join(STORE, f"{k}_conc.vspec.arrow"))

    # ---------------- TRFP paths
    def raw_ready(name, f):
        evict(f)
        if a.warm_trfp:
            timed(f"prefetch_{name}", ["dd", f"if={f}", "of=/dev/null", "bs=16M", "status=none"], raw="lustre_cold", out="devnull", stage="prefetch")
    if "trfp_mzml" not in skip:
        for k in [k for k in ("01A", "01B") if k in files]:
            f = files[k]
            raw_ready(f"trfp_mzml_{k}", f)
            timed(f"trfp_mzml_{k}", [TRFP, "-i", f, "-o", STAGE, "-f", "2", "-l", "3"], file=k, path="trfp_mzml",
                  stage="raw_to_mzml", raw="lustre_cold", out="tmp")
            mz = os.path.join(STAGE, os.path.basename(f)[:-4] + ".mzML")
            timed(f"mzml_parse_{k}", convert_cmd(mz, os.path.join(STORE, f"{k}_mzml.vspec.arrow"), "--run-name",
                                                 os.path.basename(f)[:-4]),
                  file=k, path="trfp_mzml", stage="mzml_to_store", out="tmp", mzml_bytes=os.path.getsize(mz))
            if k == "01A":
                timed(f"pyteomics_iter_{k}", [PY, os.path.join(ROOT, "scripts/ingest/pyteomics_iter.py"), mz], file=k,
                      path="trfp_mzml", stage="mzml_parse_pyteomics_reference", out="none")
            if k != "01A":
                os.remove(mz)
    if "trfp_parquet" not in skip:
        f = files["01A"]
        raw_ready("trfp_parquet_01A", f)
        timed("trfp_parquet_01A", [TRFP, "-i", f, "-o", STAGE, "-f", "3", "-l", "3"], file="01A", path="trfp_parquet",
              stage="raw_to_parquet", raw="lustre_cold", out="tmp")
        pq = os.path.join(STAGE, os.path.basename(f)[:-4] + ".mzparquet")
        timed("parquet_read_01A", [PY, os.path.join(ROOT, "scripts/ingest/parquet_peaks.py"), pq], file="01A",
              path="trfp_parquet", stage="parquet_to_numpy", out="none", parquet_bytes=os.path.getsize(pq))
        rm(pq)
    if "trfp_none" not in skip:
        f = files["01A"]
        raw_ready("trfp_none_01A", f)
        timed("trfp_none_01A", [TRFP, "-i", f, "-o", STAGE, "-f", "4", "-l", "3"], file="01A", path="trfp_none",
              stage="raw_read_only", raw="lustre_cold", out="none")
    if "trfp_conc" not in skip and "01B" in files:
        os.makedirs(os.path.join(STAGE, "c"), exist_ok=True)
        for k in ("01A", "01B"):
            raw_ready(f"trfp_conc2_{k}", files[k])
        concurrent("trfp_mzml_conc2", [[TRFP, "-i", files[k], "-o", os.path.join(STAGE, "c"), "-f", "2", "-l", "3"]
                                       for k in ("01A", "01B")], files=["01A", "01B"], path="trfp_mzml",
                   stage="raw_to_mzml", raw="lustre_prefetched" if a.warm_trfp else "lustre_cold", out="tmp")
        rm(os.path.join(STAGE, "c"))

    # ---------------- Bruker
    if "bruker" not in skip:
        dd = os.path.join(STAGE, "bruker")
        rm(dd)
        os.makedirs(dd, exist_ok=True)
        evict(BRUKER_ZIP)
        timed("bruker_unzip", [PY, "-c", f"import zipfile; zipfile.ZipFile({BRUKER_ZIP!r}).extractall({dd!r})"],
              file="A9_R1", path="bruker", stage="zip_to_d", raw="lustre_cold", out="tmp")
        d = os.path.join(dd, "A9_G_DIA_nLC_tTOF_R1.d")
        for mode, out in (("frames", "tmp"), ("frames", "lustre"), ("frames_pyo3", "tmp"), ("spectra", "tmp")):
            outp = os.path.join(STORE if out == "tmp" else LUSTRE_OUT, f"A9_R1_{mode}.vspec.arrow")
            timed(f"bruker_{mode}_{out}", convert_cmd(d, outp, "--bruker-mode", mode), file="A9_R1", path=f"bruker_{mode}",
                  raw="tmp_staged_warm", out=out)

    # ---------------- read-back
    if "readback" not in skip:
        rb = os.path.join(ROOT, "scripts/ingest/readback.py")
        targets = [("01A_direct_tmp", os.path.join(STORE, "01A_direct.vspec.arrow")),
                   ("01A_direct_lz4_tmp", os.path.join(STORE, "01A_direct_lz4.vspec.arrow")),
                   ("01A_direct_lustre", os.path.join(LUSTRE_OUT, "01A_direct.vspec.arrow")),
                   ("A9_R1_frames_tmp", os.path.join(STORE, "A9_R1_frames.vspec.arrow")),
                   ("A9_R1_frames_lustre", os.path.join(LUSTRE_OUT, "A9_R1_frames.vspec.arrow"))]
        for name, p in targets:
            if not os.path.exists(p):
                continue
            timed(f"readback_{name}_cold", [PY, rb, p, "--cold"], store=name, cache="cold")
            subprocess.run(["cat", p], stdout=subprocess.DEVNULL)  # warm the page cache
            timed(f"readback_{name}_warm", [PY, rb, p], store=name, cache="warm")
        rm(LUSTRE_OUT)

    with open(os.path.join(OUT, "metrics_raw.json"), "w") as fh:
        json.dump(dict(started=t_start, finished=time.time(), files=files, runs=RUNS), fh, indent=1, default=str)
    print("[E1] done", flush=True)


if __name__ == "__main__":
    main()
