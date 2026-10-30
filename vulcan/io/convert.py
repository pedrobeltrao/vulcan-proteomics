"""CLI: python -m vulcan.io.convert <in> <out> [options]

<in>:  Thermo .raw  -> direct in-process reader (default) or --via mzml (TRFP -> mzML -> parser)
       .mzML        -> streaming mzML parser
       Bruker .d    -> timsrust (--bruker-mode frames|spectra)
<out>: a .vspec.arrow file, or a directory (file named <run>.vspec.arrow)

Prints the run metadata JSON (counts, timings) on stdout.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

from . import schema as S
from .store import configure_threads

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _out_path(out: str, inp: str) -> str:
    if out.endswith(S.FILE_SUFFIX):
        return out
    base = os.path.basename(inp.rstrip("/"))
    for ext in (".raw", ".RAW", ".mzML", ".mzml", ".d"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return os.path.join(out, base + S.FILE_SUFFIX)


def trfp(raw: str, outdir: str, fmt: str = "2", extra=()) -> float:
    exe = shutil.which("ThermoRawFileParser") or os.path.join(ROOT, "conda_envs", "trfp", "bin", "ThermoRawFileParser")
    t = time.perf_counter()
    subprocess.run([exe, "-i", raw, "-o", outdir, "-f", fmt, "-l", "3", *extra], check=True)
    return time.perf_counter() - t


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m vulcan.io.convert", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--via", choices=["direct", "mzml"], default="direct", help="Thermo path (default direct)")
    ap.add_argument("--compression", choices=["none", "lz4", "zstd"], default="none")
    ap.add_argument("--run-name")
    ap.add_argument("--run-id", type=int)
    ap.add_argument("--bruker-mode", choices=["frames", "frames_pyo3", "spectra"], default="frames")
    ap.add_argument("--stage-dir", default=None,
                    help="copy the raw input here first (local disk), delete afterwards")
    ap.add_argument("--no-prefetch", action="store_true", help="Thermo direct: skip the sequential page-cache prefetch")
    ap.add_argument("--no-filter-string", action="store_true")
    ap.add_argument("--keep-intermediate", action="store_true")
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args(argv)
    configure_threads(a.threads)
    comp = None if a.compression == "none" else a.compression
    inp = a.input.rstrip("/")
    out = _out_path(a.output, inp)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    staged = None
    t_stage = 0.0
    if a.stage_dir:
        os.makedirs(a.stage_dir, exist_ok=True)
        staged = os.path.join(a.stage_dir, os.path.basename(inp))
        t = time.perf_counter()
        (shutil.copytree if os.path.isdir(inp) else shutil.copyfile)(inp, staged)
        t_stage = time.perf_counter() - t
        inp = staged
    low = inp.lower()
    try:
        if low.endswith(".raw"):
            if a.via == "direct":
                from . import thermo_direct
                meta = thermo_direct.convert(inp, out, run_name=a.run_name, run_id=a.run_id, compression=comp,
                                             include_filter=not a.no_filter_string, do_prefetch=not a.no_prefetch)
            else:
                from . import mzml
                tmpdir = a.stage_dir or os.path.dirname(os.path.abspath(out))
                t_trfp = trfp(inp, tmpdir)
                mz = os.path.join(tmpdir, os.path.basename(inp)[:-4] + ".mzML")
                meta = mzml.convert(mz, out, run_name=a.run_name, run_id=a.run_id, compression=comp, source_raw=inp)
                meta["timing_s"]["trfp"] = t_trfp
                if not a.keep_intermediate:
                    os.remove(mz)
        elif low.endswith(".mzml"):
            from . import mzml
            meta = mzml.convert(inp, out, run_name=a.run_name, run_id=a.run_id, compression=comp)
        elif low.endswith(".d") and os.path.isdir(inp):
            from . import bruker
            meta = bruker.convert(inp, out, mode=a.bruker_mode, run_name=a.run_name, run_id=a.run_id, compression=comp)
        elif low.endswith(".d.zip"):
            ap.error("unzip the .d.zip first (e.g. to a local staging dir); zip members are read-only inputs")
        else:
            ap.error(f"unrecognised input: {inp}")
    finally:
        if staged:
            (shutil.rmtree if os.path.isdir(staged) else os.remove)(staged)
    meta.setdefault("timing_s", {})["stage_copy"] = t_stage
    meta["store_path"] = os.path.abspath(out)
    meta["store_bytes"] = os.path.getsize(out)
    meta.pop("im_lut", None)
    json.dump(meta, sys.stdout, default=str)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
