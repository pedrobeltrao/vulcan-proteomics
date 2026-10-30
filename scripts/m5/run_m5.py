#!/usr/bin/env python3
"""M5a driver: cleavable XL (DSSO, D3b) and non-cleavable XL (DSS, D3) on the GPU index, with stage timings.
Dev only: D3b = PXD029252 mainlib_rep2 (not T3), D3 = PXD014337 DSS R1; nothing under data/sealed or
data/truth/sealed is read (the D3b library table/FASTAs are re-fetched public copies in /tmp/vulcan_m5/truth,
sha256-identical to the manifest).

    flock /tmp/vulcan_bench.lock python scripts/m5/run_m5.py --which D3b D3 --out /tmp/vulcan_m5/out

Stages timed per dataset (wall + process core-s; GPU ms from CUDA events inside the search):
  ingest (RAW on NVMe -> store, vulcan.io.convert, separate process), L1 index build, GPU index build,
  load store, XL search (k0, doublets/alpha K1, K2 linear-like / alpha+beta, evidence, pairing, ranking,
  linear/mono/loop competitors, CSM table), link-site localization (CPU), L4 FDR + truth (separate step).
The search is repeated `--repeat` times in-process (first = cold kernels); outputs are written from pass 1 and
checked identical across passes.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time

t_start = time.perf_counter()
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
WORK = "/tmp/vulcan_m5"

DATASETS = {
    "D3b": {"raw": "data/raw/PXD029252/20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2.raw", "linker": "DSSO",
            "fasta": f"{WORK}/fasta/D3b_ribo_human671.fasta"},
    "D3": {"raw": "data/raw/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw", "linker": "DSS",
           "fasta": f"{WORK}/fasta/D3_cas9_contam.fasta"},
}
INDEX_PARAMS = {"enzyme": {"cleave_at": "KR", "restrict": "P", "c_terminal": True, "missed_cleavages": 3,
                           "min_len": 5, "max_len": 40}, "peptide_min_mass": 400.0, "peptide_max_mass": 6000.0,
                "max_variable_mods": 2}
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
SEARCH_PARAMS = {"max_fragment_charge": 3, "max_peaks": 200}


def core_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime, c.ru_utime + c.ru_stime


def timed(fn):
    t0, (c0, _) = time.perf_counter(), core_s()
    r = fn()
    return r, time.perf_counter() - t0, core_s()[0] - c0


def ingest(raw, out_dir):
    """RAW -> store in a child process (as in M2b: copy to NVMe first, never decode on Lustre)."""
    stage = os.path.join(WORK, "rawstage")
    os.makedirs(stage, exist_ok=True)
    loc = os.path.join(stage, os.path.basename(raw))
    t0 = time.perf_counter()
    _, ch0 = core_s()
    shutil.copy(os.path.join(ROOT, raw), loc)
    t_copy = time.perf_counter() - t0
    subprocess.run([os.path.join(ROOT, "scripts/m2b/vpy.sh"), "-m", "vulcan.io.convert", loc, out_dir,
                    "--threads", "4"], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)
    os.remove(loc)
    _, ch1 = core_s()
    return {"copy_to_nvme_s": t_copy, "wall_s": time.perf_counter() - t0, "child_core_s": ch1 - ch0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", nargs="+", default=["D3b", "D3"])
    ap.add_argument("--out", default=f"{WORK}/out")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--skip-ingest", action="store_true")
    a = ap.parse_args()
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu import engine as E
    from vulcan.hw import fingerprint
    from vulcan.index.build import L1Index, build
    from vulcan.xl import linkers as LK
    from vulcan.xl import search as XS
    from vulcan.xl import sites as S
    t_import = time.perf_counter() - t_start
    t0 = time.perf_counter()
    torch.zeros(1, device="cuda")
    t_cuda = time.perf_counter() - t0
    os.makedirs(a.out, exist_ok=True)
    report = {"process": {"python_import_s": t_import, "cuda_init_s": t_cuda}, "datasets": {}}
    p = E.SearchParams(**SEARCH_PARAMS)
    for which in a.which:
        cfg = DATASETS[which]
        linker = LK.LINKERS[cfg["linker"]]
        xp = XS.XLParams()
        R = {"linker": linker.name, "cleavable": linker.cleavable, "search_params": p.to_dict(),
             "xl_params": xp.to_dict(), "index_params": INDEX_PARAMS, "fasta": cfg["fasta"], "raw": cfg["raw"]}
        store_dir = os.path.join(WORK, "store_timed")
        run = os.path.basename(cfg["raw"])[:-4]
        store = os.path.join(store_dir, run + ".vspec.arrow")
        if not a.skip_ingest:
            R["ingest"] = ingest(cfg["raw"], store_dir)
            R["ingest"]["store_meta"] = json.load(open(store + ".json"))
        idx_dir = os.path.join(WORK, "index_timed", which)
        meta, w, c = timed(lambda: build([cfg["fasta"]], idx_dir, INDEX_PARAMS))
        R["index_build"] = {"wall_s": w, "core_s": c, "n_peptidoforms": meta["n_peptidoforms"],
                            "n_targets": meta["n_targets"], "content_hash": meta["content_hash"]}
        l1, w, c = timed(lambda: L1Index(idx_dir))
        R["l1_load"] = {"wall_s": w, "core_s": c}

        def mk_gi():
            g = E.GpuIndex(l1, p, STATIC, VAR)
            torch.cuda.synchronize()
            return g
        gi, w, c = timed(mk_gi)
        R["gpu_index_build"] = {"wall_s": w, "core_s": c, "device_bytes": gi.device_bytes,
                                "n_postings": gi.n_postings}
        sp, w, c = timed(lambda: E.load_store(store))
        R["load_store"] = {"wall_s": w, "core_s": c, "n_ms2": len(sp)}
        passes, first = [], None
        for rep in range(a.repeat):
            def srch():
                out = XS.search(sp, gi, l1, p, xp, linker)
                torch.cuda.synchronize()
                return out
            (tbl, lin, st), w, c = timed(srch)
            passes.append({"wall_s": w, "core_s": c, "stage_wall_s": st.wall_s, "stage_gpu_ms": st.gpu_ms,
                           "counts": st.counts})
            if first is None:
                first = (tbl, lin)
            else:
                passes[-1]["identical_to_pass1"] = bool(tbl.equals(first[0]))
            print(which, rep, json.dumps({"wall": w, "core": c, "csm_rows": tbl.num_rows}), flush=True)
        R["search_passes"] = passes
        tbl, lin = first
        m = (tbl.column("rank").to_numpy() == 1) & tbl.column("xl_wins").to_numpy(zero_copy_only=False)
        t1 = tbl.filter(pa.array(m))
        (sa, sb, na, nb), w, c = timed(lambda: S.localize_csms(t1, sp, l1, linker, STATIC, VAR))
        R["sites"] = {"wall_s": w, "core_s": c, "n_csms": t1.num_rows}
        t1 = t1.append_column("site_a", pa.array(sa)).append_column("site_b", pa.array(sb))
        od = os.path.join(a.out, which)
        os.makedirs(od, exist_ok=True)
        pq.write_table(tbl, os.path.join(od, f"{run}.xl_all.parquet"))
        pq.write_table(t1, os.path.join(od, f"{run}.csm.parquet"))
        if lin is not None:
            pq.write_table(lin, os.path.join(od, f"{run}.competitors.psms.parquet"))
        R["n_rank1_xl_wins"] = int(m.sum())
        R["n_rank1"] = int((tbl.column("rank").to_numpy() == 1).sum())
        R["competitor_kind_of_lost_spectra"] = {
            str(k): int(v) for k, v in zip(*np.unique(
                tbl.column("competitor_kind").to_numpy()[(tbl.column("rank").to_numpy() == 1) & ~m],
                return_counts=True))}
        R["gpu_max_mem_bytes"] = torch.cuda.max_memory_allocated()
        report["datasets"][which] = R
        del gi
        torch.cuda.empty_cache()
    report["process"]["wall_s"] = time.perf_counter() - t_start
    report["process"]["core_s"], report["process"]["child_core_s"] = core_s()
    report["hw"] = fingerprint()
    json.dump(report, open(os.path.join(a.out, "run_m5.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
