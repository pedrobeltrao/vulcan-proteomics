"""Derive E1 metrics (core-s per 1e5 MS2 per path/stage, throughput, bytes/spectrum, read-back) from runs.jsonl.

usage: python scripts/ingest/e1_summarize.py results/2026-10-26_E1_ingest [more runs.jsonl ...]
Writes <dir>/metrics.json and prints a compact table.
"""
import json
import os
import sys

THRESH = 30.0  # DESIGN §10 E1: core-s per 1e5 MS2


def load(paths):
    runs = []
    for p in paths:
        with open(p) as fh:
            runs += [json.loads(line) for line in fh if line.strip()]
    by = {}
    for r in runs:
        name, k = r["name"], 2
        while name in by:  # keep replicates (e.g. direct_conc2 measured in the main and a follow-up run)
            name = f"{r['name']}_rep{k}"
            k += 1
        r["name"] = name
        by[name] = r
    return runs, by


def per1e5(core, n_ms2):
    return round(core / n_ms2 * 1e5, 2) if n_ms2 else None


def main():
    out_dir = sys.argv[1]
    import glob
    paths = [os.path.join(out_dir, "runs.jsonl")] + sorted(glob.glob(os.path.join(out_dir, "runs_followup*.jsonl"))) + sys.argv[2:]
    runs, by = load(paths)
    M = {"experiment": "E1 ingest", "date": "2026-10-26", "threshold_core_s_per_1e5_ms2": THRESH,
         "hardware": "4 pinned cores AMD EPYC 7742 (shared node, cgroup cpus 112-115), local /tmp NVMe-backed sda2, Lustre shared-fs",
         "gpu_s": 0.0, "h2d_bytes": 0, "entrapment": False, "measured_from": "raw", "thermo": {}, "bruker": {},
         "concurrency": {}, "lustre_vs_tmp": {}, "readback": {}, "notes": []}
    nms2 = {}
    # ---- Thermo direct
    for k in ("01A", "01B", "XL", "S4"):
        r = by.get(f"direct_{k}")
        if not r or not r.get("result"):
            continue
        res = r["result"]
        nms2[k] = res["n_ms2"]
        cpu = res.get("cpu_s", {})
        stage = {s: round(v, 2) for s, v in cpu.items()}
        stage["process_overhead(imports,.NET init,exit)"] = round(r["core_s"] - sum(cpu.values()), 2)
        M["thermo"][f"direct_{k}"] = dict(
            file=k, n_spectra=res["n_spectra"], n_ms1=res["n_ms1"], n_ms2=res["n_ms2"], n_peaks=res["n_peaks"],
            raw_bytes=res["source_size"], store_bytes=res["store_bytes"], wall_s=r["wall_s"], core_s=round(r["core_s"], 2),
            user_s=r["user_s"], sys_s=r["sys_s"], max_rss_mb=round(r["max_rss_kb"] / 1024),
            core_s_per_1e5_ms2=per1e5(r["core_s"], res["n_ms2"]),
            core_s_per_1e5_ms2_excl_raw_read=per1e5(r["core_s"] - cpu.get("prefetch", 0), res["n_ms2"]),
            core_s_per_1e5_spectra=per1e5(r["core_s"], res["n_spectra"]),
            spectra_per_s_wall=round(res["n_spectra"] / r["wall_s"]), spectra_per_core_s=round(res["n_spectra"] / r["core_s"]),
            store_bytes_per_spectrum=round(res["store_bytes"] / res["n_spectra"]),
            store_bytes_per_peak=round(res["store_bytes"] / res["n_peaks"], 2),
            stage_cpu_s=stage, stage_wall_s={s: round(v, 2) for s, v in res.get("timing_s", {}).items()},
            cgroup_cpu_stall_s=r.get("cgroup_cpu_stall_s"))
    for name in ("direct_01A_noprefetch_cold", "direct_01A_from_tmp", "direct_01A_lustre_out", "direct_01A_lz4", "direct_01A_nofilter", "stage_cp_01A"):
        r = by.get(name)
        if not r:
            continue
        res = r.get("result") or {}
        M["thermo"][name] = dict(wall_s=r["wall_s"], core_s=round(r["core_s"], 2), user_s=r["user_s"], sys_s=r["sys_s"],
                                 core_s_per_1e5_ms2=per1e5(r["core_s"], nms2.get("01A")),
                                 store_bytes=res.get("store_bytes"), stage_cpu_s=res.get("cpu_s"))
    if "direct_01A_from_tmp" in M["thermo"] and "stage_cp_01A" in M["thermo"]:
        c = M["thermo"]["direct_01A_from_tmp"]["core_s"] + M["thermo"]["stage_cp_01A"]["core_s"]
        M["thermo"]["staged_total_01A"] = dict(core_s=round(c, 2), core_s_per_1e5_ms2=per1e5(c, nms2.get("01A")),
                                               wall_s=M["thermo"]["direct_01A_from_tmp"]["wall_s"] + M["thermo"]["stage_cp_01A"]["wall_s"])
    # ---- TRFP paths
    for k in ("01A", "01B"):
        a, b = by.get(f"trfp_mzml_{k}"), by.get(f"mzml_parse_{k}")
        if a and b and b.get("result"):
            res = b["result"]
            n2 = nms2.get(k, res["n_ms2"])
            tot = a["core_s"] + b["core_s"]
            M["thermo"][f"trfp_mzml_{k}"] = dict(
                n_ms2=n2, mzml_bytes=b.get("mzml_bytes"), raw_to_mzml=dict(wall_s=a["wall_s"], core_s=round(a["core_s"], 2), max_rss_mb=round(a["max_rss_kb"] / 1024)),
                mzml_to_store=dict(wall_s=b["wall_s"], core_s=round(b["core_s"], 2), stage_cpu_s=res.get("cpu_s")),
                total_core_s=round(tot, 2), total_wall_s=round(a["wall_s"] + b["wall_s"], 1),
                core_s_per_1e5_ms2=per1e5(tot, n2), parse_core_s_per_1e5_ms2=per1e5(b["core_s"], n2),
                store_bytes=res["store_bytes"], parse_mb_per_core_s=round(b.get("mzml_bytes", 0) / 1e6 / b["core_s"], 1))
    r = by.get("pyteomics_iter_01A")
    if r and r.get("result"):
        res = r["result"]
        M["thermo"]["pyteomics_reference_01A"] = dict(n_spectra=res["n_spectra"], cpu_s=round(res["cpu_s"], 2),
                                                      spectra_per_core_s=round(res["n_spectra"] / res["cpu_s"]),
                                                      note="first N spectra of the 01A mzML, pyteomics MzML iterator with decoding")
    a, b = by.get("trfp_parquet_01A"), by.get("parquet_read_01A")
    if a and b:
        tot = a["core_s"] + b["core_s"]
        M["thermo"]["trfp_parquet_01A"] = dict(raw_to_parquet=dict(wall_s=a["wall_s"], core_s=round(a["core_s"], 2)),
                                               parquet_to_numpy=dict(wall_s=b["wall_s"], core_s=round(b["core_s"], 2)),
                                               parquet_bytes=b.get("parquet_bytes"), total_core_s=round(tot, 2),
                                               core_s_per_1e5_ms2=per1e5(tot, nms2.get("01A")),
                                               metadata_complete=False)
    r = by.get("trfp_none_01A")
    if r and r["core_s"] < 1:
        M["notes"].append("TRFP -f 4 (no spectra output) exits in 0.3 s without reading spectra, so it gives no raw-read floor; dropped.")
    elif r:
        M["thermo"]["trfp_none_01A"] = dict(wall_s=r["wall_s"], core_s=round(r["core_s"], 2),
                                            core_s_per_1e5_ms2=per1e5(r["core_s"], nms2.get("01A")),
                                            note="TRFP -f 4 (no output): RawFileReader read + TRFP per-scan work, floor of any TRFP path")
    # ---- concurrency
    for name in ("direct_conc2", "direct_conc2_rep2", "direct_conc3", "trfp_mzml_conc2", "trfp_mzml_conc2_rep2"):
        r = by.get(name)
        if not r:
            continue
        if r.get("core_s", 0) < 5:  # failed launch (e.g. missing output dir)
            continue
        files = r.get("files", [])
        n2 = sum(nms2.get(k, 0) for k in files)
        singles = [by.get(f"direct_{k}") if name.startswith("direct") else by.get(f"trfp_mzml_{k}") for k in files]
        s_core = sum(s["core_s"] for s in singles if s)
        s_wall = sum(s["wall_s"] for s in singles if s)
        M["concurrency"][name] = dict(files=files, batch_wall_s=round(r["batch_wall_s"], 1), core_s=round(r["core_s"], 2),
                                      core_s_per_1e5_ms2=per1e5(r["core_s"], n2) if name.startswith("direct") else
                                      per1e5(r["core_s"], n2),
                                      sum_single_core_s=round(s_core, 2), sum_single_wall_s=round(s_wall, 1),
                                      core_s_inflation=round(r["core_s"] / s_core, 3) if s_core else None,
                                      wall_speedup_vs_serial=round(s_wall / r["batch_wall_s"], 2) if s_wall else None,
                                      runs_per_hour=round(len(files) * 3600 / r["batch_wall_s"], 1),
                                      cgroup_cpu_stall_s=r.get("cgroup_cpu_stall_s"),
                                      members=[dict(wall_s=m["wall_s"], core_s=round(m["user_s"] + m["sys_s"], 2)) for m in r["members"]])
    # ---- Bruker
    u = by.get("bruker_unzip")
    if u:
        M["bruker"]["unzip_A9_R1"] = dict(wall_s=u["wall_s"], core_s=round(u["core_s"], 2))
    for name in ("bruker_frames_tmp", "bruker_frames_lustre", "bruker_frames_pyo3_tmp", "bruker_spectra_tmp"):
        r = by.get(name)
        if not r or not r.get("result"):
            continue
        res = r["result"]
        M["bruker"][name] = dict(
            n_rows=res["n_spectra"], n_ms1=res["n_ms1"], n_ms2=res["n_ms2"], n_peaks=res["n_peaks"],
            events_outside_windows=res.get("events_outside_windows"), wall_s=r["wall_s"], core_s=round(r["core_s"], 2),
            max_rss_mb=round(r["max_rss_kb"] / 1024), core_s_per_1e5_ms2_rows=per1e5(r["core_s"], res["n_ms2"]),
            core_s_per_1e8_events=round(r["core_s"] / res["n_peaks"] * 1e8, 2) if res["n_peaks"] else None,
            store_bytes=res["store_bytes"], store_bytes_per_row=round(res["store_bytes"] / max(res["n_spectra"], 1)),
            store_bytes_per_peak=round(res["store_bytes"] / max(res["n_peaks"], 1), 2), stage_cpu_s=res.get("cpu_s"),
            converter=res.get("converter"))
    # ---- Lustre vs /tmp
    a, b = M["thermo"].get("direct_01A"), M["thermo"].get("direct_01A_lustre_out")
    if a and b:
        M["lustre_vs_tmp"]["thermo_direct_01A_store_out"] = dict(tmp=dict(wall_s=a["wall_s"], core_s=a["core_s"]),
                                                                 lustre=dict(wall_s=b["wall_s"], core_s=b["core_s"]))
    a, b = M["bruker"].get("bruker_frames_tmp"), M["bruker"].get("bruker_frames_lustre")
    if a and b:
        M["lustre_vs_tmp"]["bruker_frames_store_out"] = dict(tmp=dict(wall_s=a["wall_s"], core_s=a["core_s"]),
                                                             lustre=dict(wall_s=b["wall_s"], core_s=b["core_s"]))
    a, b, c = by.get("direct_01A"), by.get("direct_01A_noprefetch_cold"), M["thermo"].get("staged_total_01A")
    if a and b:
        M["lustre_vs_tmp"]["thermo_raw_input_01A"] = dict(
            lustre_prefetch=dict(wall_s=a["wall_s"], core_s=round(a["core_s"], 2)),
            lustre_mmap_cold_no_prefetch=dict(wall_s=b["wall_s"], core_s=round(b["core_s"], 2)),
            tmp_staged_incl_copy=c)
    # ---- read-back
    for r in runs:
        if r["name"].startswith("readback_") and r.get("result"):
            res = r["result"]
            M["readback"][r["name"][9:]] = dict(cache=r.get("cache"), process_wall_s=r["wall_s"], process_core_s=round(r["core_s"], 2),
                                                 read_wall_s=round(res["total_wall"], 3), read_cpu_s=round(res["total_cpu"], 3),
                                                 gb_per_s=round(res["gb_per_s"], 3), spectra_per_s=round(res["spectra_per_s"]),
                                                 bytes=res["bytes"], n_spectra=res["n_spectra"], mmap=res["mmap"])
    # ---- decision rule
    cands = {}
    for k, v in M["thermo"].items():
        if k in ("direct_01A", "direct_01B") and isinstance(v, dict) and v.get("n_ms2"):
            cands.setdefault("thermo_direct", []).append(v["core_s_per_1e5_ms2"])
        if k.startswith("trfp_mzml_") and v.get("core_s_per_1e5_ms2"):
            cands.setdefault("thermo_trfp_mzml", []).append(v["core_s_per_1e5_ms2"])
    if "trfp_parquet_01A" in M["thermo"]:
        cands["thermo_trfp_parquet(incomplete metadata)"] = [M["thermo"]["trfp_parquet_01A"]["core_s_per_1e5_ms2"]]
    for k in ("bruker_frames_tmp", "bruker_frames_pyo3_tmp", "bruker_spectra_tmp"):
        if k in M["bruker"]:
            cands[k + "(per 1e5 MS2 rows)"] = [M["bruker"][k]["core_s_per_1e5_ms2_rows"]]
    path_mean = {k: round(sum(v) / len(v), 2) for k, v in cands.items()}
    best = min((k for k in path_mean if k.startswith("thermo") and "incomplete" not in k), key=path_mean.get, default=None)
    M["decision"] = dict(core_s_per_1e5_ms2_by_path=path_mean, lowest_core_s_thermo_path=best,
                         lowest_value=path_mean.get(best), threshold=THRESH,
                         ingest_bound=(path_mean.get(best, 1e9) > THRESH) if best else None)
    with open(os.path.join(out_dir, "metrics.json"), "w") as fh:
        json.dump(M, fh, indent=1)
    print(json.dumps(M["decision"], indent=1))


if __name__ == "__main__":
    main()
