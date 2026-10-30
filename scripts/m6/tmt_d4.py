#!/usr/bin/env python3
"""M6 TMT on D4 (PXD016766, TMT11 SPS-MS3, IT-CID MS2, 12 fractions of ONE plex, 11 human cell lines).

Search: GPU closed search, human+contaminants, CAM + TMT-K static, TMT peptide-N-term static (added at load time:
mass32 + TMT on every row and mask bit 63 with the kernel N-term delta, i.e. exactly Sage's static `^`), Met-ox
variable, 10 ppm precursor, +-600 ppm fragment (the GPU kernels take ppm only; ~0.4 Da at m/z 700 for the IT-CID
MS2).  Scoring: plain hyperscore paired TDC (no IT-CID static model exists, M3), pooled over the 12 fractions.
Quant: vulcan.quant.tmt (SPS-MS3 linkage, per-plex tolerance, MS1-feature isolation purity, impurity NNLS).

    python scripts/m6/tmt_d4.py search|quant
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

D4 = os.path.join(C.WORK, "d4")
PARAMS = {"report_psms": 1, "fragment_ppm": [-600.0, 600.0], "precursor_ppm": [-10.0, 10.0]}


def stores():
    return sorted(glob.glob(os.path.join(D4, "store", "*.vspec.arrow")))


def tmt_l1():
    from vulcan.index.build import L1Index
    l1 = L1Index(C.D4_INDEX)
    l1.mass32 = (l1.mass32.astype(np.float32) + np.float32(C.TMT)).astype(np.float32)
    assert (np.diff(l1.mass32) >= 0).all()
    l1.mask = l1.mask | np.uint64(1 << 63)
    return l1


def do_search():
    import torch
    from vulcan.gpu import engine as E
    t0 = time.perf_counter()
    p = E.SearchParams(**{k: tuple(v) if isinstance(v, list) else v for k, v in PARAMS.items()})
    torch.zeros(1, device="cuda")
    l1 = tmt_l1()
    gi = E.GpuIndex(l1, p, dict(C.STATIC, K=C.TMT), {"M": C.VAR["M"], "^": C.TMT})
    setup = time.perf_counter() - t0
    out = os.path.join(D4, "psms")
    os.makedirs(out, exist_ok=True)
    runs = []
    for s in stores():
        t1 = time.perf_counter()
        sp = E.load_store(s)
        tbl, st, _ = E.search(sp, gi, l1, p)
        pq.write_table(tbl, os.path.join(out, f"{sp.run}.psms.parquet"))
        runs.append({"run": sp.run, "wall_s": time.perf_counter() - t1, "search_wall_s": st.wall_s["total"],
                     "n_spectra": st.counts["n_spectra"], "cpu_core_s": st.counts["cpu_core_s"]})
        print(runs[-1], flush=True)
    json.dump({"setup_s": setup, "params": PARAMS, "runs": runs}, open(os.path.join(out, "stats.json"), "w"), indent=1)


def do_quant():
    from vulcan.io import store as ST
    from vulcan.quant import tmt as TM
    t_all = time.perf_counter()
    ids = pd.read_parquet(os.path.join(D4, "ids", "ids.parquet"))
    pep = pd.read_parquet(os.path.join(D4, "ids", "peptides.parquet")).set_index("pf_id")
    acc = ids[ids["pass_0.01"] & ~ids["is_decoy"]].copy()
    A = TM.impurity_matrix()
    rows, calib, timing = [], {}, {"extract_s": 0.0, "link_s": 0.0, "purity_s": 0.0, "correct_s": 0.0}
    for s in stores():
        run = os.path.basename(s)[:-len(".vspec.arrow")]
        t0 = time.perf_counter()
        t = ST.read_table(s)
        meta = t.select(["index", "scan", "ms_level", "master_scan", "filter_string"]).to_pandas()
        lk = TM.link_ms3(meta)
        ms2_index = meta.loc[meta["ms_level"] == 2, "index"].to_numpy()
        timing["link_s"] += time.perf_counter() - t0
        a = acc[acc["run"] == run].copy()
        a["store_index"] = ms2_index[a["spec_idx"].to_numpy()]
        a = a.merge(lk[["index", "rep_index", "link_ok"]], left_on="store_index", right_on="index", how="left")
        t0 = time.perf_counter()
        pk = ST.peaks_numpy(t)
        ok = a["link_ok"].fillna(False).to_numpy().astype(bool)
        rep = a["rep_index"].fillna(-1).to_numpy().astype(np.int64)
        # first pass (+-20 ppm) -> per-plex offset/tolerance -> final extraction
        sel = np.flatnonzero(ok)
        I0, E0 = TM.extract(pk.mz, pk.intensity, pk.offsets, rep[sel], 20.0)
        off, tol, cal = TM.calibrate_tolerance(E0, I0)
        I, Er = TM.extract(pk.mz, pk.intensity, pk.offsets, rep[sel], tol, off)
        timing["extract_s"] += time.perf_counter() - t0
        calib[run] = cal
        t0 = time.perf_counter()
        Icor = TM.correct(I, A)
        timing["correct_s"] += time.perf_counter() - t0
        # isolation purity: MS1 feature of the MS2 precursor (hypothesis table)
        t0 = time.perf_counter()
        import pyarrow as pa
        import pyarrow.ipc as ipc
        with pa.memory_map(os.path.join(D4, "ms1", run + ".hypotheses.arrow")) as src:
            h = ipc.open_file(src).read_all().select(["spectrum_id", "source", "charge", "neutral_mass",
                                                      "isolation_fraction", "feature_window_intensity"]).to_pandas()
        h = h[h["source"] == 0]
        run_id = int(t.column("run_id")[0].as_py())
        a["spectrum_id"] = (np.uint64(run_id) << np.uint64(32)) | a["store_index"].to_numpy().astype(np.uint64)
        a["charge"] = a["charge"].astype(np.int8)
        h["charge"] = h["charge"].astype(np.int8)
        m = a[["spectrum_id", "charge", "calcmass"]].reset_index().merge(h, on=["spectrum_id", "charge"])
        m["ppm"] = (m["neutral_mass"] - m["calcmass"]) / m["calcmass"] * 1e6
        # isotope-error tolerant: the hypothesis mass already subtracts e * NEUTRON; keep |ppm| <= 15
        m = m[np.abs(m["ppm"]) <= 15].sort_values(["index", "feature_window_intensity"], ascending=[True, False])
        m = m.drop_duplicates("index")
        pur = np.full(len(a), np.nan)
        pur[m["index"].to_numpy()] = m["isolation_fraction"].to_numpy()
        timing["purity_s"] += time.perf_counter() - t0
        R = pd.DataFrame(np.full((len(a), 11), np.nan), columns=TM.CHANNELS)
        Rc = R.copy()
        R.iloc[sel] = I
        Rc.iloc[sel] = Icor
        out = a[["run", "spec_idx", "scan", "pf_id", "charge", "score", "rep_index", "link_ok"]].copy()
        out["purity"] = pur
        out = pd.concat([out.reset_index(drop=True), R.add_prefix("raw_"), Rc.add_prefix("cor_")], axis=1)
        rows.append(out)
    Q = pd.concat(rows, ignore_index=True)
    Q["modseq"] = pep["modseq"].reindex(Q["pf_id"].to_numpy()).to_numpy()
    Q["proteins"] = pep["proteins"].reindex(Q["pf_id"].to_numpy()).to_numpy()
    Q.to_parquet(os.path.join(D4, "tmt_psms.parquet"))
    raw = Q[[f"raw_{c}" for c in TM.CHANNELS]].to_numpy()
    cor = Q[[f"cor_{c}" for c in TM.CHANNELS]].to_numpy()
    linked = Q["link_ok"].fillna(False).to_numpy().astype(bool)
    allc = linked & (raw > 0).all(axis=1)
    met = {"n_psms_1pct": int(len(Q)), "n_peptidoforms": int(Q["pf_id"].nunique()),
           "ms3_linked_frac": float(linked.mean()),
           "reporter_complete_frac (all 11 > 0, of accepted PSMs)": float(allc.mean()),
           "reporter_any_frac": float((linked & (raw > 0).any(axis=1)).mean()),
           "per_channel_nonzero_frac": {c: float((raw[linked, i] > 0).mean()) for i, c in enumerate(TM.CHANNELS)},
           "tolerance_per_plex": calib,
           "purity": {"n_with_feature": int(np.isfinite(Q["purity"]).sum()),
                      "median": float(np.nanmedian(Q["purity"])),
                      "frac_ge_0.5": float((Q["purity"] >= 0.5).mean()),
                      "frac_ge_0.7": float((Q["purity"] >= 0.7).mean())},
           "impurity_matrix": "OpenMS TMT11 default (lot-specific sheet not deposited)"}
    # technical reproducibility of within-plex channel profiles (no replicate channels in this design)
    share = cor / np.where(cor.sum(axis=1, keepdims=True) > 0, cor.sum(axis=1, keepdims=True), np.nan)
    S = pd.DataFrame(share, columns=TM.CHANNELS)
    S["modseq"] = Q["modseq"].to_numpy()
    S["run"] = Q["run"].to_numpy()
    S["purity"] = Q["purity"].to_numpy()
    for tag, mask in (("all", allc), ("purity_ge_0.5", allc & (Q["purity"].to_numpy() >= 0.5))):
        X = S[mask]
        g = X.groupby("modseq")[TM.CHANNELS]
        n = g.size()
        sd = g.std(ddof=1)[n >= 2]
        mu = g.mean()[n >= 2]
        cv = (sd / mu).to_numpy().ravel()
        # across fractions: peptide channel profile per fraction (sum of PSM shares), then CV across fractions
        F = X.groupby(["modseq", "run"])[TM.CHANNELS].mean()
        nf = F.groupby(level=0).size()
        cvf = (F.groupby(level=0).std(ddof=1) / F.groupby(level=0).mean())[nf >= 2].to_numpy().ravel()
        met[f"cv_{tag}"] = {"n_psms": int(mask.sum()), "peptides_ge2_psms": int((n >= 2).sum()),
                            "median_cv_channel_share_psm_level": float(np.nanmedian(cv)),
                            "peptides_in_ge2_fractions": int((nf >= 2).sum()),
                            "median_cv_channel_share_across_fractions": float(np.nanmedian(cvf)) if len(cvf) else None}
    # ratio compression vs purity: spread (max/min channel share) for pure vs impure spectra
    spread = np.log2(np.nanmax(share, axis=1) / np.clip(np.nanmin(share, axis=1), 1e-6, None))
    pq_ = Q["purity"].to_numpy()
    met["log2_maxmin_channel_spread_by_purity"] = {
        f"{lo}-{hi}": float(np.nanmedian(spread[allc & (pq_ >= lo) & (pq_ < hi)]))
        for lo, hi in ((0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01))}
    met["timing_s"] = dict(timing, total_s=time.perf_counter() - t_all)
    json.dump(met, open(os.path.join(D4, "tmt_metrics.json"), "w"), indent=1, default=float)
    print(json.dumps(met, indent=1, default=float)[:3000])


if __name__ == "__main__":
    {"search": do_search, "quant": do_quant}[sys.argv[1]]()
