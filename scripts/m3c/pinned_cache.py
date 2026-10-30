#!/usr/bin/env python3
"""3c follow-up, item 1: the pinned, content-addressed Prosit/Koina prediction store (vulcan.predict.pinned) and the
dev replicate features computed FROM it.

    python scripts/m3c/pinned_cache.py features [--datasets D4 D10cf] [--kinds prod E2001 E2002 E2003 C] [--mode fill]
        IT-CID L3 features for every dev replicate PSM directory of 3c ($VULCAN_3C_WORK/psms/<ds>_<kind>, search output
        unchanged) with every prediction read from the pinned store; written to $VULCAN_3C_FU/psms/<ds>_<kind>/
        (<run>.psms.parquet symlinked, <run>.feat.parquet new).  --mode fill queries Koina once for keys not yet in
        the store (this is how the store is built); --mode frozen refuses any Koina call (re-runs).
    python scripts/m3c/pinned_cache.py verify     shard sha256s + per-table content digests vs MANIFEST.json
    python scripts/m3c/pinned_cache.py rebuild --out DIR
        Rebuild: re-query Koina for every key of the pinned store into a NEW store at DIR (same tables, same keys,
        single fill per table) and report how many entries differ from the pinned floats.  The pinned store is never
        overwritten by this command (Koina is not bit-deterministic; the pinned snapshot is the reference).
    python scripts/m3c/pinned_cache.py manifest --dest FILE   copy MANIFEST.json + snapshot id (for git)
Store root: $VULCAN_PINNED (default data/cache/koina_pinned).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

import numpy as np  # noqa: E402

FU_WORK = os.environ.get("VULCAN_3C_FU", "/tmp/vulcan_3c_fu")
FU_RES = os.path.join(C.RES, "followup")
PINNED = os.environ.get("VULCAN_PINNED", os.path.join(C.ROOT, "data", "cache", "koina_pinned"))


def fu_psm_dir(tag):
    return os.path.join(FU_WORK, "psms", tag)


_stores = {}


def store(mode="frozen", root=PINNED):
    from vulcan.predict.pinned import PinnedStore
    k = (root, mode)
    if k not in _stores:
        _stores[k] = PinnedStore(root, mode)
    return _stores[k]


def builder(family, mode="frozen", root=PINNED):
    """ItcidFeatureBuilder whose cache factory returns a PinnedPredCache over the pinned store."""
    from vulcan.predict.pinned import PinnedPredCache
    from vulcan.rescore.itcid import ItcidFeatureBuilder
    st, var = C.mods(family)
    kind = "tmt" if C.is_tmt(family) else "lf"
    s = store(mode, root)
    caches = {}

    def factory(nce):
        if nce not in caches:
            caches[nce] = PinnedPredCache(s, kind, st, var, ce=nce)
        return caches[nce]
    return ItcidFeatureBuilder(factory, st, var)


def features(datasets, kinds, mode):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from vulcan.gpu.engine import load_store
    sys.path.insert(0, HERE)
    import features as F3
    for ds in datasets:
        fam = C.DEV[ds]["family"]
        fb = builder(fam, mode)
        for k in kinds:
            src = C.psm_dir(f"{ds}_{k}")
            dst = fu_psm_dir(f"{ds}_{k}")
            os.makedirs(dst, exist_ok=True)
            l1 = C.load_l1(fam, k)
            ip = os.path.join(dst, "feat_info.json")
            infos = json.load(open(ip)) if os.path.exists(ip) else {"store": PINNED, "runs": {}}
            for f in sorted(glob.glob(os.path.join(src, "*.psms.parquet"))):
                run = os.path.basename(f)[:-len(".psms.parquet")]
                lp = os.path.join(dst, run + ".psms.parquet")
                if not os.path.lexists(lp):
                    os.symlink(f, lp)
                outp = os.path.join(dst, run + ".feat.parquet")
                if os.path.exists(outp):
                    continue
                t = pq.read_table(f)
                psm = {c: t.column(c).to_numpy() for c in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
                sp = load_store(F3.find_store(run))
                t0 = time.perf_counter()
                feats, info = fb.run(l1, psm, sp)
                info["wall_s"] = time.perf_counter() - t0
                pq.write_table(pa.table(feats), outp + ".partial")
                os.replace(outp + ".partial", outp)
                infos["runs"][run] = info
                json.dump(infos, open(ip, "w"), indent=1, default=float)
                s = store(mode)
                print(ds, k, run, f"{info['wall_s']:.1f}s", "filled", s.stats["filled"], "koina_s",
                      round(s.stats["koina_s"], 1), flush=True)
            open(os.path.join(dst, "stats.json"), "w").write(json.dumps({"src": src, "store": PINNED}))


def compare_with_3c(datasets, kinds):
    """Feature-level comparison of the pinned-store features with the 3c features (fresh Koina, pf_id cache)."""
    import pyarrow.parquet as pq
    out = {}
    for ds in datasets:
        for k in kinds:
            tag = f"{ds}_{k}"
            n = nd = 0
            mx = {}
            for f in sorted(glob.glob(os.path.join(fu_psm_dir(tag), "*.feat.parquet"))):
                a = pq.read_table(f)
                b = pq.read_table(os.path.join(C.psm_dir(tag), os.path.basename(f)))
                for c in a.column_names:
                    x, y = a.column(c).to_numpy(), b.column(c).to_numpy()
                    d = ~((x == y) | (np.isnan(x) & np.isnan(y)))
                    nd += int(d.sum())
                    n += int(d.size)
                    if d.any():
                        mx[c] = max(mx.get(c, 0.0), float(np.nanmax(np.abs(x[d].astype(float) - y[d].astype(float)))))
            out[tag] = {"n_values": n, "n_values_differ_from_3c": nd, "max_abs_diff_by_feature": mx}
            print(tag, out[tag], flush=True)
    return out


def rebuild(out_root):
    from vulcan.predict.koina import KoinaProsit
    from vulcan.predict.pinned import PinnedStore
    ref = PinnedStore(PINNED, "frozen")
    new = PinnedStore(out_root, "fill")
    rep = {"reference_snapshot": ref.snapshot_id(), "tables": {}}
    for name in sorted(ref.manifest["tables"]):
        model = name.split("|")[0]
        z = name.split("z=")[1]
        frag = name.split("frag=")[1].split("|")[0]
        ce = name.split("ce=")[1].split("|")[0]
        kind = "tmt" if "TMT" in model else "lf"
        kp = KoinaProsit(kind, {}, {}, ce=float(ce) if ce != "na" else 35.0)
        tab = ref.table(name)
        peps = [bytes(p).decode() for p in tab.pep]
        d = tab.val.shape[1] if tab.val.ndim == 2 else 1
        out_name = "irt" if z == "na" else "intensities"

        def fetch(ps):
            got = KoinaProsit._run(kp, model, list(ps), None if z == "na" else [int(z)] * len(ps))
            v = np.zeros((len(ps), d), np.float32)
            for bi, o in got:
                v[bi] = np.asarray(o[out_name]["data"], np.float32).reshape(len(bi), d)
            return v
        _ = frag
        t0 = time.perf_counter()
        v = new.get(name, peps, fetch).reshape(len(peps), -1)
        r = tab.val.reshape(len(peps), -1)
        diff = ~np.all(v == r, axis=1)
        rep["tables"][name] = {"n": len(peps), "n_entries_differ": int(diff.sum()),
                               "max_abs_diff": float(np.abs(v - r).max()) if len(peps) else 0.0,
                               "wall_s": time.perf_counter() - t0}
        print(name, rep["tables"][name], flush=True)
    new.rehash()
    rep["rebuilt_snapshot"] = PinnedStore(out_root, "frozen").snapshot_id()
    rep["bit_identical_to_pinned"] = rep["rebuilt_snapshot"] == rep["reference_snapshot"]
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["features", "verify", "rebuild", "manifest", "compare"])
    ap.add_argument("--datasets", nargs="+", default=["D10cf", "D4"])
    ap.add_argument("--kinds", nargs="+", default=list(C.KINDS))
    ap.add_argument("--mode", default="frozen", choices=["frozen", "fill"])
    ap.add_argument("--out")
    ap.add_argument("--dest")
    a = ap.parse_args()
    os.makedirs(FU_RES, exist_ok=True)
    if a.cmd == "features":
        features(a.datasets, a.kinds, a.mode)
    elif a.cmd == "compare":
        r = compare_with_3c(a.datasets, a.kinds)
        json.dump(r, open(a.out or os.path.join(FU_RES, "pinned_vs_3c_features.json"), "w"), indent=1)
    elif a.cmd == "verify":
        s = store("frozen")
        r = {"root": PINNED, "snapshot_id": s.snapshot_id(), "tables": s.verify()}
        r["all_ok"] = all(v["content_ok"] for v in r["tables"].values())
        print(json.dumps(r, indent=1))
        json.dump(r, open(a.out or os.path.join(FU_RES, "pinned_verify.json"), "w"), indent=1)
    elif a.cmd == "rebuild":
        r = rebuild(a.out)
        json.dump(r, open(os.path.join(FU_RES, "pinned_rebuild.json"), "w"), indent=1)
    elif a.cmd == "manifest":
        s = store("frozen")
        shutil.copyfile(os.path.join(PINNED, "MANIFEST.json"), a.dest)
        print(s.snapshot_id())


if __name__ == "__main__":
    main()
