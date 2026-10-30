"""Resident, consolidated view of the content-hashed prediction cache (vulcan.predict.cache) for production scoring.

Same on-disk cache (<root>/<key>/<kind>/shard_*.npz, key = sha256(model, instrument, NCE, mods)) and the SAME lookup
semantics as `PredCache`: a peptidoform's value is the first copy in (disk shards in file-name order, then in-process
appends in order).  What changes is speed only:

* `consolidate(cache_dir)` merges every shard of a key once (stable, first copy wins = exactly what
  `PredCache._Table._merge` + `searchsorted(side='left')` returns) into `<key>/consolidated/<kind>/*.npy`
  (ids sorted unique u64; rt f32 | off i64 + inten f16) plus a manifest of the shards it covers.  A process then
  memory-maps it (no 8 s npz merge per NCE per process; shared page cache across processes).  Shards written
  after consolidation are merged at load and searched after it (file-name order is creation order).
* lookups are vectorized (numba gather of the variable-length fragment rows); MS2 misses of all precursor charges
  of one call are predicted in ONE predictor call (per-row charge) instead of one call per charge.  Misses are
  appended to disk as ordinary shards (readable by `PredCache`).

Predicted values for misses depend (at the 1-ulp fp16 level) on the predictor batch composition, as they already do
in `PredCache` (results/2026-10-27_scoring_speed/NOTES.md); hits return the stored bits unchanged.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import threading
import time

import numpy as np

from vulcan.predict.cache import _Table, cache_key

try:
    import numba as _nb

    @_nb.njit(cache=True, nogil=True)
    def _copy_rows_nb(src, src_start, ln, dst_start, out):  # pragma: no cover - jitted
        for i in range(ln.shape[0]):
            s = src_start[i]
            d = dst_start[i]
            for j in range(ln[i]):
                out[d + j] = src[s + j]
except Exception:  # noqa: BLE001
    _copy_rows_nb = None


def _copy_rows(src, src_start, ln, dst_start, out):
    """out[dst_start[i] + j] = src[src_start[i] + j] for j < ln[i] (pure data movement; u16 views of f16)."""
    src_start = np.ascontiguousarray(src_start, np.int64)
    ln = np.ascontiguousarray(ln, np.int64)
    dst_start = np.ascontiguousarray(dst_start, np.int64)
    if _copy_rows_nb is not None:
        _copy_rows_nb(src, src_start, ln, dst_start, out)
        return
    rel = np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln)
    out[np.repeat(dst_start, ln) + rel] = src[np.repeat(src_start, ln) + rel]


def _gather(src_u16, src_off, pos, dst_off, out_u16):
    """Rows `pos` of a (src_off, src) ragged table into the ragged layout dst_off."""
    _copy_rows(src_u16, src_off[pos], dst_off[1:] - dst_off[:-1], dst_off[:-1], out_u16)


def _shards(d):
    return sorted(glob.glob(os.path.join(d, "shard_*.npz")))


def consolidate(cache_dir: str, kinds=None, log=print) -> dict:
    """Write <cache_dir>/consolidated/<kind>/ for every table kind (rt, ms2_z*).  Idempotent; returns manifest."""
    kinds = kinds or sorted(k for k in os.listdir(cache_dir) if k == "rt" or k.startswith("ms2_z"))
    man = {}
    for kind in kinds:
        src = os.path.join(cache_dir, kind)
        files = _shards(src)
        if not files:
            continue
        out = os.path.join(cache_dir, "consolidated", kind)
        mf = os.path.join(out, "manifest.json")
        names = [os.path.basename(f) for f in files]
        if os.path.exists(mf) and json.load(open(mf))["shards"] == names:
            man[kind] = json.load(open(mf))
            continue
        t0 = time.perf_counter()
        t = _Table(src, "rt" if kind == "rt" else "ms2")          # the PredCache merge (stable, file order)
        sg = t.segs[0]
        ids = sg["ids"]
        first = np.r_[True, ids[1:] != ids[:-1]] if len(ids) else np.zeros(0, bool)
        tmp = out + ".partial"
        os.makedirs(tmp, exist_ok=True)
        np.save(os.path.join(tmp, "ids.npy"), ids[first])
        if kind == "rt":
            np.save(os.path.join(tmp, "rt.npy"), sg["rt"][first])
        else:
            keep = np.flatnonzero(first)
            ln = sg["off"][keep + 1] - sg["off"][keep]
            off = np.zeros(len(keep) + 1, np.int64)
            np.cumsum(ln, out=off[1:])
            inten = np.empty(off[-1], np.float16)
            _gather(sg["inten"].view(np.uint16), sg["off"], keep.astype(np.int64), off, inten.view(np.uint16))
            np.save(os.path.join(tmp, "off.npy"), off)
            np.save(os.path.join(tmp, "inten.npy"), inten)
        m = {"shards": names, "n_ids": int(first.sum()), "n_duplicates_dropped": int((~first).sum()),
             "built_s": time.perf_counter() - t0}
        json.dump(m, open(os.path.join(tmp, "manifest.json"), "w"), indent=1)
        if os.path.exists(out):
            import shutil
            shutil.rmtree(out)
        os.replace(tmp, out)
        man[kind] = m
        log(f"[library] consolidated {cache_dir}/{kind}: {m['n_ids']} ids, {len(names)} shards, {m['built_s']:.1f}s")
    return man


class _LibTable:
    """Segments: [consolidated (mmap)] + [merge of shards newer than the consolidation] + in-process appends."""

    def __init__(self, d, kind, cdir, resident=False):
        self.d, self.kind = d, kind
        mm = None if resident else "r"
        self.segs = []
        t0 = time.perf_counter()
        covered = set()
        mf = os.path.join(cdir, "manifest.json")
        if os.path.exists(mf):
            covered = set(json.load(open(mf))["shards"])
            seg = {"ids": np.load(os.path.join(cdir, "ids.npy"), mmap_mode=mm)}
            if kind == "rt":
                seg["rt"] = np.load(os.path.join(cdir, "rt.npy"), mmap_mode=mm)
            else:
                seg["off"] = np.load(os.path.join(cdir, "off.npy"), mmap_mode=mm)
                seg["inten"] = np.load(os.path.join(cdir, "inten.npy"), mmap_mode=mm)
            self.segs.append(seg)
        rest = [f for f in _shards(d) if os.path.basename(f) not in covered]
        if rest:
            if covered and min(os.path.basename(f) for f in rest) < max(covered):
                # an older-named shard appeared after consolidation (concurrent writer): keep PredCache order
                # exactly by not using the consolidation for this table
                self.segs = list(_Table(d, kind).segs)
            else:
                self.segs.append(_merge(kind, [dict(np.load(f)) for f in rest]))
        self.load_s = time.perf_counter() - t0

    def find(self, ids):
        seg = np.full(len(ids), -1, np.int64)
        pos = np.zeros(len(ids), np.int64)
        if not len(ids):
            return seg, pos
        o = np.argsort(ids, kind="stable")              # sorted needles: cache-friendly binary search
        sid = ids[o]
        for k, sg in enumerate(self.segs):
            todo = np.flatnonzero(seg[o] < 0)
            n = len(sg["ids"])
            if not len(todo) or not n:
                continue
            p = np.minimum(np.searchsorted(sg["ids"], sid[todo]), n - 1)
            hit = np.asarray(sg["ids"][p]) == sid[todo]
            seg[o[todo[hit]]] = k
            pos[o[todo[hit]]] = p[hit]
        return seg, pos

    def append(self, ids, **arrays):
        os.makedirs(self.d, exist_ok=True)
        name = f"shard_{time.time_ns():020d}_{os.getpid()}"
        tmp = os.path.join(self.d, "." + name + ".npz")
        np.savez(tmp, pf_id=ids, **arrays)
        os.replace(tmp, os.path.join(self.d, name + ".npz"))
        self.segs.append(_merge(self.kind, [{"pf_id": ids, **arrays}]))


def _merge(kind, parts):
    t = _Table.__new__(_Table)
    t.kind = "rt" if kind == "rt" else "ms2"
    return t._merge(parts)


class LibCache:
    """Drop-in for PredCache (same key, stats, lookup_rt / lookup_ms2 results for hits) over a consolidated cache."""

    def __init__(self, root: str, predictor, nce: float, instrument: str = "QE", static_mods=None, var_mods=None,
                 model_sha: str | None = None, resident: bool = False, device=None):
        """resident: read the consolidated tables into process memory (a persistent worker; immune to page-cache
        eviction) instead of memory-mapping them.  device (e.g. 'cuda'): additionally mirror every table on the GPU
        and serve `lookup_ms2_dev` / `lookup_rt_dev` there (pure data movement: the same bits); host memory then
        only holds the memory map."""
        self.pred = predictor
        self.resident = resident
        self.device = device
        self.gsegs = {}
        self.nce, self.instrument = float(nce), instrument
        self.static = static_mods or {"C": 57.021464}
        self.var = var_mods or {"M": 15.9949}
        sha = model_sha or predictor.version["model_zip_sha256"]
        self.key = cache_key(sha, instrument, nce, self.static, self.var)
        self.dir = os.path.join(root, self.key)
        os.makedirs(self.dir, exist_ok=True)
        kj = os.path.join(self.dir, "key.json")
        if not os.path.exists(kj):
            with open(kj, "w") as fh:
                json.dump({"predictor": predictor.version, "instrument": instrument, "nce": self.nce,
                           "static": self.static, "var": self.var}, fh, indent=1)
        self.tables = {}
        self.stats = {"rt_hits": 0, "rt_miss": 0, "ms2_hits": 0, "ms2_miss": 0, "predict_s": 0.0,
                      "n_predicted_rt": 0, "n_predicted_ms2": 0, "library_load_s": 0.0}
        self.lock = threading.RLock()          # several cohort GPU-stage threads may share one cache

    def _tab(self, kind):
        if kind not in self.tables:
            t = _LibTable(os.path.join(self.dir, kind), "rt" if kind == "rt" else "ms2",
                          os.path.join(self.dir, "consolidated", kind), resident=self.resident)
            self.tables[kind] = t
            self.stats["library_load_s"] += t.load_s
        return self.tables[kind]

    def warm(self, kinds=("rt", "ms2_z2", "ms2_z3", "ms2_z4")):
        """Open (resident: load; device: upload) the tables; returns bytes held (on the device if any)."""
        nb = 0
        for k in kinds:
            if self.device is not None:
                nb += sum(v.numel() * v.element_size() for d in self._gseg(k) for v in d.values())
            else:
                nb += sum(v.nbytes for sg in self._tab(k).segs for v in sg.values())
        return nb

    # ------------------------------------------------------------------ device mirror
    def _gseg(self, kind):
        import torch
        t = self._tab(kind)
        g = self.gsegs.setdefault(kind, [])
        import warnings
        while len(g) < len(t.segs):
            sg = t.segs[len(g)]
            d = {"ids": torch.from_numpy(_okey(np.asarray(sg["ids"]))).to(self.device)}
            with warnings.catch_warnings():             # read-only memory map: copied to the device, never written
                warnings.simplefilter("ignore", UserWarning)
                for k in ("rt", "off", "inten"):
                    if k in sg:
                        d[k] = torch.from_numpy(np.ascontiguousarray(sg[k])).to(self.device)
            g.append(d)
        return g

    def _gfind(self, kind, ids):
        import torch
        g = self._gseg(kind)
        q = torch.from_numpy(_okey(ids)).to(self.device)
        seg = torch.full((len(ids),), -1, dtype=torch.int64, device=self.device)
        pos = torch.zeros(len(ids), dtype=torch.int64, device=self.device)
        for k, d in enumerate(g):
            n = d["ids"].numel()
            if not n or not len(ids):
                continue
            p = torch.searchsorted(d["ids"], q).clamp_max(n - 1)
            hit = (d["ids"][p] == q) & (seg < 0)
            seg = torch.where(hit, k, seg)
            pos = torch.where(hit, p, pos)
        return seg, pos

    def lookup_rt_dev(self, l1, rows):
        """As lookup_rt, served from the device mirror; returns a float32 tensor on the device."""
        import torch
        with self.lock:
            ids = l1.pf_id[rows]
            seg, pos = self._gfind("rt", ids)
            miss = (seg < 0).cpu().numpy()
            if miss.any():                          # host path predicts + appends the misses (and counts them)
                self._lookup_rt(l1, np.asarray(rows)[miss])
                seg, pos = self._gfind("rt", ids)
            self.stats["rt_hits"] += int((~miss).sum())
            out = torch.zeros(len(ids), dtype=torch.float32, device=self.device)
            for k, d in enumerate(self._gseg("rt")):
                m = seg == k
                out[m] = d["rt"][pos[m]]
            return out

    def lookup_ms2_dev(self, l1, rows, charge):
        """As lookup_ms2, served from the device mirror: (off [n+1] numpy, flat float16 tensor on the device)."""
        import torch
        charge = np.asarray(charge, np.int64)
        rows = np.asarray(rows)
        with self.lock:
            n = len(rows)
            lens = (l1.length[rows].astype(np.int64) - 1) * 4
            off = np.zeros(n + 1, np.int64)
            np.cumsum(lens, out=off[1:])
            dev = self.device
            out = torch.zeros(int(off[-1]), dtype=torch.float16, device=dev)
            off_t = torch.from_numpy(off).to(dev)
            lens_t = torch.from_numpy(lens).to(dev)
            groups = []
            miss_idx = []
            for z in np.unique(charge):
                m = np.flatnonzero(charge == z)
                seg, pos = self._gfind(f"ms2_z{int(z)}", l1.pf_id[rows[m]])
                groups.append((int(z), m, seg, pos))
                mm = (seg < 0).cpu().numpy()
                if mm.any():
                    miss_idx.append(m[mm])
            n_miss_rows = 0
            if miss_idx:                     # host path on the missing rows only: predicts + appends (counts misses)
                mi = np.concatenate(miss_idx)
                n_miss_rows = len(mi)
                self._lookup_ms2(l1, rows[mi], charge[mi])
                groups = [(z, m) + self._gfind(f"ms2_z{z}", l1.pf_id[rows[m]]) for z, m, _, _ in groups]
            self.stats["ms2_hits"] += n - n_miss_rows
            for z, m, seg, pos in groups:
                m_t = torch.from_numpy(m).to(dev)
                for k, d in enumerate(self._gseg(f"ms2_z{z}")):
                    mk = torch.nonzero(seg == k).squeeze(1)
                    if not mk.numel():
                        continue
                    p = pos[mk]
                    src = d["off"][p]
                    ln = d["off"][p + 1] - src
                    rr = m_t[mk]
                    if not torch.equal(ln, lens_t[rr]):
                        raise AssertionError("library row length mismatch")
                    tot = int(ln.sum())
                    rep = torch.repeat_interleave(torch.arange(len(mk), device=dev), ln)
                    rel = torch.arange(tot, device=dev) - torch.repeat_interleave(torch.cumsum(ln, 0) - ln, ln)
                    out[off_t[rr][rep] + rel] = d["inten"][src[rep] + rel]
            return off, out

    def lookup_rt(self, l1, rows):
        with self.lock:
            return self._lookup_rt(l1, rows)

    def lookup_ms2(self, l1, rows, charge):
        """(off [n+1], flat f16) in row order, fragment layout of vulcan.predict.apd (as PredCache.lookup_ms2)."""
        with self.lock:
            return self._lookup_ms2(l1, rows, charge)

    def _lookup_rt(self, l1, rows):
        ids = l1.pf_id[rows]
        t = self._tab("rt")
        seg, pos = t.find(ids)
        miss = seg < 0
        self.stats["rt_hits"] += int((~miss).sum())
        if miss.any():
            u, first = np.unique(ids[miss], return_index=True)
            r = rows[miss][first]
            t0 = time.perf_counter()
            o = self.pred.predict(l1.residues, l1.seq_off[r], l1.length[r], l1.mask[r].view(np.int64), None,
                                  self.nce, self.instrument, want_ms2=False, want_rt=True)
            self.stats["predict_s"] += time.perf_counter() - t0
            self.stats["n_predicted_rt"] += len(u)
            self.stats["rt_miss"] += int(miss.sum())
            t.append(u, rt=o["rt"])
            seg, pos = t.find(ids)
            assert (seg >= 0).all()
        out = np.zeros(len(ids), np.float32)
        for k, sg in enumerate(t.segs):
            m = seg == k
            if m.any():
                out[m] = sg["rt"][pos[m]]
        return out

    def _lookup_ms2(self, l1, rows, charge):
        charge = np.asarray(charge, np.int64)
        rows = np.asarray(rows)
        n = len(rows)
        lens = (l1.length[rows].astype(np.int64) - 1) * 4
        off = np.zeros(n + 1, np.int64)
        np.cumsum(lens, out=off[1:])
        out = np.zeros(off[-1], np.float16)
        zs = np.unique(charge)
        found = {}
        miss_rows, miss_z = [], []
        for z in zs:
            m = np.flatnonzero(charge == z)
            t = self._tab(f"ms2_z{int(z)}")
            ids = l1.pf_id[rows[m]]
            seg, pos = t.find(ids)
            miss = seg < 0
            self.stats["ms2_hits"] += int((~miss).sum())
            found[int(z)] = (m, ids, seg, pos)
            if miss.any():
                u, first = np.unique(ids[miss], return_index=True)
                miss_rows.append(rows[m][miss][first])
                miss_z.append(np.full(len(u), z))
                self.stats["ms2_miss"] += int(miss.sum())
        if miss_rows:                                   # one predictor call for every charge's misses
            r = np.concatenate(miss_rows)
            zz = np.concatenate(miss_z)
            t0 = time.perf_counter()
            o = self.pred.predict(l1.residues, l1.seq_off[r], l1.length[r], l1.mask[r].view(np.int64), zz,
                                  self.nce, self.instrument, want_ms2=True, want_rt=False)
            self.stats["predict_s"] += time.perf_counter() - t0
            self.stats["n_predicted_ms2"] += len(r)
            for z in np.unique(zz):
                k = np.flatnonzero(zz == z)
                ln = o["ms2_off"][k + 1] - o["ms2_off"][k]
                so = np.zeros(len(k) + 1, np.int64)
                np.cumsum(ln, out=so[1:])
                buf = np.empty(so[-1], np.float16)
                _gather(o["ms2"].view(np.uint16), o["ms2_off"], k.astype(np.int64), so, buf.view(np.uint16))
                ids_z = l1.pf_id[r[k]]
                t = self._tab(f"ms2_z{int(z)}")
                t.append(ids_z, off=so, inten=buf)
                m, ids, _, _ = found[int(z)]
                seg, pos = t.find(ids)
                assert (seg >= 0).all()
                found[int(z)] = (m, ids, seg, pos)
        ov = out.view(np.uint16)
        for z, (m, ids, seg, pos) in found.items():
            t = self._tab(f"ms2_z{z}")
            for k, sg in enumerate(t.segs):
                mk = np.flatnonzero(seg == k)
                if not len(mk):
                    continue
                p = pos[mk]
                so = np.asarray(sg["off"])
                ln = so[p + 1] - so[p]
                assert np.array_equal(ln, lens[m][mk])
                _copy_rows(np.asarray(sg["inten"]).view(np.uint16), so[p], ln, off[m][mk], ov)
        return off, out


def _okey(ids):
    """u64 ids -> int64 with the same ORDER (sign bit flipped), for torch.searchsorted."""
    return (np.ascontiguousarray(ids, np.uint64) ^ np.uint64(1 << 63)).view(np.int64)


def library_digest(root: str, keys) -> dict:
    """sha256 over every consolidated manifest (shard lists) of the given keys: identifies the library state."""
    h = hashlib.sha256()
    out = {}
    for k in sorted(keys):
        for mf in sorted(glob.glob(os.path.join(root, k, "consolidated", "*", "manifest.json"))):
            b = open(mf, "rb").read()
            h.update(mf[len(root):].encode())
            h.update(b)
            out[os.path.relpath(mf, root)] = json.loads(b)["n_ids"]
    return {"sha256": h.hexdigest(), "n_ids": out}
