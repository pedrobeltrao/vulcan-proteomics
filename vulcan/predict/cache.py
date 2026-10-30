"""Content-hashed prediction library/cache (DESIGN §4 L3: "tryptic libraries are precomputed once per
(proteome, instrument class)"; lazy prediction for anything else, with a measured hit rate).

Keyed by peptidoform id (xxhash64 of the canonical peptidoform string, vulcan.index) -> content-defined, so one
library serves every index built from the same peptides (prod, entrapment arms, decoys).  Layout under
<root>/<key>/, key = sha256(model sha, instrument, NCE, mods)[:16]:
  rt/shard_*.npz      pf_id (u64), rt (f32)                      -- charge independent
  ms2_z<z>/shard_*.npz pf_id (u64), off (i64 n+1), inten (f16)    -- per precursor charge
Shards are append-only; `lookup_*` loads all shards once (sorted merge) and predicts + appends misses.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import time

import numpy as np


def cache_key(model_sha: str, instrument: str, nce: float, static_mods: dict, var_mods: dict) -> str:
    s = json.dumps({"model": model_sha, "instrument": instrument, "nce": round(float(nce), 2),
                    "static": sorted(static_mods.items()), "var": sorted(var_mods.items())}, sort_keys=True)
    return hashlib.sha256(s.encode()).hexdigest()[:16]


class _Table:
    """In-memory view of one (rt | ms2_z) table as a list of sorted segments: the shards on disk (merged once, at
    first use) plus one segment per in-process append (no reload of the library after a lazy append)."""

    def __init__(self, d, kind):
        self.d, self.kind = d, kind
        self.segs = []
        self.load_s = 0.0
        t0 = time.perf_counter()
        files = sorted(glob.glob(os.path.join(self.d, "shard_*.npz")))
        if files:
            self.segs.append(self._merge([dict(np.load(f)) for f in files]))
        self.load_s = time.perf_counter() - t0

    def _merge(self, parts):
        ids = np.concatenate([p["pf_id"] for p in parts])
        o = np.argsort(ids, kind="stable")
        if self.kind == "rt":
            return {"ids": ids[o], "rt": np.concatenate([p["rt"] for p in parts])[o]}
        lens = np.concatenate([np.diff(p["off"]) for p in parts])
        flat = np.concatenate([p["inten"] for p in parts])
        src_off = np.zeros(len(lens) + 1, np.int64)
        np.cumsum(lens, out=src_off[1:])
        l2 = lens[o]
        off = np.zeros(len(l2) + 1, np.int64)
        np.cumsum(l2, out=off[1:])
        gather = np.repeat(src_off[:-1][o] - off[:-1], l2) + np.arange(off[-1])
        return {"ids": ids[o], "off": off, "inten": flat[gather]}

    def find(self, ids):
        """-> (segment index or -1, position) per id."""
        seg = np.full(len(ids), -1, np.int64)
        pos = np.zeros(len(ids), np.int64)
        for k, sg in enumerate(self.segs):
            todo = np.flatnonzero(seg < 0)
            if not len(todo) or not len(sg["ids"]):
                continue
            p = np.minimum(np.searchsorted(sg["ids"], ids[todo]), len(sg["ids"]) - 1)
            hit = sg["ids"][p] == ids[todo]
            seg[todo[hit]] = k
            pos[todo[hit]] = p[hit]
        return seg, pos

    def append(self, ids, **arrays):
        os.makedirs(self.d, exist_ok=True)
        # unique shard name (several processes may append to one cache concurrently; entries are content-defined,
        # so duplicates across shards are harmless: lookups take any copy)
        name = f"shard_{time.time_ns():020d}_{os.getpid()}"
        tmp = os.path.join(self.d, "." + name + ".npz")
        np.savez(tmp, pf_id=ids, **arrays)
        os.replace(tmp, os.path.join(self.d, name + ".npz"))
        self.segs.append(self._merge([{"pf_id": ids, **arrays}]))


class PredCache:
    def __init__(self, root: str, predictor, nce: float, instrument: str = "QE", static_mods=None, var_mods=None):
        self.pred = predictor
        self.nce, self.instrument = float(nce), instrument
        self.static = static_mods or {"C": 57.021464}
        self.var = var_mods or {"M": 15.9949}
        self.key = cache_key(predictor.version["model_zip_sha256"], instrument, nce, self.static, self.var)
        self.dir = os.path.join(root, self.key)
        os.makedirs(self.dir, exist_ok=True)
        with open(os.path.join(self.dir, "key.json"), "w") as fh:
            json.dump({"predictor": predictor.version, "instrument": instrument, "nce": self.nce,
                       "static": self.static, "var": self.var}, fh, indent=1)
        self.tables = {}
        self.stats = {"rt_hits": 0, "rt_miss": 0, "ms2_hits": 0, "ms2_miss": 0, "predict_s": 0.0,
                      "n_predicted_rt": 0, "n_predicted_ms2": 0,
                      "library_load_s": 0.0}

    def _tab(self, kind):
        if kind not in self.tables:
            self.tables[kind] = _Table(os.path.join(self.dir, kind), "rt" if kind == "rt" else "ms2")
            self.stats["library_load_s"] += self.tables[kind].load_s
        return self.tables[kind]

    def lookup_rt(self, l1, rows):
        """rows: index rows (of L1Index l1). Returns rt per row (predicting misses)."""
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
            out[m] = sg["rt"][pos[m]]
        return out

    def lookup_ms2(self, l1, rows, charge):
        """Returns (off [n+1], flat f16 intensities) in row order: fragment layout of vulcan.predict.apd."""
        charge = np.asarray(charge, np.int64)
        n = len(rows)
        lens = (l1.length[rows].astype(np.int64) - 1) * 4
        off = np.zeros(n + 1, np.int64)
        np.cumsum(lens, out=off[1:])
        out = np.zeros(off[-1], np.float16)
        for z in np.unique(charge):
            m = np.flatnonzero(charge == z)
            t = self._tab(f"ms2_z{int(z)}")
            ids = l1.pf_id[rows[m]]
            seg, pos = t.find(ids)
            miss = seg < 0
            self.stats["ms2_hits"] += int((~miss).sum())
            if miss.any():
                u, first = np.unique(ids[miss], return_index=True)
                r = rows[m][miss][first]
                t0 = time.perf_counter()
                o = self.pred.predict(l1.residues, l1.seq_off[r], l1.length[r], l1.mask[r].view(np.int64),
                                      np.full(len(r), z), self.nce, self.instrument, want_ms2=True, want_rt=False)
                self.stats["predict_s"] += time.perf_counter() - t0
                self.stats["n_predicted_ms2"] += len(u)
                self.stats["ms2_miss"] += int(miss.sum())
                t.append(u, off=o["ms2_off"], inten=o["ms2"])
                seg, pos = t.find(ids)
                assert (seg >= 0).all()
            for k, sg in enumerate(t.segs):
                mk = seg == k
                if not mk.any():
                    continue
                src = sg["off"][pos[mk]]
                ln = sg["off"][pos[mk] + 1] - src
                assert np.array_equal(ln, lens[m][mk])
                dst = off[m][mk]
                rel = np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln)
                out[np.repeat(dst, ln) + rel] = sg["inten"][np.repeat(src, ln) + rel]
        return off, out

    def precompute(self, l1, charges=(2, 3), batch_rows=2_000_000):
        """Library precompute for every row of an index (RT + MS2 at `charges`)."""
        rows = np.arange(len(l1))
        t0 = time.perf_counter()
        for b in range(0, len(rows), batch_rows):
            r = rows[b:b + batch_rows]
            self.lookup_rt(l1, r)
            for z in charges:
                self.lookup_ms2(l1, r, np.full(len(r), z))
        return time.perf_counter() - t0
