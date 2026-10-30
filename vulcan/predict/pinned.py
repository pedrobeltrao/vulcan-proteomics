"""Pinned, content-addressed local store of Prosit (Koina) predictions (tree 3c follow-up; LESSONS: remote
predictors are not bit-deterministic, so every prediction is fetched ONCE and every consumer reads the stored floats).

Key of one entry  = (Koina model, fragmentation type, collision energy, precursor charge, peptide string), where the
peptide string is exactly the ProForma-like string sent to Koina (vulcan.predict.koina.KoinaProsit._strings: UNIMOD
tokens for CAM, Met-ox, TMT).  Models without a CE / fragmentation / charge input carry "na" in that field, so one
table = one Koina call signature:
    ms2  "<model>|frag=<CID|na>|ce=<35.0|na>|z=<z>"   value = the raw float32 'intensities' output (29 x 2 x 3 = 174)
    rt   "<model>|frag=na|ce=na|z=na"                  value = the raw float32 'irt' output (1)
Values are stored exactly as Koina returned them (float32, incl. the -1 "impossible ion" markers); the conversion to
the APD fragment layout is KoinaProsit.predict's own code, so it is a deterministic function of the stored floats.

Layout under <root> (default data/cache/koina_pinned; persistent, gitignored; rebuilt by scripts/m3c/pinned_cache.py):
    MANIFEST.json                 {"tables": {table: {"shards": [{"file", "sha256", "n"}], "n", "content_sha256"}}}
    shards/<sha256[:32]>.npz      key (u64, blake2b-64 of the peptide string), pep (bytes), val (float32 [n, d])
Shards are immutable and named by the sha256 of their bytes; a reader loads exactly the shards listed in MANIFEST.json
and verifies every sha256 (so the manifest pins the snapshot).  content_sha256 of a table = sha256 over its entries
sorted by key (key bytes, peptide bytes, value bytes): identical for any two stores holding the same predictions, so
it identifies the snapshot independently of how the shards were cut.  Appends (mode "fill") take an exclusive flock,
re-read the manifest and never add a key that is already present (a key is predicted at most once, ever).
Mode "frozen" refuses to query Koina: a missing key raises.
"""
from __future__ import annotations

import fcntl
import hashlib
import io
import json
import os
import time

import numpy as np

from .koina import KoinaProsit

ROOT_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                            "data", "cache", "koina_pinned")
NA = "na"


def pep_key(peps) -> np.ndarray:
    return np.array([int.from_bytes(hashlib.blake2b(p.encode() if isinstance(p, str) else p, digest_size=8).digest(),
                                    "little") for p in peps], dtype=np.uint64)


def table_name(model: str, frag=NA, ce=NA, z=NA) -> str:
    ce = NA if ce in (None, NA) else f"{float(ce):.1f}"
    return f"{model}|frag={frag or NA}|ce={ce}|z={NA if z in (None, NA) else int(z)}"


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


class _Tab:
    def __init__(self, key, pep, val):
        o = np.argsort(key, kind="stable")
        self.key, self.pep, self.val = key[o], pep[o], val[o]
        if len(self.key) > 1 and (np.diff(self.key) == 0).any():
            raise RuntimeError("duplicate key in pinned table")

    def find(self, keys):
        if not len(self.key):
            return np.full(len(keys), -1, np.int64)
        p = np.minimum(np.searchsorted(self.key, keys), len(self.key) - 1)
        hit = self.key[p] == keys
        return np.where(hit, p, -1)


class PinnedStore:
    def __init__(self, root: str = ROOT_DEFAULT, mode: str = "frozen"):
        assert mode in ("frozen", "fill")
        self.root, self.mode = root, mode
        os.makedirs(os.path.join(root, "shards"), exist_ok=True)
        self.tabs: dict = {}
        self.loaded: dict = {}
        self.manifest = self._read_manifest()
        self.stats = {"hits": 0, "filled": 0, "koina_s": 0.0, "load_s": 0.0}

    # -------------------------------------------------------------------------------------------- manifest / shards
    def _read_manifest(self):
        p = os.path.join(self.root, "MANIFEST.json")
        return json.load(open(p)) if os.path.exists(p) else {"version": 1, "tables": {}}

    def snapshot_id(self) -> str:
        """sha256 over the per-table content digests (requires rehash() after the last fill)."""
        t = self.manifest["tables"]
        if any(t[k]["content_sha256"] is None for k in t):
            raise RuntimeError("pinned store has un-hashed tables: run PinnedStore(root, 'fill').rehash()")
        return _sha(json.dumps({k: t[k]["content_sha256"] for k in sorted(t)}, sort_keys=True).encode())

    def _load_shard(self, ent):
        b = open(os.path.join(self.root, "shards", ent["file"]), "rb").read()
        if _sha(b) != ent["sha256"]:
            raise RuntimeError(f"pinned shard {ent['file']} sha256 mismatch")
        z = np.load(io.BytesIO(b))
        return z["key"], z["pep"], z["val"]

    def table(self, name):
        if name not in self.tabs:
            t0 = time.perf_counter()
            ents = self.manifest["tables"].get(name, {}).get("shards", [])
            parts = [self._load_shard(e) for e in ents]
            if parts:
                w = max(p[1].dtype.itemsize for p in parts)
                tab = _Tab(np.concatenate([p[0] for p in parts]), np.concatenate([p[1].astype(f"S{w}") for p in parts]),
                           np.concatenate([p[2] for p in parts]))
            else:
                tab = _Tab(np.zeros(0, np.uint64), np.zeros(0, "S1"), np.zeros((0,), np.float32))
            self.tabs[name] = tab
            self.loaded[name] = [e["file"] for e in ents]
            self.stats["load_s"] += time.perf_counter() - t0
        return self.tabs[name]

    @staticmethod
    def content_sha256(tab: _Tab) -> str:
        """sha256 over the entries sorted by key: all keys (u64 LE), then the peptide strings joined by newlines,
        then all values (float32 LE, row-major).  Independent of shard cut and of the fixed-width string dtype."""
        h = hashlib.sha256()
        h.update(np.ascontiguousarray(tab.key, "<u8").tobytes())
        h.update(b"\n".join(bytes(p) for p in tab.pep))
        h.update(np.ascontiguousarray(tab.val, "<f4").tobytes())
        return h.hexdigest()

    def rehash(self):
        """Recompute every table's content_sha256 (under the lock)."""
        with open(os.path.join(self.root, ".lock"), "a") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            try:
                self.manifest = self._read_manifest()
                self.tabs.clear()
                for name, t in self.manifest["tables"].items():
                    tab = self.table(name)
                    t["n"], t["content_sha256"] = int(len(tab.key)), self.content_sha256(tab)
                self._write_manifest()
            finally:
                fcntl.flock(lk, fcntl.LOCK_UN)

    def _write_manifest(self):
        mp = os.path.join(self.root, "MANIFEST.json")
        json.dump(self.manifest, open(mp + ".partial", "w"), indent=1, sort_keys=True)
        os.replace(mp + ".partial", mp)

    def verify(self) -> dict:
        out = {}
        for name in sorted(self.manifest["tables"]):
            self.tabs.pop(name, None)
            tab = self.table(name)                       # verifies every shard sha256
            c = self.content_sha256(tab)
            out[name] = {"n": int(len(tab.key)), "content_sha256": c,
                         "content_ok": c == self.manifest["tables"][name]["content_sha256"]}
        return out

    # -------------------------------------------------------------------------------------------- lookup / fill
    def get(self, name, peps, fetch=None):
        """peps: list[str] -> float32 [n, d].  Missing keys: mode fill -> fetch(missing peps) -> values (queried once,
        appended under the lock); mode frozen -> KeyError."""
        keys = pep_key(peps)
        tab = self.table(name)
        pos = tab.find(keys)
        miss = np.flatnonzero(pos < 0)
        if len(miss):
            if self.mode != "fill" or fetch is None:
                raise KeyError(f"{len(miss)} peptides not in the pinned store table {name} (frozen); e.g. "
                               f"{peps[int(miss[0])]}")
            self._fill(name, [peps[i] for i in miss], fetch)
            tab = self.table(name)
            pos = tab.find(keys)
            assert (pos >= 0).all()
        hit = tab.pep[pos]
        if not all(bytes(h) == p.encode() for h, p in zip(hit, peps)):
            raise RuntimeError("blake2b-64 key collision in pinned store")
        self.stats["hits"] += int(len(peps) - len(miss))
        return tab.val[pos]

    def _fill(self, name, peps, fetch):
        lockp = os.path.join(self.root, ".lock")
        with open(lockp, "a") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            try:
                self.manifest = self._read_manifest()      # another process may have appended meanwhile
                files = [e["file"] for e in self.manifest["tables"].get(name, {}).get("shards", [])]
                if self.loaded.get(name) != files:
                    self.tabs.pop(name, None)
                tab = self.table(name)
                u = sorted(set(peps))
                keys = pep_key(u)
                keep = tab.find(keys) < 0
                new = [p for p, k in zip(u, keep) if k]
                if not new:
                    return
                t0 = time.perf_counter()
                val = np.asarray(fetch(new), np.float32)
                self.stats["koina_s"] += time.perf_counter() - t0
                assert len(val) == len(new)
                key = pep_key(new)
                o = np.argsort(key, kind="stable")
                pep = np.array([p.encode() for p in new], dtype=bytes)[o]
                buf = io.BytesIO()
                np.savez_compressed(buf, key=key[o], pep=pep, val=val[o])
                b = buf.getvalue()
                sha = _sha(b)
                fn = sha[:32] + ".npz"
                tmp = os.path.join(self.root, "shards", "." + fn)
                open(tmp, "wb").write(b)
                os.replace(tmp, os.path.join(self.root, "shards", fn))
                t = self.manifest["tables"].setdefault(name, {"shards": [], "n": 0, "content_sha256": None})
                t["shards"].append({"file": fn, "sha256": sha, "n": len(new),
                                    "added": time.strftime("%Y-%m-%dT%H:%M:%S")})
                w = max(tab.pep.dtype.itemsize, pep.dtype.itemsize)
                vv = val[o] if tab.val.size else val[o]
                if tab.val.size:
                    vv = np.concatenate([tab.val.reshape(len(tab.key), -1), val[o].reshape(len(new), -1)])
                tab = _Tab(np.concatenate([tab.key, key[o]]), np.concatenate([tab.pep.astype(f"S{w}"),
                                                                              pep.astype(f"S{w}")]), vv)
                self.tabs[name] = tab
                self.loaded[name] = files + [fn]
                t["n"] = int(len(tab.key))
                t["content_sha256"] = None       # recomputed by rehash() (O(n); not on every append)
                self._write_manifest()
                self.stats["filled"] += len(new)
            finally:
                fcntl.flock(lk, fcntl.LOCK_UN)


class PinnedKoina(KoinaProsit):
    """KoinaProsit whose remote calls are served from a PinnedStore (Koina is queried only to fill missing keys in
    mode "fill").  predict() / _strings() / the APD conversion are inherited unchanged."""

    def __init__(self, kind, static_mods, var_mods, store: PinnedStore, threads=8, ce=35.0):
        super().__init__(kind, static_mods, var_mods, threads=threads, ce=ce)
        self.store = store
        self.version = dict(self.version, pinned_store=store.root)

    def _sig(self, model, ms2):
        if not ms2:
            return dict(model=model)
        if self.kind == "tmt":
            return dict(model=model, frag="CID", ce=self.ce)
        return dict(model=model)

    def _run(self, model, strs, charges=None):
        idx = [i for i, s in enumerate(strs) if s is not None]
        if not idx:
            return []
        ms2 = charges is not None
        out_name = "intensities" if ms2 else "irt"
        d = 174 if ms2 else 1
        res = np.zeros((len(idx), d), np.float32)
        if ms2:
            zz = np.array([int(min(max(charges[i], 1), 6)) for i in idx])
            groups = {int(z): np.flatnonzero(zz == z) for z in np.unique(zz)}
        else:
            groups = {None: np.arange(len(idx))}
        for z, g in groups.items():
            name = table_name(z=z, **self._sig(model, ms2))
            peps = [strs[idx[j]] for j in g]

            def fetch(new, z=z):
                ch = None if z is None else [z] * len(new)
                got = KoinaProsit._run(self, model, list(new), ch)
                v = np.zeros((len(new), d), np.float32)
                for bi, o in got:
                    v[bi] = np.asarray(o[out_name]["data"], np.float32).reshape(len(bi), d)
                return v
            res[g] = self.store.get(name, peps, fetch).reshape(len(g), d)
        return [(idx, {out_name: {"data": res.reshape(-1)}})]


class PinnedPredCache:
    """Drop-in for vulcan.predict.cache.PredCache (lookup_rt / lookup_ms2 / key / stats) that reads every prediction
    from the pinned store (no pf_id layer, no node-local state).  One instance per consumer is fine: all instances
    over one store return identical floats."""

    def __init__(self, store: PinnedStore, kind: str, static_mods: dict, var_mods: dict, ce: float = 35.0):
        self.pred = PinnedKoina(kind, static_mods, var_mods, store, ce=ce)
        self.store = store
        self.key = "pinned:" + kind
        self.stats = {"predict_s": 0.0, "rt_rows": 0, "ms2_rows": 0}

    # Per-instance memo keyed by (index object, row[, charge]): values come from the pinned store, so memoizing
    # cannot change any float; it only avoids re-converting the same peptidoform within / across runs.
    def _memo(self, l1, kind):
        m = self.__dict__.setdefault("_memos", {})
        k = (id(l1), kind)
        if k not in m:
            m[k] = (l1, {})                       # keep l1 referenced so id() stays unique
        return m[k][1]

    def lookup_rt(self, l1, rows):
        rows = np.asarray(rows)
        t0 = time.perf_counter()
        memo = self._memo(l1, "rt")
        u = np.unique(rows)
        todo = np.array([r for r in u.tolist() if r not in memo], np.int64)
        if len(todo):
            o = self.pred.predict(l1.residues, l1.seq_off[todo], l1.length[todo], l1.mask[todo].view(np.int64), None,
                                  self.pred.ce, "IT", want_ms2=False, want_rt=True)
            memo.update(zip(todo.tolist(), o["rt"].tolist()))
        out = np.array([memo[r] for r in rows.tolist()], np.float32)
        self.stats["predict_s"] += time.perf_counter() - t0
        self.stats["rt_rows"] += len(rows)
        return out

    def lookup_ms2(self, l1, rows, charge):
        rows = np.asarray(rows)
        charge = np.asarray(charge, np.int64)
        t0 = time.perf_counter()
        memo = self._memo(l1, "ms2")
        keys = list(zip(rows.tolist(), charge.tolist()))
        todo = sorted({k for k in keys if k not in memo})
        if todo:
            tr = np.array([k[0] for k in todo], np.int64)
            tz = np.array([k[1] for k in todo], np.int64)
            o = self.pred.predict(l1.residues, l1.seq_off[tr], l1.length[tr], l1.mask[tr].view(np.int64), tz,
                                  self.pred.ce, "IT", want_ms2=True, want_rt=False)
            fo, fl = o["ms2_off"], o["ms2"]
            for i, k in enumerate(todo):
                memo[k] = fl[fo[i]:fo[i + 1]]
        lens = (l1.length[rows].astype(np.int64) - 1) * 4
        off = np.zeros(len(rows) + 1, np.int64)
        np.cumsum(lens, out=off[1:])
        flat = np.concatenate([memo[k] for k in keys]) if keys else np.zeros(0, np.float16)
        assert len(flat) == off[-1]
        self.stats["predict_s"] += time.perf_counter() - t0
        self.stats["ms2_rows"] += len(rows)
        return off, flat
