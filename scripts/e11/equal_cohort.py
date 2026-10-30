#!/usr/bin/env python3
"""E11 cohort equivalence (dev only): every run of a GPU-MS1 chimeric cohort vs its base run in a CPU-MS1 reference
cohort (E7b base12, /tmp/vulcan_e7b_<uid>/out/base12, written by the E7b default CPU path).

Per run (<base> or <base>_rNN): the MS1 feature table and hypothesis table (.arrow: every column bitwise + schema
metadata), the chimeric PSM table, the peptide table and the protein table (.parquet: every column bitwise, `run`
excluded because copies are renamed by construction, as in scripts/e7b/equal_cohort.py copies mode).
Any column that differs is listed; the run counts as ok only if nothing differs.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "e7b"))

from equal_cohort import pq_np, same  # noqa: E402
from equal_ms1 import read_arrow, table_equal  # noqa: E402


def identity_only(run_prefix: str, base_prefix: str, rec: dict) -> bool:
    """A renamed copy (<base>_rNN) gets its own run_id (store meta), hence its own spectrum_id = (run_id << 32) | index
    (vulcan.io.schema) and its own run_id in the hypothesis-table metadata.  True iff those are the ONLY differences:
    run_id is one constant per table, spectrum_id's low 32 bits equal the base's and its high bits are that run_id,
    the metadata equals the base's after removing run_id, and no other column differs."""
    import numpy as np
    if set(rec[".features.arrow"]) - {"run_id"} or set(rec[".hypotheses.arrow"]) - {"spectrum_id", "<schema metadata>"}:
        return False
    if any(rec[s] for s in (".chim.psms.parquet", ".peptides.parquet", ".proteins.parquet")):
        return False
    fx, fy = read_arrow(run_prefix + ".features.arrow"), read_arrow(base_prefix + ".features.arrow")
    rid = np.unique(fx.column("run_id").to_numpy())
    if len(rid) > 1 or fx.drop(["run_id"]).equals(fy.drop(["run_id"])) is False:
        return False
    hx, hy = read_arrow(run_prefix + ".hypotheses.arrow"), read_arrow(base_prefix + ".hypotheses.arrow")
    sx = hx.column("spectrum_id").to_numpy()
    sy = hy.column("spectrum_id").to_numpy()
    if len(sx) != len(sy) or not np.array_equal(sx & 0xFFFFFFFF, sy & 0xFFFFFFFF):
        return False
    if len(sx) and len(rid) and not np.all((sx >> np.uint64(32)) == np.uint64(rid[0])):
        return False
    if not hx.drop(["spectrum_id"]).equals(hy.drop(["spectrum_id"])):
        return False
    mx = json.loads(hx.schema.metadata[b"vulcan_ms1_hypotheses"]) if b"vulcan_ms1_hypotheses" in (hx.schema.metadata or {}) \
        else {k.decode(): v.decode() for k, v in (hx.schema.metadata or {}).items()}
    my = json.loads(hy.schema.metadata[b"vulcan_ms1_hypotheses"]) if b"vulcan_ms1_hypotheses" in (hy.schema.metadata or {}) \
        else {k.decode(): v.decode() for k, v in (hy.schema.metadata or {}).items()}
    return _drop_run_id(mx) == _drop_run_id(my)


def _drop_run_id(m):
    """Remove run_id at any depth of a (possibly JSON-string-valued) metadata dict."""
    if isinstance(m, str):
        try:
            m = json.loads(m)
        except ValueError:
            return m
    if isinstance(m, dict):
        return {k: _drop_run_id(v) for k, v in m.items() if k != "run_id"}
    if isinstance(m, list):
        return [_drop_run_id(v) for v in m]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="GPU-MS1 cohort output dir")
    ap.add_argument("--ref", default=f"/tmp/vulcan_e7b_{os.getuid()}/out/base12")
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    res = {"cohort": a.out, "reference_cohort": a.ref, "runs": {}}
    for f in sorted(os.listdir(a.out)):
        if not f.endswith(".chim.psms.parquet"):
            continue
        run = f[:-len(".chim.psms.parquet")]
        base = re.sub(r"_r\d\d$", "", run)
        rec = {"base": base}
        for suf in (".features.arrow", ".hypotheses.arrow"):
            t = table_equal(read_arrow(os.path.join(a.out, run + suf)), read_arrow(os.path.join(a.ref, base + suf)))
            rec[suf] = sorted(t["differing_columns"]) + ([] if t["metadata_equal"] else ["<schema metadata>"])
        for suf in (".chim.psms.parquet", ".peptides.parquet", ".proteins.parquet"):
            rec[suf] = same(pq_np(os.path.join(a.out, run + suf), drop=("run",)),
                            pq_np(os.path.join(a.ref, base + suf), drop=("run",)))
        rec["ok"] = not any(v for k, v in rec.items() if k.startswith("."))
        if not rec["ok"] and run != base:
            rec["identity_only"] = identity_only(os.path.join(a.out, run), os.path.join(a.ref, base), rec)
        res["runs"][run] = rec
    res["n_runs"] = len(res["runs"])
    res["n_bases"] = len({r["base"] for r in res["runs"].values()})
    res["all_ok"] = bool(res["runs"]) and all(r["ok"] for r in res["runs"].values())
    res["n_base_named_runs_bitwise_equal"] = sum(r["ok"] for k, r in res["runs"].items() if k == r["base"])
    res["n_copies_identity_only"] = sum(bool(r.get("identity_only")) for r in res["runs"].values())
    res["all_ok_up_to_run_identity"] = bool(res["runs"]) and all(r["ok"] or r.get("identity_only")
                                                                 for r in res["runs"].values())
    bad = {}
    for r in res["runs"].values():
        for k, v in r.items():
            if k.startswith(".") and v:
                bad.setdefault(k, set()).update(v)
    res["differing_columns_any_run"] = {k: sorted(v) for k, v in bad.items()}
    json.dump(res, open(a.json, "w"), indent=1)
    print("runs", res["n_runs"], "bases", res["n_bases"], "all_ok", res["all_ok"], "base-named equal",
          res["n_base_named_runs_bitwise_equal"], "copies identity-only", res["n_copies_identity_only"],
          "all_ok_up_to_run_identity", res["all_ok_up_to_run_identity"], res["differing_columns_any_run"])


if __name__ == "__main__":
    main()
