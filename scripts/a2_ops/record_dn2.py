#!/usr/bin/env python3
"""Operator recording for A2.9.4 steps 2-4 (label-free only).

  record_dn2.py select   -> results/<R>/selection_log.tsv (every walk decision + every downloaded file and checksum)
  record_dn2.py admit    -> results/<R>/admission_counts.json, results/<R>/membership.sha256 (per-project membership
                            sha256 of the admit codes + sha256 of every admission file), swap test copied in
  record_dn2.py persist DIR...  copy $W/<DIR> trees to data/cache/a2/<DIR> and write results/<R>/cache_<name>.sha256
Reads only selector/admission outputs (codes, label-free counts).  Never opens a gate-index label file.
"""
from __future__ import annotations

import csv
import glob
import hashlib
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MAIN = os.environ.get("A2_MAIN_ROOT", ROOT)                 # main checkout (results / cache are written there)
W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
R = os.path.join(MAIN, "results", "2026-10-27_A2_dn2")
CACHE = os.path.join(MAIN, "data", "cache", "a2")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def events(run_label):
    return [e for e in map(json.loads, open(os.path.join(W, "walk_events.jsonl"))) if str(e.get("run")) == run_label]


def select(run_label="2r"):
    log = json.load(open(os.path.join(MAIN, "data", "manifests", "dn2_select.log.json")))
    dl = {}
    for e in events(run_label):
        if e["ev"] == "downloaded":
            dl.setdefault(e["accession"], []).append(e)
    with open(os.path.join(R, "selection_log.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["walk_pos", "accession", "organism", "result", "reason_or_slot", "upid", "file", "size", "sha256",
                    "sha1", "pride_sha1_source", "sha1_status"])
        for i, x in enumerate(log["log"]):
            why = x.get("reason") or x.get("slot") or (",".join(x.get("trigger_members", [])))
            files = dl.get(x["accession"], []) if x["result"] in ("admitted", "marker screen exclude") else []
            if not files:
                w.writerow([i, x["accession"], x["organism"], x["result"], why, "", "", "", "", "", "", ""])
            for f in files:
                upid = next((e["upid"] for e in events(run_label) if e["ev"] == "marker_admit"
                             and e["accession"] == x["accession"]), "")
                w.writerow([i, x["accession"], x["organism"], x["result"], why, upid, f["filename"],
                            f["expected_size"], f["sha256"], f["sha1"], f["pride_sha1_source"], f["sha1_status"]])


def admit():
    log = json.load(open(os.path.join(MAIN, "data", "manifests", "dn2_select.log.json")))
    per, lines = {}, []
    for acc in log["admitted"]:
        d = os.path.join(W, "admission", acc)
        c = json.load(open(os.path.join(d, "counts.json")))
        per[acc] = {"upid": c["upid"], "runs": c["runs"], "counts": c["counts"],
                    "membership_sha256": c["membership_sha256"]}
        lines.append(f"{c['membership_sha256']}  membership:{acc}")
        for f in sorted(glob.glob(os.path.join(d, "*"))):
            lines.append(f"{sha256(f)}  admission/{acc}/{os.path.basename(f)}")
    tot = {k: sum(p["counts"].get(k, 0) for p in per.values())
           for k in ("n_spectra", "n_certified", "n_A", "n_B_candidates", "n_A_with_candidate", "tail_proxy",
                     "removed_b_f", "removed_b_hamming", "removed_d_family", "removed_c_pre")}
    swap = os.path.join(W, "swap_test.json")
    out = {"admitted": log["admitted"], "per_project": per, "totals": tot, "final_state": log["final_state"],
           "fs_donors_equal_n_A": tot["n_A"],
           "host_threshold": json.load(open(os.path.join(W, "host_threshold.json"))),
           "swap_test": json.load(open(swap)) if os.path.exists(swap) else None,
           "note": "label-free counts only (A2.9.4); B candidates are final only after (c-final) in step 5"}
    json.dump(out, open(os.path.join(R, "admission_counts.json"), "w"), indent=1)
    open(os.path.join(R, "membership.sha256"), "w").write("\n".join(lines) + "\n")
    print(json.dumps({"admitted": log["admitted"], "totals": tot, "final_state": log["final_state"]}, indent=1))


def persist(dirs):
    for d in dirs:
        src, dst = os.path.join(W, d), os.path.join(CACHE, d)
        lines = []
        for p in sorted(glob.glob(os.path.join(src, "**", "*"), recursive=True)):
            if os.path.isfile(p):
                rel = os.path.relpath(p, W)
                q = os.path.join(CACHE, rel)
                os.makedirs(os.path.dirname(q), exist_ok=True)
                shutil.copyfile(p, q + ".partial")
                os.replace(q + ".partial", q)
                h = sha256(q)
                if h != sha256(p):
                    raise SystemExit(f"copy mismatch {rel}")
                lines.append(f"{h}  {rel}")
        name = d.replace("/", "_")
        open(os.path.join(R, f"cache_{name}.sha256"), "w").write("\n".join(lines) + "\n")
        print(d, len(lines), "files ->", dst)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "select":
        select(*(sys.argv[2:3] or []))
    elif cmd == "admit":
        admit()
    elif cmd == "persist":
        persist(sys.argv[2:])
