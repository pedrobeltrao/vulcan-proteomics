#!/usr/bin/env python3
"""C9 pilot: METADATA-ONLY count of the human OT-HCD DDA label-free pool a removed-peptide null set could draw from.

Frame rule = the E9/D8 frame (scripts/e9/select_e9.py: eligible, Homo sapiens only, DDA, LFQ, no enrichment, Thermo
single-instrument QE / QE Plus / HF / HF-X / Exploris 240/480, BAD_KW title/keyword screen).  Exclusions (stricter
than E9): every accession in data/manifests/*.tsv except t4_frame.tsv (dev, sealed T1-T5 sets, T4 selection, tier S,
train, train_itcid, train_m7, predictor training, DN2, DNIT, D8/E9 incl. every project E9 tried), every T4 draw found
under results/**/t4_selector/**/t4_selection.tsv (T4-reachable), and the C9 pilot's own dev projects.
Clusters: the crawl 'lab' field, normalised order-insensitively (as reviews/checks/a4_r3_walk.py).
No file listing, no download, no spectrum.  Writes results/2026-10-28_C9_pilot/pool_count.json.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "e9"))
import select_e9 as S  # noqa: E402

MAN = os.path.join(ROOT, "data", "manifests")
OUT = os.path.join(ROOT, "results", "2026-10-28_C9_pilot", "pool_count.json")


def norm_lab(s):
    toks = sorted(t for t in re.findall(r"[a-z0-9]+", (s or "").lower()) if len(t) > 1)
    return " ".join(toks)


def main():
    excl = {}
    for f in sorted(glob.glob(os.path.join(MAN, "*.tsv"))):
        if os.path.basename(f) == "t4_frame.tsv":
            continue
        try:
            rows = list(csv.DictReader(open(f), delimiter="\t"))
        except Exception:  # noqa: BLE001
            continue
        for r in rows:
            for k in ("accession", "Accession", "project"):
                a = r.get(k)
                if a and a.startswith("PXD"):
                    excl.setdefault(a, set()).add(os.path.basename(f))
    for f in glob.glob(os.path.join(ROOT, "results", "**", "t4_selector", "**", "t4_selection.tsv"), recursive=True):
        for r in csv.DictReader(open(f), delimiter="\t"):
            excl.setdefault(r["accession"], set()).add("T4-reachable:" + os.path.relpath(f, ROOT))
    e9log = json.load(open(os.path.join(MAN, "e9_d8.log.json")))
    tried = e9log.get("tried", [])
    n_usable = []                               # E9 listing: usable RAW files per tried human project (label-free)
    for t in tried:
        a, st = (t[0], t[1]) if isinstance(t, (list, tuple)) else (t, "")
        excl.setdefault(a, set()).add("e9_tried")
        m = re.match(r"selected \d+ of (\d+)", st)
        n_usable.append(int(m.group(1)) if m else 0)
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] == "Homo sapiens (human)"
            and r["acquisition"] == "DDA" and r["label"] == "LFQ" and r["enrichment"] == "none"
            and r["vendor"] == "thermo" and S.INSTR.match(r["instruments"])
            and not S.BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])]
    n_frame = len(cand)
    reasons = {}
    for r in cand:
        if r["accession"] in excl:
            for x in excl[r["accession"]]:
                reasons[x.split(":")[0]] = reasons.get(x.split(":")[0], 0) + 1
    pool = [r for r in cand if r["accession"] not in excl]
    labs = {}
    for r in pool:
        labs.setdefault(norm_lab(r["lab"]) or ("?" + r["accession"]), []).append(r["accession"])
    sizes = sorted((len(v) for v in labs.values()), reverse=True)
    by_instr = {}
    for r in pool:
        by_instr[r["instruments"]] = by_instr.get(r["instruments"], 0) + 1
    out = {"frame_rule": "E9/D8 frame (scripts/e9/select_e9.py)", "frame_candidates": n_frame,
           "excluded": n_frame - len(pool), "exclusion_reasons_hits": reasons, "pool_projects": len(pool),
           "pool_labs": len(labs), "lab_size_top10": sizes[:10],
           "labs_with_1_project": sum(1 for s in sizes if s == 1), "pool_by_instrument": by_instr,
           "e9_tried_listed": len(tried),
           "e9_usable_raw_per_project": {"n": len(n_usable), "zero": sum(x == 0 for x in n_usable),
                                         "quantiles_p10_25_50_75_90": [sorted(n_usable)[int(q * (len(n_usable) - 1))]
                                                                       for q in (0.1, 0.25, 0.5, 0.75, 0.9)],
                                         "values": n_usable},
           "note": "metadata only; usable RAW files per project unknown before listing (E9: 142 of 166 tried projects "
                   "had >= 1 usable RAW in 0.1-2 GB, mean 1.08 GB)"}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
