#!/usr/bin/env python3
"""E9 run stream (DESIGN §7 E9, tree 3d): a seeded, nested sample of human DDA Orbitrap-HCD runs from the PRIDE crawl
frame (data/manifests/t4_frame.tsv), recorded as dev-role set D8.  Stdlib only; metadata-only (no RAW is read here).

Hard exclusions (project level):
  * every T4 *selection* accession (data/manifests/t4_selection.tsv; sealed T4 sample),
  * every dev / sealed accession (tier_s.tsv, dev_sealed.tsv, dev_sealed_extra.tsv),
  * the 13 M3 static-scorer training projects (train.tsv, DATA.md role=train),
  * every predictor-training accession (predictor_training_accessions.tsv, any confidence/role),
  * keywords/titles suggesting HLA / XL / glyco / SILAC / DIA / phospho / TMT / single-cell / AP-MS etc. (as M3 train).
Frame filter: eligible, organism human, DDA, label LFQ, enrichment none, Thermo, instrument class OT-HCD
(Q Exactive family / Exploris 240/480, single-instrument projects) = the instrument class of the default scorer (D7).
Files: RAW category, .raw, 0.1-2.0 GB, names not suggesting QC/blank/wash/standard/DIA/phospho/TMT; <= 8 per project.

Sampling (seed 11, random.Random(11)): projects sorted by accession and shuffled; files of each project sorted by name
and sampled with random.Random(11 + accession digits).  Global run order = round-robin over the seeded project order
(round 1 = first run of every project, round 2 = second run, ...), so every prefix N spans as many projects as possible
and N = 1 / 10 / 100 / 1000 are NESTED prefixes of one list.  Projects are drawn until the list holds --n-runs runs.
Projects whose /files call fails or that have < --min-files usable files are skipped (logged).

Output: data/manifests/e9_d8.tsv (order, id, accession, filename, url, expected_size, pride_sha1, instrument, year,
lab, round, role=dev, set_id=D8) and data/manifests/e9_d8.log.json.  Downloaded sha256s are appended by the E9 streamer
(data/manifests/e9_d8.sha256).

    python3 scripts/e9/select_e9.py [--n-runs 1000] [--per-project 8]
"""
import argparse
import csv
import json
import os
import random
import re
import sys
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100&page=%d"
SEED = 11
INSTR = re.compile(r"^(Q Exactive( HF| HF-X| Plus)?|Orbitrap Exploris (480|240))$")
BAD_KW = re.compile(r"HLA|MHC|immunopept|cross-?link|XL-MS|glyco|SILAC|DIA\b|SWATH|data.independent|phospho|single.cell|"
                    r"TMT|iTRAQ|ubiquitin|acetyl|N-?termin|top-down|crosslink|metaproteom|interactome|AP-MS|BioID|"
                    r"pull.?down|immunoprecip|labell?ing|dimethyl|PRM|SRM|targeted", re.I)
BAD_FILE = re.compile(r"phos|TMT|blank|wash|QC|std|standard|DIA|SWATH|ptm|enrich|\bIP\b|_IP_|hela_?digest|BSA|calib", re.I)


def accessions(path, col="accession"):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {r[col] for r in csv.DictReader(fh, delimiter="\t") if r.get(col)}


def pride_files(acc):
    out, page = [], 0
    while True:
        for attempt in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(PRIDE % (acc, page), headers={"User-Agent": "vulcan-e9"}),
                                            timeout=60) as r:
                    got = json.load(r)
                break
            except Exception:  # noqa: BLE001
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt + 1))
        out += got
        if len(got) < 100:
            return out
        page += 1


def ftp_url(f):
    for loc in f.get("publicFileLocations", []):
        if loc.get("value", "").startswith("ftp://"):
            return loc["value"].replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
    return None


def exclusions():
    ex = {}
    for name in ("t4_selection.tsv", "tier_s.tsv", "dev_sealed.tsv", "dev_sealed_extra.tsv", "train.tsv",
                 "predictor_training_accessions.tsv"):
        for a in accessions(os.path.join(MAN, name)):
            ex.setdefault(a, []).append(name)
    return ex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-runs", type=int, default=1000)
    ap.add_argument("--per-project", type=int, default=8)
    ap.add_argument("--min-files", type=int, default=1)
    ap.add_argument("--min-bytes", type=float, default=1e8)
    ap.add_argument("--max-bytes", type=float, default=2e9)
    ap.add_argument("--out", default=os.path.join(MAN, "e9_d8.tsv"))
    a = ap.parse_args()
    excl = exclusions()
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] == "Homo sapiens (human)"
            and r["acquisition"] == "DDA" and r["label"] == "LFQ" and r["enrichment"] == "none"
            and r["vendor"] == "thermo" and INSTR.match(r["instruments"])
            and not BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])]
    n_frame = len(cand)
    n_excl = sum(r["accession"] in excl for r in cand)
    cand = [r for r in cand if r["accession"] not in excl]
    cand.sort(key=lambda r: r["accession"])
    random.Random(SEED).shuffle(cand)
    log = {"seed": SEED, "frame_candidates": n_frame, "excluded_by_manifest": n_excl,
           "exclusion_manifest_sizes": {k: len(accessions(os.path.join(MAN, k))) for k in
                                        ("t4_selection.tsv", "tier_s.tsv", "dev_sealed.tsv", "dev_sealed_extra.tsv",
                                         "train.tsv", "predictor_training_accessions.tsv")},
           "eligible_after_exclusion": len(cand), "tried": []}
    projects, total = [], 0
    for r in cand:
        if total >= a.n_runs:
            break
        acc = r["accession"]
        try:
            files = pride_files(acc)
        except Exception as e:  # noqa: BLE001
            log["tried"].append([acc, f"api error {e}"])
            continue
        raws = sorted([f for f in files if f.get("fileCategory", {}).get("value") == "RAW"
                       and f["fileName"].lower().endswith(".raw") and not BAD_FILE.search(f["fileName"])
                       and a.min_bytes <= f["fileSizeBytes"] <= a.max_bytes and ftp_url(f)], key=lambda f: f["fileName"])
        if len(raws) < a.min_files:
            log["tried"].append([acc, f"{len(raws)} usable raw files"])
            continue
        rng = random.Random(SEED + int(re.sub(r"\D", "", acc)))
        pick = rng.sample(raws, min(a.per_project, len(raws)))
        projects.append((r, pick))
        total += len(pick)
        log["tried"].append([acc, f"selected {len(pick)} of {len(raws)}"])
        print(acc, len(pick), total, flush=True)
    # round-robin global order over the seeded project order
    rows, rnd = [], 0
    while len(rows) < a.n_runs and any(rnd < len(p) for _, p in projects):
        for r, pick in projects:
            if rnd < len(pick) and len(rows) < a.n_runs:
                f = pick[rnd]
                rows.append({"order": str(len(rows)), "id": f"E9_{len(rows):04d}", "accession": r["accession"],
                             "filename": f["fileName"], "url": ftp_url(f), "expected_size": str(f["fileSizeBytes"]),
                             "pride_sha1": f.get("checksum", ""), "instrument": r["instruments"], "year": r["year"],
                             "lab": r["lab"].replace("\t", " "), "round": str(rnd), "role": "dev", "set_id": "D8"})
        rnd += 1
    for row in rows:  # paranoia: hard exclusion re-check
        assert row["accession"] not in excl, row
    cols = ["order", "id", "accession", "filename", "url", "expected_size", "pride_sha1", "instrument", "year", "lab",
            "round", "role", "set_id"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in cols) + "\n")
    log["n_runs"] = len(rows)
    log["n_projects"] = len({r["accession"] for r in rows})
    log["bytes"] = sum(int(r["expected_size"]) for r in rows)
    for n in (1, 10, 100, 1000):
        sub = rows[:n]
        log[f"N{n}"] = {"projects": len({r["accession"] for r in sub}), "bytes": sum(int(r["expected_size"]) for r in sub)}
    json.dump(log, open(a.out.replace(".tsv", ".log.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in log.items() if k != "tried"}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
