#!/usr/bin/env python3
"""Select the project-disjoint TRAINING runs for the static-first scorer (DESIGN §4 L3, M3a).  Stdlib only.

Frame: data/manifests/t4_frame.tsv (eligible PRIDE DDA projects, crawled 2026-10-26).  Exclusions (hard):
  * every accession in t4_selection.tsv (sealed T4 sample), tier_s.tsv / dev_sealed*.tsv (D*/T*),
  * every accession in predictor_training_accessions.tsv (AlphaPeptDeep/DeepLC/MS2PIP/Prosit training),
  * yeast (held out entirely, memorization arm M / T5) and Arabidopsis (arm-C entrapment proteome),
  * enrichment != none, label != LFQ, keywords suggesting HLA / XL / glyco / SILAC / DIA / phospho / single-cell.
Instrument class OT-HCD: Q Exactive family and Exploris only (Fusion/Lumos/Eclipse may use IT-CID MS2).
Seeded (pre-registered project-fold seed 11, DESIGN §6.1): projects shuffled per species, raw files 0.3–1.3 GB chosen at random.
Output: data/manifests/train.tsv (id, accession, filename, url, expected_size, species, instrument, fold, role=train)
Fold map: project -> fold (0..2), assigned round-robin over the seeded project order (committed with the manifest).

    python3 scripts/m3/select_train.py [--human 9 --mouse 2 --ecoli 2 --runs-per-project 1]
"""
import argparse, csv, json, os, random, re, sys, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100&page=%d"
SPECIES = {"human": "Homo sapiens (human)", "mouse": "Mus musculus (mouse)", "ecoli": "Escherichia coli"}
INSTR = re.compile(r"^(Q Exactive( HF| HF-X| Plus)?|Orbitrap Exploris (480|240))$")
BAD_KW = re.compile(r"HLA|MHC|immunopept|cross-?link|XL-MS|glyco|SILAC|DIA\b|SWATH|data.independent|phospho|single.cell|"
                    r"TMT|iTRAQ|ubiquitin|acetyl|N-?termin|top-down|crosslink|metaproteom|ribosom|interactome|AP-MS|BioID|"
                    r"pull.?down|immunoprecip", re.I)
SEED = 11
BAD_FILE = re.compile(r"phos|TMT|blank|wash|QC|std|standard|DIA|ptm|enrich|IP", re.I)


def accessions(path, col):
    with open(path) as fh:
        return {r[col] for r in csv.DictReader(fh, delimiter="\t") if r.get(col)}


def pride_files(acc):
    out, page = [], 0
    while True:
        with urllib.request.urlopen(urllib.request.Request(PRIDE % (acc, page), headers={"User-Agent": "vulcan-m3"}),
                                    timeout=60) as r:
            got = json.load(r)
        out += got
        if len(got) < 100:
            return out
        page += 1


def ftp_url(f):
    for loc in f.get("publicFileLocations", []):
        if loc.get("value", "").startswith("ftp://"):
            return loc["value"].replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--human", type=int, default=9)
    ap.add_argument("--mouse", type=int, default=2)
    ap.add_argument("--ecoli", type=int, default=2)
    ap.add_argument("--runs-per-project", type=int, default=1)
    ap.add_argument("--min-bytes", type=float, default=3e8)
    ap.add_argument("--max-bytes", type=float, default=1.3e9)
    a = ap.parse_args()
    excl = set()
    excl |= accessions(os.path.join(MAN, "t4_selection.tsv"), "accession")
    for m in ("tier_s.tsv", "dev_sealed.tsv", "dev_sealed_extra.tsv"):
        excl |= accessions(os.path.join(MAN, m), "accession")
    excl |= accessions(os.path.join(MAN, "predictor_training_accessions.tsv"), "accession")
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    rng = random.Random(SEED)
    rows, log = [], {"excluded_accessions": len(excl)}
    for sp, n in (("human", a.human), ("mouse", a.mouse), ("ecoli", a.ecoli)):
        cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] == SPECIES[sp] and r["acquisition"] == "DDA"
                and r["label"] == "LFQ" and r["enrichment"] == "none" and r["vendor"] == "thermo"
                and INSTR.match(r["instruments"]) and r["accession"] not in excl and r["year"] >= "2017"
                and not BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])]
        cand.sort(key=lambda r: r["accession"])
        rng.shuffle(cand)
        log[sp] = {"candidates": len(cand), "tried": []}
        got = 0
        for r in cand:
            if got >= n:
                break
            try:
                files = pride_files(r["accession"])
            except Exception as e:  # noqa: BLE001
                log[sp]["tried"].append([r["accession"], f"api error {e}"]); continue
            raws = sorted([f for f in files if f["fileName"].lower().endswith(".raw") and not BAD_FILE.search(f["fileName"])
                           and a.min_bytes <= f["fileSizeBytes"] <= a.max_bytes and ftp_url(f)], key=lambda f: f["fileName"])
            if len(raws) < a.runs_per_project:
                log[sp]["tried"].append([r["accession"], f"{len(raws)} raw files in size range"]); continue
            pick = rng.sample(raws, a.runs_per_project)
            for f in pick:
                rows.append({"id": f"TR{len(rows):02d}", "accession": r["accession"], "filename": f["fileName"],
                             "url": ftp_url(f), "expected_size": str(f["fileSizeBytes"]), "species": sp,
                             "instrument": r["instruments"], "year": r["year"], "lab": r["lab"], "role": "train"})
            log[sp]["tried"].append([r["accession"], "selected"])
            got += 1
    # fold map: round-robin over the seeded project order (project-disjoint folds)
    projs = []
    for r in rows:
        if r["accession"] not in projs:
            projs.append(r["accession"])
    fold = {p: i % 3 for i, p in enumerate(projs)}
    for r in rows:
        r["fold"] = str(fold[r["accession"]])
    cols = ["id", "accession", "filename", "url", "expected_size", "species", "instrument", "year", "lab", "fold", "role"]
    with open(os.path.join(MAN, "train.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c].replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(MAN, "train_folds.tsv"), "w") as fh:
        fh.write("accession\tfold\tspecies\n")
        for p in projs:
            fh.write(f"{p}\t{fold[p]}\t{next(r['species'] for r in rows if r['accession'] == p)}\n")
    json.dump(log, open(os.path.join(MAN, "train_select.log.json"), "w"), indent=1)
    print(json.dumps(log, indent=1))
    print(f"{len(rows)} runs from {len(projs)} projects, {sum(int(r['expected_size']) for r in rows)/1e9:.1f} GB")


if __name__ == "__main__":
    main()
