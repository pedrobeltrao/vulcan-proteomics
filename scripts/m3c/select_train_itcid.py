#!/usr/bin/env python3
"""Select the project-disjoint TRAINING runs for the IT-CID static-first scorer (tree 3c; DESIGN §4 L3, same recipe
as scripts/m3/select_train.py for OT-HCD).  Stdlib only.  Written and committed BEFORE any download.

Frame: data/manifests/t4_frame.tsv (eligible PRIDE DDA projects, crawled 2026-10-26).
Hard exclusions:
  * the union of the `accession` column of EVERY data/manifests/*.tsv except t4_frame.tsv (dev, sealed, tier S, T4
    selection, M3 OT-HCD training (train.tsv), predictor-training accessions (AlphaPeptDeep/DeepLC/MS2PIP/Prosit,
    incl. the Prosit 2020 CID/TMT training projects), DN (g2null), D8/E9 stream, D7b, D9-D12, D10cf, ...);
  * species other than human / mouse (so: no yeast (arm M / T5), no Arabidopsis (arm C), and no free-living microbe,
    i.e. nothing from the A2 DN2 pool; the A2 selector in scripts/a2 is NOT run);
  * enrichment != none; keywords suggesting HLA / XL / glyco / SILAC / DIA / phospho / single-cell / AP-MS / ...
Strata (instrument class IT-CID = ion-trap CID MS2 with Orbitrap MS1):
  TMT  label == TMT, instrument in {Orbitrap Fusion, Fusion Lumos, Eclipse}   (SPS-MS3 expected)      -> 6 projects
  LF   label == LFQ, instrument in {LTQ Orbitrap (XL/Velos/Elite), Orbitrap Fusion, Fusion Lumos}     -> 7 projects
Order: per stratum, candidates sorted by accession and shuffled with random.Random(11) (project-fold seed, §6.1);
one .raw per project, 0.15-1.5 GB, file name not matching BAD_FILE, drawn with the same rng.
Metadata cannot see the MS2 analyzer (LESSONS), so the candidate list is longer than the target; the fetcher
(fetch_train_itcid.py) walks it IN ORDER and accepts a run only if, after ingest (label-free, pre-specified here):
  * >= 95% of MS2 spectra are ion-trap (analyzer IT) + CID;
  * TMT stratum: MS3 spectra >= 50% of MS2 count (SPS-MS3), and the label is TMT 6/10/11 (229.163 Da): MS3 reporter
    peaks at 126.1277..131.1445 present and TMTpro-only reporters (132.14/133.14/134.15/135.15 +-0.003) absent
    in >= 90% of MS3 with any reporter;
  * LF stratum: no MS3 spectra (an MS3 would indicate an isobaric-label method).
Rejected runs are logged with the reason and replaced by the next candidate of the same stratum.
Fold map: 3 project-disjoint folds, assigned round-robin within each stratum in acceptance order (so every fold
has TMT and LF runs); written by the fetcher with the final manifest.

    python3 scripts/m3c/select_train_itcid.py          # -> data/manifests/train_itcid_candidates.tsv + log
"""
import argparse, csv, glob, json, os, random, re, sys, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100&page=%d"
SPECIES = {"Homo sapiens (human)": "human", "Mus musculus (mouse)": "mouse"}
INSTR = {"TMT": re.compile(r"^Orbitrap (Fusion|Fusion Lumos|Eclipse)$"),
         "LF": re.compile(r"^(LTQ Orbitrap( XL| Velos| Elite)?|Orbitrap Fusion|Orbitrap Fusion Lumos)$")}
LABEL = {"TMT": "TMT", "LF": "LFQ"}
TARGET = {"TMT": 6, "LF": 7}
N_CAND = {"TMT": 30, "LF": 30}
BAD_KW = re.compile(r"HLA|MHC|immunopept|cross-?link|XL-MS|glyco|SILAC|DIA\b|SWATH|data.independent|phospho|single.cell|"
                    r"iTRAQ|ubiquitin|acetyl|N-?termin|top-down|crosslink|metaproteom|ribosom|interactome|AP-MS|BioID|"
                    r"pull.?down|immunoprecip|TMTpro|16.?plex|18.?plex|ETD|EThcD|UVPD", re.I)
BAD_FILE = {"TMT": re.compile(r"phos|blank|wash|QC|std|standard|DIA|ptm|enrich|IP|pro16|tmtpro", re.I),
            "LF": re.compile(r"phos|TMT|blank|wash|QC|std|standard|DIA|ptm|enrich|IP", re.I)}
SEED = 11


def excluded_accessions():
    excl = set()
    for f in sorted(glob.glob(os.path.join(MAN, "*.tsv"))):
        if os.path.basename(f) in ("t4_frame.tsv", "train_itcid_candidates.tsv", "train_itcid.tsv",
                                   "train_itcid_folds.tsv"):
            continue
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if rd.fieldnames and "accession" in rd.fieldnames:
                excl |= {r["accession"] for r in rd if r.get("accession")}
    return excl


def pride_files(acc):
    out, page = [], 0
    while True:
        with urllib.request.urlopen(urllib.request.Request(PRIDE % (acc, page), headers={"User-Agent": "vulcan-3c"}),
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
    ap.add_argument("--min-bytes", type=float, default=1.5e8)
    ap.add_argument("--max-bytes", type=float, default=1.5e9)
    a = ap.parse_args()
    excl = excluded_accessions()
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    rng = random.Random(SEED)
    rows, log = [], {"excluded_accessions": len(excl), "seed": SEED, "target": TARGET}
    for stratum in ("TMT", "LF"):
        cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] in SPECIES and r["acquisition"] == "DDA"
                and r["label"] == LABEL[stratum] and r["enrichment"] == "none" and r["vendor"] == "thermo"
                and INSTR[stratum].match(r["instruments"]) and r["accession"] not in excl
                and not BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])]
        cand.sort(key=lambda r: r["accession"])
        rng.shuffle(cand)
        log[stratum] = {"candidates": len(cand), "tried": []}
        got = 0
        for r in cand:
            if got >= N_CAND[stratum]:
                break
            try:
                files = pride_files(r["accession"])
            except Exception as e:  # noqa: BLE001
                log[stratum]["tried"].append([r["accession"], f"api error {e}"]); continue
            raws = sorted([f for f in files if f["fileName"].lower().endswith(".raw")
                           and not BAD_FILE[stratum].search(f["fileName"])
                           and a.min_bytes <= f["fileSizeBytes"] <= a.max_bytes and ftp_url(f)],
                          key=lambda f: f["fileName"])
            if not raws:
                log[stratum]["tried"].append([r["accession"], "0 raw files in size range"]); continue
            f = rng.choice(raws)
            rows.append({"cand": f"{stratum}{got:02d}", "stratum": stratum, "accession": r["accession"],
                         "filename": f["fileName"], "url": ftp_url(f), "expected_size": str(f["fileSizeBytes"]),
                         "species": SPECIES[r["organisms"]], "instrument": r["instruments"], "year": r["year"],
                         "lab": r["lab"]})
            log[stratum]["tried"].append([r["accession"], "candidate"])
            got += 1
            print(stratum, got, r["accession"], r["instruments"], f["fileName"], flush=True)
    cols = ["cand", "stratum", "accession", "filename", "url", "expected_size", "species", "instrument", "year", "lab"]
    with open(os.path.join(MAN, "train_itcid_candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    json.dump(log, open(os.path.join(MAN, "train_itcid_select.log.json"), "w"), indent=1)
    print(f"{len(rows)} candidates")


if __name__ == "__main__":
    sys.exit(main())


# ------------------------------------------------------------------------------------------------ DN2 guard
def dn2_guard(accessions, organisms=()):
    """PI instruction (2026-10-27, before any null download): never use an A2 DN2 candidate for any purpose.
    Mechanically refuse (raise) if any accession appears in data/manifests/dn2*.tsv / dn2_select.log.json or in
    results/2026-10-27_A2_dn2/{selection_log.tsv, walk_events.jsonl} (every project the walk has seen), or if an
    organism could be reached by the DN2 ALLOW-list walk (genus on config/a2.yaml `allow`, not excluded there).
    Read-only use of the A2 files; scripts/a2 is not run.  Returns the number of DN2 accessions checked against."""
    seen = set()
    files = glob.glob(os.path.join(MAN, "dn2*")) + [os.path.join(ROOT, "results", "2026-10-27_A2_dn2", f)
                                                    for f in ("selection_log.tsv", "walk_events.jsonl")]
    for f in files:
        if os.path.exists(f):
            seen |= set(re.findall(r"PXD\d{6}", open(f).read()))
    bad = sorted(set(accessions) & seen)
    if bad:
        raise SystemExit(f"DN2 guard: {bad} are DN2 candidates; refusing")
    allow = re.findall(r"\{genus: (\w+)", open(os.path.join(ROOT, "config", "a2.yaml")).read())
    prev = ["pyrococcus furiosus", "saccharomyces cerevisiae", "trypanosoma brucei", "staphylococcus aureus",
            "synechococcus", "bacillus cereus", "homo sapiens", "mus musculus", "escherichia coli", "arabidopsis"]
    for o in organisms:
        g = o.split()[0] if o.split() else ""
        if g in allow and not any(o.lower().startswith(p) for p in prev):
            raise SystemExit(f"DN2 guard: organism {o} is reachable by the DN2 ALLOW walk; refusing")
    return len(seen)
