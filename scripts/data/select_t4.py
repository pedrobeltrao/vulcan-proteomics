#!/usr/bin/env python3
"""T4 selector (DESIGN.md 6.6), metadata-only dry run.  No raw file is ever downloaded.

Usage:  python3 scripts/data/select_t4.py [--stage crawl|frame|select|all] [--refresh] [--out-dir DIR] [--check [TSV]]
        (--exclusion-mode draw|frame, default draw: see item 7; --out-dir DIR writes t4_frame.tsv, t4_selection.tsv, summary.md, all_projects_classified.tsv to DIR instead)
Outputs (project-root relative):
  data/manifests/t4_frame.tsv       full metadata-eligible frame + classification fields
  data/manifests/t4_selection.tsv   selected projects x files with sizes
  results/2026-10-26_T4_dryrun/summary.md
Cache: data/cache/pride_meta/ (projects_p<N>.json pages, files_<ACC>_p<N>.json, uniprot_reference_proteomes.tsv).
Re-running without --refresh uses only the cache, so the frame is reproducible from a frozen crawl.

SOURCES (PRIDE Archive API v3, https://www.ebi.ac.uk/pride/ws/archive/v3/)
  /projects?pageSize=100&page=N   (page size is capped at 100; ~418 pages; sequential, 0.5 s pause, retries)
  /projects/<ACC>/files           (ONLY for selection candidates, to confirm raw files and read sizes)
  UniProt /proteomes/stream?query=reference:true  (organism -> reference proteome check, by species name in lineage)

FRAME (project-level, all from metadata; every rule is in classify_project())
  1. vendor   : instruments all Thermo or all Bruker (mixed-vendor and other vendors, incl. SCIEX TripleTOF, Waters,
                Agilent, dropped). Thermo ion-trap-only (LTQ/LCQ/Velos-only without Orbitrap) dropped: no HR MS2.
  2. acquisition: DIA/SWATH/diaPASEF/PRM/SRM evidence -> dropped (see DIA RULES). Otherwise DDA.
  3. label    : TMT (TMT/iTRAQ/tandem mass tag) or LFQ (label-free/LFQ/iBAQ/spectral-count in structured fields or text),
                SILAC/dimethyl/metabolic/15N -> dropped; TMT+LFQ both -> TMT.
  4. enrichment: phospho = phosphoproteomics/phosphopeptide + an enrichment term (TiO2, IMAC, Fe-NTA, ...) in
                title/keywords/tags/protocols; else none. Phospho projects are in the frame even with label 'unspecified'.
     frame membership: label in {LFQ,TMT} or enrichment == phospho.
  5. XL (crosslink, DSS, DSSO, BS3, XL-MS ...), top-down, MALDI, metaproteomics, multi-organism -> dropped.
  6. organism : single organism whose species name occurs in a UniProt *reference* proteome lineage.
  7. exclusions (mechanical, review M7; each adds its own exclude_reasons tag):
                dev_sealed:<set>       hard-coded DEV_SEALED accessions below;
                train                  data/manifests/train.tsv (required) + DATA.md table rows with role `train`;
                dev_manifest:<file>    every accession (all roles) of dev_sealed, tier_s, g2null (required) and
                                       dev_sealed_extra, d9, d4_ref (optional: warn if missing);
                e9_d8_stream           data/manifests/e9_d8.tsv, the D8/E9 stream (required);
                predictor_training     data/manifests/predictor_training_accessions.tsv (required), confidence in
                                       {high, medium} and role train/train_or_test/unknown/derived; eval/holdout/low rows only flagged.
                --exclusion-mode draw (default) allocates and shuffles on the frame without the train/dev_manifest/
                e9_d8_stream exclusions and skips those projects while drawing, so the seed-11 selection is unchanged
                when none of them would have been drawn; --exclusion-mode frame removes them from the frame first.
                assert_disjoint(rows) checks a selection against all of these; `--check [TSV]` runs it from the CLI and
                the selection stage runs it before writing t4_selection.tsv.
  8. raw files present: checked at selection time through /files (RAW category, right extension for the vendor: Thermo
                .raw; Bruker .d / .d.zip / .d.tar*). Frame size is therefore 'metadata-eligible'; candidates that fail
                the check are replaced by the next seeded candidate in the same cell and the failure rate is reported.
  Year bins by submissionDate: <=2014, 2015-2017, 2018-2020, 2021-2023, >=2024.
  Instrument family: orbitrap (Orbitrap/Q Exactive/Exploris/Fusion/Lumos/Eclipse/Velos/Elite/Exactive),
                astral (Orbitrap Astral), tims (timsTOF), bruker_qtof (maXis, impact, compact, others Bruker).
  Projects with both Orbitrap and Astral -> astral. Cell = year_bin x family x label x enrichment (label 'UNK' for
  phospho projects with no TMT/LFQ declaration).

DIA RULES (error-prone; see summary.md)
  strong DIA -> drop: PRIDE experimentType or keyword/title matches DIA|SWATH|diaPASEF|data-independent|DIA-NN|
                Spectronaut|EncyclopeDIA|OpenSWATH|Skyline-targeted (PRM/SRM/MRM) in structured fields or software.
  weak DIA (protocol text only): drop unless the text also mentions data-dependent/DDA/TopN (then 'mixed', dropped
                and counted as ambiguity); no DIA and no DDA evidence -> 'dda_assumed' (kept, lower confidence).
  Files whose names contain DIA/SWATH/diaPASEF are dropped at file level.

SAMPLING (seed 11, random.Random(11) everywhere; target N_PROJECTS = 40, allowed 30-50)
  Allocation = systematic PPS over cells with weight sqrt(n_cell) (square-root damping so large cells do not
  dominate but every large cell is certain to be hit): cells sorted lexicographically by key; interval =
  sum(w)/N; hits at start+k*interval with start ~ U(0, interval) from the seeded RNG; a cell with weight > interval
  gets several hits, each capped at its eligible count. Within a cell, candidates are shuffled (seeded) and taken
  in order; a candidate failing the raw check (no raw file of the right type / <3 usable files) is skipped and the
  next one taken. One project is never taken twice; each submitter lab (first labPI name) at most 2 projects.
FILES: usable = RAW category, vendor extension, >= 20 MB, no DIA name; <= 8 chosen per project with
  random.Random(11 + int(acc digits)). Global caps: <= 300 runs (trim from projects with most files, one file at a
  time, seeded) and <= 1 TB (drop the largest remaining file of the project with the largest selected bytes while
  keeping >= 1 file per project).
INGEST ESTIMATE: spectra per GB assumed 40k (Thermo) / 25k (Bruker .d.zip), a rough planning figure, divided by 2k and
  5k spectra/s.
"""
import argparse
import collections
import csv
import gzip
import json
import math
import os
import random
import re
import sys
import time
import urllib.request
import urllib.error

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CACHE = os.path.join(ROOT, "data", "cache", "pride_meta")
MAN = os.path.join(ROOT, "data", "manifests")
RES = os.path.join(ROOT, "results", "2026-10-26_T4_dryrun")
API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
SEED = 11
N_PROJECTS = 40
MAX_FILES_PER_PROJECT = 8
MAX_RUNS = 300
MAX_BYTES = 10 ** 12
MIN_FILE_BYTES = 20 * 10 ** 6
MAX_PER_LAB = 2
PAUSE = 0.5
CRAWL_THREADS = 4

# --- exclusions (DESIGN.md 6.6 + DATA.md) ------------------------------------------------------------------
DEV_SEALED = {
    "PXD001468": "D0", "PXD001077": "D1", "PXD000138": "D2/T2", "PXD009449": "D2/T2", "PXD007058": "D2",
    "PXD014337": "D3", "PXD029252": "D3b/T3", "PXD016766": "D4", "PXD059608": "D6",
    "PXD049412": "T7", "PXD015446": "T8",
    # ProteoBench: DDA-QExactive module (Van Puyvelde et al., PXD028735), DIA-Astral (PXD046444),
    # DIA-plasma PYE (jPOST JPST003358, no PXD), diaPASEF module (PXD026600): excluded by accession.
    "PXD028735": "T1/ProteoBench-DDA-QExactive", "PXD046444": "T6/ProteoBench-DIA-Astral",
    "PXD026600": "ProteoBench-diaPASEF",
    # ProteomeTools (Prosit training; also the PXD009449 relatives) - excluded irrespective of Task B
    "PXD004732": "ProteomeTools", "PXD010595": "ProteomeTools", "PXD021013": "ProteomeTools",
}
PRED_TSV = os.path.join(MAN, "predictor_training_accessions.tsv")
DATA_MD = os.path.join(ROOT, "DATA.md")
# Manifest-driven exclusions (review M7). (file, reason, required). Every accession in the file (all roles: dev, sealed,
# train) is excluded from the frame with the given reason. A missing required manifest is a hard error: the sealed T4 frame
# must be reproducible from tracked inputs. Optional manifests only warn.
TRAIN_TSV = "train.tsv"
E9_D8_TSV = "e9_d8.tsv"
DEV_MANIFESTS = [("dev_sealed.tsv", True), ("dev_sealed_extra.tsv", False), ("d9.tsv", False), ("d4_ref.tsv", False),
                 ("tier_s.tsv", True), ("g2null.tsv", True),
                 # Amendment A2.9.1: the fresh null set DN2 (written by scripts/a2/select_dn2.py after prereg-v0.4.2).
                 # Optional until it exists; DN2 is drawn from outside the T4 selection, so draw mode is unchanged.
                 ("dn2.tsv", False)]

THERMO_RE = re.compile(r"orbitrap|q exactive|exactive|exploris|fusion|lumos|eclipse|astral|ltq|lcq|velos|elite|"
                       r"thermo|tsq|ion trap|tribrid|ltq-ft|ltq ft", re.I)
BRUKER_RE = re.compile(r"timstof|tims tof|maxis|impact|compact|bruker|amazon|hct|ultraflex|autoflex|solarix|"
                       r"apex|micrOTOF|esquire", re.I)
OTHER_RE = re.compile(r"triple\s?tof|qtrap|sciex|waters|synapt|xevo|q-tof|agilent|6\d\d\d|ab sciex|shimadzu|jeol|"
                      r"maldi|tof/tof|ltq-?orbitrap", re.I)
# 'ltq orbitrap' is Thermo; handled before OTHER (see vendor_of)
DIA_STRUCT = re.compile(r"\bDIA\b|swath|dia-?pasef|data[- ]independent|dia-?nn|spectronaut|encyclopedia|openswath|"
                        r"\bPRM\b|\bSRM\b|\bMRM\b|parallel reaction|selected reaction|multiple reaction|"
                        r"\bdia-umpire|msplit|\bBoxCar\b|\bnDIA\b|wide-window", re.I)
DIA_TEXT = re.compile(r"\bDIA\b|swath|dia-?pasef|data[- ]independent", re.I)
DDA_TEXT = re.compile(r"data[- ]dependent|\bDDA\b|\bTop[- ]?\d+\b|\btop[- ]?n\b|shotgun|\bDDA-PASEF\b|"
                      r"most intense (precursor|ion)|\bTopN\b", re.I)
TMT_RE = re.compile(r"\bTMT\w*|tandem mass tag|\biTRAQ\b|\bTMTpro\b", re.I)
LFQ_RE = re.compile(r"label[- ]?free|\bLFQ\b|\biBAQ\b|spectral count|\bemPAI\b|\bNSAF\b|\bTop3\b|MaxLFQ|"
                    r"unlabell?ed", re.I)
SILAC_RE = re.compile(r"\bSILAC\b|dimethyl|metabolic labell?ing|\b15N\b|\bSILAM\b|\bpSILAC\b|\bNeuCode\b|"
                      r"\bPulsed SILAC|18O", re.I)
XL_RE = re.compile(r"cross-?link|\bXL-?MS\b|\bDSSO?\b|\bBS3\b|\bDSBU\b|\bDSG\b|\bEDC\b|\bXLMS\b|\bCX-?MS\b|"
                   r"disuccinimid", re.I)
MISC_DROP = re.compile(r"top[- ]down|\bMALDI\b|metaproteom|\bimaging\b|intact protein|native ms|"
                       r"\bHDX\b|hydrogen[- ]deuterium", re.I)
PHOS_A = re.compile(r"phosphoproteom|phosphopeptide|phospho-?enrich|phosphorylation", re.I)
PHOS_B = re.compile(r"TiO2|titanium|IMAC|Fe-?NTA|Ti-?IMAC|Zr-?IMAC|zirconium|enrich|\bMOAC\b|Ga-?IMAC|Ni-?NTA", re.I)


def _http(url, tries=5):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "vulcan-t4-dryrun/0.1"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            last = e
            if e.code == 404:
                raise
        except Exception as e:  # noqa
            last = e
        time.sleep(2 * (i + 1))
    raise RuntimeError("failed %s: %s" % (url, last))


def cached_json(name, url, refresh=False, allow_network=True):
    p = os.path.join(CACHE, name)
    if os.path.exists(p) and not refresh:
        with gzip.open(p, "rt") as f:
            return json.load(f)
    if not allow_network:
        return None
    data = _http(url)
    time.sleep(PAUSE)
    os.makedirs(CACHE, exist_ok=True)
    with gzip.open(p, "wt") as f:
        f.write(data.decode())
    return json.loads(data)


# --- stage 1: crawl -------------------------------------------------------------------------------------------
def crawl(refresh):
    os.makedirs(CACHE, exist_ok=True)
    req = urllib.request.Request(API + "/projects?pageSize=1&page=0", headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        total = int(r.headers.get("total_records"))
    pages = math.ceil(total / 100)
    print("PRIDE projects:", total, "pages:", pages, flush=True)
    from concurrent.futures import ThreadPoolExecutor
    projs = {}

    def get(pg):
        return cached_json("projects_p%04d.json.gz" % pg, "%s/projects?pageSize=100&page=%d" % (API, pg), refresh)
    with ThreadPoolExecutor(CRAWL_THREADS) as ex:  # 4 concurrent requests, 0.5 s pause each: polite
        for pg, d in enumerate(ex.map(get, range(pages))):
            for p in d:
                projs[p["accession"]] = p
            if pg % 20 == 0:
                print("page", pg, len(projs), flush=True)
    with open(os.path.join(CACHE, "crawl_info.json"), "w") as f:
        json.dump({"total_records_header": total, "pages": pages, "unique_accessions": len(projs),
                   "crawl_date": time.strftime("%Y-%m-%d")}, f)
    return projs, total


def load_cached_projects():
    projs = {}
    pg = 0
    while True:
        p = os.path.join(CACHE, "projects_p%04d.json.gz" % pg)
        if not os.path.exists(p):
            break
        with gzip.open(p, "rt") as f:
            for x in json.load(f):
                projs[x["accession"]] = x
        pg += 1
    info = {}
    ip = os.path.join(CACHE, "crawl_info.json")
    if os.path.exists(ip):
        info = json.load(open(ip))
    return projs, info


def uniprot_ref_species(refresh=False):
    p = os.path.join(CACHE, "uniprot_reference_proteomes.tsv")
    if not os.path.exists(p) or refresh:
        url = ("https://rest.uniprot.org/proteomes/stream?query=reference:true&fields=upid,organism,organism_id,"
               "lineage&format=tsv")
        data = urllib.request.urlopen(urllib.request.Request(url), timeout=300).read()
        open(p, "wb").write(data)
    names = set()
    with open(p) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            for n in c[3].split(", "):
                names.add(n.strip().lower())
            names.add(re.sub(r"\s*\(.*", "", c[1]).strip().lower())
    return names


# --- stage 2: classification ------------------------------------------------------------------------------------
def names(lst, key="name"):
    return [x.get(key, "") if isinstance(x, dict) else str(x) for x in (lst or [])]


def vendor_of(inst_names):
    vs = set()
    fam = set()
    for n in inst_names:
        nl = n.lower()
        if "astral" in nl:
            vs.add("thermo"); fam.add("astral")
        elif "timstof" in nl.replace(" ", "") or "tims tof" in nl:
            vs.add("bruker"); fam.add("tims")
        elif re.search(r"orbitrap|q exactive|exactive|exploris|fusion|lumos|eclipse|tribrid|ltq|lcq|velos|elite", nl):
            vs.add("thermo")
            fam.add("orbitrap" if re.search(r"orbitrap|exactive|exploris|fusion|lumos|eclipse|tribrid|elite", nl) else "trap")
        elif re.search(r"maxis|impact|compact|bruker|amazon|hct|esquire|solarix|apex", nl):
            vs.add("bruker"); fam.add("bruker_qtof")
        else:
            vs.add("other:" + n)
    return vs, fam


def classify_project(p):
    r = {"accession": p["accession"]}
    inst = names(p.get("instruments"))
    qm = names(p.get("quantificationMethods"))
    et = names(p.get("experimentTypes"))
    kw = list(p.get("keywords") or []) + list(p.get("projectTags") or [])
    sw = names(p.get("softwares"))
    ptms = names(p.get("identifiedPTMStrings"))
    orgs = names(p.get("organisms"))
    title = p.get("title") or ""
    desc = p.get("projectDescription") or ""
    sp = p.get("sampleProcessingProtocol") or ""
    dp = p.get("dataProcessingProtocol") or ""
    struct = " | ".join([title] + kw + qm + et)
    prot = " ".join([desc, sp, dp])
    allt = struct + " " + prot + " " + " ".join(sw)
    year = int((p.get("submissionDate") or "0000")[:4])
    r.update(title=title[:160], submission_date=(p.get("submissionDate") or "")[:10], year=year,
             instruments=";".join(inst), organisms=";".join(orgs), n_organisms=len(orgs),
             quant_methods=";".join(qm), experiment_types=";".join(et), keywords=";".join(kw)[:200],
             software=";".join(sw), license=p.get("license", ""),
             lab=(p.get("labPIs") or [{}])[0].get("name", "") if p.get("labPIs") else "")
    reasons = []
    # vendor
    vs, fam = vendor_of(inst)
    if not inst:
        reasons.append("no_instrument")
    elif any(v.startswith("other") for v in vs):
        reasons.append("other_vendor")
    elif len(vs) > 1:
        reasons.append("mixed_vendor")
    elif fam == {"trap"}:
        reasons.append("ion_trap_only")
    vendor = next(iter(vs)) if len(vs) == 1 and not any(v.startswith("other") for v in vs) else ""
    family = ("astral" if "astral" in fam else "orbitrap" if "orbitrap" in fam else "tims" if "tims" in fam
              else "bruker_qtof" if "bruker_qtof" in fam else "trap" if "trap" in fam else "")
    r.update(vendor=vendor, family=family)
    # acquisition
    dia_s = bool(DIA_STRUCT.search(struct)) or any(re.search(r"dia-?nn|spectronaut|encyclopedia|openswath|dia-umpire", s, re.I) for s in sw)
    dia_w = bool(DIA_TEXT.search(prot))
    dda_t = bool(DDA_TEXT.search(allt))
    if dia_s:
        acq, conf = "DIA", "structured"
        reasons.append("dia")
    elif dia_w and dda_t:
        acq, conf = "mixed", "text"
        reasons.append("dda_dia_mixed")
    elif dia_w:
        acq, conf = "DIA", "text"
        reasons.append("dia")
    elif dda_t:
        acq, conf = "DDA", "text"
    else:
        acq, conf = "DDA", "assumed"
    r.update(acquisition=acq, acq_confidence=conf)
    # label
    sil = SILAC_RE.search(struct) or len(SILAC_RE.findall(prot)) >= 2
    tmt_s = TMT_RE.search(struct) is not None
    tmt_t = len(TMT_RE.findall(prot)) >= 2
    lfq_s = LFQ_RE.search(struct) is not None
    lfq_t = LFQ_RE.search(prot) is not None
    if sil:
        label, lconf = "SILAC", "structured" if SILAC_RE.search(struct) else "text"
        reasons.append("silac")
    elif tmt_s or tmt_t:
        label, lconf = "TMT", "structured" if tmt_s else "text"
    elif lfq_s or lfq_t:
        label, lconf = "LFQ", "structured" if lfq_s else "text"
    else:
        label, lconf = "UNK", ""
    # enrichment
    phos_ptm = any("phospho" in x.lower() for x in ptms)
    phos = bool(PHOS_A.search(struct)) and (bool(PHOS_B.search(allt)))
    if not phos and re.search(r"phosphoproteom|phosphopeptide enrich", prot, re.I) and PHOS_B.search(prot):
        phos = True
    enrich = "phospho" if phos else "none"
    r.update(label=label, label_confidence=lconf, enrichment=enrich, phospho_ptm_reported=int(phos_ptm))
    if XL_RE.search(struct) or len(XL_RE.findall(prot)) >= 2:
        reasons.append("xl")
    if MISC_DROP.search(struct) or len(MISC_DROP.findall(prot)) >= 2:
        reasons.append("topdown_maldi_meta")
    if label == "UNK" and enrich != "phospho":
        reasons.append("no_quant_or_phospho_evidence")
    r["year_bin"] = ("<=2014" if year <= 2014 else "2015-2017" if year <= 2017 else "2018-2020" if year <= 2020
                     else "2021-2023" if year <= 2023 else ">=2024")
    r["reasons"] = ",".join(sorted(set(reasons)))
    return r


def organism_ok(row, species):
    orgs = [o for o in row["organisms"].split(";") if o]
    if len(orgs) != 1:
        return False, "organism_count_%d" % len(orgs)
    nm = re.sub(r"\s*\(.*", "", orgs[0]).strip().lower()
    if nm in species:
        return True, ""
    return False, "no_reference_proteome"


def load_predictor_exclusions():
    ex = {}
    flag = {}
    if os.path.exists(PRED_TSV):
        with open(PRED_TSV) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                acc = r.get("accession", "")
                if not re.match(r"PXD\d+", acc):
                    continue
                is_train = r.get("role", "train") in ("train", "train_or_test", "unknown_role", "train_and_eval", "derived")
                tgt = ex if (r.get("confidence", "").lower() in ("high", "medium") and is_train) else flag
                tgt.setdefault(acc, set()).add(r["predictor"])
    return ex, flag


def _manifest_accessions(fname, required):
    path = os.path.join(MAN, fname)
    if not os.path.exists(path):
        if required:
            raise SystemExit("select_t4: required exclusion manifest missing: %s" % path)
        print("select_t4: WARNING optional exclusion manifest missing: %s" % path, file=sys.stderr)
        return set()
    with open(path, newline="") as f:
        accs = {r["accession"].strip() for r in csv.DictReader(f, delimiter="\t") if (r.get("accession") or "").strip()}
    if required and not accs:
        raise SystemExit("select_t4: required exclusion manifest has no accessions: %s" % path)
    return accs


def _datamd_train_accessions(path=None):
    """PXD accessions of DATA.md table rows whose 2nd column (role) is `train`."""
    path = path or DATA_MD
    out = set()
    if not os.path.exists(path):
        print("select_t4: WARNING DATA.md missing: %s" % path, file=sys.stderr)
        return out
    for line in open(path):
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(c) > 2 and c[1] == "train":
            out.update(re.findall(r"PXD\d{6}", line))
    return out


def load_exclusions():
    """acc -> set of exclusion reasons, from DEV_SEALED, train.tsv (+DATA.md role=train), the dev/sealed manifests,
    the D8/E9 stream manifest and the predictor-training manifest (high/medium train rows)."""
    ex = collections.defaultdict(set)
    for acc, tag in DEV_SEALED.items():
        ex[acc].add("dev_sealed:" + tag)
    train = _manifest_accessions(TRAIN_TSV, True)
    md_train = _datamd_train_accessions()
    if md_train - train:
        print("select_t4: WARNING DATA.md role=train accessions not in train.tsv: %s" % sorted(md_train - train), file=sys.stderr)
    for acc in train | md_train:
        ex[acc].add("train")
    for fname, req in DEV_MANIFESTS:
        for acc in _manifest_accessions(fname, req):
            ex[acc].add("dev_manifest:" + fname[:-4])
    for acc in _manifest_accessions(E9_D8_TSV, True):
        ex[acc].add("e9_d8_stream")
    if not os.path.exists(PRED_TSV):
        raise SystemExit("select_t4: required exclusion manifest missing: %s" % PRED_TSV)
    pred_ex, _ = load_predictor_exclusions()
    for acc in pred_ex:
        ex[acc].add("predictor_training")
    return dict(ex)


NEW_EXCLUSIONS = ("train", "dev_manifest:", "e9_d8_stream")  # added by review M7


def is_new_exclusion(reason):
    return reason == "train" or reason == "e9_d8_stream" or reason.startswith("dev_manifest:")


def assert_disjoint(selection_rows, exclusions=None):
    """Assert selection ∩ (train ∪ dev/sealed ∪ D8 ∪ predictor) = ∅.  selection_rows: iterable of accessions or dicts
    with an 'accession' key (e.g. csv.DictReader rows of t4_selection.tsv).  Returns the number of distinct accessions
    checked; raises AssertionError listing every overlap and its reasons."""
    ex = load_exclusions() if exclusions is None else exclusions
    accs = {(r["accession"] if isinstance(r, dict) else r).strip() for r in selection_rows}
    hit = {a: sorted(ex[a]) for a in sorted(accs) if a in ex}
    assert not hit, "T4 selection intersects excluded sets: %s" % hit
    return len(accs)


def assert_disjoint_file(path=None):
    path = path or os.path.join(MAN, "t4_selection.tsv")
    with open(path, newline="") as f:
        return assert_disjoint(csv.DictReader(f, delimiter="\t"))


def build_frame(projs):
    species = uniprot_ref_species()
    pred_ex, pred_flag = load_predictor_exclusions()
    excl = load_exclusions()
    rows, drops = [], collections.Counter()
    dropped_first = collections.Counter()
    for acc, p in projs.items():
        r = classify_project(p)
        ok, why = organism_ok(r, species)
        reasons = [x for x in r["reasons"].split(",") if x]
        if not ok:
            reasons.append(why)
        reasons.extend(sorted(excl.get(acc, ())))  # dev_sealed:*, train, dev_manifest:*, e9_d8_stream, predictor_training
        r["predictor_flag_lowconf"] = ";".join(sorted(pred_flag.get(acc, [])))
        r["exclude_reasons"] = ",".join(sorted(set(reasons)))
        r["eligible"] = int(not reasons)
        # pre-M7 eligibility (only the original rules + dev_sealed dict + predictor_training): the draw-mode frame
        r["eligible_base"] = int(not [x for x in reasons if not is_new_exclusion(x)])
        r["cell"] = "|".join([r["year_bin"], r["family"], r["label"], r["enrichment"]])
        for x in set(reasons):
            drops[x if x.startswith("dev_manifest:") else x.split(":")[0]] += 1
        if reasons:
            dropped_first[reasons[0].split(":")[0]] += 1
        rows.append(r)
    return rows, drops, excl


FRAME_COLS = ["accession", "title", "submission_date", "year", "year_bin", "vendor", "family", "instruments",
              "organisms", "acquisition", "acq_confidence", "label", "label_confidence", "enrichment",
              "phospho_ptm_reported", "quant_methods", "experiment_types", "keywords", "software", "license", "lab",
              "predictor_flag_lowconf", "cell", "eligible", "exclude_reasons"]


def write_tsv(path, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols])


# --- stage 3: selection ------------------------------------------------------------------------------------------
def raw_ok(vendor, name):
    n = name.lower()
    if re.search(r"\bdia\b|swath|diapasef|[_-]dia[_-]", n):
        return False
    if vendor == "thermo":
        return n.endswith(".raw")
    return bool(re.search(r"\.d(\.zip|\.tar(\.gz)?|\.7z|\.rar)?$", n)) or (n.endswith(".zip") and ".d" in n)


def project_files(acc, vendor, refresh):
    files, pg = [], 0
    while pg < 60:
        d = cached_json("files_%s_p%03d.json.gz" % (acc, pg), "%s/projects/%s/files?pageSize=100&page=%d" % (API, acc, pg), refresh)
        if not d:
            break
        files.extend(d)
        if len(d) < 100:
            break
        pg += 1
    out = []
    for x in files:
        if (x.get("fileCategory") or {}).get("value") != "RAW":
            continue
        nm, sz = x["fileName"], x.get("fileSizeBytes") or 0
        if raw_ok(vendor, nm) and sz >= MIN_FILE_BYTES:
            url = ""
            for loc in x.get("publicFileLocations", []):
                if "FTP" in loc.get("name", ""):
                    url = loc["value"]
            out.append({"file": nm, "bytes": sz, "url": url})
    return out, len(files)


def pps_allocate(cells, rng, n):
    keys = sorted(cells)
    w = {k: math.sqrt(len(cells[k])) for k in keys}
    interval = sum(w.values()) / n
    start = rng.random() * interval
    alloc, cum, nxt = collections.Counter(), 0.0, start
    for k in keys:
        cum += w[k]
        while nxt < cum - 1e-12:
            alloc[k] += 1
            nxt += interval
    return {k: min(v, len(cells[k])) for k, v in alloc.items()}


def select(rows, refresh, skip=None):
    """rows: the frame used for allocation and shuffling.  skip: acc -> reasons; such candidates are skipped while
    drawing (like a raw-check failure) without changing allocation or shuffle order (exclusion mode 'draw')."""
    skip = skip or {}
    rng = random.Random(SEED)
    elig = list(rows)
    cells = collections.defaultdict(list)
    for r in elig:
        cells[r["cell"]].append(r)
    alloc = pps_allocate(cells, rng, N_PROJECTS)
    chosen, lab_count, tried, failed = [], collections.Counter(), 0, collections.Counter()
    for k in sorted(alloc):
        cand = cells[k][:]
        random.Random(SEED + sum(map(ord, k))).shuffle(cand)
        got = 0
        for r in cand:
            if got >= alloc[k]:
                break
            if r["accession"] in skip:
                for x in skip[r["accession"]]:
                    failed["skipped_excluded:" + x] += 1
                continue
            if r["lab"] and lab_count[r["lab"]] >= MAX_PER_LAB:
                failed["lab_cap"] += 1
                continue
            tried += 1
            try:
                files, nall = project_files(r["accession"], r["vendor"], refresh)
            except Exception as e:  # noqa
                failed["api_error"] += 1
                continue
            if len(files) < 3:
                failed["fewer_than_3_usable_raw_files"] += 1
                r["raw_check"] = "fail:%d_usable_of_%d" % (len(files), nall)
                continue
            r["raw_check"] = "ok:%d_usable_of_%d" % (len(files), nall)
            frng = random.Random(SEED + int(re.sub(r"\D", "", r["accession"])))
            pick = frng.sample(files, min(MAX_FILES_PER_PROJECT, len(files)))
            chosen.append({"row": r, "files": sorted(pick, key=lambda x: x["file"]), "n_usable": len(files),
                           "n_files_all": nall, "alloc_cell": k})
            if r["lab"]:
                lab_count[r["lab"]] += 1
            got += 1
        if got < alloc[k]:
            failed["cell_underfilled:" + k] += alloc[k] - got
    # global caps
    trim = random.Random(SEED + 1)
    total = lambda: sum(len(c["files"]) for c in chosen)
    while total() > MAX_RUNS:
        mx = max(len(c["files"]) for c in chosen)
        c = trim.choice([c for c in chosen if len(c["files"]) == mx])
        c["files"].pop(trim.randrange(len(c["files"])))
    nbytes = lambda: sum(f["bytes"] for c in chosen for f in c["files"])
    dropped_for_bytes = 0
    while nbytes() > MAX_BYTES:
        c = max((c for c in chosen if len(c["files"]) > 1), key=lambda c: sum(f["bytes"] for f in c["files"]))
        c["files"].remove(max(c["files"], key=lambda f: f["bytes"]))
        dropped_for_bytes += 1
    return chosen, alloc, cells, failed, tried, dropped_for_bytes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["crawl", "frame", "select", "all"])
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--out-dir", default=None,
                    help="write t4_frame.tsv, t4_selection.tsv, summary.md and all_projects_classified.tsv here "
                         "(default: data/manifests/ + results/2026-10-26_T4_dryrun/ + the cache, as before)")
    ap.add_argument("--exclusion-mode", choices=["draw", "frame"], default="draw",
                    help="draw (default): allocate/shuffle on the pre-M7 frame and skip M7-excluded candidates while "
                         "drawing (reproduces the committed selection); frame: remove them from the frame first")
    ap.add_argument("--check", nargs="?", const=os.path.join(MAN, "t4_selection.tsv"), default=None,
                    help="only run assert_disjoint on a selection TSV (default data/manifests/t4_selection.tsv) and exit")
    a = ap.parse_args()
    if a.check:
        n = assert_disjoint_file(a.check)
        print("assert_disjoint OK: %d accessions in %s, none excluded" % (n, a.check))
        return
    man_out = a.out_dir or MAN
    res_out = a.out_dir or RES
    cls_out = a.out_dir or CACHE
    if a.stage in ("crawl", "all"):
        crawl(a.refresh)
        if a.stage == "crawl":
            return
    projs, info = load_cached_projects()
    if not projs:
        sys.exit("empty cache; run --stage crawl")
    rows, drops, excl = build_frame(projs)
    elig = [r for r in rows if r["eligible"]]
    base = [r for r in rows if r["eligible_base"]]
    print("frame without new (M7) exclusions:", len(base), "with:", len(elig), "exclusion mode:", a.exclusion_mode)
    write_tsv(os.path.join(man_out, "t4_frame.tsv"), FRAME_COLS, elig)
    write_tsv(os.path.join(cls_out, "all_projects_classified.tsv"), FRAME_COLS, rows)
    print("projects", len(rows), "eligible", len(elig))
    if a.stage == "frame":
        return
    if a.exclusion_mode == "frame":
        chosen, alloc, cells, failed, tried, dropped_bytes = select(elig, a.refresh)
    else:
        skip = {acc: sorted(x for x in rs if is_new_exclusion(x)) for acc, rs in excl.items()
                if any(is_new_exclusion(x) for x in rs)}
        chosen, alloc, cells, failed, tried, dropped_bytes = select(base, a.refresh, skip)
    sel_rows = []
    for c in chosen:
        r = c["row"]
        for f in c["files"]:
            sel_rows.append({"accession": r["accession"], "cell": r["cell"], "year_bin": r["year_bin"],
                             "vendor": r["vendor"], "family": r["family"], "label": r["label"],
                             "enrichment": r["enrichment"], "organism": r["organisms"], "lab": r["lab"],
                             "acq_confidence": r["acq_confidence"], "label_confidence": r["label_confidence"],
                             "file": f["file"], "bytes": f["bytes"], "url": f["url"],
                             "n_usable_raw_in_project": c["n_usable"]})
    assert_disjoint(sel_rows, excl)
    write_tsv(os.path.join(man_out, "t4_selection.tsv"),
              ["accession", "cell", "year_bin", "vendor", "family", "label", "enrichment", "organism", "lab",
               "acq_confidence", "label_confidence", "file", "bytes", "url", "n_usable_raw_in_project"], sel_rows)
    write_summary(rows, drops, elig, chosen, alloc, cells, failed, tried, dropped_bytes, info, res_out,
                  frame_path=os.path.relpath(os.path.join(man_out, "t4_frame.tsv"), ROOT), n_base=len(base),
                  mode=a.exclusion_mode)


def write_summary(rows, drops, elig, chosen, alloc, cells, failed, tried, dropped_bytes, info, res_dir=RES,
                  frame_path="data/manifests/t4_frame.tsv", n_base=None, mode="draw"):
    os.makedirs(res_dir, exist_ok=True)
    nb = sum(f["bytes"] for c in chosen for f in c["files"])
    nr = sum(len(c["files"]) for c in chosen)
    gb = nb / 1e9
    spg = {"thermo": 40000, "bruker": 25000}
    spectra = sum(f["bytes"] / 1e9 * spg[c["row"]["vendor"]] for c in chosen for f in c["files"])
    cnt = lambda key, rs: collections.Counter(r[key] for r in rs)
    L = []
    L.append("# T4 selector dry run (metadata only), %s\n" % time.strftime("%Y-%m-%d"))
    L.append("No raw file downloaded. Script: `scripts/data/select_t4.py` (seed 11). Cache: `data/cache/pride_meta/`.\n")
    L.append("## Frame\n")
    L.append("- PRIDE projects crawled: %d unique accessions (API header total_records: %s; crawl %s)." % (
        len(rows), info.get("total_records_header"), info.get("crawl_date")))
    L.append("- **Metadata-eligible frame: %d projects** (`%s`). Raw-file presence is verified only "
             "for selection candidates: %d candidates checked, failures/skips %s." % (len(elig), frame_path, tried,
                                                                                    dict(failed) or "none"))
    L.append("- Exclusion mode `%s`. Frame without the M7 exclusions (train, dev_manifest:*, e9_d8_stream): %s projects; "
             "with them: %d. In `draw` mode allocation and within-cell shuffles use the former and M7-excluded "
             "candidates are skipped while drawing; in `frame` mode they use the latter." % (mode, n_base, len(elig)))
    L.append("- Exclusion reasons (a project can have several; counts per reason over all %d projects):\n" % len(rows))
    for k, v in drops.most_common():
        L.append("  - %s: %d" % (k, v))
    L.append("\n### Frame composition\n")
    for key in ("year_bin", "family", "label", "enrichment", "acq_confidence", "label_confidence"):
        L.append("- %s: %s" % (key, dict(sorted(cnt(key, elig).items()))))
    L.append("- cells (year x family x label x enrichment) non-empty: %d" % len(cells))
    L.append("\n## Selection\n")
    L.append("- Projects: %d (target 40, allowed 30-50); runs: %d (cap 300); bytes: %.1f GB = %.3f TB (cap 1 TB); "
             "files dropped for byte cap: %d." % (len(chosen), nr, gb, gb / 1000, dropped_bytes))
    L.append("- Cells hit: %d. Allocation rule: sqrt-weighted systematic PPS (see script header)." % len(set(c["alloc_cell"] for c in chosen)))
    sel = [c["row"] for c in chosen]
    for key in ("year_bin", "family", "label", "enrichment", "vendor"):
        L.append("- selected by %s: %s" % (key, dict(sorted(cnt(key, sel).items()))))
    L.append("- selected acq_confidence: %s (non-'text' rows rely on defaults; 'assumed' = no DIA or DDA wording found)" %
             dict(cnt("acq_confidence", sel)))
    L.append("\n### Estimated ingest time\n")
    L.append("Assumes %d spectra (MS2) per GB raw in total (40k/GB Thermo, 25k/GB Bruker; planning figure) -> %.2e spectra." %
             (round(spectra / max(gb, 1e-9)), spectra))
    for rate in (2000, 5000):
        L.append("- at %d spectra/s: %.1f h (single process)" % (rate, spectra / rate / 3600))
    L.append("\n## DDA vs DIA classification and error modes\n")
    L.append("""PRIDE has no machine-readable acquisition-mode field. Rules (in order): (1) strong DIA = DIA/SWATH/diaPASEF/
PRM/SRM terms or DIA software (DIA-NN, Spectronaut, OpenSWATH, EncyclopeDIA) in title, keywords, tags, experiment types, quantification
methods or software -> dropped; (2) DIA wording only in protocol text -> dropped unless DDA wording is also present, in which
case labelled `mixed` and dropped; (3) DDA wording (data-dependent, DDA, TopN, shotgun) -> DDA/text; (4) none of these -> DDA/assumed.
Error modes: (a) DDA-labelled projects that are DDA library + DIA quantification with the DIA only described in files (kept as
`assumed`/`text` DDA; file-name filter removes DIA/SWATH names only); (b) projects mentioning DIA as a comparison or in the
Introduction-style description (dropped = false negatives of the frame, harmless for validity but shrink it); (c) software lists
rarely record DIA tools so structured recall is low; (d) 'shotgun proteomics' experiment type is applied to DIA too, so it is not used as
DDA evidence; (e) projects with both DDA and DIA runs pass only if no DIA wording exists, hence residual DIA contamination possible
in `assumed` projects (quantified above via acq_confidence). Label errors: TMT/LFQ mentions in protocols can describe a comparison or
earlier work (TMT needs >=2 text mentions or one structured one; SILAC needs >=2); LFQ is often undeclared, so many true LFQ projects are
absent from the frame (`no_quant_or_phospho_evidence`). Phospho needs phospho wording and an enrichment term. Instrument lists are
submitter-chosen, may be incomplete (e.g. Astral listed as Exploris) and may mix vendors (dropped). Organism is matched by species
name to UniProt reference-proteome lineages (strain/synonym mismatches lose a few projects; multi-organism projects dropped).
Year is submissionDate, not acquisition date.

Open uncertainties: spectra per GB; whether the Bruker `.d.zip` ingest path exists; ground-truth DDA/DIA per selected project
should be confirmed from raw-file headers at ingest time (a cheap check: Thermo scan filter / tdf acquisition mode) and a project
failing it replaced from the same cell. Predictor-training exclusions: see `data/manifests/predictor_training_accessions.tsv`
(re-run `--stage frame`/`all` after it changes; ProteomeTools accessions are already excluded).
""")
    with open(os.path.join(res_dir, "summary.md"), "w") as f:
        f.write("\n".join(L))
    print("\n".join(L[:30]))


if __name__ == "__main__":
    main()
