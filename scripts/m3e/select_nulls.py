#!/usr/bin/env python3
"""A1.3 null-source selection (DESIGN §6.0 A1.3, §6.4 G2): wrong-species OT-HCD DDA runs for the G2 / tail-G2 / G1 /
G3 / G5 mixtures, recorded as dev set DN (data/manifests/g2null.tsv).  Stdlib only; metadata-only.

Rule (fixed before any of these runs is searched):
  * frame: data/manifests/t4_frame.tsv; single organism from ALLOW (prokaryotes / unicellular eukaryotes, far from
    human); DDA, LFQ, enrichment none, Thermo, single instrument matching INSTR (Q Exactive family / Exploris = OT-HCD);
    title/keywords not matching BAD_KW;
  * hard exclusions: T4 selection, every dev/sealed/tier-S accession, the 13 M3 training projects, every predictor-
    training accession, the E9/D8 stream (same manifests as scripts/e9/select_e9.py + e9_d8.tsv);
  * already-used null organisms excluded (P. furiosus = D1, S. cerevisiae = D7) and training species (human, mouse,
    E. coli), Arabidopsis (arm C);
  * projects sorted by accession, shuffled with random.Random(11); walk in order; keep a project if its organism is
    not yet taken and it has >= RUNS usable RAW files (.raw, 0.3-1.5 GB, name not matching BAD_FILE); pick RUNS files
    with random.Random(11 + accession digits); stop at N_PROJ projects.
  * reference proteome: first UniProt proteome hit (free-text organism query) whose organism name starts with the
    frame's genus+species; recorded with the row.  Used ONLY by the homolog filter (A1.3).

    python3 scripts/m3e/select_nulls.py [--n-proj 4] [--runs 3]
"""
import argparse, csv, json, os, random, re, sys, time, urllib.parse, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
sys.path.insert(0, os.path.join(ROOT, "scripts", "e9"))
SEED = 11
ALLOW = ["Staphylococcus aureus", "Bacillus cereus", "Pseudomonas aeruginosa", "Synechocystis", "Mycobacterium tuberculosis",
         "Streptomyces", "Myxococcus xanthus", "Klebsiella pneumoniae", "Francisella tularensis", "Neisseria gonorrhoeae",
         "Synechococcus", "Cupriavidus necator", "Trypanosoma brucei", "Plasmodium falciparum", "Toxoplasma gondii",
         "Dictyostelium discoideum", "Chlamydomonas reinhardtii", "Leishmania donovani", "Tetrahymena thermophila"]
INSTR = re.compile(r"^(Q Exactive( HF| HF-X| Plus)?|Orbitrap Exploris (480|240))$")
BAD_KW = re.compile(r"HLA|MHC|immunopept|cross-?link|XL-MS|glyco|SILAC|DIA\b|SWATH|data.independent|phospho|single.cell|"
                    r"TMT|iTRAQ|ubiquitin|acetyl|N-?termin|top-down|crosslink|metaproteom|interactome|AP-MS|BioID|"
                    r"pull.?down|immunoprecip|labell?ing|dimethyl|PRM|SRM|targeted|infect|host|co-?culture", re.I)
BAD_FILE = re.compile(r"phos|TMT|blank|wash|QC|std|standard|DIA|SWATH|ptm|enrich|\bIP\b|_IP_|hela|BSA|calib|human|mouse", re.I)
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100&page=%d"


def get_json(url):
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "vulcan-m3e"}), timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            if k == 3:
                raise
            time.sleep(3 * (k + 1))


def pride_files(acc):
    out, page = [], 0
    while True:
        got = json.loads(get_json(PRIDE % (acc, page)))
        out += got
        if len(got) < 100:
            return out
        page += 1


def ftp_url(f):
    for loc in f.get("publicFileLocations", []):
        if loc.get("value", "").startswith("ftp://"):
            return loc["value"].replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
    return None


def accessions(name):
    p = os.path.join(MAN, name)
    if not os.path.exists(p):
        return set()
    return {r["accession"] for r in csv.DictReader(open(p), delimiter="\t") if r.get("accession")}


def proteome(org):
    gs = " ".join(org.split()[:2])
    q = urllib.parse.quote_plus(gs)
    txt = get_json(f"https://rest.uniprot.org/proteomes/search?query={q}&fields=upid,organism,organism_id,protein_count"
                   f"&format=tsv&size=25").decode()
    for line in txt.splitlines()[1:]:
        upid, name, taxid, n = line.split("\t")
        if name.lower().startswith(gs.lower()):
            return {"upid": upid, "organism": name, "taxid": taxid, "n_proteins": int(n)}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-proj", type=int, default=4)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(MAN, "g2null.tsv"))
    a = ap.parse_args()
    excl = set()
    for n in ("t4_selection.tsv", "tier_s.tsv", "dev_sealed.tsv", "dev_sealed_extra.tsv", "train.tsv",
              "predictor_training_accessions.tsv", "e9_d8.tsv"):
        excl |= accessions(n)
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    cand = [r for r in frame if r["organisms"] and ";" not in r["organisms"]
            and any(r["organisms"].startswith(x) for x in ALLOW) and r["acquisition"] == "DDA" and r["label"] == "LFQ"
            and r["enrichment"] == "none" and r["vendor"] == "thermo" and INSTR.match(r["instruments"])
            and not BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])
            and r["accession"] not in excl]
    cand.sort(key=lambda r: r["accession"])
    random.Random(SEED).shuffle(cand)
    log = {"seed": SEED, "allow": ALLOW, "n_candidates": len(cand), "tried": []}
    taken, rows = set(), []
    for r in cand:
        if len({x["accession"] for x in rows}) >= a.n_proj:
            break
        gs = " ".join(r["organisms"].split()[:2])
        if gs in taken:
            log["tried"].append([r["accession"], gs, "organism taken"])
            continue
        try:
            files = pride_files(r["accession"])
        except Exception as e:  # noqa: BLE001
            log["tried"].append([r["accession"], gs, f"api error {e}"])
            continue
        raws = sorted([f for f in files if f.get("fileCategory", {}).get("value") == "RAW"
                       and f["fileName"].lower().endswith(".raw") and not BAD_FILE.search(f["fileName"])
                       and 3e8 <= f["fileSizeBytes"] <= 1.5e9 and ftp_url(f)], key=lambda f: f["fileName"])
        if len(raws) < a.runs:
            log["tried"].append([r["accession"], gs, f"{len(raws)} usable raw"])
            continue
        prot = proteome(r["organisms"])
        if prot is None:
            log["tried"].append([r["accession"], gs, "no UniProt proteome"])
            continue
        pick = random.Random(SEED + int(re.sub(r"\D", "", r["accession"]))).sample(raws, a.runs)
        taken.add(gs)
        log["tried"].append([r["accession"], gs, f"selected {a.runs} of {len(raws)}", prot])
        for f in sorted(pick, key=lambda f: f["fileName"]):
            rows.append({"id": f"DN{len(rows):02d}", "accession": r["accession"], "filename": f["fileName"],
                         "url": ftp_url(f), "expected_size": str(f["fileSizeBytes"]),
                         "purpose": f"A1.3 wrong-species G2 null source: {r['organisms']} ({r['instruments']}); "
                                    f"homolog filter proteome {prot['upid']} (taxid {prot['taxid']})",
                         "role": "dev", "set_id": "DN", "member": ""})
        print(r["accession"], gs, prot["upid"], flush=True)
    cols = ["id", "accession", "filename", "url", "expected_size", "purpose", "role", "set_id", "member"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in cols) + "\n")
    json.dump(log, open(os.path.join(MAN, "g2null_select.log.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
