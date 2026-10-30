#!/usr/bin/env python3
"""Select IT-CID wrong-species NULL runs for the IT-CID G2-original / G5 gates (tree 3c; DESIGN §6.4 G2: "strata:
instrument class only"; the OT-HCD null sets DN / D7 / D1 are not IT-CID).  Stdlib only; written and committed
BEFORE any download.  These runs are gate nulls only (role dev, set DN3c): never trained on, never used for
selection, only searched against the human DB inside G2 mixtures.

Species: S. cerevisiae (eukaryote with a reviewed reference proteome already in data/fasta, as the M3 G2 yeast
nulls; not a free-living microbe of the A2 DN2 ALLOW list, and A2 excludes S. cerevisiae from DN2 anyway).
Frame rows: single organism S. cerevisiae, DDA, LFQ, enrichment none, Thermo, instrument in {LTQ Orbitrap (XL/Velos/
Elite), Orbitrap Fusion, Fusion Lumos}; excluded: every accession of every data/manifests/*.tsv except t4_frame.tsv
(incl. D7 PXD031616, DN, the 3c training candidates); BAD_KW as the training selection.
Order: sorted by accession, shuffled with random.Random(23); per project 2 raw files 0.15-1.5 GB (rng choice).
Acceptance after ingest: >= 95% of MS2 IT+CID, no MS3.  Walk in order until 3 projects (6 runs) are accepted.
Homolog filter, mixture design and test: exactly M3's G2 (scripts/m3/g2_nulls.py, gates.py), hosts = D10cf.

    python3 scripts/m3c/select_nulls_itcid.py  -> data/manifests/null_itcid_candidates.tsv
"""
import csv, json, os, random, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import select_train_itcid as S  # noqa: E402

TARGET_PROJECTS = 3
N_CAND = 12
INSTR = re.compile(r"^(LTQ Orbitrap( XL| Velos| Elite)?|Orbitrap Fusion|Orbitrap Fusion Lumos)$")
SEED = 23


def main():
    excl = S.excluded_accessions()
    excl |= {r["accession"] for r in csv.DictReader(open(os.path.join(S.MAN, "train_itcid_candidates.tsv")),
                                                    delimiter="\t")}
    frame = list(csv.DictReader(open(os.path.join(S.MAN, "t4_frame.tsv")), delimiter="\t"))
    cand = [r for r in frame if r["eligible"] == "1" and r["organisms"].startswith("Saccharomyces cerevisiae")
            and ";" not in r["organisms"] and r["acquisition"] == "DDA" and r["label"] == "LFQ"
            and r["enrichment"] == "none" and r["vendor"] == "thermo" and INSTR.match(r["instruments"])
            and r["accession"] not in excl and not S.BAD_KW.search(r["title"] + " " + r["keywords"] + " "
                                                                  + r["experiment_types"])]
    cand.sort(key=lambda r: r["accession"])
    rng = random.Random(SEED)
    rng.shuffle(cand)
    rows, log = [], {"candidates": len(cand), "tried": []}
    for r in cand:
        if len({x["accession"] for x in rows}) >= N_CAND:
            break
        try:
            files = S.pride_files(r["accession"])
        except Exception as e:  # noqa: BLE001
            log["tried"].append([r["accession"], f"api error {e}"]); continue
        raws = sorted([f for f in files if f["fileName"].lower().endswith(".raw") and not S.BAD_FILE["LF"].search(
            f["fileName"]) and 1.5e8 <= f["fileSizeBytes"] <= 1.5e9 and S.ftp_url(f)], key=lambda f: f["fileName"])
        if len(raws) < 2:
            log["tried"].append([r["accession"], f"{len(raws)} raw files in size range"]); continue
        for k, f in enumerate(rng.sample(raws, 2)):
            rows.append({"cand": f"N{len({x['accession'] for x in rows}):02d}{k}", "accession": r["accession"],
                         "filename": f["fileName"], "url": S.ftp_url(f), "expected_size": str(f["fileSizeBytes"]),
                         "instrument": r["instruments"], "year": r["year"], "lab": r["lab"]})
        log["tried"].append([r["accession"], "candidate"])
        print(r["accession"], r["instruments"], flush=True)
    cols = ["cand", "accession", "filename", "url", "expected_size", "instrument", "year", "lab"]
    with open(os.path.join(S.MAN, "null_itcid_candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    json.dump(log, open(os.path.join(S.MAN, "null_itcid_select.log.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
