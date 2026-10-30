#!/usr/bin/env python3
"""M7: select ADDITIONAL project-disjoint OT-HCD training runs for the encoder arm (DESIGN §4 L3 cross-encoder).
Stdlib only.  Same frame and criteria as scripts/m3/select_train.py (M3 training set), with a NEW seed (23 =
SEEDS[1]; M3 used 11) and a wider exclusion list fixed before any download:

  * every accession in any data/manifests/*.tsv except t4_frame.tsv (t4_selection, tier_s, dev_sealed*, d7b, d9-d12,
    d10cf, d4_ref, dn2, g2null (DN), e9_d8 (D8 incl. the holdout), train_itcid*, null_itcid_candidates,
    predictor_training_accessions, and the existing M3 train.tsv -- the latter are re-used as-is, not re-drawn),
  * every accession appearing in the A2 DN2 walk (data/manifests/dn2_select.log.json, results/.../selection_log.tsv),
  * yeast (arm M / T5) and Arabidopsis (arm C) by species restriction (human / mouse / E. coli only).

Output: data/manifests/train_m7.tsv (M3 rows + new rows; id TRxx kept for M3 rows, TMxx for new), train_m7_folds.tsv
(project -> fold 0..4, round-robin over the seeded order; M3 projects keep their M3 fold mod 5), train_m7_select.log.json.

    python3 scripts/m7/select_train_m7.py [--human 24 --mouse 3 --ecoli 3]
"""
import argparse, csv, glob, json, os, random, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import select_train as ST  # noqa: E402  (frame filters, PRIDE helpers)

SEED = 23
N_FOLDS = 5


def all_excluded():
    excl = set()
    for f in glob.glob(os.path.join(MAN, "*.tsv")):
        if os.path.basename(f) == "t4_frame.tsv":
            continue
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if "accession" in (rd.fieldnames or []):
                excl |= {r["accession"] for r in rd if r.get("accession")}
    log = json.load(open(os.path.join(MAN, "dn2_select.log.json")))
    excl |= {e.get("accession") for e in log.get("log", []) if e.get("accession")}
    sl = os.path.join(ROOT, "results", "2026-10-27_A2_dn2", "selection_log.tsv")
    if os.path.exists(sl):
        excl |= ST.accessions(sl, "accession")
    return excl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--human", type=int, default=24)
    ap.add_argument("--mouse", type=int, default=3)
    ap.add_argument("--ecoli", type=int, default=3)
    a = ap.parse_args()
    m3 = list(csv.DictReader(open(os.path.join(MAN, "train.tsv")), delimiter="\t"))
    excl = all_excluded()
    frame = list(csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t"))
    rng = random.Random(SEED)
    rows, log = [], {"seed": SEED, "excluded_accessions": len(excl)}
    for sp, n in (("human", a.human), ("mouse", a.mouse), ("ecoli", a.ecoli)):
        cand = [r for r in frame if r["eligible"] == "1" and r["organisms"] == ST.SPECIES[sp] and r["acquisition"] == "DDA"
                and r["label"] == "LFQ" and r["enrichment"] == "none" and r["vendor"] == "thermo"
                and ST.INSTR.match(r["instruments"]) and r["accession"] not in excl and r["year"] >= "2017"
                and not ST.BAD_KW.search(r["title"] + " " + r["keywords"] + " " + r["experiment_types"])]
        cand.sort(key=lambda r: r["accession"])
        rng.shuffle(cand)
        log[sp] = {"candidates": len(cand), "tried": []}
        got = 0
        for r in cand:
            if got >= n:
                break
            try:
                files = ST.pride_files(r["accession"])
            except Exception as e:  # noqa: BLE001
                log[sp]["tried"].append([r["accession"], f"api error {e}"]); continue
            raws = sorted([f for f in files if f["fileName"].lower().endswith(".raw") and not ST.BAD_FILE.search(f["fileName"])
                           and 3e8 <= f["fileSizeBytes"] <= 1.3e9 and ST.ftp_url(f)], key=lambda f: f["fileName"])
            if not raws:
                log[sp]["tried"].append([r["accession"], "0 raw files in size range"]); continue
            f = rng.choice(raws)
            sha1 = next((c.get("value") for c in f.get("checksum", []) if isinstance(c, dict)), None) \
                if isinstance(f.get("checksum"), list) else f.get("checksum")
            rows.append({"id": f"TM{len(rows):02d}", "accession": r["accession"], "filename": f["fileName"],
                         "url": ST.ftp_url(f), "expected_size": str(f["fileSizeBytes"]), "species": sp,
                         "instrument": r["instruments"], "year": r["year"], "lab": r["lab"], "role": "train",
                         "pride_sha1": str(sha1 or "")})
            log[sp]["tried"].append([r["accession"], "selected"])
            got += 1
    projs = []
    for r in rows:
        if r["accession"] not in projs:
            projs.append(r["accession"])
    fold = {p: i % N_FOLDS for i, p in enumerate(projs)}
    for r in m3:
        fold[r["accession"]] = int(r["fold"]) % N_FOLDS
    out = []
    for r in m3:
        out.append({**r, "pride_sha1": "", "source": "m3_train"})
    for r in rows:
        out.append({**r, "source": "m7_new"})
    for r in out:
        r["fold"] = str(fold[r["accession"]])
    cols = ["id", "accession", "filename", "url", "expected_size", "species", "instrument", "year", "lab", "fold",
            "role", "source", "pride_sha1"]
    with open(os.path.join(MAN, "train_m7.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(MAN, "train_m7_folds.tsv"), "w") as fh:
        fh.write("accession\tfold\tspecies\tsource\n")
        for r in out:
            fh.write(f"{r['accession']}\t{r['fold']}\t{r['species']}\t{r['source']}\n")
    json.dump(log, open(os.path.join(MAN, "train_m7_select.log.json"), "w"), indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else {"candidates": v["candidates"],
                                                             "selected": sum(t[1] == 'selected' for t in v["tried"])})
                      for k, v in log.items()}, indent=1))
    print(f"{len(rows)} new runs from {len(projs)} projects, {sum(int(r['expected_size']) for r in rows)/1e9:.1f} GB")


if __name__ == "__main__":
    main()
