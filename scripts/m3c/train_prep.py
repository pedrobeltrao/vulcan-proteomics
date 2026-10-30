#!/usr/bin/env python3
"""Search (IT-CID production settings, family prod index) + L3 features for every ingested IT-CID training
candidate store (idempotent; only candidates that pass validation are used by train_static_itcid.py).

    python scripts/m3c/train_prep.py [--only-accepted]
"""
import argparse, csv, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-accepted", action="store_true")
    a = ap.parse_args()
    cands = list(csv.DictReader(open(os.path.join(C.ROOT, "data", "manifests", "train_itcid_candidates.tsv")),
                                delimiter="\t"))
    log = json.load(open(os.path.join(C.ROOT, "data", "manifests", "train_itcid_fetch.log.json")))
    out = C.psm_dir("train")
    by_fam = {}
    for r in cands:
        st = C.train_store(r)
        if not os.path.exists(st):
            continue
        if a.only_accepted and not log.get(r["cand"], {}).get("accept"):
            continue
        if log.get(r["cand"]) and not log[r["cand"]]["accept"]:
            continue
        run = os.path.splitext(r["filename"])[0]
        if os.path.exists(os.path.join(out, run + ".psms.parquet")):
            continue
        by_fam.setdefault(C.train_family(r), []).append(st)
    for fam, stores in by_fam.items():
        subprocess.run([sys.executable, os.path.join(C.ROOT, "scripts", "m3c", "search.py"), "--family", fam,
                        "--stores", *stores, "--out", os.path.join(out, "_tmp_" + fam)], check=True)
        for f in os.listdir(os.path.join(out, "_tmp_" + fam)):
            if f.endswith(".psms.parquet"):
                os.replace(os.path.join(out, "_tmp_" + fam, f), os.path.join(out, f))
        os.replace(os.path.join(out, "_tmp_" + fam, "stats.json"), os.path.join(out, f"search_stats_{fam}_"
                                                                                 f"{len(os.listdir(out))}.json"))
        os.rmdir(os.path.join(out, "_tmp_" + fam))
    # features per family (features.py skips runs that already have them)
    fams = {}
    for r in cands:
        run = os.path.splitext(r["filename"])[0]
        if os.path.exists(os.path.join(out, run + ".psms.parquet")):
            fams.setdefault(C.train_family(r), []).append(run)
    import features as F
    import shutil
    for fam, runs in fams.items():
        d = os.path.join(out, "_feat_" + fam)
        os.makedirs(d, exist_ok=True)
        todo = [x for x in runs if not os.path.exists(os.path.join(out, x + ".feat.parquet"))]
        for x in todo:
            shutil.copyfile(os.path.join(out, x + ".psms.parquet"), os.path.join(d, x + ".psms.parquet"))
        if todo:
            F.run_dir(d, fam)
            for x in todo:
                os.replace(os.path.join(d, x + ".feat.parquet"), os.path.join(out, x + ".feat.parquet"))
            os.replace(os.path.join(d, "feat_info.json"), os.path.join(out, f"feat_info_{fam}_{len(todo)}_{todo[0]}.json"))
        shutil.rmtree(d)
        print(fam, "features", len(todo), flush=True)


if __name__ == "__main__":
    main()
