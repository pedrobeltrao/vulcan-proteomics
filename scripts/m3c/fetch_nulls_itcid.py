#!/usr/bin/env python3
"""Fetch + ingest + validate the IT-CID yeast null candidates (null_itcid_candidates.tsv; rule in
select_nulls_itcid.py): projects in seeded order, a project is accepted if both its runs are >= 95% IT-CID MS2 with
no MS3; stop at 3 accepted projects.  Raw -> data/raw/null_itcid/<acc>/, stores -> $VULCAN_3C_WORK/store_null.
Writes data/manifests/null_itcid.tsv (role dev, set DN3c), null_itcid.sha256, null_itcid_fetch.log.json.
"""
import csv, json, os, sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fetch_train_itcid as F  # noqa: E402

TARGET = 3


def dest(r):
    return os.path.join(F.DATA, "raw", "null_itcid", r["accession"], r["filename"])


def store_path(r):
    return os.path.join(F.WORK, "store_null", os.path.splitext(r["filename"])[0] + ".vspec.arrow")


F.dest = dest
F.store_path = store_path


def proc(r):
    try:
        F.download(r)
        v = F.validate(dict(r, stratum="LF"), F.ingest(r))
    except Exception as e:  # noqa: BLE001
        v = {"accept": False, "reason": f"error: {e}"}
    print(r["cand"], r["accession"], v["reason"], flush=True)
    return v


def main():
    cands = list(csv.DictReader(open(os.path.join(F.MAN, "null_itcid_candidates.tsv")), delimiter="\t"))
    import select_train_itcid as S
    n_dn2 = S.dn2_guard([r["accession"] for r in cands], ["Saccharomyces cerevisiae"])
    print("DN2 guard passed (checked against", n_dn2, "DN2-walk accessions)", flush=True)
    accs = []
    for r in cands:
        if r["accession"] not in accs:
            accs.append(r["accession"])
    log, acc_ok = {}, []
    i = 0
    while len(acc_ok) < TARGET and i < len(accs):
        batch = accs[i:i + TARGET - len(acc_ok) + 1]
        i += len(batch)
        rows = [r for r in cands if r["accession"] in batch]
        with ThreadPoolExecutor(6) as ex:
            for r, v in zip(rows, ex.map(proc, rows)):
                log[r["cand"]] = dict(v, accession=r["accession"], filename=r["filename"])
        for a in batch:
            if all(log[r["cand"]]["accept"] for r in cands if r["accession"] == a):
                acc_ok.append(a)
    acc_ok = [a for a in accs if a in acc_ok][:TARGET]
    json.dump({"accepted_projects": acc_ok, "runs": log}, open(os.path.join(F.MAN, "null_itcid_fetch.log.json"), "w"),
              indent=1)
    rows = [r for r in cands if r["accession"] in acc_ok]
    cols = ["id", "accession", "filename", "url", "expected_size", "instrument", "year", "lab", "role", "set_id"]
    with open(os.path.join(F.MAN, "null_itcid.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k, r in enumerate(rows):
            fh.write("\t".join([f"DN3c_{k:02d}"] + [r[c] for c in cols[1:8]] + ["dev", "DN3c"]) + "\n")
    lines = [f"{F.sha(dest(r))}  {os.path.relpath(dest(r), F.DATA)}" for r in rows]
    open(os.path.join(F.MAN, "null_itcid.sha256"), "w").write("\n".join(lines) + "\n")
    print("accepted", acc_ok)


if __name__ == "__main__":
    main()
