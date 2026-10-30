#!/usr/bin/env python3
"""Fetch + ingest + validate the IT-CID training candidates (data/manifests/train_itcid_candidates.tsv, written by
select_train_itcid.py) and write the accepted training manifest, fold map and sha256.

Acceptance rule: exactly as pre-registered in select_train_itcid.py (docstring).  Deviation fixed before the first
download (logged): a TMT-stratum candidate whose FILE NAME declares TMTpro (TMT16/TMT18/16plex/18plex/TMTpro) is
skipped without download, because it would fail the pre-specified TMT 6/10/11 label check anyway (saves bandwidth;
cannot change which runs are accepted).
Per stratum the candidates are walked in seeded order; the first TARGET runs that pass are accepted.

Raw files -> data/raw/train_itcid/<acc>/ (resumable, size-verified); a copy on node-local NVMe is decoded with
vulcan.io.convert (LESSONS: never let RawFileReader mmap Lustre) into $VULCAN_3C_WORK/store_train.
Outputs: data/manifests/train_itcid.tsv (role=train), train_itcid_folds.tsv, train_itcid.sha256,
train_itcid_fetch.log.json.

    python scripts/m3c/fetch_train_itcid.py [--jobs 6]
"""
import argparse, csv, hashlib, json, os, re, shutil, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, "data")
MAN = os.path.join(DATA, "manifests")
WORK = os.environ.get("VULCAN_3C_WORK", "/tmp/vulcan_3c")
UA = {"User-Agent": "vulcan-3c-train-fetch (research)"}
TARGET = {"TMT": 6, "LF": 7}
TMTPRO_NAME = re.compile(r"TMT ?1[68]|1[68] ?plex|TMTpro", re.I)
TMT_REP = np.array([126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
                    130.141145, 131.138180, 131.144500])
# TMTpro-only marker reporters.  BUG FIX (2026-10-27, after the first 3 TMT candidates were rejected, before any
# run was accepted; no scorer/dev data involved): the first rule also used 132N/132C/133N/133C, but those coincide
# with the 13C isotope peaks of the TMT 131N/131C (and 132) reporters (131.1382 + 1.00335 = 132.1415), so every
# TMT10/11 run "had TMTpro reporters".  Now only 134N/134C/135N (no TMT 6/10/11 reporter or isotope within 0.003 Th)
# count, and only with intensity >= 5% of the most intense reporter-region peak of that MS3.
PRO_ONLY = np.array([134.148245, 134.154565, 135.151600])


def dest(r):
    return os.path.join(DATA, "raw", "train_itcid", r["accession"], r["filename"])


def store_path(r):
    return os.path.join(WORK, "store_train", os.path.splitext(r["filename"])[0] + ".vspec.arrow")


def download(r):
    path, size = dest(r), int(r["expected_size"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) == size:
        return "skip"
    part = path + ".part"
    for attempt in range(8):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            req = urllib.request.Request(r["url"], headers={**UA, **({"Range": f"bytes={have}-"} if have else {})})
            with urllib.request.urlopen(req, timeout=120) as resp, open(part, "ab" if have else "wb") as fh:
                if have and resp.status != 206:
                    fh.truncate(0)
                while True:
                    b = resp.read(1 << 22)
                    if not b:
                        break
                    fh.write(b)
            if os.path.getsize(part) == size:
                os.replace(part, path)
                return "ok"
        except Exception as e:  # noqa: BLE001
            print(f"  retry {r['cand']} {attempt}: {e}", file=sys.stderr, flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{r['cand']} download failed")


def ingest(r):
    out = store_path(r)
    if os.path.exists(out):
        return out
    stage = os.path.join(WORK, "rawstage")
    os.makedirs(stage, exist_ok=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    src = os.path.join(stage, os.path.basename(dest(r)))
    shutil.copyfile(dest(r), src)
    try:
        subprocess.run([sys.executable, "-m", "vulcan.io.convert", src, out, "--threads", "4"], check=True, cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": ROOT}, capture_output=True, timeout=3600)
    finally:
        os.remove(src)
    return out


def validate(r, store):
    from vulcan.io import store as ST
    t = ST.read_table(store, columns=["ms_level", "activation", "analyzer"])
    lvl = t.column("ms_level").to_numpy()
    m2 = lvl == 2
    act = t.column("activation").to_numpy(zero_copy_only=False)[m2]
    ana = t.column("analyzer").to_numpy(zero_copy_only=False)[m2]
    n2, n3 = int(m2.sum()), int((lvl == 3).sum())
    frac = float(np.mean((act == 1) & (ana == 2))) if n2 else 0.0
    rec = {"n_ms2": n2, "n_ms3": n3, "it_cid_fraction": frac}
    ok, why = frac >= 0.95, [] if frac >= 0.95 else [f"IT-CID fraction {frac:.3f} < 0.95"]
    if r["stratum"] == "LF":
        if n3 > 0:
            ok = False; why.append(f"{n3} MS3 spectra (isobaric-label method?)")
    else:
        if n3 < 0.5 * n2:
            ok = False; why.append(f"MS3/MS2 = {n3 / max(n2, 1):.2f} < 0.5 (no SPS-MS3)")
        else:
            tt = ST.read_table(store)
            tt = tt.filter(tt.column("ms_level").to_numpy() == 3)
            pk = ST.peaks_numpy(tt)
            n_any = n_pro = 0
            for i in range(min(len(pk.offsets) - 1, 5000)):
                mz = pk.mz[pk.offsets[i]:pk.offsets[i + 1]]
                it = pk.intensity[pk.offsets[i]:pk.offsets[i + 1]]
                sel = (mz > 125.5) & (mz < 135.5)
                lo, li = mz[sel], it[sel]
                if not len(lo):
                    continue
                d = np.abs(lo[:, None] - TMT_REP[None, :]).min(0) <= 0.003
                strong = lo[li >= 0.05 * li.max()]
                p = (np.abs(strong[:, None] - PRO_ONLY[None, :]).min(0) <= 0.003) if len(strong) else np.zeros(1, bool)
                if d.any() or p.any():
                    n_any += 1
                    n_pro += bool(p.any())
            rec.update(ms3_with_reporter=n_any, ms3_with_tmtpro_only_reporter=n_pro)
            if n_any == 0:
                ok = False; why.append("no TMT reporter in MS3")
            elif n_pro > 0.1 * n_any:
                ok = False; why.append(f"TMTpro reporters in {n_pro}/{n_any} MS3 (not TMT 6/10/11)")
    rec["accept"] = bool(ok)
    rec["reason"] = "; ".join(why) if why else "pass"
    return rec


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def process(r):
    try:
        st = download(r)
        s = ingest(r)
        v = validate(r, s)
        v["download"] = st
    except Exception as e:  # noqa: BLE001
        v = {"accept": False, "reason": f"error: {e}"}
    print(r["cand"], r["accession"], v["reason"], flush=True)
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    cands = list(csv.DictReader(open(os.path.join(MAN, "train_itcid_candidates.tsv")), delimiter="\t"))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import select_train_itcid as S
    S.dn2_guard([r["accession"] for r in cands], ["Homo sapiens", "Mus musculus"])
    logp = os.path.join(MAN, "train_itcid_fetch.log.json")
    log = json.load(open(logp)) if os.path.exists(logp) else {}
    accepted = {}
    for stratum in ("TMT", "LF"):
        cs = [r for r in cands if r["stratum"] == stratum]
        for r in cs:
            if stratum == "TMT" and TMTPRO_NAME.search(r["filename"]) and r["cand"] not in log:
                log[r["cand"]] = {"accept": False, "reason": "skipped without download: file name declares TMTpro"}
        i = 0
        while True:
            acc = [r for r in cs if log.get(r["cand"], {}).get("accept")]
            # accepted = first TARGET passing in seeded order, provided every earlier candidate is resolved
            done = []
            for r in cs:
                if r["cand"] not in log:
                    break
                if log[r["cand"]]["accept"]:
                    done.append(r)
            if len(done) >= TARGET[stratum]:
                accepted[stratum] = done[:TARGET[stratum]]
                break
            todo = [r for r in cs if r["cand"] not in log]
            if not todo:
                accepted[stratum] = done
                print("WARNING: candidates exhausted for", stratum, file=sys.stderr)
                break
            need = TARGET[stratum] - len(acc)
            batch = todo[:max(a.jobs, need + 2)]
            with ThreadPoolExecutor(a.jobs) as ex:
                for r, v in zip(batch, ex.map(process, batch)):
                    log[r["cand"]] = dict(v, accession=r["accession"], filename=r["filename"])
                    json.dump(log, open(logp, "w"), indent=1)   # persist per result (a reader segfault lost a batch)
            json.dump(log, open(logp, "w"), indent=1)
            i += 1
    json.dump(log, open(logp, "w"), indent=1)
    # manifest + fold map (round-robin within stratum, acceptance order)
    rows = []
    for stratum in ("TMT", "LF"):
        for k, r in enumerate(accepted[stratum]):
            rows.append(dict(r, id=f"TI{len(rows):02d}", fold=str(k % 3), role="train",
                             label="TMT" if stratum == "TMT" else "none"))
    cols = ["id", "accession", "filename", "url", "expected_size", "species", "instrument", "year", "lab", "stratum",
            "label", "fold", "role"]
    with open(os.path.join(MAN, "train_itcid.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(MAN, "train_itcid_folds.tsv"), "w") as fh:
        fh.write("accession\tfold\tstratum\tspecies\n")
        for r in rows:
            fh.write(f"{r['accession']}\t{r['fold']}\t{r['stratum']}\t{r['species']}\n")
    lines = [f"{sha(dest(r))}  {os.path.relpath(dest(r), DATA)}" for r in rows]
    open(os.path.join(MAN, "train_itcid.sha256"), "w").write("\n".join(lines) + "\n")
    print(json.dumps({s: [r["accession"] for r in v] for s, v in accepted.items()}))


if __name__ == "__main__":
    main()
