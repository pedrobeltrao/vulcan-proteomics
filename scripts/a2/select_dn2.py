#!/usr/bin/env python3
"""Amendment A2.9.1-A2.9.3 + A2.3: fresh-null (DN2) selector.  DO NOT RUN before the tag prereg-v0.4.2 exists and
contains this file (enforced by require_tag()); runs exactly once (P4).

Procedure (every constant from config/a2.yaml "selection"; order fixed in selection.step_order):
  A2.9.1 exclusions   union of the `accession` columns of every data/manifests/*.tsv except t4_frame.tsv (read at run
                      time) ∪ config/a2_excluded_accessions.tsv (the 300 accessions at adoption; the run-time set must
                      contain it); non-PXD frame rows dropped; previously used null species and training species
                      dropped by organism prefix.
  A2.3 step 1         frame metadata filter (single organism on ALLOW, DDA, Thermo, LFQ, enrichment none, OT-HCD
                      instrument regex of A1; title/keywords/experiment types not matching BAD_KW = A1 list + A2 list).
  A2.9.3 walk         candidates sorted by accession, shuffled with random.Random(23); for each candidate in order:
                        1. genus slot: skip if an ADMITTED project already holds its genus (Sulfolobus/Saccharolobus
                           share one slot);
                        2. A2.3 step 2 protocol text + SDRF (PRIDE project record; metadata only);
                        3. list PRIDE files; >= 3 usable RAW (.raw, 0.3-1.5 GB, name not matching BAD_FILE) else skip;
                        4. pick 3 with random.Random(23 + int(accession digits)) (sample from the name-sorted list);
                        5. UniProt reference proteome (rule in config); none -> skip;
                        6. download the 3 RAWs (sha256 recorded);
                        7. own-species search + marker screen (A2.2a, A2.3 step 3) -> exclude (logged with counts) or
                        8. admit.py (A2.2) -> label-free counts; the project is admitted and takes its genus slot;
                        9. stopping rule (A2.9.3) on the label-free counts of admitted projects; stop or continue.
                      The walk ends when the rule is met, the cap (16) is reached, or candidates run out.
Outputs: data/manifests/dn2.tsv (+ dn2.sha256), data/manifests/dn2_select.log.json (every candidate tried, reason,
label-free counts only).  No per-label (target/decoy) quantity is computed, printed or stored (A2.9.4).

    python scripts/a2/select_dn2.py            # after the tag only
The dry-run on a synthetic project table is tests/a2/test_select_dryrun.py (no network, no files).
"""
from __future__ import annotations

import csv
import glob
import hashlib
import json
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Callable

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates import a2_config  # noqa: E402
from vulcan.gates.a2_admission import stopping_state  # noqa: E402

MAN = os.path.join(ROOT, "data", "manifests")


# ------------------------------------------------------------------------------------------------ A2.9.1
def manifest_accessions(man_dir: str = MAN, exclude=("t4_frame.tsv",)) -> set:
    out = set()
    for p in sorted(glob.glob(os.path.join(man_dir, "*.tsv"))):
        if os.path.basename(p) in exclude:
            continue
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                a = (r.get("accession") or "").strip()
                if a:
                    out.add(a)
    return out


def adoption_accessions(cfg: dict, root: str = ROOT) -> set:
    e = cfg["selection"]["exclusions"]
    with open(os.path.join(root, e["at_adoption_file"]), newline="") as fh:
        s = {r["accession"].strip() for r in csv.DictReader(fh, delimiter="\t")}
    if len(s) != e["at_adoption_count"]:
        raise RuntimeError("adoption exclusion list count mismatch")
    if hashlib.sha256("\n".join(sorted(s)).encode()).hexdigest() != e["at_adoption_sha256"]:
        raise RuntimeError("adoption exclusion list sha256 mismatch")
    return s


def load_exclusions(cfg: dict, man_dir: str = MAN, root: str = ROOT) -> set:
    """Run-time manifest union ∪ the adoption list (untracked manifests such as d7b/d10cf/e9_d8 may be absent on
    another checkout; the union keeps them, so the run-time exclusion set always contains the adoption set)."""
    run_time = manifest_accessions(man_dir, tuple(cfg["selection"]["exclusions"]["except"]))
    return run_time | adoption_accessions(cfg, root)


# ------------------------------------------------------------------------------------------------ ALLOW / organism
def genus_slot(organism: str, cfg: dict):
    """-> genus slot name if the (single) organism is on the ALLOW list, else None."""
    w = organism.split()
    if not w:
        return None
    genus = w[0]
    species = w[1].lower() if len(w) > 1 else ""
    for a in cfg["selection"]["allow"]:
        if a["genus"] != genus:
            continue
        if "species" in a and species not in [s.lower() for s in a["species"]]:
            return None
        if species in [s.lower() for s in a.get("deny_species", [])]:
            return None
        return a.get("genus_group", genus)
    return None


def organism_excluded(organism: str, cfg: dict) -> str | None:
    e = cfg["selection"]["exclusions"]
    low = organism.lower()
    for s in e["previous_null_species"] + e["training_species"]:
        if low.startswith(s.lower()):
            return f"excluded species {s}"
    return None


def regexes(cfg: dict) -> dict:
    s = cfg["selection"]
    bad = re.compile(s["bad_kw_a1"] + "|" + s["bad_kw_a2_add"], re.I)
    prot = re.compile(s["bad_kw_a1"] + "|" + s["bad_kw_a2_add"] + "|" + s["protocol_kw_add"], re.I)
    return {"bad_kw": bad, "protocol": prot, "instr": re.compile(s["frame_filter"]["instruments_regex"]),
            "bad_file": re.compile(s["bad_file_regex"], re.I)}


def frame_candidates(frame_rows, cfg: dict, excl: set) -> tuple[list, dict]:
    """A2.3 step 1 + A2.9.1 + A2.9.2 on the frame.  Returns (candidates in accession order, drop reason counts)."""
    ff = cfg["selection"]["frame_filter"]
    rx = regexes(cfg)
    out, drops = [], {}

    def drop(why):
        drops[why] = drops.get(why, 0) + 1

    for r in frame_rows:
        acc = r["accession"].strip()
        org = (r.get("organisms") or "").strip()
        if not re.fullmatch(r"PXD\d+", acc):
            drop("non-PXD")
        elif acc in excl:
            drop("excluded accession")
        elif not org or ";" in org:
            drop("not single organism")
        elif organism_excluded(org, cfg):
            drop("excluded species")
        elif genus_slot(org, cfg) is None:
            drop("not on ALLOW")
        elif r.get("acquisition") != ff["acquisition"] or r.get("vendor") != ff["vendor"] \
                or r.get("label") != ff["label"] or r.get("enrichment") != ff["enrichment"]:
            drop("frame fields")
        elif not rx["instr"].match(r.get("instruments", "")):
            drop("instrument")
        elif rx["bad_kw"].search(" ".join([r.get("title", ""), r.get("keywords", ""), r.get("experiment_types", "")])):
            drop("BAD_KW")
        else:
            out.append(r)
    out.sort(key=lambda r: r["accession"])
    return out, drops


def shuffled(cands: list, cfg: dict) -> list:
    c = sorted(cands, key=lambda r: r["accession"])
    random.Random(cfg["seeds"]["walk"]).shuffle(c)
    return c


# ------------------------------------------------------------------------------------------------ A2.3 step 2
def protocol_reject(record: dict, sdrf_rows, cfg: dict) -> str | None:
    """None if the protocol text / SDRF do not exclude the project; else the reason.  Missing text/SDRF never
    excludes."""
    rx = regexes(cfg)["protocol"]
    for f in cfg["selection"]["protocol_fields"]:
        v = record.get(f)
        if isinstance(v, list):
            v = " ".join(map(str, v))
        if v and rx.search(str(v)):
            return f"protocol text ({f}): {rx.search(str(v)).group(0)}"
    if sdrf_rows:
        orgs = {(r.get("characteristics[organism]") or "").strip().lower() for r in sdrf_rows}
        orgs.discard("")
        if len(orgs) > cfg["selection"]["sdrf"]["max_organisms"]:
            return f"SDRF organisms {sorted(orgs)}"
        allowed = {x.lower() for x in cfg["selection"]["sdrf"]["allowed_organism_parts"]}
        parts = {(r.get("characteristics[organism part]") or "").strip().lower() for r in sdrf_rows}
        parts.discard("")
        bad = parts - allowed
        if bad:
            return f"SDRF organism part {sorted(bad)}"
    return None


def usable_raws(files: list[dict], cfg: dict) -> list[dict]:
    lo, hi = cfg["selection"]["raw_size_bytes"]
    bf = regexes(cfg)["bad_file"]
    return sorted([f for f in files if f["fileName"].lower().endswith(".raw") and not bf.search(f["fileName"])
                   and lo <= f["fileSizeBytes"] <= hi and f.get("url")], key=lambda f: f["fileName"])


def pick_files(raws: list[dict], acc: str, cfg: dict) -> list[dict]:
    k = cfg["selection"]["runs_per_project"]
    pick = random.Random(cfg["seeds"]["file_pick_base"] + int(re.sub(r"\D", "", acc))).sample(raws, k)
    return sorted(pick, key=lambda f: f["fileName"])


# ------------------------------------------------------------------------------------------------ walk
@dataclass
class WalkIO:
    """Injected I/O.  Production: PrideIO; dry-run tests: fakes.  Every callable sees one project at a time."""
    project_record: Callable[[str], dict]
    sdrf: Callable[[str], list]
    list_files: Callable[[str], list]
    proteome: Callable[[dict], dict | None]
    download: Callable[[str, list], list]                      # -> rows with sha256
    marker_and_admit: Callable[[str, list, dict], dict]       # -> {"marker": {...}, "counts": {...} | None}
    log: list = field(default_factory=list)


def walk(cands: list, cfg: dict, io: WalkIO) -> dict:
    taken, admitted, rows = set(), [], []
    state = stopping_state([], cfg["selection"]["stopping"])
    for r in cands:
        acc, org = r["accession"], r["organisms"].strip()
        slot = genus_slot(org, cfg)
        if slot in taken:
            io.log.append({"accession": acc, "organism": org, "result": "genus slot taken", "slot": slot})
            continue
        why = protocol_reject(io.project_record(acc), io.sdrf(acc), cfg)
        if why:
            io.log.append({"accession": acc, "organism": org, "result": "protocol", "reason": why})
            continue
        raws = usable_raws(io.list_files(acc), cfg)
        if len(raws) < cfg["selection"]["runs_per_project"]:
            io.log.append({"accession": acc, "organism": org, "result": f"{len(raws)} usable raw"})
            continue
        pick = pick_files(raws, acc, cfg)
        prot = io.proteome(r)
        if prot is None:
            io.log.append({"accession": acc, "organism": org, "result": "no reference proteome"})
            continue
        got = io.download(acc, pick)
        res = io.marker_and_admit(acc, got, prot)
        mk = res["marker"]
        if mk["exclude"]:
            io.log.append({"accession": acc, "organism": org, "result": "marker screen exclude",
                           "trigger_members": mk["trigger_members"],
                           "n_trigger_members_detected": mk["n_trigger_members_detected"],
                           "n_panel_members_detected": mk["n_panel_members_detected"]})
            continue
        counts = res["counts"]
        taken.add(slot)
        admitted.append(counts)
        for f in got:
            rows.append({"accession": acc, "organism": org, "upid": prot["upid"], **f})
        state = stopping_state(admitted, cfg["selection"]["stopping"])
        io.log.append({"accession": acc, "organism": org, "result": "admitted", "slot": slot, "counts": counts,
                       "marker": {k: mk[k] for k in ("n_panel_members_detected", "n_trigger_members_detected",
                                                     "detected_members")},
                       "stopping": state})
        if state["stop"]:
            break
    return {"admitted": [x["accession"] for x in io.log if x["result"] == "admitted"], "rows": rows,
            "final_state": state, "log": io.log}


# ------------------------------------------------------------------------------------------------ guard
def require_tag(tag: str, paths=("scripts/a2", "vulcan/gates", "config/a2.yaml", "config/a2_panels.tsv",
                                  "config/a2_excluded_accessions.tsv")):
    """The selector may only run on a checkout that contains the tag, with the A2 code and config unmodified."""
    def git(*a):
        return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True)
    if git("rev-parse", "-q", "--verify", f"refs/tags/{tag}").returncode != 0:
        raise SystemExit(f"select_dn2: tag {tag} does not exist: A2 is not adopted; do not list DN2")
    if git("merge-base", "--is-ancestor", tag, "HEAD").returncode != 0:
        raise SystemExit(f"select_dn2: HEAD does not contain {tag}")
    d = git("diff", "--quiet", tag, "--", *paths)
    st = git("status", "--porcelain", "--", *paths).stdout.strip()
    if d.returncode != 0 or st:
        raise SystemExit("select_dn2: A2 code/config differ from the tag; refusing to run")


# ------------------------------------------------------------------------------------------------ reference proteome
# D16 (RESEARCH_TREE; results/2026-10-27_A2_dn2/DEVIATIONS.md items 11-13): two implementation fixes.
#  (1) UniProt REST no longer matches `proteome_type:1`; the reference filter is `reference:true` (as select_t4.py).
#  (2) The match is at species level: a genus-only organism (one word, or a placeholder epithet sp./sp/spp./spp)
#      has no reference proteome (stricter reading).
PLACEHOLDER_EPITHETS = {"sp", "sp.", "spp", "spp."}


def genus_species(organism: str) -> str | None:
    w = organism.split()
    if len(w) < 2 or w[1].lower() in PLACEHOLDER_EPITHETS:
        return None
    return " ".join(w[:2])


def proteome_query(gs: str) -> str:
    return f'(organism_name:"{gs}") AND (reference:true)'


def match_reference_proteome(organism: str, hits: list) -> dict | None:
    """hits: [upid, organism name, taxid, protein count] rows of the UniProt reference-proteome search.  First hit whose
    name starts with the full frame organism, else with "genus species"; genus-only organisms -> None."""
    gs = genus_species(organism)
    if gs is None:
        return None
    for key in (organism.strip().lower(), gs.lower()):
        for upid, name, taxid, n in hits:
            if name.lower().startswith(key):
                return {"upid": upid, "organism": name, "taxid": taxid, "n_proteins": int(n)}
    return None


# ------------------------------------------------------------------------------------------------ production I/O
class PrideIO:  # pragma: no cover  (network; runs once, after the tag)
    PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3"

    def __init__(self, cfg, raw_dir):
        self.cfg, self.raw_dir = cfg, raw_dir

    def _get(self, url):
        import time
        import urllib.request
        for k in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "vulcan-a2"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    return r.read()
            except Exception:  # noqa: BLE001
                if k == 3:
                    raise
                time.sleep(3 * (k + 1))

    def project_record(self, acc):
        return json.loads(self._get(f"{self.PRIDE}/projects/{acc}"))

    def _files(self, acc):
        out, page = [], 0
        while True:
            got = json.loads(self._get(f"{self.PRIDE}/projects/{acc}/files?pageSize=100&page={page}"))
            out += got
            if len(got) < 100:
                return out
            page += 1

    @staticmethod
    def _url(f):
        for loc in f.get("publicFileLocations", []):
            if loc.get("value", "").startswith("ftp://"):
                return loc["value"].replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
        return None

    def sdrf(self, acc):
        fs = [f for f in self._files(acc) if f["fileName"].lower().endswith(".sdrf.tsv") and self._url(f)]
        rows = []
        for f in fs:
            txt = self._get(self._url(f)).decode(errors="replace")
            rows += list(csv.DictReader(txt.splitlines(), delimiter="\t"))
        return rows

    def list_files(self, acc):
        return [{"fileName": f["fileName"], "fileSizeBytes": f["fileSizeBytes"], "url": self._url(f)}
                for f in self._files(acc) if f.get("fileCategory", {}).get("value") == "RAW"]

    def proteome(self, row):
        import urllib.parse
        org = row["organisms"].strip()
        gs = genus_species(org)
        if gs is None:                      # D16: genus-only organism -> no reference proteome
            return None
        q = urllib.parse.quote_plus(proteome_query(gs))
        txt = self._get(f"https://rest.uniprot.org/proteomes/search?query={q}&fields=upid,organism,organism_id,"
                        f"protein_count&format=tsv&size=50").decode()
        hits = [line.split("\t") for line in txt.splitlines()[1:] if line.strip()]
        return match_reference_proteome(org, hits)

    def download(self, acc, files):
        out = []
        d = os.path.join(self.raw_dir, acc)
        os.makedirs(d, exist_ok=True)
        for f in files:
            p = os.path.join(d, f["fileName"])
            subprocess.run(["curl", "-sSf", "-C", "-", "-o", p, f["url"]], check=True)
            if os.path.getsize(p) != f["fileSizeBytes"]:
                raise RuntimeError(f"size mismatch {p}")
            out.append({"filename": f["fileName"], "url": f["url"], "expected_size": f["fileSizeBytes"],
                        "sha256": a2_config.sha256_file(p), "path": p})
        return out

    def marker_and_admit(self, acc, files, prot):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import admit
        return admit.run_project(acc, files, prot, self.cfg)


def main():  # pragma: no cover
    cfg = a2_config.load()
    require_tag(cfg["tag_required"])
    a2_config.check_panel_files(cfg)
    a2_config.check_seqclf(cfg)
    excl = load_exclusions(cfg)
    with open(os.path.join(ROOT, cfg["selection"]["frame"]), newline="") as fh:
        frame = list(csv.DictReader(fh, delimiter="\t"))
    cands, drops = frame_candidates(frame, cfg, excl)
    io_ = PrideIO(cfg, os.path.join(ROOT, "data", "raw", "dn2"))
    w = walk(shuffled(cands, cfg), cfg, WalkIO(io_.project_record, io_.sdrf, io_.list_files, io_.proteome,
                                                io_.download, io_.marker_and_admit))
    cols = ["id", "accession", "filename", "url", "expected_size", "purpose", "role", "set_id", "member"]
    with open(os.path.join(MAN, "dn2.tsv"), "w") as fh, open(os.path.join(MAN, "dn2.sha256"), "w") as sh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(w["rows"]):
            fh.write("\t".join([f"DN2_{i:02d}", r["accession"], r["filename"], r["url"], str(r["expected_size"]),
                                f"A2 fresh null (DN2): {r['organism']}; own proteome {r['upid']}", "dev", "DN2",
                                ""]) + "\n")
            sh.write(f"{r['sha256']}  {r['accession']}/{r['filename']}\n")
    json.dump({"seed": cfg["seeds"]["walk"], "n_frame_candidates": len(cands), "frame_drops": drops,
               "n_exclusions": len(excl), **{k: w[k] for k in ("admitted", "final_state", "log")}},
              open(os.path.join(MAN, "dn2_select.log.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
