"""T2 driver datasets: raw files, search libraries and truth, per accession.

D2   = the M4 E8 dev set exactly (PXD000138 dev pools 58/60/66/72, PXD009449 Ymod_Phospho, PXD007058 pools 1-4); used to
       reproduce results/2026-10-27_M4_ptm/e8.json.
D2T  = the T2-shaped dev analogue: PXD000138 dev pools + the four dev PXD009449 pools (Ymod_Phospho, Kmod_Dimethyl,
       Kmod_Formyl, Rmod_Citrullin), i.e. the same two components as T2 including non-phospho kit pools.
T2   = sealed: PXD000138 pools 76/25/24/61 + PXD009449 Kmod_Acetyl/Succinyl/Ubiquitinyl/Crotonyl; RAW only from
       $SEALED_FILES; libraries and truth from data/truth/sealed/ (read only in sealed mode, sha256-verified first).

Truth rules are the M4 §0 declarations (scripts/m4/truth.py logic, parameterised by file paths, outputs identical on D2),
plus the literal §6.5 pool rule (a pool containing positional isomers of one sequence is excluded as a whole).
"""
from __future__ import annotations

import csv
import glob
import hashlib
import os
import sys
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
M4 = os.path.join(ROOT, "scripts", "m4")
if M4 not in sys.path:
    sys.path.insert(0, M4)

HUMAN = ["data/fasta/human_sp_canonical.fasta", "data/fasta/contaminants_frankenfield2022.fasta"]
SEALED_DIRS = (os.path.join(ROOT, "data", "sealed"), os.path.join(ROOT, "data", "truth", "sealed"))

DEV_POOL_WELL = {"58": "E10", "60": "E12", "66": "F6", "72": "F12"}
T2_POOL_WELL = {"76": "G4", "25": "C1", "24": "B12", "61": "F1"}
PHOSPHO_KITS = {"Ymod_Phospho"}        # PXD009449 pools whose modified peptides carry phospho (single-Y kit truth)


def sealed_raw_names() -> list[tuple[str, str]]:
    """(accession, filename) of the T2 sealed RAW files from data/manifests/dev_sealed.tsv (manifest only)."""
    out = []
    for r in csv.DictReader(open(_p("data/manifests/dev_sealed.tsv")), delimiter="\t"):
        if r["role"] == "sealed" and r["set_id"] == "T2" and r["filename"].endswith(".raw"):
            out.append((r["accession"], r["filename"]))
    return out


def pxd009449_pool(raw: str) -> str:
    """'02330a_GG1_3990_07_PTM_TrainKit_Kmod_Acetyl_200fmol_3xHCD_R1.raw' -> 'Kmod_Acetyl'."""
    b = os.path.basename(raw)
    return b.split("TrainKit_", 1)[1].split("_200fmol", 1)[0]


def run_name(raw: str) -> str:
    return os.path.basename(raw)[:-4]


def _p(rel):
    return os.path.join(ROOT, rel)


def dataset(name: str, sealed_files: list[str] | None = None) -> dict:
    """Return {'accessions': {acc: {'raws': [...], 'kind': ...}}, 'role': 'dev'|'sealed', ...}."""
    if name == "D2":
        acc = {
            "PXD000138": {"raws": sorted(glob.glob(_p("data/raw/PXD000138/*.raw"))), "kind": "marx",
                          "pool_well": DEV_POOL_WELL, "pool_fasta_dir": _p("data/truth/PXD000138")},
            "PXD009449": {"raws": sorted(glob.glob(_p("data/raw/PXD009449/*Ymod_Phospho*3xHCD*.raw"))), "kind": "ptmkit",
                          "evidence_dir": _p("data/truth/PXD009449"), "lib_tag": "Ymod"},
            "PXD007058": {"raws": sorted(glob.glob(_p("data/raw/PXD007058/*.raw"))), "kind": "ferries"},
        }
        role = "dev"
    elif name == "D2T":
        acc = {
            "PXD000138": {"raws": sorted(glob.glob(_p("data/raw/PXD000138/*.raw"))), "kind": "marx",
                          "pool_well": DEV_POOL_WELL, "pool_fasta_dir": _p("data/truth/PXD000138")},
            "PXD009449": {"raws": sorted(glob.glob(_p("data/raw/PXD009449/*3xHCD*.raw"))), "kind": "ptmkit",
                          "evidence_dir": _p("data/truth/PXD009449"), "lib_tag": "dev4"},
        }
        role = "dev"
    elif name == "T2":
        if not sealed_files:
            raise SystemExit("t2: T2 needs $SEALED_FILES (run through bench/sealed_run.py --dataset T2)")
        sf = sorted(sealed_files)
        acc = {
            "PXD000138": {"raws": [f for f in sf if os.sep + "PXD000138" + os.sep in f], "kind": "marx",
                          "pool_well": T2_POOL_WELL, "pool_fasta_dir": _p("data/truth/sealed/PXD000138")},
            "PXD009449": {"raws": [f for f in sf if os.sep + "PXD009449" + os.sep in f], "kind": "ptmkit",
                          "evidence_dir": _p("data/truth/sealed/PXD009449"), "lib_tag": "T2Kmod"},
        }
        if len(acc["PXD000138"]["raws"]) != 4 or len(acc["PXD009449"]["raws"]) != 4 or len(sf) != 8:
            raise SystemExit("t2: T2 must have exactly 4 PXD000138 + 4 PXD009449 sealed RAW files")
        role = "sealed"
    else:
        raise SystemExit(f"t2: unknown dataset {name}")
    for a, d in acc.items():
        if not d["raws"]:
            raise SystemExit(f"t2: no raw files for {name}/{a}")
        if role == "dev":
            for f in d["raws"] + [d.get("pool_fasta_dir", ""), d.get("evidence_dir", "")]:
                if f and os.path.realpath(f).startswith(SEALED_DIRS):
                    raise SystemExit("t2: dev dataset refers to a sealed path")
        if d["kind"] == "marx":
            got = sorted(run_name(r) for r in d["raws"])
            if got != sorted(d["pool_well"]):
                raise SystemExit(f"t2: PXD000138 runs {got} do not match the pool map")
    return {"name": name, "role": role, "accessions": acc}


# ------------------------------------------------------------------------------------------------ input verification
def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def manifest_hashes() -> dict:
    out = {}
    for m in ("dev_sealed.sha256", "derived_truth.sha256", "fasta.sha256", "tier_s.sha256"):
        p = _p(os.path.join("data/manifests", m))
        if not os.path.exists(p):
            continue
        for l in open(p):
            if l.strip():
                h, rel = l.rstrip("\n").split("  ", 1)
                out[os.path.normpath(os.path.join(ROOT, "data", rel))] = h
    return out


def library_inputs(ds: dict) -> list[str]:
    """Every truth/library file the dataset reads (hashed before use; contents never printed)."""
    files = [_p(f) for f in HUMAN]
    for a, d in ds["accessions"].items():
        if d["kind"] == "marx":
            files += [os.path.join(d["pool_fasta_dir"], f"pool_{w}.fasta") for w in sorted(d["pool_well"].values())]
        elif d["kind"] == "ptmkit":
            files += [os.path.join(d["evidence_dir"], f"{pxd009449_pool(r)}_evidence.txt") for r in d["raws"]]
        elif d["kind"] == "ferries":
            files += sorted(glob.glob(_p("data/truth/PXD007058/*.zip")))
    return files


def verify_inputs(files: list[str], strict: bool) -> dict:
    """sha256 of every input against the committed manifests. strict (sealed): any file without a manifest hash or with
    a mismatch aborts."""
    want = manifest_hashes()
    rep = {}
    for f in files:
        h = sha256(f)
        w = want.get(os.path.normpath(f))
        ok = (w == h) if w else None
        rep[os.path.relpath(f, ROOT)] = {"sha256": h, "manifest": ok}
        if strict and ok is not True:
            raise SystemExit("t2: REFUSED input hash mismatch")
    return rep


# ------------------------------------------------------------------------------------------------ libraries / FASTA
def evidence_sequences(paths: list[str]) -> list[str]:
    s = set()
    for p in paths:
        for x in csv.DictReader(open(p), delimiter="\t"):
            s.add(x["Sequence"])
    return sorted(s)


def ferries_sequences() -> list[str]:
    import common as C  # scripts/m4 (dev files only)
    return C.ferries_sequences()


def write_fasta(ds: dict, out: str) -> str:
    """human SwissProt canonical + contaminants + the dataset's libraries, in the M4 d2_fasta order (D2: identical
    bytes to /tmp/vulcan_m4/d2_phospho.fasta)."""
    with open(out + ".partial", "w") as fh:
        for f in HUMAN:
            fh.write(open(_p(f)).read().rstrip("\n") + "\n")
        a = ds["accessions"]
        if "PXD000138" in a:
            d = a["PXD000138"]
            for f in sorted(glob.glob(os.path.join(d["pool_fasta_dir"], "pool_*.fasta"))):
                if os.path.basename(f)[5:-6] in d["pool_well"].values():
                    fh.write(open(f).read().rstrip("\n") + "\n")
        if "PXD009449" in a:
            d = a["PXD009449"]
            ev = [os.path.join(d["evidence_dir"], f"{pxd009449_pool(r)}_evidence.txt") for r in d["raws"]]
            for i, s in enumerate(evidence_sequences(ev)):
                fh.write(f">lib_PXD009449_{d['lib_tag']}_{i}\n{s}\n")
        if "PXD007058" in a:
            for i, s in enumerate(ferries_sequences()):
                fh.write(f">lib_PXD007058_{i}\n{s}\n")
    os.replace(out + ".partial", out)
    return out


# ------------------------------------------------------------------------------------------------ truth (M4 §0)
def truth_marx(d: dict) -> dict:
    from truth import infer_design  # scripts/m4/truth.py (unchanged)
    truth, rep = {}, {}
    for pool, well in sorted(d["pool_well"].items()):
        f = os.path.join(d["pool_fasta_dir"], f"pool_{well}.fasta")
        seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
        P, var, units = infer_design(seq)
        sets = {j: sorted({u[j] for u in units}) for j in var}
        slots = [j for j, s in sets.items() if set(s) == set("STY")]
        r = {"well": well, "unit_len": P, "n_units": len(units), "n_variable_positions": len(var),
             "n_acceptor_slots": len(slots)}
        if len(slots) != 1:
            r["dropped"] = "not exactly one {S,T,Y} slot"
            rep[pool] = r
            continue
        slot = slots[0]
        seen = defaultdict(set)
        for u in units:
            seen[u].add(slot)
        iso = [u for u, s in seen.items() if len(s) > 1]
        if iso:                                    # §6.5 literal: the whole pool is excluded
            r["dropped"] = "pool contains positional isomers (§6.5)"
            r["n_positional_isomers"] = len(iso)
            rep[pool] = r
            continue
        n = 0
        for u in seen:
            truth[(pool, u)] = slot
            n += 1
        r.update(n_truth_peptides=n, n_positional_isomers_dropped=0)
        rep[pool] = r
    return {"truth": truth, "report": rep}


def truth_ptmkit(d: dict) -> dict:
    truth, rep = {}, {}
    for raw in d["raws"]:
        pool, run = pxd009449_pool(raw), run_name(raw)
        if pool not in PHOSPHO_KITS:
            rep[pool] = {"run": run, "phospho_truth": "none (not a phospho kit; E8 not evaluable on this pool)"}
            continue
        seqs = evidence_sequences([os.path.join(d["evidence_dir"], f"{pool}_evidence.txt")])
        ny = Counter(s.count("Y") for s in seqs)
        n = 0
        for s in seqs:
            if s.count("Y") == 1:
                truth[(run, s)] = s.index("Y")
                n += 1
        rep[pool] = {"run": run, "n_sequences": len(seqs), "by_Y_count": {str(k): v for k, v in sorted(ny.items())},
                     "n_truth_peptides": n}
    return {"truth": truth, "report": rep}


def truth_ferries() -> dict:
    import truth as TR  # scripts/m4/truth.py (dev only)
    return TR.pxd007058()


def truth_for(ds: dict) -> dict:
    out = {}
    for a, d in ds["accessions"].items():
        out[a] = {"marx": truth_marx, "ptmkit": truth_ptmkit}.get(d["kind"], lambda _d: truth_ferries())(d)
    return out
