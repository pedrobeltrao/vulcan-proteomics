"""E8 truth tables per the definitions committed in results/2026-10-27_M4_ptm/NOTES.md §0 (commit 97daf98, before
any localization result). Dev files only; nothing under data/truth/sealed/ is read.

Each builder returns {"truth": {(pool_or_run_key, stripped_seq): 0-based true site position}, "report": {...}}.
"""
from __future__ import annotations

import csv
import glob
import os
import re
from collections import Counter, defaultdict

import common as C


# ------------------------------------------------------------------------------------------------ PXD000138
def infer_design(seq: str, pmin: int = 8, pmax: int = 30):
    """Unit length P (minimal number of variable positions; ties -> smallest P) and the variable positions."""
    best = None
    for P in range(pmin, pmax + 1):
        if len(seq) % P:
            continue
        units = [seq[i:i + P] for i in range(0, len(seq), P)]
        if len(units) < 3:
            continue
        var = [j for j in range(P) if len({u[j] for u in units}) > 1]
        if not var:
            continue
        if best is None or len(var) < len(best[1]):
            best = (P, var, units)
    return best


def pxd000138():
    truth, rep = {}, {}
    for f in sorted(glob.glob(os.path.join(C.ROOT, "data/truth/PXD000138/pool_*.fasta"))):
        well = os.path.basename(f)[5:-6]
        seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
        P, var, units = infer_design(seq)
        sets = {j: sorted({u[j] for u in units}) for j in var}
        slots = [j for j, s in sets.items() if set(s) == set("STY")]
        pool = [k for k, v in C.POOL_WELL.items() if v == well][0]
        r = {"well": well, "unit_len": P, "n_units": len(units), "variable_positions": {str(j): "".join(s) for j, s
             in sets.items()}, "acceptor_slots": slots}
        if len(slots) != 1:
            r["dropped"] = "not exactly one {S,T,Y} slot"
            rep[pool] = r
            continue
        slot = slots[0]
        seen = defaultdict(set)
        for u in units:
            seen[u].add(slot)
        iso = [u for u, s in seen.items() if len(s) > 1]
        for u, s in seen.items():
            if len(s) == 1:
                truth[(pool, u)] = slot
        r.update(n_truth_peptides=sum(1 for k in truth if k[0] == pool), n_positional_isomers_dropped=len(iso))
        rep[pool] = r
    return {"truth": truth, "report": rep}


# ------------------------------------------------------------------------------------------------ PXD009449
def pxd009449():
    seqs = C.ymod_sequences()
    run = C.run_name(C.D2_RAW["PXD009449"][0])
    truth = {}
    ny = Counter()
    for s in seqs:
        n = s.count("Y")
        ny[n] += 1
        if n == 1:
            truth[(run, s)] = s.index("Y")
    return {"truth": truth, "report": {"n_sequences": len(seqs), "by_Y_count": dict(ny),
                                       "n_truth_peptides": len(truth)}}


# ------------------------------------------------------------------------------------------------ PXD007058
_PROB = re.compile(r"([A-Z])\(([\d.]+)\)")


def _site_probs(s: str):
    """'AQS(0.99)FPDNR' -> (stripped, {pos: prob})."""
    out, seq, i = {}, [], 0
    for m in re.finditer(r"([A-Z])(\(([\d.]+)\))?", s):
        seq.append(m.group(1))
        if m.group(3) is not None:
            out[len(seq) - 1] = float(m.group(3))
    return "".join(seq), out


def pxd007058():
    truth, rep = {}, {}
    tables = C.ferries_tables()
    for pool, t in tables.items():
        run = f"SF_200217_pPeptideLibrary_pool{pool}_HCDOT_rep1"
        per = defaultdict(list)          # seq -> [(argmax pos, max prob, probs)]
        prot = {}
        n_rows = 0
        for r in t["msms"]:
            if r.get("Reverse", "") == "+" or "CON__" in r.get("Proteins", ""):
                continue
            if r.get("Phospho (STY)", "0") != "1":
                continue
            ps = r.get("Phospho (STY) Probabilities", "")
            if not ps:
                continue
            seq, pr = _site_probs(ps)
            if not pr:
                continue
            n_rows += 1
            pos = max(pr, key=lambda k: (pr[k], -k))
            per[seq].append((pos, pr[pos], pr))
            prot[seq] = r.get("Proteins", "")
        n_elig = n_seal = 0
        drop = Counter()
        for seq, rows in per.items():
            wells = set(re.findall(r"_([A-H]\d+)(?:;|$)", prot[seq]))
            if wells and wells <= C.SEALED_T2_WELLS:
                n_seal += 1
                continue
            if not any(mp >= 0.99 for _, mp, _ in rows):
                drop["a_no_psm_ge_0.99"] += 1
                continue
            conf = [pos for pos, mp, _ in rows if mp >= 0.75]
            top = Counter(conf).most_common(1)[0]
            if top[1] < 0.9 * len(conf):
                drop["b_disagreement"] += 1
                continue
            hi = {pos for pos, mp, _ in rows if mp >= 0.99}
            if len(hi) > 1:
                drop["c_positional_isomer"] += 1
                continue
            truth[(run, seq)] = top[0]
            n_elig += 1
        rep[pool] = {"run": run, "n_mono_phospho_msms_rows": n_rows, "n_sequences": len(per),
                     "n_truth_peptides": n_elig, "excluded_sealed_wells": n_seal, "dropped": dict(drop)}
    return {"truth": truth, "report": rep}


def all_truth():
    return {"PXD000138": pxd000138(), "PXD009449": pxd009449(), "PXD007058": pxd007058()}


if __name__ == "__main__":
    import json
    out = all_truth()
    print(json.dumps({k: v["report"] for k, v in out.items()}, indent=1))
