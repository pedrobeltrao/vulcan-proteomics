#!/usr/bin/env python3
"""M1e analysis: IDs (peptide + protein group), entrapment FDP (LB / combined / paired), §6.3 verdicts, speed.

    python scripts/baselines/analyze.py            # -> results/2026-10-27_M1e_baselines/metrics.json + tables/*.tsv

Estimators and the decision rule come from vulcan.stats (imported, not re-implemented):
  entrapment.{paired_counts, combined_counts, protein_paired_counts, evaluate_arm, combine_arms}, tdc.tdc_qvalues,
  protein.picked_group_fdr.

Units
  * peptide = stripped sequence, I->L collapsed (all peptidoforms pooled); accepted at alpha if the tool's own
    peptide-level q <= alpha (Sage: LDA peptide_q; mokapot: peptide q-value); score = best PSM score.
    Extra column `peptides_vulcan_tdc`: vulcan.stats.tdc TDC over the same per-peptide best scores (production only).
  * protein group = vulcan.stats.protein.picked_group_fdr (subset grouping, rescue on, seed 11) on the tool's per-peptide
    best scores of target AND decoy peptides; counted at picked-group q <= alpha (non-obsolete target groups).
    Augmented searches: entrapment proteins are built from the peptide-level pairing (arm P: entrapment protein
    "E|<acc>" = the paired entrapment peptides of acc's peptides; shared peptides are therefore swapped consistently;
    arm C: foreign proteins "F|<acc>" over the filtered foreign peptides).  Decoy proteins "rev_<protein>".
    A group with >= 1 original target protein is a target (§6.3).  Arm-C protein-level r = (# foreign proteins with
    >= 1 retained in-space peptide) / (# target proteins with >= 1 in-space peptide).
  * PSM counts are "not FDR-controlled" (DESIGN §6.8) and come from the production search.
Pipelines: sage (built-in LDA), sage+mokapot (Sage PIN), comet+mokapot (Comet PIN).
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "entrapment"))
from vulcan.stats import entrapment as ent  # noqa: E402
from vulcan.stats import protein as prot  # noqa: E402
from vulcan.stats import tdc  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from run_searches import DATASETS, SEEDS, ENT, RES, TIMINGS, WORK  # noqa: E402
from peptides import DECOY_SEEDS, STD_AA, accession, digest, il, in_space, read_fasta  # noqa: E402

CACHE = os.path.join(WORK, "cache")          # the full-budget (affinity threads) pass defines the metrics
ALPHAS = [0.001, 0.005, 0.01, 0.05]
PIPELINES = ["sage", "sage+mokapot", "comet+mokapot"]
MODEL_CLASS = "static_low_capacity"   # LDA / linear SVM over hand-crafted features: required arms P and C
N_BOOT = 2000
PICK_SEED = 11
FA = os.path.join(ROOT, "data", "fasta")
FOREIGN_FASTA = {"human_contam": "arath_UP000006548.fasta", "pfu_contam": "human_sp_canonical.fasta"}

_mod = re.compile(r"\[[^\]]*\]|\([^)]*\)")


def strip(p: str) -> str:
    if p.count(".") >= 2 and re.match(r"^[A-Za-z\-n]?\.", p):
        p = p[p.index(".") + 1: p.rindex(".")]
    p = _mod.sub("", p)
    return re.sub(r"[^A-Z]", "", p.replace("n", "")).replace("I", "L")


def read_tsv(path):
    with open(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def read_ragged(path, last="Proteins"):
    """PIN / mokapot-style TSV whose last named column (`Proteins`) may spill into extra tab-separated fields."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        n = len(hdr)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = dict(zip(hdr[:n - 1], f[:n - 1]))
            d[hdr[n - 1]] = "\t".join(f[n - 1:]) if hdr[n - 1] == last else f[n - 1]
            yield d


def split_prots(s: str):
    return {x.strip() for x in re.split(r"[;,\t]", s) if x.strip()}


# ------------------------------------------------------------------------- tool outputs -> peptide tables

class Res:
    """Target-side peptides {seq: [q, score]}, decoy-side peptides {seq: score}, protein lists {(side, seq): set}."""

    def __init__(self):
        self.psm_q = []
        self.peps = {}
        self.decs = {}
        self.pep_prots = defaultdict(set)

    def add_t(self, k, q, s, prots=()):
        cur = self.peps.get(k)
        self.peps[k] = [min(q, cur[0]), max(s, cur[1])] if cur else [q, s]
        self.pep_prots[("T", k)] |= set(prots)

    def add_d(self, k, s, prots=()):
        self.decs[k] = max(s, self.decs.get(k, -np.inf))
        self.pep_prots[("D", k)] |= set(prots)


def load(pipeline, ds, db, cache=CACHE):
    r = Res()
    if pipeline == "sage":
        f = os.path.join(cache, "sage", ds, db, "results.sage.tsv")
        if not os.path.exists(f):
            return None
        for row in read_tsv(f):
            s = float(row["sage_discriminant_score"])
            k = strip(row["peptide"])
            if row["label"] == "1":
                r.psm_q.append(float(row["spectrum_q"]))
                r.add_t(k, float(row["peptide_q"]), s, split_prots(row["proteins"]))
            else:
                r.add_d(k, s, split_prots(row["proteins"]))
        r.psm_q = np.array(r.psm_q)
        return r
    tool = pipeline.split("+")[0]
    d = os.path.join(cache, f"{tool}_mokapot", ds, db)
    fpep = os.path.join(d, "mokapot.peptides.txt")
    if not os.path.exists(fpep):
        return None
    r.psm_q = np.array([float(x["mokapot q-value"]) for x in read_ragged(os.path.join(d, "mokapot.psms.txt"))])
    for x in read_ragged(fpep):
        r.add_t(strip(x["Peptide"]), float(x["mokapot q-value"]), float(x["mokapot score"]), split_prots(x["Proteins"]))
    fd = os.path.join(d, "mokapot.decoy.peptides.txt")      # mokapot 0.10.0 spelling
    if os.path.exists(fd):
        for x in read_ragged(fd):
            r.add_d(strip(x["Peptide"]), float(x["mokapot score"]), split_prots(x["Proteins"]))
    return r


_dbmaps = {}


def peptide_map(ds, db):
    """kind {seq: T/E}, pair key {seq: target seq} for forward peptides of a peptide-level DB; plus decoy_of."""
    fam = DATASETS[ds]["family"]
    if (fam, db) not in _dbmaps:
        kind, key, dec_of = {}, {}, {}
        for r in read_tsv(os.path.join(ENT, "db", f"{fam}.{db}.peptides.tsv")):
            if r["kind"] in ("T", "E"):
                kind[r["seq"]] = r["kind"]
                key[r["seq"]] = r["pair_key"]
            else:
                dec_of[r["seq"]] = (r["kind"], r["decoy_of"])
        _dbmaps[(fam, db)] = (kind, key, dec_of)
    return _dbmaps[(fam, db)]


def db_meta(ds, db):
    fam = DATASETS[ds]["family"]
    with open(os.path.join(ENT, "db", f"{fam}.{db}.json")) as fh:
        return json.load(fh)


# ------------------------------------------------------------------------------------------- peptide level

def id_counts(res, with_tdc=False):
    pq = np.array([v[0] for v in res.peps.values()])
    out = {f"{a:g}": {"peptides": int((pq <= a).sum()), "psms_not_fdr_controlled": int((res.psm_q <= a).sum())}
           for a in ALPHAS}
    if with_tdc and res.decs:
        sc = np.array([v[1] for v in res.peps.values()] + list(res.decs.values()))
        isd = np.r_[np.zeros(len(res.peps), bool), np.ones(len(res.decs), bool)]
        q = tdc.tdc_qvalues(sc, isd)
        for a in ALPHAS:
            out[f"{a:g}"]["peptides_vulcan_tdc"] = tdc.count_discoveries(q, isd, a)
    return out


def entrapment_rows(res, kind, key):
    seqs = [s for s in res.peps if s in kind]
    is_e = np.array([kind[s] == "E" for s in seqs])
    q = np.array([res.peps[s][0] for s in seqs])
    score = np.array([res.peps[s][1] for s in seqs])
    pk = np.array([key[s] or s for s in seqs], dtype=object)
    return seqs, is_e, q, score, pk, len(res.peps) - len(seqs)


def consistent(ev):
    """DESIGN §6.3: 'consistent with control' (weaker, labelled) if the lower CI of the FDP-hat is <= alpha."""
    return bool(ev["ci95"][0] <= ev["alpha"])


# ------------------------------------------------------------------------------------------- protein level

_tmaps = {}


def target_prot_map(fam):
    """{I/L peptide: set(acc)} for the family's target FASTA (human/pfu + contaminants), digested like the DBs."""
    if fam not in _tmaps:
        m = defaultdict(set)
        for h, s in read_fasta(os.path.join(ENT, "inputs", f"{fam}.fasta")):
            a = accession(h)
            for p in digest(s):
                q = il(p)
                if set(q) <= STD_AA:
                    m[q].add(a)
        _tmaps[fam] = m
    return _tmaps[fam]


_fmaps = {}


def foreign_prot_map(fam):
    if fam not in _fmaps:
        m = defaultdict(set)
        for h, s in read_fasta(os.path.join(FA, FOREIGN_FASTA[fam])):
            a = "F|" + accession(h)
            for p in digest(s):
                q = il(p)
                if set(q) <= STD_AA:
                    m[q].add(a)
        _fmaps[fam] = m
    return _fmaps[fam]


def pair_of_prefix(p):
    for pre in ("rev_", "DECOY_"):
        if p.startswith(pre):
            return p[len(pre):]
    return p


def run_picked(pep2prot, pep_score, is_dec):
    pair_of = {}
    for ps in pep2prot.values():
        for p in ps:
            pair_of[p] = pair_of_prefix(p)
    is_decoy_prot = {p: (pair_of[p] != p) for p in pair_of}
    return prot.picked_group_fdr(pep2prot, pep_score, pair_of, is_decoy_prot, seed=PICK_SEED)


def protein_ids_production(res):
    """Picked protein-group IDs at q<=alpha using the tool's own protein lists (target + decoy peptides)."""
    pep2prot, score = {}, {}
    for (side, k), ps in res.pep_prots.items():
        if not ps:
            continue
        pep2prot[(side, k)] = sorted(ps)
        score[(side, k)] = res.peps[k][1] if side == "T" else res.decs[k]
    r = run_picked(pep2prot, score, None)
    keep = (~r.is_decoy) & (~r.is_obsolete)
    return {f"{a:g}": int(((r.q <= a) & keep).sum()) for a in ALPHAS}


def augmented_pep2prot(ds, db, res):
    """pep2prot / scores for an augmented (peptide-level) search, plus entrapment-protein metadata."""
    fam = DATASETS[ds]["family"]
    kind, key, dec_of = peptide_map(ds, db)
    tmap = target_prot_map(fam)
    fmap = foreign_prot_map(fam) if db == "C" else None
    ent_prot, partner = {}, {}

    def fwd_prots(seq):
        if kind.get(seq) == "T":
            return sorted(tmap[seq])
        if kind.get(seq) == "E":
            if fmap is not None:
                return sorted(fmap[seq])
            t = key[seq]
            if not t:
                return []
            out = ["E|" + a for a in sorted(tmap[t])]
            for e in out:
                partner[e] = e[2:]
            return out
        return []
    pep2prot, score = {}, {}
    for k, (q, s) in res.peps.items():
        ps = fwd_prots(k)
        if ps:
            pep2prot[("T", k)], score[("T", k)] = ps, s
    for k, s in res.decs.items():
        if k not in dec_of:
            continue
        _, fseq = dec_of[k]
        ps = fwd_prots(fseq)
        if ps:
            pep2prot[("D", k)], score[("D", k)] = ["rev_" + p for p in ps], s
    return pep2prot, score, partner


def protein_r(fam, kind_map):
    """Protein-level r for arm C: foreign proteins with >=1 retained peptide / target proteins with >=1 in-space peptide."""
    tmap, fmap = target_prot_map(fam), foreign_prot_map(fam)
    tp = {a for p, accs in tmap.items() if in_space(p) for a in accs}
    fp = {a for p, accs in fmap.items() if kind_map.get(p) == "E" for a in accs}
    return len(fp) / max(len(tp), 1), len(tp), len(fp)


def eval_protein(ds, pipeline):
    """Protein-group level: IDs, arm P (3 entrapment seeds) and arm C evaluation."""
    fam = DATASETS[ds]["family"]
    out = {}
    reps = defaultdict(list)
    seeds_done = []
    for seed in SEEDS:
        db = f"P{seed}"
        res = load(pipeline, ds, db)
        if res is None:
            continue
        p2p, sc, partner = augmented_pep2prot(ds, db, res)
        r = run_picked(p2p, sc, None)
        seeds_done.append(seed)
        isE = defaultdict(bool)
        for g in r.groups:
            for p in g:
                isE[p] = p.startswith("E|") or p.startswith("rev_E|")
        keep = (~r.is_decoy) & (~r.is_obsolete)
        for a in ALPHAS:
            sel = np.flatnonzero((r.q <= a) & keep)
            groups = [r.groups[i] for i in sel]
            c = ent.protein_paired_counts(groups, r.score[sel], np.ones(len(sel), bool), isE, partner)
            reps[a].append(ent.Replicate(counts=c, entrap_seed=seed))
    if seeds_done:
        out["arm_P"] = {"seeds": seeds_done, "alpha": {}}
        for a in ALPHAS:
            ev = ent.evaluate_arm(reps[a], a, estimator="paired", r=1.0, n_boot=N_BOOT, seed=PICK_SEED)
            ev["consistent_with_control"] = consistent(ev)
            out["arm_P"]["alpha"][f"{a:g}"] = ev
    res = load(pipeline, ds, "C")
    if res is not None:
        kind, key, _ = peptide_map(ds, "C")
        p2p, sc, partner = augmented_pep2prot(ds, "C", res)
        r_count, n_tp, n_fp = protein_r(fam, kind)
        rr = db_meta(ds, "C")["r"]    # spec 6.3: r on the effective peptide search space; protein-count ratio = sensitivity only
        r = run_picked(p2p, sc, None)
        keep = (~r.is_decoy) & (~r.is_obsolete)
        out["arm_C"] = {"r_used_peptide_space": rr, "r_protein_count_ratio_sensitivity": r_count, "n_target_proteins_in_space": n_tp, "n_foreign_proteins": n_fp, "alpha": {}}
        for a in ALPHAS:
            sel = np.flatnonzero((r.q <= a) & keep)
            is_t = ent.protein_group_is_target([r.groups[i] for i in sel], defaultdict(bool, {
                p: p.startswith("F|") or p.startswith("rev_F|") for i in sel for p in r.groups[i]}))
            ids = np.array([str(i) for i in sel], dtype=object)
            c = ent.combined_counts(ids, ~is_t, np.ones(len(sel), bool))
            ev = ent.evaluate_arm([ent.Replicate(counts=c)], a, estimator="combined", r=rr, n_boot=N_BOOT, seed=PICK_SEED)
            ev["N_E"] = int(c.e.sum())
            ev["consistent_with_control"] = consistent(ev)
            alt = ent.evaluate_arm([ent.Replicate(counts=c)], a, estimator="combined", r=r_count, n_boot=N_BOOT, seed=PICK_SEED)
            ev["sensitivity_r_protein_count"] = {"r": r_count, "mean_fdp": alt["mean_fdp"], "ci95": alt["ci95"],
                                                 "verdict": alt["verdict"]}
            out["arm_C"]["alpha"][f"{a:g}"] = ev
    return out


# ------------------------------------------------------------------------------------------- entrapment (peptide)

def eval_P(pipeline, ds):
    """Arm P: replicates = 3 entrapment seeds (pseudo-reverse decoy only in dev)."""
    reps_by_alpha = defaultdict(list)
    per_seed = {}
    for seed in SEEDS:
        db = f"P{seed}"
        res = load(pipeline, ds, db)
        if res is None:
            continue
        kind, key, _ = peptide_map(ds, db)
        seqs, is_e, q, score, pk, unknown = entrapment_rows(res, kind, key)
        per_seed[seed] = {"unknown_peptides": unknown, "ids": id_counts(res)}
        for a in ALPHAS:
            acc = q <= a
            c = ent.paired_counts(pk[acc].astype(str), is_e[acc], score[acc], acc[acc])
            reps_by_alpha[a].append(ent.Replicate(counts=c, entrap_seed=seed, decoy_seed=None))
            n_t, n_e, n_eo, n_egt = c.totals()
            per_seed[seed][f"{a:g}"] = {"N_T": int(n_t), "N_E": int(n_e), "N_E_ge_s_gt_T": int(n_eo),
                                        "N_E_gt_T_ge_s": int(n_egt),
                                        "fdp_lb": float(ent.lower_bound(n_t, n_e)),
                                        "fdp_combined_r1": float(ent.combined(n_t, n_e, 1.0)),
                                        "fdp_paired": float(ent.paired(n_t, n_e, n_eo, n_egt))}
    if not per_seed:
        return None
    out = {"seeds": sorted(per_seed), "per_seed": per_seed, "alpha": {}}
    for a in ALPHAS:
        reps = reps_by_alpha[a]
        p = ent.evaluate_arm(reps, a, estimator="paired", r=1.0, n_boot=N_BOOT, seed=PICK_SEED)
        p["consistent_with_control"] = consistent(p)
        c = ent.evaluate_arm(reps, a, estimator="combined", r=1.0, n_boot=N_BOOT, seed=PICK_SEED)
        out["alpha"][f"{a:g}"] = {"paired": p, "combined_r1": c}
    return out


def eval_C(pipeline, ds):
    db = "C"
    res = load(pipeline, ds, db)
    if res is None:
        return None
    meta = db_meta(ds, db)
    r = meta["r"]
    kind, key, _ = peptide_map(ds, db)
    seqs, is_e, q, score, pk, unknown = entrapment_rows(res, kind, key)
    out = {"r": r, "r_per_mass_decile": meta["r_per_mass_decile"], "unknown_peptides": unknown,
           "ids": id_counts(res), "alpha": {}}
    for a in ALPHAS:
        acc = q <= a
        c = ent.combined_counts(np.array(seqs, dtype=object)[acc].astype(str), is_e[acc], acc[acc])
        ev = ent.evaluate_arm([ent.Replicate(counts=c)], a, estimator="combined", r=r, n_boot=N_BOOT, seed=PICK_SEED)
        ev["N_E"] = int(c.e.sum())
        ev["consistent_with_control"] = consistent(ev)
        out["alpha"][f"{a:g}"] = ev
    return out


def eval_S(pipeline, ds):
    """Decoy-variability arm (DESIGN 6.2): production-style search against T + 5 external shuffled-decoy databases
    (DECOY_SEEDS 1001-1005; Sage generate_decoys=false / Comet decoy_search=0). IDs only (no entrapment)."""
    per = {}
    for seed in DECOY_SEEDS:
        db = f"S{seed}"
        res = load(pipeline, ds, db)
        if res is None:
            continue
        p2p, sc, _ = augmented_pep2prot(ds, db, res)
        r = run_picked(p2p, sc, None)
        keep = (~r.is_decoy) & (~r.is_obsolete)
        ids = id_counts(res, with_tdc=True)
        per[seed] = {f"{a:g}": {"peptides": ids[f"{a:g}"]["peptides"],
                                "peptides_vulcan_tdc": ids[f"{a:g}"].get("peptides_vulcan_tdc"),
                                "protein_groups_picked": int(((r.q <= a) & keep).sum())} for a in ALPHAS}
    if not per:
        return None
    out = {"seeds": sorted(per), "per_seed": per, "mean_sd": {}}
    for a in ALPHAS:
        k = f"{a:g}"
        for key in ("peptides", "protein_groups_picked"):
            v = np.array([per[sd][k][key] for sd in per], float)
            out["mean_sd"].setdefault(k, {})[key] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0}
    return out


def verdicts_for(P, C):
    v = {}
    for a in ALPHAS:
        k = f"{a:g}"
        arms = {}
        if P:
            arms["P"] = P["alpha"][k]["paired"]["verdict"] if "paired" in P["alpha"][k] else P["alpha"][k]["verdict"]
        if C:
            arms["C"] = C["alpha"][k]["verdict"]
        v[k] = {"arms": arms, "overall": ent.combine_arms(arms, MODEL_CLASS) if arms else None,
                "gated": a in ent.GATED_ALPHAS}
    return v


# ------------------------------------------------------------------------------------------- speed

def n_ms2(ds):
    tot = 0
    for f in DATASETS[ds]["files"]:
        m = json.load(open(os.path.join(WORK, "meta", f"{f}-metadata.json")))
        tot += int(next(e["value"] for e in m["MsData"] if e["name"] == "Number of MS2 spectra"))
    return tot


def speed():
    recs = [json.loads(x) for x in open(TIMINGS)] if os.path.exists(TIMINGS) else []
    last = {}
    for r in recs:     # keep the last successful record per label
        if r.get("rc") == 0:
            last[r["label"]] = r
    out = {}
    for ds in DATASETS:
        try:
            nms2 = n_ms2(ds)
        except (FileNotFoundError, StopIteration):
            continue
        files = DATASETS[ds]["files"]

        def agg(labels):
            rs = [last[x] for x in labels if x in last]
            if len(rs) != len(labels):
                return None
            return {"wall_s": round(sum(r["wall_s"] for r in rs), 2), "core_s": round(sum(r["core_s"] for r in rs), 2),
                    "max_rss_mb": round(max(r["max_rss_mb"] for r in rs), 1), "locked": all(r["locked"] for r in rs),
                    "n_runs": len(rs), "threads": rs[0].get("threads_env"), "start": rs[0].get("start"),
                    "hw": rs[0].get("hw")}
        conv = agg([f"trfp_mzml:{ds}:{f}" for f in files])
        if conv:
            conv.pop("hw", None)
        dsd = {"n_ms2": nms2, "conversion_trfp_threads1_serial": conv, "runs": {}}
        for T in (4, 32):
            tag = f"@t{T}"
            for db in ("prod", "target", "P2001", "P2002", "P2003", "C"):
                d = {"sage_search": agg([f"sage:{ds}:{db}{tag}"]),
                     "sage_mokapot": agg([f"mokapot:sage:{ds}:{db}{tag}"]),
                     "comet_search": agg([f"comet:{ds}:{db}:{f}{tag}" for f in files]),
                     "comet_mokapot": agg([f"mokapot:comet:{ds}:{db}{tag}"])}
                if not any(d.values()):
                    continue
                for k in ("sage_search", "comet_search"):
                    if d[k]:
                        d[k]["spectra_per_s"] = round(nms2 / d[k]["wall_s"], 1)
                        d[k]["core_s_per_1e5_ms2"] = round(d[k]["core_s"] / nms2 * 1e5, 1)
                pip = {}
                for name, parts in {"sage": ["sage_search"], "sage+mokapot": ["sage_search", "sage_mokapot"],
                                    "comet+mokapot": ["comet_search", "comet_mokapot"]}.items():
                    ps = [d[p] for p in parts]
                    if all(ps):
                        pip[name] = {"wall_s": round(sum(p["wall_s"] for p in ps), 1),
                                     "core_s": round(sum(p["core_s"] for p in ps), 1),
                                     "max_rss_mb": max(p["max_rss_mb"] for p in ps),
                                     "all_locked": all(p["locked"] for p in ps),
                                     "excludes_conversion": True}
                d["pipeline_search_plus_rescore"] = pip
                dsd["runs"].setdefault(f"t{T}", {})[db] = d
        out[ds] = dsd
    return out


def main():
    metrics = {"experiment": "M1e baselines (dev data only)", "date": "2026-10-27", "alphas": ALPHAS,
               "pipelines": PIPELINES, "model_class_for_required_arms": MODEL_CLASS, "n_boot": N_BOOT,
               "pick_seed": PICK_SEED, "entrapment_seeds": SEEDS,
               "estimators": "vulcan.stats.{entrapment,tdc,protein}", "datasets": {}}
    for ds in DATASETS:
        dsd = {"accession": DATASETS[ds]["accession"], "instrument": DATASETS[ds]["instrument"],
               "files": DATASETS[ds]["files"], "pipelines": {}}
        for pl in PIPELINES:
            prod = load(pl, ds, "prod")
            tgt = load(pl, ds, "target")
            P = eval_P(pl, ds)
            C = eval_C(pl, ds)
            if prod is None and P is None and C is None:   # (S arm alone does not create a pipeline entry)
                continue
            PR = eval_protein(ds, pl) if (P or C) else {}
            dsd["pipelines"][pl] = {
                "ids_production": id_counts(prod, with_tdc=True) if prod else None,
                "protein_ids_production_picked_group": protein_ids_production(prod) if prod else None,
                "ids_peptide_level_target_only": id_counts(tgt) if tgt else None,
                "arm_P": P, "arm_C": C,
                "decoy_variability_S_external_shuffled_1001_1005": eval_S(pl, ds),
                "verdicts_peptide": verdicts_for(P, C),
                "protein": PR,
                "verdicts_protein": verdicts_for(PR.get("arm_P"), PR.get("arm_C")) if PR else None}
            print(ds, pl, "done", flush=True)
        for db in ("prod", "target", "P2001", "P2002", "P2003", "C"):
            p = os.path.join(ENT, "db", f"{DATASETS[ds]['family']}.{db}.json")
            if os.path.exists(p):
                dsd.setdefault("databases", {})[db] = {k: v for k, v in json.load(open(p)).items()
                                                       if not k.endswith("_examples")}
        metrics["datasets"][ds] = dsd
    metrics["speed"] = speed()
    with open(os.path.join(RES, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("wrote", os.path.join(RES, "metrics.json"))


if __name__ == "__main__":
    main()
