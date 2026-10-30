#!/usr/bin/env python3
"""Render results/2026-10-27_M1e_baselines/tables/*.tsv and tables/main_tables.md from metrics.json."""
from __future__ import annotations

import csv
import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "2026-10-27_M1e_baselines")
TAB = os.path.join(RES, "tables")
PL = ["sage", "sage+mokapot", "comet+mokapot"]
SPEED_KEY = {"sage": "sage", "sage+mokapot": "sage+mokapot", "comet+mokapot": "comet+mokapot"}


def pct(x):
    return "n/a" if x is None else f"{100 * x:.2f}"


def write(name, header, rows):
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def verdict_label(ev):
    v = ev["verdict"]
    if v == "inconclusive" and ev.get("consistent_with_control"):
        return "inconclusive*"
    return v


def main():
    m = json.load(open(os.path.join(RES, "metrics.json")))
    ids, entp, entpr, spd = [], [], [], []
    md = []
    for ds, d in m["datasets"].items():
        for pl, x in d["pipelines"].items():
            prod, pp = x["ids_production"], x["protein_ids_production_picked_group"]
            for a, v in (prod or {}).items():
                ids.append([ds, pl, a, v["peptides"], v.get("peptides_vulcan_tdc", ""), v["psms_not_fdr_controlled"],
                            (pp or {}).get(a, "")])
            for tgt, level in ((entp, x), (entpr, x["protein"])):
                if not level:
                    continue
                P, C = level.get("arm_P"), level.get("arm_C")
                for a in [f"{q:g}" for q in m["alphas"]]:
                    if P:
                        ev = P["alpha"][a]["paired"] if "paired" in P["alpha"][a] else P["alpha"][a]
                        tgt.append([ds, pl, "P", a, ev["mean_n_t"], ev["mean_lb"], ev["mean_fdp"], ev["ci95"][0],
                                    ev["ci95"][1], ev["verdict"], ev["consistent_with_control"], ev["n_replicates"]])
                    if C:
                        ev = C["alpha"][a]
                        tgt.append([ds, pl, "C", a, ev["mean_n_t"], ev["mean_lb"], ev["mean_fdp"], ev["ci95"][0],
                                    ev["ci95"][1], ev["verdict"], ev["consistent_with_control"], 1])
        for tag, runs in d.get("_", {}).items():
            pass
        for tag, dbs in m["speed"].get(ds, {}).get("runs", {}).items():
            for db, r in dbs.items():
                for stage, key in (("search", "sage_search"), ("rescore", "sage_mokapot")):
                    if r.get(key):
                        spd.append([ds, "sage" if stage == "search" else "sage+mokapot", db, tag, stage] + speedrow(r[key]))
                for stage, key in (("search", "comet_search"), ("rescore", "comet_mokapot")):
                    if r.get(key):
                        spd.append([ds, "comet+mokapot", db, tag, stage] + speedrow(r[key]))
        conv = m["speed"].get(ds, {}).get("conversion_trfp_threads1_serial")
        if conv:
            spd.append([ds, "trfp", "-", "t1", "convert"] + speedrow(conv))
        # main table
        md.append(f"#### {ds} ({d['accession']}, {d['instrument']})\n")
        md.append("| tool | peptides@1% | proteins@1% | verdict P (pep) | verdict C (pep) | verdict P (prot) | verdict C (prot) | "
                  "wall s / core-s @32 thr | wall s / core-s @4 thr |")
        md.append("|---|---|---|---|---|---|---|---|---|")
        for pl in PL:
            x = d["pipelines"].get(pl)
            if not x:
                continue
            def vl(level, arm):
                a = (level or {}).get("arm_" + arm)
                if not a:
                    return "n/a"
                ev = a["alpha"]["0.01"]
                return verdict_label(ev["paired"] if "paired" in ev else ev)
            vp = {"P": vl(x, "P"), "C": vl(x, "C")}
            vq = {"P": vl(x["protein"], "P"), "C": vl(x["protein"], "C")}
            cells = [pl, str((x["ids_production"] or {}).get("0.01", {}).get("peptides", "n/a")),
                     str((x["protein_ids_production_picked_group"] or {}).get("0.01", "n/a")),
                     vp.get("P", "n/a"), vp.get("C", "n/a"), vq.get("P", "n/a"), vq.get("C", "n/a")]
            for tag in ("t32", "t4"):
                p = m["speed"].get(ds, {}).get("runs", {}).get(tag, {}).get("prod", {}).get("pipeline_search_plus_rescore", {}).get(pl)
                cells.append("n/a" if not p else f"{p['wall_s']:.0f} / {p['core_s']:.0f}")
            md.append("| " + " | ".join(cells) + " |")
        conv = m["speed"].get(ds, {}).get("conversion_trfp_threads1_serial")
        if conv:
            md.append(f"\nRaw -> mzML (TRFP, counted separately; both Sage and Comet need it): wall {conv['wall_s']:.0f} s, "
                      f"core-s {conv['core_s']:.0f}, max RSS {conv['max_rss_mb']:.0f} MB.\n")
    write("ids.tsv", ["dataset", "pipeline", "alpha", "peptides_tool_q", "peptides_vulcan_tdc", "psms_not_fdr_controlled",
                      "protein_groups_picked"], ids)
    hdr = ["dataset", "pipeline", "arm", "alpha", "mean_N_T", "mean_LB", "mean_FDPhat", "ci95_lo", "ci95_hi", "verdict",
           "consistent_with_control", "n_replicates"]
    write("entrapment_peptide.tsv", hdr, entp)
    write("entrapment_protein.tsv", hdr, entpr)
    write("speed.tsv", ["dataset", "pipeline", "db", "threads", "stage", "wall_s", "core_s", "max_rss_mb", "locked",
                        "n_runs", "start"], spd)
    with open(os.path.join(TAB, "main_tables.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(md))


def speedrow(r):
    return [r["wall_s"], r["core_s"], r["max_rss_mb"], r["locked"], r["n_runs"], r.get("start")]


if __name__ == "__main__":
    main()
