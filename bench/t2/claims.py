#!/usr/bin/env python3
"""Deterministic T2 claim table from a driver metrics.json (rules and labels: config/t2.yaml).

    python bench/t2/claims.py <out>/metrics.json [--out DIR]

Every row is printed whatever its outcome.  The completeness rule (A5.5.4 pattern) is checked here: a missing
component, decoy, tool, run or bootstrap makes the look "incomplete" (the driver then exits non-zero), and a missing or
failed comparator makes every comparison row "not evaluable -> not beaten".  Output depends only on metrics.json and
config/t2.yaml (no clock, no environment).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GATE_COMPONENTS = ("PXD000138", "PXD009449")       # PXD007058 (dev proxy truth) never gates (M4 NOTES §0)
EXPECTED = {"T2": {"PXD000138": 4, "PXD009449": 4}, "D2T": {"PXD000138": 4, "PXD009449": 4},
            "D2": {"PXD000138": 4, "PXD009449": 1, "PXD007058": 4}}
ALPHAS = ("0.01", "0.05")


def f(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def pct(x):
    return "NA" if x is None else f"{100 * x:.2f}%"


def bnd(d, path, kind="pct"):
    """' [run-to-run min–max over N runs]' from the component's t2_band (PI D35); '' if absent."""
    b = d.get("t2_band", {}).get(path)
    if not b or "min" not in b:
        return " [band NA]"
    fm = pct if kind == "pct" else (lambda x: f(x) if isinstance(x, float) else str(x))
    return f" [{fm(b['min'])}–{fm(b['max'])} over {b['n']} runs]"


def vband(d, path):
    b = d.get("t2_band", {}).get(path)
    if not b or "counts" not in b:
        return ""
    if len(b["counts"]) == 1:
        return f" (same verdict in all {b['n']} runs)"
    return " (verdict NOT stable across runs: " + ", ".join(f"{k} {v}/{b['n']}" for k, v in sorted(b["counts"].items())) + ")"


def completeness(m) -> list[str]:
    miss = []
    dsn = m["meta"]["dataset"]
    for acc, nrun in EXPECTED[dsn].items():
        d = m["datasets"].get(acc)
        if d is None:
            miss.append(f"component {acc}")
            continue
        if len(d.get("runs", [])) != nrun:
            miss.append(f"{acc} runs")
        for tool in ("ours", "pyascore"):
            for dec in ("A", "L"):
                k = f"{tool}_{dec}"
                if k not in d or "error" in d[k] or "eval" not in d[k]:
                    miss.append(f"{acc} {k}")
                elif d["n_truth_eligible_psms"] and k not in d.get("t2", {}):
                    miss.append(f"{acc} t2 {k}")
        if d["n_truth_eligible_psms"]:
            for dec in ("A", "L"):
                for a in ALPHAS:
                    c = d.get("t2", {}).get(f"ours_{dec}", {}).get(a, {})
                    if "verdict" not in c:
                        miss.append(f"{acc} verdict {dec} {a}")
                    elif c.get("n_accepted_eligible") and "n_boot_nan" not in c:
                        miss.append(f"{acc} bootstrap {dec} {a}")
                    elif f"ours_{dec}/{a}/true_flr" not in d.get("t2_band", {}):
                        miss.append(f"{acc} repeat band {dec} {a}")
    n_runs = sum(EXPECTED[dsn].values())
    if len(m.get("runtime_per_run", {})) != n_runs:
        miss.append("runtime_per_run")
    return miss


def build(m, cfg):
    L = cfg["labels"]
    rows = []
    dsn = m["meta"]["dataset"]
    miss = completeness(m)
    psm_label = ["L-PSM", "L-D-A5", "L-TRUTH", "L-TUNE"]

    # ---------------------------------------------------------------- E8 cells and gate
    gate_cells = []
    for acc in sorted(m["datasets"]):
        d = m["datasets"][acc]
        has_truth = bool(d.get("n_truth_eligible_psms"))
        for dec in ("A", "L"):
            for a in ALPHAS:
                c = d.get("t2", {}).get(f"ours_{dec}", {}).get(a, {})
                v = c.get("verdict", "not evaluable (no truth)" if not has_truth else "missing")
                if dec == "A" and acc in GATE_COMPONENTS and has_truth:
                    gate_cells.append(v)
                rows.append({"id": f"E8-{acc}-{dec}-{a}", "kind": "E8 cell" if dec == "A" else "E8 cell (Leu, secondary; reported, not a gate input)",
                             "component": acc, "decoy": dec, "alpha": a,
                             "statistic": f"true FLR {pct(c.get('true_flr'))}{bnd(d, f'ours_{dec}/{a}/true_flr')} "
                                          f"(n {f(c.get('n_accepted_eligible'))}{bnd(d, f'ours_{dec}/{a}/n_accepted_eligible', 'n')} accepted eligible STY sites, "
                                          f"{f(c.get('k_wrong'))} wrong); one-sided 95% UCB {pct(c.get('ucb'))}{bnd(d, f'ours_{dec}/{a}/ucb')} "
                                          f"(CP {pct(c.get('ucb_cp'))}, cluster bootstrap {pct(c.get('ucb_boot'))})",
                             "verdict": v + (vband(d, f"ours_{dec}/{a}/verdict") if has_truth else ""), "labels": psm_label})
    if not gate_cells:
        gate = "not evaluable (no gate component with truth)"
    elif any(v == "FAILED" for v in gate_cells):
        gate = "FAILED"
    elif all(v == "PASS" for v in gate_cells):
        gate = "PASS"
    else:
        gate = "not passed (inconclusive or not evaluable cells)"
    no_truth = [acc for acc in GATE_COMPONENTS if acc in m["datasets"] and not m["datasets"][acc].get("n_truth_eligible_psms")]
    rows.insert(0, {"id": "G-E8", "kind": "E8 gate (DESIGN §6.5), Ala, mono-phospho, PSM level (D8)",
                    "statistic": f"cells {len(gate_cells)}: " + "; ".join(sorted(set(gate_cells))) +
                                 (f"; components without phospho truth (excluded): {', '.join(no_truth)}" if no_truth else ""),
                    "verdict": gate, "labels": psm_label})

    # ---------------------------------------------------------------- site recall vs pyAscore
    for acc in sorted(m["datasets"]):
        d = m["datasets"][acc]
        t2 = d.get("t2", {})
        for dec in ("A", "L"):
            dd = t2.get(f"diff_{dec}", {})
            for kind in ("true", "estimated"):
                for a in ALPHAS:
                    c = dd.get(f"{kind}_{a}")
                    rid = f"R-{kind}-{acc}-{dec}-{a}"
                    title = (f"site recall at matched true FLR {a}" if kind == "true" else f"site recall at estimated FLR̂ {a}")
                    if dec == "L":
                        title += " (Leu, secondary)"
                    if not d.get("n_truth_eligible_psms"):
                        v, stat = "not evaluable (no truth)", "NA"
                    elif c is None:
                        v, stat = "not evaluable (comparator missing) -> not beaten", "NA"
                    else:
                        v = c["verdict"]
                        if v == "not evaluable":
                            v = "not evaluable -> not beaten"
                        pth = f"diff_{dec}/{kind}_{a}"
                        stat = (f"ours {f(c['ours'])}{bnd(d, pth + '/ours', 'n')} vs pyAscore {f(c['pyascore'])}; "
                                f"diff {f(c['diff'])}{bnd(d, pth + '/diff', 'n')} "
                                f"[{f(c['ci95'][0])}, {f(c['ci95'][1])}] (paired sequence-cluster bootstrap, run 0)")
                        v = v + vband(d, pth + "/verdict")
                        if kind == "estimated" and gate != "PASS":
                            v = L["L-E8-VOID"] + f" (computed: {c['verdict']})"
                    rows.append({"id": rid, "kind": title, "component": acc, "decoy": dec, "alpha": a,
                                 "statistic": stat, "verdict": v, "labels": psm_label})

    # ---------------------------------------------------------------- descriptive rows
    for acc in sorted(m["datasets"]):
        d = m["datasets"][acc]
        rows.append({"id": f"D-ids-{acc}", "kind": "phospho PSMs localized (q <= 1%)", "component": acc,
                     "statistic": f"rank-1 PSMs {d.get('n_rank1_psms')}; peptidoforms q<=1% {d.get('n_peptidoforms_q01')}; "
                                  f"phospho PSMs {d.get('n_phospho_psms_localized')} (mono {d.get('n_mono')}, multi {d.get('n_multi')}); "
                                  f"truth-eligible {d.get('n_truth_eligible_psms')}; per run "
                                  + ", ".join(f"{r}: {v['n_phospho_psms_localized']}" for r, v in sorted(d.get('per_run', {}).items())),
                     "verdict": L["L-DESC"], "labels": psm_label})
        for dec in ("A", "L"):
            e = d.get(f"ours_{dec}", {}).get("eval", {})
            mu = e.get("multi")
            rows.append({"id": f"D-multi-{acc}-{dec}", "kind": "multi-phospho (reported separately, §6.5)", "component": acc,
                         "statistic": "none" if not mu else
                         f"sites {mu['n_sites']}, STY {mu['n_sty']}, decoy {mu['n_decoy']}, STY at FLR̂ 1% {mu.get('n_sty_at_flrhat_0.01')}, 5% {mu.get('n_sty_at_flrhat_0.05')}",
                         "verdict": L["L-DESC"], "labels": psm_label})
            t = d.get("t2", {}).get(f"ours_{dec}", {})
            for a in ALPHAS:
                s = t.get(f"strata_{a}")
                if s:
                    rows.append({"id": f"D-strata-{acc}-{dec}-{a}", "kind": "acceptor within 2 residues (§6.5 strata)", "component": acc,
                                 "statistic": "; ".join(f"{k}: n_elig {v['n_eligible_at_global_thr']} true FLR {pct(v['true_flr_at_global_thr'])} "
                                                        f"(own-stratum FLR̂: n {v['n_eligible_at_own_stratum_thr']}, true {pct(v['true_flr_at_own_stratum_thr'])})"
                                                        for k, v in sorted(s.items())),
                                 "verdict": L["L-DESC"], "labels": psm_label})
            sl = t.get("site_level")
            if sl:
                rows.append({"id": f"D-site-{acc}-{dec}", "kind": "site-level FLR̂ (D8 secondary)", "component": acc,
                             "statistic": "; ".join(f"FLR̂ {a}: STY sites {sl[a]['n_sty_sites']}, eligible {sl[a]['n_eligible']}, true {pct(sl[a]['true_flr'])}"
                                                    for a in ALPHAS if a in sl),
                             "verdict": L["L-DESC"], "labels": psm_label})
        tr = m["truth_reports"].get(acc, {})
        nophos = [p for p, v in tr.items() if isinstance(v, dict) and str(v.get("phospho_truth", "")).startswith("none")]
        if nophos:
            pr = d.get("per_run", {})
            rows.append({"id": f"D-nonphospho-{acc}", "kind": "phospho PSMs at q<=1% in pools without phospho (kit design)",
                         "component": acc,
                         "statistic": ", ".join(f"{p}: {pr.get(tr[p]['run'], {}).get('n_phospho_psms_localized')}" for p in sorted(nophos)),
                         "verdict": L["L-DESC"], "labels": psm_label})
    rt = m.get("runtime_per_run", {})
    for r in sorted(rt):
        v = rt[r]
        rows.append({"id": f"S-run-{r}", "kind": "runtime per run", "component": r,
                     "statistic": f"convert {f(v.get('convert_wall_s'), 1)} s (MS2 {v.get('n_ms2')}), search {f(v.get('search_wall_s'), 2)} s",
                     "verdict": L["L-DESC"], "labels": []})
    for acc in sorted(m["datasets"]):
        d = m["datasets"][acc]
        rows.append({"id": f"S-loc-{acc}", "kind": "localization throughput", "component": acc,
                     "statistic": "; ".join(f"{k} {f(d[k].get('time_s'), 1)} s ({f(d[k].get('psms_per_s'), 0)} PSM/s)"
                                            for k in ("ours_A", "ours_L", "pyascore_A", "pyascore_L") if k in d and "time_s" in d[k]),
                     "verdict": L["L-DESC"], "labels": []})
    tm = m.get("timing", {})
    rows.append({"id": "S-total", "kind": "driver total", "statistic":
                 f"wall {f(tm.get('total_wall_s'), 0)} s, core-s {f(tm.get('total_core_s'), 0)}; stages " +
                 ", ".join(f"{k} {f(v.get('wall_s'), 0)} s" for k, v in sorted(tm.get("stages", {}).items())),
                 "verdict": L["L-DESC"], "labels": []})
    return {"dataset": dsn, "look_id": m["meta"]["look_id"], "git_sha": m["meta"]["git_sha"], "e8_gate": gate,
            "incomplete": bool(miss), "missing": miss, "rows": rows,
            "labels": {k: L[k] for k in sorted(L)}}


def render(rep) -> str:
    out = [f"# T2 claim table ({rep['dataset']}, look {rep['look_id']}, git {rep['git_sha'][:12]})", ""]
    if rep["incomplete"]:
        out += ["**look failed: incomplete** — missing: " + "; ".join(rep["missing"]), ""]
    out += ["Run 0 of our localizer is the result of record; [a–b over N runs] is the empirical min–max over N repeated "
            "localizations of the same PSMs (frozen path, CUDA atomic sums break exact ties at random; PI D35); it is a range, not a confidence bound (on dev, M4's own earlier draw fell up to 2 sites outside an 11-run range). "
            "pyAscore is deterministic.", ""]
    out += [f"**E8 gate: {rep['e8_gate']}**", "",
            "Labels (apply to every row that lists them): " + "; ".join(f"**{k}**: {v}" for k, v in rep["labels"].items()), "",
            "| id | claim | statistic | verdict | labels |", "|---|---|---|---|---|"]
    for r in rep["rows"]:
        out.append(f"| {r['id']} | {r['kind']} | {r['statistic']} | {r['verdict']} | {', '.join(r['labels'])} |")
    return "\n".join(out) + "\n"


def write_claims(metrics_path, cfg, out_dir):
    m = json.load(open(metrics_path))
    rep = build(m, cfg)
    rep["metrics_sha256"] = hashlib.sha256(open(metrics_path, "rb").read()).hexdigest()
    with open(os.path.join(out_dir, "claims.json"), "w") as fh:
        json.dump(rep, fh, indent=1, sort_keys=True)
    with open(os.path.join(out_dir, "claims.md"), "w") as fh:
        fh.write(render(rep))
    return rep


def main():
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("metrics")
    ap.add_argument("--out", default=None)
    ap.add_argument("--config", default=os.path.join(ROOT, "config", "t2.yaml"))
    a = ap.parse_args()
    rep = write_claims(a.metrics, yaml.safe_load(open(a.config)), a.out or os.path.dirname(os.path.abspath(a.metrics)))
    print("look failed: incomplete" if rep["incomplete"] else f"E8 gate: {rep['e8_gate']}")
    sys.exit(1 if rep["incomplete"] else 0)


if __name__ == "__main__":
    main()
