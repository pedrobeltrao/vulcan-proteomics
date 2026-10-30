"""Markdown tables from the T1 driver metrics.json."""
from __future__ import annotations

LEVELS = ("peptidoform", "sequence", "protein")


def _f(x, n=3):
    return "-" if x is None else f"{x:.{n}f}"


def tables(m: dict) -> str:
    out = [f"# T1 driver: {m['meta']['dataset']} ({m['meta']['role']}), runs={len(m['meta']['runs'])}, "
           f"decoy seeds {m['meta']['decoy_seeds']}, entrapment seeds {m['meta']['entrap_seeds']}", ""]
    out.append("## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)")
    out.append("")
    out.append("| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |")
    out.append("|---|---|---|---|---|---|---|")
    for meth, rec in m["methods"].items():
        if "ids_dataset_level" not in rec:
            out.append(f"| {meth} | missing {rec.get('missing_replicates')} | | | | | |")
            continue
        for lv in LEVELS:
            ids = rec["ids_dataset_level"][lv]
            v1 = rec["verdicts"][lv]["0.01"]["overall"]
            v5 = rec["verdicts"][lv]["0.05"]["overall"]
            out.append(f"| {meth} | {lv} | {ids['0.001']['mean']:.0f} | {ids['0.01']['mean']:.0f} ± {ids['0.01']['sd']:.0f} | "
                       f"{ids['0.05']['mean']:.0f} | {v1.get('a')} / {v1.get('b', '-')} | {v5.get('a')} / {v5.get('b', '-')} |")
    out += ["", "## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)", "",
            "| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |", "|---|---|---|---|---|---|---|---|"]
    for meth, rec in m["methods"].items():
        if "verdicts" not in rec:
            continue
        for lv in LEVELS:
            v = rec["verdicts"][lv]["0.01"]
            P = v["P"]
            C = v.get("C")
            out.append(f"| {meth} | {lv} | {_f(P['a']['mean_fdp'], 4)} | [{_f(P['a']['ci95'][0], 4)}, {_f(P['a']['ci95'][1], 4)}] | "
                       f"[{_f(P.get('b', P['a'])['ci95'][0], 4)}, {_f(P.get('b', P['a'])['ci95'][1], 4)}] | "
                       f"{_f(C['a']['mean_fdp'], 4) if C else '-'} | "
                       f"{'[' + _f(C['a']['ci95'][0], 4) + ', ' + _f(C['a']['ci95'][1], 4) + ']' if C else '-'} | {_f(P['a']['mean_lb'], 4)} |")
    out += ["", "## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; "
            "95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)", "",
            "| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |",
            "|---|---|---|---|---|---|"]
    for comp, rec in m.get("noninferiority", {}).items():
        for lv in LEVELS:
            d = rec[lv]
            o = d["own_q"]["0.01"]
            m1, m5 = d["matched_fdp"]["0.01"], d["matched_fdp"]["0.05"]
            pf = lambda x: (f"{_f(x['ratio'])} [{_f(x['ci95'][0])}, {_f(x['ci95'][1])}]"  # noqa: E731
                            + (" PASS" if x["lower_ci_ge_0.97"] else " fail"))
            out.append(f"| {comp} | {lv} | {pf(o)} | {'yes' if o['counts_for_noninferiority'] else 'no'} "
                       f"({o['comparator_verdict_deciding']}) | {pf(m1)} | {pf(m5)} |")
    if m.get("sage_rt_check"):
        out += ["", f"Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {m['sage_rt_check']}"]
    if m.get("sage_predict_rt_control"):
        out += ["", f"Sage RT on vs off (prod, first decoy seed): {m['sage_predict_rt_control']}"]
    t = m.get("driver_timing", {})
    if t:
        out += ["", f"## Driver timing: total {t['wall_s']:.0f} s wall, {t['core_s_self_and_children']:.0f} core-s", "",
                "| stage | wall s | core-s (self+children) | loadavg start/end |", "|---|---|---|---|"]
        for k, v in t["stages"].items():
            out.append(f"| {k} | {v['wall_s']:.0f} | {v['core_s_self_and_children']:.0f} | "
                       f"{v['loadavg_start']:.0f}/{v['loadavg_end']:.0f} |")
    return "\n".join(out) + "\n"
