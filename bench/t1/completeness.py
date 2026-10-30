"""Completeness rule of the A5 look (DESIGN §6.0d; reviews/amendment_A5_draft.md A5.5.4, normative).

A look is FAILED unless all of the following hold for the T1 metrics.json:
  - every method of A3_METHODS has all 25 replicates (evaluate.py writes a method record with `verdicts` only when no
    replicate is missing) and no `missing_replicates` key;
  - all runs are present: metrics["runs"] (and meta.runs) equal the expected run list, and every leave-one-run-out
    list of a3_secondary has one entry per run;
  - meta.a3 lists every A3.3.3 comparator;
  - every matched-FDP ratio has `fdp_reached` recorded;
  - every bootstrap's `n_boot_nan` is recorded (matched-FDP ratios: the only bootstrap that can produce NaN draws);
  - A5 additions (metrics["a5"]) exist and report complete (bench/t1/a5.py).
Checked twice: by the driver in look mode (exit non-zero) and by bench/t1/claims_a5.py ("look failed: incomplete").
`check` returns the list of violations; empty means complete.  It never looks at a value, only at presence.
"""
from __future__ import annotations

from . import a3 as A3

N_REPLICATES = 25            # 5 production + 15 arm P + 5 arm C (DESIGN §6.3 milestone design)


def _ratio_entries(ni: dict):
    for comp, byl in (ni or {}).items():
        for lv, d in (byl or {}).items():
            for kind in ("own_q", "matched_fdp"):
                for thr, v in ((d or {}).get(kind) or {}).items():
                    yield comp, lv, kind, thr, v


def check(metrics: dict, runs=None, methods=None, require_a5: bool = True) -> list:
    """Violations of the A5.5.4 completeness rule (empty list = complete).  `runs`: expected run names (the look's
    sealed files); `methods`: defaults to A3_METHODS."""
    bad = []
    methods = list(methods or A3.A3_METHODS)
    if not isinstance(metrics, dict):
        return ["metrics is not a JSON object"]
    mm = metrics.get("methods") or {}
    for m in methods:
        rec = mm.get(m)
        if rec is None:
            bad.append(f"method {m}: absent")
        elif "missing_replicates" in rec:
            bad.append(f"method {m}: missing_replicates {rec['missing_replicates']}")
        elif "verdicts" not in rec:
            bad.append(f"method {m}: no verdicts")
    got_runs = metrics.get("runs")
    if runs is not None:
        if sorted(got_runs or []) != sorted(runs):
            bad.append(f"runs: {got_runs} != expected {sorted(runs)}")
    if not got_runs:
        bad.append("runs: empty")
    if len(metrics.get("decoy_seeds") or []) * (1 + len(metrics.get("entrap_seeds") or [])) + \
            len(metrics.get("decoy_seeds") or []) != N_REPLICATES:
        bad.append(f"replicate design: decoy {metrics.get('decoy_seeds')} x entrap {metrics.get('entrap_seeds')} "
                   f"!= {N_REPLICATES} replicates")
    meta_a3 = ((metrics.get("meta") or {}).get("a3") or {})
    for c in A3.A3_COMPARATORS:
        if c not in (meta_a3.get("methods") or []):
            bad.append(f"meta.a3: comparator {c} not listed")
    nis = {"noninferiority": metrics.get("noninferiority")}
    for o, ni in (metrics.get("noninferiority_by_numerator") or {}).items():
        nis[f"noninferiority_by_numerator.{o}"] = ni
    if not metrics.get("noninferiority"):
        bad.append("noninferiority: absent")
    for num in (A3.PRIMARY, A3.SECONDARY):
        ni = (metrics.get("noninferiority_by_numerator") or {}).get(num)
        if not ni:
            bad.append(f"noninferiority_by_numerator.{num}: absent")
            continue
        for c in A3.A3_COMPARATORS + ([A3.SECONDARY] if num == A3.PRIMARY else []):
            if c not in ni:
                bad.append(f"noninferiority_by_numerator.{num}: comparator {c} absent")
    for name, ni in nis.items():
        for comp, lv, kind, thr, v in _ratio_entries(ni):
            if kind == "matched_fdp":
                if "fdp_reached" not in v:
                    bad.append(f"{name}.{comp}.{lv}.matched_fdp.{thr}: fdp_reached not recorded")
                if "n_boot_nan" not in v:
                    bad.append(f"{name}.{comp}.{lv}.matched_fdp.{thr}: n_boot_nan not recorded")
    a3s = metrics.get("a3_secondary")
    if not a3s:
        bad.append("a3_secondary: absent")
    else:
        R = len(got_runs or [])
        for num, byc in (a3s.get("loo") or {}).items():
            for comp, byl in byc.items():
                for lv, d in byl.items():
                    for kind, byt in d.items():
                        for thr, v in byt.items():
                            if R > 1 and len(v.get("loo") or []) != R:
                                bad.append(f"a3_secondary.loo.{num}.{comp}.{lv}.{kind}.{thr}: "
                                           f"{len(v.get('loo') or [])} LOO entries != {R} runs")
        if not a3s.get("loo"):
            bad.append("a3_secondary.loo: absent")
    if require_a5:
        a5 = metrics.get("a5")
        if not a5:
            bad.append("a5: absent")
        elif not a5.get("complete", False):
            bad.append(f"a5: incomplete {a5.get('missing')}")
    return bad
