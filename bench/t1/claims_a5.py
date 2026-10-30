#!/usr/bin/env python3
"""A5 claim script (DESIGN §6.0d; reviews/amendment_A5_draft.md A5.11.4, normative): maps the look's outputs to the
A5.9.2 claim table, the A5.7 voids, family S (A3.3.5 as amended by A5.8) and the A5.9.1 labels.

    python bench/t1/claims_a5.py --t1 <T1 out>/metrics.json --t5 <T5 out>/armM.json --out <dir> \
        [--ledger-t1-manifest-sha <sha256>] [--ledger-t5-manifest-sha <sha256>] [--base-commit <sha>] [--config config/a5.yaml]

Pure function of: T1 metrics.json (+ its manifest.sha256), T5 armM.json (+ manifest), config/a5.yaml (A5.4 gate
tables, C9 record, labels, A6 hook, comparator non-reproducibility), D23's committed A3.1.8 record, and the A6 status
file AS COMMITTED AT THE LOOK'S BASE COMMIT.  Writes claims.json, claims.md and claims.sha256 (sha256 of both).
Deterministic: sorted keys, no timestamps, no randomness.  It decides nothing that is not written in A5 or A3; where
the text leaves a reading open the stricter reading is taken (listed in STRICTER and in config/a5.yaml deviations).

Exit codes: 0 = written; 2 = "look failed: incomplete" (claims.md says so and lists what is missing); 3 = an input
hash does not match the ledger result line.

Decision rules (all from the text; references in brackets):
  * completeness [A5.5.4]: bench/t1/completeness.check + T5 armM completeness + (sealed) 6 T1 runs, 2 T5 runs,
    15 arm-M replicates.  Any violation -> look failed: incomplete.
  * a ratio is NOT EVALUABLE if it is missing, its point ratio is None, its CI is not finite, fdp_reached is not
    True for both sides (matched FDP), or n_boot_nan > 5% of the 2000 draws [A5.5.4].  A not-evaluable comparator
    counts as NOT BEATEN in family S and in P1 [A5.5.4, A5.8.4].
  * "beaten" = evaluable AND LB (lower 95% bootstrap bound) clears the bar AND, where the LOO backstop applies
    (family S, P1(i)), the leave-one-run-out minimum clears it.  A LOO minimum that is None (one run, or a LOO ratio
    undefined) does not clear any bar [A3.3.4; stricter].
  * v1 void [A5.7.1]: any "violated" (convention a) at 1%/5% in arm P (peptidoform, sequence, protein) or arm C
    (peptidoform, sequence) on T1; any "violated" T5 arm-M v1 verdict (peptidoform, 1%/5%, convention a); A6 FAILED.
    Claim-level: each matched-FDP ratio with v1 as numerator at peptidoform/sequence, f in {1%, 5%}, is void if the
    lower CI of v1's mean arm-C LB at v1's matched threshold is > f.
  * v2_sep void [A5.7.2]: any "violated" (convention a) in any arm of the overall or any merged-stratum peptidoform
    list at 1%/5% (recomputed here from the verdicts AND taken from a3_secondary.void_flags; union); any "violated"
    T5 arm-M verdict in the v2_sep overall or any merged-stratum list (peptidoform, 1%/5%, convention a); claim-level
    arm-C void as for v1 on the overall list.  A protein-group "violated" (arm P, 1% or 5%) stops family S at S3.
  * family S [A3.3.5 + A5.8]: S1..S5 in order on the full v2_sep list, each against EVERY A3.3.3 comparator (IUT), by
    LB and LOO minimum; tested only if every earlier one was rejected.
  * P1 [A5.9.2]: (i) peptidoform matched 1%, LB >= 0.97 and LOO min >= 0.97 vs every A1 comparator; (ii) peptidoform
    own q 1%, LB >= 0.97 vs every A1 comparator whose own-q verdict is "supported"; (ii) is "no supported comparator"
    if there is none, and P1 then rests on (i).
  * headline comparator [A5.9.4]: the A1 comparator with the lowest P1(i) point ratio (a None ratio counts as lowest).
  * A5.9.3 state: C9 record (config c9.verdict) and, in the GO branch, the A6 status file (A6.2 hook).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

VERSION = "claims-a5-1"
A6_STATUS_PATH = "results/2026-10-28_A6_hn1/a6_status.json"
A6_ENUM = ("PASS", "FAILED", "not evaluable", "invalid pool", "not completed", "label-free stop", "not tagged",
           "stopped: node time limit")
A6_ABSENT = "stopped: node time limit"
A6_PASS_DEADLINE = "2026-10-28T17:55:00+01:00"     # A6.8.2 step 5 / D31.3 (config a6.pass_deadline)

V1, V2 = "vulcan", "vulcan_v2_sep"
A1_COMPARATORS = ["sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
CHIM_COMPARATORS = ["sage_chim2", "sage_chim2+pred+mokapot"]
A3_COMPARATORS = A1_COMPARATORS + CHIM_COMPARATORS        # == bench/t1/a3.A3_COMPARATORS (asserted in tests)
NI = 0.97
N_BOOT = 2000
NAN_FRAC = 0.05
VOID_ALPHAS = ("0.01", "0.05")
LEVELS = ("peptidoform", "sequence", "protein")
OWN_ALPHAS = ("0.001", "0.01", "0.05")
FMATCH = ("0.01", "0.05")
SEALED_T1_RUNS, SEALED_T5_RUNS, SEALED_ARMM_REPS = 6, 2, 15
LQ, RQ = "“", "”"

# stricter readings taken where the text leaves a choice (also recorded as dated deviations in config/a5.yaml)
STRICTER = [
    "LOO minimum None (one run or an undefined LOO ratio) does not clear any bar (family S, P1(i)).",
    "Matched-FDP ratio not evaluable also when fdp_reached is not True on both sides or the CI is not finite "
    "(A5.5.4 names only >5% NaN draws).",
    "v2_sep void: union of a3_secondary.void_flags and a recomputation from the per-arm verdicts.",
    "Family S S1 with no 'supported' own-q comparator: S1's own-q part is not met, so S1 is not rejected and the "
    "sequence stops (A3.3.5 gives the 'not tested' rule for S5 only; A5.9.2's 'rests on (i)' is written for P1 only).",
    "A6 FAILED keeps the FAILED state even if its gate output cannot be verified (it can only void); this departs "
    "from A6.2's literal 'unverifiable -> stopped' and needs PI confirmation as a deviation. PASS needs a "
    "verifiable gate output and a step-5 commit that is an ancestor of the look base commit, else 'stopped: node "
    "time limit'.",
    "A6 PASS also needs its step-5 commit dated no later than the A6 clock's step-5 deadline (config a6."
    "pass_deadline, 17:55 CET 2026-10-28; A6.8.2/D31.3), else 'stopped: node time limit'.",
    "A6 is v1-only (D32): v2_sep wording (L-Q, L-R, title) is always the no-A6 wording; its C9 clause says so; the "
    "A5.7.2 'v2_sep FAILED a completed null-spectrum test' trigger has no input and cannot fire (stated in claims.md).",
    "In the A6 PASS state the designation is not changed by this script (A5.10.6 is a separate act); the L-D label "
    "is printed unchanged with that note.",
    "Arm M (T5) verdicts enter L-E for v1/v2_sep lists at peptidoform and sequence; 'not applicable to this unit' "
    "at protein level.  Arm M has no paired control, so '(consistent with control)' is never attached to it.",
    "Sealed run counts (6 T1, 2 T5) and 15 arm-M replicates are part of the completeness rule; sealed inputs "
    "require the ledger manifest hashes (run() refuses without them).",
    "Claim-level arm-C void also when the arm-C LB lower CI is missing or not finite.",
    "Headline comparator: the lowest FINITE P1(i) point ratio; a None ratio (0 comparator IDs) never headlines; "
    "ties by comparator name; with no finite ratio the first A1 comparator by name is named.",
    "A5.9.3 PASS wording: the '(…)' of L-S's pass clause is filled with '(A1, 2026-10-27: tail-G2, G1, G3)', the "
    "gates named in A5.4.",
    "Own-q ratios carry no n_boot_nan: evaluate.py floors the comparator denominator at 1e-9, so an own-q bootstrap "
    "draw cannot be NaN; A5.5.4's 'every bootstrap's n_boot_nan' is read as every matched-FDP bootstrap.",
    "A void (scorer-level or claim-level) is reported as 'VOID', never as 'not met'/'not rejected'; voided numbers "
    "are printed in full under captions titled 'VOID under §6.0d A5.7 (reason)'.",
]


# ------------------------------------------------------------------------------------------------ small helpers
def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_yaml(p):
    import yaml
    with open(p) as fh:
        return yaml.safe_load(fh)


def fmt(x, nd=4):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, (int,)):
        return str(x)
    if isinstance(x, float) and not math.isfinite(x):
        return "nan"
    return f"{x:.{nd}g}"


def pct(x):
    return "—" if x is None else f"{100 * x:.3g}%"


def word(v):
    """§6.3 verdict words only (code spelling 'not-evaluable' -> 'not evaluable')."""
    if v is None:
        return "not computed"
    return {"not-evaluable": "not evaluable"}.get(v, v)


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def g(d, *ks, default=None):
    for k in ks:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def q(s):
    return LQ + s + RQ


# ------------------------------------------------------------------------------------------------ A6 hook
def a6_state(root: str, base_commit: str | None, cfg: dict | None = None) -> dict:
    """The A6 status as committed at the look base commit (A6.2).  -> {"status", "step5_commit", "reason"}."""
    path = g(cfg or {}, "a6", "status_path", default=A6_STATUS_PATH)
    enum = tuple(g(cfg or {}, "a6", "enum", default=A6_ENUM))
    absent = {"status": A6_ABSENT, "step5_commit": None}
    if not base_commit:
        return {**absent, "reason": "no look base commit given"}
    r = subprocess.run(["git", "-C", root, "show", f"{base_commit}:{path}"], capture_output=True, text=True)
    if r.returncode:
        return {**absent, "reason": "status file not committed at the look base commit"}
    try:
        d = json.loads(r.stdout)
    except ValueError:
        return {**absent, "reason": "status file unparseable"}
    st = d.get("status") if isinstance(d, dict) else None
    if st not in enum:
        return {**absent, "reason": f"status outside the closed enum: {st!r}"}
    c5, gsha, gpath = d.get("step5_commit"), d.get("gate_output_sha256"), d.get("gate_output_path")
    malformed = (c5 is not None and not (isinstance(c5, str) and re.fullmatch(r"[0-9a-f]{7,40}", c5))) or \
                (gsha is not None and not (isinstance(gsha, str) and re.fullmatch(r"[0-9a-f]{64}", gsha)))
    if malformed and st != "FAILED":
        return {**absent, "reason": "malformed sha fields"}
    if st == "PASS":
        ok = bool(c5 and gsha and gpath)
        if ok:
            ok = subprocess.run(["git", "-C", root, "merge-base", "--is-ancestor", c5, base_commit]).returncode == 0
        if ok:
            r2 = subprocess.run(["git", "-C", root, "show", f"{base_commit}:{gpath}"], capture_output=True)
            ok = r2.returncode == 0 and hashlib.sha256(r2.stdout).hexdigest() == gsha
        if not ok:
            return {**absent, "reason": "PASS not verifiable (step-5 commit / gate output)"}
        dl = g(cfg or {}, "a6", "pass_deadline", default=A6_PASS_DEADLINE)
        t5 = subprocess.run(["git", "-C", root, "log", "-1", "--format=%cI", c5], capture_output=True,
                            text=True).stdout.strip()
        import datetime as _dt
        if not t5 or _dt.datetime.fromisoformat(t5) > _dt.datetime.fromisoformat(dl):
            return {**absent, "reason": f"PASS step-5 commit dated {t5 or '?'} after the A6 deadline {dl}"}
    return {"status": st, "step5_commit": c5 if isinstance(c5, str) else None, "reason": "read at base commit"}


def a59_state(cfg: dict, a6: dict) -> str:
    """A5.9.3 row: 'NO-GO' | 'GO, A6 PASS' | 'GO, A6 FAILED' | 'GO, A6 other'."""
    if g(cfg, "c9", "verdict") != "GO":
        return "NO-GO"
    return {"PASS": "GO, A6 PASS", "FAILED": "GO, A6 FAILED"}.get(a6["status"], "GO, A6 other")


# ------------------------------------------------------------------------------------------------ labels
def render_labels(cfg: dict, state: str, a6: dict) -> dict:
    L = {k: v["text"] for k, v in cfg["labels"].items()}
    c9 = cfg["c9"]
    if c9.get("verdict") == "GO":
        st = a6["status"]
        c9v1 = L["C9_clause_a6"].replace("PASS | FAILED | not evaluable | not completed (time box)", st) \
                                .replace("(commit c)", f"(commit {a6.get('step5_commit') or 'none'})")
    elif c9.get("verdict") == "NO-GO":
        c9v1 = L["C9_clause_nogo"].replace("= x <", f"= {c9.get('p_v1_minimum', 'x')} <") \
                                  .replace("commit c", f"commit {c9.get('commit')}")
    else:
        c9v1 = L["C9_clause_timebox"]
    c9v2 = (f"; C9 pilot: GO (commit {c9.get('commit')}); A6 gates v1 only (D32) and was not run for v2_sep"
            if c9.get("verdict") == "GO" else "; " + c9v1)
    out = {
        "L-G_v1": L["L-G_v1"].replace("{C9 clause}", c9v1),
        "L-G_v2_sep": L["L-G_v2_sep"].replace(" {C9 clause for v2_sep}", c9v2),
        "L-2": L["L-2"], "L-D": L["L-D"],
        "v2_sep_title": L["v2_sep_title"], "stratum1_title": L["stratum1_title"], "R4": L["R4"],
        "narrow_tolerance": L["narrow_tolerance"].replace("**", ""), "E7b_rank_caveat": L["E7b_rank_caveat"],
        "D14_limitation": L["D14_limitation"], "void_power": L["void_power"].replace("**", ""),
        "A1.8_statement": L["A1.8_statement"],
    }
    ls = L["L-S"]
    lq = L["L-Q"]
    lr_template = L["L-R"]
    last_ls = "not validated by null-spectrum tests: the only completed test (A1) FAILED"
    if state == "GO, A6 FAILED":
        ls = ls.replace("the only completed test (A1) FAILED",
                        L["L-S_last_clause_a6_verdict"].replace("<verdict>", "FAILED"))
    elif state == "GO, A6 PASS":
        ls = ls.replace(last_ls, L["L-S_pass_last_clause"].replace("(…)", "(A1, 2026-10-27: tail-G2, G1, G3)"))
        lq = lq.replace("estimator validity for this scorer is not established by null-spectrum tests (A1 FAILED)",
                        L["L-Q_pass_last_clause"])
        lr_template = None
    out["L-S"] = ls
    # A6 gates v1 only (D32): only v1's L-S / L-Q / L-R change with the A6 state; v2_sep keeps the no-A6 texts.
    out["L-Q_v1"] = lq
    out["L-Q_v2_sep"] = L["L-Q"]
    out["L-R_template"] = lr_template
    out["L-R_template_v2_sep"] = L["L-R"]
    out["L-R_pass"] = L["L-R_pass"]
    out["L-M_template"] = L["L-M"]
    out["L-E_template"] = L["L-E"]
    return out


def render_lr(lab: dict, p, qq, scorer="v1"):
    t = lab["L-R_template"] if scorer == "v1" else lab["L-R_template_v2_sep"]
    if t is None:
        return lab["L-R_pass"]
    ps = pct(p) if p is not None else "— (not defined for this pipeline)"
    return t.replace("p%", ps).replace("q%", pct(qq))


def render_lm(lab: dict, f: str, level: str, armc: dict | None):
    t = lab["L-M_template"].replace("FDP̂ = f;", f"FDP̂ = {pct(float(f))};")
    if level == "protein":
        return t.replace("arm C for the same list at the same threshold: mean FDP̂ = x%, LB lower CI = y% "
                         "[protein level: arm C not computed]", "arm C for the same list at the same threshold: "
                         "not computed at protein level")
    x = g(armc or {}, "mean_fdp_c")
    y = (g(armc or {}, "lb_ci95") or [None])[0]
    return t.replace("mean FDP̂ = x%", f"mean FDP̂ = {pct(x)}").replace("LB lower CI = y%", f"LB lower CI = {pct(y)}") \
            .replace(" [protein level: arm C not computed]", "")


def render_le(lab: dict, unit: str, alpha: str, p, c, m):
    return lab["L-E_template"].replace("arm P: <verdict>", f"arm P: {p}") \
        .replace("arm C: <verdict | not computed at protein level>", f"arm C: {c}") \
        .replace("arm M (T5): <verdict | not applicable to this unit>", f"arm M (T5): {m}") \
        .replace("(<unit>, α = <α>)", f"({unit}, α = {pct(float(alpha))})")


# ------------------------------------------------------------------------------------------------ verdict access
def arm_word(rec_verdicts: dict, level, alpha, arm, conv="a", with_control=True):
    v = g(rec_verdicts, level, alpha)
    if v is None:
        return "not computed"
    if arm == "C" and level == "protein":
        return "not computed at protein level"
    w = word(g(v, arm, conv, "verdict"))
    if with_control and arm == "P" and w == "inconclusive" and g(v, "consistent_with_control_P", conv):
        w += " (consistent with control)"
    return w


def arm_m_word(L: dict | None, unit, alpha, conv="a"):
    if unit == "protein":
        return "not applicable to this unit"
    v = g(L or {}, unit, alpha, conv)
    if v is None:
        return "not computed"
    return word(v.get("verdict"))           # arm M has no paired control: no "consistent with control" suffix


def violated(v):
    return v == "violated"


# ------------------------------------------------------------------------------------------------ ratio access
def ratio_rec(metrics, num, comp, level, kind, thr):
    return g(metrics, "noninferiority_by_numerator", num, comp, level, kind, thr)


def loo_min(metrics, num, comp, level, kind, thr):
    return g(metrics, "a3_secondary", "loo", num, comp, level, kind, thr, "loo_min")


def evaluable(rec: dict | None, kind: str):
    if not rec:
        return False, "missing"
    if rec.get("ratio") is None:
        return False, "no ratio (comparator has 0 IDs)"
    ci = rec.get("ci95") or []
    if len(ci) != 2 or not all(finite(x) for x in ci):
        return False, "CI not finite"
    if kind == "matched_fdp":
        if rec.get("fdp_reached") != [True, True]:
            return False, f"fdp_reached {rec.get('fdp_reached')}"
        nn = rec.get("n_boot_nan")
        if nn is None or nn > NAN_FRAC * N_BOOT:
            return False, f"n_boot_nan {nn} > {NAN_FRAC:.0%} of {N_BOOT}"
    return True, ""


def clears(x, bar, strict):
    return finite(x) and (x > bar if strict else x >= bar)


def beaten(metrics, num, comp, level, kind, thr, bar, strict, use_loo, ni_source=None, loo_source=None):
    rec = (ni_source or {}).get(comp, {}).get(level, {}).get(kind, {}).get(thr) if ni_source is not None \
        else ratio_rec(metrics, num, comp, level, kind, thr)
    ok, why = evaluable(rec, kind)
    lm = (g(loo_source or {}, comp, level, kind, thr, "loo_min") if loo_source is not None
          else loo_min(metrics, num, comp, level, kind, thr))
    row = {"comparator": comp, "level": level, "kind": kind, "thr": thr, "bar": (">" if strict else "≥") + str(bar),
           "ratio": g(rec or {}, "ratio"), "lb": (g(rec or {}, "ci95") or [None])[0],
           "ci95": g(rec or {}, "ci95"), "loo_min": lm if use_loo else None, "evaluable": ok,
           "not_evaluable_reason": why or None}
    if not ok:
        row["beaten"] = False
        row["why"] = f"not evaluable ({why}) = not beaten"
        return row
    b = clears(row["lb"], bar, strict) and (clears(lm, bar, strict) if use_loo else True)
    row["beaten"] = bool(b)
    row["why"] = "" if b else ("LB" if not clears(row["lb"], bar, strict) else "LOO minimum " +
                               ("undefined" if lm is None else "below bar"))
    return row


def own_q_supported(metrics, comp, alpha="0.01", level="peptidoform"):
    """§1: own-q comparisons use comparators whose own-q verdict (deciding convention) is 'supported'."""
    return g(metrics, "methods", comp, "verdicts", level, alpha, "deciding") == "supported"


# ------------------------------------------------------------------------------------------------ completeness
def completeness(metrics, armm, sealed: bool):
    bad = []
    try:
        from bench.t1 import completeness as CP
        bad += CP.check(metrics)
    except ImportError:                                    # pragma: no cover (module committed with A5)
        bad.append("bench.t1.completeness not importable")
    a5 = metrics.get("a5") or {}
    for k in ("stratum1", "armC_matched", "rank_change", "screened_armC"):
        if k not in a5:
            bad.append(f"a5.{k}: absent")
    for m in [V1, V2] + A3_COMPARATORS:
        for lv in ("peptidoform", "sequence"):
            for f in FMATCH:
                if g(a5, "armC_matched", m, lv, f) is None:
                    bad.append(f"a5.armC_matched.{m}.{lv}.{f}: absent")
    if not g(metrics, "a3_secondary", "strata", "merged"):
        bad.append("a3_secondary.strata.merged: absent")
    if g(metrics, "a3_secondary", "void_flags") is None:
        bad.append("a3_secondary.void_flags: absent")
    if not isinstance(armm, dict):
        bad.append("T5 armM.json: absent")
    else:
        if not armm.get("complete", False):
            bad.append(f"T5 armM: incomplete {armm.get('missing')}")
        for path in (("v1", "overall"), ("v2_sep", "overall"), ("v2_sep", "merged")):
            if not g(armm, *path):
                bad.append(f"T5 armM.{'.'.join(path)}: absent")
        for path in (("v1", "overall"), ("v2_sep", "overall")):
            for u in ("peptidoform", "sequence"):
                for al in VOID_ALPHAS:
                    if g(armm, *path, u, al, "a", "verdict") is None:
                        bad.append(f"T5 armM.{'.'.join(path)}.{u}.{al}: no verdict")
        for blk, L in (g(armm, "v2_sep", "merged") or {}).items():
            for al in VOID_ALPHAS:
                if g(L, "peptidoform", al, "a", "verdict") is None:
                    bad.append(f"T5 armM.v2_sep.merged.{blk}.peptidoform.{al}: no verdict")
    if sealed:
        if len(metrics.get("runs") or []) != SEALED_T1_RUNS:
            bad.append(f"T1 runs: {len(metrics.get('runs') or [])} != {SEALED_T1_RUNS}")
        if isinstance(armm, dict):
            if len(armm.get("runs") or []) != SEALED_T5_RUNS:
                bad.append(f"T5 runs: {len(armm.get('runs') or [])} != {SEALED_T5_RUNS}")
            if armm.get("n_replicates") != SEALED_ARMM_REPS:
                bad.append(f"T5 arm-M replicates: {armm.get('n_replicates')} != {SEALED_ARMM_REPS}")
            if armm.get("mode") != "sealed":
                bad.append(f"T5 armM mode {armm.get('mode')!r} != 'sealed'")
    return bad


# ------------------------------------------------------------------------------------------------ voids
def voids(metrics, armm, a6) -> dict:
    out = {"v1": [], "v2_sep": [], "claim_level": {}}
    vv = g(metrics, "methods", V1, "verdicts") or {}
    for al in VOID_ALPHAS:
        for lv in LEVELS:
            if violated(g(vv, lv, al, "P", "a", "verdict")):
                out["v1"].append(f"T1 arm P violated ({lv}, α = {pct(float(al))}, convention a)")
        for lv in ("peptidoform", "sequence"):
            if violated(g(vv, lv, al, "C", "a", "verdict")):
                out["v1"].append(f"T1 arm C violated ({lv}, α = {pct(float(al))}, convention a)")
        if violated(g(armm, "v1", "overall", "peptidoform", al, "a", "verdict")):
            out["v1"].append(f"T5 arm M violated (v1 overall, peptidoform, α = {pct(float(al))}, convention a)")
    if a6["status"] == "FAILED":
        out["v1"].append("A6 FAILED (completed null-spectrum test before the look; A5.7.1.3)")
    v2 = g(metrics, "methods", V2, "verdicts") or {}
    strata = g(metrics, "a3_secondary", "strata", "merged") or {}
    for al in VOID_ALPHAS:
        for arm in ("P", "C"):
            if violated(g(v2, "peptidoform", al, arm, "a", "verdict")):
                out["v2_sep"].append(f"T1 arm {arm} violated (v2_sep overall peptidoform, α = {pct(float(al))})")
            for k in sorted(strata):
                if violated(g(strata, k, "verdicts", al, arm, "a", "verdict")):
                    out["v2_sep"].append(f"T1 arm {arm} violated (v2_sep merged block containing stratum {k}, "
                                         f"peptidoform, α = {pct(float(al))})")
        if violated(g(armm, "v2_sep", "overall", "peptidoform", al, "a", "verdict")):
            out["v2_sep"].append(f"T5 arm M violated (v2_sep overall, peptidoform, α = {pct(float(al))})")
        for blk in sorted(g(armm, "v2_sep", "merged") or {}):
            if violated(g(armm, "v2_sep", "merged", blk, "peptidoform", al, "a", "verdict")):
                out["v2_sep"].append(f"T5 arm M violated (v2_sep merged {blk}, peptidoform, α = {pct(float(al))})")
    for h in g(metrics, "a3_secondary", "void_flags", "a", "peptidoform_violated") or []:
        r = f"T1 a3_secondary void flag {h} (convention a)"
        out["v2_sep"].append(r)
    out["v2_sep_protein_violated"] = sorted(set(
        [al for al in VOID_ALPHAS if violated(g(v2, "protein", al, "P", "a", "verdict"))] +
        list(g(metrics, "a3_secondary", "void_flags", "a", "protein_group_overall_violated") or [])))
    # claim-level arm-C void (A5.7.1 (a)-(f)); numerators v1 (all lists) and v2_sep (overall list)
    for num in (V1, V2):
        out["claim_level"][num] = {}
        for lv in ("peptidoform", "sequence"):
            for f in FMATCH:
                rec = g(metrics, "a5", "armC_matched", num, lv, f)
                lo = (g(rec or {}, "lb_ci95") or [None])[0]
                void = (not finite(lo)) or lo > float(f)
                out["claim_level"][num][f"{lv}/{f}"] = {
                    "void": bool(void), "lb_lower_ci": lo, "mean_lb": g(rec or {}, "mean_lb"),
                    "mean_fdp_c": g(rec or {}, "mean_fdp_c"), "ours_q_at_match": g(rec or {}, "ours_q_at_match"),
                    "reason": ("arm-C LB lower CI not computed" if not finite(lo) else
                               f"arm-C mean LB lower CI {pct(lo)} > f = {pct(float(f))}" if void else "")}
    out["v1_void"] = bool(out["v1"])
    out["v2_sep_void"] = bool(out["v2_sep"])
    return out


# ------------------------------------------------------------------------------------------------ P1 / family S
def headline_choice(rows_i):
    """A5.9.4: the A1 comparator with the lowest P1(i) point ratio.  Lowest FINITE ratio; a None ratio (comparator
    with 0 IDs) never headlines; ties broken by comparator name; if no ratio is finite, the first A1 comparator by
    name is named and its ratio printed as not computed (recorded in STRICTER)."""
    fin = sorted((r["ratio"], r["comparator"]) for r in rows_i if finite(r["ratio"]))
    return fin[0][1] if fin else sorted(r["comparator"] for r in rows_i)[0]


def p1(metrics, vd):
    rows_i = [beaten(metrics, V1, c, "peptidoform", "matched_fdp", "0.01", NI, False, True) for c in A1_COMPARATORS]
    sup = [c for c in A1_COMPARATORS if own_q_supported(metrics, c)]
    rows_ii = [beaten(metrics, V1, c, "peptidoform", "own_q", "0.01", NI, False, False) for c in sup]
    cl = vd["claim_level"][V1]["peptidoform/0.01"]
    i_ok = all(r["beaten"] for r in rows_i) and not cl["void"]
    ii_ok = all(r["beaten"] for r in rows_ii) if sup else None
    met = i_ok and (ii_ok is not False)
    head = headline_choice(rows_i)
    if vd["v1_void"]:
        outcome, void_reason = "VOID", "; ".join(vd["v1"])
    elif cl["void"]:
        outcome = "VOID (claim-level, A5.7.1: peptidoform/1%)"
        void_reason = "claim-level, A5.7.1: peptidoform/1%: " + cl["reason"]
    else:
        outcome, void_reason = ("met" if met else "not met"), None
    return {"i": rows_i, "ii": rows_ii, "ii_supported": sup,
            "ii_note": None if sup else "no supported comparator; P1 rests on (i)",
            "claim_level_void_i": cl["void"], "i_met": i_ok, "ii_met": ii_ok,
            "outcome": outcome, "void_reason": void_reason,
            "headline_comparator": head,
            "headline_row": next(r for r in rows_i if r["comparator"] == head)}


S_SPEC = [
    ("S1", "non-inferiority", [("peptidoform", "matched_fdp", "0.01", NI, False, "all"),
                               ("peptidoform", "own_q", "0.01", NI, False, "supported")]),
    ("S2", "more IDs", [("peptidoform", "matched_fdp", "0.01", 1.0, True, "all")]),
    ("S3", "more protein groups", [("protein", "matched_fdp", "0.01", 1.0, True, "all")]),
    ("S4", "more IDs at 5%", [("peptidoform", "matched_fdp", "0.05", 1.0, True, "all")]),
    ("S5", "more IDs at own q", [("peptidoform", "own_q", "0.01", 1.0, True, "supported")]),
]


def _s_rows(metrics, parts):
    rows, notes, none_supported = [], [], False
    for lv, kind, thr, bar, strict, which in parts:
        comps = A3_COMPARATORS if which == "all" else [c for c in A3_COMPARATORS if own_q_supported(metrics, c)]
        if which == "supported" and not comps:
            notes.append("own q: no supported comparator")
            none_supported = True
            continue
        rows += [beaten(metrics, V2, c, lv, kind, thr, bar, strict, True) for c in comps]
    return rows, notes, none_supported


def family_s(metrics, vd):
    """A3.3.5 fixed sequence on the full v2_sep list (A5.2 F1).  -> (rows, void title or None).
    A void v2_sep still has every S row computed and printed in full, with outcome VOID (A5.7.3)."""
    if vd["v2_sep_void"]:
        title = "VOID under §6.0d A5.7 (" + "; ".join(vd["v2_sep"]) + ")"
        out = []
        for sid, hyp, parts in S_SPEC:
            rows, notes, _ = _s_rows(metrics, parts)
            out.append({"id": sid, "hypothesis": hyp, "tested": False, "outcome": "VOID", "rows": rows,
                        "notes": notes})
        return out, title
    out = []
    go = True
    stop = None
    for sid, hyp, parts in S_SPEC:
        if not go:
            out.append({"id": sid, "hypothesis": hyp, "tested": False, "outcome": "not tested (" + stop + ")",
                        "rows": []})
            continue
        if sid == "S3" and vd["v2_sep_protein_violated"]:
            stop = "stopped at S3: v2_sep protein-group list 'violated'"
            out.append({"id": sid, "hypothesis": hyp, "tested": False, "outcome": "not tested (" + stop + ")",
                        "rows": []})
            go = False
            continue
        rows, notes, none_sup = _s_rows(metrics, parts)
        if sid == "S5" and none_sup:                       # A3.3.5: "If none is supported, S5 is not tested"
            stop = "S5 not tested: no supported comparator"
            out.append({"id": sid, "hypothesis": hyp, "tested": False, "outcome": "not tested (no supported "
                        "comparator)", "rows": [], "notes": notes})
            go = False
            continue
        cl_void = []
        for lv, kind, thr, *_ in parts:
            if kind == "matched_fdp" and lv in ("peptidoform", "sequence") and \
                    vd["claim_level"][V2][f"{lv}/{thr}"]["void"]:
                cl_void.append(f"{lv}/{thr}: " + vd["claim_level"][V2][f"{lv}/{thr}"]["reason"])
        # S1 with no supported own-q comparator: the own-q part cannot be met (stricter reading; A3.3.5 gives the
        # "not tested" rule for S5 only, and A5.9.2's "rests on (i)" is written for P1 only)
        s1_unmet = sid == "S1" and none_sup
        if s1_unmet:
            notes.append("S1 own-q part not met: no supported comparator (stricter reading)")
        rej = bool(rows) and all(r["beaten"] for r in rows) and not cl_void and not s1_unmet
        out.append({"id": sid, "hypothesis": hyp, "tested": True, "rows": rows, "notes": notes,
                    "claim_level_void": cl_void,
                    "outcome": ("VOID (claim-level, A5.7.2: " + "; ".join(cl_void) + ")") if cl_void else
                               ("rejected (claim holds)" if rej else "not rejected")})
        if not rej:
            go = False
            stop = f"{sid} not rejected" if not cl_void else f"{sid} VOID (claim-level)"
    return out, None


# ------------------------------------------------------------------------------------------------ markdown
def ratio_table_rows(metrics, num, comps, vd, lab, armc_all):
    lines = ["| comparator | unit | threshold | ratio | 95% CI | LOO min | evaluable | claim-level void | "
             "label | arm C (num) mean FDP̂ / LB low | arm C (comp) mean FDP̂ / LB low | comparator non-reproducibility |",
             "|" + "---|" * 12]
    nonrep = lab.get("_nonrep", {})
    for c in comps:
        for lv in LEVELS:
            for kind, thrs in (("matched_fdp", FMATCH), ("own_q", OWN_ALPHAS)):
                for t in thrs:
                    rec = ratio_rec(metrics, num, c, lv, kind, t)
                    ok, why = evaluable(rec, kind)
                    lm = loo_min(metrics, num, c, lv, kind, t)
                    cv = ""
                    if kind == "matched_fdp" and lv != "protein":
                        cv = "VOID" if vd["claim_level"].get(num, {}).get(f"{lv}/{t}", {}).get("void") else "no"
                    an, ac = g(armc_all, num, lv, t), g(armc_all, c, lv, t)
                    fa = (lambda a: "—" if lv == "protein" or kind != "matched_fdp" else
                          f"{pct(g(a or {}, 'mean_fdp_c'))} / {pct((g(a or {}, 'lb_ci95') or [None])[0])}")
                    nr = nonrep.get(c)
                    nrs = ", ".join(f"{k} ≤ {v:g}" for k, v in sorted(nr.items())) if nr else "—"
                    lbl = "[L-M]" if kind == "matched_fdp" else "[L-Q]"
                    lines.append(f"| {c} | {lv} | {'matched FDP̂' if kind == 'matched_fdp' else 'own q'} "
                                 f"{pct(float(t))} | {fmt(g(rec or {}, 'ratio'))} | "
                                 f"[{', '.join(fmt(x) for x in (g(rec or {}, 'ci95') or []))}] | {fmt(lm)} | "
                                 f"{'yes' if ok else 'no: ' + why} | {cv or '—'} | {lbl} | {fa(an)} | {fa(ac)} | {nrs} |")
    return lines


def le_all(res, key):
    return " ".join(f"[L-E] {q(res['le'][key][lv][a])}" for lv in LEVELS for a in VOID_ALPHAS)


def build_md(res, lab, metrics, armm, cfg):
    L = []
    w = L.append
    st = res["state"]
    p = res["P1"]
    gate_tags = f"[L-G v1] {q(lab['L-G_v1'])}"
    head = p["headline_row"]
    w("# A5 claim table (bench/t1/claims_a5.py; DESIGN §6.0d)")
    w("")
    w(f"Generated by `{VERSION}` from the look outputs; copied into REPORT.md verbatim (A5.2 F9). "
      "Nothing in this file may be edited by hand.")
    w("")
    w("## Headline")
    w("")
    if p["outcome"].startswith("VOID"):
        p1_word = f"P1 (v1 no loss) is VOID under §6.0d A5.7 ({p['void_reason']})"
        void_ratio = " (VOID under §6.0d A5.7; printed in full, not a claim)"
    else:
        p1_word = {"met": "P1 (v1 no loss) is met", "not met": "P1 (v1 no loss) is not met"}[p["outcome"]]
        void_ratio = ""
    w(f"[L-S] {q(lab['L-S'])} {p1_word}; against the A1 comparator with the lowest P1(i) point ratio, "
      f"{head['comparator']}, v1's peptidoform ratio at matched arm-P FDP̂ 1%{void_ratio} is {fmt(head['ratio'])} "
      f"(95% CI [{', '.join(fmt(x) for x in (head['ci95'] or []))}], leave-one-run-out minimum {fmt(head['loo_min'])}) "
      f"[L-M] {q(res['lm_head'])}")
    w("")
    w("## Abstract")
    w("")
    w(f"[L-S] {q(lab['L-S'])}")
    w("")
    w(f"[L-D] {q(lab['L-D'])}")
    w("")
    w("## State (A5.9.3)")
    w("")
    w(f"- C9 record: {cfg['c9'].get('verdict')} (commit {cfg['c9'].get('commit')}, "
      f"{cfg['c9'].get('verdict_file')} sha256 {cfg['c9'].get('verdict_file_sha256')})")
    w(f"- A6 status (read at the look base commit {res['inputs'].get('base_commit') or 'none'}): "
      f"{res['a6']['status']} ({res['a6']['reason']})")
    w(f"- A5.9.3 row: **{st}**")
    w("- A3.1.8 (i)–(iii) hold (D23; record verified by sha256: "
      f"{res['inputs'].get('a3_1_8_sha256')}).")
    w("- A5.7.2 'v2_sep FAILED a completed null-spectrum test before the look': no such test exists for v2_sep "
      "(A6 gates v1 only, D32), so this trigger has no input and does not fire.")
    if st == "GO, A6 PASS":
        w("- Designation: unchanged by this script; a change happens only under A5.10.6, as a separate dated act.")
    w("")
    w("## Gate status (G; A5.4, verbatim)")
    for sc in ("v1", "v2_sep"):
        w("")
        w(f"**Table G-{sc}: {cfg['gates'][sc]['title']}** [no FDR-filtered numbers]")
        w("| Gate | Data | Result | Source |")
        w("|---|---|---|---|")
        for r in cfg["gates"][sc]["rows"]:
            w(f"| {r['gate']} | {r['data']} | {r['result']} | {r['source']} |")
    w("")
    w(f"[L-G v1] {q(lab['L-G_v1'])}")
    w("")
    w(f"[L-G v2_sep] {q(lab['L-G_v2_sep'])}")
    w("")
    w("## Claim table (A5.9.2; every row printed whatever its outcome)")
    w("")
    w("**Table A5.9.2: claim rows** [no FDR-filtered numbers]")
    w("| # | claim | outcome | labels | void |")
    w("|---|---|---|---|---|")
    for r in res["rows"]:
        w(f"| {r['id']} | {r['claim']} | {r['outcome']} | {r['labels']} | {r['void']} |")
    w("")
    # ---------------- v1
    w("## v1 primary (static_gbm_pred), reported first and in full")
    w("")
    w(f"Every number in this section carries [L-R] {q(res['lr_v1'])}")
    w("")
    vt1 = f" — {lab['void_title_v1']}" if res["voids"]["v1_void"] else ""
    cl1 = [k for k, d in sorted(res["voids"]["claim_level"][V1].items()) if d["void"]]
    vt1_i = vt1 or (f" — VOID under §6.0d A5.7 (claim-level, A5.7.1: peptidoform/1%: "
                    f"{res['voids']['claim_level'][V1]['peptidoform/0.01']['reason']})"
                    if res["voids"]["claim_level"][V1]["peptidoform/0.01"]["void"] else "")
    vt1_p2 = vt1 or (f" — VOID under §6.0d A5.7 (claim-level, A5.7.1: matched ratios at {', '.join(cl1)}; "
                     f"rows marked VOID)" if cl1 else "")
    if vt1:
        w(f"**{lab['void_title_v1']}**: every v1 FDR-based number below is printed in full and is not a claim.")
        w("")
    w("### P1: v1 no loss (§1)")
    w("")
    w(f"**Table P1(i){vt1_i}: peptidoform, matched arm-P FDP̂ 1%, bar LB ≥ 0.97 and LOO minimum ≥ 0.97, every A1 comparator** "
      f"{gate_tags} [L-E] {q(res['le']['v1']['peptidoform']['0.01'])} [L-M] {q(res['lm_head'])} "
      f"[L-R] {q(res['lr_v1'])}")
    w("| comparator | ratio | 95% CI | LOO min | evaluable | beaten | why |")
    w("|---|---|---|---|---|---|---|")
    for r in p["i"]:
        w(f"| {r['comparator']} | {fmt(r['ratio'])} | [{', '.join(fmt(x) for x in (r['ci95'] or []))}] | "
          f"{fmt(r['loo_min'])} | {r['evaluable']} | {r['beaten']} | {r['why'] or '—'} |")
    w("")
    w(f"Claim-level arm-C void at peptidoform/1%: {res['voids']['claim_level'][V1]['peptidoform/0.01']['void']} "
      f"({res['voids']['claim_level'][V1]['peptidoform/0.01']['reason'] or 'lower CI ≤ f'}).")
    w("")
    if p["ii_supported"]:
        w(f"**Table P1(ii){vt1}: peptidoform, own q 1%, bar LB ≥ 0.97, A1 comparators with own-q verdict 'supported'** "
          f"{gate_tags} [L-E] {q(res['le']['v1']['peptidoform']['0.01'])} [L-Q] {q(lab['L-Q_v1'])} "
          f"[L-R] {q(res['lr_v1'])}")
        w("| comparator | ratio | 95% CI | evaluable | beaten | why |")
        w("|---|---|---|---|---|---|")
        for r in p["ii"]:
            w(f"| {r['comparator']} | {fmt(r['ratio'])} | [{', '.join(fmt(x) for x in (r['ci95'] or []))}] | "
              f"{r['evaluable']} | {r['beaten']} | {r['why'] or '—'} |")
    else:
        w("P1(ii): no supported comparator; P1 rests on (i).")
    w("")
    w(f"P1 outcome: **{p['outcome']}** (i: {p['i_met']}; ii: {p['ii_met'] if p['ii_supported'] else 'no supported comparator'}).")
    w("")
    w("### P2/P4: v1 ratios vs every comparator (all units, own q and matched FDP̂)")
    w("")
    w(f"P4 rows (sage_chim2, sage_chim2+pred+mokapot) are reported in the primary section (A3 Q2), not a claim: "
      f"{q(lab['narrow_tolerance'])}; {q(lab['E7b_rank_caveat'])}.")
    w("")
    w(f"**Table P2{vt1_p2}: v1 / comparator** {gate_tags} {le_all(res, 'v1')} [L-M] {q(res['lm_head'])} (matched rows; "
      f"arm C columns per row) [L-Q] {q(lab['L-Q_v1'])} (own-q rows) [L-R] {q(res['lr_v1'])}")
    L += ratio_table_rows(metrics, V1, A3_COMPARATORS + [V2], res["voids"], lab, g(metrics, "a5", "armC_matched"))
    w("")
    w("### P3: v1 IDs at own q")
    w("")
    w(f"**Table P3{vt1}: v1 IDs (dataset level, mean over decoy shuffles)** {gate_tags} {le_all(res, 'v1')} "
      f"[L-Q] {q(lab['L-Q_v1'])} [L-R] {q(res['lr_v1'])}")
    w("| unit | own q | IDs (mean) | per decoy seed |")
    w("|---|---|---|---|")
    for lv in LEVELS:
        for a in OWN_ALPHAS:
            d = g(metrics, "methods", V1, "ids_dataset_level", lv, a) or {}
            w(f"| {lv} | {pct(float(a))} | {fmt(d.get('mean'))} | {d.get('per_decoy_seed')} |")
    w("")
    w("### P5/P6: v1 entrapment verdicts, per arm (literal §6.3 words; convention (a) decides, (b) reported)")
    w("")
    w("**Table P5: v1 T1 arms P and C; T5 arm M** [no FDR-filtered numbers]")
    w("| unit | α | arm P (a) | arm P (b) | arm C (a) | arm C (b) | arm M T5 (a) | arm M T5 (b) |")
    w("|---|---|---|---|---|---|---|---|")
    vv = g(metrics, "methods", V1, "verdicts") or {}
    for lv in LEVELS:
        for a in sorted(vv.get(lv, {}), key=float):
            w(f"| {lv} | {pct(float(a))} | {arm_word(vv, lv, a, 'P', 'a')} | {arm_word(vv, lv, a, 'P', 'b')} | "
              f"{arm_word(vv, lv, a, 'C', 'a')} | {arm_word(vv, lv, a, 'C', 'b')} | "
              f"{arm_m_word(g(armm, 'v1', 'overall'), lv, a, 'a')} | {arm_m_word(g(armm, 'v1', 'overall'), lv, a, 'b')} |")
    w("")
    w("Hyperscore reference through the v1 pipeline on T5 (reported only, never used to excuse a violation):")
    w("")
    w("**Table P6-ref: hyperscore (v1 pipeline) T5 arm M** [no FDR-filtered numbers]")
    w("| unit | α | arm M (a) | arm M (b) |")
    w("|---|---|---|---|")
    for lv in ("peptidoform", "sequence"):
        for a in OWN_ALPHAS:
            w(f"| {lv} | {pct(float(a))} | {arm_m_word(g(armm, 'hyperscore_v1', 'overall'), lv, a, 'a')} | "
              f"{arm_m_word(g(armm, 'hyperscore_v1', 'overall'), lv, a, 'b')} |")
    w("")
    w("### P7: v1 quant; P8: speed")
    w("")
    w(f"- P7: {res['p7']}")
    w(f"- P8: {res['p8']}")
    w("")
    # ---------------- v2_sep
    w(f"## {lab['v2_sep_title']} (v2_sep, family S)")
    w("")
    w(f"[L-2] {q(lab['L-2'])}")
    w("")
    w(f"[L-G v2_sep] {q(lab['L-G_v2_sep'])}")
    w("")
    vt2 = f" — {res['S_void']}" if res["S_void"] else ""
    cl2 = [k for k, d in sorted(res["voids"]["claim_level"][V2].items()) if d["void"]]
    vt2_r3 = vt2 or (f" — VOID under §6.0d A5.7 (claim-level, A5.7.2: matched ratios at {', '.join(cl2)}; "
                     f"rows marked VOID)" if cl2 else "")
    vt2_s = vt2 or ("".join(f" — VOID under §6.0d A5.7 ({s_['id']} claim-level, A5.7.2: "
                            f"{'; '.join(s_['claim_level_void'])})" for s_ in res["S"] if s_.get("claim_level_void")))
    if vt2:
        w(f"**{res['S_void']}**: every v2_sep FDR-based number below is printed in full and is not a claim.")
        w("")
    w(f"**Table S{vt2_s}: {lab['v2_sep_title']}** [L-G v2_sep] {q(lab['L-G_v2_sep'])} "
      f"[L-E] {q(res['le']['v2_sep']['peptidoform']['0.01'])} (S1, S2, S5) "
      f"[L-E] {q(res['le']['v2_sep']['protein']['0.01'])} (S3) "
      f"[L-E] {q(res['le']['v2_sep']['peptidoform']['0.05'])} (S4) "
      f"[L-M] {q(res['lm_v2'])} (S1, S2) [L-M] {q(res['lm_v2_prot'])} (S3) [L-M] {q(res['lm_v2_05'])} (S4) "
      f"[L-Q] {q(lab['L-Q_v2_sep'])} (S1 own q, S5) [L-R] {q(res['lr_v2'])} [L-2] {q(lab['L-2'])}")
    w("| S | hypothesis | outcome | comparator | unit / threshold | bar | ratio | LB | LOO min | beaten | why |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for s in res["S"]:
        if not s["rows"]:
            w(f"| {s['id']} | {s['hypothesis']} | {s['outcome']} | — | — | — | — | — | — | — | "
              f"{'; '.join(s.get('notes') or []) or '—'} |")
        for r in s["rows"]:
            w(f"| {s['id']} | {s['hypothesis']} | {s['outcome']} | {r['comparator']} | {r['level']} "
              f"{'matched FDP̂' if r['kind'] == 'matched_fdp' else 'own q'} {pct(float(r['thr']))} | {r['bar']} | "
              f"{fmt(r['ratio'])} | {fmt(r['lb'])} | {fmt(r['loo_min'])} | {r['beaten']} | {r['why'] or '—'} |")
    w("")
    w(f"S3 dev expectation (D25), quoted whatever T1 shows: {q(cfg['labels']['D25_s3_dev']['text'])}. "
      f"Comparisons with sage_chim2 / C2: {q(lab['narrow_tolerance'])}; {q(lab['E7b_rank_caveat'])}.")
    w("")
    w("### R2: v2_sep per-merged-stratum entrapment verdicts (A3.3.6; literal; merged block named)")
    w("")
    w(f"Merge maps (production replicates): {json.dumps(res['merge_maps'], sort_keys=True)}")
    w("")
    w("**Table R2: v2_sep T1 arms P and C per merged block and raw stratum; T5 arm M** [no FDR-filtered numbers]")
    w("| list | stratum | α | arm P (a) | arm C (a) | arm P (b) | arm C (b) |")
    w("|---|---|---|---|---|---|---|")
    sv = g(metrics, "a3_secondary", "strata") or {}
    for kd in ("merged", "raw"):
        for k in sorted(sv.get(kd, {})):
            vv2 = {"x": sv[kd][k]["verdicts"]}
            for a in VOID_ALPHAS:
                pa = word(g(vv2, "x", a, "P", "a", "verdict"))
                lab_ne = " (not evaluable at α; no stratum-specific statement)" if pa == "not evaluable" else ""
                w(f"| {kd}{' (reported)' if kd == 'raw' else ''} | {res['block_names'].get(k, k) if kd == 'merged' else k} | "
                  f"{pct(float(a))} | {pa}{lab_ne} | {word(g(vv2, 'x', a, 'C', 'a', 'verdict'))} | "
                  f"{word(g(vv2, 'x', a, 'P', 'b', 'verdict'))} | {word(g(vv2, 'x', a, 'C', 'b', 'verdict'))} |")
    w("")
    w("**Table R2-T1-overall: v2_sep overall list, T1 arms P and C, T5 arm M** [no FDR-filtered numbers]")
    w("| unit | α | arm P (a) | arm C (a) | arm M T5 (a) | arm P (b) | arm C (b) | arm M T5 (b) |")
    w("|---|---|---|---|---|---|---|---|")
    v2 = g(metrics, "methods", V2, "verdicts") or {}
    for lv in LEVELS:
        for a in sorted(v2.get(lv, {}), key=float):
            w(f"| {lv} | {pct(float(a))} | {arm_word(v2, lv, a, 'P', 'a')} | {arm_word(v2, lv, a, 'C', 'a')} | "
              f"{arm_m_word(g(armm, 'v2_sep', 'overall'), lv, a, 'a')} | {arm_word(v2, lv, a, 'P', 'b')} | "
              f"{arm_word(v2, lv, a, 'C', 'b')} | {arm_m_word(g(armm, 'v2_sep', 'overall'), lv, a, 'b')} |")
    w("")
    w("**Table R2-T5: v2_sep T5 arm M per merged block (deciding) and raw stratum (reported); hyperscore reference** "
      "[no FDR-filtered numbers]")
    w("| list | block | unit | α | arm M (a) | arm M (b) |")
    w("|---|---|---|---|---|---|")
    for kd, key in (("merged", ("v2_sep", "merged")), ("raw (reported)", ("v2_sep", "raw")),
                    ("hyperscore v2 pipeline (reference)", ("hyperscore_v2",))):
        for blk, Lm in sorted((g(armm, *key) or {}).items()):
            for u in ("peptidoform", "sequence"):
                for a in OWN_ALPHAS:
                    w(f"| {kd} | {blk} | {u} | {pct(float(a))} | {arm_m_word(Lm, u, a, 'a')} | {arm_m_word(Lm, u, a, 'b')} |")
    w("")
    w("### R3: v2_sep / v1 and v2_sep vs every comparator, all units (reported)")
    w("")
    w(f"**Table R3{vt2_r3}: v2_sep / comparator** [L-G v2_sep] {q(lab['L-G_v2_sep'])} {le_all(res, 'v2_sep')} "
      f"[L-M] {q(res['lm_v2'])} (matched rows) [L-Q] {q(lab['L-Q_v2_sep'])} (own-q rows) [L-R] {q(res['lr_v2'])} "
      f"[L-2] {q(lab['L-2'])}")
    L += ratio_table_rows(metrics, V2, [V1] + A3_COMPARATORS, res["voids"], lab, g(metrics, "a5", "armC_matched"))
    w("")
    w(f"**Table R3-strata{vt2}: v2_sep production IDs per merged block** [L-G v2_sep] "
      f"{q(lab['L-G_v2_sep'])} {le_all(res, 'v2_sep')} (overall list; per block: Table R2) [L-Q] {q(lab['L-Q_v2_sep'])} [L-R] {q(res['lr_v2'])} [L-2] {q(lab['L-2'])}")
    w("| block | own q | IDs (mean over decoy shuffles) |")
    w("|---|---|---|")
    for k in sorted(sv.get("merged", {})):
        for a, d in sorted((sv["merged"][k].get("ids_prod_dataset_level") or {}).items()):
            w(f"| {res['block_names'].get(k, k)} | {pct(float(a))} | {fmt(d.get('mean'))} |")
    w("")
    w(f"### R1: {lab['stratum1_title']}")
    w("")
    w(f"**Table R1: {lab['stratum1_title']}** [L-G v2_sep] {q(lab['L-G_v2_sep'])} {le_all(res, 'v2_sep')} (overall list; per block: Table R2) "
      f"[L-M] {q(res['lm_v2'])} (matched rows) [L-Q] {q(lab['L-Q_v2_sep'])} (own-q rows) [L-R] {q(res['lr_v2'])} "
      f"[L-2] {q(lab['L-2'])}")
    w("| comparator | unit / threshold | ratio | 95% CI | LOO min |")
    w("|---|---|---|---|---|")
    s1n = g(metrics, "a5", "stratum1", "noninferiority") or {}
    s1l = g(metrics, "a5", "stratum1", "loo") or {}
    for c in sorted(s1n):
        for lv in LEVELS:
            for kind in ("matched_fdp", "own_q"):
                for t, rec in sorted((g(s1n, c, lv, kind) or {}).items()):
                    w(f"| {c} | {lv} {'matched FDP̂' if kind == 'matched_fdp' else 'own q'} {pct(float(t))} | "
                      f"{fmt(rec.get('ratio'))} | [{', '.join(fmt(x) for x in (rec.get('ci95') or []))}] | "
                      f"{fmt(g(s1l, c, lv, kind, t, 'loo_min'))} |")
    w("")
    w(f"### R4: {lab['R4']}")
    w("")
    w(f"Entrapment peptides excluded by the D18 screen: {g(metrics, 'a5', 'screened_armC', 'n_excluded')}. "
      "It never changes a verdict, a void or a claim.")
    w("")
    w(f"**Table R4: {lab['R4']}** [no FDR-filtered numbers]")
    w("| method | unit | α | screened arm C (a) | mean FDP̂ | screened arm C (b) |")
    w("|---|---|---|---|---|---|")
    for m, byl in sorted((g(metrics, "a5", "screened_armC", "methods") or {}).items()):
        for lv, bya in sorted(byl.items()):
            for a, d in sorted(bya.items(), key=lambda x: float(x[0])):
                w(f"| {m} | {lv} | {pct(float(a))} | {word(g(d, 'a', 'verdict'))} | {pct(g(d, 'a', 'mean_fdp'))} | "
                  f"{word(g(d, 'b', 'verdict'))} |")
    w("")
    w("### Comparator entrapment verdicts (their own; comparator rows are never void)")
    w("")
    w("**Table C: comparators, T1 arms P and C (convention a)** [no FDR-filtered numbers]")
    w("| method | unit | α | arm P | arm C |")
    w("|---|---|---|---|---|")
    for m in A3_COMPARATORS:
        vm = g(metrics, "methods", m, "verdicts") or {}
        for lv in LEVELS:
            for a in sorted(vm.get(lv, {}), key=float):
                w(f"| {m} | {lv} | {pct(float(a))} | {arm_word(vm, lv, a, 'P', 'a')} | {arm_word(vm, lv, a, 'C', 'a')} |")
    w("")
    w("## VOID record (A5.7)")
    w("")
    for who, key in (("v1", "v1"), ("v2_sep", "v2_sep")):
        rs = res["voids"][key]
        w(f"- {who}: " + (f"**VOID under §6.0d A5.7 ({'; '.join(rs)})**; the voided numbers are printed in full above"
                         if rs else "no void trigger"))
    for num in (V1, V2):
        for k, d in sorted(res["voids"]["claim_level"][num].items()):
            if d["void"]:
                w(f"- claim-level: {num} matched ratios at {k}: VOID under §6.0d A5.7 ({d['reason']})")
    w("")
    w("## Limitations (stated with the results)")
    w("")
    w(f"- {q(lab['void_power'])}")
    w(f"- {q(lab['D14_limitation'])}")
    w("- Speed is claims-grade only if the node was otherwise idle during the look (A5.5.1.7, A5.13.9).")
    w("- Comparator non-reproducibility (D33; reported, not a decision input): " +
      "; ".join(f"{m}: " + ", ".join(f"{k} ≤ {v:g}" for k, v in sorted(d.items()))
                for m, d in sorted(cfg.get("comparator_nonreproducibility", {}).get("values", {}).items())) + ".")
    w("")
    w("## Stricter readings taken by this script")
    w("")
    for s in STRICTER:
        w(f"- {s}")
    w("")
    return "\n".join(L)


# ------------------------------------------------------------------------------------------------ main
def _verify_manifest(path, want_sha, name):
    d = os.path.dirname(os.path.abspath(path))
    man = os.path.join(d, "manifest.sha256")
    if not os.path.exists(man):
        return f"{name}: manifest.sha256 missing next to {os.path.basename(path)}"
    got = sha256_file(man)
    if got != want_sha:
        return f"{name}: manifest.sha256 sha256 {got} != ledger {want_sha}"
    listed = {}
    for l in open(man):
        if l.strip():
            h, p = l.rstrip("\n").split("  ", 1)
            listed[os.path.normpath(p)] = h
    rel = os.path.normpath(os.path.relpath(os.path.abspath(path), d))
    if listed.get(rel) != sha256_file(path):
        return f"{name}: {rel} not in manifest.sha256 or sha256 differs"
    return None


def run(t1_metrics_path, t5_armM_path, a5_cfg_path, out_dir, ledger_results=None, base_commit=None, root=ROOT):
    """Write claims.json / claims.md / claims.sha256 into out_dir.  Returns (exit code, result dict)."""
    os.makedirs(out_dir, exist_ok=True)
    cfg = load_yaml(a5_cfg_path)
    inputs = {"t1_metrics_sha256": sha256_file(t1_metrics_path) if os.path.exists(t1_metrics_path) else None,
              "t5_armM_sha256": sha256_file(t5_armM_path) if os.path.exists(t5_armM_path) else None,
              "config_sha256": sha256_file(a5_cfg_path), "base_commit": base_commit}
    rec = cfg.get("a3_1_8_records") or {}
    if rec.get("path"):
        p = os.path.join(ROOT, rec["path"])          # committed in this repository (D23)
        inputs["a3_1_8_sha256"] = sha256_file(p) if os.path.exists(p) else None
        if inputs["a3_1_8_sha256"] != rec.get("sha256"):
            return _fail(out_dir, 3, ["A3.1.8 record sha256 mismatch"], inputs)
    if ledger_results:
        errs = [e for e in (_verify_manifest(t1_metrics_path, ledger_results.get("T1"), "T1"),
                            _verify_manifest(t5_armM_path, ledger_results.get("T5"), "T5")) if e]
        if errs:
            return _fail(out_dir, 3, errs, inputs)
    metrics = json.load(open(t1_metrics_path)) if os.path.exists(t1_metrics_path) else {}
    armm = json.load(open(t5_armM_path)) if os.path.exists(t5_armM_path) else None
    sealed = g(metrics, "meta", "role") == "sealed"
    missing = completeness(metrics, armm, sealed)
    if missing:
        return _fail(out_dir, 2, missing, inputs)
    if sealed and not ledger_results:
        return _fail(out_dir, 3, ["sealed inputs require the ledger manifest sha256s (--look-id)"], inputs)
    a6 = a6_state(root, base_commit, cfg)
    state = a59_state(cfg, a6)
    lab = render_labels(cfg, state, a6)
    lab["_nonrep"] = (cfg.get("comparator_nonreproducibility") or {}).get("values", {})
    vd = voids(metrics, armm, a6)
    lab["void_title_v1"] = ("VOID under §6.0d A5.7 (" + "; ".join(vd["v1"]) + ")") if vd["v1"] else ""
    P1 = p1(metrics, vd)
    S, S_void = family_s(metrics, vd)
    rc = g(metrics, "a5", "rank_change") or {}
    res = {"version": VERSION, "inputs": inputs, "state": state, "a6": a6, "voids": vd, "P1": P1, "S": S,
           "S_void": S_void, "stricter_readings": STRICTER}
    res["lr_v1"] = render_lr(lab, g(rc, V1, "frac_spectra_rank_gt1"), g(rc, V1, "frac_accepted_units_rank_gt1"))
    res["lr_v2"] = render_lr(lab, g(rc, V2, "frac_spectra_rank_gt1"), g(rc, V2, "frac_accepted_units_rank_gt1"),
                             "v2_sep")
    armc = g(metrics, "a5", "armC_matched") or {}
    res["lm_head"] = render_lm(lab, "0.01", "peptidoform", g(armc, V1, "peptidoform", "0.01"))
    res["lm_v2"] = render_lm(lab, "0.01", "peptidoform", g(armc, V2, "peptidoform", "0.01"))
    res["lm_v2_05"] = render_lm(lab, "0.05", "peptidoform", g(armc, V2, "peptidoform", "0.05"))
    res["lm_v2_prot"] = render_lm(lab, "0.01", "protein", None)
    res["le"] = {}
    for key, m, armkey in (("v1", V1, "v1"), ("v2_sep", V2, "v2_sep")):
        vv = g(metrics, "methods", m, "verdicts") or {}
        res["le"][key] = {lv: {a: render_le(lab, lv, a, arm_word(vv, lv, a, "P"), arm_word(vv, lv, a, "C"),
                                            arm_m_word(g(armm, armkey, "overall"), lv, a))
                               for a in VOID_ALPHAS} for lv in LEVELS}
    mm = g(metrics, "a3_secondary", "strata", "merge_maps_prod") or {}
    res["merge_maps"] = {r: g(v, "merge") for r, v in sorted(mm.items())}
    blocks = {}
    for r, v in sorted(mm.items()):
        for blk, members in (g(v, "merge") or {}).items():
            for k in members:
                kk = str(k) if int(k) < 3 else ">=3"
                blocks.setdefault(kk, set()).add("{" + ", ".join("≥3" if int(x) >= 3 else str(x) for x in members) + "}")
    res["block_names"] = {k: "merged stratum " + " / ".join(sorted(v)) for k, v in blocks.items()}
    res["p7"] = ("quant (ε, CV, DE pAUC) on q-filtered lists: " + (
        "see metrics.json key 'quant' [L-Q] " + q(lab["L-Q_v1"]) if metrics.get("quant") else
        "not recorded by the T1 driver (no 'quant' key in metrics.json); no quant claim"))
    dt = metrics.get("driver_timing") or {}
    res["p8"] = (f"driver wall {fmt(dt.get('wall_s'))} s, core-s {fmt(dt.get('core_s_self_and_children'))} "
                 f"(all stages, all methods, {len(metrics.get('runs') or [])} runs); per-stage timing in metrics.json "
                 "driver_timing; claims-grade only if the node was idle (A5.5.1.7)")
    vvoid = "VOID (" + "; ".join(vd["v1"]) + ")" if vd["v1_void"] else "—"
    v2void = "VOID (" + "; ".join(vd["v2_sep"]) + ")" if vd["v2_sep_void"] else "—"
    cl1 = [k for k, d in vd["claim_level"][V1].items() if d["void"]]
    cl2 = [k for k, d in vd["claim_level"][V2].items() if d["void"]]
    rows = [
        {"id": "P1", "claim": "v1 no loss (§1): (i) matched FDP̂ 1% LB and LOO min ≥ 0.97 vs every A1 comparator; "
         "(ii) own q 1% LB ≥ 0.97 vs supported A1 comparators", "outcome": P1["outcome"] +
         ("" if P1["ii_supported"] else "; (ii) no supported comparator, P1 rests on (i)"),
         "labels": "L-S, L-M (i), L-Q (ii), L-E, L-R, L-G", "void": vvoid if vd["v1_void"] else
         ("claim-level: peptidoform/1%" if P1["claim_level_void_i"] else "—")},
        {"id": "P2", "claim": "v1 vs each comparator at 5%, sequence and protein level, 0.1% own q",
         "outcome": "reported (Table P2)", "labels": "L-M or L-Q, L-E, L-R", "void": vvoid if vd["v1_void"] else
         ("claim-level: " + ", ".join(cl1) if cl1 else "—")},
        {"id": "P3", "claim": "v1 IDs at own q {0.1, 1, 5}%", "outcome": "reported (Table P3)",
         "labels": "L-Q, L-E, L-R", "void": vvoid},
        {"id": "P4", "claim": "v1 vs C2 and sage_chim2 (A3 Q2; not a claim)", "outcome": "reported (Table P2)",
         "labels": "L-M, narrow tolerance, E7b rank caveat", "void": vvoid},
        {"id": "P5", "claim": "v1 entrapment verdicts", "outcome": "literal (Table P5)", "labels": "literal",
         "void": "—"},
        {"id": "P6", "claim": "v1 arm M (T5)", "outcome": "literal (Table P5)", "labels": "literal", "void": "—"},
        {"id": "P7", "claim": "v1 quant on q-filtered lists", "outcome": res["p7"].split(":", 1)[1].strip()[:80],
         "labels": "L-Q", "void": vvoid},
        {"id": "P8", "claim": "speed (core-s, runs/h)", "outcome": "reported", "labels": "none", "void": "—"},
    ]
    for s in S:
        rows.append({"id": s["id"], "claim": f"v2_sep {s['hypothesis']} (family S; Secondary arm; not gate-validated)",
                     "outcome": s["outcome"], "labels": "title, L-2, L-M/L-Q, L-E, L-R, L-G (v2_sep)",
                     "void": v2void if vd["v2_sep_void"] else ("claim-level: " + "; ".join(s.get("claim_level_void") or [])
                                                               if s.get("claim_level_void") else "—")})
    rows += [
        {"id": "R1", "claim": "v2_sep stratum-1-only list", "outcome": "descriptive; not family S (Table R1)",
         "labels": "descriptive; not family S, L-2", "void": "—"},
        {"id": "R2", "claim": "v2_sep per-merged-stratum verdicts", "outcome": "literal (Table R2)",
         "labels": "literal; merged block named; raw ≥3 shown", "void": "—"},
        {"id": "R3", "claim": "v2_sep/v1 ratio, per-stratum IDs, sequence level", "outcome": "reported (Table R3)",
         "labels": "L-M or L-Q, L-2", "void": v2void if vd["v2_sep_void"] else
         ("claim-level: " + ", ".join(cl2) if cl2 else "—")},
        {"id": "R4", "claim": "D18 screened arm C", "outcome": "reported only (Table R4)",
         "labels": "secondary arm C (D18 screen); reported only", "void": "—"},
        {"id": "D", "claim": "designations", "outcome": "fixed (F3)", "labels": "L-D", "void": "—"},
        {"id": "G", "claim": "gate-status tables", "outcome": "A5.4 verbatim", "labels": "verbatim", "void": "—"},
    ]
    res["rows"] = rows
    res["labels_rendered"] = {k: v for k, v in lab.items() if isinstance(v, str) and not k.startswith("_")
                              and "template" not in k}
    res["labels_rendered"].update({"L-R v1": res["lr_v1"], "L-R v2_sep": res["lr_v2"], "L-M head": res["lm_head"],
                                   "L-M v2": res["lm_v2"],
                                   "L-M v2 5%": res["lm_v2_05"], "L-M v2 protein": res["lm_v2_prot"]})
    for key in res["le"]:
        for lv in res["le"][key]:
            for a, t in res["le"][key][lv].items():
                res["labels_rendered"][f"L-E {key} {lv} {a}"] = t
    md = build_md(res, lab, metrics, armm, cfg)
    for lbl in [cfg["labels"]["D25_s3_dev"]["text"]]:
        res["labels_rendered"]["D25_s3_dev"] = lbl
    res["gates"] = cfg["gates"]
    return _write(out_dir, res, md, 0)


def _write(out_dir, res, md, code):
    js = json.dumps(res, indent=1, sort_keys=True, default=str, allow_nan=True)
    open(os.path.join(out_dir, "claims.json"), "w").write(js + "\n")
    open(os.path.join(out_dir, "claims.md"), "w").write(md + ("\n" if not md.endswith("\n") else ""))
    with open(os.path.join(out_dir, "claims.sha256"), "w") as fh:
        for n in ("claims.json", "claims.md"):
            fh.write(f"{sha256_file(os.path.join(out_dir, n))}  {n}\n")
    return code, res


def _fail(out_dir, code, items, inputs):
    head = "look failed: incomplete" if code == 2 else "look failed: input hash mismatch"
    md = head + "\n\n" + "\n".join(f"- {x}" for x in items) + "\n"
    return _write(out_dir, {"version": VERSION, "failed": head, "items": items, "inputs": inputs}, md, code)


def ledger_look(root, look_id):
    """(ledger_results {"T1": manifest sha256, "T5": ...}, base commit) from the ledger lines of one look
    (bench/sealed_run.py format).  Refuses unless both datasets have an exit-0 result line with a manifest hash."""
    import re
    import subprocess
    sys.path.insert(0, os.path.join(root, "bench"))
    import sealed_run as SR
    _, looks = SR.look_state(root)
    L = looks.get(look_id)
    if not L or not SR.look_succeeded(L):
        sys.exit(f"claims_a5: look {look_id} is not a successful look in the ledger (both exit-0 result lines and a "
                 "passed leak test are required)")
    lr = {}
    for ds in ("T1", "T5"):
        m = re.search(r"manifest_sha256=([0-9a-f]{64})", L["results"][ds])
        if not m:
            sys.exit(f"claims_a5: no manifest sha256 in the {ds} result line of look {look_id}")
        lr[ds] = m.group(1)
    short = None
    for c in SR.ledger_lines(root):
        if c[3] in ("T1", "T5") and f"look={look_id}" in c[5]:
            short = c[1]
    if not short:
        sys.exit(f"claims_a5: no ledger line for look {look_id}")
    base = subprocess.run(["git", "-C", root, "rev-parse", short], capture_output=True, text=True).stdout.strip()
    if not base:
        sys.exit(f"claims_a5: base commit {short} of look {look_id} not found")
    return lr, base


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1", required=True)
    ap.add_argument("--t5", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default=os.path.join(ROOT, "config", "a5.yaml"))
    ap.add_argument("--ledger-t1-manifest-sha", default=None)
    ap.add_argument("--ledger-t5-manifest-sha", default=None)
    ap.add_argument("--base-commit", default=None)
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--look-id", default=None, help="read the manifest sha256s and the base commit from the ledger "
                    "lines of this look (required for sealed inputs)")
    a = ap.parse_args(argv)
    lr = None
    if a.look_id:
        lr, a.base_commit = ledger_look(a.root, a.look_id)
    elif json.load(open(a.t1)).get("meta", {}).get("role") == "sealed" if os.path.exists(a.t1) else False:
        sys.exit("claims_a5: sealed inputs need --look-id (manifest sha256s and base commit come from the ledger)")
    if a.ledger_t1_manifest_sha or a.ledger_t5_manifest_sha:
        if a.look_id:
            sys.exit("claims_a5: give --look-id or the --ledger-*-manifest-sha options, not both")
        lr = {"T1": a.ledger_t1_manifest_sha, "T5": a.ledger_t5_manifest_sha}
    code, res = run(a.t1, a.t5, a.config, a.out, lr, a.base_commit, a.root)
    print(res.get("failed") or f"claims written: state {res['state']}; see {a.out}/claims.md")
    return code


if __name__ == "__main__":
    sys.exit(main())
