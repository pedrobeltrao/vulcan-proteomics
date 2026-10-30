"""bench/t1/claims_a5.py (DESIGN §6.0d, A5.11.4.4): every branch on synthetic inputs, the A6 status-file hook in a
temp git repo, incomplete inputs, the wording lint on the generated claims.md, and an end-to-end run on the committed
D1 dry-run metrics.

Synthetic inputs are built from the committed D1 dry-run metrics (results/2026-10-28_A3_freeze/D1/with_a3_metrics.json,
real driver structure) with: 6 fake run names and 6-entry LOO lists (loo_min set per test), a synthetic `a5` block
(bench/t1/a5.py schema) and a synthetic T5 armM.json (bench/t5_armM.py schema).  D1 itself has one run, so its
LOO minimum is undefined; the end-to-end test checks that this yields "not beaten" (never a claim by omission)."""
import copy
import json
import os
import subprocess

import pytest

from bench.t1 import a3 as A3
from bench.t1 import claims_a5 as C
from bench.t1 import wording_lint as WL

ROOT = C.ROOT
D1 = os.path.join(ROOT, "results/2026-10-28_A3_freeze/D1/with_a3_metrics.json")
CFG = os.path.join(ROOT, "config/a5.yaml")
RUNS6 = [f"run{i}" for i in range(6)]
LEVELS = ("peptidoform", "sequence")


def _L(verdict="supported", alphas=("0.001", "0.01", "0.05"), units=LEVELS):
    rec = lambda: {"verdict": verdict, "mean_fdp": 0.009, "ci95": [0.007, 0.011], "mean_lb": 0.005,
                   "lb_ci95": [0.004, 0.006], "mean_n_t": 5000.0}
    return {u: {a: {"a": rec(), "b": rec()} for a in alphas} for u in units}


def make_armm(mode="dev", runs=("d7a",), n=3):
    return {"mode": mode, "runs": list(runs), "entrap_seeds": [2001, 2002, 2003], "decoy_seeds": [1001],
            "n_replicates": n, "complete": True, "missing": [],
            "v1": {"overall": _L()}, "v2_sep": {"overall": _L(), "merged": {"{1}": _L(), "{2, >=3}": _L()},
                                                "raw": {"1": _L(), "2": _L(), ">=3": _L()}},
            "hyperscore_v1": {"overall": _L()}, "hyperscore_v2": {"overall": _L()}}


def _armc(lo=0.003):
    return {"ours_q_at_match": 0.011, "fdp_reached": True, "lb_reps": [0.004] * 5, "mean_lb": 0.004,
            "lb_ci95": [lo, 0.005], "lb_ci95_b": [lo, 0.005], "mean_fdp_c": 0.008, "r": 13.5}


def make_metrics(loo_scale=1.0):
    m = json.load(open(D1))
    m["runs"] = list(RUNS6)
    m["meta"]["runs"] = list(RUNS6)
    for num, byc in m["a3_secondary"]["loo"].items():
        for comp, byl in byc.items():
            for lv, d in byl.items():
                for kind, byt in d.items():
                    for t, v in byt.items():
                        r = v["full"]
                        v["loo"] = [None if r is None else r * loo_scale] * 6
                        v["loo_min"] = None if r is None else r * loo_scale
    meths = A3.A3_METHODS
    m["a5"] = {"version": "a5-1",
               "stratum1": {"noninferiority": copy.deepcopy(m["noninferiority_by_numerator"]["vulcan_v2_sep"]),
                            "loo": copy.deepcopy(m["a3_secondary"]["loo"]["vulcan_v2_sep"]), "ids": {}},
               "armC_matched": {x: {lv: {f: _armc() for f in ("0.01", "0.05")} for lv in LEVELS} for x in meths},
               "rank_change": {"vulcan": {"frac_spectra_rank_gt1": 0.12, "frac_accepted_units_rank_gt1": 0.03},
                               "vulcan_v2_sep": {"frac_accepted_units_rank_gt1": 0.05}},
               "screened_armC": {"n_excluded": 17, "methods": {"vulcan": {"peptidoform": {"0.01": {
                   "a": {"verdict": "supported", "mean_fdp": 0.008, "ci95": [0.006, 0.01], "mean_lb": 0.007,
                         "lb_ci95": [0.006, 0.008]}, "b": {"verdict": "supported"}}}}}},
               "complete": True, "missing": []}
    return m


def set_ratio(m, num, comp, lv, kind, thr, ratio=None, lb=None, loo=None, **kw):
    r = m["noninferiority_by_numerator"][num][comp][lv][kind][thr]
    if ratio is not None:
        r["ratio"] = ratio
    if lb is not None:
        r["ci95"] = [lb, max(lb, r["ci95"][1])]
    r.update(kw)
    if loo is not None:
        e = m["a3_secondary"]["loo"][num][comp][lv][kind][thr]
        e["loo_min"] = loo
        e["loo"] = [loo] * 6


def all_ratios(m, num, ratio, lb, loo):
    for comp, byl in m["noninferiority_by_numerator"][num].items():
        for lv, d in byl.items():
            for kind, byt in d.items():
                for t in byt:
                    kw = {"fdp_reached": [True, True], "n_boot_nan": 0} if kind == "matched_fdp" else {}
                    set_ratio(m, num, comp, lv, kind, t, ratio, lb, loo, **kw)


def set_verdict(m, method, lv, al, arm, v, conv="a"):
    m["methods"][method]["verdicts"][lv][al][arm][conv]["verdict"] = v


@pytest.fixture
def run(tmp_path):
    def _run(metrics, armm, base_commit=None, root=ROOT, cfg=CFG):
        t1 = tmp_path / "t1"
        t5 = tmp_path / "t5"
        t1.mkdir(exist_ok=True)
        t5.mkdir(exist_ok=True)
        json.dump(metrics, open(t1 / "metrics.json", "w"))
        json.dump(armm, open(t5 / "armM.json", "w"))
        out = tmp_path / "out"
        code, res = C.run(str(t1 / "metrics.json"), str(t5 / "armM.json"), cfg, str(out), None, base_commit, root)
        return code, res, open(out / "claims.md").read()
    return _run


def lint_ok(res, md):
    labels = {"rendered": list(res["labels_rendered"].values()), "L-S": [res["labels_rendered"]["L-S"]]}
    v = WL.lint(md, "claims", labels)
    assert not v, v[:10]


# ------------------------------------------------------------------------------------------------ basics
def test_comparator_set_is_a3s():
    assert C.A3_COMPARATORS == A3.A3_COMPARATORS
    assert C.A6_ENUM == tuple(C.load_yaml(CFG)["a6"]["enum"])
    assert C.A6_STATUS_PATH == C.load_yaml(CFG)["a6"]["status_path"]


def test_all_pass_branch_and_lint(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    code, res, md = run(m, make_armm())
    assert code == 0
    assert res["P1"]["outcome"] == "met"
    assert [s["outcome"] for s in res["S"]][:4] == ["rejected (claim holds)"] * 4
    assert not res["voids"]["v1_void"] and not res["voids"]["v2_sep_void"]
    for rid in ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "S1", "S2", "S3", "S4", "S5", "R1", "R2", "R3",
                "R4", "D", "G"]:
        assert any(r["id"] == rid for r in res["rows"]), rid
    assert "FAILED" in res["labels_rendered"]["L-S"]
    assert "Secondary arm; not gate-validated" in md
    lint_ok(res, md)


def test_deterministic(run, tmp_path):
    m = make_metrics()
    _, _, md1 = run(m, make_armm())
    j1 = open(tmp_path / "out" / "claims.json").read()
    _, _, md2 = run(m, make_armm())
    assert md1 == md2 and j1 == open(tmp_path / "out" / "claims.json").read()


# ------------------------------------------------------------------------------------------------ v1 voids
@pytest.mark.parametrize("lv,arm", [("peptidoform", "P"), ("sequence", "P"), ("protein", "P"),
                                    ("peptidoform", "C"), ("sequence", "C")])
@pytest.mark.parametrize("al", ["0.01", "0.05"])
def test_v1_void_t1(run, lv, arm, al):
    m = make_metrics()
    set_verdict(m, "vulcan", lv, al, arm, "violated")
    code, res, md = run(m, make_armm())
    assert res["voids"]["v1_void"] and res["P1"]["outcome"] == "VOID"
    for t in ("P1(i)", "P1(ii)", "P2", "P3"):
        assert f"**Table {t} — VOID under §6.0d A5.7 (T1 arm {arm} violated" in md, t
    assert "is VOID under §6.0d A5.7 (T1 arm" in md.split("## Headline", 1)[1].split("## Abstract")[0]
    lint_ok(res, md)


def test_v1_void_convention_b_does_not_void(run):
    m = make_metrics()
    set_verdict(m, "vulcan", "peptidoform", "0.01", "P", "violated", conv="b")
    _, res, _ = run(m, make_armm())
    assert not res["voids"]["v1_void"]


def test_v1_void_at_0001_does_not_void(run):
    m = make_metrics()
    set_verdict(m, "vulcan", "peptidoform", "0.001", "P", "violated")
    _, res, _ = run(m, make_armm())
    assert not res["voids"]["v1_void"]


@pytest.mark.parametrize("al", ["0.01", "0.05"])
def test_v1_void_t5_armm(run, al):
    a = make_armm()
    a["v1"]["overall"]["peptidoform"][al]["a"]["verdict"] = "violated"
    _, res, _ = run(make_metrics(), a)
    assert res["voids"]["v1_void"]
    assert not res["voids"]["v2_sep_void"]


def test_hyperscore_violation_never_voids(run):
    a = make_armm()
    a["hyperscore_v1"]["overall"]["peptidoform"]["0.01"]["a"]["verdict"] = "violated"
    _, res, _ = run(make_metrics(), a)
    assert not res["voids"]["v1_void"]


def test_claim_level_void_boundary(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    m["a5"]["armC_matched"]["vulcan"]["peptidoform"]["0.01"]["lb_ci95"] = [0.01, 0.02]    # == f: not void
    _, res, _ = run(m, make_armm())
    assert not res["voids"]["claim_level"]["vulcan"]["peptidoform/0.01"]["void"]
    assert res["P1"]["outcome"] == "met"
    m["a5"]["armC_matched"]["vulcan"]["peptidoform"]["0.01"]["lb_ci95"] = [0.0100001, 0.02]   # > f: void
    _, res, md = run(m, make_armm())
    assert res["voids"]["claim_level"]["vulcan"]["peptidoform/0.01"]["void"]
    assert res["P1"]["outcome"] == "VOID (claim-level, A5.7.1: peptidoform/1%)"
    assert not res["voids"]["v1_void"]          # claim-level only
    head = md.split("## Headline\n\n", 1)[1].split("\n")[0]
    assert "is VOID under §6.0d A5.7 (claim-level" in head and "(VOID under §6.0d A5.7; printed in full" in head
    assert "**Table P1(i) — VOID under §6.0d A5.7 (claim-level" in md
    assert "**Table P2 — VOID under §6.0d A5.7 (claim-level" in md
    lint_ok(res, md)


def test_comparator_armc_never_voids(run):
    m = make_metrics()
    m["a5"]["armC_matched"]["sage"]["peptidoform"]["0.01"]["lb_ci95"] = [0.5, 0.6]
    _, res, _ = run(m, make_armm())
    assert not any(d["void"] for d in res["voids"]["claim_level"]["vulcan"].values())


# ------------------------------------------------------------------------------------------------ v2_sep voids / S
def test_v2_void_stratum(run):
    m = make_metrics()
    k = sorted(m["a3_secondary"]["strata"]["merged"])[-1]
    m["a3_secondary"]["strata"]["merged"][k]["verdicts"]["0.05"]["C"]["a"]["verdict"] = "violated"
    _, res, md = run(m, make_armm())
    assert res["voids"]["v2_sep_void"] and not res["voids"]["v1_void"]
    assert all(s["outcome"] == "VOID" for s in res["S"])
    assert all(s["rows"] for s in res["S"] if s["id"] != "S5" or s["rows"])          # numbers kept (A5.7.3)
    assert len(res["S"][0]["rows"]) >= len(C.A3_COMPARATORS)
    for t in ("S", "R3", "R3-strata"):
        assert f"**Table {t} — VOID under §6.0d A5.7 (" in md, t
    lint_ok(res, md)


def test_v2_void_t5_merged(run):
    a = make_armm()
    a["v2_sep"]["merged"]["{2, >=3}"]["peptidoform"]["0.01"]["a"]["verdict"] = "violated"
    _, res, _ = run(make_metrics(), a)
    assert res["voids"]["v2_sep_void"]


def test_v2_raw_stratum_t5_does_not_void(run):
    a = make_armm()
    a["v2_sep"]["raw"][">=3"]["peptidoform"]["0.01"]["a"]["verdict"] = "violated"
    _, res, _ = run(make_metrics(), a)
    assert not res["voids"]["v2_sep_void"]


def test_v2_claim_level_void(run):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    m["a5"]["armC_matched"]["vulcan_v2_sep"]["peptidoform"]["0.01"]["lb_ci95"] = [0.02, 0.03]
    _, res, md = run(m, make_armm())
    assert res["S"][0]["outcome"].startswith("VOID (claim-level, A5.7.2: peptidoform/0.01")
    assert all(not s["tested"] for s in res["S"][1:])
    assert "**Table S — VOID under §6.0d A5.7 (S1 claim-level" in md
    lint_ok(res, md)


@pytest.mark.parametrize("stop", [0, 1, 2, 3, 4])
def test_s_stops_at_each(run, stop):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    spec = C.S_SPEC[stop]
    lv, kind, thr = spec[2][0][:3]
    comp = "sage" if spec[2][0][5] == "all" else next(c for c in C.A3_COMPARATORS if C.own_q_supported(m, c))
    set_ratio(m, "vulcan_v2_sep", comp, lv, kind, thr, lb=0.98 if stop in (1, 4) else 0.5)   # S2/S5 share S1's rows
    _, res, _ = run(m, make_armm())
    outs = [s["outcome"] for s in res["S"]]
    assert outs[:stop] == ["rejected (claim holds)"] * stop
    assert outs[stop] == "not rejected"
    assert all(o.startswith("not tested") for o in outs[stop + 1:])


def test_s_loo_backstop(run):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    set_ratio(m, "vulcan_v2_sep", "comet+mokapot", "peptidoform", "matched_fdp", "0.01", loo=0.99)
    _, res, _ = run(m, make_armm())
    assert res["S"][0]["outcome"] == "rejected (claim holds)"      # S1 bar 0.97
    assert res["S"][1]["outcome"] == "not rejected"                 # S2 bar > 1.00 fails on LOO min


def test_s3_skip_on_protein_violated(run):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    set_verdict(m, "vulcan_v2_sep", "protein", "0.01", "P", "violated")
    _, res, md = run(m, make_armm())
    assert not res["voids"]["v2_sep_void"]
    outs = [s["outcome"] for s in res["S"]]
    assert outs[:2] == ["rejected (claim holds)"] * 2
    assert "violated" in outs[2] and outs[3].startswith("not tested") and outs[4].startswith("not tested")


def test_not_evaluable_comparator_is_not_beaten(run):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    set_ratio(m, "vulcan_v2_sep", "sage", "peptidoform", "matched_fdp", "0.01", n_boot_nan=101)
    set_ratio(m, "vulcan", "sage", "peptidoform", "matched_fdp", "0.01", n_boot_nan=101)
    _, res, _ = run(m, make_armm())
    assert res["S"][0]["outcome"] == "not rejected"
    assert res["P1"]["outcome"] == "not met"
    row = next(r for r in res["P1"]["i"] if r["comparator"] == "sage")
    assert not row["evaluable"] and "not beaten" in row["why"]
    set_ratio(m, "vulcan", "sage", "peptidoform", "matched_fdp", "0.01", n_boot_nan=100)    # exactly 5%: evaluable
    _, res, _ = run(m, make_armm())
    assert res["P1"]["outcome"] == "met"


def test_fdp_not_reached_is_not_evaluable(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    set_ratio(m, "vulcan", "sage+mokapot", "peptidoform", "matched_fdp", "0.01", fdp_reached=[True, False])
    _, res, _ = run(m, make_armm())
    assert res["P1"]["outcome"] == "not met"


def test_no_supported_comparator(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    for c in C.A3_COMPARATORS:
        m["methods"][c]["verdicts"]["peptidoform"]["0.01"]["deciding"] = "inconclusive"
    _, res, md = run(m, make_armm())
    assert res["P1"]["ii_supported"] == [] and res["P1"]["outcome"] == "met"
    assert "no supported comparator" in md
    assert res["S"][0]["outcome"] == "not rejected"           # stricter: S1's own-q part cannot be met (M1)
    assert all(s["outcome"].startswith("not tested") for s in res["S"][1:])
    lint_ok(res, md)


def test_headline_lowest_ratio_and_loss(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    set_ratio(m, "vulcan", "sage+pred+mokapot", "peptidoform", "matched_fdp", "0.01", ratio=0.99, lb=0.96, loo=0.98)
    _, res, md = run(m, make_armm())
    assert res["P1"]["headline_comparator"] == "sage+pred+mokapot"
    assert res["P1"]["outcome"] == "not met"
    head = md.split("## Headline\n\n", 1)[1].split("\n")[0]
    assert head.startswith("[L-S] ") and "sage+pred+mokapot" in head and "0.99" in head


def test_not_evaluable_stratum_labelled(run):
    m = make_metrics()
    k = sorted(m["a3_secondary"]["strata"]["merged"])[-1]
    m["a3_secondary"]["strata"]["merged"][k]["verdicts"]["0.01"]["P"]["a"]["verdict"] = "not-evaluable"
    _, res, md = run(m, make_armm())
    assert "not evaluable at α; no stratum-specific statement" in md
    assert not res["voids"]["v2_sep_void"]


# ------------------------------------------------------------------------------------------------ A5.9.3 states (A6 hook)
def _repo(tmp_path, status, date="2026-10-28T17:00:00+01:00"):
    r = tmp_path / "repo"
    r.mkdir()
    env = {**os.environ, "GIT_COMMITTER_DATE": date, "GIT_AUTHOR_DATE": date}
    g = lambda *a: subprocess.run(["git", "-C", str(r), *a], check=True, capture_output=True, text=True,
                                  env=env).stdout.strip()
    g("init", "-q")
    g("config", "user.email", "t@t")
    g("config", "user.name", "t")
    (r / "README").write_text("x")
    g("add", "README")
    g("commit", "-qm", "init")
    if status is not None:
        p = r / C.A6_STATUS_PATH
        p.parent.mkdir(parents=True)
        out = p.parent / "gates.json"
        out.write_text("{}")
        g("add", str(out))
        g("commit", "-qm", "step5")
        c5 = g("rev-parse", "HEAD")
        import hashlib
        d = {"status": status, "step5_commit": c5,
             "gate_output_sha256": hashlib.sha256(b"{}").hexdigest(),
             "gate_output_path": "results/2026-10-28_A6_hn1/gates.json"}
        if isinstance(status, dict):
            d = status
        p.write_text(json.dumps(d))
        g("add", str(p))
        g("commit", "-qm", "status")
    return str(r), g("rev-parse", "HEAD")


@pytest.mark.parametrize("status,state", [
    ("PASS", "GO, A6 PASS"), ("FAILED", "GO, A6 FAILED"), ("not evaluable", "GO, A6 other"),
    ("invalid pool", "GO, A6 other"), ("not completed", "GO, A6 other"), ("label-free stop", "GO, A6 other"),
    ("not tagged", "GO, A6 other"), ("stopped: node time limit", "GO, A6 other"), ("pass", "GO, A6 other"),
    (None, "GO, A6 other")])
def test_a6_states(run, tmp_path, status, state):
    root, base = _repo(tmp_path, status)
    a6 = C.a6_state(root, base, C.load_yaml(CFG))
    if status is None or status == "pass":
        assert a6["status"] == "stopped: node time limit"
    assert C.a59_state(C.load_yaml(CFG), a6) == state
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    code, res, md = run(m, make_armm(), base_commit=base, root=root)
    assert code == 0 and res["state"] == state
    ls = res["labels_rendered"]["L-S"]
    assert "FAILED" in ls
    if state == "GO, A6 FAILED":
        assert res["voids"]["v1_void"] and ls.endswith("A1 FAILED; A6 FAILED.")
        assert res["P1"]["outcome"] == "VOID"
    elif state == "GO, A6 PASS":
        assert "passed null-spectrum gates at its second completed test (A6, α/2) before T1" in ls
        assert "A1.8's precondition met at test 2 (A6)" in res["lr_v1"]
        assert "validated by A6 (test 2), after A1 FAILED" in res["labels_rendered"]["L-Q_v1"]
    else:
        assert ls == C.load_yaml(CFG)["labels"]["L-S"]["text"]
        st = a6["status"]
        assert f"A6 (removed-peptide human nulls, test 2 at α/2): {st} (commit" in res["labels_rendered"]["L-G_v1"]
    assert "A6 gates v1 only (D32)" in res["labels_rendered"]["L-G_v2_sep"]
    # B1: v2_sep never carries A6 wording (A6 gates v1 only, D32)
    assert res["labels_rendered"]["L-Q_v2_sep"] == C.load_yaml(CFG)["labels"]["L-Q"]["text"]
    assert "precondition met" not in res["lr_v2"] and "is used although" in res["lr_v2"]
    for t in ("**Table S", "**Table R3:", "**Table R3-strata", "**Table R1"):
        cap = md.split(t, 1)[1].split("\n")[0]
        assert "validated by A6" not in cap and "precondition met" not in cap, t
    if state == "GO, A6 PASS":
        for t in ("**Table P1(ii)", "**Table P2", "**Table P3"):
            if t in md:
                assert "validated by A6" in md.split(t, 1)[1].split("\n")[0], t
    lint_ok(res, md)


def test_a6_not_committed_at_base(tmp_path):
    root, base = _repo(tmp_path, None)
    p = os.path.join(root, C.A6_STATUS_PATH)
    os.makedirs(os.path.dirname(p))
    open(p, "w").write(json.dumps({"status": "PASS"}))        # present on disk, not committed
    assert C.a6_state(root, base, None)["status"] == "stopped: node time limit"


def test_a6_pass_unverifiable(tmp_path):
    root, base = _repo(tmp_path, {"status": "PASS", "step5_commit": None, "gate_output_sha256": None})
    assert C.a6_state(root, base, None)["status"] == "stopped: node time limit"


def test_a6_failed_unverifiable_stays_failed(tmp_path):
    root, base = _repo(tmp_path, {"status": "FAILED", "step5_commit": None, "gate_output_sha256": None})
    assert C.a6_state(root, base, None)["status"] == "FAILED"


def test_a6_no_base_commit():
    assert C.a6_state(ROOT, None, None)["status"] == "stopped: node time limit"


def test_nogo_state_wording(run, tmp_path):
    import yaml
    cfg = C.load_yaml(CFG)
    cfg["c9"]["verdict"] = "NO-GO"
    cfg["c9"]["p_v1_minimum"] = 0.5
    p = tmp_path / "a5.yaml"
    p.write_text(yaml.safe_dump(cfg, allow_unicode=True))
    _, res, md = run(make_metrics(), make_armm(), cfg=str(p))
    assert res["state"] == "NO-GO"
    assert "C9 removed-peptide pilot: NO-GO (P(v1 minimum) = 0.5 < 0.8, commit 3097ef5); not run." in \
        res["labels_rendered"]["L-G_v1"]
    lint_ok(res, md)


# ------------------------------------------------------------------------------------------------ incomplete / hashes
@pytest.mark.parametrize("breaker", [
    lambda m, a: m["methods"].__setitem__("sage", {"missing_replicates": ["P_1001_2001"]}),
    lambda m, a: m.pop("a5"),
    lambda m, a: m["a5"].__setitem__("complete", False),
    lambda m, a: m["a5"]["armC_matched"].pop("sage"),
    lambda m, a: m["meta"]["a3"]["methods"].remove("sage_chim2"),
    lambda m, a: m["noninferiority_by_numerator"]["vulcan"]["sage"]["peptidoform"]["matched_fdp"]["0.01"].pop("n_boot_nan"),
    lambda m, a: m["noninferiority_by_numerator"]["vulcan"]["sage"]["peptidoform"]["matched_fdp"]["0.01"].pop("fdp_reached"),
    lambda m, a: a.__setitem__("complete", False),
    lambda m, a: a["v2_sep"].pop("merged"),
    lambda m, a: a["v1"]["overall"]["peptidoform"].pop("0.05"),
    lambda m, a: m["a3_secondary"]["loo"]["vulcan"]["sage"]["peptidoform"]["own_q"]["0.01"].__setitem__("loo", [1.0] * 5),
])
def test_incomplete_refused(run, breaker):
    m, a = make_metrics(), make_armm()
    breaker(m, a)
    code, res, md = run(m, a)
    assert code == 2
    assert md.startswith("look failed: incomplete")


def test_sealed_counts(run):
    m = make_metrics()
    m["meta"]["role"] = "sealed"
    code, _, md = run(m, make_armm())
    assert code == 2 and "T5 runs" in md and "replicates" in md
    code, _, md = run(m, make_armm(mode="sealed", runs=("a", "b"), n=15))
    assert code == 3 and "ledger" in md            # complete, but sealed inputs need the ledger hashes
    m["runs"] = RUNS6[:5]
    code, _, md = run(m, make_armm(mode="sealed", runs=("a", "b"), n=15))
    assert code == 2


def test_ledger_hash_check(tmp_path):
    m = make_metrics()
    t1, t5 = tmp_path / "T1", tmp_path / "T5"
    t1.mkdir()
    t5.mkdir()
    json.dump(m, open(t1 / "metrics.json", "w"))
    json.dump(make_armm(), open(t5 / "armM.json", "w"))
    for d, f in ((t1, "metrics.json"), (t5, "armM.json")):
        (d / "manifest.sha256").write_text(f"{C.sha256_file(str(d / f))}  {f}\n")
    good = {"T1": C.sha256_file(str(t1 / "manifest.sha256")), "T5": C.sha256_file(str(t5 / "manifest.sha256"))}
    code, _ = C.run(str(t1 / "metrics.json"), str(t5 / "armM.json"), CFG, str(tmp_path / "o1"), good)
    assert code == 0
    code, _ = C.run(str(t1 / "metrics.json"), str(t5 / "armM.json"), CFG, str(tmp_path / "o2"), {**good, "T5": "0" * 64})
    assert code == 3
    json.dump({**m, "x": 1}, open(t1 / "metrics.json", "w"))          # metrics edited after the manifest
    code, _ = C.run(str(t1 / "metrics.json"), str(t5 / "armM.json"), CFG, str(tmp_path / "o3"), good)
    assert code == 3


# ------------------------------------------------------------------------------------------------ end to end on D1
def test_end_to_end_d1(run):
    """Committed D1 dry-run metrics (one run) + synthetic a5 block and D7-like armM: must not crash; with one run the
    LOO minimum is undefined, so nothing that needs it can be claimed."""
    m = json.load(open(D1))
    m["a5"] = make_metrics()["a5"]
    code, res, md = run(m, make_armm())
    assert code == 0
    assert res["P1"]["outcome"] == "not met"
    assert all(not r["beaten"] for r in res["P1"]["i"])
    assert res["S"][0]["outcome"] == "not rejected"
    lint_ok(res, md)


def test_a6_pass_after_deadline(tmp_path):
    root, base = _repo(tmp_path, "PASS", date="2026-10-28T17:56:00+01:00")
    a6 = C.a6_state(root, base, C.load_yaml(CFG))
    assert a6["status"] == "stopped: node time limit" and "deadline" in a6["reason"]


def test_p1_ii_fails_while_i_passes(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    sup = [c for c in C.A1_COMPARATORS if C.own_q_supported(m, c)]
    assert sup
    set_ratio(m, "vulcan", sup[0], "peptidoform", "own_q", "0.01", lb=0.9)
    _, res, _ = run(m, make_armm())
    assert res["P1"]["i_met"] and res["P1"]["ii_met"] is False and res["P1"]["outcome"] == "not met"


@pytest.mark.parametrize("al", ["0.01", "0.05"])
def test_v2_void_t5_overall(run, al):
    a = make_armm()
    a["v2_sep"]["overall"]["peptidoform"][al]["a"]["verdict"] = "violated"
    _, res, _ = run(make_metrics(), a)
    assert res["voids"]["v2_sep_void"] and not res["voids"]["v1_void"]


def test_stricter_list_equals_config():
    got = [x.split(": ", 1)[1] for x in C.load_yaml(CFG)["claim_script_stricter_readings"]]
    assert got == C.STRICTER


def test_no_stratum1_decision_arithmetic(run, tmp_path):
    m = make_metrics()
    all_ratios(m, "vulcan_v2_sep", 1.2, 1.1, 1.1)
    _, res, md = run(m, make_armm())
    js = open(tmp_path / "out" / "claims.json").read()
    assert "R1_descriptive" not in js
    r1 = md.split("### R1:", 1)[1].split("### R4")[0]
    assert "rejected" not in r1 and "claim holds" not in r1


def test_headline_ignores_none_ratio(run):
    m = make_metrics()
    all_ratios(m, "vulcan", 1.2, 1.1, 1.1)
    m["noninferiority_by_numerator"]["vulcan"]["sage"]["peptidoform"]["matched_fdp"]["0.01"]["ratio"] = None
    set_ratio(m, "vulcan", "comet+mokapot", "peptidoform", "matched_fdp", "0.01", ratio=1.05)
    _, res, _ = run(m, make_armm())
    assert res["P1"]["headline_comparator"] == "comet+mokapot"


def test_sealed_requires_ledger_hashes(tmp_path):
    m = make_metrics()
    m["meta"]["role"] = "sealed"
    t1, t5 = tmp_path / "a", tmp_path / "b"
    t1.mkdir()
    t5.mkdir()
    json.dump(m, open(t1 / "metrics.json", "w"))
    json.dump(make_armm(mode="sealed", runs=("a", "b"), n=15), open(t5 / "armM.json", "w"))
    code, _ = C.run(str(t1 / "metrics.json"), str(t5 / "armM.json"), CFG, str(tmp_path / "o"))
    assert code == 3


def test_end_to_end_d0(run):
    """Committed D0 dry-run metrics (2 runs, dev) + synthetic a5 block and armM (A5.11.4.4)."""
    m = json.load(open(os.path.join(ROOT, "results/2026-10-28_A3_freeze/D0/with_a3_metrics.json")))
    m["a5"] = make_metrics()["a5"]
    code, res, md = run(m, make_armm())
    assert code == 0, md[:500]
    assert len(m["runs"]) == 2
    assert {r["id"] for r in res["rows"]} >= {"P1", "S1", "R1", "G"}
    lint_ok(res, md)


def _dev_armm_files():
    import glob
    out = []
    for p in sorted(glob.glob(os.path.join(ROOT, "results", "**", "*armM.json"), recursive=True)):
        try:
            d = json.load(open(p))
        except ValueError:
            continue
        if isinstance(d, dict) and d.get("mode") == "dev" and "v1" in d and "v2_sep" in d:
            if subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", p], capture_output=True).returncode == 0:
                out.append(p)
    return out


def test_dev_armm_output_committed():
    """Review R1: the A5.6 dev arm-M output of bench/t5_armM.py must be committed (the end-to-end test is not skipped)."""
    assert _dev_armm_files(), "no committed bench/t5_armM.py dev armM.json"


@pytest.mark.parametrize("path", _dev_armm_files() or ["none"])
def test_end_to_end_d7_armm(run, path):
    """End to end on the committed D7 dev arm-M output of bench/t5_armM.py (A5.11.4.4) with the D1 metrics."""
    m = make_metrics()
    code, res, md = run(m, json.load(open(path)))
    assert code == 0, md[:500]
    lint_ok(res, md)
