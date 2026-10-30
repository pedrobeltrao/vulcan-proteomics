"""bench/t5_armM.py (Amendment A5, A5.6 / A5.11.3): sealed-mode refusals, dev-mode guards, the combined arm-M estimator
on synthetic winners, schema / completeness, manifest and work clean-up.  No sealed data, no GPU."""
import json
import os
import types

import numpy as np
import pytest

from bench import t5_armM as T5
from vulcan.stats._util import DECOY_SEEDS, ENTRAP_SEEDS_SEALED


def _args(**kw):
    d = dict(mode="sealed", runs=None, entrap_seeds=None, decoy_seeds=None, out=None, work=None, stores=None,
             targets=None, db_recipe=None, workers=None)
    d.update(kw)
    return types.SimpleNamespace(**d)


@pytest.fixture
def look(tmp_path):
    root = tmp_path / "look" / "T5"
    root.mkdir(parents=True)
    env = {"SEALED_DATASET": "T5", "SEALED_FILES": "/x/data/sealed/PXD052597/b.raw /x/data/sealed/PXD052597/a.raw",
           "SEALED_LOOK_ROOT": str(root)}
    return root, env


def test_sealed_ok_defaults(look):
    root, env = look
    st = T5.check_args(_args(out=str(root / "out"), work=str(root / "work")), env)
    assert st["entrap"] == list(ENTRAP_SEEDS_SEALED) and st["decoys"] == list(DECOY_SEEDS)
    assert len(st["entrap"]) * len(st["decoys"]) == 15
    assert st["recipe"] == "shuffle" and st["raws"] == sorted(env["SEALED_FILES"].split())
    assert [os.path.basename(t) for t in st["targets"]] == ["yeast_UP000002311.fasta",
                                                             "contaminants_frankenfield2022.fasta"]


@pytest.mark.parametrize("kw", [
    dict(entrap_seeds=[2001, 2002, 2003]),
    dict(decoy_seeds=[1001]),
    dict(stores="/tmp/x"),
    dict(targets=["a.fasta"]),
    dict(db_recipe="dev_foreign"),
    dict(runs=["/other.raw"]),
])
def test_sealed_refuses_nondefault(look, kw):
    root, env = look
    with pytest.raises(SystemExit):
        T5.check_args(_args(out=str(root / "out"), work=str(root / "work"), **kw), env)


def test_sealed_refuses_outside_root_and_nonempty(look, tmp_path):
    root, env = look
    with pytest.raises(SystemExit):
        T5.check_args(_args(out=str(tmp_path / "elsewhere"), work=str(root / "work")), env)
    (root / "work").mkdir()
    (root / "work" / "f").write_text("x")
    with pytest.raises(SystemExit):
        T5.check_args(_args(out=str(root / "out"), work=str(root / "work")), env)


def test_sealed_needs_sealed_run_env(look):
    root, env = look
    for drop in ("SEALED_DATASET", "SEALED_FILES", "SEALED_LOOK_ROOT"):
        e = {k: v for k, v in env.items() if k != drop}
        with pytest.raises(SystemExit):
            T5.check_args(_args(out=str(root / "out"), work=str(root / "work")), e)
    with pytest.raises(SystemExit):
        T5.check_args(_args(out=str(root / "out"), work=str(root / "work")), {**env, "SEALED_DATASET": "T1"})


def test_dev_guards(tmp_path):
    ok = dict(mode="dev", runs=["data/raw/PXD031616/LM2486_DDA.raw"], out=str(tmp_path / "o"), work=str(tmp_path / "w"))
    st = T5.check_args(_args(**ok, decoy_seeds=[1001]), {})
    assert st["recipe"] == "dev_foreign" and st["entrap"] == [2001, 2002, 2003]
    with pytest.raises(SystemExit):
        T5.check_args(_args(**ok, entrap_seeds=[2004]), {})
    with pytest.raises(SystemExit):
        T5.check_args(_args(**{**ok, "runs": ["data/sealed/PXD052597/a.raw"]}), {})
    with pytest.raises(SystemExit):
        T5.check_args(_args(**ok), {"SEALED_FILES": "x"})


def _winners(rng, n, ent_rate, lists=True):
    is_dec = rng.random(n) < 0.3
    ent = (~is_dec) & (rng.random(n) < ent_rate)
    q = np.sort(rng.random(n) * 0.1)
    x = {"is_decoy": is_dec, "q": q, "ent": ent, "unit": np.arange(n), "seq": np.array([f"S{i // 2}" for i in range(n)]),
         "lists": {"overall": None}}
    if lists:
        k = rng.integers(1, 4, n)
        for s in (1, 2, 3):
            x["lists"][("raw", s)] = k == s
            x["lists"][("merged", s)] = k == s
    return x


def test_arm_m_combined_estimator():
    rng = np.random.default_rng(0)
    reps = {T5.rep_name(j, s): {"decoy_seed": j, "entrap_seed": s} for j in (1001, 1002) for s in (2001, 2002, 2003)}
    per = {r: _winners(rng, 4000, 0.01) for r in reps}
    L = T5.arm_m(per, reps, None, 1.0)
    assert set(L) == {"peptidoform", "sequence"} and set(L["peptidoform"]) == {"0.001", "0.01", "0.05"}
    for a in ("0.001", "0.01", "0.05"):
        rec = L["peptidoform"][a]
        assert set(rec) == {"a", "b", "counts"} and set(rec["a"]) == set(T5.ARM_KEYS)
        # combined estimator with r = 1: FDP_M = 2 N_E / (N_T + N_E), per replicate, mean over replicates
        f = [2 * c["N_E"] / (c["N_T"] + c["N_E"]) for c in rec["counts"].values()]
        assert rec["a"]["mean_fdp"] == pytest.approx(np.mean(f), rel=1e-12)
        assert rec["a"]["mean_lb"] == pytest.approx(np.mean(f) / 2, rel=1e-12)
        assert rec["a"]["verdict"] in ("supported", "supported-noisy", "violated", "not-evaluable", "inconclusive")
    # a gross violation (true FDP ~ 2 x 30%) is 'violated' at 5%
    bad = {r: _winners(rng, 4000, 0.3) for r in reps}
    assert T5.arm_m(bad, reps, None, 1.0)["peptidoform"]["0.05"]["a"]["verdict"] == "violated"
    # list restriction: counts restricted to the stratum mask
    Ls = T5.arm_m(per, reps, ("raw", 1), 1.0)
    r0 = next(iter(reps))
    x = per[r0]
    m = (~x["is_decoy"]) & x["lists"][("raw", 1)] & (x["q"] <= 0.01)
    assert Ls["peptidoform"]["0.01"]["counts"][r0]["N_T"] == int((m & ~x["ent"]).sum())


def test_evaluate_schema_and_completeness(tmp_path, monkeypatch):
    rng = np.random.default_rng(1)
    reps = {T5.rep_name(1001, s): {"index": "i", "decoy_seed": 1001, "entrap_seed": s} for s in (2001, 2002, 2003)}
    dbinfo = {"seeds": {str(s): {"r_peptide_eff": 1.0} for s in (2001, 2002, 2003)}}
    st = {"mode": "dev", "entrap": [2001, 2002, 2003], "decoys": [1001]}
    drop = {"rep": None}

    def closed(rep, info, runs, work, missing):
        if rep == drop["rep"]:
            missing.append(f"v1/{rep}:LM1")
            return None
        return {"v1": _winners(rng, 3000, 0.01, False), "hyperscore_v1": _winners(rng, 3000, 0.01, False)}

    def chim(rep, info, runs, work, missing, mm):
        mm[rep] = {"1": [1], "2": [2, 3], "3": [2, 3]}
        w = _winners(rng, 3000, 0.01)
        w["lists"][("merged", 2)] = w["lists"][("merged", 3)] = w["lists"][("raw", 2)] | w["lists"][("raw", 3)]
        return {"v2_sep": w, "hyperscore_v2": _winners(rng, 3000, 0.01, False)}
    monkeypatch.setattr(T5, "winners_closed", closed)
    monkeypatch.setattr(T5, "winners_chim", chim)
    res = T5.evaluate(st, reps, ["LM1", "LM2"], str(tmp_path), dbinfo, print)
    assert res["complete"] and not res["missing"] and res["n_replicates"] == 3
    for k in ("mode", "runs", "entrap_seeds", "decoy_seeds", "n_replicates", "complete", "missing"):
        assert k in res
    assert set(res["v2_sep"]) == {"overall", "merged", "raw"}
    assert set(res["v2_sep"]["raw"]) == {"1", "2", ">=3"} and set(res["v2_sep"]["merged"]) == {"1", "2_3"}
    for k in ("v1", "hyperscore_v1", "hyperscore_v2"):
        assert set(res[k]) == {"overall"}
    json.dumps(res, default=float)
    assert set(T5.violated(res)) == {"v1", "v2_sep"}
    drop["rep"] = T5.rep_name(1001, 2002)
    res = T5.evaluate(st, reps, ["LM1", "LM2"], str(tmp_path), dbinfo, print)
    assert not res["complete"] and res["missing"] and "v1" not in res


def test_read_runs_flags_missing_run(tmp_path):
    import pandas as pd
    d = tmp_path / "rep"
    d.mkdir()
    pd.DataFrame({"spec_idx": [1], "pf_row": [0]}).to_parquet(d / "A.scored.parquet")
    missing = []
    df = T5._read_runs(str(d), ["A", "B"], ".scored.parquet", missing, "v1/rep")
    assert len(df) == 1 and missing == ["v1/rep:B"]


def test_manifest_and_clean_work(tmp_path):
    out = tmp_path / "out"
    (out / "sub").mkdir(parents=True)
    (out / "armM.json").write_text("{}")
    (out / "sub" / "x").write_text("y")
    m = T5.write_manifest(str(out))
    lines = open(m).read().splitlines()
    assert [l.split("  ")[1] for l in lines] == ["armM.json", "sub/x"]
    assert lines[0].split("  ")[0] == T5._sha(str(out / "armM.json"))
    w = tmp_path / "work"
    for d in ("db", "db_2004", "stage", "ours", "ms1", "predcache_t5_v1"):
        (w / d).mkdir(parents=True)
        (w / d / "f").write_text("z")
    (w / "t5_armM.log").write_text("l")
    T5.clean_work(str(w))
    assert sorted(os.listdir(w)) == ["db", "db_2004"]


def test_sealed_work_in_node_scratch(look, tmp_path):
    root, env = look
    sc = tmp_path / "vulcan_look_x"
    env = {**env, "SEALED_LOOK_SCRATCH": str(sc)}
    T5.check_args(_args(out=str(root / "out"), work=str(sc / "T5" / "work")), env)
    with pytest.raises(SystemExit):                       # --out must be under the look root, not the scratch
        T5.check_args(_args(out=str(sc / "T5" / "out"), work=str(sc / "T5" / "work")), env)
    with pytest.raises(SystemExit):                       # another dataset's scratch
        T5.check_args(_args(out=str(root / "out"), work=str(sc / "T1" / "work")), env)


def test_console_whitelist():
    ok = ["t5_armM: stage db_armc", "t5_armM: stage db_seed_2004", "t5_armM: stage v2_score", "t5_armM: complete",
          "t5_armM: wrote armM.json and manifest.sha256", "t5_armM: stage incomplete"]
    bad = ["t5_armM: stage done complete=True", "t5_armM: scored M_1001_2004: 13.5s", "t5_armM: v1 0.01 violated"]
    assert all(T5.CONSOLE_RE.match(x) for x in ok) and not any(T5.CONSOLE_RE.match(x) for x in bad)


def test_single_bootstrap_equals_two_calls_and_pool():
    """Convention (b) read off the (a) call equals a separate ci_convention='b' call; the fork pool changes nothing."""
    from vulcan.stats import entrapment as ENT
    rng = np.random.default_rng(3)
    reps = {T5.rep_name(j, s): {"decoy_seed": j, "entrap_seed": s} for j in (1001, 1002) for s in (2001, 2002)}
    per = {r: _winners(rng, 3000, 0.02) for r in reps}
    L1 = T5.arm_m(per, reps, None, 1.0, workers=1)
    L4 = T5.arm_m(per, reps, None, 1.0, workers=4)
    assert L1 == L4
    for unit in T5.UNITS:
        for a in T5.ALPHAS:
            rs = [ENT.Replicate(T5._counts(per[r], None, unit, a)[0], reps[r]["entrap_seed"], reps[r]["decoy_seed"])
                  for r in reps]
            for cv in ("a", "b"):
                e = ENT.evaluate_arm(rs, a, "combined", r=1.0, ci_convention=cv)
                got = L1[unit][str(a)][cv]
                for k in T5.ARM_KEYS:
                    assert (tuple(got[k]) if isinstance(got[k], (list, tuple)) else got[k]) == \
                        (tuple(e[k]) if isinstance(e[k], (list, tuple)) else e[k]), (unit, a, cv, k)
