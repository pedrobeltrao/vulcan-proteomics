"""gates_hn1 (step 5) on synthetic per-null-spectrum arrays: the evaluation path end to end through a6_stats, and the
guards refuse before the tag.  No data, no labels file."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "a6"))
import gates_hn1 as G  # noqa: E402


def synth(n=6000, J=12, seed=1):
    rng = np.random.default_rng(seed)
    proj = np.array([f"PXD{rng.integers(J * 2):06d}" for _ in range(n)])
    clus = np.array([str(hash(p) % J) for p in proj])
    pool = np.where(rng.random(n) < 0.8, "A", "FS")

    def pipe():
        s = rng.normal(size=n)
        return {"win_t": rng.random(n) < 0.5, "score": s, "tail1": s > 1.6, "logit": rng.normal(size=n)}
    g3 = {k: rng.normal(size=n) for k in ("bt", "bd", "kt", "kd")}
    return {"pool": pool, "project": proj, "cluster": clus, "instrument": np.full(n, "QE HF"),
            "run": np.full(n, "r"), "fixed": pipe(), "v1": pipe(), "v1_noRT": pipe(), "g3": g3}


def test_evaluate_returns_closed_status():
    D = synth()
    res = G.evaluate(D, np.random.default_rng(2).random(len(D["pool"])) < 0.9)
    st = res["v1"]["verdict"]["status"]
    assert st in ("PASS", "FAILED", "not evaluable", "invalid pool")
    # n far below G1 n_min / tail n_min -> cannot pass
    assert st != "PASS"
    assert set(res["reported"]) == {"noRT_capped", "noRT_uncapped"}


def test_guards_refuse_without_tag(tmp_path):
    if os.popen(f"git -C {G.ROOT} tag -l {G.TAG}").read().strip():
        pytest.skip("tag exists")
    with pytest.raises(SystemExit):
        G.guards(str(tmp_path / "rec.json"))


def test_guards_refuse_after_hard_stop(tmp_path):
    import datetime as dt
    with pytest.raises(SystemExit, match="hard stop"):
        G.guards(str(tmp_path / "rec.json"), now=G.HARD_STOP + dt.timedelta(seconds=1))
