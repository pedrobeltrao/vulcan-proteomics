"""MC power checks of the implemented A2 gate code vs reviews/checks/a2_power.json (synthetic only).
Quick version here (reduced replicates, 4-SE agreement); the full run is scripts/a2/power_checks.py
(results/2026-10-27_A2_implementation/power_checks.json) and the --runslow test below."""
import pytest

import power_checks as PC


def test_power_quick():
    r = PC.run(scale=0.15)
    bad = [x for x in r["S3"] + r["S5"]["table"] + r["S6"] + r["S7"] if not x["agree"]]
    assert not bad, bad
    assert r["S5"]["n_min"]["config_ge_draft"]
    assert r["tail_rule"]["a2_never_passes_where_a1_fails"]
    # A1 margins kept: a +0.10 tail bias fails (>= 0.95 at this replicate count) from n_tail 400
    for g in r["tail_rule"]["grid"]:
        if g["n_tail_mean"] >= 400 and g["p_learned"] >= 0.60:
            assert g["pass"] <= 0.10, g


@pytest.mark.slow
def test_power_full():
    r = PC.run(scale=1.0)
    assert r["all_agree"], [x for x in r["S3"] + r["S5"]["table"] + r["S6"] if not x["agree"]]
