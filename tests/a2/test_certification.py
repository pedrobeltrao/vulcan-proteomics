"""A2.2(a) certification (a2_admission.certify_run) and A2.3 step 3 marker screen on a synthetic own search."""
import numpy as np

from vulcan.gates import a2_admission as AD
from vulcan.stats._util import hash_tiebreak


def own_search(seed=0, n_spec=4000):
    rng = np.random.default_rng(seed)
    spec, unit, pair, score, dec, pep = [], [], [], [], [], []
    for s in range(n_spec):
        good = rng.random() < 0.6
        for k in range(3):
            u = int(rng.integers(0, 20000))
            d = bool(rng.random() < 0.5) if not (good and k == 0) else False
            spec.append(s)
            unit.append(2 * u + d)
            pair.append(u)
            score.append(rng.gumbel() + (12.0 if (good and k == 0) else 0.0))
            dec.append(d)
            pep.append(f"PEP{u % 997}K" if u % 50 else f"CONTAM{u % 7}K")
    return (np.array(spec), np.array(unit, np.uint64), np.array(pair), np.array(score), np.array(dec),
            np.array(pep, dtype=object))


def test_certify_run_rules():
    spec, unit, pair, score, dec, pep = own_search()
    sources = {p: ({"CRAP|x"} if p.startswith("CONTAM") else {"OWN|P1"}) for p in set(pep)}
    sources["PEP5K"] = {"OWN|P1", "PNL|H|DSP|9606|P15924"}           # also in a panel -> must not certify
    near = lambda ps: np.array([p.startswith("CONTAM") or p == "PEP5K" for p in ps])  # noqa: E731
    tb = hash_tiebreak(unit, 11)
    r = AD.certify_run(spec, unit, pair, score, dec, pep, sources, near, tb)
    rows = r["rows"]
    assert r["n_certified"] == rows.size > 1000
    assert not dec[rows].any()                                          # targets only
    assert not any(pep[i].startswith("CONTAM") or pep[i] == "PEP5K" for i in rows)
    # rank 1: no other candidate of the same spectrum scores higher
    for i in rows[:200]:
        assert score[i] >= score[spec == spec[i]].max()
    # stricter q never certifies more
    r2 = AD.certify_run(spec, unit, pair, score, dec, pep, sources, near, tb, q_psm=0.0001)
    assert r2["n_certified"] <= r["n_certified"]


def test_marker_screen():
    sources = {"AAAK": {"PNL|H|DSP|9606|P1"}, "CCCK": {"PNL|H|DSP|9606|P1"}, "DDDK": {"PNL|H|DCD|9606|P2"},
               "EEEK": {"PNL|H|DCD|9606|P2", "PNL|B|ALB|9913|P3"}, "FFFK": {"PNL|B|ALB|9913|P3"},
               "GGGK": {"PNL|K|KRT1|9606|P4"}, "HHHK": {"PNL|K|KRT1|9606|P4"}, "OWNK": {"OWN|Q1", "PNL|H|LYZ|9606|P5"},
               "LYZK": {"PNL|H|LYZ|9606|P5"}}
    near_own = lambda ps: np.array([p == "OWNK" for p in ps])  # noqa: E731
    acc = {"r1": {"AAAK", "DDDK", "GGGK", "OWNK", "LYZK"}, "r2": {"CCCK", "EEEK", "HHHK"}, "r3": {"FFFK"}}
    m = AD.marker_screen(acc, sources, near_own)
    # union over runs: DSP {AAAK, CCCK}, DCD {DDDK, EEEK}, ALB {EEEK, FFFK} (shared EEEK counts for both),
    # KRT1 {GGGK, HHHK} detected but never triggers; LYZ has 1 peptide after the own-near OWNK is dropped
    assert m["detected_members"] == ["B|ALB|9913", "H|DCD|9606", "H|DSP|9606", "K|KRT1|9606"]
    assert m["n_trigger_members_detected"] == 3 and m["exclude"]
    m1 = AD.marker_screen({"r1": {"GGGK", "HHHK", "AAAK", "CCCK"}}, sources, near_own)
    assert m1["n_trigger_members_detected"] == 1 and not m1["exclude"]


def test_stopping_state():
    cfg = {"min_projects": 10, "min_A": 120000, "min_A_with_candidate": 72000, "min_tail_proxy": 400,
           "cap_projects": 16}
    p = {"n_A": 12000, "n_A_with_candidate": 8000, "tail_proxy": 40}
    assert not AD.stopping_state([p] * 9, cfg)["stop"]
    assert AD.stopping_state([p] * 10, cfg)["reason"] == "criteria met"
    q = {"n_A": 12000, "n_A_with_candidate": 8000, "tail_proxy": 10}
    assert not AD.stopping_state([q] * 15, cfg)["stop"]
    assert AD.stopping_state([q] * 16, cfg)["reason"] == "cap reached"
