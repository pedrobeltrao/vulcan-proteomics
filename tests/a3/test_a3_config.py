"""config/a3.yaml pinned to the A3 text, to the committed files and to reviews/checks/a3_power.json."""
import json
import os

import pytest

from vulcan.gates import a2_config, a3_config

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CFG = a3_config.load()


def test_text_and_tag():
    assert CFG["amendment_text"] == "reviews/amendment_A3_draft.md" and CFG["tag_required"] == "prereg-v0.4.3"
    assert os.path.exists(os.path.join(ROOT, CFG["amendment_text"]))


def test_model_and_table_hashes():
    got = a3_config.check_models(CFG)
    assert got["v2_sep.model"].startswith("3f2c69cc") and got["v2_sep.preregistration"].startswith("de018e95")
    assert got["ce_s.model"].startswith("aa19c8f3") and got["v1_static_gbm_pred.model"].startswith("31bfe685")
    assert CFG["scorers"]["v1_static_gbm_pred"]["sha256"] == a2_config.load()["gates"]["frozen_scorer"]["sha256"]
    assert CFG["scorers"]["seq_classifier"]["sha256"] == a2_config.load()["gates"]["seq_classifier"]["sha256"]
    assert a3_config.check_offsets(CFG) == CFG["asyn"]["offsets_sha256"]
    fz = CFG["freeze"]
    assert a2_config.sha256_file(os.path.join(ROOT, fz["alphapeptdeep_weights"]["path"])) == fz["alphapeptdeep_weights"]["sha256"]
    assert a2_config.sha256_file(os.path.join(ROOT, fz["sage_binary"]["path"])) == fz["sage_binary"]["sha256"]


def test_asyn_and_rows_match_text():
    a = CFG["asyn"]
    assert a["seed"] == 47 and a["ranks"] == [2, 3] and a["max_tries"] == 32
    assert a["envelope"]["ppm"] == 20.0 and a["envelope"]["k"] == [-1, 0, 1, 2, 3] and a["envelope"]["zmax"] == 6
    assert a["window_screens"]["identity_f_threshold"] == 0.20 and a["in_competition"] is False
    r = CFG["rows"]
    assert r["V1"]["n_min"] == 200 and r["V1"]["margin_abs"] == 0.60 and r["V1"]["margin_diff"] == 0.10
    assert r["V1"]["bootstrap_B"] == 10000 and r["V1"]["bootstrap_seed"] == 83
    assert r["V2"]["n_min_top"] == 7200 and r["V14"]["n_min_top"] == 7200
    assert r["V4"]["n_min"] == 2000 and r["V5"]["n_min"] == 2000 and r["V13"]["n_min"] == 200
    assert r["V8"]["n_min"] == "g2_sample_size" and r["V9"]["alarm"] == 0.55
    assert r["V9"]["parts"]["c"]["min_winners_per_label"] == 200 and r["V11"]["unused"] is True
    t = CFG["tiers"]
    assert t["tier1"] == ["V1", "V2", "V3", "V4", "V6.class1", "V7", "V9.a", "V9.b", "V10", "V12"]
    assert t["tier2"] == ["V5", "V6.class2", "V8", "V9.c", "V13", "V14"] and t["fs_rows"] == ["V10"]


def test_n_min_equal_power_output():
    p = json.load(open(os.path.join(ROOT, CFG["power"]["output"])))
    assert CFG["rows"]["V6"]["n_min_per_label_per_bin"] == p["n_min"]["V6"]["n_min"]
    assert CFG["rows"]["V12"]["n_min_per_label_per_bin"] == p["n_min"]["V12"]["n_min"]
    assert p["seed"] == CFG["power"]["seed"] == 131
    dev = json.load(open(os.path.join(ROOT, CFG["rows"]["V12"]["dev_record"])))
    assert abs(dev["max_variance_decile_miss_rate"] - CFG["rows"]["V12"]["dev_miss_rate"]) < 5e-5


def test_freeze_not_recorded_refuses():
    if CFG["freeze"]["commit"] != a3_config.TODO:
        pytest.skip("freeze recorded")
    st = a3_config.freeze_status(CFG)
    assert st["ok"] is False and "TODO" in st["reason"]
    with pytest.raises(SystemExit):
        a3_config.require_freeze_commit(CFG)


def test_freeze_status_detects_wrong_commit():
    import copy
    c = copy.deepcopy(CFG)
    c["freeze"]["commit"] = "prereg-v0.4.2"
    c["freeze"]["tree_hashes"] = {t: "0" * 40 for t in a3_config.FROZEN_TREES}   # A3.1.3 + A5.11.6 tree list
    st = a3_config.freeze_status(c)
    assert st["ok"] is False and st["tree_hashes_match"] is False


def test_row_tol_conf_equal_statistics_constants():
    from vulcan.gates import a3_stats as S
    for k in ("V6", "V12"):
        assert CFG["rows"][k]["tol"] == S.TOL and CFG["rows"][k]["conf"] == S.CONF
