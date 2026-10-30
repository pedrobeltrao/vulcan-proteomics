"""config/a2.yaml pins every A2 threshold and seed to the draft text (reviews/amendment_A2_draft.md) and the A1
margins; panel members equal the Appendix; recorded hashes match the files present."""
import csv
import hashlib
import os

import pytest

from vulcan.gates import a2_config, a2_panels, a2_stats

CFG = a2_config.load()
ROOT = a2_config.ROOT


def test_seeds():
    s = CFG["seeds"]
    assert s["walk"] == 23 and s["file_pick_base"] == 23 and s["mixture_permutation"] == 23
    assert s["cluster_bootstrap"] == 83 and s["synthetic_checks"] == 131


def test_thresholds_match_draft():
    st = CFG["selection"]["stopping"]
    assert (st["min_projects"], st["min_A"], st["min_A_with_candidate"], st["min_tail_proxy"], st["cap_projects"]) == \
        (10, 120000, 72000, 400, 16)
    c = CFG["own_search"]["certification"]
    assert (c["q_psm"], c["q_peptide"]) == (0.001, 0.01)
    m = CFG["marker_screen"]
    assert (m["q_peptide"], m["min_distinct_peptides"], m["exclude_if_trigger_members_at_least"]) == (0.01, 2, 2)
    a = CFG["admission"]
    assert (a["identity_f_threshold"], a["identity_frag_ppm"], a["identity_prec_ppm"], a["hamming_max"]) == \
        (0.20, 20.0, 10.0, 2)
    assert abs(CFG["fs"]["delta_th"] - 2 / 7) < 1e-15
    g = CFG["gates"]
    t = g["tail_g2"]
    assert (t["n_min"], t["margin_abs"], t["margin_diff"], t["conf_one_sided"], t["bootstrap_B"], t["bootstrap_seed"],
            t["validity_wilson_conf_two_sided"]) == (200, 0.60, 0.10, 0.95, 10000, 83, 0.95)
    assert (g["g1"]["n_min_top"], g["g1"]["bound"], g["g1"]["conf_two_sided"]) == (7200, 0.05, 0.99)
    assert (g["g3"]["n_min_per_label_per_bin"], g["g3"]["tol"]) == (2000, 0.01)
    assert g["g2_original_g5"]["g5_alarm"] == 0.55
    sel = CFG["selection"]
    assert sel["runs_per_project"] == 3 and sel["raw_size_bytes"] == [3.0e8, 1.5e9]


def test_config_and_code_agree():
    t = CFG["gates"]
    assert a2_stats.BOOT_B == t["tail_g2"]["bootstrap_B"] and a2_stats.BOOT_SEED == t["tail_g2"]["bootstrap_seed"]
    assert a2_stats.G1_N_MIN == t["g1"]["n_min_top"] and a2_stats.G3_N_MIN == t["g3"]["n_min_per_label_per_bin"]
    from vulcan.gates import a2_admission as AD, a2_identity as I
    assert (AD.Q_PSM_CERT, AD.Q_PEP_CERT) == (0.001, 0.01)
    assert I.F_THRESHOLD == CFG["admission"]["identity_f_threshold"]
    assert list(I.ISOTOPE_ERRORS) == CFG["admission"]["identity_isotope_errors"]


def test_allow_list_matches_draft():
    genera = {a["genus"] for a in CFG["selection"]["allow"]}
    draft = {"Chlamydomonas", "Synechocystis", "Sordaria", "Phaeodactylum", "Streptomyces", "Myxococcus",
             "Cupriavidus", "Caballeronia", "Frankia", "Abortiporus", "Bacillus", "Komagataella", "Fervidobacterium",
             "Corynebacterium", "Emiliania", "Propionibacterium", "Deinococcus", "Chlorella", "Aeropyrum",
             "Acaryochloris", "Aspergillus", "Amycolatopsis", "Sorangium", "Neurospora", "Pseudomonas", "Rhodobacter",
             "Haloferax", "Halobacterium", "Shewanella", "Geobacter", "Trichoderma", "Schizosaccharomyces",
             "Caulobacter", "Sulfolobus", "Saccharolobus", "Thermococcus", "Methanosarcina", "Thermus", "Zymomonas",
             "Sinorhizobium", "Thalassiosira"}
    assert genera == draft


def test_panels_match_appendix():
    m = a2_panels.PANEL_MEMBERS
    assert len([x for x in m if x[0] == "H"]) == 38 + 16 + 17
    assert len([x for x in m if x[0] == "K"]) == 10
    assert len([x for x in m if x[0] == "B"]) == 13 + 4 + 1
    rows = list(csv.DictReader(open(os.path.join(ROOT, CFG["panels"]["table"])), delimiter="\t"))
    assert [(r["panel"], r["symbol"], int(r["taxid"])) for r in rows] == m
    assert all(r["accessions"] for r in rows)


def test_recorded_hashes():
    p = os.path.join(ROOT, CFG["panels"]["fasta"])
    if not os.path.exists(p):
        pytest.skip("panel FASTA not on this checkout (gitignored); re-run scripts/a2/resolve_panels.py")
    a2_config.check_panel_files(CFG)
    e = CFG["selection"]["exclusions"]
    accs = [r["accession"] for r in csv.DictReader(open(os.path.join(ROOT, e["at_adoption_file"])), delimiter="\t")]
    assert len(accs) == e["at_adoption_count"]
    assert hashlib.sha256("\n".join(sorted(accs)).encode()).hexdigest() == e["at_adoption_sha256"]


def test_seqclf_hash_if_present():
    if not os.path.exists(CFG["gates"]["seq_classifier"]["path"]):
        pytest.skip("classifier pickle not on this node")
    a2_config.check_seqclf(CFG)
