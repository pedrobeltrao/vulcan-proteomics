"""Dry run of the DN2 selector (scripts/a2/select_dn2.py) on a SYNTHETIC project table with fake I/O.  No network, no
PRIDE, no files: the real t4_frame is not read here."""
import random

import pytest

import select_dn2 as SEL
from vulcan.gates import a2_config

CFG = a2_config.load()


def row(acc, org, **kw):
    r = {"accession": acc, "organisms": org, "acquisition": "DDA", "vendor": "thermo", "label": "LFQ",
         "enrichment": "none", "instruments": "Q Exactive HF", "title": "whole-cell proteome study", "keywords": "",
         "experiment_types": "Shotgun proteomics"}
    r.update(kw)
    return r


GENERA = ["Chlamydomonas reinhardtii", "Synechocystis sp. PCC 6803", "Streptomyces coelicolor", "Myxococcus xanthus",
          "Cupriavidus necator", "Komagataella phaffii", "Deinococcus radiodurans", "Neurospora crassa",
          "Pseudomonas putida KT2440", "Haloferax volcanii", "Shewanella oneidensis", "Geobacter sulfurreducens",
          "Caulobacter crescentus", "Thermococcus kodakarensis", "Thermus thermophilus", "Zymomonas mobilis",
          "Sinorhizobium meliloti", "Thalassiosira pseudonana", "Schizosaccharomyces pombe", "Trichoderma reesei"]


def synthetic_frame():
    rows, k = [], 100000
    for g in GENERA:
        for _ in range(2):
            k += 7
            rows.append(row(f"PXD{k}", g))
    k += 1
    rows += [row(f"PXD{k + 1}", "Bacillus cereus ATCC 14579"),           # cereus group: denied
             row(f"PXD{k + 2}", "Bacillus subtilis 168"),                # allowed
             row(f"PXD{k + 3}", "Corynebacterium diphtheriae"),          # species restriction
             row(f"PXD{k + 4}", "Corynebacterium glutamicum"),
             row(f"PXD{k + 5}", "Sulfolobus acidocaldarius"),            # shares a slot with Saccharolobus
             row(f"PXD{k + 6}", "Saccharolobus solfataricus"),
             row(f"PXD{k + 7}", "Staphylococcus aureus"),                # previous null species, not on ALLOW
             row(f"PXD{k + 8}", "Homo sapiens"),
             row(f"PXD{k + 9}", "Streptomyces coelicolor;Homo sapiens"),  # not single organism
             row(f"PXD{k + 10}", "Neurospora crassa", title="Neurospora skin of the fungus"),   # BAD_KW 'skin'
             row(f"PXD{k + 11}", "Neurospora crassa", instruments="LTQ Orbitrap"),
             row(f"PXD{k + 12}", "Neurospora crassa", label="TMT"),
             row(f"PXD{k + 13}", "Thermus thermophilus", keywords="HeLa spike-in"),
             row("PXD023811", "Chlamydomonas reinhardtii"),               # excluded accession (DN)
             row("MSV000099999", "Chlamydomonas reinhardtii")]           # non-PXD
    return rows


class FakeIO:
    def __init__(self, seed=0, bad_protocol=(), marker_exclude=(), few_raw=(), no_proteome=()):
        self.r = random.Random(seed)
        self.bad_protocol, self.marker_exclude = set(bad_protocol), set(marker_exclude)
        self.few_raw, self.no_proteome = set(few_raw), set(no_proteome)
        self.calls = []

    def project_record(self, acc):
        self.calls.append(("record", acc))
        return {"description": "grown in brain heart infusion" if acc in self.bad_protocol else "minimal medium"}

    def sdrf(self, acc):
        return []

    def list_files(self, acc):
        self.calls.append(("list", acc))
        n = 2 if acc in self.few_raw else 6
        fs = [{"fileName": f"{acc}_run{i}.raw", "fileSizeBytes": 5e8, "url": f"x://{acc}/{i}"} for i in range(n)]
        fs += [{"fileName": f"{acc}_QC_{i}.raw", "fileSizeBytes": 5e8, "url": "x"} for i in range(3)]
        fs += [{"fileName": f"{acc}_big.raw", "fileSizeBytes": 2e9, "url": "x"}]
        return fs

    def proteome(self, r):
        return None if r["accession"] in self.no_proteome else {"upid": "UP0001", "taxid": "1"}

    def download(self, acc, files):
        self.calls.append(("download", acc))
        return [{"filename": f["fileName"], "url": f["url"], "expected_size": f["fileSizeBytes"], "sha256": "0" * 64}
                for f in files]

    def marker_and_admit(self, acc, files, prot):
        self.calls.append(("admit", acc))
        excl = acc in self.marker_exclude
        mk = {"exclude": excl, "trigger_members": ["H|DSP|9606", "H|DCD|9606"] if excl else [],
              "n_trigger_members_detected": 2 if excl else 0, "n_panel_members_detected": 2 if excl else 0,
              "detected_members": []}
        c = None if excl else {"n_A": 13000, "n_A_with_candidate": 8000, "tail_proxy": 45, "n_spectra": 90000}
        return {"marker": mk, "counts": c}


def wio(f):
    return SEL.WalkIO(f.project_record, f.sdrf, f.list_files, f.proteome, f.download, f.marker_and_admit)


def test_frame_filter_and_exclusions():
    excl = SEL.load_exclusions(CFG)
    assert len(excl) >= 300 and "PXD023811" in excl
    cands, drops = SEL.frame_candidates(synthetic_frame(), CFG, excl)
    accs = {r["accession"] for r in cands}
    orgs = {r["organisms"] for r in cands}
    assert "PXD023811" not in accs and "MSV000099999" not in accs
    assert "Bacillus cereus ATCC 14579" not in orgs and "Bacillus subtilis 168" in orgs
    assert "Corynebacterium diphtheriae" not in orgs and "Corynebacterium glutamicum" in orgs
    assert not any(o.startswith(("Staphylococcus", "Homo")) or ";" in o for o in orgs)
    assert drops["BAD_KW"] == 2 and drops["instrument"] == 1 and drops["frame fields"] == 1
    assert [r["accession"] for r in cands] == sorted(accs)


def test_shuffle_is_random23_on_sorted():
    cands, _ = SEL.frame_candidates(synthetic_frame(), CFG, SEL.load_exclusions(CFG))
    s = SEL.shuffled(list(reversed(cands)), CFG)
    ref = sorted(cands, key=lambda r: r["accession"])
    random.Random(23).shuffle(ref)
    assert [r["accession"] for r in s] == [r["accession"] for r in ref]


def test_walk_stopping_one_per_genus_replacement():
    cands, _ = SEL.frame_candidates(synthetic_frame(), CFG, SEL.load_exclusions(CFG))
    order = SEL.shuffled(cands, CFG)
    accs = [r["accession"] for r in order]
    io = FakeIO(bad_protocol=accs[:1], marker_exclude=accs[1:3], few_raw=accs[3:4], no_proteome=accs[4:5])
    w = SEL.walk(order, CFG, wio(io))
    adm = w["admitted"]
    # 13,000 A per project -> |A| >= 120,000 after 10 projects; proxy 45 x 10 >= 400; candidates 80,000 >= 72,000
    assert len(adm) == 10 and w["final_state"]["reason"] == "criteria met"
    slots = [SEL.genus_slot(next(r["organisms"] for r in order if r["accession"] == a), CFG) for a in adm]
    assert len(set(slots)) == len(slots)
    res = {x["accession"]: x["result"] for x in io_log(w)}
    assert res[accs[0]] == "protocol" and res[accs[1]] == "marker screen exclude"
    assert res[accs[3]].endswith("usable raw") and res[accs[4]] == "no reference proteome"
    # metadata-rejected projects are never listed / downloaded
    assert ("list", accs[0]) not in io.calls and ("download", accs[3]) not in io.calls
    # files: 3 per admitted project, never QC / oversized
    assert len(w["rows"]) == 30 and not any("QC" in r["filename"] or "big" in r["filename"] for r in w["rows"])
    # deterministic
    w2 = SEL.walk(order, CFG, wio(FakeIO(bad_protocol=accs[:1], marker_exclude=accs[1:3], few_raw=accs[3:4],
                                         no_proteome=accs[4:5])))
    assert w2["rows"] == w["rows"] and w2["admitted"] == adm


def io_log(w):
    return w["log"]


def test_walk_cap_16():
    cands, _ = SEL.frame_candidates(synthetic_frame(), CFG, SEL.load_exclusions(CFG))

    class Small(FakeIO):
        def marker_and_admit(self, acc, files, prot):
            r = super().marker_and_admit(acc, files, prot)
            r["counts"] = {"n_A": 1000, "n_A_with_candidate": 500, "tail_proxy": 3, "n_spectra": 9000}
            return r
    w = SEL.walk(SEL.shuffled(cands, CFG), CFG, wio(Small()))
    assert len(w["admitted"]) == 16 and w["final_state"]["reason"] == "cap reached"


def test_pick_files_seed():
    raws = [{"fileName": f"f{i}.raw"} for i in range(8)]
    p = SEL.pick_files(raws, "PXD012345", CFG)
    ref = sorted(random.Random(23 + 12345).sample(raws, 3), key=lambda f: f["fileName"])
    assert p == ref


def test_protocol_and_sdrf_rules():
    assert SEL.protocol_reject({"sampleProcessingProtocol": "cells grown on blood agar"}, [], CFG)
    assert SEL.protocol_reject({"description": "searched against the human proteome"}, [], CFG)
    assert SEL.protocol_reject({}, [], CFG) is None                                   # missing text never excludes
    sd = [{"characteristics[organism]": "Neurospora crassa", "characteristics[organism part]": "whole organism"}]
    assert SEL.protocol_reject({}, sd, CFG) is None
    sd2 = sd + [{"characteristics[organism]": "Homo sapiens", "characteristics[organism part]": "whole organism"}]
    assert "SDRF organisms" in SEL.protocol_reject({}, sd2, CFG)
    sd3 = [{"characteristics[organism]": "Neurospora crassa", "characteristics[organism part]": "mycelium tip"}]
    assert "organism part" in SEL.protocol_reject({}, sd3, CFG)


def test_selector_refuses_without_tag():
    import subprocess
    have = subprocess.run(["git", "-C", SEL.ROOT, "rev-parse", "-q", "--verify", "refs/tags/prereg-v0.4.2"],
                          capture_output=True).returncode == 0
    if have:
        pytest.skip("tag exists")
    with pytest.raises(SystemExit):
        SEL.require_tag("prereg-v0.4.2")


def test_a1_bad_kw_substring_quirks_kept_verbatim():
    # A1's BAD_KW is reused verbatim (A2.3 step 1).  Without word boundaries it also matches inside words; the
    # stricter reading keeps that (documented in IMPLEMENTATION_NOTES.md).
    rx = SEL.regexes(CFG)["bad_kw"]
    assert rx.search("Chlamydomonas proteome")            # 'HLA'
    assert rx.search("growth media")                      # 'DIA\b'
    assert rx.search("protein groups")                    # 'UPS'
