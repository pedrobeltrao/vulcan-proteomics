"""P1 label-swap invariance of the A2.2 admission (scripts/a2/admit.admit_project) on synthetic data: swapping the T
and D lists (and permuting their order) leaves the admitted set bit-identical.  Also: the screens are not vacuous, the
gate DB can be read from an L1-style directory that has no labels.arrow, and outputs carry label-free keys only."""
import hashlib
import json
import os
import random
import re

import numpy as np
import pyarrow as pa

import admit as ADM
from vulcan.gates import a2_admission as AD
from vulcan.gates import a2_config
from vulcan.gates import a2_identity as I
from vulcan.gates.homolog import HomologFilter

AA = "ACDEFGHIKLMNPQRSTVWY"


def rev(s):
    return s[-2::-1] + s[-1]


def _pep(rng, lo=8, hi=18):
    return "".join(rng.choice(AA.replace("K", "").replace("R", "")) for _ in range(rng.randint(lo, hi))) + \
        rng.choice("KR")


def _world(seed=1):
    rng = random.Random(seed)
    T = sorted({_pep(rng) for _ in range(3000)})
    T = [(s, 0) for s in T]
    D = [(rev(s), 0) for s, _ in T]
    own_prot = ["".join(rng.choice(AA) for _ in range(400)) for _ in range(60)]
    own_peps = []
    for prot in own_prot:
        for m in re.finditer(r"[^KR]{7,20}[KR]", prot):
            own_peps.append(m.group(0))
    cert = {}
    for run in ("r1", "r2"):
        n_spec = 3000
        spec = sorted(rng.sample(range(n_spec), 500))
        ps = []
        for i in range(len(spec)):
            u = rng.random()
            if u < 0.08:                                    # shares a b/y ladder with a T peptide (same mass)
                t = rng.choice(T)[0]
                p = t[:-3] + t[-2] + t[-3] + t[-1]
            elif u < 0.14:                                  # one substitution from a D peptide
                d = list(rng.choice(D)[0])
                d[rng.randrange(len(d) - 1)] = rng.choice(AA)
                p = "".join(d)
            else:
                p = rng.choice(own_peps)
            ps.append(p)
        z = np.array([rng.choice([2, 3]) for _ in spec])
        mz = np.array([(I.peptide_mass(p) + zz * I.PROTON) / zz for p, zz in zip(ps, z)])
        cert[run] = {"spec": np.array(spec), "p_seq": np.array(ps, dtype=object), "p_mask": np.zeros(len(ps), np.uint64),
                     "z": z, "mz_mono": mz, "hyperscore_p": np.array([rng.uniform(15, 40) for _ in ps]),
                     "n_spectra": n_spec}
    return T, D, own_prot, cert


def _hash01(s):
    return int(hashlib.blake2b(s.encode(), digest_size=8).hexdigest(), 16) / 2 ** 64


def fake_exhaustive(db):
    """Label-free stand-in for the GPU exhaustive search: candidate score = 30 f(c) + 30 h(seq); rank 1 by score then
    sequence; partner = pseudo-reverse (a symmetric pairing)."""
    def run(_run, c):
        mx, top, par = [], [], []
        for p, z, mz in zip(c["p_seq"], c["z"], c["mz_mono"]):
            pi = I.ion_ladder(p, 0, first=2)
            best = None
            for r in db.window_rows(float(mz), int(z)):
                s = db.seq(r)
                sc = 30 * I.f_overlap(pi, I.ion_ladder(s, 0)) + 30 * _hash01(s)
                if best is None or (sc, s) > best:
                    best = (sc, s)
            mx.append(-np.inf if best is None else best[0])
            top.append(None if best is None else best[1])
            par.append(None if best is None else rev(best[1]))
        return {"max_hs": np.array(mx), "top_seq": top, "partner_seq": par}
    return run


def _admit(T, D, own_prot, cert):
    cfg = a2_config.load()
    db = I.GateDB.from_lists(T, D)
    ph = I.PeptideHamming(db.peptides())
    near = HomologFilter(own_prot).near
    return ADM.admit_project(cert, db, ph, fake_exhaustive(db), near, 25.0, cfg)


def _keys(o, acc):
    if isinstance(o, dict):
        for k, v in o.items():
            acc.add(str(k))
            _keys(v, acc)
    return acc


def test_swap_T_and_D_bit_identical():
    T, D, own, cert = _world()
    a = _admit(T, D, own, cert)
    b = _admit(D, T, own, cert)
    rng = random.Random(9)
    T2, D2 = T[:], D[:]
    rng.shuffle(T2)
    rng.shuffle(D2)
    c = _admit(D2, T2, own, cert)
    for x in (b, c):
        assert x["membership_sha256"] == a["membership_sha256"]
        assert x["counts"] == a["counts"]
        for run in a["codes"]:
            assert a["codes"][run].tobytes() == x["codes"][run].tobytes()
            for k in ("spec", "mz_mono", "z"):
                assert np.asarray(a["A"][run][k]).tobytes() == np.asarray(x["A"][run][k]).tobytes()


def test_screens_not_vacuous_and_label_free_output():
    T, D, own, cert = _world(2)
    a = _admit(T, D, own, cert)
    allc = np.concatenate(list(a["codes"].values()))
    for code in (AD.A, AD.REM_B_F, AD.REM_B_HAM, AD.REM_D):
        assert (allc == code).sum() > 0, AD.CODE_NAME[code]
    keys = _keys(a["counts"], set())
    assert keys <= set(AD.ADMISSION_COUNT_KEYS)
    assert not [k for k in keys if re.search(r"target|decoy|label", k, re.I)]


def test_gate_db_from_l1_needs_no_labels(tmp_path):
    T = [("ELVISLIVESK", 0), ("PEPTIDEKR", 0)]
    D = [(rev(s), m) for s, m in T]
    ref = I.GateDB.from_lists(T, D)
    seqs = [ref.seq(r) for r in range(len(ref))]
    res = np.frombuffer("".join(seqs).encode(), np.uint8)
    off = np.r_[0, np.cumsum([len(s) for s in seqs])[:-1]]
    t = pa.table({"mass32": pa.array(ref.mass.astype(np.float32)), "mask": pa.array(ref.mask),
                  "seq_off": pa.array(off.astype(np.int64)), "length": pa.array([len(s) for s in seqs], pa.uint8())})
    with pa.OSFile(str(tmp_path / "peptidoforms.arrow"), "wb") as sink, pa.ipc.new_file(sink, t.schema) as w:
        w.write_table(t)
    np.save(tmp_path / "residues.npy", res)
    json.dump({"params": {"variable_mods": {"M": 15.9949}}}, open(tmp_path / "meta.json", "w"))
    assert not os.path.exists(tmp_path / "labels.arrow")
    db = I.GateDB.from_l1(str(tmp_path))
    assert db.fingerprint() == ref.fingerprint()


def test_own_labels_refuses_gate_index(tmp_path):
    json.dump({"params": {}}, open(tmp_path / "meta.json", "w"))
    import pytest
    with pytest.raises(PermissionError):
        AD.own_labels(str(tmp_path))


def test_membership_hash_and_combine_order():
    codes = AD.combine(10, [1, 3, 5, 7], [True, False, False, False], [True, True, False, False],
                       [False, True, True, False], [False, False, True, True])
    assert codes.tolist() == [0, 2, 0, 3, 0, 4, 0, 5, 0, 0]          # first failing step wins
    assert AD.membership_hash({"b": codes, "a": codes}) == AD.membership_hash({"a": codes, "b": codes})
