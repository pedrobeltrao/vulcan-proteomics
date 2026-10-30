"""vulcan.quant unit tests (synthetic): impurity correction, MaxLFQ, alignment, transfer q-values, SPS-MS3 linkage,
metrics."""
import numpy as np
import pandas as pd

from vulcan.quant import lfq as Q
from vulcan.quant import metrics as MET
from vulcan.quant import tmt as TM


def test_impurity_matrix_and_correction():
    A = TM.impurity_matrix()
    assert A.shape == (11, 11)
    assert (A.sum(axis=0) <= 1 + 1e-12).all() and (np.diag(A) > 0.8).all()   # mass leaving the 11 channels is lost
    rng = np.random.default_rng(0)
    X = rng.uniform(1e3, 1e5, size=(50, 11))
    Y = X @ A.T
    Xc = TM.correct(Y, A)
    assert np.allclose(Xc, X, rtol=1e-6)


def test_link_ms3():
    rows = pd.DataFrame({
        "index": [0, 1, 2, 3, 4], "scan": [1, 2, 3, 4, 5], "ms_level": [1, 2, 3, 2, 3],
        "master_scan": [None, 1, 2, 1, 4],
        "filter_string": ["FTMS + c NSI Full ms", "ITMS + c NSI r d Full ms2 783.7729@cid35.00 [211-1578]",
                          "FTMS + c NSI d Full ms3 783.7729@cid35.00 468.2771@hcd65.00 [100-500]",
                          "ITMS + c NSI r d Full ms2 570.7941@cid35.00 [153-1152]",
                          "FTMS + c NSI d Full ms3 570.7000@cid35.00 442.9741@hcd65.00 [100-500]"]})
    m = TM.link_ms3(rows).set_index("index")
    assert m.loc[1, "rep_index"] == 2 and bool(m.loc[1, "link_ok"])
    assert not bool(m.loc[3, "link_ok"])                      # parent m/z mismatch -> not linked


def test_extract_reporters():
    mz = np.array([126.127726, 127.124761 * (1 + 2e-6), 127.131081, 131.1445 * (1 - 30e-6)], np.float64)
    it = np.array([10.0, 20.0, 30.0, 40.0])
    I, E = TM.extract(mz, it, np.array([0, 4]), np.array([0]), 10.0)
    assert I[0, 0] == 10 and I[0, 1] == 20 and I[0, 2] == 30 and I[0, 10] == 0   # 131C is 30 ppm off
    assert abs(E[0, 1] - 2.0) < 1e-3


def test_maxlfq_recovers_ratios():
    rng = np.random.default_rng(1)
    runs = [f"r{i}" for i in range(6)]
    true = np.array([0, 0.5, 1.0, -0.5, 0.2, 0.0])
    rows, prot = [], []
    for p in range(20):
        base = rng.uniform(10, 20, size=8)
        for k in range(8):
            v = 2 ** (base[k] + true + rng.normal(0, 0.01, 6))
            v[rng.random(6) < 0.2] = np.nan
            rows.append(v)
            prot.append(f"P{p}")
    M = pd.DataFrame(rows, columns=runs, index=[f"x{i}" for i in range(len(rows))])
    P = Q.maxlfq(M, pd.Series(prot, index=M.index))
    L = np.log2(P)
    d = (L.sub(L["r0"], axis=0)).median()
    assert np.allclose(d.to_numpy(), true - true[0], atol=0.03)


def test_alignment_recovers_shift():
    rng = np.random.default_rng(2)
    n = 2000
    rt0 = np.sort(rng.uniform(1, 60, n))
    tabs = {}
    for run, f in (("a", lambda x: x), ("b", lambda x: x * 1.05 + 0.5), ("c", lambda x: x - 0.3 + 0.002 * x ** 1.5)):
        tabs[run] = pd.DataFrame({"pf_id": np.arange(n, dtype=np.uint64), "charge": 2,
                                  "rt_apex": f(rt0) + rng.normal(0, 0.02, n)})
    al = Q.Alignment(tabs)
    for run in ("a", "b", "c"):
        back = al.to_run(run, al.to_consensus(run, tabs[run]["rt_apex"].to_numpy()))
        assert np.median(np.abs(back - tabs[run]["rt_apex"].to_numpy())) < 1e-6
    x = al.to_consensus("b", tabs["b"]["rt_apex"].to_numpy()) - al.to_consensus("a", tabs["a"]["rt_apex"].to_numpy())
    assert np.median(np.abs(x)) < 0.05


def test_transfer_q_null_calibration():
    """Pure null candidates: target and decoy windows exchangeable -> few accepted at 1%."""
    rng = np.random.default_rng(3)
    n = 20000
    T = pd.DataFrame({"St": rng.normal(size=n), "Sd_rt": rng.normal(size=n), "Sd_mz": rng.normal(size=n),
                      "is_decoy_donor": np.zeros(n, bool)})
    q, qp = Q.transfer_qvalues(T)
    assert (q <= 0.01).sum() <= 10 and np.all(np.isnan(q) | (q >= qp - 1e-12))
    # half true transfers with high scores -> accepted, and decoy-donor rows carry no q
    T2 = T.copy()
    T2.loc[: n // 2, "St"] += 6
    T2.loc[n - 10:, "is_decoy_donor"] = True
    q2, _ = Q.transfer_qvalues(T2)
    assert (q2 <= 0.01).sum() > n // 4 and np.isnan(q2[-10:]).all()


def test_epsilon_and_cv():
    runs = ["x_Condition_A_1", "x_Condition_A_2", "x_Condition_A_3", "x_Condition_B_1", "x_Condition_B_2", "x_Condition_B_3"]
    A, B = runs[:3], runs[3:]
    M = pd.DataFrame([[100, 100, 100, 100, 100, 100], [200, 200, 200, 100, 100, 100], [25, 25, 25, 100, 100, 100]],
                     columns=runs, index=["h", "y", "e"], dtype=float)
    sp = pd.Series({"h": "HUMAN", "y": "YEAST", "e": "ECOLI"})
    e = MET.epsilon(M, sp, A, B, {"HUMAN": 0.0, "YEAST": 1.0, "ECOLI": -2.0})
    assert e["n"] == 3 and e["median_abs_eps"] < 1e-12
    assert MET.cv_stats(M, A, B)["median_cv"] == 0.0


def test_combine_isotopes_avgproj():
    """M6d: the averagine projection equals the plain sum for an ideal envelope and down-weights a censored M+2."""
    from vulcan.ms1.averagine import envelope
    from vulcan.quant.lfq import combine_isotopes
    m = np.array([1500.0, 2500.0])
    env = np.stack([envelope(x, 3) for x in m]) * 1e6
    s = combine_isotopes(env, m, "avgproj")
    assert np.allclose(s, env.sum(1), rtol=2e-3)
    cens = env.copy()
    cens[:, 2] = 0.0                       # M+2 lost below the centroiding threshold
    loss_sum = 1 - combine_isotopes(cens, m, "sum") / env.sum(1)
    loss_proj = 1 - combine_isotopes(cens, m, "avgproj") / env.sum(1)
    assert np.all(loss_proj < loss_sum)


def test_xic_iso_window():
    """M6d kernel: a Gaussian elution of a 3-isotope envelope integrates to the analytic area per isotope."""
    from vulcan.quant.lfq import NEUTRON, _xic_iso
    rt = np.linspace(0, 2, 401)
    z, mz0 = 2.0, 600.0
    rel = np.array([1.0, 0.6, 0.25])
    prof = np.exp(-0.5 * ((rt - 1.0) / 0.05) ** 2)
    mz, it, off = [], [], [0]
    for p in prof:
        for j in range(3):
            mz.append(mz0 + j * NEUTRON / z)
            it.append(1e6 * rel[j] * p)
        off.append(len(mz))
    a, nsc = _xic_iso(rt, np.array(off, np.int64), np.array(mz), np.array(it), np.array([mz0]), np.array([z]),
                      np.array([0.5]), np.array([1.5]), 5.0, 3)
    exp = 1e6 * rel * 0.05 * np.sqrt(2 * np.pi)
    assert np.allclose(a[0], exp, rtol=1e-3)
    assert nsc[0] == 201
