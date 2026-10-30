"""vulcan.predict.apd vectorized input path == peptdeep's own featurization (same model, same outputs); D4 guard.
Runs only in conda_envs/pred (peptdeep installed) with CUDA."""
import numpy as np
import pytest

pytest.importorskip("peptdeep")
torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("needs CUDA", allow_module_level=True)


def test_d4_guard():
    from vulcan.predict.apd import check_allowed
    for bad in ("D2", "T2"):
        with pytest.raises(PermissionError):
            check_allowed(bad)
    with pytest.raises(PermissionError):
        check_allowed(accession="PXD009449")
    check_allowed("D0")


def test_fast_path_equals_peptdeep_featurization():
    import pandas as pd
    from peptdeep.model.featurize import get_batch_aa_indices, get_batch_mod_feature, parse_instrument_indices
    from vulcan.predict.apd import APD
    seqs = ["PEPTCIDEK", "MLMAARK", "ACDEFGHLKMNPQR", "LLLLLLLK"]
    masks = [0, 0b101, 1 << 9, 0]
    m = APD()
    res = np.frombuffer("".join(seqs).encode(), np.uint8).copy()
    ln = np.array([len(s) for s in seqs])
    off = np.r_[0, np.cumsum(ln)[:-1]]
    z = np.array([2, 3, 2, 4])
    o = m.predict(res, off, ln, np.array(masks, np.int64), z, 28.0)
    mdl = m.mm.ms2_model
    for i, (s, mk) in enumerate(zip(seqs, masks)):
        mods, sites = [], []
        for k, c in enumerate(s):
            if c == "C":
                mods.append("Carbamidomethyl@C"); sites.append(str(k + 1))
            if c == "M" and (mk >> k) & 1:
                mods.append("Oxidation@M"); sites.append(str(k + 1))
        df = pd.DataFrame({"sequence": [s], "mods": [";".join(mods)], "mod_sites": [";".join(sites)], "nAA": [len(s)]})
        with torch.inference_mode():
            aa = torch.tensor(get_batch_aa_indices(df.sequence.values.astype("U")), dtype=torch.long, device="cuda")
            mx = torch.tensor(get_batch_mod_feature(df), dtype=torch.float32, device="cuda")
            rt = m.rt(aa, mx).item()
            p = mdl.model(aa, mx, torch.tensor([[z[i] * mdl.charge_factor]], device="cuda", dtype=torch.float32),
                          torch.tensor([[28.0 * mdl.NCE_factor]], device="cuda"),
                          torch.tensor(parse_instrument_indices(["QE"]), device="cuda")).cpu().numpy()[0][:, m.col]
        p = p / max(p.max(), 1e-30)
        p[p < mdl.min_inten] = 0
        assert abs(rt - o["rt"][i]) < 1e-6
        got = o["ms2"][o["ms2_off"][i]:o["ms2_off"][i + 1]].reshape(-1, 4).astype(np.float32)
        assert np.abs(got - p).max() < 1e-3
