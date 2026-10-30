#!/usr/bin/env python3
"""Diagnostic (not a gate): what are the offset-null spectra that reach the 1% tail?  For mixtures 00/01, each offset
copy is mapped back to its source D0 spectrum; its winning candidate (frozen scorer and hyperscore) is compared with
the source spectrum's own production rank-1 peptide (hyperscore search of the unshifted spectrum, same mixture)."""
import json, os, pickle, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3")); sys.path.insert(0, HERE)
import m3_common as C, m3_data as MD  # noqa: E402
import numpy as np  # noqa: E402
from gates_a1 import _load, neutral_rt  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402
from vulcan.index.build import L1Index  # noqa: E402
from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition, psm_competition  # noqa: E402
W = "/tmp/vulcan_m3e"
idx = C.index_dir("human_contam", "prod"); lab = MD.labels(idx); l1 = L1Index(idx)
model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
out = []
for name in ("g2mix00_h0", "g2mix01_h1"):
    run, psm, feat = _load(os.path.join(W, "psms", "G2mix"), [name])[0]
    prov = np.load(os.path.join(W, "null_mix", name + ".prov.npy"))
    sp = pickle.load(open(os.path.join(W, "null_mix", name + ".spectra.pkl"), "rb"))
    nhost = int((prov == 0).sum()); known = np.flatnonzero(sp.charge[:nhost] > 0); nk = known.size
    rows, spec = psm["pf_row"], psm["spec_idx"]; dec = lab["is_decoy"][rows]; tb = hash_tiebreak(psm["pf_id"], 11)
    for pname, sc in (("hyperscore", psm["hyperscore"].astype(float)), ("gbm", model.apply(psm, feat))):
        w = psm_competition(spec, sc, tiebreak=tb)
        r = double_competition(spec.astype(np.int64), psm["pf_id"], lab["tdc_pair"][rows], sc, dec, seed=11, tiebreak=tb)
        thr = GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01)
        best = {int(spec[i]): i for i in w}
        for i in w:
            s = int(spec[i])
            if prov[s] != 1 or sc[i] < thr:
                continue
            src = int(known[(s - nhost) % nk])
            j = best.get(src)
            out.append({"mix": name, "pipe": pname, "cand": l1.sequence(int(rows[i])), "decoy": bool(dec[i]),
                        "cand_charge": int(psm["charge"][i]), "iso": int(psm["isotope_error"][i]),
                        "delta_mass_ppm": float(psm["delta_mass_ppm"][i]),
                        "src_rank1": None if j is None else l1.sequence(int(rows[j])),
                        "src_rank1_decoy": None if j is None else bool(dec[j]),
                        "precursor_shift_da": float((sp.precursor_mz[s] - sp.precursor_mz[src]) * sp.charge[s])})
json.dump(out, open(os.path.join(W, "diag_offset_tail.json"), "w"), indent=1)
for o in out[:40]:
    print(o["pipe"], o["cand"], o["decoy"], o["cand_charge"], o["iso"], round(o["delta_mass_ppm"], 1), "| src:", o["src_rank1"], round(o["precursor_shift_da"], 2))
