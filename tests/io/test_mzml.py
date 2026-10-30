"""mzML parser vs pyteomics on synthetic mzML (TRFP and msconvert attribute orders)."""
import base64
import zlib

import numpy as np
import pytest

from vulcan.io import mzml
from vulcan.io import store as St


def _arr(x, bits=64, comp=True):
    raw = np.asarray(x, dtype=np.float64 if bits == 64 else np.float32).tobytes()
    if comp:
        raw = zlib.compress(raw)
    return base64.b64encode(raw).decode()


def _cv(acc, name, value="", order="trfp"):
    if order == "trfp":
        return f'<cvParam cvRef="MS" accession="{acc}" value="{value}" name="{name}" />'
    return f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value="{value}"/>'


def _bda(arr, kind, bits, comp, order):
    acc, nm = ("MS:1000514", "m/z array") if kind == "mz" else ("MS:1000515", "intensity array")
    enc = _arr(arr, bits, comp)
    return (f'<binaryDataArray encodedLength="{len(enc)}">' + _cv(acc, nm, order=order)
            + _cv("MS:1000523" if bits == 64 else "MS:1000521", f"{bits}-bit float", order=order)
            + _cv("MS:1000574" if comp else "MS:1000576", "zlib compression" if comp else "no compression", order=order)
            + f"<binary>{enc}</binary></binaryDataArray>")


SPECS = [
    dict(scan=1, level=1, rt=0.5, mz=[400.1, 500.2, 600.3], it=[10, 20, 30]),
    dict(scan=2, level=2, rt=0.51, mz=[110.07, 245.1], it=[5, 6], prec=500.25, iso=500.2, z=2, ce=28, act="MS:1000422",
         filt="FTMS + c NSI d Full ms2 500.25@hcd28.00 [100.0000-1500.0000]", inj=22.5),
    dict(scan=3, level=2, rt=0.52, mz=[], it=[], prec=650.33, iso=650.3, z=0, ce=25, act="MS:1000598", supp="MS:1003181",
         filt="FTMS + c NSI d Full ms2 650.33@etd25.00@hcd20.00 [120.0000-1800.0000]", inj=50),
]


def _mzml(order="trfp", bits=64, comp=True):
    parts = ['<?xml version="1.0" encoding="utf-8"?><indexedmzML><mzML>',
             '<referenceableParamGroupList count="1"><referenceableParamGroup id="g">'
             + _cv("MS:1002877", "Q Exactive HF-X", order=order) + "</referenceableParamGroup></referenceableParamGroupList>",
             '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1"><componentList count="1">'
             '<analyzer order="2">' + _cv("MS:1000484", "orbitrap", order=order) + "</analyzer></componentList>"
             "</instrumentConfiguration></instrumentConfigurationList>",
             '<run id="synthetic" defaultInstrumentConfigurationRef="IC1"><spectrumList count="3">']
    for i, s in enumerate(SPECS):
        sid = f"controllerType=0 controllerNumber=1 scan={s['scan']}"
        x = [f'<spectrum id="{sid}" index="{i}" defaultArrayLength="{len(s["mz"])}">',
             _cv("MS:1000511", "ms level", s["level"], order), _cv("MS:1000130", "positive scan", order=order),
             '<scanList count="1"><scan instrumentConfigurationRef="IC1">',
             _cv("MS:1000016", "scan start time", s["rt"], order)]
        if "filt" in s:
            x.append(_cv("MS:1000512", "filter string", s["filt"], order))
            x.append(_cv("MS:1000927", "ion injection time", s["inj"], order))
        x.append("</scan></scanList>")
        if s["level"] == 2:
            x.append('<precursorList count="1"><precursor spectrumRef="controllerType=0 controllerNumber=1 scan=1"><isolationWindow>'
                     + _cv("MS:1000827", "isolation window target m/z", s["iso"], order)
                     + _cv("MS:1000828", "isolation window lower offset", 0.7, order)
                     + _cv("MS:1000829", "isolation window upper offset", 0.7, order)
                     + "</isolationWindow><selectedIonList count=\"1\"><selectedIon>"
                     + _cv("MS:1000744", "selected ion m/z", s["prec"], order)
                     + (_cv("MS:1000041", "charge state", s["z"], order) if s["z"] else "")
                     + "</selectedIon></selectedIonList><activation>"
                     + _cv("MS:1000045", "collision energy", s["ce"], order) + _cv(s["act"], "act", order=order)
                     + (_cv(s["supp"], "supp", order=order) if "supp" in s else "")
                     + "</activation></precursor></precursorList>")
        x.append('<binaryDataArrayList count="2">' + _bda(s["mz"], "mz", bits, comp, order)
                 + _bda(s["it"], "it", bits, comp, order) + "</binaryDataArrayList></spectrum>")
        parts.append("".join(x))
    parts.append("</spectrumList></run></mzML></indexedmzML>")
    return "\n".join(parts)


@pytest.mark.parametrize("order,bits,comp", [("trfp", 64, True), ("msconvert", 32, False), ("msconvert", 64, True)])
def test_mzml_parser_matches_pyteomics(tmp_path, order, bits, comp):
    from pyteomics import mzml as pmz
    p = tmp_path / "s.mzML"
    p.write_text(_mzml(order, bits, comp))
    out = str(tmp_path / "s.vspec.arrow")
    meta = mzml.convert(str(p), out, block=256)  # tiny blocks exercise the block cutting
    assert meta["n_spectra"] == 3 and meta["n_ms1"] == 1 and meta["n_ms2"] == 2
    t = St.read_table(out)
    P = St.peaks_numpy(t)
    ref = list(pmz.MzML(str(p)))
    for i, r in enumerate(ref):
        a, b = P.offsets[i], P.offsets[i + 1]
        assert np.array_equal(P.mz[a:b], np.asarray(r["m/z array"], np.float32))
        assert np.array_equal(P.intensity[a:b], np.asarray(r["intensity array"], np.float32))
        assert t.column("ms_level")[i].as_py() == r["ms level"]
        assert abs(t.column("rt_min")[i].as_py() - float(r["scanList"]["scan"][0]["scan start time"])) < 1e-9
    row = t.slice(1, 1).to_pylist()[0]
    assert row["precursor_mz"] == pytest.approx(500.25) and row["precursor_charge"] == 2
    assert row["activation"] == 2 and row["nce"] == pytest.approx(28) and row["master_scan"] == 1
    assert row["isolation_target_mz"] == pytest.approx(500.2) and row["isolation_lower_offset"] == pytest.approx(0.7)
    assert row["filter_string"] == SPECS[1]["filt"] and row["injection_time_ms"] == pytest.approx(22.5)
    assert row["analyzer"] == 1
    row = t.slice(2, 1).to_pylist()[0]
    assert row["precursor_charge"] is None and row["activation"] == 3 and row["supplemental_activation"] == 2
    assert t.slice(0, 1).to_pylist()[0]["precursor_mz"] is None
    assert meta["instrument_model"] == "Q Exactive HF-X"


def test_numpress_rejected(tmp_path):
    s = _mzml().replace("MS:1000574", "MS:1002312", 1)
    p = tmp_path / "n.mzML"
    p.write_text(s)
    with pytest.raises(NotImplementedError):
        mzml.convert(str(p), str(tmp_path / "n.vspec.arrow"))
