#!/usr/bin/env python3
"""Independent review check (M3 static scorer): training-set independence + label-related feature AUCs.
Dev/training data only.  Usage: conda_envs/pred/bin/python reviews/checks/train_independence.py"""
import csv, json, os, sys, glob
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAN = os.path.join(ROOT, "data", "manifests")
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

def col(path, c):
    with open(path) as fh:
        return [r[c] for r in csv.DictReader(fh, delimiter="\t")]

out = {}
train = list(csv.DictReader(open(os.path.join(MAN, "train.tsv")), delimiter="\t"))
tacc = {r["accession"] for r in train}
sets = {"t4_selection": set(col(os.path.join(MAN, "t4_selection.tsv"), "accession")),
        "tier_s": set(col(os.path.join(MAN, "tier_s.tsv"), "accession")),
        "dev_sealed": set(col(os.path.join(MAN, "dev_sealed.tsv"), "accession")) | set(col(os.path.join(MAN, "dev_sealed_extra.tsv"), "accession")),
        "predictor_training": set(col(os.path.join(MAN, "predictor_training_accessions.tsv"), "accession")),
        "known_T1_source_PXD028735": {"PXD028735"}}
out["overlap"] = {k: sorted(tacc & v) for k, v in sets.items()}
out["train_species"] = sorted({r["species"] for r in train})
# labs: training labs vs dev/sealed labs in frame
frame = {r["accession"]: r for r in csv.DictReader(open(os.path.join(MAN, "t4_frame.tsv")), delimiter="\t")}
devlabs = {frame[a]["lab"] for a in sets["tier_s"] | sets["dev_sealed"] if a in frame}
out["train_labs_in_dev_sealed_labs"] = sorted({r["lab"] for r in train} & devlabs)
out["t4_selection_labs_overlap_train_labs"] = sorted({r["lab"] for r in train} &
    set(col(os.path.join(MAN, "t4_selection.tsv"), "lab")))
# raw sha overlap
def shas(p):
    return {l.split()[0] for l in open(p) if l.strip()}
tsha = shas(os.path.join(MAN, "train.sha256"))
out["sha_overlap"] = {f: len(tsha & shas(os.path.join(MAN, f))) for f in ("tier_s.sha256", "dev_sealed.sha256")}
# T4 selector enforcement of training exclusion
src = open(os.path.join(ROOT, "scripts", "data", "select_t4.py")).read()
out["select_t4_reads_train_tsv"] = "train.tsv" in src
# folds
folds = {}
for r in train:
    folds.setdefault(r["accession"], set()).add(r["fold"])
out["fold_project_disjoint"] = all(len(v) == 1 for v in folds.values())
# feature AUCs, training rank-1 PSMs, target vs decoy
import m3_common as C, m3_data as MD
from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, feature_matrix
from sklearn.metrics import roc_auc_score
names = BASE_FEATURES + PRED_FEATURES
Xa, da, hsa = [], [], []
for sp in ("human", "mouse", "ecoli"):
    lab = MD.labels(C.index_dir(f"{sp}_contam", "prod"))
    for run, psm, feat in MD.load_runs(C.psm_dir(f"train_{sp}")):
        r1 = psm["rank"] == 1
        X = feature_matrix(psm, feat, names)[r1]
        Xa.append(X); da.append(lab["is_decoy"][psm["pf_row"]][r1]); hsa.append(psm["hyperscore"][r1])
X, d, hs = np.vstack(Xa), np.concatenate(da), np.concatenate(hsa)
y = ~d
# null-like region: bottom 50% hyperscore (mostly wrong PSMs for both labels)
low = hs <= np.quantile(hs, 0.5)
aucs = {}
for i, n in enumerate(names):
    try:
        aucs[n] = {"all": round(roc_auc_score(y, X[:, i]), 4), "low_half": round(roc_auc_score(y[low], X[low, i]), 4)}
    except ValueError:
        aucs[n] = None
out["n_rank1"] = int(len(y)); out["decoy_frac_rank1"] = float(d.mean()); out["decoy_frac_low_half"] = float(d[low].mean())
out["feature_auc_target_vs_decoy"] = aucs
# pure-sequence features specifically
out["sequence_only_features_low_half"] = {n: aucs[n]["low_half"] for n in ("peptide_len", "missed_cleavages", "n_var") if aucs.get(n)}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "train_independence.json"), "w"), indent=1)

# ---- per-feature target-vs-decoy AUC of rank-1 winners on G2 null-provenance spectra (D0 mixture; no homolog filter,
# so slightly conservative). Under a label-blind feature set every AUC should be ~0.5.
try:
    lab = MD.labels(C.index_dir("human_contam", "prod"))
    Xn, dn = [], []
    for run, psm, feat in MD.load_runs(C.psm_dir("G2_mix")):
        prov = np.load(os.path.join(C.WORK, "store_null", run.replace("_g2mix", "") + ".provenance.npy")) \
            if os.path.exists(os.path.join(C.WORK, "store_null", run.replace("_g2mix", "") + ".provenance.npy")) \
            else np.load(os.path.join(C.WORK, "store_null", run + ".provenance.npy"))
        r1 = (psm["rank"] == 1) & (prov[psm["spec_idx"]] > 0)
        Xn.append(feature_matrix(psm, feat, names)[r1]); dn.append(lab["is_decoy"][psm["pf_row"]][r1])
    Xn, dn = np.vstack(Xn), np.concatenate(dn)
    out["null_n_rank1"] = int(len(dn)); out["null_target_win_rate"] = float((~dn).mean())
    out["null_feature_auc"] = {n: round(roc_auc_score(~dn, Xn[:, i]), 4) for i, n in enumerate(names)}
except Exception as e:  # noqa: BLE001
    out["null_feature_auc_error"] = repr(e)
print(json.dumps({k: out[k] for k in out if k.startswith("null")}, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "train_independence.json"), "w"), indent=1)

# ---- spectrum-level dedup: hash(precursor_mz rounded 1e-4, top-5 peak m/z rounded 1e-3) of MS2 in training stores vs D0/D1
import pyarrow.ipc as ipc, hashlib
def spec_hashes(path):
    t = ipc.open_file(path).read_all()
    names_ = t.schema.names
    mzc = "mz" if "mz" in names_ else [n for n in names_ if n.endswith("mz") and "precursor" not in n and "sps" not in n][-1]
    ic = "intensity" if "intensity" in names_ else [n for n in names_ if "intens" in n and "precursor" not in n][-1]
    ms = t.column("ms_level").to_numpy() == 2
    pm = t.column("precursor_mz").to_numpy(zero_copy_only=False)
    mz = t.column(mzc).to_pylist(); it = t.column(ic).to_pylist()
    hs = set()
    for i in np.nonzero(ms)[0]:
        m, w = np.asarray(mz[i]), np.asarray(it[i])
        if m.size == 0: continue
        top = np.sort(np.round(m[np.argsort(w)[-5:]], 3))
        hs.add(hashlib.blake2b(f"{pm[i]:.4f}|{top.tolist()}".encode(), digest_size=8).hexdigest())
    return hs
try:
    dev = set()
    for p in glob.glob("/tmp/vulcan_e2/store/*.vspec.arrow"):
        dev |= spec_hashes(p)
    tr = set()
    for p in glob.glob(os.path.join(C.WORK, "store_train", "*.vspec.arrow")):
        tr |= spec_hashes(p)
    out["spectrum_hash_dedup"] = {"n_dev_D0D1": len(dev), "n_train": len(tr), "overlap": len(dev & tr)}
except Exception as e:  # noqa: BLE001
    out["spectrum_hash_dedup_error"] = repr(e)
print(json.dumps({k: out[k] for k in out if k.startswith("spectrum")}, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "train_independence.json"), "w"), indent=1)
