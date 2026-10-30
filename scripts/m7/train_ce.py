#!/usr/bin/env python3
"""Train the split-then-joint encoder (size fixed by E2b) on the M7 training set (scripts/m7/build_dataset.py).

* Project-disjoint folds (train_m7_folds.tsv): fold 4 = validation (early stopping on validation loss, LR choice);
  folds 0-3 = training.  No dev/eval data is touched.
* Loss = BCE(candidate logit, label) over every candidate of the batch (positives + incorrect-target negatives)
        + 0.5 * symmetric InfoNCE between pooled spectrum and positive-peptide embeddings (in-batch negatives).
* bf16 autocast, AdamW, linear warm-up + cosine decay; seed 11 (SEEDS[0]).
Output: results/2026-10-28_M7_encoders/models/ce_<size>.pt (+ .json with sha256, config, training log).

    python scripts/m7/train_ce.py --size s [--epochs 12 --lr 5e-4]
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402

from vulcan.encoders.model import SIZES, SplitJointEncoder, n_params  # noqa: E402

W = os.environ.get("VULCAN_M7_WORK", "/tmp/vulcan_m7")
RES = os.path.join(ROOT, "results", "2026-10-28_M7_encoders")
VAL_FOLD = int(os.environ.get("M7_VAL_FOLD", "4"))   # fixed: 4 (env override only for smoke tests)


class Data:
    def __init__(self, path, dev):
        z = np.load(path)
        self.n_spec = len(z["spec_fold"])
        self.fold = z["spec_fold"]
        self.mz = torch.from_numpy(z["spec_mz"]).to(dev)
        self.int = torch.from_numpy(z["spec_int"]).to(dev)
        self.n = torch.from_numpy(z["spec_n"].astype(np.int64)).to(dev)
        self.sch = torch.from_numpy(z["spec_charge"].astype(np.int64)).to(dev)
        sid = z["cand_sid"]
        o = np.argsort(sid, kind="stable")
        self.sid = sid[o]
        self.tok = torch.from_numpy(z["cand_tok"][o]).to(dev)
        self.cch = torch.from_numpy(z["cand_charge"][o].astype(np.int64)).to(dev)
        self.lab = torch.from_numpy(z["cand_label"][o].astype(np.float32)).to(dev)
        self.kind = z["cand_kind"][o]
        self.start = np.searchsorted(self.sid, np.arange(self.n_spec), "left")
        self.end = np.searchsorted(self.sid, np.arange(self.n_spec), "right")
        self.dev = dev

    def batch(self, spec_ids):
        spec_ids = np.asarray(spec_ids)
        cnt = self.end[spec_ids] - self.start[spec_ids]
        idx = np.concatenate([np.arange(a, b) for a, b in zip(self.start[spec_ids], self.end[spec_ids])])
        local = np.repeat(np.arange(len(spec_ids)), cnt)
        it = torch.from_numpy(idx).to(self.dev)
        tok = self.tok[it].long()
        L = int((tok > 0).sum(1).max())
        tok = tok[:, :L]
        sidt = torch.from_numpy(spec_ids).to(self.dev)
        spec = {"mz": self.mz[sidt], "int": self.int[sidt], "n": self.n[sidt], "charge": self.sch[sidt]}
        pep = {"tokens": tok, "charge": self.cch[it]}
        return spec, pep, torch.from_numpy(local).to(self.dev), self.lab[it], idx


def loss_fn(model, d: Data, ids, w_nce=0.5):
    spec, pep, local, y, idx = d.batch(ids)
    logit, se, pe = model(spec, pep, local)
    bce = F.binary_cross_entropy_with_logits(logit, y)
    posm = y > 0.5
    # one positive per spectrum (by construction); InfoNCE spectrum <-> positive peptide
    ps = local[posm]
    sc = model.logit_scale.exp().clamp(max=100)
    lg = sc * se[ps] @ pe[posm].T
    tgt = torch.arange(len(ps), device=lg.device)
    nce = 0.5 * (F.cross_entropy(lg, tgt) + F.cross_entropy(lg.T, tgt))
    return bce + w_nce * nce, bce, nce, logit, y


@torch.no_grad()
def evaluate(model, d: Data, ids, bs):
    model.eval()
    tot, tb, tn, nb = 0.0, 0.0, 0.0, 0
    top1, nsp = 0, 0
    allp, ally = [], []
    for b in range(0, len(ids), bs):
        bi = ids[b:b + bs]
        with torch.autocast("cuda", dtype=torch.bfloat16):
            l, bce, nce, logit, y = loss_fn(model, d, bi)
        tot += float(l); tb += float(bce); tn += float(nce); nb += 1
        # top-1 among each spectrum's candidates
        spec, pep, local, yy, _ = d.batch(bi)
        lg = logit.float()
        best = torch.full((len(bi),), -1e9, device=lg.device).scatter_reduce(0, local, lg, "amax")
        top1 += int(((lg >= best[local]) & (yy > 0.5)).sum())
        nsp += len(bi)
        allp.append(lg.cpu().numpy()); ally.append(yy.cpu().numpy())
    p, y = np.concatenate(allp), np.concatenate(ally)
    from sklearn.metrics import roc_auc_score
    model.train()
    return {"loss": tot / nb, "bce": tb / nb, "nce": tn / nb, "top1_acc": top1 / nsp, "auc": float(roc_auc_score(y, p))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", default="s")
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--bs", type=int, default=384, help="spectra per batch")
    ap.add_argument("--patience", type=int, default=2)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    torch.manual_seed(11)
    np.random.seed(11)
    rng = np.random.default_rng(11)
    dev = torch.device("cuda")
    d = Data(os.path.join(W, "dataset", "train.npz"), dev)
    tr = np.flatnonzero(d.fold != VAL_FOLD)
    va = np.flatnonzero(d.fold == VAL_FOLD)
    cfg = SIZES[a.size]
    model = SplitJointEncoder(cfg).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=0.01)
    steps_ep = math.ceil(len(tr) / a.bs)
    total = steps_ep * a.epochs
    warm = min(500, total // 10)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / warm) * 0.5 * (1 + math.cos(math.pi * min(s, total) / total)))
    log = {"size": a.size, "params": n_params(model), "config": cfg.to_dict(), "n_train_spectra": int(len(tr)),
           "n_val_spectra": int(len(va)), "lr": a.lr, "bs": a.bs, "epochs": [], "seed": 11}
    print(json.dumps({k: v for k, v in log.items() if k != "epochs"}), flush=True)
    best, best_state, bad = 1e9, None, 0
    t0 = time.time()
    for ep in range(a.epochs):
        perm = rng.permutation(tr)
        te = time.time()
        run_l = 0.0
        for i, b in enumerate(range(0, len(perm), a.bs)):
            with torch.autocast("cuda", dtype=torch.bfloat16):
                l, bce, nce, _, _ = loss_fn(model, d, perm[b:b + a.bs])
            opt.zero_grad(set_to_none=True)
            l.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
            run_l += float(l)
            if i % 500 == 0:
                print(f"ep {ep} step {i}/{steps_ep} loss {float(l):.4f} bce {float(bce):.4f} nce {float(nce):.4f}", flush=True)
        v = evaluate(model, d, va, 1024)
        rec = {"epoch": ep, "train_loss": run_l / steps_ep, "val": v, "wall_s": time.time() - te}
        log["epochs"].append(rec)
        print(json.dumps(rec), flush=True)
        if v["loss"] < best - 1e-4:
            best, bad = v["loss"], 0
            best_state = {k: t.detach().cpu().clone() for k, t in model.state_dict().items()}
            log["best_epoch"] = ep
        else:
            bad += 1
            if bad >= a.patience:
                break
    log["train_wall_s"] = time.time() - t0
    os.makedirs(os.path.join(RES, "models"), exist_ok=True)
    name = f"ce_{a.size}{a.tag}"
    buf = io.BytesIO()
    torch.save({"config": cfg.to_dict(), "state_dict": best_state}, buf)
    blob = buf.getvalue()
    open(os.path.join(RES, "models", name + ".pt"), "wb").write(blob)
    log["sha256"] = hashlib.sha256(blob).hexdigest()
    json.dump(log, open(os.path.join(RES, "models", name + ".json"), "w"), indent=1)
    print("saved", name, log["sha256"], "best epoch", log.get("best_epoch"), flush=True)


if __name__ == "__main__":
    main()
