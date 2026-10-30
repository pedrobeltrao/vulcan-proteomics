#!/usr/bin/env python3
"""3g step 1 (PLAN.md): homology measures of foreign pair base sequences E (and of the actual decoy sequences) against
the human+contaminant target proteome (I/L collapsed), for the E9 arm-C tail tables and a random sample of 20,000
arm-C foreign pairs (seed 11, background removal rates).

Per base sequence E:  ham2 (Hamming <= 2 to a same-length human target peptide), k6 / k5 (shares a 6-/5-mer with the
human proteome), sw (best Smith-Waterman local score, BLOSUM62, gap open 11 / extend 1, vs every human protein;
exhaustive numba DP) and its Karlin-Altschul E-value (lambda 0.267, K 0.041, n = proteome residues, m = len(E)).
Per sequence x in {E, dec(E)}:  h_y (longest C-terminal suffix of x found in human ending at a tryptic site or protein
C-term; capped at 12) and h_b (longest N-terminal prefix found in human starting after a tryptic site or at a protein
N-term / after Met1; capped at 12).
    conda_envs/pe/bin/python scripts/decoy_realism/homology.py
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numba as nb  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common as C  # noqa: E402

HUMAN = os.path.join(C.ROOT, "data", "cache", "M1e", "entrapment", "inputs", "human_contam.fasta")
KMAX = 12
AA = "ARNDCQEGHILKMFPSTWYV"
B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4"""
M = np.array([[int(x) for x in ln.split()] for ln in B62.strip().split("\n")], np.int16)
CODE = np.full(256, 20, np.int8)          # 20 = unknown / separator
for i, a in enumerate(AA):
    CODE[ord(a)] = i
SUB = np.full((21, 21), -4, np.int16)
SUB[:20, :20] = M
SUB[20, :] = -100
SUB[:, 20] = -100


def il(s):
    return s.replace("I", "L")


def read_fasta(p):
    seqs, cur = [], []
    for ln in open(p):
        if ln.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    return [il(s) for s in seqs]


@nb.njit(cache=True)
def sw_one(q, db, sub, go, ge):
    """Gotoh local alignment, query q vs concatenated db (separator code 20 scores -100). Returns the best score."""
    m = len(q)
    H = np.zeros(m + 1, np.int32)
    E = np.zeros(m + 1, np.int32)
    best = 0
    for j in range(len(db)):
        c = db[j]
        diag = 0
        F = 0
        Hprev_row_im1 = 0      # H[i-1] of the current column
        for i in range(1, m + 1):
            # E: gap in db direction (horizontal), from previous column
            e = max(E[i] - ge, H[i] - go - ge)
            E[i] = e
            F = max(F - ge, Hprev_row_im1 - go - ge)
            h = diag + sub[q[i - 1], c]
            if e > h:
                h = e
            if F > h:
                h = F
            if h < 0:
                h = 0
            diag = H[i]
            H[i] = h
            Hprev_row_im1 = h
            if h > best:
                best = h
    return best


@nb.njit(parallel=True, cache=True)
def sw_many(qs, qoff, qlen, db, sub, go, ge):
    out = np.zeros(len(qoff), np.int32)
    for k in nb.prange(len(qoff)):
        out[k] = sw_one(qs[qoff[k]:qoff[k] + qlen[k]], db, sub, go, ge)
    return out


@nb.njit(parallel=True, cache=True)
def ham_many(qs, qoff, qlen, refs_by_len, ref_off, ref_n, maxd):
    """min(Hamming, maxd+1) of each query vs all same-length refs (refs_by_len: concatenated blocks per length)."""
    out = np.full(len(qoff), maxd + 1, np.int32)
    for k in nb.prange(len(qoff)):
        L = qlen[k]
        if L >= len(ref_off):
            continue
        o, n = ref_off[L], ref_n[L]
        q = qs[qoff[k]:qoff[k] + L]
        best = maxd + 1
        for r in range(n):
            base = o + r * L
            d = 0
            for i in range(L):
                if refs_by_len[base + i] != q[i]:
                    d += 1
                    if d >= best:
                        break
            if d < best:
                best = d
                if best == 0:
                    break
        out[k] = best
    return out


def encode(seqs):
    lens = np.array([len(s) for s in seqs], np.int64)
    off = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
    buf = CODE[np.frombuffer("".join(seqs).encode(), np.uint8)]
    return buf, off, lens


def tryptic_sets(prots):
    suf, pre = set(), set()
    for s in prots:
        n = len(s)
        ends = [j for j in range(n) if (s[j] in "KR" and (j + 1 >= n or s[j + 1] != "P"))] + [n - 1]
        for j in ends:
            for k in range(2, KMAX + 1):
                if j - k + 1 < 0:
                    break
                suf.add(s[j - k + 1:j + 1])
        starts = [0, 1] + [j + 1 for j in ends if j + 1 < n]
        for st in starts:
            for k in range(2, KMAX + 1):
                if st + k > n:
                    break
                pre.add(s[st:st + k])
    return suf, pre


def h_y(x, suf):
    best = 0
    for k in range(2, min(KMAX, len(x)) + 1):
        if x[-k:] in suf:
            best = k
        else:
            break
    return best


def h_b(x, pre):
    best = 0
    for k in range(2, min(KMAX, len(x)) + 1):
        if x[:k] in pre:
            best = k
        else:
            break
    return best


def main():
    t0 = time.time()
    prots = read_fasta(HUMAN)
    n_res = sum(len(p) for p in prots)
    k6 = {p[i:i + 6] for p in prots for i in range(len(p) - 5)}
    k5 = {p[i:i + 5] for p in prots for i in range(len(p) - 4)}
    suf, pre = tryptic_sets(prots)
    print("sets", len(k6), len(k5), len(suf), len(pre), time.time() - t0, flush=True)
    labC = C.labels(C.IDX["C"])
    hum = np.unique(labC["target_seq"][labC["family"] == 0])
    hum = [il(s) for s in hum]
    # Hamming refs grouped by length
    maxL = 51
    blocks, ref_off, ref_n = [], np.zeros(maxL, np.int64), np.zeros(maxL, np.int64)
    pos = 0
    byL = {}
    for s in hum:
        byL.setdefault(len(s), []).append(s)
    for L in range(maxL):
        ss = byL.get(L, [])
        ref_off[L], ref_n[L] = pos, len(ss)
        if ss:
            blocks.append(CODE[np.frombuffer("".join(ss).encode(), np.uint8)])
            pos += L * len(ss)
    refs = np.concatenate(blocks)
    db = CODE[np.frombuffer(("*".join(prots)).encode(), np.uint8)]
    # query sets: tail tables (base seqs) + random 20k arm-C foreign pairs
    rng = np.random.default_rng(C.SEED)
    eidx = np.flatnonzero(labC["family"] == 1)
    samp = rng.choice(eidx, 20000, replace=False)
    rnd_base = np.array([il(s) for s in labC["target_seq"][samp]])
    tabs = {}
    for k in ("C_score", "C_hyperscore"):
        p = os.path.join(C.SCR, f"e9tail_{k}.parquet")
        if os.path.exists(p):
            tabs[k] = pd.read_parquet(p)
    allq = sorted(set(rnd_base.tolist()).union(*[set(t.base_seq.map(il)) for t in tabs.values()]))
    print("queries", len(allq), flush=True)
    qs, qoff, qlen = encode(allq)
    t1 = time.time()
    ham = ham_many(qs, qoff, qlen, refs, ref_off, ref_n, 2)
    print("ham", time.time() - t1, flush=True)
    t1 = time.time()
    sw = sw_many(qs, qoff, qlen, db, SUB, 11, 1)
    t_sw = time.time() - t1
    print("sw", t_sw, flush=True)
    lam, K = 0.267, 0.041
    ev = K * qlen * n_res * np.exp(-lam * sw.astype(float))
    selfsc = np.array([sum(int(M[CODE[ord(c)], CODE[ord(c)]]) for c in s) for s in allq])
    H = pd.DataFrame({"seq": allq, "ham": ham, "k6": [any(s[i:i + 6] in k6 for i in range(len(s) - 5)) for s in allq],
                      "k5": [any(s[i:i + 5] in k5 for i in range(len(s) - 4)) for s in allq],
                      "sw": sw, "sw_evalue": ev, "sw_norm": sw / selfsc,
                      "hy": [h_y(s, suf) for s in allq], "hb": [h_b(s, pre) for s in allq]})
    H.to_parquet(os.path.join(C.SCR, "homology_base.parquet"))
    # decoy-side ladders
    decs = sorted(set().union(*[set(t.dec_seq.dropna().map(il)) for t in tabs.values()]))
    D = pd.DataFrame({"seq": decs, "hy": [h_y(s, suf) for s in decs], "hb": [h_b(s, pre) for s in decs]})
    D.to_parquet(os.path.join(C.SCR, "homology_dec.parquet"))
    pd.DataFrame({"seq": rnd_base}).to_parquet(os.path.join(C.SCR, "random_base.parquet"))
    json.dump({"n_res": n_res, "n_queries": len(allq), "sw_s": t_sw, "total_s": time.time() - t0},
              open(os.path.join(C.SCR, "homology_meta.json"), "w"))


if __name__ == "__main__":
    main()
