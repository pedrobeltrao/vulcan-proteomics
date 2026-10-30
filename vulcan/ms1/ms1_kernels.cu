// E11: GPU MS1 feature finder (hill detection + isotope grouping), opt-in backend `ms1_backend=gpu`.
//
// Reproduces vulcan/ms1/hills.py (detect_hills, gap=0 fused path) and vulcan/ms1/isotopes.py (group_isotopes fast
// path, up to the greedy acceptance) BITWISE.  The CPU numba path stays the reference; tests/ms1/test_gpu_exact.py
// and scripts/e11/equal_gpu.py check every array.
//
// How exactness is kept:
//   * every floating-point reduction is a sequential loop inside ONE thread, in the CPU loop order (per hill: split
//     smoothing, summaries, sums of squares; per lookup: trace-cosine dot product; per chain: averagine cosine);
//   * numba's typing is replicated: f32*f32 and f32+f32 round to f32 (__fmul_rn/__fadd_rn) before widening to f64,
//     f32 op f64 widens first; nvcc -fmad=false (no FMA contraction, numba/LLVM does not contract by default);
//     double sqrt/div are IEEE (correctly rounded) on both sides;
//   * every order-dependent integer step is order-exact: the one-to-one link claim is (min distance, then min index)
//     via two atomicMin passes = the CPU's ascending strict-< scan; peaks are placed by list ranking (chain depth =
//     index order within a hill = the stable counting sort); the m/z argsort and the scan-block bucketing use CUB's
//     stable radix sort; candidate order is the row-major nonzero order of has1 (exclusive scan); the priority order
//     is a total order (ties end at (mono hill, charge), unique), so any correct sort gives numpy's lexsort result;
//   * the isotope lookup keeps the explicit tie rule of _find_iso_binned (best cosine, then lowest global m/z rank).
// The sequential greedy acceptance and the feature table stay on the CPU (same numba/numpy code).
//
// Host API (ctypes, vulcan/ms1/gpu.py): a session holds device buffers between calls; all calls are synchronous on
// the session's own stream and return 0 or a cudaError_t code.  No torch dependency (static cudart).
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define CK(x) do { int _e = (int)(x); if (_e != 0) return _e; } while (0)
#define TPB 256
#define MAXI 16      // max_iso supported on the GPU (the CPU default is 10)

static inline unsigned nblk(int64_t n) { return (unsigned)((n + TPB - 1) / TPB); }

struct Session {
  cudaStream_t st;
  cudaEvent_t ev[8];
  cudaEvent_t evs;       // blocking-sync event: host waits sleep instead of spinning (core-s)
  float gpu_ms[4];       // hills, isotopes, peak_map, (unused)
  // hills (device)
  int64_t n_h, n_p, n_scans;
  int64_t* off; int32_t* pscan; double* pmz; float* pint; int64_t* psrc;
  double* hmz; int32_t* hstart; int32_t* hend; int32_t* hapex_scan; float* hapex_int; float* harea; int32_t* hlen;
  // isotopes (device)
  int64_t n_cand, m;
  int max_iso;
  int64_t* ci; int64_t* cz; int64_t* hills; float* tcos; float* acos; int32_t* nn; int64_t* order;
};

template <class T> static int dalloc(Session* S, T** p, size_t n) {
  *p = nullptr;
  if (n == 0) n = 1;
  return (int)cudaMallocAsync((void**)p, n * sizeof(T), S->st);
}
template <class T> static void dfree(Session* S, T** p) {
  if (*p) cudaFreeAsync(*p, S->st);
  *p = nullptr;
}

// ------------------------------------------------------------------------------------------------ CUB wrappers
template <class TI, class TO> static int ex_sum(Session* S, const TI* in, TO* out, int64_t n) {
  size_t tb = 0; void* tmp = nullptr;
  CK(cub::DeviceScan::ExclusiveSum(nullptr, tb, in, out, n, S->st));
  CK(cudaMallocAsync(&tmp, tb ? tb : 1, S->st));
  CK(cub::DeviceScan::ExclusiveSum(tmp, tb, in, out, n, S->st));
  CK(cudaFreeAsync(tmp, S->st));
  return 0;
}
template <class TI, class TO> static int in_sum(Session* S, const TI* in, TO* out, int64_t n) {
  size_t tb = 0; void* tmp = nullptr;
  CK(cub::DeviceScan::InclusiveSum(nullptr, tb, in, out, n, S->st));
  CK(cudaMallocAsync(&tmp, tb ? tb : 1, S->st));
  CK(cub::DeviceScan::InclusiveSum(tmp, tb, in, out, n, S->st));
  CK(cudaFreeAsync(tmp, S->st));
  return 0;
}
template <class K, class V> static int sort_pairs(Session* S, const K* kin, K* kout, const V* vin, V* vout, int64_t n,
                                                  int end_bit) {
  size_t tb = 0; void* tmp = nullptr;
  CK(cub::DeviceRadixSort::SortPairs(nullptr, tb, kin, kout, vin, vout, n, 0, end_bit, S->st));
  CK(cudaMallocAsync(&tmp, tb ? tb : 1, S->st));
  CK(cub::DeviceRadixSort::SortPairs(tmp, tb, kin, kout, vin, vout, n, 0, end_bit, S->st));
  CK(cudaFreeAsync(tmp, S->st));
  return 0;
}
static int ssync(Session* S) {   // wait for the stream without spinning a core
  CK(cudaEventRecord(S->evs, S->st));
  CK(cudaEventSynchronize(S->evs));
  return 0;
}
template <class T> static int d2h_scalar(Session* S, const T* d, T* h) {
  CK(ssync(S));                   // stream idle first: the pageable copy below then returns at once
  CK(cudaMemcpyAsync(h, d, sizeof(T), cudaMemcpyDeviceToHost, S->st));
  CK(ssync(S));
  return 0;
}

// ------------------------------------------------------------------------------------------------ hills kernels
__global__ void k_prep(int64_t n, int64_t n_scans, const int64_t* soff, const void* mz_in, int mz_f64,
                       const void* int_in, int int_f64, int32_t* scan, double* mz, float* inten) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  // scan = last s with soff[s] <= i (soff = searchsorted(scan, arange(n_scans+1)), scan ascending)
  int64_t lo = 0, hi = n_scans;           // soff[0] <= i always; find largest s in [0, n_scans) with soff[s] <= i
  while (hi - lo > 1) {
    int64_t mid = (lo + hi) >> 1;
    if (soff[mid] <= i) lo = mid; else hi = mid;
  }
  scan[i] = (int32_t)lo;
  mz[i] = mz_f64 ? ((const double*)mz_in)[i] : (double)((const float*)mz_in)[i];
  inten[i] = int_f64 ? __double2float_rn(((const double*)int_in)[i]) : ((const float*)int_in)[i];
}

// _link_adjacent: nearest peak of the previous MS1 scan (candidates j, j+1; strict < keeps j on ties)
__global__ void k_link_cand(int64_t n, const int32_t* scan, const int64_t* soff, const double* mz, double htol,
                            int32_t* cand, double* dist) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  cand[i] = -1;
  int s = scan[i];
  if (s < 1) return;
  int64_t pa = soff[s - 1], pb = soff[s];
  if (pb <= pa) return;
  double x = mz[i];
  // j = max(pa, (first index in [pa, pb) with mz > x) - 1): the CPU's forward pointer (mz ascending in a scan)
  int64_t lo = pa, hi = pb;
  while (lo < hi) {
    int64_t mid = (lo + hi) >> 1;
    if (mz[mid] <= x) lo = mid + 1; else hi = mid;
  }
  int64_t j = lo - 1;
  if (j < pa) j = pa;
  double bd = INFINITY;
  int64_t bc = -1;
  for (int64_t c = j; c <= j + 1; ++c) {
    if (c < pb) {
      double d = __dmul_rn(__ddiv_rn(fabs(__dsub_rn(mz[c], x)), x), 1e6);
      if (d < bd) { bd = d; bc = c; }
    }
  }
  if (bd <= htol) { cand[i] = (int32_t)bc; dist[i] = bd; }
}

__global__ void k_claim_d(int64_t n, const int32_t* cand, const double* dist, unsigned long long* claim) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n || cand[i] < 0) return;
  atomicMin(&claim[cand[i]], (unsigned long long)__double_as_longlong(dist[i]));   // dist >= +0: bits are monotone
}
__global__ void k_claim_q(int64_t n, const int32_t* cand, const double* dist, const unsigned long long* claim,
                          int32_t* claim_q) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n || cand[i] < 0) return;
  if ((unsigned long long)__double_as_longlong(dist[i]) == claim[cand[i]]) atomicMin(&claim_q[cand[i]], (int32_t)i);
}
// anc/rank init for list ranking: anc = pred (or self for a hill start), rank = 1 if linked
__global__ void k_pred(int64_t n, const int32_t* claim_q, int32_t* anc, int32_t* rk) {
  int64_t c = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (c >= n) return;
  int32_t q = claim_q[c];
  if (q != INT32_MAX) { anc[q] = (int32_t)c; rk[q] = 1; }
}
__global__ void k_init_self(int64_t n, int32_t* anc, int32_t* rk) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  anc[i] = (int32_t)i;
  rk[i] = 0;
}
__global__ void k_is_root(int64_t n, const int32_t* rk, int32_t* isr) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  isr[i] = rk[i] == 0;     // before ranking: rk = 1 iff linked
}
__global__ void k_jump(int64_t n, const int32_t* anc, const int32_t* rk, int32_t* anc2, int32_t* rk2, int* changed) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  int32_t a = anc[i];
  int32_t aa = anc[a];
  anc2[i] = aa;
  rk2[i] = rk[i] + rk[a];
  if (aa != a) *changed = 1;
}
__global__ void k_hill_of(int64_t n, const int32_t* anc, const int32_t* hid_root, int32_t* hill, int32_t* cnt) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  int32_t h = hid_root[anc[i]];
  hill[i] = h;
  atomicAdd(&cnt[h], 1);
}
__global__ void k_gather(int64_t n, const int32_t* hill, const int32_t* rk, const int64_t* off0, const int32_t* scan,
                         const double* mz, const float* inten, const int64_t* src, int32_t* g_scan, double* g_mz,
                         float* g_int, int64_t* g_src) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  int64_t p = off0[hill[i]] + rk[i];
  g_scan[p] = scan[i];
  g_mz[p] = mz[i];
  g_int[p] = inten[i];
  g_src[p] = src ? src[i] : i;
}

// _split_labels (per hill, sequential, same arithmetic)
__global__ void k_split(int64_t nh, const int64_t* off, const float* pint, double hvf, int minlh, int32_t* start,
                        double* sbuf, double* rbuf) {
  int64_t h = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (h >= nh) return;
  int64_t a = off[h], b = off[h + 1];
  int64_t L = b - a;
  if (L == 0) return;
  start[a] = 1;
  if (L < 2 * (int64_t)minlh + 1 || L < 5) return;
  double* s = sbuf + a;
  double* r = rbuf + a;
  for (int64_t i = 0; i < L; ++i) {
    int64_t lo = i - 1 > 0 ? i - 1 : 0;
    int64_t hi = i + 1 < L - 1 ? i + 1 : L - 1;
    double acc = 0.0;
    for (int64_t j = lo; j <= hi; ++j) acc = __dadd_rn(acc, (double)pint[a + j]);
    s[i] = __ddiv_rn(acc, (double)(hi - lo + 1));
  }
  double rm = 0.0;                        // rmax[L] = 0; rmax[i] = max(rmax[i+1], s[i])
  for (int64_t i = L - 1; i >= 0; --i) {
    if (s[i] > rm) rm = s[i];
    r[i] = rm;
  }
  int64_t seg = 0;
  double lmax = s[0];
  for (int64_t i = 1; i < L; ++i) {
    if (i - seg >= minlh && L - i >= minlh && i < L - 1) {
      if (s[i] < s[i - 1] && s[i] <= s[i + 1]) {
        double sh = __dmul_rn(s[i], hvf);
        if (sh < lmax && sh < r[i + 1]) {
          start[a + i] = 1;
          seg = i;
          lmax = s[i];
          continue;
        }
      }
    }
    if (s[i] > lmax) lmax = s[i];
  }
}
__global__ void k_lab_start(int64_t n, const int32_t* start, const int32_t* lab, int64_t* ls) {
  int64_t p = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (p >= n) return;
  if (start[p]) ls[lab[p]] = p;
}
__global__ void k_keep(int64_t nl, const int64_t* ls, int minlh, int32_t* keep, int64_t* kcnt) {
  int64_t l = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (l >= nl) return;
  int64_t c = ls[l + 1] - ls[l];
  int k = c >= minlh;
  keep[l] = k;
  kcnt[l] = k ? c : 0;
}
__global__ void k_compact(int64_t n, const int32_t* lab, const int32_t* keep, const int64_t* ls, const int64_t* noff,
                          const int32_t* g_scan, const double* g_mz, const float* g_int, const int64_t* g_src,
                          int32_t* o_scan, double* o_mz, float* o_int, int64_t* o_src) {
  int64_t p = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (p >= n) return;
  int32_t l = lab[p];
  if (!keep[l]) return;
  int64_t o = noff[l] + (p - ls[l]);
  o_scan[o] = g_scan[p];
  o_mz[o] = g_mz[p];
  o_int[o] = g_int[p];
  o_src[o] = g_src[p];
}
__global__ void k_new_off(int64_t nl, const int32_t* keep, const int32_t* nid, const int64_t* noff, int64_t* off) {
  int64_t l = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (l >= nl) return;
  if (keep[l]) off[nid[l]] = noff[l];
}
// _summaries
__global__ void k_summaries(int64_t nh, const int64_t* off, const int32_t* pscan, const double* pmz, const float* pint,
                            const double* rt, double* hmz, int32_t* hs, int32_t* he, int32_t* hap, float* hai,
                            float* harea, int32_t* hlen) {
  int64_t h = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (h >= nh) return;
  int64_t a = off[h], b = off[h + 1];
  double sw = 0.0, smz = 0.0, best = -1.0, area = 0.0;
  int64_t bi = a;
  for (int64_t i = a; i < b; ++i) {
    float w = pint[i];
    sw = __dadd_rn(sw, (double)w);
    smz = __dadd_rn(smz, __dmul_rn((double)w, pmz[i]));
    if ((double)w > best) { best = (double)w; bi = i; }
    if (i > a) {
      double s2 = (double)__fadd_rn(pint[i], pint[i - 1]);                 // numba: f32 + f32 -> f32
      double dt = __dsub_rn(rt[pscan[i]], rt[pscan[i - 1]]);
      area = __dadd_rn(area, __dmul_rn(__dmul_rn(0.5, s2), dt));
    }
  }
  hmz[h] = sw > 0 ? __ddiv_rn(smz, sw) : pmz[a];
  hs[h] = pscan[a];
  he[h] = pscan[b - 1];
  hap[h] = pscan[bi];
  hai[h] = __double2float_rn(best);
  harea[h] = __double2float_rn(area);
  hlen[h] = (int32_t)(b - a);
}

// ------------------------------------------------------------------------------------------------ isotope kernels
__global__ void k_iota32(int64_t n, int32_t* v) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i < n) v[i] = (int32_t)i;
}
__global__ void k_norm2(int64_t nh, const int64_t* off, const float* pint, double* out) {
  int64_t h = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (h >= nh) return;
  double s = 0.0;
  for (int64_t q = off[h]; q < off[h + 1]; ++q) s = __dadd_rn(s, (double)__fmul_rn(pint[q], pint[q]));
  out[h] = s;
}
__global__ void k_nblocks(int64_t nh, const int32_t* sidx, const int32_t* hs, const int32_t* he, int bs, int64_t* nb) {
  int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (q >= nh) return;
  int32_t h = sidx[q];
  nb[q] = he[h] / bs - hs[h] / bs + 1;
}
__global__ void k_expand(int64_t nh, const int32_t* sidx, const int32_t* hs, const int64_t* eoff, int bs, int32_t* key,
                         int32_t* val, int32_t* bcnt) {
  int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (q >= nh) return;
  int32_t b0 = hs[sidx[q]] / bs;
  for (int64_t e = eoff[q], k = 0; e < eoff[q + 1]; ++e, ++k) {
    key[e] = b0 + (int32_t)k;
    val[e] = (int32_t)q;
    atomicAdd(&bcnt[b0 + k], 1);
  }
}
__global__ void k_bin_fill(int64_t ne, const int32_t* bq, const double* smz, const int32_t* sidx, const int32_t* hs,
                           const int32_t* he, double* b_mz, int32_t* b_h, int32_t* b_lo, int32_t* b_hi) {
  int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (e >= ne) return;
  int32_t q = bq[e];
  int32_t h = sidx[q];
  b_mz[e] = smz[q];
  b_h[e] = h;
  b_lo[e] = hs[h];
  b_hi[e] = he[h];
}

struct Iso {
  const int64_t* off; const int32_t* hs; const int32_t* he; const float* pint; const double* norm2;
  const double* hmz; const float* hai;
  const int64_t* boff; const double* b_mz; const int32_t* b_pos; const int32_t* b_h; const int32_t* b_lo;
  const int32_t* b_hi;
  int bs;
  double itol, mtc, min_avg_cos, min_expl;
  int min_ov;
  const double* table; int rows, cols;
};

__device__ __forceinline__ void trace_cos(const Iso& I, int32_t a, int32_t b, double& c, int64_t& ov) {
  int32_t s0 = max(I.hs[a], I.hs[b]);
  int32_t s1 = min(I.he[a], I.he[b]);
  if (s1 < s0) { c = -1.0; ov = 0; return; }
  int64_t i = I.off[a] + (s0 - I.hs[a]);
  int64_t j = I.off[b] + (s0 - I.hs[b]);
  double sab = 0.0;
  for (int64_t k = 0; k < (int64_t)(s1 - s0) + 1; ++k)
    sab = __dadd_rn(sab, (double)__fmul_rn(I.pint[i + k], I.pint[j + k]));
  double saa = I.norm2[a], sbb = I.norm2[b];
  ov = (int64_t)(s1 - s0) + 1;
  if (saa <= 0 || sbb <= 0) { c = 0.0; return; }
  c = __ddiv_rn(sab, __dsqrt_rn(__dmul_rn(saa, sbb)));
}

// _find_iso_binned
__device__ int32_t find_iso(const Iso& I, int32_t mono, double target, bool first_only, double& bestc_out) {
  double lo = __dmul_rn(target, __dsub_rn(1.0, __dmul_rn(I.itol, 1e-6)));
  double hi = __dmul_rn(target, __dadd_rn(1.0, __dmul_rn(I.itol, 1e-6)));
  int32_t m0 = I.hs[mono], m1 = I.he[mono];
  int32_t b0 = m0 / I.bs;
  int32_t best = -1;
  double bestc = -2.0;
  int64_t bestp = (int64_t)1 << 62;
  for (int32_t b = b0; b <= m1 / I.bs; ++b) {
    int64_t a = I.boff[b], e = I.boff[b + 1];
    int64_t l = a, r = e;                  // searchsorted(b_mz[a:e], lo, 'left')
    while (l < r) {
      int64_t mid = (l + r) >> 1;
      if (I.b_mz[mid] < lo) l = mid + 1; else r = mid;
    }
    for (int64_t p = l; p < e && I.b_mz[p] <= hi; ++p) {
      if (I.b_lo[p] > m1 || I.b_hi[p] < m0) continue;
      int32_t fb = I.b_lo[p] / I.bs;
      if ((fb > b0 ? fb : b0) != b) continue;
      int32_t h = I.b_h[p];
      double c;
      int64_t ov;
      trace_cos(I, mono, h, c, ov);
      if (ov >= 1) {
        if (ov < I.min_ov) c = c > I.mtc ? c : I.mtc;
        int64_t bp = I.b_pos[p];
        if (c >= I.mtc && (c > bestc || (c == bestc && bp < bestp))) {
          bestc = c;
          best = h;
          bestp = bp;
          if (first_only) { bestc_out = bestc; return best; }
        }
      }
    }
  }
  bestc_out = bestc;
  return best;
}

__device__ void avg_cos(const Iso& I, double mass, const double* Iv, int L, double& c, double& expl) {
  double r = __ddiv_rn(mass, 10.0);
  int64_t r0 = (int64_t)r;
  double fr;
  if (r0 >= I.rows - 1) { r0 = I.rows - 2; fr = 1.0; } else fr = __dsub_rn(r, (double)r0);
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  expl = 0.0;
  for (int k = 0; k < L; ++k) {
    double t = k < I.cols ? __dadd_rn(__dmul_rn(__dsub_rn(1.0, fr), I.table[r0 * I.cols + k]),
                                      __dmul_rn(fr, I.table[(r0 + 1) * I.cols + k])) : 0.0;
    sab = __dadd_rn(sab, __dmul_rn(t, Iv[k]));
    saa = __dadd_rn(saa, __dmul_rn(Iv[k], Iv[k]));
    sbb = __dadd_rn(sbb, __dmul_rn(t, t));
    expl = __dadd_rn(expl, t);
  }
  if (saa <= 0 || sbb <= 0) { c = 0.0; return; }
  c = __ddiv_rn(sab, __dsqrt_rn(__dmul_rn(saa, sbb)));
}
__device__ int best_prefix(const Iso& I, double mass, const double* Iv, int n, double& cout_) {
  for (int L = n; L > 1; --L) {
    double c, e;
    avg_cos(I, mass, Iv, L, c, e);
    if (c >= I.min_avg_cos && e >= I.min_expl) { cout_ = c; return L; }
  }
  cout_ = 0.0;
  return 1;
}

__global__ void k_count(Iso I, int64_t nh, const int32_t* sidx, const int64_t* charges, int nz, uint8_t* has1) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= nh * nz) return;
  int64_t q = t / nz;
  int zi = (int)(t - q * nz);
  int32_t i = sidx[q];                         // m/z order: neighbouring threads touch neighbouring bins
  double target = __dadd_rn(I.hmz[i], __ddiv_rn(1.0033548378, (double)charges[zi]));
  double c;
  has1[(int64_t)i * nz + zi] = find_iso(I, i, target, true, c) >= 0;
}
__global__ void k_cand(int64_t nh, int nz, const uint8_t* has1, const int32_t* cidx, const int64_t* charges,
                       int64_t* ci, int64_t* cz) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= nh * nz) return;
  if (!has1[t]) return;
  int64_t k = cidx[t];
  ci[k] = t / nz;
  cz[k] = charges[t % nz];
}
__global__ void k_fill(Iso I, int64_t nh, const int32_t* sidx, int nz, const uint8_t* has1, const int32_t* cidx,
                       const int64_t* charges, int max_iso, int64_t* hills, float* tcos, float* acos, int32_t* nn) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= nh * nz) return;
  int64_t q = t / nz;
  int zi = (int)(t - q * nz);
  int32_t i = sidx[q];
  int64_t f = (int64_t)i * nz + zi;
  if (!has1[f]) return;
  int64_t cq = cidx[f];
  int64_t z = charges[zi];
  int32_t hb[MAXI];
  double cb[MAXI];
  hb[0] = i;
  cb[0] = 1.0;
  int n = 1;
  for (int k = 1; k < max_iso; ++k) {
    double target = __dadd_rn(I.hmz[i], __ddiv_rn(__dmul_rn((double)k, 1.0033548378), (double)z));
    double c;
    int32_t j = find_iso(I, i, target, false, c);
    if (j < 0) break;
    hb[n] = j;
    cb[n] = c;
    ++n;
  }
  int L = 1;
  double ac = 0.0;
  if (n >= 2) {
    double Iv[MAXI];
    for (int k = 0; k < n; ++k) Iv[k] = (double)I.hai[hb[k]];
    double mass = __dmul_rn(__dsub_rn(I.hmz[i], 1.007276466812), (double)z);
    L = best_prefix(I, mass, Iv, n, ac);
  }
  nn[cq] = L;
  acos[cq] = __double2float_rn(ac);
  for (int k = 0; k < L; ++k) {
    hills[cq * max_iso + k] = hb[k];
    tcos[cq * max_iso + k] = __double2float_rn(cb[k]);
  }
}
__global__ void k_fill_init(int64_t n, int64_t* hills, float* tcos) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= n) return;
  hills[t] = -1;
  tcos[t] = 0.0f;
}
__global__ void k_flag_ok(int64_t m, const int32_t* nn, int32_t* flag) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t < m) flag[t] = nn[t] >= 2;
}
__global__ void k_grid(int64_t m, const int32_t* flag, const int64_t* ci, const int64_t* cz, const int32_t* nn, int zw,
                       int32_t* grid) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= m || !flag[t]) return;
  grid[ci[t] * zw + cz[t]] = nn[t];             // (mono, z) unique: max == the value
}
// charge_ban on the nn >= 2 candidates; flag2 = ok && kept
__global__ void k_ban(int64_t m, const int32_t* flag, const int64_t* ci, const int64_t* cz, const int32_t* nn, int zw,
                      const int32_t* grid, int use_ban, int32_t* flag2) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= m) return;
  int keep = flag[t];
  if (keep && use_ban) {
    int64_t z = cz[t];
    int zmax = zw - 1;
    if (2 * z <= zmax) {
      int32_t hi = 0;
      for (int64_t kz = 2 * z; kz <= zmax; kz += z) {
        int32_t g = grid[ci[t] * zw + kz];
        hi = g > hi ? g : hi;
      }
      if (hi >= nn[t]) keep = 0;
    }
  }
  flag2[t] = keep;
}
__global__ void k_compact_cand(int64_t m, const int32_t* flag, const int64_t* pos, int max_iso, const int64_t* ci,
                               const int64_t* cz, const int64_t* hills, const float* tcos, const float* acos,
                               const int32_t* nn, int64_t* o_ci, int64_t* o_cz, int64_t* o_hills, float* o_tcos,
                               float* o_acos, int32_t* o_nn) {
  int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (t >= m || !flag[t]) return;
  int64_t o = pos[t];
  o_ci[o] = ci[t];
  o_cz[o] = cz[t];
  o_acos[o] = acos[t];
  o_nn[o] = nn[t];
  for (int k = 0; k < max_iso; ++k) {
    o_hills[o * max_iso + k] = hills[t * max_iso + k];
    o_tcos[o * max_iso + k] = tcos[t * max_iso + k];
  }
}
struct PrioLess {   // np.lexsort((cz, ci, -apex[ci], -acos, -nn)): nn desc, acos desc, apex desc, ci asc, cz asc
  const int32_t* nn; const float* acos; const float* hai; const int64_t* ci; const int64_t* cz;
  __device__ bool operator()(int64_t a, int64_t b) const {
    if (nn[a] != nn[b]) return nn[a] > nn[b];
    if (acos[a] != acos[b]) return acos[a] > acos[b];
    float xa = hai[ci[a]], xb = hai[ci[b]];
    if (xa != xb) return xa > xb;
    if (ci[a] != ci[b]) return ci[a] < ci[b];
    return cz[a] < cz[b];
  }
};
__global__ void k_widen(int64_t n, const uint8_t* a, int32_t* w) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i < n) w[i] = a[i];
}
__global__ void k_dec(int64_t n, int32_t* a) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i < n) a[i] -= 1;
}
__global__ void k_iota64(int64_t n, int64_t* v) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i < n) v[i] = i;
}

// peak -> feature / isotope index (features.py _feature_table tail)
__global__ void k_peak_map(int64_t nh, const int64_t* off, const int64_t* psrc, const int32_t* hill_feat,
                           const int8_t* hill_iso, int32_t* peak_feature, int8_t* peak_iso) {
  int64_t h = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (h >= nh) return;
  int32_t f = hill_feat[h];
  int8_t k = hill_iso[h];
  for (int64_t p = off[h]; p < off[h + 1]; ++p) {
    peak_feature[psrc[p]] = f;
    peak_iso[psrc[p]] = k;
  }
}
__global__ void k_fill_i32(int64_t n, int32_t* a, int32_t v) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i < n) a[i] = v;
}

// ================================================================================================= host API
extern "C" {

int ms1g_create(void** out) {
  Session* S = new Session();
  memset(S, 0, sizeof(Session));
  cudaError_t e = cudaStreamCreateWithFlags(&S->st, cudaStreamNonBlocking);
  if (e != cudaSuccess) { delete S; return (int)e; }
  for (int k = 0; k < 8; ++k) CK(cudaEventCreate(&S->ev[k]));
  CK(cudaEventCreateWithFlags(&S->evs, cudaEventBlockingSync | cudaEventDisableTiming));
  *out = S;
  return 0;
}

static void free_hills(Session* S) {
  dfree(S, &S->off); dfree(S, &S->pscan); dfree(S, &S->pmz); dfree(S, &S->pint); dfree(S, &S->psrc);
  dfree(S, &S->hmz); dfree(S, &S->hstart); dfree(S, &S->hend); dfree(S, &S->hapex_scan); dfree(S, &S->hapex_int);
  dfree(S, &S->harea); dfree(S, &S->hlen);
}
static void free_iso(Session* S) {
  dfree(S, &S->ci); dfree(S, &S->cz); dfree(S, &S->hills); dfree(S, &S->tcos); dfree(S, &S->acos); dfree(S, &S->nn);
  dfree(S, &S->order);
}

int ms1g_release(void* h) {    // free the run's device buffers (session stays usable)
  Session* S = (Session*)h;
  free_hills(S);
  free_iso(S);
  CK(ssync(S));
  return 0;
}
int ms1g_destroy(void* h) {
  Session* S = (Session*)h;
  ms1g_release(h);
  for (int k = 0; k < 8; ++k) cudaEventDestroy(S->ev[k]);
  cudaEventDestroy(S->evs);
  cudaStreamDestroy(S->st);
  delete S;
  return 0;
}
int ms1g_set_pool_threshold(uint64_t bytes) {   // keep up to `bytes` of freed memory cached in the default pool
  int dev;
  CK(cudaGetDevice(&dev));
  cudaMemPool_t pool;
  CK(cudaDeviceGetDefaultMemPool(&pool, dev));
  CK(cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &bytes));
  return 0;
}
int ms1g_mem_info(uint64_t* free_b, uint64_t* total_b) {
  size_t f = 0, t = 0;
  CK(cudaMemGetInfo(&f, &t));
  *free_b = f;
  *total_b = t;
  return 0;
}
int ms1g_trim(void) {          // return the default pool's cached memory to the driver (idle decode workers)
  int dev;
  CK(cudaGetDevice(&dev));
  cudaMemPool_t pool;
  CK(cudaDeviceGetDefaultMemPool(&pool, dev));
  CK(cudaDeviceSynchronize());
  CK(cudaMemPoolTrimTo(pool, 0));
  return 0;
}
void ms1g_gpu_ms(void* h, float* out) {
  Session* S = (Session*)h;
  for (int k = 0; k < 4; ++k) out[k] = S->gpu_ms[k];
}

// Hills from flat MS1 peaks (scan ascending, mz ascending within a scan; soff = searchsorted(scan, arange(n_scans+1))).
// mz/inten: float32 or float64 host arrays (flags); src: int64 host array or NULL (= arange).  Returns n_h, n_p.
int ms1g_hills(void* h, int64_t n, int64_t n_scans, const int64_t* soff_h, const void* mz_h, int mz_f64,
               const void* int_h, int int_f64, const int64_t* src_h, const double* rt_h, int64_t n_rt, double htol,
               int minlh, double hvf, int64_t* out_nh, int64_t* out_np) {
  Session* S = (Session*)h;
  free_hills(S);
  free_iso(S);
  if (n >= INT32_MAX) return -2;
  S->n_scans = n_scans;
  CK(cudaEventRecord(S->ev[0], S->st));
  int64_t *soff, *src = nullptr;
  void *mz_in, *int_in;
  double* rt;
  CK(dalloc(S, &soff, n_scans + 1));
  CK(dalloc(S, (char**)&mz_in, n * (mz_f64 ? 8 : 4)));
  CK(dalloc(S, (char**)&int_in, n * (int_f64 ? 8 : 4)));
  CK(dalloc(S, &rt, n_rt));
  CK(cudaMemcpyAsync(soff, soff_h, (n_scans + 1) * 8, cudaMemcpyHostToDevice, S->st));
  CK(cudaMemcpyAsync(mz_in, mz_h, n * (mz_f64 ? 8 : 4), cudaMemcpyHostToDevice, S->st));
  CK(cudaMemcpyAsync(int_in, int_h, n * (int_f64 ? 8 : 4), cudaMemcpyHostToDevice, S->st));
  CK(cudaMemcpyAsync(rt, rt_h, n_rt * 8, cudaMemcpyHostToDevice, S->st));
  if (src_h) {
    CK(dalloc(S, &src, n));
    CK(cudaMemcpyAsync(src, src_h, n * 8, cudaMemcpyHostToDevice, S->st));
  }
  int32_t* scan; double* mz; float* inten;
  CK(dalloc(S, &scan, n)); CK(dalloc(S, &mz, n)); CK(dalloc(S, &inten, n));
  if (n) k_prep<<<nblk(n), TPB, 0, S->st>>>(n, n_scans, soff, mz_in, mz_f64, int_in, int_f64, scan, mz, inten);
  dfree(S, (char**)&mz_in); dfree(S, (char**)&int_in);
  // ---- link
  int32_t *cand, *claim_q, *anc, *rk, *anc2, *rk2, *isr, *hid_root;
  double* dist;
  unsigned long long* claim;
  CK(dalloc(S, &cand, n)); CK(dalloc(S, &dist, n)); CK(dalloc(S, &claim, n)); CK(dalloc(S, &claim_q, n));
  CK(cudaMemsetAsync(claim, 0xff, n * 8, S->st));
  if (n) {
    k_link_cand<<<nblk(n), TPB, 0, S->st>>>(n, scan, soff, mz, htol, cand, dist);
    k_claim_d<<<nblk(n), TPB, 0, S->st>>>(n, cand, dist, claim);
    k_fill_i32<<<nblk(n), TPB, 0, S->st>>>(n, claim_q, INT32_MAX);
    k_claim_q<<<nblk(n), TPB, 0, S->st>>>(n, cand, dist, claim, claim_q);
  }
  dfree(S, &cand); dfree(S, &dist); dfree(S, &claim);
  CK(dalloc(S, &anc, n)); CK(dalloc(S, &rk, n)); CK(dalloc(S, &anc2, n)); CK(dalloc(S, &rk2, n));
  CK(dalloc(S, &isr, n)); CK(dalloc(S, &hid_root, n + 1));
  if (n) {
    k_init_self<<<nblk(n), TPB, 0, S->st>>>(n, anc, rk);
    k_pred<<<nblk(n), TPB, 0, S->st>>>(n, claim_q, anc, rk);
    k_is_root<<<nblk(n), TPB, 0, S->st>>>(n, rk, isr);
  }
  dfree(S, &claim_q);
  // ---- list ranking (Wyllie): anc -> chain root, rk -> depth (= rank within the hill, index order)
  int* d_changed;
  CK(dalloc(S, &d_changed, 1));
  for (int it = 0; it < 64 && n; ++it) {
    CK(cudaMemsetAsync(d_changed, 0, sizeof(int), S->st));
    k_jump<<<nblk(n), TPB, 0, S->st>>>(n, anc, rk, anc2, rk2, d_changed);
    int32_t* t;
    t = anc; anc = anc2; anc2 = t;
    t = rk; rk = rk2; rk2 = t;
    int changed;
    CK(d2h_scalar(S, d_changed, &changed));
    if (!changed) break;
  }
  dfree(S, &d_changed); dfree(S, &anc2); dfree(S, &rk2);
  // ---- hill ids (rank of the root by index), counting-sort placement
  CK(ex_sum(S, isr, hid_root, n));
  int32_t last_isr = 0, last_hid = 0;
  if (n) { CK(d2h_scalar(S, isr + n - 1, &last_isr)); CK(d2h_scalar(S, hid_root + n - 1, &last_hid)); }
  int64_t nh0 = n ? (int64_t)last_hid + last_isr : 0;
  int32_t *hill, *cnt;
  int64_t* off0;
  CK(dalloc(S, &hill, n)); CK(dalloc(S, &cnt, nh0 + 1)); CK(dalloc(S, &off0, nh0 + 1));
  CK(cudaMemsetAsync(cnt, 0, (nh0 + 1) * 4, S->st));
  if (n) k_hill_of<<<nblk(n), TPB, 0, S->st>>>(n, anc, hid_root, hill, cnt);
  CK(ex_sum(S, cnt, off0, nh0 + 1));
  dfree(S, &isr); dfree(S, &hid_root); dfree(S, &anc); dfree(S, &cnt);
  int32_t* g_scan; double* g_mz; float* g_int; int64_t* g_src;
  CK(dalloc(S, &g_scan, n)); CK(dalloc(S, &g_mz, n)); CK(dalloc(S, &g_int, n)); CK(dalloc(S, &g_src, n));
  if (n) k_gather<<<nblk(n), TPB, 0, S->st>>>(n, hill, rk, off0, scan, mz, inten, src, g_scan, g_mz, g_int, g_src);
  dfree(S, &hill); dfree(S, &rk); dfree(S, &scan); dfree(S, &mz); dfree(S, &inten); dfree(S, &src);
  dfree(S, &soff);
  // ---- valley split -> labels
  int32_t *start, *lab;
  double *sbuf, *rbuf;
  CK(dalloc(S, &start, n)); CK(dalloc(S, &lab, n)); CK(dalloc(S, &sbuf, n)); CK(dalloc(S, &rbuf, n));
  CK(cudaMemsetAsync(start, 0, n * 4, S->st));
  if (nh0) k_split<<<nblk(nh0), TPB, 0, S->st>>>(nh0, off0, g_int, hvf, minlh, start, sbuf, rbuf);
  dfree(S, &sbuf); dfree(S, &rbuf); dfree(S, &off0);
  CK(in_sum(S, start, lab, n));                       // lab = cumsum(start) - 1
  if (n) k_dec<<<nblk(n), TPB, 0, S->st>>>(n, lab);
  int32_t last_lab = -1;
  if (n) CK(d2h_scalar(S, lab + n - 1, &last_lab));
  int64_t nl = (int64_t)last_lab + 1;
  // ---- drop hills shorter than minlh, compact
  int64_t *ls, *kcnt, *noff;
  int32_t *keep, *nid;
  CK(dalloc(S, &ls, nl + 1)); CK(dalloc(S, &kcnt, nl + 1)); CK(dalloc(S, &noff, nl + 1));
  CK(dalloc(S, &keep, nl + 1)); CK(dalloc(S, &nid, nl + 1));
  if (n) k_lab_start<<<nblk(n), TPB, 0, S->st>>>(n, start, lab, ls);
  CK(cudaMemcpyAsync(ls + nl, &n, 8, cudaMemcpyHostToDevice, S->st));
  CK(cudaMemsetAsync(keep + nl, 0, 4, S->st));
  CK(cudaMemsetAsync(kcnt + nl, 0, 8, S->st));
  if (nl) k_keep<<<nblk(nl), TPB, 0, S->st>>>(nl, ls, minlh, keep, kcnt);
  CK(ex_sum(S, keep, nid, nl + 1));
  CK(ex_sum(S, kcnt, noff, nl + 1));
  int32_t nh_i;
  int64_t np_;
  CK(d2h_scalar(S, nid + nl, &nh_i));
  CK(d2h_scalar(S, noff + nl, &np_));
  int64_t nh = nh_i;
  CK(dalloc(S, &S->off, nh + 1)); CK(dalloc(S, &S->pscan, np_)); CK(dalloc(S, &S->pmz, np_));
  CK(dalloc(S, &S->pint, np_)); CK(dalloc(S, &S->psrc, np_));
  if (n) k_compact<<<nblk(n), TPB, 0, S->st>>>(n, lab, keep, ls, noff, g_scan, g_mz, g_int, g_src, S->pscan, S->pmz,
                                               S->pint, S->psrc);
  if (nl) k_new_off<<<nblk(nl), TPB, 0, S->st>>>(nl, keep, nid, noff, S->off);
  CK(cudaMemcpyAsync(S->off + nh, &np_, 8, cudaMemcpyHostToDevice, S->st));
  CK(ssync(S));     // &np_ is a host stack value
  dfree(S, &ls); dfree(S, &kcnt); dfree(S, &noff); dfree(S, &keep); dfree(S, &nid); dfree(S, &start); dfree(S, &lab);
  dfree(S, &g_scan); dfree(S, &g_mz); dfree(S, &g_int); dfree(S, &g_src);
  // ---- summaries
  CK(dalloc(S, &S->hmz, nh)); CK(dalloc(S, &S->hstart, nh)); CK(dalloc(S, &S->hend, nh));
  CK(dalloc(S, &S->hapex_scan, nh)); CK(dalloc(S, &S->hapex_int, nh)); CK(dalloc(S, &S->harea, nh));
  CK(dalloc(S, &S->hlen, nh));
  if (nh) k_summaries<<<nblk(nh), TPB, 0, S->st>>>(nh, S->off, S->pscan, S->pmz, S->pint, rt, S->hmz, S->hstart,
                                                   S->hend, S->hapex_scan, S->hapex_int, S->harea, S->hlen);
  dfree(S, &rt);
  CK(cudaEventRecord(S->ev[1], S->st));
  CK(ssync(S));
  CK(cudaGetLastError());
  cudaEventElapsedTime(&S->gpu_ms[0], S->ev[0], S->ev[1]);
  S->n_h = nh;
  S->n_p = np_;
  *out_nh = nh;
  *out_np = np_;
  return 0;
}

// Copy the hills to host arrays (NULL pointers are skipped).
int ms1g_hills_get(void* h, int64_t* off, int32_t* pscan, double* pmz, float* pint, int64_t* psrc, double* hmz,
                   int32_t* hs, int32_t* he, int32_t* hap, float* hai, float* harea, int32_t* hlen) {
  Session* S = (Session*)h;
  int64_t nh = S->n_h, np_ = S->n_p;
#define G(dst, src, cnt) if (dst && (cnt)) CK(cudaMemcpyAsync(dst, src, (cnt) * sizeof(*(dst)), cudaMemcpyDeviceToHost, S->st))
  G(off, S->off, nh + 1); G(pscan, S->pscan, np_); G(pmz, S->pmz, np_); G(pint, S->pint, np_); G(psrc, S->psrc, np_);
  G(hmz, S->hmz, nh); G(hs, S->hstart, nh); G(he, S->hend, nh); G(hap, S->hapex_scan, nh); G(hai, S->hapex_int, nh);
  G(harea, S->harea, nh); G(hlen, S->hlen, nh);
#undef G
  CK(ssync(S));
  return 0;
}

// Isotope grouping up to (and excluding) the greedy pass.  Output m candidates (nn >= 2, after the charge ban) in
// np.nonzero(has1) order + the priority order.  Requires contiguous hills (checked by the caller).
int ms1g_isotopes(void* h, const int64_t* charges_h, int nz, double itol, int max_iso, double mtc, int min_ov,
                  const double* table_h, int rows, int cols, double min_avg_cos, double min_expl, int use_ban,
                  int scan_bin, int64_t* out_ncand, int64_t* out_m) {
  Session* S = (Session*)h;
  free_iso(S);
  if (max_iso > MAXI || max_iso < 1) return -3;
  int64_t nh = S->n_h;
  CK(cudaEventRecord(S->ev[2], S->st));
  int64_t* charges;
  double* table;
  CK(dalloc(S, &charges, nz)); CK(dalloc(S, &table, (int64_t)rows * cols));
  CK(cudaMemcpyAsync(charges, charges_h, nz * 8, cudaMemcpyHostToDevice, S->st));
  CK(cudaMemcpyAsync(table, table_h, (int64_t)rows * cols * 8, cudaMemcpyHostToDevice, S->st));
  int64_t zmax = 0;
  for (int k = 0; k < nz; ++k) zmax = charges_h[k] > zmax ? charges_h[k] : zmax;
  // ---- m/z order (stable)
  int32_t *iota, *sidx;
  double* smz;
  CK(dalloc(S, &iota, nh)); CK(dalloc(S, &sidx, nh)); CK(dalloc(S, &smz, nh));
  if (nh) k_iota32<<<nblk(nh), TPB, 0, S->st>>>(nh, iota);
  CK(sort_pairs(S, S->hmz, smz, iota, sidx, nh, 64));
  dfree(S, &iota);
  double* norm2;
  CK(dalloc(S, &norm2, nh));
  if (nh) k_norm2<<<nblk(nh), TPB, 0, S->st>>>(nh, S->off, S->pint, norm2);
  // ---- scan-block bins (stable by block, m/z rank inside)
  int bs = scan_bin;
  int64_t nbin = S->n_scans / bs + 1;
  int64_t *nbq, *eoff;
  CK(dalloc(S, &nbq, nh + 1)); CK(dalloc(S, &eoff, nh + 1));
  CK(cudaMemsetAsync(nbq + nh, 0, 8, S->st));
  if (nh) k_nblocks<<<nblk(nh), TPB, 0, S->st>>>(nh, sidx, S->hstart, S->hend, bs, nbq);
  CK(ex_sum(S, nbq, eoff, nh + 1));
  int64_t ne;
  CK(d2h_scalar(S, eoff + nh, &ne));
  if (ne >= INT32_MAX) return -4;
  int32_t *ekey, *eval, *skey, *bq, *bcnt;
  int64_t* boff;
  CK(dalloc(S, &ekey, ne)); CK(dalloc(S, &eval, ne)); CK(dalloc(S, &skey, ne)); CK(dalloc(S, &bq, ne));
  CK(dalloc(S, &bcnt, nbin + 1)); CK(dalloc(S, &boff, nbin + 1));
  CK(cudaMemsetAsync(bcnt, 0, (nbin + 1) * 4, S->st));
  if (nh) k_expand<<<nblk(nh), TPB, 0, S->st>>>(nh, sidx, S->hstart, eoff, bs, ekey, eval, bcnt);
  int end_bit = 1;
  while (((int64_t)1 << end_bit) <= nbin) ++end_bit;
  CK(sort_pairs(S, ekey, skey, eval, bq, ne, end_bit));
  CK(ex_sum(S, bcnt, boff, nbin + 1));
  dfree(S, &ekey); dfree(S, &eval); dfree(S, &skey); dfree(S, &nbq); dfree(S, &eoff); dfree(S, &bcnt);
  double* b_mz;
  int32_t *b_h, *b_lo, *b_hi;
  CK(dalloc(S, &b_mz, ne)); CK(dalloc(S, &b_h, ne)); CK(dalloc(S, &b_lo, ne)); CK(dalloc(S, &b_hi, ne));
  if (ne) k_bin_fill<<<nblk(ne), TPB, 0, S->st>>>(ne, bq, smz, sidx, S->hstart, S->hend, b_mz, b_h, b_lo, b_hi);
  Iso I;
  I.off = S->off; I.hs = S->hstart; I.he = S->hend; I.pint = S->pint; I.norm2 = norm2; I.hmz = S->hmz;
  I.hai = S->hapex_int; I.boff = boff; I.b_mz = b_mz; I.b_pos = bq; I.b_h = b_h; I.b_lo = b_lo; I.b_hi = b_hi;
  I.bs = bs; I.itol = itol; I.mtc = mtc; I.min_avg_cos = min_avg_cos; I.min_expl = min_expl; I.min_ov = min_ov;
  I.table = table; I.rows = rows; I.cols = cols;
  // ---- first-isotope pre-screen
  int64_t nt = nh * nz;
  uint8_t* has1;
  int32_t* cidx;
  CK(dalloc(S, &has1, nt)); CK(dalloc(S, &cidx, nt + 1));
  if (nt) k_count<<<nblk(nt), TPB, 0, S->st>>>(I, nh, sidx, charges, nz, has1);
  {
    // exclusive scan of has1 (widened to int32) -> candidate index in row-major nonzero order
    int32_t* w;
    CK(dalloc(S, &w, nt));
    if (nt) k_widen<<<nblk(nt), TPB, 0, S->st>>>(nt, has1, w);
    CK(ex_sum(S, w, cidx, nt));
    dfree(S, &w);
  }
  int32_t last_c = 0;
  uint8_t last_h = 0;
  if (nt) { CK(d2h_scalar(S, cidx + nt - 1, &last_c)); CK(d2h_scalar(S, has1 + nt - 1, &last_h)); }
  int64_t m0 = nt ? (int64_t)last_c + last_h : 0;
  *out_ncand = m0;
  S->n_cand = m0;
  // ---- chains for every candidate
  int64_t *ci0_, *cz0, *hills0;
  float *tcos0, *acos0;
  int32_t* nn0;
  CK(dalloc(S, &ci0_, m0)); CK(dalloc(S, &cz0, m0)); CK(dalloc(S, &hills0, m0 * max_iso));
  CK(dalloc(S, &tcos0, m0 * max_iso)); CK(dalloc(S, &acos0, m0)); CK(dalloc(S, &nn0, m0));
  if (m0) k_fill_init<<<nblk(m0 * max_iso), TPB, 0, S->st>>>(m0 * max_iso, hills0, tcos0);
  if (nt) {
    k_cand<<<nblk(nt), TPB, 0, S->st>>>(nh, nz, has1, cidx, charges, ci0_, cz0);
    k_fill<<<nblk(nt), TPB, 0, S->st>>>(I, nh, sidx, nz, has1, cidx, charges, max_iso, hills0, tcos0, acos0, nn0);
  }
  dfree(S, &has1); dfree(S, &cidx); dfree(S, &sidx); dfree(S, &smz); dfree(S, &norm2); dfree(S, &boff);
  dfree(S, &b_mz); dfree(S, &b_h); dfree(S, &b_lo); dfree(S, &b_hi); dfree(S, &bq); dfree(S, &table);
  dfree(S, &charges);
  // ---- nn >= 2, charge ban, compaction (candidate order kept)
  int zw = (int)zmax + 1;
  int32_t *flag, *flag2, *grid;
  int64_t* pos;
  CK(dalloc(S, &flag, m0)); CK(dalloc(S, &flag2, m0 + 1)); CK(dalloc(S, &pos, m0 + 1));
  CK(dalloc(S, &grid, nh * zw));
  CK(cudaMemsetAsync(grid, 0, nh * zw * 4, S->st));
  CK(cudaMemsetAsync(flag2 + m0, 0, 4, S->st));
  if (m0) {
    k_flag_ok<<<nblk(m0), TPB, 0, S->st>>>(m0, nn0, flag);
    k_grid<<<nblk(m0), TPB, 0, S->st>>>(m0, flag, ci0_, cz0, nn0, zw, grid);
    k_ban<<<nblk(m0), TPB, 0, S->st>>>(m0, flag, ci0_, cz0, nn0, zw, grid, use_ban, flag2);
  }
  CK(ex_sum(S, flag2, pos, m0 + 1));
  int64_t m;
  CK(d2h_scalar(S, pos + m0, &m));
  CK(dalloc(S, &S->ci, m)); CK(dalloc(S, &S->cz, m)); CK(dalloc(S, &S->hills, m * max_iso));
  CK(dalloc(S, &S->tcos, m * max_iso)); CK(dalloc(S, &S->acos, m)); CK(dalloc(S, &S->nn, m));
  CK(dalloc(S, &S->order, m));
  if (m0) k_compact_cand<<<nblk(m0), TPB, 0, S->st>>>(m0, flag2, pos, max_iso, ci0_, cz0, hills0, tcos0, acos0, nn0,
                                                      S->ci, S->cz, S->hills, S->tcos, S->acos, S->nn);
  dfree(S, &flag); dfree(S, &flag2); dfree(S, &pos); dfree(S, &grid); dfree(S, &ci0_); dfree(S, &cz0);
  dfree(S, &hills0); dfree(S, &tcos0); dfree(S, &acos0); dfree(S, &nn0);
  // ---- priority order (total order: any correct sort == numpy's lexsort)
  if (m) {
    k_iota64<<<nblk(m), TPB, 0, S->st>>>(m, S->order);
    PrioLess cmp{S->nn, S->acos, S->hapex_int, S->ci, S->cz};
    thrust::sort(thrust::cuda::par.on(S->st), S->order, S->order + m, cmp);
  }
  CK(cudaEventRecord(S->ev[3], S->st));
  CK(ssync(S));
  CK(cudaGetLastError());
  cudaEventElapsedTime(&S->gpu_ms[1], S->ev[2], S->ev[3]);
  S->m = m;
  S->max_iso = max_iso;
  *out_m = m;
  return 0;
}

int ms1g_isotopes_get(void* h, int64_t* ci, int64_t* cz, int64_t* hills, float* tcos, float* acos, int32_t* nn,
                      int64_t* order) {
  Session* S = (Session*)h;
  int64_t m = S->m, mi = S->m * S->max_iso;
#define G(dst, src, cnt) if (dst && (cnt)) CK(cudaMemcpyAsync(dst, src, (cnt) * sizeof(*(dst)), cudaMemcpyDeviceToHost, S->st))
  G(ci, S->ci, m); G(cz, S->cz, m); G(hills, S->hills, mi); G(tcos, S->tcos, mi); G(acos, S->acos, m);
  G(nn, S->nn, m); G(order, S->order, m);
#undef G
  CK(ssync(S));
  return 0;
}

// peak_feature / peak_iso over all n_all flat MS1 peaks (-1 where the peak is in no feature hill)
int ms1g_peak_map(void* h, const int32_t* hill_feat_h, const int8_t* hill_iso_h, int64_t n_all, int32_t* pf_h,
                  int8_t* pi_h) {
  Session* S = (Session*)h;
  int64_t nh = S->n_h;
  CK(cudaEventRecord(S->ev[4], S->st));
  int32_t *hf, *pf;
  int8_t *hi, *pi;
  CK(dalloc(S, &hf, nh)); CK(dalloc(S, &hi, nh)); CK(dalloc(S, &pf, n_all)); CK(dalloc(S, &pi, n_all));
  if (nh) {
    CK(cudaMemcpyAsync(hf, hill_feat_h, nh * 4, cudaMemcpyHostToDevice, S->st));
    CK(cudaMemcpyAsync(hi, hill_iso_h, nh, cudaMemcpyHostToDevice, S->st));
  }
  CK(cudaMemsetAsync(pf, 0xff, n_all * 4, S->st));
  CK(cudaMemsetAsync(pi, 0xff, n_all, S->st));
  if (nh) k_peak_map<<<nblk(nh), TPB, 0, S->st>>>(nh, S->off, S->psrc, hf, hi, pf, pi);
  CK(cudaEventRecord(S->ev[5], S->st));
  CK(ssync(S));
  if (n_all) {
    CK(cudaMemcpyAsync(pf_h, pf, n_all * 4, cudaMemcpyDeviceToHost, S->st));
    CK(cudaMemcpyAsync(pi_h, pi, n_all, cudaMemcpyDeviceToHost, S->st));
  }
  dfree(S, &hf); dfree(S, &hi); dfree(S, &pf); dfree(S, &pi);
  CK(ssync(S));
  CK(cudaGetLastError());
  cudaEventElapsedTime(&S->gpu_ms[2], S->ev[4], S->ev[5]);
  return 0;
}

}  // extern "C"
