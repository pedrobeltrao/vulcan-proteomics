// Vulcan GPU kernels K1 (tiled inverted-index candidate generation) and K2 (exact rescoring),
// plus the fragment emitter used to build the K1 index.  DESIGN.md §4 L1-L3, §10 E2.
//
// Built by scripts/m2/build_kernels.sh (nvcc, sm_80 + sm_86 + PTX, -fmad=false so f32 arithmetic
// rounds exactly like Sage's Rust code).  Loaded with ctypes by vulcan/gpu/_lib.py; every launcher
// takes raw device pointers and a cudaStream_t (torch's current stream).
//
// LABEL BLINDNESS: nothing here receives a target/decoy/entrapment flag.  Candidates are identified by
// their row in the content-sorted peptidoform table (sorted by Sage f32 mass, then hash id).
//
// Fragment arithmetic replicates sage-core 99407db6 (ion_series.rs, scoring.rs, mass.rs) in f32:
//   residue mass = MONO[r] + mod[r] (f32), b = running sum from 0, y = running (peptide mass - sum).
#include <cuda_runtime.h>
#include <stdint.h>
#include <math.h>

#define FULL 0xffffffffu

struct FragParams {
  float mono[26];      // Sage MONOISOTOPIC_MASSES (f32)
  float stat[26];      // static mod per residue (f32, 0 = none)
  float var[26];       // variable mod per residue (f32, 0 = none)
};

__constant__ FragParams c_fp;
// M4: variable peptide-N-term modification (Sage `^`), carried in mask bit 63; 0 when the index has none.
// Sage ion_series.rs: b starts at nterm, y starts at (monoisotopic - nterm).  With the bit unset both are
// exact no-ops (0.f + x, x - 0.f), so closed-search outputs are bit-identical to E2/M2b.
__constant__ float c_nterm;
#define NTERM_BIT (1ull << 63)
__device__ __forceinline__ float nterm_of(uint64_t mask) { return (mask & NTERM_BIT) ? c_nterm : 0.f; }
// 3c: fragment tolerance unit.  0 (default) = ppm: bound = m + m*tol/1e6 (Sage Tolerance::Ppm, the E2 arithmetic,
// unchanged op for op, so OT-HCD outputs are bit-identical); 1 = Da: bound = m + tol (Sage Tolerance::Da), for
// low-resolution ion-trap CID MS2.  Applied to K1 (index lookup), the K1 brute-force reference and K2.
__constant__ int c_tol_da;
__device__ __forceinline__ float tol_bound(float m, float tol) {
  return c_tol_da ? __fadd_rn(m, tol) : __fadd_rn(m, __fdiv_rn(__fmul_rn(m, tol), 1000000.0f));
}

// ------------------------------------------------------------------------------------ index build
// key layout (int64, non-negative): chunk[50..62] | bin[32..49] | local[16..31] | fine[0..15]
// q(m) = floor(log(m / fmin) * inv_logw * 65536); bin = q >> 16, fine = q & 0xffff.
// Matching uses floor() on both the fragment and the tolerance bounds: a quantum (bin/65536 = 3e-10 relative at
// 20 ppm bins) is far finer than the f32 spacing (6e-8 relative), so distinct f32 values never share a quantum and
// q(lo) <= q(f) <= q(hi)  <=>  lo <= f <= hi  exactly as Sage's f32 comparisons (incl. f == lo ties).

__device__ __forceinline__ float resmass(uint8_t r, int k, uint64_t mask) {
  int i = r - 'A';
  float m = ((mask >> k) & 1ull) ? c_fp.var[i] : c_fp.stat[i];
  return __fadd_rn(c_fp.mono[i], m);
}

extern "C" __global__ void frag_count_kernel(int64_t n, const int64_t* seq_off, const uint8_t* len,
                                             const uint64_t* mask, const float* mass32, const uint8_t* res,
                                             int min_ion_index, float fmin, float fmax, int32_t* count) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  const uint8_t* s = res + seq_off[i];
  int L = len[i];
  uint64_t mk = mask[i];
  float b = nterm_of(mk), y = __fsub_rn(mass32[i], nterm_of(mk));
  int c = 0;
  for (int k = 0; k < L - 1; ++k) {
    float rm = resmass(s[k], k, mk);
    b = __fadd_rn(b, rm);
    y = __fsub_rn(y, rm);
    if (k + 1 > min_ion_index && b >= fmin && b <= fmax) ++c;
    if (L - 1 - k > min_ion_index && y >= fmin && y <= fmax) ++c;
  }
  count[i] = c;
}

extern "C" __global__ void frag_emit_kernel(int64_t n, const int64_t* seq_off, const uint8_t* len,
                                            const uint64_t* mask, const float* mass32, const uint8_t* res,
                                            int min_ion_index, float fmin, float fmax, double inv_logw,
                                            int chunk_size, int64_t row0, const int64_t* out_off,
                                            int64_t* keys) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  const uint8_t* s = res + seq_off[i];
  int L = len[i];
  uint64_t mk = mask[i];
  float b = nterm_of(mk), y = __fsub_rn(mass32[i], nterm_of(mk));
  int64_t g = row0 + i;
  int64_t chunk = g / chunk_size, local = g % chunk_size;
  int64_t o = out_off[i];
  double lf = log((double)fmin);
  for (int k = 0; k < L - 1; ++k) {
    float rm = resmass(s[k], k, mk);
    b = __fadd_rn(b, rm);
    y = __fsub_rn(y, rm);
    if (k + 1 > min_ion_index && b >= fmin && b <= fmax) {
      int64_t q = (int64_t)floor((log((double)b) - lf) * inv_logw * 65536.0);
      keys[o++] = (chunk << 50) | ((q >> 16) << 32) | (local << 16) | (q & 0xffff);
    }
    if (L - 1 - k > min_ion_index && y >= fmin && y <= fmax) {
      int64_t q = (int64_t)floor((log((double)y) - lf) * inv_logw * 65536.0);
      keys[o++] = (chunk << 50) | ((q >> 16) << 32) | (local << 16) | (q & 0xffff);
    }
  }
}

// ------------------------------------------------------------------------------------ K1
// One warp per work item (hypothesis, chunk, [lo, hi) local-id tile).  Shared-memory u32 count per
// candidate in the tile.  Lanes stride over the spectrum's peaks (x fragment charges); for each peak the
// q-range of Sage's tolerance bounds is looked up in the (chunk, bin) posting lists, each list binary-
// searched for local >= lo.  Then lanes scan the tile and append every candidate with count >= thr
// (identical threshold for every candidate) to the global buffer.

#define K1_WARPS 4

extern "C" __global__ void k1_kernel(int n_items, const int32_t* it_hyp, const int32_t* it_chunk,
                                     const int32_t* it_lo, const int32_t* it_hi,
                                     const int32_t* hyp_spec, const int8_t* hyp_maxfz,
                                     const int64_t* pk_off, const float* pk_mass,
                                     const uint32_t* offsets, const uint32_t* entries,
                                     int n_bins, int chunk0, int chunk_size, int tile,
                                     float tol_lo, float tol_hi, float fmin, double inv_logw,
                                     int thr, int64_t cap, unsigned long long* n_out,
                                     int32_t* out_hyp, int32_t* out_row, int32_t* out_cnt,
                                     unsigned int* hyp_scored, unsigned long long* hyp_matched) {
  extern __shared__ uint32_t smem[];
  int warp = threadIdx.x >> 5, lane = threadIdx.x & 31;
  int item = blockIdx.x * K1_WARPS + warp;
  if (item >= n_items) return;
  uint32_t* cnt = smem + warp * tile;
  int h = it_hyp[item], c = it_chunk[item], lo = it_lo[item], hi = it_hi[item];
  int n = hi - lo;
  for (int j = lane; j < n; j += 32) cnt[j] = 0;
  __syncwarp();
  int s = hyp_spec[h];
  int maxfz = hyp_maxfz[h];
  int64_t p0 = pk_off[s], p1 = pk_off[s + 1];
  const uint32_t* off = offsets + (int64_t)(c - chunk0) * n_bins;
  double lf = log((double)fmin);
  uint32_t lo_key = (uint32_t)lo << 16;
  uint32_t hi_key = (uint32_t)hi << 16;
  for (int64_t p = p0 + lane; p < p1; p += 32) {
    float pm = pk_mass[p];
    for (int z = 1; z < maxfz; ++z) {
      float m = __fmul_rn(pm, (float)z);
      float a = tol_bound(m, tol_lo);
      float bnd = tol_bound(m, tol_hi);
      if (bnd < fmin) continue;
      int64_t qlo = (int64_t)floor((log((double)a) - lf) * inv_logw * 65536.0);
      int64_t qhi = (int64_t)floor((log((double)bnd) - lf) * inv_logw * 65536.0);
      if (qlo < 0) qlo = 0;
      int64_t blo = qlo >> 16, bhi = qhi >> 16;
      if (bhi >= n_bins) bhi = n_bins - 1;
      for (int64_t bb = blo; bb <= bhi; ++bb) {
        uint32_t e0 = off[bb], e1 = off[bb + 1];
        // lower bound of lo_key in entries[e0, e1)
        uint32_t L = e0, R = e1;
        while (L < R) {
          uint32_t M = (L + R) >> 1;
          if (entries[M] < lo_key) L = M + 1; else R = M;
        }
        for (uint32_t e = L; e < e1; ++e) {
          uint32_t v = entries[e];
          if (v >= hi_key) break;
          int64_t q = (bb << 16) | (v & 0xffff);
          if (q >= qlo && q <= qhi) atomicAdd(&cnt[(v >> 16) - lo], 1u);
        }
      }
    }
  }
  __syncwarp();
  unsigned int scored = 0;
  unsigned long long matched = 0;
  for (int j = lane; j < n; j += 32) {
    uint32_t k = cnt[j];
    if (k > 0) { ++scored; matched += k; }
    if ((int)k >= thr) {
      unsigned long long w = atomicAdd(n_out, 1ull);
      if ((int64_t)w < cap) {
        out_hyp[w] = h;
        out_row[w] = c * chunk_size + lo + j;
        out_cnt[w] = (int32_t)k;
      }
    }
  }
  for (int o = 16; o > 0; o >>= 1) {
    scored += __shfl_down_sync(FULL, scored, o);
    matched += __shfl_down_sync(FULL, matched, o);
  }
  if (lane == 0) {
    if (scored) atomicAdd(&hyp_scored[h], scored);
    if (matched) atomicAdd(&hyp_matched[h], matched);
  }
}

// Brute-force preliminary counts (reference for the index): one thread per (hypothesis, candidate row)
// pair, fragments generated on the fly, every peak checked with Sage's page_search f32 tolerance semantics
// (theoretical fragment in [m(1+lo), m(1+hi)] of the peak mass m), index ion filter applied.
extern "C" __global__ void k1_brute_kernel(int64_t n, const int32_t* c_hyp, const int32_t* c_row,
                                           const int32_t* hyp_spec, const int8_t* hyp_maxfz,
                                           const int64_t* pk_off, const float* pk_mass,
                                           const int64_t* seq_off, const uint8_t* len, const uint64_t* mask,
                                           const float* mass32, const uint8_t* res, int min_ion_index,
                                           float fmin, float fmax, float tol_lo, float tol_hi,
                                           double inv_logw, int32_t* out_cnt) {
  int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
  if (i >= n) return;
  int h = c_hyp[i], r = c_row[i];
  int s = hyp_spec[h], maxfz = hyp_maxfz[h];
  int64_t p0 = pk_off[s], p1 = pk_off[s + 1];
  const uint8_t* sq = res + seq_off[r];
  int L = len[r];
  uint64_t mk = mask[r];
  float b = nterm_of(mk), y = __fsub_rn(mass32[r], nterm_of(mk));
  int cnt = 0;
  for (int k = 0; k < L - 1; ++k) {
    float rm = resmass(sq[k], k, mk);
    b = __fadd_rn(b, rm);
    y = __fsub_rn(y, rm);
    for (int t = 0; t < 2; ++t) {
      float f = t == 0 ? b : y;
      bool ok = t == 0 ? (k + 1 > min_ion_index) : (L - 1 - k > min_ion_index);
      if (!ok || f < fmin || f > fmax) continue;
      // plain f32 comparisons (Sage page_search semantics): independent of the index's log-ppm quantization
      for (int64_t p = p0; p < p1; ++p) {
        for (int z = 1; z < maxfz; ++z) {
          float m = __fmul_rn(pk_mass[p], (float)z);
          float a = tol_bound(m, tol_lo);
          float bd = tol_bound(m, tol_hi);
          if (f >= a && f <= bd) ++cnt;
        }
      }
    }
  }
  out_cnt[i] = cnt;
}

// ------------------------------------------------------------------------------------ K2
// Exact rescoring, Sage scoring.rs::score_candidate semantics.  One block per spectrum: its processed peaks
// (mass-sorted) go to shared memory; threads stride over the spectrum's candidates; each thread generates all
// b and y ions (no index filter), for fragment charges 1..maxfz-1, and picks the most intense peak within
// tolerance (ties -> the later peak, as Sage's `>=`).

struct K2Out {
  double hyperscore;
  float summed_b, summed_y, ppm;
  int32_t matched_b, matched_y, longest_b, longest_y;
};

__device__ __forceinline__ double lnfact(int n) {
  if (n == 0) return 1.0;  // sic (sage scoring.rs)
  double x = (double)n;
  return x * log(x) - x + 0.5 * log(x) + 0.5 * log(2.0 * 3.141592653589793 * x);
}

struct Run { int start, length, last, longest; };
__device__ __forceinline__ void run_matched(Run& r, int idx) {
  if (r.last == idx) return;
  if (r.start + r.length == idx) { r.length += 1; r.longest = max(r.longest, r.length); }
  else { r.start = idx; r.length = 1; r.longest = max(r.longest, r.length); }
  r.last = idx;
}

#define K2_MAXPK 512
#define K2_MAXLEN 64   // peptide length bound for the localization-aware shifted scoring (index max_len 50)

extern "C" __global__ void k2_kernel(int n_spec, const int64_t* cand_off, const int32_t* c_row,
                                     const int8_t* c_maxfz, const int32_t* spec_ids,
                                     const int64_t* pk_off, const float* pk_mass, const float* pk_int,
                                     const int64_t* seq_off, const uint8_t* len, const uint64_t* mask,
                                     const float* mass32, const uint8_t* res, float tol_lo, float tol_hi,
                                     const float* c_shift,
                                     double* o_hs, float* o_sb, float* o_sy, float* o_ppm,
                                     int16_t* o_mb, int16_t* o_my, int16_t* o_lb, int16_t* o_ly) {
  __shared__ float sm[K2_MAXPK];
  __shared__ float si[K2_MAXPK];
  int bs = blockIdx.x;
  if (bs >= n_spec) return;
  int s = spec_ids[bs];
  int64_t p0 = pk_off[s];
  int np = (int)(pk_off[s + 1] - p0);
  if (np > K2_MAXPK) np = K2_MAXPK;
  for (int j = threadIdx.x; j < np; j += blockDim.x) { sm[j] = pk_mass[p0 + j]; si[j] = pk_int[p0 + j]; }
  __syncthreads();
  int64_t c0 = cand_off[bs], c1 = cand_off[bs + 1];
  for (int64_t c = c0 + threadIdx.x; c < c1; c += blockDim.x) {
    int r = c_row[c];
    int maxfz = c_maxfz[c];
    const uint8_t* sq = res + seq_off[r];
    int L = len[r];
    uint64_t mk = mask[r];
    float sb = 0.f, sy = 0.f, ppm = 0.f;
    int mb = 0, my = 0;
    Run rb = {0, 0, 0, 0}, ry = {0, 0, 0, 0};
    // M4 open / mass-offset mode: c_shift[c] = observed - candidate mass (unlocalized delta).  Every b and y ion
    // is looked up unshifted AND shifted by delta; the more intense peak of the two windows counts once
    // (shifted + unshifted ions, MSFragger/pFind style).  c_shift == nullptr: closed search, unchanged code path.
    float shift = c_shift ? c_shift[c] : 0.f;
    int n_alt = (c_shift && shift != 0.f) ? 2 : 1;
    // Localization-aware shifted scoring (n_alt == 2): per ion index k keep the unshifted and the shifted match
    // separately (summed over fragment charges); the delta is then placed on the single residue s that maximises
    // the hyperscore (b_k carries it iff k >= s, y_k iff s > k), as MSFragger's localization-aware open search.
    // ppm and the longest ladders are taken over the union of both windows (features only).
    float iu[2][K2_MAXLEN], is_[2][K2_MAXLEN];
    unsigned char nu[2][K2_MAXLEN], ns[2][K2_MAXLEN];
    if (n_alt == 2) {
      for (int k = 0; k < K2_MAXLEN; ++k) { iu[0][k] = iu[1][k] = is_[0][k] = is_[1][k] = 0.f;
                                            nu[0][k] = nu[1][k] = ns[0][k] = ns[1][k] = 0; }
    }
    for (int kind = 0; kind < 2; ++kind) {
      float cum = kind == 0 ? nterm_of(mk) : __fsub_rn(mass32[r], nterm_of(mk));
      for (int k = 0; k < L - 1; ++k) {
        float rm = resmass(sq[k], k, mk);
        cum = kind == 0 ? __fadd_rn(cum, rm) : __fsub_rn(cum, rm);
        for (int z = 1; z < maxfz; ++z) {
          int best = -1;
          float mx = 0.f, mz = 0.f, best_mz = 0.f;
          for (int alt = 0; alt < n_alt; ++alt) {
            int best_a = -1;
            float mx_a = 0.f;
            float base = alt == 0 ? cum : __fadd_rn(cum, shift);
            mz = __fdiv_rn(base, (float)z);
            float a = tol_bound(mz, tol_lo);
            float bd = tol_bound(mz, tol_hi);
            int L2 = 0, R2 = np;
            while (L2 < R2) { int M = (L2 + R2) >> 1; if (sm[M] < a) L2 = M + 1; else R2 = M; }
            for (int j = L2; j < np && sm[j] <= bd; ++j) {
              if (si[j] >= mx) { mx = si[j]; best = j; best_mz = mz; }
              if (si[j] >= mx_a) { mx_a = si[j]; best_a = j; }
            }
            if (n_alt == 2 && best_a >= 0 && k < K2_MAXLEN) {
              if (alt == 0) { iu[kind][k] = __fadd_rn(iu[kind][k], si[best_a]); nu[kind][k] += 1; }
              else { is_[kind][k] = __fadd_rn(is_[kind][k], si[best_a]); ns[kind][k] += 1; }
            }
          }
          mz = best_mz;
          if (best >= 0) {
            float pmass = sm[best], pint = si[best];
            ppm = __fadd_rn(ppm, __fdiv_rn(__fmul_rn(__fmul_rn(pint, fabsf(__fsub_rn(mz, pmass))), 2.0e6f),
                                           __fadd_rn(mz, pmass)));
            if (kind == 0) { mb += 1; sb = __fadd_rn(sb, pint); run_matched(rb, k); }
            else { my += 1; sy = __fadd_rn(sy, pint); run_matched(ry, k); }
          }
        }
      }
    }
    if (n_alt == 2) {
      // site s in [0, L): b_k (k = 0..L-2, residues 0..k) shifted iff k >= s; y_k (residues k+1..L-1) iff k < s
      int Lk = min(L - 1, K2_MAXLEN);
      float b_u = 0.f, y_s = 0.f, b_s = 0.f, y_u = 0.f;     // running sums for s = 0: all b shifted, all y unshifted
      int nb_u = 0, ny_s = 0, nb_s = 0, ny_u = 0;
      for (int k = 0; k < Lk; ++k) { b_s = __fadd_rn(b_s, is_[0][k]); nb_s += ns[0][k];
                                     y_u = __fadd_rn(y_u, iu[1][k]); ny_u += nu[1][k]; }
      double best_hs = -1e300;
      for (int st = 0; st < L; ++st) {
        if (st > 0 && st - 1 < Lk) {        // ion index st-1 moves: b unshifted, y shifted
          int k = st - 1;
          b_u = __fadd_rn(b_u, iu[0][k]); nb_u += nu[0][k]; b_s = __fsub_rn(b_s, is_[0][k]); nb_s -= ns[0][k];
          y_s = __fadd_rn(y_s, is_[1][k]); ny_s += ns[1][k]; y_u = __fsub_rn(y_u, iu[1][k]); ny_u -= nu[1][k];
        }
        float sb2 = __fadd_rn(b_u, b_s), sy2 = __fadd_rn(y_u, y_s);
        int mb2 = nb_u + nb_s, my2 = ny_u + ny_s;
        double h2 = log((double)__fadd_rn(fmaxf(sb2, 0.f), 1.0f) * (double)__fadd_rn(fmaxf(sy2, 0.f), 1.0f))
                    + lnfact(mb2) + lnfact(my2);
        if (h2 > best_hs) { best_hs = h2; mb = mb2; my = my2; sb = fmaxf(sb2, 0.f); sy = fmaxf(sy2, 0.f); }
      }
    }
    double i2 = (double)__fadd_rn(sb, 1.0f) * (double)__fadd_rn(sy, 1.0f);
    double hs = log(i2) + lnfact(mb) + lnfact(my);
    if (!isfinite(hs)) hs = 255.0;
    o_hs[c] = hs;
    o_sb[c] = sb;
    o_sy[c] = sy;
    o_ppm[c] = __fdiv_rn(ppm, __fadd_rn(sb, sy));
    o_mb[c] = (int16_t)mb;
    o_my[c] = (int16_t)my;
    o_lb[c] = (int16_t)rb.longest;
    o_ly[c] = (int16_t)ry.longest;
  }
}

// ------------------------------------------------------------------------------------ launchers
extern "C" {

int vk_set_params(const float* mono, const float* stat, const float* var) {
  FragParams fp;
  for (int i = 0; i < 26; ++i) { fp.mono[i] = mono[i]; fp.stat[i] = stat[i]; fp.var[i] = var[i]; }
  return (int)cudaMemcpyToSymbol(c_fp, &fp, sizeof(FragParams));
}

int vk_set_nterm(float nterm) { return (int)cudaMemcpyToSymbol(c_nterm, &nterm, sizeof(float)); }

int vk_set_tol_da(int da) { return (int)cudaMemcpyToSymbol(c_tol_da, &da, sizeof(int)); }

int vk_frag_count(int64_t n, const int64_t* seq_off, const uint8_t* len, const uint64_t* mask,
                  const float* mass32, const uint8_t* res, int min_ion_index, float fmin, float fmax,
                  int32_t* count, cudaStream_t st) {
  if (n == 0) return 0;
  int bs = 256;
  frag_count_kernel<<<(unsigned)((n + bs - 1) / bs), bs, 0, st>>>(n, seq_off, len, mask, mass32, res,
                                                                    min_ion_index, fmin, fmax, count);
  return (int)cudaGetLastError();
}

int vk_frag_emit(int64_t n, const int64_t* seq_off, const uint8_t* len, const uint64_t* mask,
                 const float* mass32, const uint8_t* res, int min_ion_index, float fmin, float fmax,
                 double inv_logw, int chunk_size, int64_t row0, const int64_t* out_off, int64_t* keys,
                 cudaStream_t st) {
  if (n == 0) return 0;
  int bs = 256;
  frag_emit_kernel<<<(unsigned)((n + bs - 1) / bs), bs, 0, st>>>(n, seq_off, len, mask, mass32, res,
                                                                   min_ion_index, fmin, fmax, inv_logw,
                                                                   chunk_size, row0, out_off, keys);
  return (int)cudaGetLastError();
}

int vk_k1(int n_items, const int32_t* it_hyp, const int32_t* it_chunk, const int32_t* it_lo,
          const int32_t* it_hi, const int32_t* hyp_spec, const int8_t* hyp_maxfz, const int64_t* pk_off,
          const float* pk_mass, const uint32_t* offsets, const uint32_t* entries, int n_bins, int chunk0,
          int chunk_size, int tile, float tol_lo, float tol_hi, float fmin, double inv_logw, int thr,
          int64_t cap, unsigned long long* n_out, int32_t* out_hyp, int32_t* out_row, int32_t* out_cnt,
          unsigned int* hyp_scored, unsigned long long* hyp_matched, cudaStream_t st) {
  if (n_items == 0) return 0;
  size_t sh = (size_t)K1_WARPS * tile * sizeof(uint32_t);
  if (sh > 48 * 1024) {
    cudaError_t e = cudaFuncSetAttribute(k1_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)sh);
    if (e != cudaSuccess) return (int)e;
  }
  unsigned grid = (unsigned)((n_items + K1_WARPS - 1) / K1_WARPS);
  k1_kernel<<<grid, K1_WARPS * 32, sh, st>>>(n_items, it_hyp, it_chunk, it_lo, it_hi, hyp_spec, hyp_maxfz,
                                              pk_off, pk_mass, offsets, entries, n_bins, chunk0, chunk_size,
                                              tile, tol_lo, tol_hi, fmin, inv_logw, thr, cap, n_out, out_hyp,
                                              out_row, out_cnt, hyp_scored, hyp_matched);
  return (int)cudaGetLastError();
}

int vk_k1_brute(int64_t n, const int32_t* c_hyp, const int32_t* c_row, const int32_t* hyp_spec,
                const int8_t* hyp_maxfz, const int64_t* pk_off, const float* pk_mass, const int64_t* seq_off,
                const uint8_t* len, const uint64_t* mask, const float* mass32, const uint8_t* res,
                int min_ion_index, float fmin, float fmax, float tol_lo, float tol_hi, double inv_logw,
                int32_t* out_cnt, cudaStream_t st) {
  if (n == 0) return 0;
  int bs = 128;
  k1_brute_kernel<<<(unsigned)((n + bs - 1) / bs), bs, 0, st>>>(n, c_hyp, c_row, hyp_spec, hyp_maxfz, pk_off,
                                                                  pk_mass, seq_off, len, mask, mass32, res,
                                                                  min_ion_index, fmin, fmax, tol_lo, tol_hi,
                                                                  inv_logw, out_cnt);
  return (int)cudaGetLastError();
}

int vk_k2(int n_spec, const int64_t* cand_off, const int32_t* c_row, const int8_t* c_maxfz,
          const int32_t* spec_ids, const int64_t* pk_off, const float* pk_mass, const float* pk_int,
          const int64_t* seq_off, const uint8_t* len, const uint64_t* mask, const float* mass32,
          const uint8_t* res, float tol_lo, float tol_hi, double* o_hs, float* o_sb, float* o_sy,
          float* o_ppm, int16_t* o_mb, int16_t* o_my, int16_t* o_lb, int16_t* o_ly, int threads,
          cudaStream_t st) {
  if (n_spec == 0) return 0;
  k2_kernel<<<n_spec, threads, 0, st>>>(n_spec, cand_off, c_row, c_maxfz, spec_ids, pk_off, pk_mass, pk_int,
                                        seq_off, len, mask, mass32, res, tol_lo, tol_hi, nullptr, o_hs, o_sb, o_sy,
                                        o_ppm, o_mb, o_my, o_lb, o_ly);
  return (int)cudaGetLastError();
}

// M4: K2 with a per-candidate precursor delta (open / mass-offset search; shifted + unshifted ions).
int vk_k2_shift(int n_spec, const int64_t* cand_off, const int32_t* c_row, const int8_t* c_maxfz,
                const int32_t* spec_ids, const int64_t* pk_off, const float* pk_mass, const float* pk_int,
                const int64_t* seq_off, const uint8_t* len, const uint64_t* mask, const float* mass32,
                const uint8_t* res, float tol_lo, float tol_hi, const float* c_shift, double* o_hs, float* o_sb,
                float* o_sy, float* o_ppm, int16_t* o_mb, int16_t* o_my, int16_t* o_lb, int16_t* o_ly,
                int threads, cudaStream_t st) {
  if (n_spec == 0) return 0;
  k2_kernel<<<n_spec, threads, 0, st>>>(n_spec, cand_off, c_row, c_maxfz, spec_ids, pk_off, pk_mass, pk_int,
                                        seq_off, len, mask, mass32, res, tol_lo, tol_hi, c_shift, o_hs, o_sb, o_sy,
                                        o_ppm, o_mb, o_my, o_lb, o_ly);
  return (int)cudaGetLastError();
}

}  // extern "C"
