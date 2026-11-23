// Micro-benchmark: scalar vs AVX2 uppass state-update loop.
//
// Isolates the hot loop in ts::uppass_node (ts_fitch.cpp:54-61), which VTune
// shows at 13.2% of DLL self-CPU on the standard-Fitch (Zhu2013, "-"->"?")
// workload.  The any_intersect reduce (line 47) is ALREADY vectorised
// (any_hit_reduce_avx2); the *state update* that follows is scalar, unlike
// the analogous fitch_downpass which uses simd::fitch_combine (line 102).
//
// uppass new_val[s] = (anc_final[s] & node_prelim[s] & has_isect)
//                   | (node_prelim[s] & no_isect)
// plus a "changed" flag (new_val != old final_) that drives the incremental
// dirty-propagation.  This bench reproduces BOTH the value and the changed
// flag, asserts the AVX2 path is bit-identical to scalar, then times them.
//
// Build (matches R's flags + AVX2):
//   g++ -O2 -mavx2 -std=c++17 -o bench_uppass.exe bench_uppass_combine.cpp
// Run: ./bench_uppass.exe

#include <cstdint>
#include <cstdio>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <immintrin.h>

static inline uint64_t popcount64(uint64_t x){ return __builtin_popcountll(x); }

// ---- horizontal OR of a 256-bit reg (matches ts_simd.h hor_or256) ----
static inline uint64_t hor_or256(__m256i v){
  __m128i lo = _mm256_castsi256_si128(v);
  __m128i hi = _mm256_extracti128_si256(v, 1);
  __m128i c  = _mm_or_si128(lo, hi);
  __m128i s  = _mm_srli_si128(c, 8);
  c = _mm_or_si128(c, s);
  return (uint64_t)_mm_cvtsi128_si64(c);
}
static inline uint64_t any_hit_reduce4(const uint64_t* a, const uint64_t* b){
  __m256i va = _mm256_loadu_si256((const __m256i*)a);
  __m256i vb = _mm256_loadu_si256((const __m256i*)b);
  return hor_or256(_mm256_and_si256(va, vb));
}

// ---- BASELINE: scalar update loop (verbatim from uppass_node) ----
static inline bool update_scalar(const uint64_t* anc_final,
                                 const uint64_t* node_prelim,
                                 uint64_t* node_final, int n_states,
                                 uint64_t has_isect, uint64_t no_isect){
  bool changed = false;
  for (int s = 0; s < n_states; ++s){
    uint64_t isect = anc_final[s] & node_prelim[s];
    uint64_t new_val = (isect & has_isect) | (node_prelim[s] & no_isect);
    if (new_val != node_final[s]) changed = true;
    node_final[s] = new_val;
  }
  return changed;
}

// ---- CANDIDATE: AVX2 update for n_states==4 (+ scalar tail) ----
static inline bool update_avx2(const uint64_t* anc_final,
                               const uint64_t* node_prelim,
                               uint64_t* node_final, int n_states,
                               uint64_t has_isect, uint64_t no_isect){
  int s = 0;
  uint64_t diff_acc = 0;
  if (n_states >= 4){
    __m256i H = _mm256_set1_epi64x((long long)has_isect);
    __m256i N = _mm256_set1_epi64x((long long)no_isect);
    __m256i diff = _mm256_setzero_si256();
    for (; s + 4 <= n_states; s += 4){
      __m256i a = _mm256_loadu_si256((const __m256i*)(anc_final + s));
      __m256i p = _mm256_loadu_si256((const __m256i*)(node_prelim + s));
      __m256i nv = _mm256_or_si256(
          _mm256_and_si256(_mm256_and_si256(a, p), H),
          _mm256_and_si256(p, N));
      __m256i old = _mm256_loadu_si256((const __m256i*)(node_final + s));
      diff = _mm256_or_si256(diff, _mm256_xor_si256(nv, old));
      _mm256_storeu_si256((__m256i*)(node_final + s), nv);
    }
    diff_acc |= hor_or256(diff);
  }
  for (; s < n_states; ++s){
    uint64_t isect = anc_final[s] & node_prelim[s];
    uint64_t new_val = (isect & has_isect) | (node_prelim[s] & no_isect);
    diff_acc |= (new_val ^ node_final[s]);
    node_final[s] = new_val;
  }
  return diff_acc != 0;
}

int main(){
  const int N = 200000;        // node-blocks
  const int NS = 4;            // states (DNA / common morph)
  const int REPS = 60;
  std::mt19937_64 rng(5813);
  // 4-state data: each word has ~half its bits set among 64 chars.
  std::vector<uint64_t> anc(N*NS), prelim(N*NS), final_a(N*NS), final_b(N*NS);
  std::vector<uint64_t> hasv(N), nov(N);
  for (int i=0;i<N*NS;++i){ anc[i]=rng(); prelim[i]=rng(); }
  for (int i=0;i<N;++i){
    uint64_t ai = any_hit_reduce4(&anc[i*NS], &prelim[i*NS]);
    hasv[i] = ai; nov[i] = ~ai;               // active_mask = all-ones here
    for (int s=0;s<NS;++s){ final_a[i*NS+s]=rng(); final_b[i*NS+s]=final_a[i*NS+s]; }
  }

  // Correctness: both kernels must agree on value AND changed flag.
  int mism=0, chg_mism=0;
  for (int i=0;i<N;++i){
    uint64_t ta[NS], tb[NS];
    for (int s=0;s<NS;++s){ ta[s]=final_a[i*NS+s]; tb[s]=final_a[i*NS+s]; }
    bool ca = update_scalar(&anc[i*NS], &prelim[i*NS], ta, NS, hasv[i], nov[i]);
    bool cb = update_avx2  (&anc[i*NS], &prelim[i*NS], tb, NS, hasv[i], nov[i]);
    if (ca!=cb) ++chg_mism;
    for (int s=0;s<NS;++s) if (ta[s]!=tb[s]) { ++mism; break; }
  }
  printf("Correctness: value mismatches=%d  changed-flag mismatches=%d  (want 0,0)\n",
         mism, chg_mism);

  auto run = [&](int which)->double{
    std::vector<double> ms;
    for (int r=0;r<REPS;++r){
      // reset finals each rep so the work (and store traffic) is identical
      std::copy(final_a.begin(), final_a.end(), final_b.begin());
      volatile uint64_t sink = 0;
      auto t0 = std::chrono::steady_clock::now();
      for (int i=0;i<N;++i){
        bool c = which==0
          ? update_scalar(&anc[i*NS], &prelim[i*NS], &final_b[i*NS], NS, hasv[i], nov[i])
          : update_avx2  (&anc[i*NS], &prelim[i*NS], &final_b[i*NS], NS, hasv[i], nov[i]);
        sink += c;
      }
      auto t1 = std::chrono::steady_clock::now();
      (void)sink;
      ms.push_back(std::chrono::duration<double,std::milli>(t1-t0).count());
    }
    std::sort(ms.begin(), ms.end());
    return ms[ms.size()/2];
  };

  double s_med = run(0);
  double a_med = run(1);
  printf("scalar : %.3f ms/pass  (%.2f ns/node-block)\n", s_med, s_med*1e6/N);
  printf("avx2   : %.3f ms/pass  (%.2f ns/node-block)\n", a_med, a_med*1e6/N);
  printf("speedup: %.2fx\n", s_med/a_med);
  return 0;
}
