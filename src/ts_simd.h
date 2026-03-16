#ifndef TS_SIMD_H
#define TS_SIMD_H

// SIMD portability layer for TreeSearch bit-parallel scoring.
//
// Provides a thin abstraction over SSE2 (x86_64) and NEON (arm64)
// intrinsics. All functions are inline and header-only. A scalar
// fallback is always available via #ifdef guards.

#include <cstdint>

// ---------- Detect available SIMD level at compile time ----------

#if defined(__AVX2__)
  #define TS_SIMD_AVX2
  #define TS_SIMD_SSE2   // AVX2 implies SSE2
  #include <immintrin.h>
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
  #define TS_SIMD_SSE2
  #include <emmintrin.h>
#elif defined(__ARM_NEON) || defined(__aarch64__)
  #define TS_SIMD_NEON
  #include <arm_neon.h>
#endif

namespace ts {
namespace simd {

// ---------- 128-bit operations (SSE2 / NEON) ----------

#if defined(TS_SIMD_SSE2)

using v128 = __m128i;

inline v128 loadu128(const uint64_t* p) {
  return _mm_loadu_si128(reinterpret_cast<const v128*>(p));
}
inline void storeu128(uint64_t* p, v128 v) {
  _mm_storeu_si128(reinterpret_cast<v128*>(p), v);
}
inline v128 and128(v128 a, v128 b) { return _mm_and_si128(a, b); }
inline v128 or128(v128 a, v128 b)  { return _mm_or_si128(a, b); }
inline v128 andnot128(v128 a, v128 b) { return _mm_andnot_si128(a, b); }
inline v128 zero128() { return _mm_setzero_si128(); }

inline v128 set1_64(uint64_t x) { return _mm_set1_epi64x(static_cast<long long>(x)); }

inline uint64_t hor_or128(v128 v) {
  uint64_t tmp[2];
  storeu128(tmp, v);
  return tmp[0] | tmp[1];
}

#elif defined(TS_SIMD_NEON)

using v128 = uint64x2_t;

inline v128 loadu128(const uint64_t* p) { return vld1q_u64(p); }
inline void storeu128(uint64_t* p, v128 v) { vst1q_u64(p, v); }
inline v128 and128(v128 a, v128 b) { return vandq_u64(a, b); }
inline v128 or128(v128 a, v128 b)  { return vorrq_u64(a, b); }
inline v128 andnot128(v128 a, v128 b) {
  // SSE2 andnot: ~a & b.  NEON equivalent: bic(b, a) = b & ~a
  return vbicq_u64(b, a);
}
inline v128 zero128() { return vdupq_n_u64(0); }

inline v128 set1_64(uint64_t x) { return vdupq_n_u64(x); }

inline uint64_t hor_or128(v128 v) {
  return vgetq_lane_u64(v, 0) | vgetq_lane_u64(v, 1);
}

#endif

// ---------- Convenience: SIMD any_hit reduction ----------
//
// Computes OR( a[s] & b[s] ) for s in [0, n_states).
// Uses SSE2/NEON when available, scalar fallback otherwise.

inline uint64_t any_hit_reduce(const uint64_t* a, const uint64_t* b,
                                int n_states) {
#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
  v128 acc = zero128();
  int s = 0;
  for (; s + 2 <= n_states; s += 2) {
    v128 va = loadu128(&a[s]);
    v128 vb = loadu128(&b[s]);
    acc = or128(acc, and128(va, vb));
  }
  uint64_t result = hor_or128(acc);
  for (; s < n_states; ++s) {
    result |= (a[s] & b[s]);
  }
  return result;
#else
  uint64_t result = 0;
  for (int s = 0; s < n_states; ++s) {
    result |= (a[s] & b[s]);
  }
  return result;
#endif
}

// Same but skips word 0 (for inapplicable blocks where state 0 = NA).
inline uint64_t any_hit_reduce_from1(const uint64_t* a, const uint64_t* b,
                                      int n_states) {
  if (n_states <= 1) return 0;
  return any_hit_reduce(a + 1, b + 1, n_states - 1);
}

// ---------- Three-operand variant: OR( clip[s] & (a[s] | b[s]) ) ----------
//
// For non-cached indirect length: vroot = final_[A] | final_[D] computed inline.

inline uint64_t any_hit_reduce3(const uint64_t* clip, const uint64_t* a,
                                 const uint64_t* b, int n_states) {
#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
  v128 acc = zero128();
  int s = 0;
  for (; s + 2 <= n_states; s += 2) {
    v128 vc = loadu128(&clip[s]);
    v128 va = loadu128(&a[s]);
    v128 vb = loadu128(&b[s]);
    acc = or128(acc, and128(vc, or128(va, vb)));
  }
  uint64_t result = hor_or128(acc);
  for (; s < n_states; ++s) {
    result |= (clip[s] & (a[s] | b[s]));
  }
  return result;
#else
  uint64_t result = 0;
  for (int s = 0; s < n_states; ++s) {
    result |= (clip[s] & (a[s] | b[s]));
  }
  return result;
#endif
}

// Three-operand variant starting from word 1 (NA skip).
inline uint64_t any_hit_reduce3_from1(const uint64_t* clip, const uint64_t* a,
                                       const uint64_t* b, int n_states) {
  if (n_states <= 1) return 0;
  return any_hit_reduce3(clip + 1, a + 1, b + 1, n_states - 1);
}

// ---------- Convenience: single-array OR reduction ----------
//
// Computes OR( a[s] ) for s in [start, n_states).

inline uint64_t or_reduce(const uint64_t* a, int n_states, int start = 0) {
#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
  v128 acc = zero128();
  int s = start;
  for (; s + 2 <= n_states; s += 2) {
    acc = or128(acc, loadu128(&a[s]));
  }
  uint64_t result = hor_or128(acc);
  for (; s < n_states; ++s) {
    result |= a[s];
  }
  return result;
#else
  uint64_t result = 0;
  for (int s = start; s < n_states; ++s) {
    result |= a[s];
  }
  return result;
#endif
}

} // namespace simd
} // namespace ts

#endif // TS_SIMD_H
