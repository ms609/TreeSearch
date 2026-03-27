#ifndef TS_SIMD_H
#define TS_SIMD_H

// SIMD portability layer for TreeSearch bit-parallel scoring.
//
// Provides a thin abstraction over SSE2 (x86_64) and NEON (arm64)
// intrinsics, with AVX2 (256-bit) runtime dispatch on x86_64.
// All functions are inline and header-only.

#include <cstdint>

// ---------- Detect available SIMD level ----------

#if defined(__x86_64__) || defined(_M_X64) || defined(_M_AMD64)
  #define TS_SIMD_X86_64
  #define TS_SIMD_SSE2
  // Include immintrin.h unconditionally on x86_64: the target("avx2")
  // attribute allows AVX2 intrinsics in annotated functions even when
  // the translation unit isn't compiled with -mavx2.
  #include <immintrin.h>
#elif defined(__ARM_NEON) || defined(__aarch64__)
  #define TS_SIMD_NEON
  #include <arm_neon.h>
#endif

// Compile-time AVX2 (set when -mavx2 is on the command line)
#if defined(__AVX2__)
  #define TS_SIMD_AVX2
#endif

// Runtime target attribute (GCC/Clang, including Rtools GCC on Windows)
#if defined(TS_SIMD_X86_64) && (defined(__GNUC__) || defined(__clang__))
  #define TS_TARGET_AVX2 __attribute__((target("avx2")))
  #define TS_HAS_AVX2_DISPATCH 1
#else
  #define TS_TARGET_AVX2
  #define TS_HAS_AVX2_DISPATCH 0
#endif

namespace ts {
namespace simd {

// ---------- Runtime AVX2 detection ----------
//
// Cached flag: evaluated once on first call.  On non-x86, always false.

inline bool cpu_has_avx2() {
#if TS_HAS_AVX2_DISPATCH
  static const bool flag = __builtin_cpu_supports("avx2");
  return flag;
#elif defined(TS_SIMD_AVX2)
  return true;   // Compiled with -mavx2; always available
#else
  return false;
#endif
}

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

// ---------- 256-bit operations (AVX2, runtime dispatch) ----------
//
// These functions are annotated with target("avx2") so the compiler
// emits AVX2 instructions even when the translation unit isn't compiled
// with -mavx2.  They are ONLY called after cpu_has_avx2() returns true.

#if TS_HAS_AVX2_DISPATCH || defined(TS_SIMD_AVX2)

using v256 = __m256i;

TS_TARGET_AVX2
inline v256 loadu256(const uint64_t* p) {
  return _mm256_loadu_si256(reinterpret_cast<const v256*>(p));
}
TS_TARGET_AVX2
inline void storeu256(uint64_t* p, v256 v) {
  _mm256_storeu_si256(reinterpret_cast<v256*>(p), v);
}
TS_TARGET_AVX2
inline v256 and256(v256 a, v256 b) { return _mm256_and_si256(a, b); }
TS_TARGET_AVX2
inline v256 or256(v256 a, v256 b)  { return _mm256_or_si256(a, b); }
TS_TARGET_AVX2
inline v256 zero256() { return _mm256_setzero_si256(); }
TS_TARGET_AVX2
inline v256 set1_64_256(uint64_t x) {
  return _mm256_set1_epi64x(static_cast<long long>(x));
}

TS_TARGET_AVX2
inline uint64_t hor_or256(v256 v) {
  __m128i lo = _mm256_castsi256_si128(v);
  __m128i hi = _mm256_extracti128_si256(v, 1);
  __m128i combined = _mm_or_si128(lo, hi);
  uint64_t tmp[2];
  _mm_storeu_si128(reinterpret_cast<__m128i*>(tmp), combined);
  return tmp[0] | tmp[1];
}

#endif // AVX2

// ---------- Convenience: SIMD any_hit reduction ----------
//
// Computes OR( a[s] & b[s] ) for s in [0, n_states).

#if TS_HAS_AVX2_DISPATCH || defined(TS_SIMD_AVX2)
TS_TARGET_AVX2
inline uint64_t any_hit_reduce_avx2(const uint64_t* a, const uint64_t* b,
                                     int n_states) {
  v256 acc = zero256();
  int s = 0;
  for (; s + 4 <= n_states; s += 4) {
    acc = or256(acc, and256(loadu256(&a[s]), loadu256(&b[s])));
  }
  uint64_t result = hor_or256(acc);
  // Scalar tail (0-3 elements)
  for (; s < n_states; ++s) {
    result |= (a[s] & b[s]);
  }
  return result;
}
#endif

inline uint64_t any_hit_reduce(const uint64_t* a, const uint64_t* b,
                                int n_states) {
#if TS_HAS_AVX2_DISPATCH
  if (cpu_has_avx2()) return any_hit_reduce_avx2(a, b, n_states);
#elif defined(TS_SIMD_AVX2)
  return any_hit_reduce_avx2(a, b, n_states);
#endif

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

#if TS_HAS_AVX2_DISPATCH || defined(TS_SIMD_AVX2)
TS_TARGET_AVX2
inline uint64_t any_hit_reduce3_avx2(const uint64_t* clip, const uint64_t* a,
                                      const uint64_t* b, int n_states) {
  v256 acc = zero256();
  int s = 0;
  for (; s + 4 <= n_states; s += 4) {
    acc = or256(acc, and256(loadu256(&clip[s]),
                            or256(loadu256(&a[s]), loadu256(&b[s]))));
  }
  uint64_t result = hor_or256(acc);
  for (; s < n_states; ++s) {
    result |= (clip[s] & (a[s] | b[s]));
  }
  return result;
}
#endif

inline uint64_t any_hit_reduce3(const uint64_t* clip, const uint64_t* a,
                                 const uint64_t* b, int n_states) {
#if TS_HAS_AVX2_DISPATCH
  if (cpu_has_avx2()) return any_hit_reduce3_avx2(clip, a, b, n_states);
#elif defined(TS_SIMD_AVX2)
  return any_hit_reduce3_avx2(clip, a, b, n_states);
#endif

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

#if TS_HAS_AVX2_DISPATCH || defined(TS_SIMD_AVX2)
TS_TARGET_AVX2
inline uint64_t or_reduce_avx2(const uint64_t* a, int n_states, int start) {
  v256 acc = zero256();
  int s = start;
  for (; s + 4 <= n_states; s += 4) {
    acc = or256(acc, loadu256(&a[s]));
  }
  uint64_t result = hor_or256(acc);
  for (; s < n_states; ++s) {
    result |= a[s];
  }
  return result;
}
#endif

inline uint64_t or_reduce(const uint64_t* a, int n_states, int start = 0) {
#if TS_HAS_AVX2_DISPATCH
  if (cpu_has_avx2()) return or_reduce_avx2(a, n_states, start);
#elif defined(TS_SIMD_AVX2)
  return or_reduce_avx2(a, n_states, start);
#endif

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

// ---------- Fitch node scoring: intersect/union with broadcast masks ----------
//
// Computes: out[s] = (left[s] & right[s] & ai) | ((left[s] | right[s]) & nu)
// where ai = broadcast(any_intersect), nu = broadcast(needs_union).
// This is the Fitch downpass node state computation.

#if TS_HAS_AVX2_DISPATCH || defined(TS_SIMD_AVX2)
TS_TARGET_AVX2
inline void fitch_combine_avx2(const uint64_t* left, const uint64_t* right,
                                uint64_t* out, int n_states,
                                uint64_t any_intersect, uint64_t needs_union) {
  v256 ai = set1_64_256(any_intersect);
  v256 nu = set1_64_256(needs_union);
  int s = 0;
  for (; s + 4 <= n_states; s += 4) {
    v256 l = loadu256(&left[s]);
    v256 r = loadu256(&right[s]);
    v256 isect = and256(l, r);
    v256 uni = or256(l, r);
    storeu256(&out[s], or256(and256(isect, ai), and256(uni, nu)));
  }
  for (; s < n_states; ++s) {
    uint64_t isect = left[s] & right[s];
    uint64_t uni = left[s] | right[s];
    out[s] = (isect & any_intersect) | (uni & needs_union);
  }
}
#endif

inline void fitch_combine(const uint64_t* left, const uint64_t* right,
                           uint64_t* out, int n_states,
                           uint64_t any_intersect, uint64_t needs_union) {
#if TS_HAS_AVX2_DISPATCH
  if (cpu_has_avx2()) {
    fitch_combine_avx2(left, right, out, n_states, any_intersect, needs_union);
    return;
  }
#elif defined(TS_SIMD_AVX2)
  fitch_combine_avx2(left, right, out, n_states, any_intersect, needs_union);
  return;
#endif

#if defined(TS_SIMD_SSE2) || defined(TS_SIMD_NEON)
  v128 ai = set1_64(any_intersect);
  v128 nu = set1_64(needs_union);
  int s = 0;
  for (; s + 2 <= n_states; s += 2) {
    v128 l = loadu128(&left[s]);
    v128 r = loadu128(&right[s]);
    v128 isect = and128(l, r);
    v128 uni = or128(l, r);
    storeu128(&out[s], or128(and128(isect, ai), and128(uni, nu)));
  }
  for (; s < n_states; ++s) {
    uint64_t isect = left[s] & right[s];
    uint64_t uni = left[s] | right[s];
    out[s] = (isect & any_intersect) | (uni & needs_union);
  }
#else
  for (int s = 0; s < n_states; ++s) {
    uint64_t isect = left[s] & right[s];
    uint64_t uni = left[s] | right[s];
    out[s] = (isect & any_intersect) | (uni & needs_union);
  }
#endif
}

} // namespace simd
} // namespace ts

#endif // TS_SIMD_H
