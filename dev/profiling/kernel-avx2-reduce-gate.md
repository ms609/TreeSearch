# AVX2 reduce-codegen: real per-candidate win, MODEST end-to-end — RESOLVED on EPYC 2026-07-17

**Verdict:** the shipped `-msse2` build pays a per-block AVX2-call-boundary on EPYC; a `-mavx2`
(or `-march=native`) build removes it for a **1.26–1.67× per-candidate** win — but that dilutes to
**~1–10% end-to-end** `MaximizeParsimony` wall, because the default `auto` search is dominated by
non-scoring machinery (ratchet/sectorial/tree-ops), not the reroot reduce. **Banked as a
`-march=native` Hamilton deployment build** (free, machine-specific, off the CRAN release). Not worth
a portable source change. Confirms the mission: the per-move kernel is not the end-to-end bottleneck.

## Per-candidate A/B on EPYC (job 17908514; `avx2ab/lib_base` vs `lib_avx2`, kern_real min-of-15)

Linux `__builtin_cpu_supports("avx2")` = **true** (unlike the MinGW box).

| dataset | blocks | `-msse2` ns | `-mavx2` ns | speedup |
|---|---|---|---|---|
| project970 (1844c) | ~29 | 142.2 | 85.0 | **1.67×** |
| project510 (2954c) | ~46 | 116.0 | 84.3 | **1.38×** |
| project2668 (1227c)| ~20 | 67.5 | 53.6 | **1.26×** |
| project5432 (189c) | 3 | 22.6 | 21.0 | 1.07× |

Scales with block count. **Mechanism:** `cpu_has_avx2()` is true on EPYC, so the `-msse2` build
calls `any_hit_reduce_avx2` (`__attribute__((target("avx2")))`, ts_simd.h) **across a target
boundary it cannot inline — once per block**, plus a `hor_or256` horizontal reduce each block. On
Zen2 that call/spill is costly and compounds over ~29 blocks × millions of candidates. `-mavx2`
compiles the scorer in an AVX2 context so the helper inlines. **Windows masked this:** there
`cpu_has_avx2()` is *false*, so the base build used inline SSE2 (no boundary) and the gap looked like
only ~9%.

## End-to-end MaximizeParsimony wall (job 17908602; same libs, fixed-seed maxReplicates=4 → identical trajectory)

| dataset | base wall | avx2 wall | end-to-end | reach |
|---|---|---|---|---|
| project970 | 3357.1 s | 3321.0 s | **1.01×** | 11704 = 11704 |
| project510 | 847.8 s | 768.3 s | **1.10×** | 18345 = 18345 |

**The per-candidate win does NOT propagate.** Amdahl backs out the scoring fraction: project970
scoring ≈ **3%** of wall (1.67× on 3% → 1.01×); project510 ≈ **34%** (1.38× on 34% → 1.10×). So
~66–97% of `auto`-search wall is non-scoring machinery (ratchet cycles, sectorial decomposition,
tree management, collapse/pool), dataset-dependent. Same dilution that sank lever-6.

Reconciles why **packing** landed 1.2–1.5× end-to-end but AVX2 only ~1–10%: packing shrinks
`total_words` *pervasively* (scoring + edge-set buffer/cache + precompute across every phase),
whereas AVX2 only speeds the reduce *compute* — a smaller slice.

## Windows note (separate, unresolved)

`cpu_has_avx2()` is almost certainly **false on MinGW/Windows** (`__builtin_cpu_supports` unreliable
there; a standalone probe would not link — `ld` quirk). So the shipped **Windows** binary likely runs
the SSE2 reduce and never its AVX2 helpers. Low priority (Windows isn't the compute platform), but a
robust CPU-detect (manual CPUID / Windows API) would let Windows use the AVX2 helpers.

## Banked

- **Hamilton:** `-march=native` deployment build recipe (`/nobackup/pjjg18/tsmarch/tsmarch_build.sh`)
  → `/nobackup/pjjg18/tsmarch/lib`; point search jobs' `lib.loc` / `R_LIBS_USER` at it for the free
  ~1–10% (larger on scoring-bound char-rich data). Machine-specific; **must not** go in the CRAN
  source (`-march=native` breaks non-AVX2 CPUs — the release keeps `-msse2` + runtime dispatch).
- **Not pursued:** the portable `target(avx2)` scorer clone. It *would* work on Linux (the boundary
  is real there), but with scoring at ~3–34% of wall the end-to-end payoff is small and
  dataset-dependent — not worth the cross-platform multi-versioning risk. Recorded, not built.

## Takeaway

Real kernel headroom exists on EPYC, but time-to-optimum is bounded by the **search machinery /
stopping**, not the per-move reduce — the next real wall lever is there (Mission A / over-search),
not further per-candidate scoring. Consistent with [[hamilton-gcc-module-noop]].
