# AVX2 reduce-codegen gate: packing reopened a ~9% per-candidate win — 2026-07-16

**Status:** FINDING recorded; **no code landed** (the portable source capture did not work locally
— reverted). Resolution deferred to a Hamilton EPYC `-msse2` vs `-mavx2` A/B (fold into the packbuild
rebuild after the 5432 job `17887591` completes — do NOT rebuild packbuild while it runs).

## The finding

The package builds `-O2 -msse2` (no `-mavx2`; `~/.R/Makevars.win` CXXFLAGS). The AVX2 reduce
(`any_hit_reduce_avx2`, `__attribute__((target("avx2")))`, ts_simd.h) is therefore a **non-inlined
call boundary** on every candidate. Packing (landed today, e8a17e60) cut hot blocks to
`n_states`=2-4, where that fixed call/setup/spill is a big fraction of a tiny block body — so
inlining the reduce now pays where it did not pre-packing.

**Whole-program `-mavx2` A/B (local, back-to-back, cutoff-tight/prefetch-off, min-of-15):**

| dataset | `-msse2` ns/cand | `-mavx2` ns/cand | speedup |
|---|---|---|---|
| project970 (1844c) | 60.57 | 55.25 | **1.088× (8.8%)** |
| project2668 (1227c) | 34.72 | 31.90 | **1.081× (8.1%)** |

Reproducible, not session drift. This is the **largest exact per-candidate lever found** post-packing.
Contrast the pre-packing prior (T-P5f whole-program codegen: −1.2%) — packing *reopened* it, exactly
as the lever-hunt `map:simd-reduce` predicted. (Injected via `src/Makevars.win` `PKG_CPPFLAGS=-mavx2`
+ `touch src/*.cpp` to force a flag-correct recompile; the env `PKG_CPPFLAGS` and a `~/.R/Makevars.win`
CXXFLAGS edit both silently did NOT take — R rebuilds the flag string.)

## Why the portable source capture did NOT work here (reverted)

Attempted the map's "hoist target(avx2) to the scorer": an `always_inline` core shared by a default
wrapper and a `TS_TARGET_AVX2` clone of `fitch_indirect_cached_flat_x4`, dispatched by
`cpu_has_avx2()`. REROOT ns was **flat** (58.9→58.8) even when the clone was force-called.

Diagnosis (by elimination — a standalone `__builtin_cpu_supports` probe won't link on this rtools
box, an unrelated `ld` quirk): **`cpu_has_avx2()` is almost certainly `false` on this MinGW/Windows
machine** (`__builtin_cpu_supports("avx2")` is known-unreliable there — the CPU-feature init may not
run). If it were true, the `-msse2` baseline would already call `any_hit_reduce_avx2` and inlining it
(the clone) would have captured the win. It was flat because the clone's *internal* `any_hit_reduce`
also took the SSE2 path (`cpu_has_avx2()` false), so the `target(avx2)` wrapper ran SSE2 code. The
whole-program `-mavx2` win is therefore SSE2→AVX2 via *unconditional* codegen (which ignores the
runtime check).

**Implication (strong hypothesis, verify on Linux/CI):** the shipped **Windows** binary likely runs
the SSE2 reduce, never its AVX2 helpers. On Linux/EPYC `__builtin_cpu_supports` works, so
`cpu_has_avx2()` is true there and the `-msse2` build already uses the AVX2 helpers (with the call
boundary) — a *different* baseline, where the target-clone (killing the boundary) is the right lever.

## Two portable capture paths (platform-dependent) — resolve on Hamilton

1. **If Windows silently runs SSE2** (`cpu_has_avx2()` false): fix CPU detection on MinGW (manual
   CPUID / Windows API instead of `__builtin_cpu_supports`) → Windows uses the AVX2 helpers. Portable,
   shippable, no `-mavx2` build needed. (Then the target-clone adds the boundary-elimination on top.)
2. **Where `cpu_has_avx2()` is true** (Linux/EPYC, most non-MinGW): a `target(avx2)` scorer clone
   dispatched by `cpu_has_avx2()` kills the per-candidate call boundary. (My local test could not
   validate this because the runtime check is false here.)

## Next step (deferred, non-blocking)

Hamilton EPYC A/B once `17887591` (5432) finishes: rebuild packbuild `-msse2` vs `-mavx2`, measure
raw-`tbr_search` wall + full `MaximizeParsimony` wall (min-of-runs) on project970/510/5432, and print
`cpu_has_avx2()` on Linux. That settles (a) the achievable end-to-end win, (b) which portable path
applies, and (c) the Windows-SSE2 hypothesis. Consistent with [[hamilton-gcc-module-noop]] (codegen
wins seen on Hamilton, deployment-gated).
