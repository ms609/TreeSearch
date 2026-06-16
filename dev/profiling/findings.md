# Profiling findings

One row per verified optimisation opportunity, in `to-do.md` paste-ready
format. A finding only lands here if an isolated `std::chrono` micro-bench
reproduces the predicted delta.

Tags:
- `[Port]` — R loop on the hot path that should move to C++.
- `[Optimise]` — C++ change with verified expected speedup.
- `[AT-LIMIT]` — function is at a hardware ceiling; record so the rotation
  skips it in future rounds.

| ID-suggest | P? | Status | Depends | Headline | Detail (% time, mechanism, verified Δ, micro-bench path) |
|------------|----|--------|---------|----------|---------------------------------------------------------|
| T-300 | P1 | DONE | — | [Optimise] `full_rescore` after accepted TBR move (ts_tbr.cpp:1138): replace with incremental rescore | LANDED (commits f531bbcd EW + 014ccdea NA dirty-set). 19.2 % of NA-path DLL CPU; 15.2 % wall speedup on Zhu2013 NA (3.88→3.29 s). |

## Round 3 (2026-06-16) — standard-Fitch TNT-parity path (Zhu2013 `-`→`?`, auto→thorough)

DLL self-CPU total 2.70 s. Names resolved via `nm` (VTune reporter shows
`func@0x…` for MinGW DWARF; image base stable so addresses map 1:1).

| ID-suggest | P? | Status | Depends | Headline | Detail (% time, mechanism, verified Δ, micro-bench path) |
|------------|----|--------|---------|----------|---------------------------------------------------------|
| T-S3a | P3 | **DONE (verified)** | — | [Optimise] Per-clip allocation churn in TBR helpers — Tier 1 | **IMPLEMENTED 2026-06-16** (uncommitted, in working tree). Converted per-clip scratch vectors to `static thread_local` + clear/assign: `fitch_incremental_uppass` `dirty` (`std::vector<bool>`→`char`, allocated EVERY clip — ts_fitch.cpp:225), `collect_main_edges`/`collect_subtree_edges` DFS stacks, `compute_from_above` preorder+stack (ts_tbr.cpp). Per-thread-safe (each search thread owns its TreeState; none re-entrant). **VERIFIED:** Zhu2013 `-`→`?` score 627 unchanged (10/10 runs); ts- suite 3061 assertions, 0 fail; wall-clock 4.00 s→3.84 s = **~4.0%** (non-overlapping medians, TS_REPS=8 TS_SEED=1, identical build flags). Remaining Tiers 2 (reusable dedup table for `seen_vp_hashes`, ~2-4%) and 3a (gate `validate_topology` under NDEBUG, ~2-3%) NOT done — see design plan; Tier 3a is a safety/policy call (removes a release-build invariant check). |
| T-S3b | — | AT-LIMIT | — | [AT-LIMIT] `simd::any_hit_reduce(3)_avx2` (21 % DLL) | Core Fitch reduce. Disasm of `hor_or256` confirms GCC already emits register-only horizontal reduce (vextracti128/vpsrldq/vpor/vmovq) — the store-reload anti-pattern is elided. Compiler-optimal at -O2. No win. |
| T-S3c | — | AT-LIMIT | — | [AT-LIMIT] `uppass_node` scalar state-update loop (13 % DLL) | The update loop (ts_fitch.cpp:54-61) is scalar vs the vectorised `fitch_combine` in downpass. Micro-bench `dev/profiling/microbench/bench_uppass_combine.cpp`: AVX2 version bit-identical (value+changed flag, 0 mismatches) but only **1.22×** at n_states=4, and the 4-wide path does NOT trigger for 2-state (binary) morphology → ~1 % wall-clock. Not worth the incremental-uppass correctness risk (see fitch-scoring memory: dirty-flag invariant is delicate). |
| — | — | NOTE | — | [Strategic] Standard-Fitch is bookkeeping/strategy-bound, not scoring-bound | Per-candidate scoring at AVX2/compiler limit. Parity levers: (a) reduce per-clip O(n) bookkeeping — `build_postorder_prealloc` 5.2 % (rebuilt per clip+accept) + incremental down/uppass 6.4 %; (b) ratchet (63 %) evaluation economy. Aligns with `.positai/plans/2026-03-21-tnt-outperformance-analysis.md` (strategy > code). |
