# TreeSearch Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- When a task is being developed in a **git worktree**, set its status to
  `WORKTREE (name)` where *name* is the worktree directory (e.g.
  `WORKTREE (TS-CID-cons)`). This distinguishes human/long-running worktree
  work from agent assignments and prevents double-claiming.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md` (see workflow in AGENTS.md).
- Tasks awaiting GHA results: `PARKED (<Letter>, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (<Letter>)`.
  S-COORD cleans these up after merge.
- Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
  completed, reset it to OPEN. Their effective priority is dynamic:
  - ≥6 OPEN specific tasks → standing tasks are P3
  - 3–5 OPEN specific tasks → standing tasks are P2
  - <3 OPEN specific tasks → standing tasks are P1

---

## Active Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-150 | P2 | WORKTREE (TS-CID-cons) | — | **CID-optimal consensus tree search** | PR #213 open to cpp-search. |
| T-204 | P2 | PR #216 (B) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | On `feature/native-search`. GHA run 23495097795. |



### Bugs

(no open bugs)




### Shiny App

(no open tasks)

### Performance Optimization — TBR per-evaluation overhead (from T-260 VTune)

T-260 VTune profiling found **37.8% of TBR time** is non-scoring overhead.
Three tasks below address the top hotspots. Combined potential: ~16–19% wall
time reduction. See `dev/benchmarks/vtune_tbr_analysis.md` for full data.

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-263 | P2 | PR #231 (F) | — | **Selective StateSnapshot save/restore in TBR.** `StateSnapshot::save()` and `::restore()` in `ts_tbr.cpp:216-246` memcpy the **entire** prelim, final_, down2, subtree_actives, local_cost, and postorder arrays before each candidate evaluation and restore on rejection. VTune: **14.6% of TBR time** (4.53s/30s). At 88 tips, each save/restore copies ~190 KB; at 180 tips, ~380+ KB. Most of this data is untouched — only nodes on the path from clip/regraft points to root are modified by `apply_tbr_move()` + `full_rescore()`. **Task:** (1) In `apply_tbr_move()`, track which node indices are modified (the clip subtree, the regraft path, and their ancestors to root). (2) Replace bulk memcpy with selective save/restore of only dirty node entries. (3) For `postorder`, either save/rebuild (cheap at ~175 ints) or track whether it changed. (4) Benchmark: must show measurable wall-time improvement on Dikow2009 (88t) and ideally mbank_X30754 (180t). **Alternative approach:** eliminate the snapshot entirely — after rejection, call `full_rescore()` to recompute from scratch. This trades snapshot cost (14.6%) for an extra `full_rescore` on rejection (~29% × rejection_rate). Only viable if T-261+T-262 make `full_rescore` very cheap. Benchmark both approaches. | Medium-high effort. Most impactful single change (~10–12% savings). Needs careful correctness testing — any node missed in the dirty set produces silent score corruption. Feature branch `feature/selective-snapshot`. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |
| T-246 | P3 | PR #233 (F) | — | **AVX2 runtime dispatch for Fitch bit ops.** Widen `ts_simd.h` from SSE2 (128-bit) to AVX2 (256-bit) with runtime detection (`__builtin_cpu_supports("avx2")`) and SSE2 fallback. Estimated 5–10% on datasets with many states or character blocks. | EPYC 7702 supports AVX2. Can be done independently of T-245. Less invasive than batching. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-252 | P3 | OPEN | — | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | Uses `benchmark_mbank_sample()` in `bench_framework.R`. |
| T-253 | P3 | OPEN | T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | T-249 complete: EW gaps are 0–7 steps (mean 2.2) across 11 hard datasets at 120s. 5 datasets optimal. Remaining blocker: T-252 (broader baseline). **NB:** always compare like-for-like scoring (Fitch vs Fitch); Brazeau scores are inherently higher. |

### Strategy Tuning (from T-251 trajectory analysis)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-257 | P3 | PARKED (F, GHA 23607823258) | — | **Post-ratchet sectorial search pass.** Add a second sectorial search pass after ratchet in the pipeline: [XSS+RSS+CSS → Ratchet → XSS+RSS+CSS → TBR] instead of [XSS+RSS+CSS → Ratchet → Drift → TBR]. TNT interleaves sectorial search throughout each replicate; this is a lightweight approximation. | T-256 found extra sectorial rounds don't improve scores, but `nodrift_3x` config was best (mean gap 4.9 vs 5.3) due to more replicates. A post-ratchet sectorial pass may still help if it's cheap enough to not cut into replicate count. Needs careful cost/benefit analysis. |


### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-26 (focus 9: Wagner & addition trees). ts_wagner.h/.cpp (595 lines) + ts_constraint.h/.cpp (880 lines) reviewed. Latent stale-reference in impose_one_pass() noted (negligible severity, mitigated by retry loops). 902 constraint tests + 80 adversarial tests pass. No bugs filed. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-26 round 29. Closed T-242 (display bug, already fixed). T-257 GHA doc mismatch (tests pass). 2 unblocked OPEN → standing at P1. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-26 by F. 6 open PRs. Updated 4 feature PRs with cpp-search base merges (all clean): #233 (AVX2, 23 behind → merged), #231 (selective-snapshot, 25 behind → merged), #216 (native-search, 36 behind → merged), #213 (CID consensus, 26 behind → merged). #210 (cpp-search→main, MERGEABLE). #178 (stale, CONFLICTING, Aug 2025 — recommend close). No review comments on any PR. CI re-triggered on all updated PRs. |

