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

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-264 | P0 | PARKED (F, GHA 23600674681) | — | **[Bug] `consensusStableReps = 3` causes catastrophic early termination.** Fix committed to cpp-search (23e9f57b). GHA validation in progress. | From T-249 analysis. Fix by F: removed consensusStableReps from presets (falls back to 0 = disabled). |
| T-242 | P1 | PARKED (C, GHA 23545987517†) | — | **[Bug?] Agnarsson2004 IW search quality regression.** 230 runs, only 5 hit best score (2% hit rate). User reports "1 trees in memory: 1 sampled, each with score 50.1872 (k = 5.62)". May indicate search regression or IW landscape difficulty. | From a.20. GHA failure is stale (pre-T-214); cpp-search passes on 23547582438. Investigation task — GHA status doesn't resolve it. |




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
| T-246 | P3 | PARKED (F, GHA 23598736410) | — | **AVX2 runtime dispatch for Fitch bit ops.** Widen `ts_simd.h` from SSE2 (128-bit) to AVX2 (256-bit) with runtime detection (`__builtin_cpu_supports("avx2")`) and SSE2 fallback. Estimated 5–10% on datasets with many states or character blocks. | EPYC 7702 supports AVX2. Can be done independently of T-245. Less invasive than batching. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-249 | P3 | HAMILTON (F, job 16596844) | — | **Rerun TNT comparison on current cpp-search HEAD.** Round 2 data (120s, 3 seeds) predates NNI warmup, ratchet tuning, biased Wagner, outer cycle loop, and SA tuning. Update `round2_hard.csv` with current engine to measure remaining gaps. **Run on Hamilton.** | Partial results (5 gap datasets, local smoke run): gaps +3 to +6 steps, essentially unchanged from R2. Benchmark infra fixed: score regex handles integers, `run_treesearch()` now uses `MaximizeParsimony()` API with strategy presets. Scripts in `dev/benchmarks/run_t249_*.R`. |
| T-252 | P3 | OPEN | — | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | Uses `benchmark_mbank_sample()` in `bench_framework.R`. |
| T-253 | P3 | OPEN | T-249, T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | Depends on T-249 (updated gaps) and T-252 (broader baseline). |

### Strategy Tuning (from T-251 trajectory analysis)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-257 | P3 | OPEN | — | **Post-ratchet sectorial search pass.** Add a second sectorial search pass after ratchet in the pipeline: [XSS+RSS+CSS → Ratchet → XSS+RSS+CSS → TBR] instead of [XSS+RSS+CSS → Ratchet → Drift → TBR]. TNT interleaves sectorial search throughout each replicate; this is a lightweight approximation. | T-256 found extra sectorial rounds don't improve scores, but `nodrift_3x` config was best (mean gap 4.9 vs 5.3) due to more replicates. A post-ratchet sectorial pass may still help if it's cheap enough to not cut into replicate count. Needs careful cost/benefit analysis. |


### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | ASSIGNED (E) | — | **Standing: Red-team review** | Focus 8: T-264 consensus-stop fix, PR #232 merge, T-255 drift removal. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-26 by E (round 23). Fixed T-248 anneal test stale assertion (annealCycles 3→1). Updated T-255 GHA. Task queue healthy. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-26 by F (S-COORD round 22). Open PRs: #216 (native-search, rebased, CI re-triggered), #213 (CID-consensus, rebased, CI re-triggered). #210 (draft cpp-search→main). #178/#106 stale+CONFLICTING — recommend closing. #230 merged. |

