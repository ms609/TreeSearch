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
| T-242 | P1 | PARKED (C, GHA 23545987517†) | — | **[Bug?] Agnarsson2004 IW search quality regression.** 230 runs, only 5 hit best score (2% hit rate). User reports "1 trees in memory: 1 sampled, each with score 50.1872 (k = 5.62)". May indicate search regression or IW landscape difficulty. | From a.20. GHA failure is stale (pre-T-214); cpp-search passes on 23547582438. Investigation task — GHA status doesn't resolve it. |




### Shiny App

(no open tasks)

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |
| T-246 | P3 | OPEN | — | **AVX2 runtime dispatch for Fitch bit ops.** Widen `ts_simd.h` from SSE2 (128-bit) to AVX2 (256-bit) with runtime detection (`__builtin_cpu_supports("avx2")`) and SSE2 fallback. Estimated 5–10% on datasets with many states or character blocks. | EPYC 7702 supports AVX2. Can be done independently of T-245. Less invasive than batching. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-249 | P3 | OPEN | — | **Rerun TNT comparison on current cpp-search HEAD.** Round 2 data (120s, 3 seeds) predates NNI warmup, ratchet tuning, biased Wagner, outer cycle loop, and SA tuning. Update `round2_hard.csv` with current engine to measure remaining gaps. **Run on Hamilton.** | Partial results (5 gap datasets, local smoke run): gaps +3 to +6 steps, essentially unchanged from R2. Benchmark infra fixed: score regex handles integers, `run_treesearch()` now uses `MaximizeParsimony()` API with strategy presets. Scripts in `dev/benchmarks/run_t249_*.R`. |
| T-252 | P3 | OPEN | — | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | Uses `benchmark_mbank_sample()` in `bench_framework.R`. |
| T-253 | P3 | OPEN | T-249, T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | Depends on T-249 (updated gaps) and T-252 (broader baseline). |

### Strategy Tuning (from T-251 trajectory analysis)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-255 | P2 | PARKED (E, GHA 23591874696) | — | **Reduce drift in default and thorough presets.** T-254 confirmed drift has zero score, MPT, or diversity benefit and costs 10–22% of replicates. Set `driftCycles=0` in default; reduce from 12 to 0 in thorough. Reallocate saved time to extra replicates. | Anneal test fix committed (7dc2ed96). GHA 23591874696 in progress. |
| T-257 | P3 | OPEN | — | **Post-ratchet sectorial search pass.** Add a second sectorial search pass after ratchet in the pipeline: [XSS+RSS+CSS → Ratchet → XSS+RSS+CSS → TBR] instead of [XSS+RSS+CSS → Ratchet → Drift → TBR]. TNT interleaves sectorial search throughout each replicate; this is a lightweight approximation. | T-256 found extra sectorial rounds don't improve scores, but `nodrift_3x` config was best (mean gap 4.9 vs 5.3) due to more replicates. A post-ratchet sectorial pass may still help if it's cheap enough to not cut into replicate count. Needs careful cost/benefit analysis. |
| T-260 | P3 | DONE (E) | — | **Per-evaluation overhead profiling.** T-251 found TreeSearch evaluates 1.5–3.6× fewer rearrangements/second than TNT despite wider SIMD. Profile the per-evaluation overhead: undo stack management (PreallocUndo grow/shrink), incremental scoring setup, collapsed-flag recomputation. Identify the top 2–3 hotspots. | Overlaps T-245 (TBR batching) and T-246 (AVX2) but targets different bottleneck (overhead, not throughput). VTune 2025.10 available locally at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest`. Use software sampling (`-knob sampling-mode=sw`). Use a TBR-only workload at 68+ tips. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-26 by E (focus 7: T-258 intra-fuse, PR #230 hot-loop-opt, T-255 drift elimination). Clean — no bugs. 624 targeted tests pass. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-26 by E (round 23). Fixed T-248 anneal test stale assertion (annealCycles 3→1). Updated T-255 GHA. Task queue healthy. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-26 by F (S-COORD round 22). Open PRs: #216 (native-search, rebased, CI re-triggered), #213 (CID-consensus, rebased, CI re-triggered). #210 (draft cpp-search→main). #178/#106 stale+CONFLICTING — recommend closing. #230 merged. |

