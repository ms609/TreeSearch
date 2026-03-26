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
| T-243 | P2 | PARKED (E, GHA 23582386358) | — | **Merge `feature/hot-loop-opt` to `cpp-search`.** FlatBlock struct, flat EW indirect functions, TBR prefetch. Confirmed 1.4% speedup at 180 tips on Hamilton (median 11.538→11.360s, p=0.001, n=10). | PR #230 open. Fixed pre-existing Rd/spelling issues (TREE's, speedup, ratchetTaper/annealCycles usage); re-dispatched GHA. |

| T-245 | P3 | OPEN | T-243 | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |
| T-246 | P3 | OPEN | T-243 | **AVX2 runtime dispatch for Fitch bit ops.** Widen `ts_simd.h` from SSE2 (128-bit) to AVX2 (256-bit) with runtime detection (`__builtin_cpu_supports("avx2")`) and SSE2 fallback. Estimated 5–10% on datasets with many states or character blocks. | EPYC 7702 supports AVX2. Can be done independently of T-245. Less invasive than batching. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-249 | P3 | ASSIGNED (F) | — | **Rerun TNT comparison on current cpp-search HEAD.** Round 2 data (120s, 3 seeds) predates NNI warmup, ratchet tuning, biased Wagner, outer cycle loop, and SA tuning. Update `round2_hard.csv` with current engine to measure remaining gaps. | Existing gaps: Geisler2001 +5, Wortley2006/Zhu2013 +3, Zanol2014 +2. TNT converges 3–5× faster on most datasets. |
| T-251 | P3 | OPEN | T-249 | **TNT trajectory analysis on gap datasets.** Run TNT with verbose logging on Geisler2001, Zhu2013, Wortley2006 (largest score gaps). Focus on **candidates-per-improvement**: how many TBR/SPR candidates does TNT evaluate per score improvement vs TreeSearch? T-250 showed TreeSearch has ~4× raw Fitch throughput, so TNT's 3–5× convergence speed must come from evaluating fewer candidates (better pruning, smarter regraft ordering, or more effective perturbation). Also compare ratchet escape magnitudes and phase composition. | Depends on T-249 to confirm which gaps persist. |
| T-252 | P3 | OPEN | T-251 | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | Uses `benchmark_mbank_sample()` in `bench_framework.R`. |
| T-253 | P3 | OPEN | T-251, T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | Depends on T-249 (updated gaps) and T-252 (broader baseline). |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-25 by F (focus 6: R↔C++ interface). Clean — no bugs. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-26 by E (round 21). Closed 4 Shiny tasks (re-validated by GHA 23547582438). Updated stale GHA notes. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-26 by E (via S-COORD). Open PRs: #230 (hot-loop-opt, GHA pending), #216 (native-search, stale GHA), #213 (CID-consensus, WORKTREE). #210 (draft cpp-search→main). #178/#106 stale+CONFLICTING — recommend closing. |

