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
| T-242 | P1 | PARKED (C, GHA 23545987517) | — | **[Bug?] Agnarsson2004 IW search quality regression.** 230 runs, only 5 hit best score (2% hit rate). User reports "1 trees in memory: 1 sampled, each with score 50.1872 (k = 5.62)". May indicate search regression or IW landscape difficulty. | From a.20. Investigate whether this is a genuine regression or expected IW behaviour. |




### Shiny App

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-232 | P2 | PARKED (D, GHA 23543699366†) | — | **[Shiny] "Tips to show" input bounces back on decrement.** Clicking "down" arrow resets to previous value (e.g. 54 for Sun dataset). | From a.013. Fix committed. GHA failed pre-T-214 fix; awaiting re-validation on current HEAD (run 23547582438 queued). |
| T-240 | P2 | PARKED (D, GHA 23544604214†) | — | **[Shiny] Pool suboptimal filter not applied when changed mid-search.** After search, changing "Keep if suboptimal by" from ≤6 to ≤2 or ≤0 doesn't filter existing pool trees. | From a.17. Fix committed. GHA failed pre-T-214; awaiting re-validation. |
| T-239 | P3 | PARKED (D, GHA 23545538742†) | — | **[Shiny] Cluster consensus: highlight edges unique to a cluster.** Heatmap colouring: emphasize "unique to cluster" vs "in 5/6 clusters". Agnarsson (6 clusters) is testbed. | From a.07. Feature committed. GHA failed pre-T-214; awaiting re-validation. |
| T-241 | P3 | PARKED (D, GHA 23545261957†) | — | **[Shiny] Show cluster assignment next to tree selector.** Add "(cluster X)" in cluster colour after "Tree to plot" label. | From a.19. Feature committed. GHA failed pre-T-214; awaiting re-validation. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-243 | P2 | PARKED (E, GHA 23581391502) | — | **Merge `feature/hot-loop-opt` to `cpp-search`.** FlatBlock struct, flat EW indirect functions, TBR prefetch. Confirmed 1.4% speedup at 180 tips on Hamilton (median 11.538→11.360s, p=0.001, n=10). | PR #230 open. Fixed pre-existing Rd/spelling issues; re-dispatched GHA. |
| T-248 | P3 | ASSIGNED (E) | — | **SA (annealing) phase tuning for large preset.** Hamilton benchmark (T-244) shows SA is 7.4% of time with 14% hit rate (0.8 steps/s vs ratchet 4.5). annealCycles=3, annealPhases=5 may be overtuned. Reducing could save ~1.2s/rep → 1 extra replicate per ~17s. | Benchmark on mbank_X30754 with annealCycles=1 or annealCycles=0 to quantify. |
| T-245 | P3 | OPEN | T-243 | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |
| T-246 | P3 | OPEN | T-243 | **AVX2 runtime dispatch for Fitch bit ops.** Widen `ts_simd.h` from SSE2 (128-bit) to AVX2 (256-bit) with runtime detection (`__builtin_cpu_supports("avx2")`) and SSE2 fallback. Estimated 5–10% on datasets with many states or character blocks. | EPYC 7702 supports AVX2. Can be done independently of T-245. Less invasive than batching. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-25 by F (focus 6: R↔C++ interface). Clean — no bugs. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-25 by F (round 20). Cleaned stale entries, updated PR refs. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-25 by F (round 20 triage). Open PRs: #216 (native-search), #213 (CID-consensus). #226/#227 merged, #215/#222 CLOSED. #210 (draft cpp-search→main). #178/#106 stale — consider closing. |

