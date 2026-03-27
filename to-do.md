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
| T-150 | P2 | PARKED (E, GHA 23636944848) | — | **CID-optimal consensus tree search** | PR #213 open to cpp-search. SPIC method added (commit 6636924c); GHA 23636944848 **FAILED** — codoc mismatch in `InfoConsensus.Rd`. Fix: regenerate Rd (`roxygen2::roxygenise(load_code=roxygen2::load_installed)`), commit, re-dispatch. |
| T-204 | P2 | PR #216 (B) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | On `feature/native-search`. GHA run 23495097795. |
| T-266 | P2 | PR #235 (A) | — | **Taxon pruning-reinsertion perturbation.** Drop ~10% of leaves, TBR-optimize backbone, Wagner-reinsert, TBR-polish. Random + instability-weighted selection. Disabled by default. | On `feature/prune-reinsert` (worktree `TS-PruneRI`). 44 tests. GHA passed (run 23636145497). |


### Bugs

(no open bugs)


### Shiny App

(no open tasks)

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-252 | P3 | ASSIGNED (F) | — | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | Uses `benchmark_mbank_sample()` in `bench_framework.R`. |
| T-253 | P3 | OPEN | T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | T-249 complete: EW gaps are 0–7 steps (mean 2.2) across 11 hard datasets at 120s. 5 datasets optimal. Remaining blocker: T-252 (broader baseline). **NB:** always compare like-for-like scoring (Fitch vs Fitch); Brazeau scores are inherently higher. |

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-269 | P3 | OPEN | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-268 | P2 | ASSIGNED (F) | — | **Branch housekeeping.** Human has merged several branches on GitHub. Pull cpp-search, refresh local branches, delete stale worktrees. | From u.2. |



### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-26 (focus 9: Wagner & addition trees). ts_wagner.h/.cpp (595 lines) + ts_constraint.h/.cpp (880 lines) reviewed. Latent stale-reference in impose_one_pass() noted (negligible severity, mitigated by retry loops). 902 constraint tests + 80 adversarial tests pass. No bugs filed. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-26 by E (round 5: 180-tip large-preset benchmarks on Hamilton HPC, T-244/T-248 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 31 by A. T-266 PR #235 opened. T-150 GHA failed (InfoConsensus.Rd codoc). T-270 completed. T-272 filed. 2 unblocked OPEN specific tasks (T-245, T-269) → standing at P1. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 by A (round 31). Merged cpp-search into #235 (prune-reinsert, clean) and #216 (native-search, clean). #213 (cid-consensus) has ts_tbr.cpp conflict (CID changes vs T-263 snapshot opt) — needs human/E resolution. #178 CLOSED (stale 2025 draft). Open: #213 (GHA failing, merge conflict), #216 (clean), #235 (clean). #210 (cpp-search→main). |
