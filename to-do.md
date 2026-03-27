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
| T-150 | P2 | PARKED (F, GHA 23646972365) | — | **CID-optimal consensus tree search** | PR #213. GHA 23645305415 failed: bare Consensus() in vignettes (not on search path). Fixed 23d93f2e. Re-dispatched GHA 23646972365. |
| T-204 | P2 | PARKED (F, GHA 23647123007) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | On `feature/native-search`. GHA 23644617599 failed: 98 deprecation warnings (Morphy/PhyDat2Morphy) in tests. Fixed 41ae291a (suppressWarnings). Re-dispatched GHA 23646124106. |


### Bugs

(no open bugs)


### Shiny App

(no open tasks)

### API / UX

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-276 | P3 | OPEN | — | **Print convergence summary to console after `MaximizeParsimony()`.** The Shiny app displays a convergence evaluation report; same info should print to the R console when running headlessly. Note: integrates with B's Chao-estimate work (T-204 native-search) once that merges. | u.571. Print n_replicates, best_score, n_MPTs, last_improved_rep, time elapsed, and convergence indicator (consensus_stable / plateau_stop / timed_out / perturb_stop). **S-RED E (2026-03-27):** `DrivenResult::perturb_stop` field now added to ts_driven.h/cpp (previously missing). Use it directly. |
| T-277 | P3 | PARKED (B, GHA 23646841244) | — | **ScoreSpectrum(): Chao1-style landscape coverage estimator** | On `feature/score-spectrum`. Exports `ScoreSpectrum()`: accepts `multiPhylo` (with new `replicate_scores` attribute) or raw numeric vector; returns Good-Turing coverage + Chao1 richness. C++ side: `DrivenResult::replicate_scores` vector (serial + parallel paths). Shiny: coverage note appended to confidence text. 8 Tier-1 tests. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-252 | P3 | PARKED (F, SLURM 16599543) | — | **Hamilton MorphoBank training-set benchmarking.** Run TreeSearch on fixed 25-matrix training sample at 30s/60s/120s budgets. Baseline current engine across size/complexity spectrum before any strategy tuning. | `t252_v2.sh` submitted (job 16599543). Previous failure was httpuv/shiny not building in fresh lib — fixed by using ts-bench/lib-baseline for deps and only installing fresh TreeSearch to lib-t252. |
| T-253 | P3 | OPEN | T-252 | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | T-249 complete: EW gaps are 0–7 steps (mean 2.2) across 11 hard datasets at 120s. 5 datasets optimal. Remaining blocker: T-252 (broader baseline). **NB:** always compare like-for-like scoring (Fitch vs Fitch); Brazeau scores are inherently higher. |

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-269 | P3 | OPEN | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. **Design note (u.005):** Full TNT-style interleaving IS architecturally supported now — setting `outerCycles = ratchetCycles` achieves one sector pass per ratchet cycle. T-257 first validated that any post-ratchet sectorial helps at all (merged). T-269 benchmarks the fine-grained variant to determine whether the per-cycle overhead is worth enabling in presets. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |


### Housekeeping

(no open tasks)



### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 4 by E (ts_driven.cpp, 1054 lines). Bugs found and fixed: (1) `unsuccessful_reps` not reset on fuse improvement — perturb-stop could fire prematurely when fusing is still productive; (2) `DrivenResult::perturb_stop` flag missing — needed by T-276; (3) stale NNI-perturb constraint comment. Also noted: consensus_constrain calls extract_consensus_splits every rep when n_unanimous=0 (performance, not correctness). Next: review ts_parallel.cpp (parallel path — not reviewed since parallelism was added). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 33 by F. T-273 completed. T-275 filed and completed by B (prune-reinsert guard). T-277 filed (ScoreSpectrum, ASSIGNED B). Unblocked OPEN: T-245, T-269, T-274, T-275→done, T-276, T-253, T-277 → 6 specific OPEN → standing at P3. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 by F (round 36). #216 (T-204 native-search): GHA 23644617599 failed (98 deprecation warns), fixed 41ae291a, re-dispatched 23646124106. #213 (T-150 cid-consensus): GHA 23644644554 failed (Consensus case), fixed 0c1f0fe9, re-dispatched 23645305415. Open PRs: #213, #216, #210 (cpp-search→main DRAFT). |
