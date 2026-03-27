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
| T-150 | P2 | PARKED (F, GHA 23648267378) | — | **CID-optimal consensus tree search** | PR #213. Fixed InfoConsensus.Rd (added treeSample, screeningTopK, method; fixed maxReplicates/targetHits/maxDrop defaults) and SearchControl.Rd (added scoreTol, plateauReps). commit 5db96f42. GHA 23648267378. |
| T-204 | P2 | PARKED (F, GHA 23648401936) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | Fixed Rd example deprecation warnings: suppressWarnings() around PhyDat2Morphy/UnloadMorphy in GapHandler, MorphyWeights, PhyDat2Morphy, RearrangeEdges, SingleCharMorphy; suppressWarnings(Morphy(...)) in donttest block. commit ec5f419f. GHA 23648401936. |


### Bugs

(no open bugs)


### Shiny App

(no open tasks)

### API / UX

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-277 | P3 | PR #236 (B) | — | **ScoreSpectrum(): Chao1-style landscape coverage estimator** | On `feature/score-spectrum`. Exports `ScoreSpectrum()`: accepts `multiPhylo` (with new `replicate_scores` attribute) or raw numeric vector; returns Good-Turing coverage + Chao1 richness. C++ side: `DrivenResult::replicate_scores` vector (serial + parallel paths). Shiny: coverage note appended to confidence text. 8 Tier-1 tests. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### TNT Comparison & Strategy Learning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-253 | P3 | OPEN | — | **Gap characterization by dataset features.** Correlate TNT-vs-TreeSearch score gaps with dataset features (ntax, nchar, missing %, homoplasy, n_blocks) to identify what *types* of problems TreeSearch is weakest on. Guide targeted strategy improvements. | T-252 complete: 25-matrix training baseline at 30/60/120s downloaded (t252_mbank_*.csv). ≤35t: all converge at 30s. 36-65t: near-optimal. 66-135t: still improving. project4284 (4062t): can't complete 1 replicate. **NB:** always compare like-for-like scoring (Fitch vs Fitch). |

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-269 | P3 | OPEN | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. **Design note (u.005):** Full TNT-style interleaving IS architecturally supported now — setting `outerCycles = ratchetCycles` achieves one sector pass per ratchet cycle. T-257 first validated that any post-ratchet sectorial helps at all (merged). T-269 benchmarks the fine-grained variant to determine whether the per-cycle overhead is worth enabling in presets. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |


### Housekeeping

(no open tasks)



### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 5 by F (ts_parallel.cpp, 589 lines). Bug found and fixed: `result.perturb_stop` never initialized (UB) and not set to `true` when perturb-stop fires in parallel path — serial path had both correct. commit 1a640b73. GHA 23648703841. Next: ts_tbr.cpp or ts_ratchet.cpp (core search modules not reviewed recently). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 33 by F. T-273 completed. T-275 filed and completed by B (prune-reinsert guard). T-277 filed (ScoreSpectrum, ASSIGNED B). Unblocked OPEN: T-245, T-269, T-274, T-275→done, T-276, T-253, T-277 → 6 specific OPEN → standing at P3. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 round 37 by E. #213 (T-150): ubuntu ✓, windows running, InfoConsensus.Rd codoc mismatch noted in task. #216 (T-204): ubuntu ✓, windows running. **ASAN**: all PRs fail gcc-ASAN (rlang≤1.1.7 uses PREXPR removed in R-devel). Workaround added to ASan.yml (pre-install r-lib/rlang from GitHub). Open PRs: #213, #216, #210 (DRAFT). |
