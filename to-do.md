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
| T-150 | P2 | PARKED (F, GHA 23644644554) | — | **CID-optimal consensus tree search** | PR #213 open to cpp-search. ASAN run 23643030700 failed: `rlang` compilation failure (infrastructure, not our code). Re-dispatched GHA 23644644554. |
| T-204 | P2 | PARKED (F, GHA 23644617599) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | On `feature/native-search`. GHA 23642888974 failed: WORDLIST fix was in eb21c588 (after merge commit fadc9d7e that triggered run). Re-dispatched GHA 23644617599. |


### Bugs

(no open bugs)


### Shiny App

(no open tasks)

### API / UX

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-276 | P3 | OPEN | — | **Print convergence summary to console after `MaximizeParsimony()`.** The Shiny app displays a convergence evaluation report; same info should print to the R console when running headlessly. Note: integrates with B's Chao-estimate work (T-204 native-search) once that merges. | u.571. Print n_replicates, best_score, n_MPTs, last_improved_rep, time elapsed, and convergence indicator (consensus_stable / plateau_stop / timed_out / perturb_stop). |
| T-277 | P3 | PARKED (B, GHA 23644927459) | — | **ScoreSpectrum(): Chao1-style landscape coverage estimator** | On `feature/score-spectrum`. Exports `ScoreSpectrum()`: accepts `multiPhylo` (with new `replicate_scores` attribute) or raw numeric vector; returns Good-Turing coverage + Chao1 richness. C++ side: `DrivenResult::replicate_scores` vector (serial + parallel paths). Shiny: coverage note appended to confidence text. 8 Tier-1 tests. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-274 | P2 | OPEN | — | **NNI-perturb cycle count at thorough-preset scale (65–88 tips).** S-PROF round 6 found NNI-perturb = 34.3% of Zhu2013 search time with only 14% hit rate (1-step improvement). Ratchet is 4–8 steps/call at comparable cost. Benchmark `nniPerturbCycles=0` vs `5` at 65–88 tips using time-adjusted expected-best metric. If no-NNI wins → reduce in thorough preset; if NNI helps at 88+ → set conditional on dataset size. | Zhu2013/Dikow2009/Giles2015, 3–5 seeds, 30s/60s budgets. Use `expected_best()` in profiling.md for comparison. |
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
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 3 by F (Ratchet/perturbation/sector/prune-reinsert). ts_prune_reinsert.cpp fully reviewed (511 lines, new). T-275 guard correct. final_[tip] init safe via load_tip_states. Sector from_above_for_sector correct. XSS/CSS adaptive early-exit correct. T-273 fix correct. No new bugs found. Next: focus on ts_driven.cpp (cross-replicate constraint tightening, outer loop) — not reviewed since T-189. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 33 by F. T-273 completed. T-275 filed and completed by B (prune-reinsert guard). T-277 filed (ScoreSpectrum, ASSIGNED B). Unblocked OPEN: T-245, T-269, T-274, T-275→done, T-276, T-253, T-277 → 6 specific OPEN → standing at P3. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 by F (round 33). #216 (T-204 native-search): GHA 23642888974 failed (WORDLIST fix eb21c588 pushed AFTER merge triggered run). Re-dispatched 23644617599. #213 (T-150 cid-consensus): GHA 23643030700 failed: rlang ASAN infrastructure issue. Re-dispatched 23644644554. Open PRs: #213, #216, #210 (cpp-search→main DRAFT). |
