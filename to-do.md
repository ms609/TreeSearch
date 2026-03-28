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
| T-150 | P2 | PR #213 (F) | — | **CID-optimal consensus tree search** | PR #213. Vignette fix (TreeTools::Consensus) commit f8bfee49. GHA 23650002703. |
| T-204 | P2 | PR #216 (F) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | GHA 23649607006 PASSED. Ready for merge. |


### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-279 | P2 | PR #237 (F) | — | **Constrained drift: stale `cd->constraint_node` after RFD rejection + redundant full_rescore.** In `ts_drift.cpp::drift_phase()`, when a constrained suboptimal move passes the constraint check (so `cd` is updated), then is rejected by the RFD test, `update_constraint()` is not called on the restored topology. Bug in both IW reject (lines 647–649) and EW reject (line 689). Same class as T-278. Also eliminates a redundant `drift_full_rescore()` call in the EW reject path (return value at line 657 was discarded, causing a second rescore at line 689). Only affects constrained drift (`driftCycles > 0`); drift disabled in all presets since T-255. | Found by S-RED focus 8 (F, 2026-03-27). PR TBD. feature/drift-constraint-fix commit e85ec84f. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### Alternative Homologies (Goloboff 2026) — `feature/alt-homology` / `TS-AltHom`

Ref: Goloboff (2026) *Cladistics* doi:10.1111/cla.70033.
Plan: `.positai/plans/2026-03-27-1415-implement-goloboff-2026-alternative-homologies-with-step-matrix-recoding.md`

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-280 | P3 | OPEN | — | **AltHom Phase 1: `AlternativeHomology` S3 class & core recoding (MVP).** Create `R/AlternativeHomology.R` (constructor, validation, print), `R/recode_alt_homology.R` (correspondence enumeration, morphotype states, cost matrix, tip assignment). Wire into `TreeLength()` for scoring on a fixed tree. Reproduce paper's Definition 1 cost matrix + Table 1 as tests. | WORKTREE (TS-AltHom). Invertible, no external constraints, two part-types only. |
| T-281 | P3 | OPEN | T-280 | **AltHom Phase 2: Constraints & options.** Non-invertible (`>`), adjacent (`>>`), restricted homology (`!`), configurable part transformation costs, adjacent-loss merging (`<`). Reproduce Definitions 2–3 and their cost matrices. | WORKTREE (TS-AltHom). |
| T-282 | P3 | OPEN | T-280 | **AltHom Phase 3: Wire into `MaximizeParsimony()` search pipeline.** Accept `AlternativeHomology` in `hierarchy` param, prepare xformArgs, end-to-end search. Also wire `Resample()` and `SuccessiveApproximations()`. | WORKTREE (TS-AltHom). |
| T-283 | P3 | OPEN | T-280 | **AltHom Phase 4: External inapplicability.** An external character can make individual characters, parts, or entire part sets inapplicable. Expand state enumeration for externally-disabled states. | WORKTREE (TS-AltHom). |
| T-284 | P3 | OPEN | T-280 | **AltHom Phase 5: Combination pruning.** Implement `xlinks&` (pairwise compatibility), `xlinks!` (observed-state-only), `xlinks@` (uninformative-state restriction) to reduce supercharacter state count. Verify same optimal trees as unpruned. | WORKTREE (TS-AltHom). |
| T-285 | P3 | OPEN | T-280 | **AltHom Phase 6: Implied weighting support.** Compute combined minimum steps across all valid alignments (not sum of per-char minima). Required for correct IW homoplasy counts. | WORKTREE (TS-AltHom). |
| T-286 | P3 | OPEN | T-280 | **AltHom Phase 7: Mixed `AlternativeHomology` + `CharacterHierarchy`.** Support datasets with both simple hierarchy blocks and alternative homology blocks in one analysis. | WORKTREE (TS-AltHom). |
| T-287 | P3 | OPEN | T-284 | **AltHom Phase 8: Static alignment fallback.** For datasets where supercharacter exceeds practical state limit, generate alternative static datasets (one per alignment) and search each. | WORKTREE (TS-AltHom). |
| T-288 | P3 | OPEN | T-282 | **AltHom Phase 9: Documentation & vignette.** `vignettes/alternative-homologies.Rmd`, roxygen docs for all new exports, `inst/REFERENCES.bib` entry. | WORKTREE (TS-AltHom). |

### TNT Comparison & Strategy Learning

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-289 | P2 | PARKED (E, SLURM 16608629) | — | **Prune-reinsert benchmark (Hamilton).** Stage 1: sweep cycles (1/3/5) × drop fraction (5/10/20/30%) × RANDOM on 5 datasets (37–180t), 5 seeds, 30s. Stage 2: best configs + INSTABILITY selection + ratchet-replacement at 30s/60s. Scripts: `bench_prune_reinsert.R`, `t289_hamilton.sh`. Goal: decide whether to enable in presets and at what settings. | F's job 16607721 failed (install exit code 3). E's job 16608629 submitted ~15:59 GMT 2026-03-27; still in dep-install phase as of 16:22. Expected results: 18:30–19:30 GMT. |
| T-269 | P3 | PARKED (F, SLURM 16607719/16607720) | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. **Design note (u.005):** Full TNT-style interleaving IS architecturally supported now — setting `outerCycles = ratchetCycles` achieves one sector pass per ratchet cycle. T-257 first validated that any post-ratchet sectorial helps at all (merged). T-269 benchmarks the fine-grained variant to determine whether the per-cycle overhead is worth enabling in presets. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |
| T-290 | P3 | OPEN | — | **Brazeau-track baseline benchmark (Hamilton).** Run default + thorough presets on the fixed 20-matrix Brazeau sample (`MBANK_BRAZEAU_SAMPLE`), EW + IW(k=10), 30s, 5 seeds. Compare phase timing distributions and strategy rankings to Fitch baselines. Includes sample-size stability analysis (Kendall's W). Script: `dev/benchmarks/bench_brazeau_baseline.R` (in TS-TNT-bench). Goal: determine whether Brazeau scoring needs different strategy presets than Fitch. | Infrastructure added in TS-TNT-bench. |


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-291 | P3 | OPEN | — | **Update `bench_framework.R` `benchmark_run()` to new `ts_driven_search` interface.** The function currently passes flat keyword args but the C++ bridge now expects structured lists (`scoringConfig`, `searchControl`, `runtimeConfig`). Needs refactor to match the current Rcpp signature. Does not block Brazeau benchmarking (which uses `MaximizeParsimony` directly). | |





### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 27 by F (ts_rcpp.cpp, 2656L — Rcpp bridge). No bugs. All correctness invariants verified (sentinels, indexing, layout). F sweep now covers: ts_fitch, ts_simd, ts_fitch_na, ts_fitch_na_incr, ts_tree, ts_splits, ts_collapsed, ts_hsj, ts_sankoff, ts_rng, ts_rcpp, ts_parallel, ts_sector, ts_pool, ts_fuse, ts_driven, ts_drift, ts_tbr, ts_ratchet, ts_nni_perturb, ts_resample, ts_prune_reinsert, ts_simplify, ts_search, ts_data, ts_wagner, ts_constraint. **All modules reviewed ≥ once**. E sweep area 4: ts_rng, ts_parallel. Next: re-review post-merge additions (ts_strategy.h from T-190, ts_pcsa from T-207) or wait for new code changes. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 40 by F. PRs merged since last S-COORD: #233 (T-246 AVX2), #234 (T-257 post-ratchet sectorial), #235 (T-266 prune-reinsert), #236 (ScoreSpectrum). Still open: #213 (T-150), #216 (T-204), #237 (T-279). F-027 WORDLIST fix (config/warmup R 4.1 CI failures since ~12:50) committed ef83e8db → GHA 23656560997 PARKED. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 round 40 by F. Open PRs: #213 (T-150, GHA 23650002703 PASS), #216 (T-204, GHA 23649607006 PASS), #237 (T-279, GHA 23650290962 PASS), #210 (DRAFT cpp-search→main). F-027 WORDLIST fix parked GHA 23656560997. |
