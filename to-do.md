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
- Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
  completed, reset it to OPEN. Their effective priority is dynamic:
  - ≥6 OPEN specific tasks → standing tasks are P3
  - 3–5 OPEN specific tasks → standing tasks are P2
  - <3 OPEN specific tasks → standing tasks are P1

---

## Active Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|


| T-150 | P2 | WORKTREE (TS-CID-cons) | — | **CID-optimal consensus tree search** | Use SPR/TBR/NNI + Ratchet to find consensus trees minimising CID to input trees. Needs: collapse/resolve moves for non-binary trees, CID batch scorer, `CIDBootstrap` (resample input trees). See `briefing-cid-consensus.md`. PoC validated in TreeDist. |




### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-177 | P1 | ASSIGNED (Human+AI) | — | **Bug fix: mid-TBR/SPR timeout.** Pass `check_timeout` callback into `tbr_search()` and `spr_search()` so they can bail out mid-pass. At 180+ tips, a single TBR convergence takes ~13s; without mid-pass timeout, `maxSeconds` can overrun significantly. | **Implemented**, building and testing. Added `std::function<bool()>` param to `tbr_search`/`spr_search`; threaded through driven, ratchet, drift calls. 1762 Tier 2 tests pass. |
| T-178 | P1 | ASSIGNED (G) | — | **NNI warmup in driven pipeline.** Add `nni_search()` before TBR in `run_single_replicate()`. Always-on (NNI costs ~1.5s at 180 tips, negligible at small sizes). Disable `sprFirst` when NNI is active (SPR intermediate step is counterproductive at 180 tips — adds time, worse scores). Validated: NNI→TBR is 2× faster and ~100 steps better than TBR alone at 180t across 3 seeds. | `nni_search()` exists in `ts_search.cpp`. May need `nniFirst` param in `SearchControl()` or just always-on. Enhancement: NNI-optimize each of the `wagnerStarts` trees before picking the best (IQ-TREE generates 100 parsimony trees, NNI-optimizes top 20). |
| T-179 | P2 | OPEN | T-177, T-178 | **Large-tree strategy preset.** Add a preset for ≥120 tips (or ≥100 tips with ≥200 patterns): NNI→SPR→TBR escalation, scaled ratchet/drift cycles, tuned sector sizes. Benchmark on 180-taxon dataset. | May need `nniFirst` param in `SearchControl()`. |
| T-180 | P2 | OPEN | — | **Warm-start benchmark infrastructure.** Create benchmark mode that seeds search with a pre-computed local optimum (from a short prior search) to measure ratchet/drift escape effectiveness in isolation, separate from initial descent quality. | `start_ptr` already supported in C++. |
| T-181 | P2 | OPEN | — | **Add 180-taxon dataset to benchmark suite.** Copy `mbank_X30754` to `inst/benchmarks/`, add to `bench_framework.R` as a "large" tier with separate timing expectations. | 180t, 425c, 374 informative patterns, 40% missing, 20% inapplicable. |
| T-182 | P3 | OPEN | — | **Adaptive ratchet perturbation probability.** Start aggressive (~40%) and taper by hit rate as pool stabilizes. Extend `adaptive_level` infrastructure to also scale `perturb_prob`. Risk: premature convergence if tapering too fast. | Also consider IQ-TREE-style tree-size scaling: 50% at ≤50 tips, 5% at ≥400 tips. Two orthogonal axes: hit-rate tapering (within run) and size scaling (across datasets). |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** Use pool consensus as backbone constraint during Wagner construction for later replicates (after ≥N diverse trees). Partially randomised addition order. Concern: run independence — mitigate by only activating late. | Constraint infrastructure already exists (`consensus_constrain`). |
| T-184 | P2 | OPEN | — | **`maxTime` triggers Morphy delegation.** `maxTime` is a natural parameter name users try, but it activates the deprecated Morphy shim instead of mapping to `maxSeconds`. Fix: treat `maxTime` as an alias for `maxSeconds` with a deprecation warning (instead of delegating to `Morphy()`). | In the old-param detection block of `MaximizeParsimony()`. |
| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop search after `nTip × K` unsuccessful perturbations (no score improvement). Dataset-adaptive convergence criterion complementary to time-based `maxSeconds` and consensus-stability stopping. IQ-TREE uses `nSeq × 100`. | From T-185 IQ-TREE review. Needs counter in `driven_search()` loop tracking perturbations since last improvement. |


### Extended Implied Weighting (XPIWE) — feature/xpiwe, worktree ../TS-Xpiwe

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-156 | P1 | OPEN | T-157 | **C++ core**: `ScoringMode::XPIWE`; `DataSet::eff_k` (per-pattern `k + global_min[p]`); unify `compute_iw` / `precompute_iw_delta` to use `eff_k[p]` | See plan `.positai/plans/2026-03-20-0923-*.md` |
| T-157 | P1 | OPEN | T-158 | **Rcpp bridge**: `bool xpiwe = false` to `make_dataset`, `ts_fitch_score`, `ts_driven_search`, `ts_resample_search`, `ts_successive_approx`; update `TreeSearch-init.c` (run `check_init.R`) | Append-only to ts_rcpp.cpp sigs |
| T-158 | P1 | OPEN | T-159 | **R API**: `extended_iw = TRUE` in `MaximizeParsimony()`, `TreeLength()` (all S3), `SearchControl()`, `Resample()`, `SuccessiveApproximations()`; silently ignore when EW / profile | — |
| T-159 | P1 | OPEN | T-160 | **Tests**: `test-ts-xpiwe.R` (Tier 2) — formula unit, 5-taxon hand-computed, `m=0` IW equivalence, end-to-end; IW `eff_k` regression in existing IW test | — |
| T-160 | P2 | OPEN | — | **Docs + NEWS**: `@param extended_iw` in Rd files; para in `vignettes/profile-scores.Rmd`; NEWS breaking-change entry; `check_man()` + `spell_check_package()` clean | — |
| T-161 | P2 | OPEN | T-158 | **Shiny GUI**: Add XPIWE option; replace default "Implied" weighting with "Implied (extended)"; keep "Implied" in list; both share the concavity slider | TS-Xpiwe worktree |
| T-162 | P3 | OPEN | T-160 | **Docs**: Add Goloboff 2014 citation to suggested references when xpiwe is selected | TS-Xpiwe worktree |

### Alternative Inapplicable-Handling Algorithms — Phase 3: Integration & polish

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|



### Shiny Search UX Improvements

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-163 | P2 | OPEN | — | **Search confidence: composite diagnostic** Replace `exp(-K)` in `SearchConfidenceText()` with tightened binomial bound `(1 - K/R)^R` (falling back to `exp(-K)` when K==R). Add pool topology count (T) and last-improvement replicate index to the diagnostic message. Report landscape ruggedness when K/R < 0.3 and flag limited independence when T==1. | See discussion in AGENTS.md "Future: search convergence" |
| T-164 | P2 | OPEN | T-163 | **Search confidence: wire pool stats** Ensure `ts_driven_search` return list includes (a) number of distinct topologies in pool at best score, (b) replicate index of last score improvement. Add to Shiny search module's result handling. | May need small C++ addition to return value |





### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-19 by B (focus 7: Shiny module wiring). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-19 by A (round 3). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-19 by Agent A (round 8). |
