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
| T-148 | P2 | ASSIGNED (B) | — | **Red-team: ParsSim vectorized convolution** | S-RED focus: log-space convolution correctness and active-range bounds in `R/ParsSim.R`. Use `.positai/expertise/red-team.md` methodology. |



### Extended Implied Weighting (XPIWE) — feature/xpiwe, worktree ../TS-Xpiwe

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-156 | P1 | OPEN | T-157 | **C++ core**: `ScoringMode::XPIWE`; `DataSet::eff_k` (per-pattern `k + global_min[p]`); unify `compute_iw` / `precompute_iw_delta` to use `eff_k[p]` | See plan `.positai/plans/2026-03-20-0923-*.md` |
| T-157 | P1 | OPEN | T-158 | **Rcpp bridge**: `bool xpiwe = false` to `make_dataset`, `ts_fitch_score`, `ts_driven_search`, `ts_resample_search`, `ts_successive_approx`; update `TreeSearch-init.c` (run `check_init.R`) | Append-only to ts_rcpp.cpp sigs |
| T-158 | P1 | OPEN | T-159 | **R API**: `extended_iw = TRUE` in `MaximizeParsimony()`, `TreeLength()` (all S3), `SearchControl()`, `Resample()`, `SuccessiveApproximations()`; silently ignore when EW / profile | — |
| T-159 | P1 | OPEN | T-160 | **Tests**: `test-ts-xpiwe.R` (Tier 2) — formula unit, 5-taxon hand-computed, `m=0` IW equivalence, end-to-end; IW `eff_k` regression in existing IW test | — |
| T-160 | P2 | OPEN | — | **Docs + NEWS**: `@param extended_iw` in Rd files; para in `vignettes/profile-scores.Rmd`; NEWS breaking-change entry; `check_man()` + `spell_check_package()` clean | — |

### Alternative Inapplicable-Handling Algorithms — Phase 3: Integration & polish

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|



### Shiny Search UX Improvements

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|




### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-19 by B (focus 7: Shiny module wiring). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-19 by A (round 3). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-19 by Agent A (round 8). |
