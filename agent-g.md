# Agent G — Progress Log

## Current Task
- **Task:** T-208 — Fix `random_topology_tree` ignoring constraints
- **Status:** IN PROGRESS
- **Branch:** `feature/fix-random-tree-constraint`
- **Worktree:** TS-FixRandCons

## Parked Task
- **Task:** T-182 — Adaptive ratchet perturbation probability taper
- **Status:** PARKED (waiting on GHA run 23505912119)
- **Branch:** `feature/adaptive-ratchet`
- **Worktree:** TS-AdaptRatch

### T-208 Bug Description
When `adaptiveStart=true` (thorough preset) AND constraints active, the bandit
can select RANDOM_TREE → `random_topology_tree()` doesn't take ConstraintData →
starting tree violates constraint → TBR blocks all constraint-relevant moves
(`cn=-1`). Could return constraint-violating trees.

Fix: fall back to WAGNER_RANDOM when `strategy == RANDOM_TREE && cd && cd->active`
in `run_single_replicate()` (ts_driven.cpp:111).

---

## Recently Completed

### T-182 Implementation Summary
Cross-replicate ratchet perturbation probability tapering. When
ratchetTaper=TRUE, perturb_prob is scaled by `max(floor, 1 - strength * hitRate)`
each replicate. Default floor=0.5, strength=0.6. Serial path only.

Enabled in thorough and large presets. Neutral on 75-88 tip datasets
(median scores identical with/without taper on Zhu2013 thorough preset,
20s budget). Expected benefit at 120+ tips where ratchet cycles are expensive.

Files changed: ts_driven.h, ts_driven.cpp, ts_rcpp.cpp, TreeSearch-init.c
(64→65 args), SearchControl.R, MaximizeParsimony.R, RcppExports.*, vignette.

4 new tests; 159 driven + 40 ratchet + 29 SearchControl + 274 parallel pass.

### Also fixed
- `gha-check-pending.sh` bug: `|| true` masked exit codes from `gha-poll.sh`,
  making all runs appear PASS. Fixed with `&& EXIT=0 || EXIT=$?` pattern.

### T-205 — Fix flaky test-pp-random-tree.R (2026-03-24)
Root cause: MWC RNG static state in `build_postorder.h` not seeded by
`set.seed()`. With `stringency = 0.005`, the five-tip test had ~1%
false-positive rate per CI run.

### S-PROF round 4 — Performance profiling (2026-03-24)
Post-pipeline-overhaul profiling. Score quality improved.

### T-203 — Simulated annealing for large trees (2026-03-24)
### T-179 — Large-tree strategy preset (2026-03-24)

### Earlier completions
See `completed-tasks.md` for full history.
EOF 2>&1
