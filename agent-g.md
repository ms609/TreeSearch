# Agent G — Progress Log

## Current Task
**Status:** IDLE

---

## Recently Completed

### T-205 — Fix flaky test-pp-random-tree.R (2026-03-24)
Root cause: MWC RNG static state in `build_postorder.h` not seeded by
`set.seed()`. With `stringency = 0.005`, the five-tip test had ~1%
false-positive rate per CI run.

Fix: widened binomial bounds (`stringency` 0.005→1e-6) and increased
`nTrees` (6000→12000, 12000→24000) across all tests. False-positive
rate now ~0.0002% per run. GHA pass: run 23501977394.

### S-PROF round 4 — Performance profiling (2026-03-24)

Post-pipeline-overhaul profiling. Score quality improved (Zhu2013 639–644
vs 648–666). NNI-perturb ~23% of cycle time (new), TBR polish 31%→1%.
No core speed regression. No new optimization tasks filed.

### T-203 — Simulated annealing for large trees (2026-03-24)
**Branch:** `feature/anneal` → merged directly to `cpp-search`

### T-179 — Large-tree strategy preset (2026-03-24)
Tuned via systematic benchmarking on mbank_X30754 (180t, 418p).

### Earlier completions
See `completed-tasks.md` for full history.
