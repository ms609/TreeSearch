# Agent G — Progress Log

## Current Task
**Status:** ACTIVE — T-205 (Flaky test-pp-random-tree.R on Windows)
**Started:** 2026-03-24

Root cause: The C-level MWC RNG (static variables `z`, `w` in
`build_postorder.h`) has no reset function. `set.seed(0)` only seeds R's
RNG. The five-tip test uses `stringency = 0.005` (1% false-positive rate),
which is too high for a Tier 1 CRAN test. After ~100 CI runs, failure is
near-certain.

Plan: Increase nTrees and decrease stringency to make tests robust.

---

## Recently Completed

### S-PROF round 4 — Performance profiling (2026-03-24)

Post-pipeline-overhaul profiling round. Major changes since round 3:
ratchet tuning, biased Wagner, adaptive bandit, outer cycles, NNI-perturb,
XPIWE, simulated annealing, MPT enumeration fix.

**Findings:**
- Score quality improved: Zhu2013 639–644 (was 648–666), Dikow2009 1611–1613
- Phase distribution shifted: NNI-perturb ~23% (new), TBR polish 31%→1%
  (NNI warmup gives near-optimal start), ratchet+drift still 60–70%
- Adaptive start (bandit) and adaptive level both confirmed working
- Outer cycle reset productive (Zhu2013 rep 1 found 5 successive improvements)
- No core speed regression (bare TBR timings consistent with machine load)
- No new optimization tasks filed

Updated `.positai/expertise/profiling.md` with full results.

### T-203 — Simulated annealing for large trees (2026-03-24)
**Branch:** `feature/anneal` → merged directly to `cpp-search`

### T-179 — Large-tree strategy preset (2026-03-24)
Tuned via systematic benchmarking on mbank_X30754 (180t, 418p).

### Earlier completions
See `completed-tasks.md` for full history.
