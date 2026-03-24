# Agent G — Progress Log

## Current Task
**Status:** IDLE

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

Implemented linear-cooling simulated annealing using stochastic TBR with
Boltzmann acceptance.  Created `ts_temper.h/.cpp` with Layer 1
(`stochastic_tbr_phase`, ported from Agent C's T-198) and Layer 3
(`anneal_search`, new).  Wired into driven pipeline between drift and
final TBR polish.

R-level: `SearchControl(annealPhases, annealTStart, annealTEnd,
annealMovesPerPhase)`.  `large` preset disables drift, enables 5 annealing
phases T=20→0.

Merge resolved conflict with `enumTimeFraction` (added to `cpp-search`
since branch point).  Updated `TreeSearch-init.c` arg count 63→64.
19 new tests, 152 driven + 219 MP feature tests pass.

### T-179 — Large-tree strategy preset (2026-03-24)
Tuned via systematic benchmarking on mbank_X30754 (180t, 418p).
Commit `fab1e52c` on `feature/parallel-temper`.

### Earlier completions
See `completed-tasks.md` for full history.
