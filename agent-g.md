# Agent G — Progress Log

## Current Task
IDLE (2026-03-24)

### Completed: T-179 (2026-03-24)
Tuned large-tree strategy preset via systematic benchmarking on
mbank_X30754 (180 taxa, 418 patterns).

**Phase timing analysis:** At 180 tips, NNI-perturbation dominated the
budget (45.5% of 60s at 5.5s/cycle). The original large preset never
reached drift. Thorough preset's outerCycles=2 interleaving helped
distribute effort but added overhead.

**Design decisions (3 rounds of benchmarking, 5 seeds each):**
- ratchetCycles: 20 → 12 (cheaper per-cycle cost dominates)
- driftCycles: 12 → 4 (moderate drift still valuable)
- nniPerturbCycles: 5 → 0 (too expensive; ratchet more escapes/second)
- outerCycles: 2 → 1 (interleaving overhead not worthwhile)
- wagnerStarts: 3 → 1 + wagnerBias=1 (biased start saves ~2.6s)
- tbrMaxHits: 3 → 1 (faster TBR passes)
- xssRounds: 5 → 3, rssRounds: 3 → 2, cssRounds: 2 → 1
- adaptiveStart: tested, regressed (not enough reps for bandit)

**Validation (median scores, mbank_X30754):**
| Budget | large | thorough | Delta |
|--------|-------|----------|-------|
| 30s | 1276 | 1283 | +7 |
| 60s | 1255 | 1259 | +4 |
| 120s | 1250 | 1250 | tied (2 vs 0-1 reps) |

Commit fab1e52c.

### Completed: T-177 (2026-03-23)
- Verified mid-TBR/SPR/NNI timeout callback implementation
- `check_timeout` polled every ~n_tip clips in tbr_search, spr_search, nni_search
- 282 targeted tests pass (0 fail); build clean

### Completed: S-COORD round 9 (2026-03-23)
- Updated coordination.md: agent status, task pipeline, observations

### Completed: T-180 (2026-03-23)
- bench_warmstart.R: warm-start benchmark infrastructure

### Completed: T-181 (2026-03-23)
- Added mbank_X30754 as large-tree benchmark tier

### Completed: T-164 (2026-03-23)
- Pool stats wired to Shiny (topology count, trajectory)

### Completed: T-178: NNI warmup in driven pipeline (2026-03-23)
- NNI always-on, SPR auto-skipped, constraint guard

### Completed: T-186: Stochastic NNI-perturbation (2026-03-23)
- New ts_nni_perturb.h/cpp, integrated in driven pipeline

### Completed: T-156–T-162 XPIWE feature (2026-03-23)
- All 7 tasks on feature/xpiwe branch

### Completed: T-185: IQ-TREE review (2026-03-23)
- Top idea: stochastic NNI-perturbation (implemented as T-186)
