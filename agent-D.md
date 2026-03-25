# Agent D Progress Log

## Current Task: S-RED focus 1
**Status:** COMPLETE
**Started:** 2026-03-25

### Findings
- **T-229 (P2 bug):** Found and fixed XFORM scoring bug in `fitch_score_ew()`.
  Missing `ScoringMode::XFORM` in EW branch caused non-hierarchy chars to be
  scored via IW path (k=0) instead of EW. MaximizeParsimony reported wrong scores
  (e.g. 3 instead of 7). 1-line fix. Committed to cpp-search.
- XPIWE formulas verified correct (Goloboff 2014 Extension 3).
- NaN guard for k=0 correct.
- Minor: obs_count=0 division risk (defensive, not filed).
- IW+hierarchy correctly blocked by R validation.
