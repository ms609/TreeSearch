# Agent C Progress Log

## Current Task: T-214
**Status:** PARKED (C, GHA 23536512228)
**Description:** [Bug] Multi-split constraints not enforced during TBR search.
**Started:** 2026-03-25

### Progress
- Claimed task
- Reproduced bug: 63/50 seeds violated constraints on 10-tip trees
- Traced violation to TBR rerooting in `tbr_search()` during ratchet perturbation
- Root cause: `classify_clip_constraints()` marks clips as UNCONSTRAINED when
  they contain ALL tips from one side of a constraint split. But TBR rerooting
  at an edge between constraint tips and extras puts them on opposite sides of
  the attachment edge, destroying the split.
- Implemented two-part fix:
  1. Post-hoc `map_constraint_nodes()` after every accepted TBR/drift move;
     reject moves that introduce violations (safety net)
  2. FORBIDDEN clip zone for clips where both clip and rest straddle a split
     (early rejection optimization)
- Also fixed broken `.ts_driven_search_raw` → `ts_driven_search` callers
  (from af7601b refactor)
- Added `test-ts-constraint-multi.R` (Tier 2, 806 assertions): 10/12/15-tip
  trees, 2-3 constraint splits, EW + IW
- Committed on cpp-search (62658709d), GHA dispatched
