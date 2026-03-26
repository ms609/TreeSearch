# Agent E — Progress Log

## Current Task
- **Status:** T-265 — RESOLVED. Scoring method confound.
- **Root cause:** T-249/T-264 compared Brazeau-scored TreeSearch to EW-scored TNT. Apparent gaps (+19 to +59) are mostly scoring method difference. Actual EW-vs-EW gap: mean +2.2 steps across 11 datasets. Wilson2003 Brazeau-optimal IS 879; EW-optimal is 860 (0 gap). R2-equiv / R2-modern / auto preset all find identical scores — no preset regression.
- **Hamilton job 16597206 running** — will provide fuller confirmation at 120s with EW scoring.
- **Stale library finding:** Old .agent-E library had broken parameter passing, causing T-249 early termination artifact. Fresh rebuild fixed this (53 reps vs 8 at 30s budget).

### S-COORD Round 25 — DONE
Updated coordination.md, task queue, S-PR. T-264 GHA ARM64 passed, Windows in progress.

### S-RED Focus 8 — DONE
### Previous: T-261+T-262 — PR #232 MERGED
### Previous: T-255 — Reduce drift — DONE
### Previous: T-260 — VTune profiling — DONE
