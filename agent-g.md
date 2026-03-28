# Agent G — Progress Log

## Current State

- **Status:** ACTIVE — T-290c wagnerStarts confirmatory benchmark
- **Date:** 2026-03-28 ~18:45 GMT

## Recent History (2026-03-28)

### G-001 — T-290: Brazeau-track phase profiling + wagnerStarts analysis (COMPLETE)

Phase profiling {Brazeau,Fitch}×{EW,IW10}×{default,thorough} on 6 datasets
(23–173t), 30s, 3 seeds. Key findings:
- Wagner 3.6–5.2× more expensive under Brazeau; ratchet 1.1–1.3×; rep rate
  ~97% of Fitch.
- wagnerStarts=3 benefit is topology quality, not cycle count — datasets
  where thorough > default complete 0 reps in 30s.
- Conclusion: Fitch-tuned presets are appropriate for Brazeau scoring.
- strategies.md + AGENTS.md updated. Results in TS-TNT-bench t290_results/.

### G-002 — S-COORD round 45 (COMPLETE, 2026-03-28 ~18:15 GMT)

- PRs #237 (T-279 drift fix) + #238 (T-245 TBR batching) merged; rows
  deleted from to-do.md; entries added to completed-tasks.md.
- T-289f Stage 5: SLURM 16622224 running (~3 min elapsed at time of check).
- Open PRs: #213 (T-150), #216 (T-204), #210 (DRAFT).

### G-003 — S-RED focus 30 (COMPLETE, 2026-03-28 ~18:45 GMT)

ts_drift.cpp (T-279 post-fix review):
- IW reject: restore → rescore → update_constraint ✓
- EW reject: "score already set above" → update_constraint ✓
- Pre-existing edge: re-apply failure in EW accept path skips update_constraint
  (unreachable in practice; not filed)
- No bugs.

ts_fitch.h + ts_fitch.cpp + ts_tbr.cpp (T-245 post-merge review):
- x4 EW computation: identical to scalar (any_hit_reduce + ~a & mask) ✓
- x4 NA computation: identical to scalar (from1/clip_ha/ba[b]/mask) ✓
- Early-exit: all-4 >= cutoff, bitwise AND — conservative and correct ✓
- Partial batch (< 4): correctly falls back to scalar functions ✓
- State updates (best_candidate, best_above, best_below, best_reroot_*):
  identical to scalar path ✓
- No bugs.

### G-004 — T-290c: wagnerStarts confirmatory benchmark — NEXT

Quick benchmark (30s, 5 seeds) comparing wagnerStarts 1 vs 3 on the
thorough preset to confirm analytical prediction from T-290b.

## Earlier History

See `completed-tasks.md` for full history. Prior work includes:
- T-208: Fix random_topology_tree ignoring constraints → PR #219 merged
- T-182: Adaptive ratchet taper (superseded by T-248)
- T-205: Fix flaky test-pp-random-tree.R
- S-PROF rounds 4–6, T-203, T-179 (large-tree preset)
AGENTEOF 2>&1
