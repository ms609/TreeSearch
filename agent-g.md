# Agent G — Progress Log

## Current Task
- **Task:** IDLE
- **Status:** Completed T-289f Stage 5 analysis

## Recently Completed

### T-289f Stage 5 — Prune-Reinsert NNI vs TBR Polish (2026-03-29)
Hamilton HPC benchmark (SLURM 16622421, 7h runtime). 5 large-tree datasets
(131-206t), 20 seeds, 60s/120s budgets, EW scoring. Three configs: baseline,
pr_nni (NNI polish), pr_tbr (TBR polish).

**Results:** pr_nni wins 7/10 conditions by expected-best. Huge benefit on
project3701 (146t, -178 median at 60s). Modest benefits at 173-180t. Slight
regression at 206t. pr_tbr harmful (1/9 wins; total starvation at 206t/60s).

**Decision:** Not enabled in large preset - benefit is dataset-dependent and
reverses at >=206t. Available via SearchControl(pruneReinsertCycles=5,
pruneReinsertNni=TRUE). strategies.md updated.

### S-COORD Round 45 (2026-03-28)
PRs #237 (T-279) and #238 (T-245) merged; rows deleted from to-do.md.

### S-RED Focus 30-31 (2026-03-28)
ts_drift.cpp (T-279): correct. ts_fitch.h/ts_tbr.cpp (T-245): correct.
ts_prune_reinsert.h/.cpp: G-006 filed (nni_search lacks ConstraintData*).

### T-290c wagnerStarts Benchmark (2026-03-28)
wagnerStarts=1 vs 3 under Brazeau scoring, 2 datasets (86-91t).
Current preset assignments confirmed correct.

---

## Earlier completions
See completed-tasks.md for full history.
