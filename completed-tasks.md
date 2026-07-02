# Closed Tasks — Decisions Worth Not Re-Litigating

This is **not** a full archive of every completed task. Routine fixes live in
git history and merged PRs; do not duplicate them here. This file keeps only
the closures whose *reasoning* is not recoverable from a commit: **not-a-bug
determinations, superseded/ruled-out designs, and negative experimental
results** — the things an agent would otherwise waste budget re-investigating.

**How to consult:** before investigating a recurring symptom or reopening a
closed `T-nnn`/`<Letter>-nnn`, **`grep` this file** (by ID or keyword) for a
prior closure. Do **not** `Read` it whole. Old task IDs referenced here remain
valid (see AGENTS.md) and need not be renamed.

**How to add a row:** only when a task closes *without a routine fix* — a
negative result, a "not a bug", or a superseded decision. One row, terminal
decision + a pointer to the write-up. Routine fixes get a one-line row in
`to-do.md`'s removal commit, not an entry here.

---

## Not a bug / scoring-method confounds

| ID | Topic | Decision |
|----|-------|----------|
| T-242 | Agnarsson2004 IW "~2% hit rate" | **Display bug only.** `ThreadSafePool::extract_into()` reset `hits_to_best` to distinct-topology count, not replicate hits. Search algorithm unaffected; real hit rate ~60–67%. Fixed `bc19667f2`; regression test in `test-ts-parallel.R`. |
| T-247 | XPIWE Vinther2008 score ≠ TNT | **Not a bug.** Discrepancy is entirely Brazeau three-pass vs standard-Fitch inapplicable handling. TreeSearch's tree (EW=79) is genuinely better under three-pass scoring. XPIWE uses `eff_k` in all paths — verified correct. |
| T-265 / T-249 / T-264 | Per-replicate "regression" vs TNT | **Scoring-method confound, not a regression.** T-249/T-264 compared Brazeau-scored TreeSearch to EW-scored TNT (apparent gap +17.8 steps; real EW-vs-EW gap +2.2, 5/11 datasets at 0). **Future TNT comparisons MUST use `fitch_mode()` for apples-to-apples.** |
| T-211 | Stale `final_` in temper candidate scoring | **Not worth fixing.** Conservative-only: stale `final_` biases Boltzmann screening but `temper_full_rescore` gates every accepted move. Fix cost (per-candidate rescore or full save/restore) exceeds negligible SA benefit. |
| T-325 / T-326 | MPT-enum collapse-dedup "30 trees" red-team finding (2026-07-02) | **Already fixed before filing — stale build, not stale source.** Root defect (enum pool deduped on full-key `add()` instead of collapsed-key `add_collapsed()`) was fixed 8 days earlier by `0daea13f` (2026-06-25): R-level `collapse=TRUE` default (`ts_collapse_pool`) + matching C++ enum change. The finder cited `ts_tbr.cpp:2345,2377`, which are that commit's **pre-fix** line numbers — the trace read a version of the file predating the fix even though HEAD had it. Verified by building both revisions fresh: `0daea13f^` reproduces 30 trees/`n_topologies=30` exactly as reported; current HEAD returns 1. The decisive fix layer is the R-level `collapse=TRUE` default, not the C++ enum change alone (confirmed by disabling the latter via `TS_ENUM_RESOLVED=1` — still 1). T-326's regression ask is already covered by `test-MaximizeParsimony-features.R:402-421` (added in the same commit; re-ran green). Residual `test-ts-collapsed.R:159-179` weak assertion is now redundant-weak, not dangerous-weak — not worth a dedicated task. Full writeup: `dev/red-team/log.md` area 11 post-hoc correction. |

## Superseded / ruled-out designs

| ID | Topic | Decision |
|----|-------|----------|
| T-183 | Pool-seeded Wagner / consensus backbone | **Superseded** by `consensusConstrain` (ts_driven.cpp), which constrains the whole replicate pipeline, not just Wagner. Marginal starting-tree value given the NNI→TBR pipeline. |
| T-198–201 | Boltzmann parallel tempering | **Ruled out** by T-199: 0% cold↔warm swap acceptance across all datasets. PCSA component salvaged as T-207/PR #227. See pt-evaluation expertise note. |
| T-185 | IQ-TREE acceleration ideas | Stochastic NNI-perturbation worth trying (→T-186, implemented). **Batch NNI not worthwhile** — see batch-nni expertise note. |

## Search-tuning experiments — settled, don't re-run

Each row records a benchmark whose conclusion fixed a default or killed an idea.
Detailed data is in the named `dev/benchmarks/` write-up; re-running wastes
Hamilton/GHA budget unless the underlying kernel has changed.

| ID | Experiment | Conclusion |
|----|-----------|------------|
| T-254 | Drift cycles (0 vs 2) | Drift gives **zero** score/MPT/diversity benefit, costs 10–22% of reps → `driftCycles=0` in default+thorough (T-255). `drift_mpt_analysis.md`. |
| T-256 | Sectorial intensity | Doubling/tripling xss+rss rounds → no score gain. Current `xss=3, rss=1` sufficient. |
| T-259 | Ratchet cycle count | Reducing 12→8/6/4 is mixed-to-worse; default **12 justified** (3-seed, directional). |
| T-274 | NNI-perturb cycles (thorough) | 59–69% overhead, ≤0.1-step benefit → `nniPerturbCycles=0` in thorough. `bench_t274_nni_perturb.R`. |
| T-248 | SA phase tuning (large) | `annealCycles=3` no significant gain over AC=1 (p>0.5) → **AC=1** in large preset (~6% faster). |
| F-029 (T-269) | Outer-cycle count | Higher `outerCycles` cuts replicate throughput with no score gain → **`outerCycles=2` optimal**. |
| PA-002 | XSS↔TBR cycling | Benefit scales with **tree size, not scoring mode**. ≤88t: pure overhead. 180t: −6.8 to −9.8 EW steps. No IW-specific treatment needed. `expt_tbr_xss_v2_results.rds`. |
| PA-003 | Targeted post-clip sector search | **NET HARMFUL** — local sector refinement after each TBR move steers into worse basins (+17 to +34 steps at 180t). Confirms XSS-as-a-separate-phase-after-TBR is correct. **Do not implement.** |
| PA-001 → F-030 | TBR clip ordering | PA-001's "tips-first falsified" was an **artifact** (clip_order reached only ~10% of TBR calls). F-030 with full propagation: `TIPS_FIRST` gives +8–13% throughput on 75–88t thorough preset. Default unchanged (`clipOrder=0`); PR #239. |
| T-289 / T-289f | Prune-reinsert polish (large) | TBR polish **catastrophic** at ≥206t (0 reps); NNI polish helps 131–180t → `pruneReinsertCycles=5, pruneReinsertNni=TRUE` in large preset, PR disabled elsewhere. `t289f_pr_nni_polish.csv`. |
| G-001 (T-290) | Brazeau-track phase profiling | Wagner is 3.6–5.2× costlier under Brazeau than Fitch; Fitch-tuned presets remain appropriate. `wagnerStarts=3` justified (better starting topology dominates when TBR convergence > budget). |
| F-006 (T-253) | Gap characterization | **ntax is the dominant predictor** of TNT gap (ρ≈0.63); nchar matters only >2000. `t253_gap_characterization.md`. |
| F-004 (T-252) | MorphoBank baseline | ≤35t converge at 30s; 66–135t still improving at 120s; project4284 (4062t) can't finish 1 replicate. CSVs in `dev/benchmarks/`. |
| T-251 | TNT trajectory analysis | Drift 30–170× less efficient than the next-worst phase; TNT spends ~67% in sectorial search vs TS's single pass. `tnt_trajectory_analysis.md`. |
| T-250 | TNT Fitch kernel disassembly | TNT is 32-bit scalar, no SIMD; TreeSearch has ~4× kernel throughput. TNT's convergence edge is **strategic, not implementation**. `tnt_disassembly_analysis.md`. |
| T-260 | VTune TBR overhead | Non-scoring overhead = 37.8% of TBR (StateSnapshot 14.6%, reset_states 9.1%) → T-261/262/263 (done). `vtune_tbr_analysis.md`. |

<!--
Full per-task fix history (everything closed with a routine commit/PR) was
purged 2026-06-16 — it lived only to duplicate git. Recover any specific row
from `git log` / merged PRs. See AGENTS.md "On task completion" for what now
warrants a row here.
-->
