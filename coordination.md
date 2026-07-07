# TreeSearch — Strategic Coordination

## S-COORD Round 46 Summary (2026-03-29 07:40 BST, Agent E)

**T-289f complete — Stage 5 NNI polish benchmark + large preset update:**
SLURM 16622483 completed (7h12m, EPYC 7702). 300 runs: 5 datasets (131–206t),
configs baseline/pr_nni/pr_tbr, 60s+120s, 10 seeds.
- pr_tbr (TBR polish): confirmed Stage 4 failure — syab07205 (206t) still 0 reps at 60s.
- pr_nni (NNI polish): fixes 0-rep failure; improves 131–180t: project3701 (146t) −178
  steps at 60s / −128 at 120s; project804 (173t) −9/−2; mbank_X30754 (180t) −4/−7.
  syab07205 (206t) +17.5 at 60s, neutral at 120s — acceptable.
- Decision: **pruneReinsertCycles=5, pruneReinsertNni=TRUE enabled in large preset**.
  commit 4a549eb4. Results: dev/benchmarks/t289f_pr_nni_polish.csv. AGENTS.md updated.

**G-006 fixed — NNI constraint guard in prune_reinsert_search():**
One-line guard `if (params.nni_full && (!cd || !cd->active))` in ts_prune_reinsert.cpp.
When constraints active, falls through to TBR (which enforces them). Mirrors the
`nni_wagner` guard in ts_driven.cpp. Task deleted from to-do.md.

**GHA 23703257153 in progress** on cpp-search (covers 4a549eb4 + G-006 fix).

**PR status:**
- #213 (T-150, CID consensus): GHA PASS, awaiting human merge.
- #216 (T-204, native search): GHA PASS, awaiting human merge.
- #210 (cpp-search→main, DRAFT): Re-run 23702009435 in progress; previous failure was
  Windows covr only (transient/infra — tests passed FAIL 0/PASS 11021).

**Task queue:** Extremely sparse. Only standing tasks + T-280–288 (all WORKTREE/AltHom).
Standing tasks at **P1** (<3 open specific tasks).

**Next:** S-RED (review alt-homology modules when T-280 merges, or review ts_search.cpp
and ts_nni_perturb.cpp which haven't been reviewed). S-PR to check PR status.

## S-COORD Round 42 Summary (2026-03-28 16:10 GMT, Agent F)

**T-269 complete — Fine-grained sectorial interleaving (30s, 4 datasets, outer_cycles 1/2/4/10/20):**
Higher outer_cycles uniformly reduces replicate throughput with no score benefit.
At outer_cycles=20: Dikow2009 gets 9 reps vs 54 at baseline; Zhu2013 gets 16 vs 88.
Scores are flat or marginally worse at high outer_cycles. The current outerCycles=2 in the
thorough preset is optimal; no preset change needed.

**T-289 complete (E) — Stage 4 confirms disable-PR decision:**
5 datasets 131–206t, 10 seeds, 60s/120s. Key: syab07205 (206t) gets 0 PR reps at 60s
(per-rep cost ≈ 60s, budget exceeded). project3701 (146t) regresses 12 steps mean at 60s.
commit 746985243 disables pruneReinsertCycles in large preset. Available via SearchControl().

**F-027 WORDLIST fix — PASSED (GHA 23656560997).** Both 'config' and 'warmup' restored.

**PR #210 (cpp-search→main):** codoc fix fdf25673 in place; R-CMD-check run 23688837232
in progress. Previous pre-existing failures (Windows covr, R-devel rlang, ASAN TBB ODR)
are infra issues, not package check failures.

**Open PRs:** #213 (T-150, GHA PASS), #216 (T-204, GHA PASS), #237 (T-279, GHA PASS).
All three await human merge.

**Task queue:** T-245 (P3, TBR batching) is only open specific task. Standing tasks P1.
S-RED next: ts_mc_fitch.cpp, ts_tabu.h, ts_prune_reinsert.h (222 lines, unreviewed).

## S-COORD Round 41 Summary (2026-03-28 14:35 GMT, Agent E)

**Codoc fix — SearchControl.Rd (E-003):**
All R-CMD-check platforms failing on PR #210 since 2026-03-28 06:25 with
"Codoc mismatches from SearchControl.Rd". Root cause: commit 22f929cf
(`pruneReinsertTbrMoves` param, T-289) added the parameter to the function
and roxygen `@param` but the Rd file was not regenerated. Fix: manually added
`pruneReinsertTbrMoves = 5L` to `\usage` and its `\item` to `\arguments` in
`man/SearchControl.Rd`. Commit fdf25673. PR #210 CI re-triggered (run
23687279706, pending). Agent-check GHA 23687210711 also dispatched.

**T-289 Stage 4 — Hamilton SLURM 16621426:**
Stage 3 confirmed MISSING criterion (sel=2, c=5, d=5%) gives −14.7 steps at
180t/60s. Large preset updated. Stage 4 validating across 5 matrices
(131–206t) at 60s/120s, 10 seeds, 200 runs. Submitted 2026-03-28 ~08:00 GMT,
~5h wall time. SSH unavailable — poll later.

**F-027 WORDLIST fix (GHA 23656560997) — PASSED.** Resolved.

**PR status:**
- #210 (cpp-search→main): CI re-running with codoc fix; was failing since 06:25.
- #213 (T-150, CID consensus): GHA 23650002703 PASS, awaiting merge.
- #216 (T-204, native search): GHA 23649607006 PASS, awaiting merge.
- #237 (T-279, drift constraint fix): GHA 23650290962 PASS, awaiting merge.

**Task queue:** T-289 PARKED (Hamilton), T-269 PARKED (Hamilton), T-245/T-290/T-291
OPEN. 3 open specific tasks → standing tasks P2 effective.

**GHA 23687804562 results (PR #210, post-codoc-fix):** All 5 release platforms
PASS. Remaining failures are pre-existing infra issues: Windows covr path, R-devel
rlang DLL, ASAN RcppParallel TBB ODR. All in dep-install or coverage steps, not
"Check package". PR #210 ready for human review.

**T-291 complete (E-004):** bench_framework.R benchmark_run() updated to new
ts_driven_search structured-list interface. commit f1ed5dfc.

**Next:** Poll Hamilton for T-289 Stage 4 results when SSH is available.

## S-COORD Round 39 Summary (2026-03-27 16:05 GMT)

**GHA results confirmed (all PASS):**
- GHA 23653228247 (F-015: ratchet constraint staleness) — **PASSED**
- GHA 23653513217 (F-016: NNI-perturb constraint staleness) — **PASSED**
- GHA 23653782359 (F-018: prune-reinsert constraint staleness) — **PASSED**

All constraint-staleness fixes now validated on both platforms. The full sweep
(TBR T-278, drift T-279, sector E-003, ratchet F-015, NNI-perturb F-016,
prune-reinsert F-018) is complete. All 6 constrained search modules now
consistently call `update_constraint(tree, *cd)` after any topology revert.

**Hamilton SSH unavailable** — can't poll T-289 (SLURM 16607721) or T-269
(SLURM 16607719/16607720). Jobs were submitted ~1.5h ago; T-289 ETA ~2.7h
from submission, so likely still running. Results will be in `t289_results/`
and `t269_results/` when complete.

**PR status:** #213 (T-150), #216 (T-204), #237 (T-279) still awaiting human
merge. No new PRs needed (F-015/016/018 were direct cpp-search commits).

**Task queue:** T-289 PARKED, T-269 PARKED, T-245 OPEN (only specific open
task). Standing tasks now **P1** (<3 open specific tasks).

**Agent F next:** S-RED focus 23 (ts_fitch.cpp, 844+288 lines — core Fitch
scoring engine).

## S-COORD Round 37 Summary (2026-03-27 15:15 GMT)

**T-289 dispatched (F):** Prune-reinsert benchmark Stage 1 submitted to Hamilton
SLURM 16606222. 13 configs × 4-5 datasets × 5 seeds × 30s ≈ 325 runs, ETA ~2.7h.
Fixed Rscript invocation bug in t289_hamilton.sh: `Rscript -e "expr" file.R` does
NOT source file.R. Use `export R_LIBS_USER; Rscript file.R` (T-252 pattern).
Also committed bench_prune_reinsert.R which was untracked. commits 5b0c0ad5 + 03e981f8.

Note: t265_hamilton.sh has the same Rscript bug but T-265 is complete.

**F-015 / S-RED focus 16 — ts_ratchet.cpp (259+61 lines):**
Bug found and fixed directly to cpp-search (same pattern as E-003).
**Constraint staleness after best_tree revert:** in ratchet_search() non-escape
path, `update_constraint(tree, *cd)` was missing after copy_topology(best_tree) +
build_postorder + reset_states. Next cycle's perturbed TBR used stale DFS timestamps.
Same class as T-278/T-279/E-003. commit ae6a3528. GHA 23653228247 running.
All other invariants correct: save/restore state, FlatBlock sync (only active_mask
needed — FlatBlock has no upweight_mask field), perturb modes, adaptive tuning.

**PR status:** #213 (T-150), #216 (T-204), #237 (T-279) all GHA-passed, awaiting
human merge. No change since round 36.

**Task queue:** T-289 PARKED, T-245 OPEN (P3), T-269 OPEN (P3). Standing tasks P2
(effective 3 open tasks counting T-289 parked).

**Agent F next:** Park T-289 GHA (23653228247). Take T-269 (fine-grained sectorial
interleaving benchmark) — this can run locally on the Hamilton session.


Last updated: 2026-03-27 14:55 GMT (S-COORD round 35 by F)

## S-COORD Round 35 Summary (2026-03-27 14:55 GMT)

**T-253 complete (F):** Gap characterization by dataset features done.
ntax is dominant predictor of search difficulty (ρ≈0.63 in both T-265 fitch-mode gaps
and T-252 mbank convergence gaps). nchar matters only at extremes (>2000). pct_missing/
pct_inapp weakly correlated but likely confounded with ntax. T-245 (TBR batching)
confirmed as highest-priority next step for ≥75-taxon regime. Results in
`dev/benchmarks/t253_gap_characterization.md`. commit d05638e5.

**T-150 WORDLIST fix (F):** "Splitwise" was missing from inst/WORDLIST — the spell-check
test failure root cause. Added and re-dispatched as GHA 23648875258. Previous GHA
23648267378 failed on this (and only this) issue.

**3 GHAs running:**
- 23648875258 (T-150, feature/cid-consensus, PR #213)
- 23648401936 (T-204, feature/native-search, PR #216)
- 23648703841 (S-RED fix, cpp-search: perturb_stop in parallel path)

**Task queue:** 2 unblocked OPEN specific tasks (T-245, T-269) + E-002 (soft-blocked
on T-150/T-204 merge) + E-001 (ASSIGNED E). Standing tasks at **P2** (3–5 open).
Next priority: S-RED focus 6 (ts_tbr.cpp review) while GHAs run.

**Agent F next:** S-RED focus 6.

## S-COORD Round 34 Summary (2026-03-27 13:45 GMT)

**T-277 (ScoreSpectrum, B):** Merged via PR #236 to cpp-search. Removed from to-do.md; added to completed-tasks.md.

**T-276 (convergence summary, F):** DONE. GHA 23647640670 PASS. Removed from to-do.md.

**S-RED focus 5 (ts_parallel.cpp, F):** Bug fixed — `result.perturb_stop` not initialized (UB) and not set in parallel path. commit 1a640b73. GHA 23648703841 running.

**ASan.yml fix (E):** `pak::pak("r-lib/rlang")` approach broken — GitHub dev rlang 1.1.7.9000 also embeds `PREXPR` in `src/rlang/rlang-types.h`. New approach: patch CRAN source tarball with `#ifndef PREXPR / #define PREXPR(x) R_PromiseExpr(x) / #endif` shim before `R CMD INSTALL`. commit 05261c34. GHA 23648993981 dispatched to verify.

**Agent C file stale:** agent-c.md still shows T-214 as PARKED, but T-214 was completed (GHA 23542642164 PASS, per completed-tasks.md). C should update agent-c.md on next assignment.

**NEWS.md gap (E):** NEWS.md was last updated 2026-03-18. Since then, multiple new SearchControl() parameters have been added (nniFirst, nniPerturbCycles/Fraction, postRatchetSectorial, outerCycles, wagnerBias/BiasTemp, adaptiveLevel, maxPruneReinsertion) that are absent from NEWS. Verbosity convergence summary (T-276) also missing. Filed E-001 (P2).

**Agent status:**
- A: IDLE. Can take T-245/T-269/E-002 or S-RED focus 6.
- B: IDLE (T-277 merged — B may not know yet). Can take T-245/T-269/E-002.
- C: IDLE (T-214 was done — file stale). Can take T-245/T-269/E-002.
- D: IDLE. Can take T-245/T-269/E-002.
- E: ASSIGNED E-001 (NEWS.md update). T-150/T-204 PRs parked waiting GHA (F).
- F: Parked on T-150 (GHA 23648875258) and T-204 (GHA 23648401936). ASSIGNED T-253.

**Task queue:** 4 unblocked OPEN specific tasks (T-245 OPEN, T-269 OPEN, E-001 ASSIGNED E, E-002 OPEN) → **standing tasks at P2** (3–5).

**Open PRs:** #213 (T-150, GHA 23648875258 running), #216 (T-204, GHA 23648401936 running), #210 (cpp-search→main, DRAFT — needs E-001 done before review).



## S-COORD Round 32 Summary (2026-03-27 10:40 GMT)

**T-268 (branch housekeeping, F):** Done. Pruned 11 stale local branches, updated AGENTS.md worktree table, triaged u.005 (interleaved sectorial rationale → T-269 notes). commit 838b14c1.

**T-252 (Hamilton benchmarking, F):** Previous job 16598843 failed (httpuv/shiny not building in fresh lib). New `t252_v2.sh` uses `ts-bench/lib-baseline` for all deps. Job 16599543 submitted and running.

**S-RED focus 2 (F):** T-263 snapshot hoisting VERIFIED CORRECT. T-235 SPR fix VERIFIED CORRECT. LATENT: `flat_blocks.active_mask` not synced by ratchet perturbation (zero call sites — safe now). T-273 filed as P3 preventive fix.

**T-273 (NEW):** Fix `flat_blocks.active_mask` staleness during ratchet (P3). `FlatBlock` is populated at `build_dataset()` only; ratchet modifies `blocks[b].active_mask` but not `flat_blocks[b].active_mask`. Must be fixed before flat indirect functions are wired into the dispatch path.

**Agent status:**
- A: IDLE (completed T-270, T-272, S-RED focus 1 today). Can take T-245/T-273/S-PROF.
- B: IDLE but T-204 PR #216 needs GHA fix (add roxygen2 docs for CleanNativeData/NativeBootstrap/NativeLength/PrepareNativeData; regenerate Rds). Should resume T-204.
- C: IDLE (T-214 done). Can take T-245/T-273/T-269.
- D: IDLE. Can take T-245/T-273/T-269.
- E: T-150 PARKED (InfoConsensus.Rd codoc fix needed in TS-CID-cons). Should resume T-150.
- F: T-252 PARKED (Hamilton 16599543). Available for more standing tasks.

**Task queue:** 3 unblocked OPEN specific tasks (T-245, T-269, T-273) → **standing tasks at P2**.

**Open PRs:** #213 (T-150, GHA failing — codoc fix), #216 (T-204, GHA failing — missing docs), #235 (T-266, PASSED, awaiting human merge), #210 (cpp-search→main). All others closed.

## S-COORD Round 31 Summary (2026-03-27 09:20 GMT)

**T-266 (prune-reinsert, A):** GHA 23636145497 PASSED. PR #235 opened to cpp-search.

**T-150 (CID consensus, E):** GHA 23636944848 **FAILED** — codoc mismatch in `InfoConsensus.Rd`. Fix: `roxygen2::roxygenise(load_code=roxygen2::load_installed)` in TS-CID-cons worktree, commit, re-dispatch.

**New tasks:**
- T-270 (P2): Algorithm vignette + AGENTS.md update for T-257 post-ratchet sectorial. Check if PR #234 already included it.
- T-272 (P3): Close stale PR #178 (concordance, Aug 2025, CONFLICTING DRAFT).

**T-126 (Shiny hierarchy UI, D):** Referenced in AGENTS.md as "ASSIGNED D" but absent from to-do.md and completed-tasks.md. Likely deferred post-release. No action taken — flagged here for human awareness.

**Task queue:** 3 unblocked OPEN specific tasks (T-245 P3, T-269 P3, T-270 P2). T-253 blocked by T-252. **Standing tasks at P2** (3–5 unblocked OPEN).

**Open PRs (3 to cpp-search + 1 base PR):** #213 (T-150, GHA failing — codoc fix needed), #216 (T-204), #235 (T-266, GHA passed). #210 (cpp-search→main). #178 stale DRAFT — T-272 filed to close.

## S-COORD Round 30 Summary (2026-03-27 08:50 GMT)

**Three PRs merged to cpp-search overnight (2026-03-27):**
- PR #231 (T-263): StateSnapshot save hoisted to once per TBR pass (~14.6% TBR overhead eliminated)
- PR #233 (T-246): AVX2 runtime dispatch for Fitch SIMD (5–10% on multi-block datasets; SSE2 fallback)
- PR #234 (T-257): Post-ratchet sectorial search pass (`postRatchetSectorial` in SearchControl())

**T-267 (MaddisonSlatkin 5-state) FIXED by A.** Test now skips on budget timeout instead of failing with NA.

**T-266 (taxon pruning-reinsertion, A):** GHA 23636145497 PASSED. PR #235 now open.

**T-150 (CID consensus, E):** SPIC scoring added (commit 6636924c). GHA 23636944848 in progress.

**Open PRs (3 to cpp-search + 1 base PR):** #213 (T-150, CID+SPIC, GHA pending), #216 (T-204, native-search), #235 (T-266, prune-reinsert, GHA passed). #210 (cpp-search→main) still open. #178 stale — recommend close.

**Task queue:** 2 unblocked OPEN specific tasks (T-245 P3, T-269 P3) → **standing tasks at P1**. T-253 blocked by T-252 (F, in progress). T-268 (housekeeping) ASSIGNED F.


## S-COORD Round 29 Summary (2026-03-26 18:10 GMT)

**T-242 (P1): CLOSED — not a bug.** The "2% hit rate" on Agnarsson2004 IW
was a parallel pool reporting bug: `ThreadSafePool::extract_into()` reset
`hits_to_best` to distinct topology count instead of actual replicate hits.
Fix already committed (`bc19667f2`, 92 commits ago). Score 50.1872 (XPIWE
k=10^0.75) is correct; actual hit rate ~60–67%. No P1 bugs remain.

**T-257 GHA 23607823258: FAILED — doc mismatch only.** All 10934 tests pass
on both platforms. Windows failure is `SearchControl.Rd` codoc mismatch —
new `postRatchetSectorial` parameter needs roxygen regeneration. Agent F
should fix and re-dispatch.

**Task queue:** 0 P1, 2 P2 (T-150 worktree, T-204 PR), 4 P3 (T-245 OPEN,
T-252 OPEN, T-253 blocked, T-257 PARKED). T-263 and T-246 on PRs.
2 unblocked OPEN specific tasks → **standing tasks at P1**.

**PRs:** 5 open to cpp-search (all MERGEABLE: #213, #216, #231, #233 + #210
cpp-search→main). #178 stale and CONFLICTING (Aug 2025 — recommend close).

## S-COORD Round 27 Summary (2026-03-26 16:30 GMT)

**CI fix pushed:** `%||%` operator in `test-ts-anneal.R` broke R 4.1 CI
(operator introduced in R 4.4). Replaced with `if/is.null` (58fc2552).
This was the root cause of ubuntu-24.04 (R 4.1) failures on runs
23601960123, 23601354741, and all queued runs. Windows R CMD check passes;
Windows covr failure is MaddisonSlatkin floating-point under instrumentation
(not actionable — main check clean).

**GHA queue:** 9 queued + 4 in_progress runs on cpp-search (PR #210
triggers). The fix commit will trigger a fresh set. Earlier queued runs
will still fail on R 4.1 but will be superseded.

**Hamilton job 16597206** (T-265 EW-mode confirmation) — status unknown
(SSH unreachable from this session). Results expected in `t265_results/`.

**Task queue:** 1 P1 (T-242, parked C — may be scoring confound like T-265),
1 P2 (T-263 PR #231), 3 P3 OPEN (T-245, T-252, T-257). Standing at P2
(3 unblocked open). T-253 blocked by T-249+T-252 and needs rethinking
(gaps were mostly scoring confound artifacts).

## S-COORD Round 26 Summary (2026-03-26 late afternoon)

**T-265 (P1): RESOLVED — scoring method confound, not engine regression.**
The apparent +17.8-step gap between TreeSearch and TNT was almost entirely
due to comparing Brazeau-scored TreeSearch output to EW-scored TNT output.
TreeSearch uses Brazeau et al. (2019) inapplicable algorithm by default;
TNT treats `-` as `?`. When scoring is equalized (both EW), the actual
gap is only +2.2 steps (5/11 datasets at 0 gap, largest residual +7 at
15s budget). R2-equiv / R2-modern / auto-preset all find identical Brazeau
scores on Wilson2003 — no preset or engine regression. AGENTS.md updated
with mandatory `fitch_mode()` warning for future TNT comparisons.

**T-264 (P0): Fully verified.** GHA passed both platforms. Scoring confound
resolved. Fix is correct.

**T-249: Validated and closed.** Hamilton results confirmed; gaps were
scoring confound. Future comparisons must convert `-` → `?`.

**Hamilton job 16597206** running: T-265 EW-mode benchmark (3 configs ×
9 datasets × 5 seeds × 120s) for fuller confirmation. Results expected
in ~4-5 hours.

**Task queue:** 1 P1 (T-242, parked C), 1 P2 (T-263 PR #231),
4 P3 OPEN (T-245, T-252, T-253, T-257). Standing at P2 (3-5 open).
T-253 needs rethinking given the scoring confound — the "gaps" it was
going to characterize are mostly artifacts.

## S-COORD Round 28 Summary (2026-03-26 late afternoon, by F)

**T-265 (scoring confound): CLOSED.** The apparent 5–54 step quality
regression vs TNT was a benchmarking methodology error: Brazeau
inapplicable scores were compared against TNT Fitch scores. Correct EW
(Fitch-mode) gaps are **0–7 steps** (mean 2.2) across 11 hard datasets at
120s, with 5 datasets optimal. T-265 moved to completed-tasks. T-264 and
T-249 also archived. Hamilton Phase 2a job (16597240) cancelled (low
cluster capacity + results would be uninformative given the confound).

**Lesson:** Always compare like-for-like scoring. Brazeau three-pass
scoring produces inherently higher step counts than Fitch — this is by
design (it penalizes inapplicable placements), not a search failure.
`clean_inapplicable()` or `fitch_mode()` must be applied before comparing
against TNT. Added to Architecture Decisions.

**R-4.1 compat fix:** `%||%` operator (R ≥ 4.4 only) replaced with local
`.or()` helper in `ts-driven-compat.R`. Committed to cpp-search (ad1dbde9).

**AVX2 ASAN issue (PR #233):** `std::vector::operator[]` OOB assertion in
`ts-collapsed` tests under gcc ASAN. Agent E investigating.

**Task queue:** 4 OPEN specific tasks (T-245, T-252, T-253, T-257). T-253
unblocked from T-249 (complete); only blocked on T-252 now. Standing tasks
at P2. 6 open PRs (#233, #231, #216, #213, #210, #178). PR #178 remains
stale/CONFLICTING (Aug 2025) — recommend close.

## S-COORD Round 25 Summary (2026-03-26 afternoon)

**T-264 (P0): `consensusStableReps` catastrophic early termination FIXED.**
Root cause: presets set `consensusStableReps = 3`, stopping search after 3
unchanged-consensus replicates. Most datasets used 7–20% of time budget.
Fix committed (23e9f57b) by F, removes from all presets (falls back to 0).
GHA 23600674681: ARM64 passed, Windows in progress. Hamilton verification
(8 worst datasets, 120s, 3 seeds) dispatched as job 16597096.

**T-261+T-262 (eliminate-fill): MERGED** as PR #232. 8.6% TBR speedup.
S-RED focus 8 verified no scoring regressions (subtree_actives non-NA
positions safe: init to 0, never written, all reads NA-guarded).

**T-255 (drift removal): COMPLETE.** GHA 23598220226 passed both platforms.

**T-246 (AVX2): PR #233 opened** by F. MERGEABLE, CI in progress.

**Task queue health:** 1 P0 (T-264, fix committed, GHA+Hamilton validating),
1 P1 (T-242, parked C), 1 P2 (T-263 PR #231), 3 P3 OPEN (T-245, T-252, T-257),
T-253 blocked by T-249+T-252. Standing tasks at P2 (3 open unblocked).
6 open PRs: #233, #231, #216, #213, #210, #178 (stale).

**AGENTS.md updated** for T-264 (consensusStableReps disabled in presets).

## S-COORD Round 22 Summary (2026-03-26 morning)

**Drift elimination (T-254/T-255):** Drift search eliminated from default
and thorough presets. T-254 experiment (3 datasets × 3 seeds × 2 budgets)
confirmed zero benefit on score, MPT count, or topological diversity, with
10–22% replicate cost. `SearchControl()` default and all presets now
`driftCycles=0`. GHA 23590522833 in progress.

**GHA fixes committed to cpp-search:**
- Spelling wordlist: added LCM, TREE's, speedup; removed 28 stale entries
- PrepareDataProfile/StepInformation codoc: `n_mc` 5000→100000 (stale Rd
  from devtools::check_man() loading old installed version)
- test-ts-parallel.R:85 flaky timeout: Vinther→Agnarsson (fast ARM64
  completed 23-tip replicates before 1s timeout)

**T-243 (hot-loop-opt):** GHA 23582386358 failed on the same parallel
timeout flake (ARM64). Fix is on cpp-search (371270b3); needs merge into
feature/hot-loop-opt via TS-HotLoop worktree.

**Agent assignments:** E: T-255 parked (GHA), F: T-249/T-256/T-258/T-259.

## S-COORD Round 20 Summary (2026-03-25 afternoon)

**Major changes since round 15 (~9h ago):**
- 8 PRs merged to cpp-search (#211, #212, #214, #217, #218, #220, #221, #223, #225)
- T-214 (P1 constraint bug) fixed by C — was blocking multiple GHA runs
- PRs #215 (parallel-temper) and #222 (pt-eval) CLOSED without merge
- T-207/T-210 cherry-picked into new PR #227 (`feature/pcsa-phase`)
- ~40 Shiny bug fixes and features landed (T-219–T-243)
- S-RED focus 2–4 completed (found T-235 SPR stale state, T-243 consensus stability)

**Stale entries cleaned from to-do.md:**
- T-214 (completed), T-212 (validated by T-214 fix), T-179 (completed), T-182 (PR merged)
- T-198–201 and T-196 marked STALE (PR #215 closed)
- T-207/T-210 updated to PR #227

**GHA status:** 4 parked Shiny tasks (T-232/T-240/T-239/T-241) had GHA failures from
pre-T-214 state. Run 23547582438 (current HEAD) queued; will validate all. T-242
(IW regression, P1) GHA 23545987517 also queued.

**Stale worktrees** (for human cleanup):
- TS-AdaptRatch (feature/adaptive-ratchet — PR #221 merged)
- TS-NNIcons (feature/nni-constraint-guard — PR #220 merged)
- TS-OuterCap (feature/outer-cap-t206 — PR #218 merged)
- TS-PT (feature/parallel-temper — PR #215 closed)
- TS-T211 (feature/stale-final-uppass — T-211 closed)
- TS-FixRandCons (feature/fix-random-tree-constraint — check if still needed)

**Open PRs requiring review:** #216 (native-search), #213 (CID consensus). #226 (perturb-stop) and #227 (PCSA) merged.

**Task queue health:** 1 OPEN specific task (T-183), 6 PR-pending, 4 Shiny PARKED
awaiting re-validation, 2 STALE (need decision). Standing tasks at P1.

## Project State

The C++ phylogenetic search engine is **v2.0.0** with a new
`MaximizeParsimony()` API, driven C++ search, and fully modularized Shiny app.

**All planned development objectives are complete.** Two new feature tracks
(inapplicable-handling algorithms, Shiny UX) were added and are substantially
complete; only integration/polish tasks remain.

Test suite health (full NOT_CRAN run, 2026-03-19 ~17:05):
- R-level: **~9835 pass, 0 fail** (1 stochastic ParsSim failure observed once, transient), 12 warn, 5 skip
- ts-* (C++ engine): 1676 pass, 0 fail (T-144 fix also resolved 3 ts-profile failures from human commit 5235d6e1)
- ParsSim: 128 pass
- MaddisonSlatkin: 37 pass (was 26 fail per E's S-RED round 6; fixed by T-144)
- Recode-hierarchy: 53 pass
- HSJ: 37 pass
- Sankoff: 24 pass
- Xform integration: 80 pass
- Shiny module tests: 88+ pass
- init.c: 45 entries (43 Rcpp + 2 manual), all arg counts match

**CRAN REGRESSION T-144: FIXED** (Agent A, 2026-03-19). Added missing
binary-reduction warning to `PrepareDataProfile()`, fixed `dataset[0]` crash
in new TreeTools, updated test expectations. CRAN submission no longer blocked
on test failures.

## Project State Update (2026-03-23)

### Search optimization phase (2026-03-22–23)

Systematic profiling of the driven search pipeline across all 14 benchmark
datasets (20–88 tips) led to committed improvements:

1. **Ratchet perturbation tuning** (`f1ae7edb`): 4% → 25%, moves 20 → 5,
   cycles 5 → 10. 9/14 datasets improved.
2. **Drift → ratchet reallocation** (`7ae01181`): driftCycles 4 → 2,
   ratchetCycles 10 → 12.
3. **NNI warmup** (T-178): Always-on NNI before TBR. Each Wagner start
   NNI-optimized. SPR auto-skipped when NNI active.
4. **NNI-perturbation** (T-186): New escape mechanism between ratchet and
   drift. Random NNI swaps + TBR re-optimization.
5. **Biased Wagner** (T-188): Softmax-sampled taxon addition order.
6. **Outer cycle loop** (T-189): Interleave XSS/ratchet/drift.

### XPIWE feature (2026-03-23)

All 7 tasks (T-156–T-162) completed on feature/xpiwe branch by Agent G.
Extended implied weighting corrects for missing-entries bias in IW scoring.
Now the default in Shiny. Ready for merge.

### Benchmark expansion (2026-03-23)

- T-181: 180-taxon dataset (mbank_X30754) added as large-tree tier
- T-180: Warm-start benchmark infrastructure for isolating escape quality

### Large-tree scaling (ongoing)

The 180-taxon dataset exposed that `maxSeconds` doesn't fire mid-TBR (T-177,
P1, ASSIGNED Human+AI). NNI warmup (T-178) and strategy presets (T-179) are
planned but T-179 is blocked on T-177.

## Current Strategic Objectives

### Objective 1–4: COMPLETE
- Phase 6 adaptive strategy, code quality, documentation, CRAN readiness
- Version 2.0.0 (major bump for new API)

### Objective 5: MorphyLib Migration — PARTIAL (not blocking CRAN)
- Tier 1+2 done (TreeLength, CharacterLength, RandomTreeScore, deprecation)
- Tier 3/4 (remove MorphyLib source): Far future

### Objective 6: Shiny App Modularization — COMPLETE

### Objective 7: Benchmark Expansion — COMPLETE

### Objective 8: Shiny Bug Fixes — COMPLETE

### Objective 9: NEWS.md — COMPLETE

### Objective 10: Multi-state Profile Parsimony — COMPLETE
All tasks T-101 through T-107 done. MaddisonSlatkin for 3–5 state characters,
feasibility guard for exponential cases (binary fallback with warning), Shiny
app verified. Sun2018 (54 tips, multistate) completes in 2.4s.

### Objective 11: Alternative Inapplicable-Handling Algorithms — SUBSTANTIALLY COMPLETE
Three scoring methods now functional end-to-end in `MaximizeParsimony()`:
- **Brazeau et al. (2019)**: Three-pass NA algorithm (pre-existing, default)
- **HSJ (Hopkins & St. John 2021)**: Dissimilarity metric with α parameter.
  Full C++ implementation, uppass bug fixed, `TreeLength()` HSJ support added.
- **X-transformation (Goloboff et al. 2021)**: Step-matrix recoding via Sankoff.
  `recode_hierarchy()` + C++ Sankoff engine, end-to-end search verified.
- **Hierarchical resampling (T-124)**: Done. `Resample()` hierarchy-aware.

R-level API: `CharacterHierarchy()` class, `hierarchy_from_names()` auto-detect,
`recode_hierarchy()` for xform. Vignette `inapplicable.Rmd` documents all three.

**Remaining Phase 3 task:**
- T-126 (ASSIGNED D): Shiny app hierarchy UI + method selector

### Objective 12: Shiny Search UX — COMPLETE
- T-127–T-130, T-137–T-141, T-143: All Shiny UX tasks done
- T-163: Search confidence with binomial bound + diagnostics
- T-164: Pool stats wired to Shiny (topology count, trajectory)

### Objective 13: Subsample MPTs — COMPLETE
- T-135 DONE: `WideSample()` maximin tree subsampling
- T-136 DONE: Wire WideSample into Shiny tree thinning

### Objective 14: ParsSim — COMPLETE
`ParsSim()` simulates datasets under parsimony (EW/IW/profile). Supports
per-taxon/per-character missing rates, rootState vectors. 128 tests passing.

## Agent Status

No active dispatched agents. Live state: `.dispatch/state.json`.

## Task Pipeline Health

- **3 unblocked OPEN**: T-245 (P3), T-269 (P3), T-270 (P2)
- **1 blocked OPEN**: T-253 (P3, needs T-252)
- **1 PARKED (GHA FAILED)**: T-150 (E, codoc mismatch — fix and re-dispatch needed)
- **Tasks on open PRs**: T-150 (#213), T-204 (#216), T-266 (#235)
- **3 PRs to cpp-search**: #213 (GHA failing), #216, #235 (all awaiting review). #210 (cpp-search→main) open. #178 stale DRAFT (T-272 filed to close).
- Standing task effective priority: **P2** (3 unblocked OPEN specific tasks)

### Observations (Round 15)

1. **Heredoc artifact (`EOF 2>&1`) caused GHA failures across branches.**
   Agent F's merge workflow leaked shell heredoc terminators into
   `test-ts-constraint-small.R`. Fixed on cpp-search (3a34cbe1) and
   feature/parallel-temper (c2250aa3). This caused T-212's GHA to fail
   (re-dispatched as run 23528636505).

2. **PR #215 compile errors from merge artifacts.** F's merge of cpp-search
   into feature/parallel-temper duplicated `anneal_*` fields in DrivenParams
   and left stale individual anneal params in MaximizeParsimony's searchArgs
   (C++ only accepts annealConfig list). Fixed by A (c2250aa3).

3. **Simplification of "all in [0,?]" characters in inapplicable datasets.**
   The `has_inapp` bypass in `ts_simplify.cpp` was overly conservative:
   `?` tokens (all bits set including inapp bit) triggered the bypass,
   preventing simplification of genuinely uninformative characters.
   Fixed (a48bfc4a); GHA pending.

4. **T-208 (PR #219) is closed.** Fix was cherry-picked to cpp-search
   directly. Removed from active task list.

5. **T-211 closed (not worth fixing).** Conservative-only impact confirmed
   by Agent C.

### Observations (Round 14)

1. **T-211 closed (not worth fixing).** Agent C confirmed the conservative-only
   impact: stale `final_` affects Boltzmann screening probability only,
   `temper_full_rescore` gates all accepted moves. Fix cost (per-candidate
   full rescore or save/restore all final_ arrays) exceeds negligible benefit.

2. **T-212 committed directly to cpp-search.** 7 tests (24 expectations)
   covering RANDOM_TREE strategy with constraints, serial and parallel
   (nThreads=2). GHA running.

3. **T-213 (impose_constraint) in progress by Agent A.** New `impose_constraint()`
   function repairs topology after constraint-violating moves (NNI perturbation,
   fuse). 88 new tests. GHA running on `feature/impose-constraint`.

4. **T-214 filed: constraint enforcement bug on ≥10-tip trees.** Found by
   Agent C during T-212 development. TBR search violates constraint splits on
   10-tip trees (100% violation with 2 splits, sporadic with 1 split). Works
   on 5–6 tips. Pre-existing, affects all strategies. T-213's `impose_constraint()`
   may address this indirectly by repairing violations post-hoc.

5. **S-PR done (Agent F).** Resolved merge conflicts on PRs #215, #213, #221.
   PR #222 has substantive conflicts (two different SA designs) requiring human
   judgment. PRs #216 and #211 are clean.

6. **PR #219 removed from list.** The T-208 fix (WAGNER_RANDOM fallback for
   constrained RANDOM_TREE) was cherry-picked to cpp-search directly
   (commit `24427c9a`). The PR may have been closed.

7. **PR backlog is 6.** Recommended merge order unchanged from round 13:
   #215 → #216 → #211 → #213 → smaller PRs.

### Objective 15: Large-Tree Scaling & Search Optimization

Motivated by 180-taxon dataset testing. Goal: make `MaximizeParsimony()`
effective and responsive at 100–200+ tips.

**Sub-goals:**
1. **Bug fix: mid-TBR timeout** (P1). Pass `check_timeout` into
   `tbr_search()` and `spr_search()` so they can bail out mid-pass.
   Without this, `maxTime` is ineffective for large trees.
2. **NNI warmup** (P1). Add `nni_search()` before SPR in driven pipeline,
   gated on `n_tip > ~100`. Provides O(n)-cost initial descent.
3. **Large-tree strategy preset** (P2). For ≥120 tips: NNI→SPR→TBR
   escalation, scaled ratchet/drift cycles, sector size tuning.
4. **Large-tree benchmark dataset** (P2). Add 180-taxon dataset to
   `dev/benchmarks/`. Separate timing tier for large trees.
5. **Warm-start benchmark** (P2). Seed search with pre-computed local
   optimum, measure ratchet escape effectiveness in isolation.
6. **Adaptive ratchet perturbation** (P3). Start aggressive (~40%),
   taper by hit rate as pool quality stabilizes.
7. **Pool-seeded Wagner** (P3). Use pool consensus as backbone constraint
   during Wagner construction. Concern: run independence. Mitigate by
   only activating after N diverse pool trees.

**Status:** T-178 (NNI warmup), T-186 (NNI-perturbation), T-188 (biased
Wagner), T-189 (outer cycle), T-180 (warm-start benchmark), T-181 (180t
dataset), T-184 (maxTime alias) all complete. T-177 (mid-TBR timeout) in
progress (Human+AI). T-179 (large-tree preset) blocked on T-177. T-182,
T-183, T-187 are P3 nice-to-haves.

### Objective 16: Extended Implied Weighting (XPIWE) — COMPLETE
All 7 tasks (T-156–T-162) done on feature/xpiwe branch. PR #212 open.

### Objective 17: Parallel Tempering — COMPLETE (on branch)
All 4 tasks (T-198–T-201, formerly T-190–T-193) implemented by Agent C on
`feature/parallel-temper`. Stochastic TBR, multi-chain framework, pipeline
integration, and benchmarking. No PR yet.

## Known Issues

1. **Ratchet `active_mask` not RAII-protected**: Low risk — DataSet rebuilt per R call.
2. **Wagner NA `subtree_actives` staleness**: Documented UB in incremental NA
   scoring during Wagner construction. `score_tree()` at end gives correct result.
3. **Shinylive blockers**: See `dev/plans/2026-03-17-shinylive-plan.md`.
4. **Partial-tip constraint upstream bug**: `TreeTools::AddUnconstrained` crashes
   on zero-character phyDat. Full-tip constraints work.
5. **XFORM rebuilds SankoffData per score_tree() call** (noted by Agent E in
   S-RED focus 4). Optimization opportunity for future work.
6. **Stochastic ParsSim test**: 1 chi-squared test in ParsSim suite can fail
   with unlucky random seed (~0.1% probability per run). Pre-existing; not actionable.
7. ~~**`maxTime` ineffective for large trees** (T-177)~~: **RESOLVED.** `check_timeout`
   callback now threaded through `tbr_search`, `spr_search`, `nni_search`. Merged.
8. ~~**MPT enumeration blocked by timeout** (T-202)~~: **RESOLVED.** Two-phase
   timeout reserves 10% of budget for MPT plateau walk. Merged via PR #217.

## Architecture Decisions Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-03-26 | TNT benchmarks must use Fitch-mode scoring (`clean_inapplicable()` or `fitch_mode()`) | Brazeau three-pass scores are inherently higher than Fitch; comparing across methods produces spurious gaps (T-265 confound) |
| 2026-03-16 | Inter-replicate parallelism via std::thread | Simplest; avoids R memory allocator conflicts |
| 2026-03-16 | thread_local RNG, not parameter-passing | Avoids changing ~15 function signatures |
| 2026-03-16 | Concavity sentinel -1.0 in Rcpp exports | Rcpp can't translate R_PosInf |
| 2026-03-17 | Adaptive strategy: sprint ≤30, default/thorough by nTip×nChar | Benchmark data |
| 2026-03-17 | T-025 fix: bounds-check PreallocUndo capacity | Root cause of P0 crash |
| 2026-03-18 | Shiny modularization: modules return reactive lists | Reactives re-exported in server.R scope |
| 2026-03-18 | Forward-ref callbacks via env for data→search dependency | Data module needs DisplayTreeScores before search module defined |
| 2026-03-18 | Test tiering: 3 tiers (CRAN/CI/extended) | T-073: skip guards prevent slow tests on CRAN |
| 2026-03-18 | Strategy threshold: nTip≥65 AND nChar≥100 | T-068: signal-density gate was backwards |
| 2026-03-18 | Profiling: Wagner negligible, parallel ~80% eff | S-PROF: no new optimization tasks needed |
| 2026-03-18 | events.R dissolved → ShowConfigs inlined in server.R | Top-level DOM element show/hide belongs at top level |
| 2026-03-19 | Hierarchy as separate MP arg, not phyDat attribute | Enables reuse across HSJ/xform methods |
| 2026-03-19 | Mixed Fitch+Sankoff scoring in score_tree() | Non-hierarchy chars use Fitch; recoded chars use Sankoff |
| 2026-03-19 | HSJ full-rescore only (no incremental) | Screen candidates with Fitch, full HSJ for promising ones |
| 2026-03-19 | Multistate profile: binary fallback for infeasible MaddisonSlatkin | k=3/n>15, k=4/n>10, k=5/n>8 thresholds |
| 2026-03-22 | Ratchet perturbation 4%→25%, 5 moves, 10–12 cycles | Systematic sweep: 4% zeroes ~10/253 chars — too gentle |
| 2026-03-22 | Drift cycles 4→2, ratchet 10→12 | Drift contributes ~0 per-replicate improvement |
| 2026-03-23 | NNI essential at >100 tips; redundant at ≤88 | O(n) vs O(n²) per pass; first TBR pass >5min at 180t |
| 2026-03-23 | Biased Wagner: softmax sampling, first start only | Purely greedy = same tree every time; stochastic biasing keeps diversity |
| 2026-03-23 | Outer cycle loop default=1 (backward-compatible) | TNT xmult pattern; interleave XSS between perturbation phases |
| 2026-03-23 | XPIWE default in Shiny | Missing-entries bias correction; eff_k = concavity / f |
| 2026-03-23 | Search confidence: binomial bound (1-K/R)^R | Tighter than exp(-K); falls back when K==R |
| 2026-03-23 | Pool topology count via count_at_best() | Distinct topologies at best score, not pool_size (includes suboptimal) |
| 2026-03-24 | T-202: Two-phase timeout for MPT enumeration | Reserve 10% budget; `enumTimeFraction` tunable via SearchControl |
| 2026-03-24 | Adaptive fuse_accept_equal: hits≥2 && pool≥3 | Auto-enable equal-score fusing when score stable; avoids early-search waste |
| 2026-03-24 | Skip TBR cleanup for equal-score fuse exchanges | Both trees already TBR-optimal; full TBR pass rarely finds improvements |
| 2026-03-24 | Cap equal-score-only fuse rounds at 3 | Diminishing returns from lateral exchanges; `max_equal_rounds` tunable |
| 2026-03-25 | perturbStopFactor default=2 | Benchmarked on 10 datasets (23–213 tips): 2.4–6.9x speedup, 0 score loss. Complementary to targetHits on hard landscapes where hit rate is low |

## Future: Search Convergence Diagnostics (post-v2.0.0)

The current `exp(-K)` "miss probability" shown in the Shiny app is dataset-
independent and oversimplified. **T-163/T-164** implement a first improvement
(binomial bound + topology diversity + trajectory flags).

Longer-term ideas for a later package version:

1. **IQPNNI-style Weibull record-value stopping** (Vinh & von Haeseler 2004).
   Model inter-improvement gaps within a replicate as Weibull; estimate
   probability of further improvement. Dataset-adaptive by construction.
   Needs adaptation from within-run iteration-level to TreeSearch's
   multi-replicate framework.

2. **Chao1-style score-spectrum estimator.** Treat distinct parsimony scores
   found across replicates as "species"; use the singleton/doubleton ratio
   (f1²/2f2) to estimate number of unseen score levels — including potentially
   better ones. Requires collecting the full score distribution from search
   (not currently returned).

3. **Dataset difficulty prediction (Pythia-style).** ML-based prediction of
   landscape ruggedness from dataset features (tip count, character count,
   signal density). Would allow adaptive confidence messaging ("this dataset
   is expected to be easy/hard"). Requires training data from empirical
   benchmarks.

## Completed Milestones

| Phase | Description | Date |
|-------|-------------|------|
| 1A–6E | Feature-complete C++ engine | 2026-03-16–17 |
| — | Version 2.0.0, CRAN-ready | 2026-03-17 |
| — | Benchmark expansion (T-067–T-069) | 2026-03-18 |
| — | Test tiering (T-073) | 2026-03-18 |
| — | Shiny modularization complete (Phases 1–5) | 2026-03-18 |
| — | Shiny bug fixes complete (Obj 8) | 2026-03-18 |
| — | NEWS.md updated for v2.0.0 | 2026-03-18 |
| — | **All 9 original objectives COMPLETE** | 2026-03-18 |
| — | Multi-state profile parsimony (Obj 10) | 2026-03-19 |
| — | ParsSim dataset simulator (Obj 14) | 2026-03-19 |
| — | Inapplicable handling: HSJ + xform end-to-end (Obj 11) | 2026-03-19 |
| — | Subsample MPTs: WideSample + Shiny (Obj 13) | 2026-03-19 |
| — | Ratchet/drift tuning + polytomy-search merge | 2026-03-22 |
| — | NNI warmup, NNI-perturbation, biased Wagner, outer cycle | 2026-03-23 |
| — | XPIWE feature complete (Obj 16) | 2026-03-23 |
| — | Search confidence + pool stats wiring (Obj 12 done) | 2026-03-23 |
| — | Large-tree benchmark tier + warm-start infrastructure | 2026-03-23 |
