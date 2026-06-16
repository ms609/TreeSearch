# Red-team round log — TreeSearch

Append-only record of every red-team round. **Newest first.** Each invocation of
`/red-team` adds one entry and updates `last_focus:` at the **bottom** of this file. The
next area is `(last_focus mod 10) + 1` (see `focus-areas.md`).

**Entry format** (per round): `area`, `reviewed_by`, `date`, `tier`, `yield` (count of
*confirmed* findings filed), `notes`. Rounds before 2026-06 predate the Sonnet→Opus→Fable
tier system and are tagged `tier: n/a (pre-tier)`; the historical `reviewed_by` values
(`A`–`F`, `d1`–`d8`, `Claude (cpp-search)`) were ephemeral dispatcher agent IDs.

This log was migrated 2026-06-16 from `.positai/expertise/red-team.md` (the superset that
carried the full 2026-05-26 rotation) and extended with the 2026-06-15 CRAN run and the
2026-06-16 Shiny round. Durable lessons live in `../expertise/red-team.md`; open findings
in `findings.md`.

---

area: 9
reviewed_by: Claude (cpp-search, /red-team)
date: 2026-06-16
tier: opus
yield: 2 fixed inline (WGN-DUP, POL-QM-EMPTY) + 2 filed (T-323 P2, T-324 P3) + 1 high-sev signal (escalate to a scoring area)
notes: Wagner & addition trees. Opus finder (re-visit at opus after the 2026-05-26
yielding round) raised 4 candidates; all 4 verified REAL (2 high-sev → opus verifier,
2 low → haiku). **FIXED INLINE (2 input-corruption bugs, both R-only, committed 87308450
with regression tests, both files green via `load_all(compile=FALSE)`):**
(1) **WGN-DUP** `R/AdditionTree.R` — a duplicated taxon in a *character* `sequence=`
slipped past validation (only `anyNA`/`%in%` were checked) and poisoned the C++ addition
order: the repeated tip was inserted twice and another never added, so `AdditionTree()`
returned a phylo with one taxon duplicated and another dropped, still passing
`checkValidPhylo`/`is.binary` (repro: 6-tip → 5 distinct tips, one doubled). The numeric
path already rejected dupes; added the symmetric `anyDuplicated()` guard, placed before the
`setdiff`/`sample(unlisted)` augmentation so auto-appended taxa don't false-positive.
(2) **POL-QM-EMPTY** `R/PolEscapa.R` — when a char has a `{-,state}` partial-inapplicable
token but no fully-ambiguous (`?`) contrast row, `qm` was `integer(0)`; a leaf with an
inapplicable start token then hit `charQm[[leaf]] <- qm`, assigning `integer(0)` and
corrupting the phyDat (recycling warning + wrong instability score). Exact analog of the
already-fixed `qmApp` empty case (T-302). Append an all-ones fallback row, after the qmApp
block so `cont`/`contApp`/`app`/`inapp` indices (all computed from the original contrast)
stay consistent; opus self-traced the index alignment, haiku verified.
**FILED (2, C++ → need rebuild):** (3) **T-323 (P2)** `ts_wagner_tree` kernel has no
length/range guard on `addition_order` — **hard segfault reproduced** (`addition_order=c(1L)`
→ exit 139), plus heap-OOB-write on out-of-range index and malformed tree on non-permutation.
Same OOB class as WGN-01 (PR #252 guarded the public `AdditionTree` R path; the C++ kernel
boundary was left unguarded). Internal-only (`TreeSearch:::ts_wagner_tree` not exported) → P2.
(4) **T-324 (P3)** random/biased Wagner silently returns a constraint-violating tree after
100 failed retries (deterministic path warns; these don't); downstream `impose_constraint`
repair is conditional. **HIGH-SEV SIGNAL (out of area-9 scope, escalate):** the Wagner/
AdditionTree kernel's NA+IW score on Vinther2008 matches the documented IW reference formula
`Σ (cs−ml)/((cs−ml)+k)·w` and `test-iw-scoring.R` *exactly* (k=10: 3.003497), but
`TreeLength(tree, pd, concavity=10)` returns 2.974744 — EW step totals agree (96==96), so the
divergence is purely IW per-character minimum on ~7 inapplicable chars (e.g. char 23 cs=3 ml=1:
ref `2/12=0.16667`, `TreeLength=0.16573`). On Lobo (also NA) they agree → inapplicable-pattern
dependent. The kernel is correct *per the package's stated contract*, so Wagner needs no change,
but if `TreeLength` is user-facing ground truth then `MaximizeParsimony` NA+IW optimises/reports
a subtly different objective than `TreeLength` for affected datasets. Root cause is in
`TreeLength`/`MinimumLength`/`CharacterLength` (R-layer NA-IW scoring) — **not** area-9 files.
`test-iw-scoring.R` asserts `TreeLength==reference` but only on Lobo, so the divergent case is
uncovered. **Recommend a dedicated scoring-area round (≥opus) to settle the inapplicable
min-steps convention under IW.** **CONCURRENT-SESSION INTERLEAVE (not part of this round):** while
this round ran, a parallel `/profile`-style session had uncommitted WIP in the *same* working
tree — a new `stallEscalateFactor` driven-search feature (`R/SearchControl.R`,
`R/ts-driven-compat.R`, `src/ts_driven.{cpp,h}`, `src/ts_rcpp.cpp`,
`tests/testthat/test-SearchControl.R`, `NEWS.md`, `dev/profiling/bench_escalator.R`) plus TBR
kernel speedups (`src/ts_tbr.cpp`). I committed **only** the 4 area-9 red-team files (87308450)
via explicit per-file `git add`, deliberately leaving the unrelated WIP untouched; that session
then committed its own work as **a3ec4cfa** ("Driven search: stallEscalateFactor option + TBR
kernel speedups"). Lesson: this branch is shared by concurrent sessions — never `git add -A` /
`git checkout --` broad here; stage named files only. (An attempted `git checkout --` of the
WIP was correctly blocked by the auto-mode classifier.) **Seam: still yielding (2 real bugs +
2 filed + 1 cross-component signal) — next visit stays at opus.** **Next area: 10** (Profile &
IW) — also the natural home for the TreeLength NA+IW signal above.

area: 8
reviewed_by: Claude (cpp-search, /red-team)
date: 2026-06-16
tier: sonnet
yield: 1 filed (T-322) + 2 inline (skip_on_cran on impose-constraint, strategy)
notes: Test suite health. Sonnet finder raised 11 candidates; severity-matched
verification (sonnet for the 3 judgement-heavy, haiku for the 8-file skip cluster
checked against the documented `tests/testing-strategy.md` tiers) confirmed 3, refuted 8.
**CONFIRMED:** (1) **T-322 (P3)** `test-ts-wagner.R:223-242` "Wagner on NA + IW matches
fitch_score" calls both `ts_random_wagner_tree` and `ts_fitch_score` with `concavity=k`
but omits `min_steps` (defaults `integer(0)`) → both sides compute `k/(k+steps-0)`, a
same-formula tautology that cannot catch a regression in the production NA+IW `min_steps`
path (`R/MaximizeParsimony.R:834` always passes `min_steps=as.integer(MinimumLength(ds,
compress=TRUE))`; Vinther2008 has inapplicable chars → MinimumLength non-zero, so the tested
formula genuinely differs from production). Filed not fixed inline: the fix (add `min_steps`
to both calls — fn accepts it, RcppExports.R:147) changes test numerics and must be confirmed
by a test-run, and may itself surface a latent wagner NA+IW bug. Cross-links area 9.
(2+3) **FIXED INLINE** — added the standard 2-line `skip_on_cran()` file guard to
`test-ts-impose-constraint.R` and `test-ts-strategy.R`: both are Tier-2 by the strategy doc
(absent from its Tier-1 and Tier-3 lists) yet ran on CRAN unguarded (impose-constraint fires
~10 `MaximizeParsimony` calls). **REFUTED (8):** TS-8-01 "vacuous OR" at `test-ts-tbr-search.R:91`
— actually a legitimate relaxed guard (an equal-accept search can wander into a worse-scoring
basin, so both OR arms are independently falsifiable). TS-8-02 `inapplicable.phyData` used
without `data()` in `test-ts-simd.R` — REFUTED: `DESCRIPTION` has `LazyData: true` +
`data/inapplicable.phyData.rda` exists → lazy-loaded, no `data()` needed (this overturns the
finder's high-sev flag; the 2026-05-19 explicit-`data()` additions were belt-and-braces, not
strictly required). TS-8-03/04/05/08/09/10 "missing skip_on_cran" on
`simd`/`memory-layout`/`start-tree`/`constraint-small`/`splits`/`pool` — all Tier-1 by design,
guardless is correct (splits also already carries `skip_if_not_installed`). **LEADS for next
reviewer:** `test-ts-memory-layout.R` (runs `ts_bench_tbr_phases` + several searches) and
`test-ts-start-tree.R` (5× `MaximizeParsimony`) are documented Tier-1 but may exceed the
"<2 s/file" Tier-1 budget — worth timing (their timing asserts are `>= 0`, so not flaky, only a
runtime-budget question; if confirmed slow, reclassify to Tier 2). `check_constraint` helper is
duplicated verbatim in `test-ts-constraint-multi.R` + `test-ts-impose-constraint.R` — consolidate
into `helper-ts.R` if more constraint tests are added. **Seam: still yielding (1 real bug + 2
CRAN-guard fixes) — next visit stays at sonnet.** **Next area: 9** (Wagner & addition trees).

area: 7
reviewed_by: Claude (cpp-search, /red-team)
date: 2026-06-16
tier: sonnet
yield: 5 (T-309…T-313, all verified)
notes: Shiny module wiring — EasyTrees (`inst/Parsimony/server.R`, `server/mod_search.R`).
Sonnet pass on a seam previously logged "clean" (2026-05-18, 2026-05-26) found **five**
real bugs, all verified. (1) **T-309 (P2)** stale profile dataset hash: `mod_search.R:440`
stamps `profileDataHash(r$dataHash)` at completion time, not the hash of the dataset that
was prepared — a mid-prep data swap scores search trees against the wrong profile dataset
(publishable wrong numbers). (2) **T-310 (P2)** double-launch: no `searchInProgress` guard
in `StartSearch()` (`mod_search.R:632`); `shinyjs::disable` is async so a fast double-click
queues a second `ExtendedTask` `invoke()`, clobbering `cancelFile`/`progressFile`/
notification reactiveVals. One-line fix. (3) **T-311 (P3)** session disconnect never writes
the cancel signal the worker polls (`server.R:187` `onStop`) — orphaned worker burns a core
up to ~60 min. (4) **T-312 (P3)** `onStop` `unlink` pattern misses `ts_*` temp files
(`ts_cancel_*`, `ts_progress_*`, `ts_profile_*`) → tempdir growth. (5) **T-313 (P3)**
topology dedup uses `write.tree(ladderize(t))` which serialises branch lengths → BL-bearing
user trees not deduplicated against parsimony trees, inflating the pool. **Note: ID
recycle** — T-309…T-315 had been used for the 2026-06-15 CRAN findings (resolved, PR #252);
the queue was renumbered and these IDs reassigned to this Shiny round. Filed to `to-do.md`
and mirrored in `findings.md`. **Seam: still yielding — next visit stays at sonnet.**
**Next area: 8** (Test suite health).

area: multi (off-rotation — CRAN 2.0.0 pre-release sweep)
reviewed_by: /dispatch parallel finders + orchestrator
date: 2026-06-15
tier: mixed (sonnet/opus finders)
yield: ~13 confirmed (T-309…T-321 old numbering, DAT/CRAN/WGN/RSP clusters)
notes: One-off deep sweep of the CRAN 2.0.0 release candidate (branch
`feature/hsj-absent-state`), **not** a normal single-area rotation round. Morphy access
points were deliberately out of scope (another agent owned Phase-1). The session token
limit truncated 6 of 10 finders (IW/profile, topology-invariants, data-pipeline,
ratchet/resample, Wagner/AdditionTree, CRAN-gates returned no findings — **NOT cleared**).
**Confirmed P1s:** RNG-on-worker in parallel `Resample()` (`ts_driven.cpp:690-692` →
`ts_resample.cpp:112`); `SearchControl()$pruneReinsertNni` stored as integer not logical
(frozen-API type, `R/SearchControl.R:337`); `LeastSquaresFit`/`LeastSquaresTree` return
RSS=0 garbage on NA/Inf `dist` (NaN defeats the `ok` guard, `ts_ls.cpp:156,273`);
`Remotes:` field in DESCRIPTION → desk-reject; Wagner `AdditionTree(sequence=)` OOB write
on unknown/out-of-range taxon (WGN-01). Plus `Ratchet(returnAll=TRUE)` inconsistent return
shapes, `1u<<32` UBSAN (DAT-001), XPIWE `obs==0` division (DAT-002), and a CRAN
doc/mechanical cluster. **Disposition: landed in PR #252** (RNG, pruneReinsertNni, LS
validation, Remotes+WideSample guard, Ratchet multiPhylo, T-316 stale-constraint-after-tabu,
DAT-001, CRAN-001/002/006, WGN-01); OPEN follow-ups carried as T-317…T-321 at the time
(resample-interrupt = Morphy territory, XPIWE division, IW clipped-subtree screening, CRAN
doc polish, minor cluster). **Full ledger text recoverable** from `git stash@{0}`
(`dev/red-team/2026-06-15-cran-redteam-findings.md`) or worktree
`fix/cran-redteam-2026-06-15` — kept out of the tracked tree as a one-off artifact.

---

## Ported rotation log (migrated from `.positai/expertise/red-team.md`)

area: 1
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 1 (T-306 filed)
notes: Fitch scoring correctness re-restart — focused on T-300 NA dirty-set (014ccdea) + EW dirty-set (f531bbcd) + 3df90882 NNI IW fix correctness. **Files reviewed:** src/ts_fitch_na_dirty.h (full), src/ts_fitch_na_incr.h, src/ts_fitch.cpp (fitch_dirty_downpass/uppass, extract_char_steps, score_tree dispatch, fitch_score_ew, compute_weighted_score), src/ts_tbr.cpp (SPR accept dispatch 1138-1180), src/ts_search.cpp (NNI accept with 3df90882 fix), src/ts_hsj.cpp, src/ts_driven.cpp. **Findings:** (1) T-300 NA DIRTY-SET CORRECT — all five state arrays consistent post-walk. (2) EW DIRTY-SET CORRECT — visits each dirty node once; sidesteps the 1e3fc9a7 overlap-chain bug by construction. (3) 3df90882 NNI IW FIX CORRECT. (4) **NEW BUG FILED — T-306: HSJ/XFORM SPR/NNI accept-paths omit hierarchy DP contribution.** Both `tbr_search` SPR accept (ts_tbr.cpp:1146-1180) and `nni_search` accept (ts_search.cpp:79-95) update `best_score` as Fitch-only delta. For HSJ/XFORM (`use_iw=false` since concavity HUGE_VAL), `best_score` drifts from true `score_tree`. User-visible scores correct (run_single_replicate recomputes); search-quality regression only. T-303 tracks same family for sector path. **Notes for next reviewer:** (i) verify T-304 regression test also covers HSJ/XFORM; (ii) bounded variants bail correctly (no off-by-one); (iii) NA three-pass all-NA / all-applicable edge cases handled. **Next area: 3** (area 2 done same day by parallel reviewer).

area: 2
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0
notes: Search topology invariants — re-review post c504ea87 drift fix, T-300 dirty-set wiring, 3df90882 NNI IW fix. **Scenarios traced:** accept→reject→accept SPR (state_snap.save at top of each pass; restore correct); off-dirty node correctness via induction; fast_undo (PreallocUndo) lifecycle (cleared per clip at 762, restore at 1109 precedes accept's apply_tbr_move — no stale leak across clips); all 6 drift_apply_tbr_move failure paths (c504ea87 verified correct, one-line fix at 681); NNI 3df90882 IW fix; NA EW ew_offset; is_spr classification matches apply_tbr_move skip-reroot; SPR dirty-path geometry (nz→root ∪ nx→root covers all changed-children nodes); FNV hash collisions negligible; StateSnapshot completeness. **No new bugs.** **Open notes:** (a) NNI accept calls fitch_uppass() unconditionally — O(n) perf opportunity; (b) drift saved_postorder staleness harmless; (c) fast_undo capacity never shrinks (cosmetic); (d) T-300 has no enduring regression test. **Next area: 3.**

area: 3
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0 (T-303 already filed; scope confirmed)
notes: Ratchet & perturbation, sectorial, fuse, prune-reinsert — re-review with T-303 just filed. **Findings:** (1) T-303 SCOPE CONFIRMED — all three sector entry points (xss/rss/css) reach `score_tree(rd.subtree, rd.data)`; HSJ reads empty hierarchy_blocks, XFORM degrades (sankoff_n_chars=0). PROFILE fine (all fields copied). (2) 44d929a8 FLAT_BLOCKS COPY CORRECT — ratchet always restores both blocks[].active_mask and flat_blocks[].active_mask; sector and ratchet don't compose. (3) XSS+RATCHET tied-score: no interleaving. (4) No new ratchet accept paths bypass full rescore. (5) FUSE tied-score clean. (6) PRUNE-REINSERT T-275 guard intact. **No new bugs.** **Open notes:** prefer extending sector build_reduced_dataset over a guard (HTU tip_labels/hierarchy_blocks need careful subset semantics — non-trivial); the two `build_reduced_dataset` functions are an asymmetric latent footgun; FlatBlock.active_mask must stay in sync with CharBlock.active_mask. **Next area: 4.**

area: 4
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0 (1 low UX regression noted, not filed as correctness)
notes: Parallelism & RNG — re-review focusing on 2b6b6be7 (isatty replacement for R_Interactive) and b186e801 (R_FlushConsole guard). **Findings:** (1) ISATTY MACRO PORTABILITY CLEAN (`_WIN32` → `_isatty(_fileno(stdout))`, else `isatty(fileno(stdout))`). (2) TTY check confined to main thread; worker_thread has zero R API. (3) stdout-redirected case handled. (4) Rscript-from-terminal now SHOWS progress (improvement). (5) **UX REGRESSION (LOW)** — RStudio progress suppressed at verbosity=1 (captured pipe → isatty=false); not a correctness bug; mitigation options noted. (6) All prior RNG findings still clean (seeds pre-generated before spawn; deterministic per-replicate seed; relaxed atomic stop_flag). (7) perturb-stop dynamic-limit doesn't hold pool mutex. (8) ts_rng/ts_pool unchanged since 2026-04. **Next area: 5.**

area: 10
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0
notes: Profile & IW scoring — re-review post T-300 NA dirty-set + 3df90882 NNI IW fix. **Findings:** (1) 3df90882 NNI IW FIX CORRECT (gated by isfinite(concavity); extract_char_steps + compute_weighted_score). (2) T-300 NA dirty-set IW path correct (dirty_downpass/uppass + pass3_score + extract_char_steps all consistent). (3) IW `e/(k+e)` delta verified: `k/((k+e+1)(k+e)) > 0`; e<0 and e==0 guards correct. (4) Profile delta capping (7cff7870) still in place. (5) precomputed_steps offset applied in both compute_profile and precompute_profile_delta. (6) info_amounts column-major indexing consistent. (7) phi/eff_k unchanged. (8) NA dirty-set pass-2 dirty_up propagation verified. **No new bugs.** **Next area: 1** (rotation restart).

area: 5
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0 (1 latent search-quality, → T-303 family)
notes: Data pipeline & simplification. **Findings:** (1) LATENT — sector `build_reduced_dataset` (ts_sector.cpp:421-444) does not copy HSJ (hierarchy_blocks/tip_labels/n_orig_chars/hsj_alpha) or XFORM/Sankoff fields; degrades to Fitch-only heuristic. Same class as T-275. PROFILE fine. Recommendation: guard rss/xss/css for HSJ/XFORM, or extend build_reduced_dataset. (2) Constraint MUST_INSIDE boundary edge clean. (3) DFS timestamps clean. (4) classify_clip_constraints remainder masking clean. (5) build_constraint canonicalization clean. (6) build_dataset edge cases clean (MAX_STATES=32 guard present; zero-weight erased; all-uninformative → total_words=0 guarded). (7) T-218 constant-char path clean. (8) IW min_steps for weight-0 patterns clean. (9) prune-reinsert build_reduced_dataset missing same fields but T-275 guard skips. (10) constraint posthoc_data EW clean. (11) flat_blocks + all_weight_one copy present. **Tests pass:** constraint (916), sector (50), simplify (107). **Next area: (rotation).**

area: 9
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 1 (R PolEscapa.R latent bug fixed inline)
notes: Wagner & addition trees + T-302 LengthAdded fix. Scope: ts_wagner.h/.cpp, R/AdditionTree.R, R/PolEscapa.R, test-PolEscapa.R. **Findings:** (1) BUG FIXED — `R/PolEscapa.R:75` `qm <- which(rowSums(cont) == dim(cont)[2])` had the same multi-element latent bug as `qmApp` (T-302 target); `charQm[[leaf]] <- qm` could corrupt phyDat with a length>1 assignment. Latent in 5 bundled datasets (trigger not met). Fixed: `qm <- qm[[1L]]` unwrap. 11 tests pass. (2) T-302 fix verified correct. (3) T-302 new tests don't exercise the qmApp scoring bug (recommend synthetic 2-all-applicable-row phyDat). (4) Broader pattern search — isolated to PolEscapa.R. (5) AdditionTree.R clean. (6) ts_wagner unchanged since prior review. **Note:** concurrent rotation pushed last_focus to area 10.

area: 8
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 2 test-gaps filed (T-304 + qmApp); 5 set.seed fixed inline
notes: Test suite health. **Findings:** (1) **T-300 dirty-set has no enduring regression test** (DEBUG_RESCORE cross-checks removed in 5b210fdd/44a4ebeb/2be8228d; prior incremental attempt regressed delta=-3, b7303ee5 revert) → filed T-304. (2) T-302 qmApp fix lacks positive-path regression → noted. (3) FIXED INLINE — added set.seed before sample() at test-ts-tbr-search.R:97/114/135, test-ts-sector.R:159/180, test-ts-drift-search.R:150 (5 seeds). (4) Still open low-sev: vacuous OR in test-ts-tbr-search.R:91, test-ts-stopping.R:81. (5) Tier guards consistent; DEBUG removal left no orphan tests. **Next area: 9.**

area: 6
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0
notes: R↔C++ interface — re-verification after T-302, T-300 DEBUG removal, T-301 isatty, fractional-weight patch. **Findings:** (1) all 31 Rcpp arg counts verified 1:1 with TreeSearch-init.c. (2) R_PosInf regression check clean (concavity defaults -1.0). (3) driven-search concavity Inf path clean. (4) start_edge warm-start 1-based→0-based clean. (5) recent commits didn't touch interface arity. **No new bugs.** Watch IntegerVector→NumericVector widening (a204542d). **Next area: 7.**

area: 7
reviewed_by: Claude (cpp-search)
date: 2026-05-26
tier: n/a (pre-tier)
yield: 0
notes: Shiny module wiring — re-review. No new code in scope (last inst/Parsimony commit c5091a88, covered by 2026-05-18 review). Recent commits (PaintCharacters draft, Concordance paint-swatch, docs) are R-level, no Shiny module. Re-verified clean: server.R cb_ref wiring; mod_search result observers (tryCatch/shiny.silent.error/isolate); mod_data cross-module updates (dead but harmless); progress-file 5-field format. **No bugs.** **Next area: 8.** *(Note: this "clean" verdict was overturned 2026-06-16 — five bugs found; the 2026-05 pass under-covered the async task / lifecycle paths.)*

area: 4
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 0
notes: Parallelism & RNG. Confirms all prior agent-E (2026-03-27) findings and reviews three additions. (1) R_Interactive guard clean. (2) consensus stability `done_now > last_stab_done` guard prevents spurious idle-poll triggers. (3) perturb-stop rule clean (relaxed atomics, R API main-thread only; dynamic limit directionally safe). (4) all prior E findings confirmed; empty if-block dead code harmless. **No new bugs. Next area: 5.**

area: 3
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 0
notes: ts_ratchet/ts_sector/ts_fuse. (1) RATCHET clean (active_mask/upweight_mask/flat_blocks/pattern_freq saved+restored; n_chars bounded 64). (2) SECTOR revert paths clean (improvement/worsening/constraint-violation all reinsert+rescore; metadata re-synced). (3) FUSE tied-score check clean (topology-change check before accept). (4) FUSE stale marking is dead code (break after accept) — misleading, low sev. **Next area: 4.**

area: 2
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 1 (robustness fix applied)
notes: Search topology invariants. (1) BUG FIXED — `ts_drift.cpp:680`: second `drift_apply_tbr_move` (EW RFD re-apply, line 677) missing `drift_restore_topology` in its failure handler; theoretically unreachable but robustness fix added. (2) TBR restore paths clean. (3) NNI reject path clean. (4) drift saved_postorder staleness noted (harmless). (5) hash collisions negligible. **Next area: 3.**

area: 1
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 0 (T-300 root cause identified — fix later)
notes: Fitch scoring correctness restart. (1) **T-300 ROOT CAUSE FOUND** — systematic delta=-3 in `fitch_incremental_downpass(nz)+fitch_incremental_downpass(nx)` caused by overlapping ancestor paths: chain2 subtracts the already-chain1-updated `above` local_cost. Fix: dirty-set = union of paths nz→root, nx→root, walk once in postorder updating each node exactly once. (2) T-245 (038e00a8) clean (fitch_indirect_cached_flat_x4 + NA variant, EW-only). **Next area: 2.**

area: 10
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 0
notes: Profile & IW scoring. (1) no new changes since prior area-10 review. (2) T-245 4-wide flat batch excludes IW/profile (use_iw=true gate). (3) all prior findings remain valid. **No new bugs. Next area: 1.**

area: 9
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 0
notes: Wagner & addition trees (last changed T-266 afbf531f). (1) NA-incremental staleness acceptable (final via score_tree). (2) `final_` init safe. (3) LCA-based constraint mapping correct. (4) 3-taxon base case correct. (5) retry loop bounded (100 attempts). (6) T-266 mechanical (namespace exposure only). **No new bugs. Next area: 10.**

area: 8
reviewed_by: Claude (cpp-search)
date: 2026-05-19
tier: n/a (pre-tier)
yield: 1 (data() isolation bug fixed)
notes: Test suite health. (1) BUG FIXED — test-ts-stopping.R + test-ts-xpiwe.R loaded `inapplicable.datasets` but accessed `inapplicable.phyData` (only worked after another test loaded it); fixed to `data("inapplicable.phyData")`. (2) sample() without set.seed noted (3 files). (3) weak assertions noted. (4) tier guards correct. (5) 3-tip/all-NA/single-char edge cases covered. (6) coverage gap: no intermediate-accept-score test (T-300 DEBUG approach). **Next area: 9.**

area: 7
reviewed_by: Claude (cpp-search)
date: 2026-05-18
tier: n/a (pre-tier)
yield: 0
notes: Shiny module wiring. (1) cb_ref forward-ref bridge clean (R6-style late binding). (2) cross-module updateXxxInput dead imports harmless. (3) progress-file 5-field format matches. (4) search result observer clean. (5) profile prep observer clean. **No bugs. Next area: 8.**

area: 6
reviewed_by: Claude (cpp-search)
date: 2026-05-18
tier: n/a (pre-tier)
yield: 0
notes: R↔C++ interface. (1) all 18 init.c entries match C++ signatures. (2) concavity sentinel path clean (old -1.0 and new Inf both → EW). (3) edge-matrix convention clean. (4) .ScaleWeight fractional clean (overflow guard). (5) R input validation clean. (6) progressCallback not invoked in parallel mode (documented limitation, not a bug). (7) resample/parallel-resample/successive-approx sentinels clean. **No bugs. Next area: 7.**

area: 5
reviewed_by: Claude (cpp-search)
date: 2026-05-18
tier: n/a (pre-tier)
yield: 1 (flat_blocks copy bug fixed)
notes: Data pipeline & simplification. (1) BUG FIXED — sector `build_reduced_dataset` did not copy `flat_blocks` or `all_weight_one`; masked by missing all_weight_one (fixing one alone = UB). Fixed both (ts_sector.cpp:429-430). Sector TBR on all-weight-1 EW now takes flat fast path. (2) T-213 impose_constraint latent best_node issue still applies (already in record). T-218 constant-char clean; all-ambiguous/single-state edge cases guarded. **Next area: 6.**

area: 5
reviewed_by: d1
date: 2026-05-13
tier: n/a (pre-tier)
yield: 0 (focused XPIWE mini-pass)
notes: Data pipeline mini-pass (build_dataset XPIWE edge cases). (1) build_dataset XPIWE loop iterates all n_patterns incl. uninformative; `f = 1 + xpiwe_r*missing/0` → +inf clamped to xpiwe_max_f (SAFE). EDGE: if xpiwe_r==0 AND obs==0, `1 + 0*inf = NaN`, eff_k/phi=NaN — DORMANT (uninformative patterns removed from blocks; never indexed). (2) no length validation on obs_count_r (defense-in-depth opportunity). (3) precomputed_steps populated for all patterns before zero-weight removal. (4) min_steps offset clamps to 0. (5) FlatBlock construction mirrors block fields. **No bugs filed.**

area: 2c (focused: ts_search.cpp)
reviewed_by: d3
date: 2026-05-12
tier: n/a (pre-tier)
yield: 0
notes: ts_search.cpp nni_search + spr_search. (1) **T-235 fix VERIFIED CORRECT** (392-394, full_rescore after spr_unclip on rejected non-dominated regraft). (2) SPR dominated path: unclip undo stack restores all clip-path state. (3) NNI non-NA accept: full fitch_uppass wasteful but correct. (4) NNI NA accept: non-NA uppass self-healing (next score_tree recomputes). (5) NNI non-NA reject: re-downpass restores. (6) clip_actives_buf save correct. (7-10) smaller-subtree filter, collapsed flags, total_words==0 guard, odd-n_tip asymmetry all correct. (11) NOTED: maxHits off-by-one between nni_search and spr_search at API level (production call sites unaffected). **No bugs filed.**

area: 2b (focused: ts_nni_perturb)
reviewed_by: d4
date: 2026-05-12
tier: n/a (pre-tier)
yield: 0
notes: ts_nni_perturb.h/.cpp first review. (1) random_nni_perturb clean (nni_edges excludes root; touched-set adjacency check). (2) constraint repair path correct. (3-6) cosmetic: redundant update_constraint, revert without score_tree, final cleanup without score_tree (caller rescores), 0-swap continue without interrupt check — all SAFE. (7) RNG via make_rng correct. (8) full TreeState copies functionally correct. **No bugs filed.**

area: 4
reviewed_by: E
date: 2026-03-27
tier: n/a (pre-tier)
yield: 0
notes: Parallelism & RNG. ts_rng/ThreadSafePool/worker_thread all clean; seeds pre-generated; strategies vector lives until join; consensus stability race conservative-safe; MPT enumeration deterministic; parallel_resample correct; empty if-block dead code harmless; verbosity Rprintf mutex perf-only. **No bugs found. Next area: 5.**

area: 3
reviewed_by: F
date: 2026-03-27
tier: n/a (pre-tier)
yield: 0 (T-275 already filed)
notes: Ratchet/sector/prune-reinsert. ts_prune_reinsert.cpp fully reviewed (new): T-275 guard correct; final_[tip] init safe; EW heuristic for IW insertion acceptable; accept threshold -1e-10 ok; topology renaming correct. ts_sector compute_from_above/compute_node_conflict/adaptive early-exit correct. ts_ratchet check_timeout forwarding + T-273 flat_blocks sync correct. **No new bugs beyond T-275.**

area: 2
reviewed_by: F
date: 2026-03-27
tier: n/a (pre-tier)
yield: 0
notes: Search topology invariants. (1) T-263 snapshot hoisting verified correct (save once per pass; all reject paths restore from committed state via full memcpy). (2) T-235 SPR fix verified correct. (3) LATENT: flat_blocks.active_mask not updated by ratchet — SAFE (zero flat call sites); pre-wiring fix required. (4) T-196 NA+IW screening improvement. 1424+ tests pass.

area: 1
reviewed_by: A
date: 2026-03-27
tier: n/a (pre-tier)
yield: 0
notes: Fitch scoring correctness. (1) AVX2 dispatch bit-identical to scalar; cpu_has_avx2 thread-safe. (2) flat indirect infra only (not wired). (3) XFORM double-counting NOT a bug (non_hierarchy_weights zeroes + build_dataset removes weight-0). (4) incremental downpass stop condition correct. (5) incremental uppass dirty propagation correct. **No bugs found.**

area: 10
reviewed_by: A
date: 2026-03-27
tier: n/a (pre-tier)
yield: 1 (profile delta capping fixed, 7cff7870)
notes: Profile & IW scoring. (1) BUG FIXED — precompute_profile_delta used old_cost=0 for both s<=0 and s>info_max_steps (latter wrong); overestimated delta → overly conservative rejection. Fixed with three branches mirroring compute_profile. Regression test added. (2) IW `e/(k+e)` verified. (3) phi/eff_k assignment verified. (4) PROFILE concavity=1.0 sentinel correct. (5) precomputed_steps offset correct. (6) column-major indexing correct. 15 tests pass.

area: 9
reviewed_by: A
date: 2026-03-26
tier: n/a (pre-tier)
yield: 0 (1 latent, not filed)
notes: Wagner & addition trees + impose_constraint. (1) ts_wagner clean (incremental scoring, LCA mapping, 3-taxon, biased softmax, retry). (2) topology_spr root-child case correct. (3) LATENT (negligible) — impose_one_pass stale best_node reference when move_out_root is child of best_node; mitigated by retry loops + safety cap + caller validation; adversarial 80/80 pass. (4) regraft_violates_constraint DFS timestamps correct. (5) classify_clip_constraints correct. 902 constraint tests pass.

area: 4
reviewed_by: D
date: 2026-03-25
tier: n/a (pre-tier)
yield: 1 (consensus stability bug fixed)
notes: Parallelism & RNG. (1) BUG FIXED — consensus stability check ran every 200ms poll not per-replicate; unchanged counter incremented on idle polls → premature termination with slow replicates. Fix: track replicates_done at last check. (2) FRAGILITY — R_CheckUserInterrupt longjmp in try/catch is ABI-dependent (OK on Windows/SJLJ, fragile Linux/DWARF); needs R_UnwindProtect (R>=3.5). (3-6) DataSet/ConstraintData deep copies, thread-local RNG, pool mutex, stop_flag relaxed ordering all correct. (7) Rf_error in ts_wagner reachable from workers if n_tip<3 (practically unreachable). (8-9) dead code + redundant fuse_round cosmetic.

area: 3
reviewed_by: D
date: 2026-03-25
tier: n/a (pre-tier)
yield: (entry truncated in source — see 2026-03-27 F review for area 3)

area: 9
reviewed_by: C
date: 2026-03-20
tier: n/a (pre-tier)
yield: 1 (boundary-edge constraint bug fixed)
notes: Wagner & addition trees. (1) BUG FIXED — wagner_edge_violates_constraint + regraft_violates_constraint used is_ancestor_or_equal(cn, below) returning true when below==cn; for MUST_OUTSIDE this wrongly rejected the boundary edge. Fix: `&& below != cn`. Search-quality improvement. (2) 2 regression tests added. (3) clarified ts_random_wagner_tree can't guarantee constraint for all addition orders without posthoc (retry loop is the guarantee). (4) dead n_added/ew_score noted. 43+18+152 tests pass.

area: 8
reviewed_by: B
date: 2026-03-20
tier: n/a (pre-tier)
yield: 1 (LogCumSumExp NaN bug fixed)
notes: Test suite health + ParsSim log-space convolution. (1) BUG FIXED — `.LogCumSumExp()` in pp_info_extra_step.r produced NaN instead of -Inf when both accumulator and new value are -Inf (IEEE -Inf-(-Inf)=NaN). Reachable if MaddisonSlatkin returns NEG_INF for interior step. Fix: guard `if (is.finite(x[k]) || is.finite(Lk[k]))`. 7 assertions added. (2) active-range bounds verified. (3) .ApproxStepInformation MC branch not jointly normalized (noted, IC clamped to 0).

area: (Shiny module wiring)
reviewed_by: (B-era)
date: 2026-03 (pre-rotation snapshot)
tier: n/a (pre-tier)
yield: 0
notes: server.R/mod_*.R full review. (1) cb_ref 4 callbacks wired correctly after modules init. (2) cross-module updateXxxInput: only mod_data.R:203 (parent_session → treespace-relators), fragile but correct. (3) no orphaned observers. (4) isolate() in result observers correct. (5) progress polling observer correctly gated (invalidateLater only when gates pass). (6) ShowConfigs top-level DOM IDs correct. (7) UpdateActiveTrees reentrancy guard via on.exit. (8) change-detection pattern prevents reactive cascades. No bugs filed.

area: 5
reviewed_by: B
date: 2026-03-19
tier: n/a (pre-tier)
yield: 2 (inapp_state + MAX_STATES guard fixed)
notes: Data pipeline & simplification. (1) FIXED — build_reduced_dataset didn't copy ds.inapp_state (harmless now; would bite HSJ). Added. (2) FIXED — no guard for n_states>32; `(1u<<s)` UB. Added Rf_error. (3) NOT FILED — build_reduced_dataset omits HSJ/Sankoff fields (not used by sectors yet). (4) simplification correctness verified. (5) constraint column-major indexing verified. (6) EW offset interaction verified. 10 assertions added; 1679 ts-* pass.

area: 4
reviewed_by: E
date: 2026-03-19
tier: n/a (pre-tier)
yield: 0
notes: Parallelism & RNG. thread-local RNG set/cleared correctly; no R API from workers; pool mutex; atomic stop flag; seeds pre-generated; DataSet deep copy; HSJ/Sankoff stateless thread-safe. PERF NOTE: XFORM rebuilds SankoffData every score_tree (not filed). init.c 45 entries match. Score verification serial=parallel. No bugs.

area: 2 (complementary)
reviewed_by: A (complementary to B focus 2)
date: 2026-03-19
tier: n/a (pre-tier)
yield: 0 (1 found, not filed — conservative)
notes: Search topology — deeper state-restoration analysis. (1) FOUND — SPR stale scoring arrays after rejected regraft (nodes on regraft-to-root path retain regrafted values); CONSERVATIVE ONLY (screening, never acceptance; final always full_rescore). Test added test-ts-spr-state-restore.R (33 assertions). Not filed (SPR secondary, self-correcting). (2) TBR all reject paths correct. (3) NNI restore correct. (4) drift saved_postorder fix verified. (5) hash collision negligible.

area: 3
reviewed_by: B
date: 2026-03-19
tier: n/a (pre-tier)
yield: 1 (pattern_freq exponential-blowup fixed)
notes: Ratchet & perturbation. (1) BUG FIXED — perturb_upweight/perturb_mixed used `ds.pattern_freq[pat] *= 2` per char; shared pattern index → `original*2^N` blowup (integer overflow risk). Fix: `+= 1` (additive). (2) PerturbSnapshot save/restore correct. (3) upweight_mask has no effect on IW/profile (redundant but harmless). (4) sectorial reinsertion revert correct. (5) XSS sectors independent. (6) CSS via TBR no revert needed. (7) fuse tied-score correct. (8) fuse ancestor stale marking dead code. (9) RSS stale subtree sizes efficiency-only. 1404 tests pass.

area: 2
reviewed_by: B
date: 2026-03-19
tier: n/a (pre-tier)
yield: 1 (drift saved_postorder fixed)
notes: Search topology invariants. (1) BUG FIXED — drift_phase sets saved_postorder once, never updates after accepted moves; stale postorder restored on last-candidate reject → subsequent tbr_search full_rescore iterates stale postorder (wrong scores). Self-healing once a move accepted. Fix: build_postorder() before return when n_accepted>0. (2) TBR save/restore correct. (3) SPR build_postorder after every clip/unclip. (4) hash dedup negligible. (5) reroot path reversal correct. 1397 tests pass.

area: 1
reviewed_by: C+A
date: 2026-03-19
tier: n/a (pre-tier)
yield: 1 (upweight_mask in indirect length fixed)
notes: Fitch scoring correctness. **C:** (1) BUG FIXED — fitch_indirect_length_bounded/cached didn't account for upweight_mask, underscoring candidates during ratchet TBR; also nx_cost in ts_tbr/ts_search/ts_drift + drift RFD. (2-7) downpass_node, stop condition, dirty propagation, NA three-pass (Brazeau et al.), extract_char_steps, union-of-finals all correct. **A (complementary):** (A1) PROVED standard Fitch incremental uppass dirty-flag correct even when downpass stops before root. (A2) FOUND theoretical NA children_app staleness (conservative; full_rescore catches). (A3) confirmed C's fixes. 1397 tests pass.

area: 2
reviewed_by: C
date: 2026-03-25
tier: n/a (pre-tier)
yield: 1 (T-235 filed)
notes: Search topology invariants — thorough review of ts_tbr/ts_drift/ts_search/ts_tree. (1) BUG FOUND (T-235) — spr_search stale state arrays after rejected regraft (full_rescore overwrites all states; spr_unclip only restores clip-to-root incremental saves → stale divided_length baselines). Conservative (final always correct). Fix: full_rescore after spr_unclip on reject path. Low impact (presets disable sprFirst). (2) TBR correct pattern. (3-4) drift saved_postorder + subtree_sizes known/harmless. (5) states_valid dead code. (6) FNV hash negligible. (7) NNI restore correct. (8) constraint enforcement post-hoc check correct.

---

last_focus: 9
