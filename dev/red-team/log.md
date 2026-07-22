# Red-team round log — TreeSearch

Append-only record of every red-team round. **Newest first.** Each invocation of
`/red-team` adds one entry and updates `last_focus:` at the **bottom** of this file. The
next area is `(last_focus mod N) + 1`, where `N` is the current row count in
`focus-areas.md` (13 as of 2026-07-03 — **not** the stale `10` this line said until then,
which made areas 11-13 mathematically unreachable by normal rotation; see RT12-01,
2026-07-03 area-12 round below). Recompute `N` whenever a row is added.

**Entry format** (per round): `area`, `reviewed_by`, `date`, `tier`, `yield` (count of
*confirmed* findings filed), `notes`. Rounds before 2026-06 predate the Sonnet→Opus→Fable
tier system and are tagged `tier: n/a (pre-tier)`; the historical `reviewed_by` values
(`A`–`F`, `d1`–`d8`, `Claude (cpp-search)`) were ephemeral dispatcher agent IDs.

This log was migrated 2026-06-16 from `.positai/expertise/red-team.md` (the superset that
carried the full 2026-05-26 rotation) and extended with the 2026-06-15 CRAN run and the
2026-06-16 Shiny round. Durable lessons live in `../expertise/red-team.md`; open findings
in `findings.md`.

---

area: 13 (constraints — assigned deep-dive, NOT a rotation round; last_focus left at 11)
reviewed_by: opus orchestrator (harness author) + advisor (2 rounds)
date: 2026-07-04
tier: opus
yield: 1 (T-333 P3) — the hypothesised P1 (std::bad_alloc from `MaximizeParsimony(constraint=)`) **REFUTED**
notes: Single assigned question: is the T-327 constraint-repair guard
`build_postorder().size()==n_internal` (`impose_one_pass`, `ts_constraint.cpp:735`,
commit 6b60f235) equivalent to structural validity, or can a **net-zero** corruption
(a node double-referenced +1, a node orphaned −1) slip it and reach the advertised
`std::bad_alloc`? **Answer: the guard is NOT a full validator (the slip is real and
reachable) — but the P1 does not exist.** Built an exhaustive standalone C++ harness
(`dev/red-team/heavy-tests/impose_validity/`): enumerates ALL
(2n-3)!! rooted binary trees for n_tip=4..8 (counts verified 15/105/945/10395/135135),
runs the **verbatim** kernel fns (`topology_spr`, `collect_edges_in/outside_subtree`,
`find_maximal_subtrees`, `compute_node_tips`) extracted from `git show
HEAD:src/ts_constraint.cpp` at build time (cannot drift) + the **real**
`build_postorder`/`init_from_edge`, and checks guard-accept ⇒ an independent
full-validity oracle. Passes: **A** geometric superset (proves guard≠validator; ~1.33M
accept-on-invalid, but on triples `collect_edges` can't emit); **C1** faithful FIRST-move
emission (exact best_node/masks/move-roots + real targets); **D** faithful model of the
WHOLE `impose_constraint` loop (≤ n_splits+1 passes, branch over RNG-reachable targets,
continue past accepted-corrupt trees, model pass boundaries / moves==0 / bail).
**Findings:** the slip is real and **TYPE-1** (genuine left/right corruption = double-ref
+ orphan; one hand-verified end-to-end — `topology_spr`'s root-child degenerate case
grafts a clip onto its own parent edge, orphaning one node and doubling another),
reachable on the FIRST repair move at n≥8 (1890) and within/across passes at n≤6 (13642).
So the guard-comment "guarantees no corrupt tree ever reaches a DFS helper" is **factually
wrong** — they do reach `compute_node_tips`/`map_constraint_nodes`. **But no P1:**
(crash) **0** root-reachable left/right cycles AND **0** `parent[]` cycles across ALL
~1.33M accepted-invalid trees — `build_postorder`'s cap unconditionally rejects any
root-reachable cycle, and the parent-ascending consumers (`tbr_search`→`reroot_at_tip`
`ts_tbr.cpp:106-107`, run on the post-impose tree BEFORE verify at
`ts_nni_perturb.cpp:105` — the vector the advisor flagged) never hit a parent-cycle ⇒ no
`std::bad_alloc`; (wrong-answer) **0** corrupt FINAL trees survive the callers'
`map_constraint_nodes` verify-and-discard (all three sites `ts_driven.cpp:1056`,
`ts_nni_perturb.cpp:119`, `ts_parallel.cpp:92`) for single-split n≤6 within probe budget ⇒
nothing malformed is scored/returned. Correctness rests ENTIRELY on defense-in-depth
(unconditional cap + caller re-verify). **Filed T-333 (P3 hardening):** replace the `:735`
size-check with a full left/right arborescence check (in-degree-1 non-root + root DFS
visits all n_internal; O(n), same cost class) → self-sufficient guard + definitively
closes the two residuals not exhaustively covered here (multi-split constraints; n≥7
continuation of the full model). Advisor caught two things mid-round: (1) a
**witness-fabrication harness bug** — the verbatim `collect_edges` fns APPEND and never
clear their output; a hoisted/reused buffer accumulated stale targets → 7709 phantom
witnesses, and C1 dropped 123→0 at n≤5 once cleared (validates the verify step); (2) the
**parent-ascension crash vector** (added `has_parent_cycle` → 0). Harness durable
(`driver.cpp` + `extract_funcs.sh` + `build_and_run.sh` + `README.md`); `bash
build_and_run.sh` reproduces (n≤7 local seconds-to-minutes; n=8 Pass-D gated off).
**FIXED 2026-07-04 (same round).** Added `structurally_valid()` to
`src/ts_constraint.cpp` (in-degree-1 non-root + root-DFS-covers-all + parent[]-inverse,
O(n)); the `try_move` guard now validates-before-rebuild and reverts on any corruption.
**Perf-scoped deliberately:** the check lives in `impose_one_pass`/`try_move`
(constraint-repair, runs ONLY under `constraint=`), NOT in `TreeState::build_postorder`
(97 hot-path call sites, left untouched) — so unconstrained search pays nothing; within a
constrained run it's the same O(n) class as the `build_postorder` already there (advisor
concurred: pushing it into `build_postorder` "to cover all callers" would be the one real
regression, ruled out). Harness **re-pointed at the real guard**: `extract_funcs.sh` now
also extracts `structurally_valid` verbatim, and `probe()` asserts it equals the
independent `full_validity()` on every tree (new `g_guard_mismatch` gate, exit 3 if it
ever diverges — the durable CI signal against a future revert-to-size-check). Re-run:
**0 disagreements** over all `(2n-3)!!` trees n=4..6 (every one of the 13642 old-guard
witnesses now rejected). Package builds clean (`-Wall -pedantic`); R integration test
added (`test-ts-impose-constraint.R` "T-333: structural guard survives root-child repair",
49 assertions green — confirms no over-rejection/crash regression through the linked
engine). to-do.md T-333 → FIXED (pending CI); findings.md row removed. GHA CI (ASan + full
check) pending on push.

---

area: 13
reviewed_by: opus finder (a819d5b1) + opus orchestrator (current-HEAD re-validation)
date: 2026-07-02
tier: opus
yield: 1 (T-327 P2) — after correcting for a STALE review build + finder P1→P2
notes: Constrained search correctness — NEVER reviewed (override of rotation area 12/meta, same CRAN-correctness rationale as area 11). **CRITICAL PROCESS NOTE: the finder's harness (rt-collapse @ d25a0f6c, 2026-06-24) was ~8 days STALE** — missing `d9a4f827` (T-213 constraint fix) and `0daea13f` (collapse). Findings were therefore re-validated against a FRESH current-HEAD build (rt-current @ 91918d1a). See memory [[redteam-verify-against-current-tip]]. **T-327 (P2) FILED** — `impose_one_pass` relocates the captured `best_node` via `topology_spr` mid-repair → corrupt/cyclic tree → unbounded `build_postorder` → catchable `std::bad_alloc`; `MaximizeParsimony(constraint + nniPerturbCycles>0)` fails on valid input. CONFIRMED on current HEAD (`impose_one_pass` byte-identical to the stale build; untouched by the fixes). Severity corrected P1→P2 (catchable error, not SIGSEGV; opt-in trigger). Confirms the long-unconfirmed stale-`best_node` hypothesis [[impose-constraint-verify-gap]] (600-seed stress missed it; real trigger = constraint + nniPerturbCycles>0). **RT13-02 (nni-perturb verify-before-capture) NOT FILED — already fixed by `d9a4f827`** (the finder re-discovered an already-fixed bug because of the stale build). **Assessed, not filed:** RT13-03 (`impose_constraint` ~40% wrong) = same root cause as T-327, a fuzzer prevalence stat on random trees×splits, not a separate user-facing bug; RT13-04 (silent constraint violation) = latent, finder could NOT force end-to-end (masked: after `update_constraint` sets `cn=-1`, `regraft_violates_constraint` rejects all regrafts); RT13-05 (bail / pass-exhaustion) folds into T-327; RT13-06 (Wagner `AdditionTree` root-edge constraint fallback, warning-not-error) overlaps filed T-324. **Ruled out / clean (current-code trace):** TBR per-clip constraint gating + resync on accept/equal/reject/tabu; fuse & parallel-fuse DO verify-before-capture (`fused_ok` discard); sector resync. Harness: rt-current worktree @ 91918d1a + `.rtlib`; repros `scratchpad/repro_loop_current.R`, `repro_seed4_raw.R`, `fuzzer.cpp`. **Open follow-up:** can the `impose_constraint` corruption crash the DEFAULT fuse path inside `build_postorder`/`map_constraint_nodes` *before* the `fused_ok` discard? (finder: not in 60 seeds). Filed to findings.md + to-do.md (T-327, P2, [m:opus e:high]).

---

area: 11
reviewed_by: opus finder (a5f8e723) + opus orchestrator (verification)
date: 2026-07-02
tier: opus
yield: 2 (T-325 P2, T-326 P3) — finder over-claim corrected during verification
notes: Zero-length collapse / MPT set — **NEVER reviewed**, default-ON since 2026-06-24 (override of rotation area 10, which overlaps the perf agents' hot `ts_fitch.cpp`). **T-325 (P2) FILED.** The same `TreePool` in run_driven_search is collapse-deduped on the main-loop path (`ts_driven.cpp:929` `add_collapsed`, collapsed-split keys) but full-keyed on the MPT-enumeration path (`ts_tbr.cpp:2345/2377` `collect_pool->add`; collapsed skipped at `:1383` when `collect_pool` set); `splits_equal` can't cross-match key types → the returned pool holds trees differing only in zero-length resolutions. `MaximizeParsimony(pd)` at **full default returns 30** where collapse-dedup should return ~1 (strict consensus = 2 supported splits → all collapse to one topology). **Finder MAGNITUDE REFUTED during verification (this is the value of the verify step):** the finder's `fullKey` over-merged RF=2 trees, so its "6 distinct / 78% exact-duplicate / 4.7×" is a measurement artifact; `RF.dist`+`TreeDist` show 27 distinct full topologies + 1 exact dup. Real defect = over-return vs the COLLAPSED baseline (≈1 vs ≈30), mechanism code-proven, reproduces at no-arg default. **T-326 (P3)** test gap: test-ts-collapsed.R:160 asserts only `pool_size>=1`, never that the returned set is collapse-deduped. **Seams RULED OUT:** C over/under-collapse (finder: 150-trial brute-force MPR oracle, 0 discrepancies; aggressive path also default-OFF `TS_COLLAPSE_AGGRESSIVE`); B hash-collision (`is_duplicate` does `splits_equal` fallback — safe). **Seams STILL OPEN (yielding, next visit stays opus):** A/F rooting-invariance near root under-tested (near-root-child geometry untested); D `local_cost`/`prelim` staleness at aggressive-flag call sites (interrupt/fuse paths unverified); G constraint-violation-via-collapse (not investigated). Harness preserved: rt-collapse worktree @ d25a0f6c + `.rtlib`; repros in scratchpad (`quantify.R`, `decide_key.R`, `mixedkey3.R`). Fix = post-enum collapse-dedup pass (maintainer call; NOT inline — perf agents active on cpp-search).

**POST-HOC CORRECTION (2026-07-02, later same day):** T-325/T-326 were **already fixed 8 days before this round ran**, by `0daea13f` (2026-06-25, "collapse zero-length branches by default; consistent MPT dedup") — the exact fix this round proposed (R-level `collapse=TRUE` default via `ts_collapse_pool`, plus the C++ enum `add_collapsed` change). The finder/orchestrator cited `ts_tbr.cpp:2345,2377`, which are `0daea13f`'s **pre-fix** line numbers (its diff hunk headers read `@@ -2342,7 +2361` / `@@ -2374,7 +2406`) — the trace read a pre-fix version of the file even though HEAD had carried the fix for over a week, i.e. a stale checkout/lib, not stale source (see memory `stale-local-treesearch-lib`). Empirically confirmed: fresh tarball build of `0daea13f^` reproduces exactly 30 trees / `n_topologies=30` on the round's own repro; fresh build of current HEAD returns 1. Isolated which layer is decisive: with the C++ enum fix disabled (`TS_ENUM_RESOLVED=1`) the result is still 1 — the R-level `collapse=TRUE` default alone fixes the user-visible symptom (the C++ enum change fixes a separate internal-pool self-consistency issue). T-326's ask is already covered by `test-MaximizeParsimony-features.R:402-421` ("collapse = TRUE contracts a soft polytomy to one collapsed tree"), added in the same commit; re-ran it green on a fresh build. Rows removed from `findings.md` and `to-do.md`; closure recorded in `completed-tasks.md`. Residual: `test-ts-collapsed.R:159-179` still only asserts `pool_size >= 1` (redundant-weak now, not dangerous-weak — the invariant is covered elsewhere) — not worth a dedicated task. **Lesson for next rotation:** before filing a finding, check whether the cited line numbers still match current `git show HEAD:<file>` — a mismatch is the tell that the build/lib under test is stale, not the source.

area: 10 (signal-resolution only — NOT a full finder round; last_focus left at 9)
reviewed_by: Claude (cpp-search, /red-team orchestrator diagnostic)
date: 2026-06-16
tier: opus (orchestrator; finder not spent — signal dissolved by direct diagnostic)
yield: 0 — the area-9 high-sev signal REFUTED with positive correctness evidence
notes: User invoked `/red-team` to pick up the area-9 high-severity signal (kernel
NA+IW score ≠ `TreeLength()` on Vinther2008). **Root-caused and REFUTED directly,
without spending a finder** (best issues-per-token outcome). The area-9 finder
mis-diagnosed it as "a different per-character minimum for inapplicable chars." It is
actually **plain-IW vs XPIWE**: `TreeLength()` defaults to `extended_iw = TRUE`
(XPIWE — Goloboff 2014 Extension-3 missing-data correction, `R/tree_length.R:144-156`:
`f = 1 + r·(nTaxa−obs)/obs`, `eff_k = k/f`, `phi = (1+eff_k)/(1+k)`,
`fit = h/(h+eff_k)`, `Σ fit·w·phi`), whereas the area-9 cross-check called the kernel
via `ts_fitch_score(..., min_steps, concavity)` — **plain IW**, no XPIWE args. Two
different objectives by construction; the "non-rational" 0.16573 the finder flagged is
the `phi`/`eff_k` scaling (`eff_k≈9.3`), not a min-steps bug. **Three diagnostics
(installed pkg, Vinther2008) settle it:** (1) `TreeLength(extended_iw=FALSE)` ==
kernel `ts_fitch_score` plain IW, **exactly** (8.566683 == 8.566683). (2)
`TreeLength(extended_iw=TRUE)` differs (8.371583) — the XPIWE correction, by design.
(3) **kernel XPIWE == R XPIWE**: a real `MaximizeParsimony(concavity=10)` run's
`attr(.,'score')` (XPIWE, since production sets `useXpiwe <- isTRUE(extended_iw) &&
is.finite(concavity) && !useProfile`, `R/MaximizeParsimony.R:813-839`) == an
independent `TreeLength()` XPIWE rescore of the returned tree, **exactly** (1.521827).
So production **optimises and reports the same objective** in BOTH plain-IW and XPIWE
modes — no optimise-vs-report mismatch. (Numbers are tree-dependent; the three
equalities are not.) **Conclusion: NOT A BUG.** The signal raised in the area-9 entry
below is resolved; do not carry it forward. **last_focus deliberately left at 9** — this
was a targeted signal check, not a rotation round, so area 10 still earns a proper
finder sweep next `/red-team`. **Residual (un-spent) area-10 surface for that future
round:** profile delta capping (7cff7870), `e/(k+e)` delta, `precomputed_steps` offset,
`info_amounts` capping, `concavity=1.0` profile sentinel, DAT-002 `obs==0` XPIWE
division reachability, and the OPEN clipped-subtree IW-screening follow-up. With fable
returning, that round is the natural escalation target (area-10 opus-class seam was dry
2026-05-19 + 2026-05-26; XPIWE consistency now positively verified, so what remains is
fable-class subtlety).

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
`feature/hsj-absent-state`), **not** a normal single-area rotation round.
The session token
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
(resample-interrupt = XPIWE division, IW clipped-subtree screening, CRAN
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

area: 9 (Wagner & addition trees)
reviewed_by: opus finder a67c314099ee52908 + orchestrator source-trace verify (opus)
date: 2026-07-02
tier: opus (override — natural next was area 1 Fitch; avoided src/ts_fitch.cpp = perf-agents' hot file → collision + stale-build risk after the T-325/T-326 phantom)
yield: 1 (T-328 filed, P2)
notes: Wagner start-tree construction (ts_wagner.cpp/.h, ts_simplify.cpp, ts_rcpp.cpp boundary). Opus finder did STATIC-READ only (final message stalled twice mid-stream — API error — recovered via SendMessage); orchestrator verified all 4 candidates by independent source trace, git-confirmed at HEAD c74ee6e6 (ts_simplify.cpp/ts_data.cpp unchanged since harness 91918d1a, tree==HEAD, no in-flight edits). (1) WGN-9-A → FILED T-328 (P2, P1-tail). ts_simplify.cpp:194 computes `token = tip_data - 1` with no range check, then `token_states[token]` (:195, :198-199) unguarded; token_states sized n_tokens (ts_data.cpp:45) → OOB (value 0 → token -1, silent misread; value > n_tokens → likely SIGSEGV). make_dataset guards vector LENGTHS / n_states, NOT tip_data VALUES → distinct sibling of T-323 (fix in one PR: harden all Wagner Rcpp-boundary indices, values as well as lengths). Finder OVERCLAIMED public reachability; corrected DOWN to internal-only (AdditionTree.R:95-98 builds tip_data from a validated phyDat) → P2 with a P1 tail (a public transform dropping a contrast row without reindexing — PrepareDataProfile / PrepareDataIW / .Recompress — would inject it on a public path; follow-up, NOT hunted this round). Logic-traced, no repro; severity by the T-323 precedent. (2) WGN-9-B (P3 perf) NOT FILED — wagner_map_constraint_nodes rebuilds full postorder per insertion (ts_wagner.cpp:253, called :561) → ~quadratic constraint-build tax; a PERF lever (perf agents' turf — steer clear), logged for them, no spawn_task. (3) WGN-9-C (P3 softmax) NOT FILED — finder self-retracted (tiny-positive temperature → exp overflow → chosen<0 fallback :715-719; judged an acceptable numerical fallback). (4) WGN-9-D (P3 robustness) NOT FILED, area-13 lead — random_constrained_tree (ts_wagner.cpp:984, via ts_driven.cpp:116) assumes laminar splits with no posthoc check, BUT overlapping/non-laminar splits ⟺ INCOMPATIBLE ⟺ an impossible constraint (.PrepareConstraint MaximizeParsimony.R:64-125 builds one split per constraint char; a satisfiable constraint is always laminar) → invalid-input-only (should error on an impossible constraint), not a valid-input correctness bug. NOTES-FOR-NEXT-REVIEWER (seam STILL YIELDING → next visit stays opus): the LIVE thread is wagner_incremental_rescore (ts_wagner.cpp:108) — finder flagged it "hardest logic, needs a harness", static-read only this round, NOT stress-tested; that seam is open.

area: 13 (Constrained search correctness)
reviewed_by: opus finder a2feb110 + opus verifier abfbba47
date: 2026-07-02
tier: opus (never formally rotated; dispatch parallel wave, out-of-rotation alongside area 11)
yield: 1 (T-329 filed, P2; high-sev OOB WGN-13-C CONFIRMED by code-trace after advisor caught the verifier's incomplete proof)
notes: Opus finder spent its whole budget on the WGN-9-D lead (impossible-constraint handling in random_constrained_tree) and did NOT reach the fuse / parallel-fuse / sector posthoc-retry / TBR clip-gating verify-before-capture questions or map_constraint_nodes resync on reject paths — those remain OPEN for a follow-up opus pass. Finder + verifier both stalled mid-emit (API), recovered via SendMessage. (1) WGN-13-C — HIGH-SEV OOB → CONFIRMED by orchestrator code-trace (the P2 driver). The opus verifier REFUTED it via an invariant proof (a binary tree on n_tip leaves has exactly n_tip-1 internal nodes → consumption laminarity-invariant) — but the ADVISOR flagged that the proof fumbled its own worked example (finder: losing split collapses to {4} = 1 tip; verifier: 0 tips — the exact 0-vs-1 boundary that decides the OOB), so I re-traced. The empty-group case IS guarded (ts_wagner.cpp:1076-1078) and the root-level collection guards split_root>=0 (:1093), BUT the child-split collection at :1072 (items.push_back(split_root[j])) is NOT guarded → a collapsed split (split_root=-1) with a strict-superset parent injects -1 as a phantom item into the parent's resolve_randomly → tree.parent[-1] OOB write (before start) + phantom leaf → +1 internal node → next_internal hits 2*n_tip-1 → tree.left/right[n_tip-1] OOB past end (arrays n_tip-1 / 2*n_tip-1, init_wagner_state:25-32). Reachable via a >=3-way-overlap impossible constraint; same class as T-323/T-328. Not yet ASan-reproduced (local ASan blocked) → GHA-ASan follow-up. LESSON: never lock a high-sev REFUTE whose author got its own key worked example wrong. (2) WGN-13-A → FILED T-329 (P2, driven by C). Impossible (non-laminar) constraint silently best-efforted by random_constrained_tree (ts_wagner.cpp:1038-1040 tightest-split tip-ownership → losing split's clade collapses, no posthoc warning) → returns a constraint-violating tree. (3) WGN-13-B (regraft_violates_constraint cn<0 silent all-reject freeze, ts_constraint.cpp:358) CONFIRMED but as a DOWNSTREAM SYMPTOM of A (only from impossible input; valid constraints keep constraint_node[s] >= 0) → folded into T-329, not filed separately. Fix: (a) validate constraint laminarity/satisfiability (pairwise split compatibility) in .PrepareConstraint and stop() cleanly (also pre-empts the impossible-constraint subset of T-324); (b) defensive — add the split_root[j]>=0 guard at ts_wagner.cpp:1072 to mirror :1093. Seam status: STILL YIELDING — the 4 core verify-before-capture questions (fuse/parallel-fuse/Wagner-retry/sector) + reject-path map_constraint_nodes resync are UNTOUCHED; next visit stays opus and should start there.

area: 11 (Zero-length-branch collapse / MPT set)
reviewed_by: opus finder aa94e05a (built isolated worktree, empirical repro) + orchestrator source-verify (opus)
date: 2026-07-02
tier: opus (NEVER REVIEWED before; dispatch parallel wave, out-of-rotation alongside area 13)
yield: 3 (T-330 P1, T-331 P3, T-332 P3)
notes: First-ever review of this never-reviewed, default-ON-since-2026-06-24 area. Opus finder did DEEP empirical work — built an isolated worktree @ HEAD c74ee6e6 (CCACHE_DISABLE=1) and reproduced bugs end-to-end (119 tool-uses, ~2h; emit stalled repeatedly — it DECLINED the write-to-file recovery as a suspicious "coordinator" instruction, a reasonable safety false-positive, and eventually emitted in chat). Orchestrator independently source-verified the P1. (1) RT11-01 → FILED T-330 (P1). compute_collapsed_flags (ts_collapsed.cpp:40-106) + _aggressive (:108-209) read ONLY ds.blocks[]; has_na fallback but NO HSJ/XFORM handling → blind to ds.hierarchy_blocks / ds.sankoff_*. Over-collapses HSJ/XFORM-supported clades. TWO manifestations: (a) final ts_collapse_pool over-collapse — finder reproduced: an HSJ-supported clade destroyed, Nnode 5→3, output identical w/ and w/o hsjConfig; (b) MPT-enum dedup at ts_tbr.cpp:2374/2419 has NO scoring_mode guard (ADVISOR caught this — the rescore paths ARE gated at :1348-1352 but the collapse-flag calls are not) → MPT set CORRUPTED during search, not display-only. Fix must span ALL call sites; fall-back-to-conservative is INSUFFICIENT (the conservative fn is equally blind). Blast radius HSJ/XFORM only (prioritization, not severity — advisor: don't launder blast-radius into severity; it's the normal path for those users). (2) RT11-03 → FILED T-331 (P3). total_words==0 early return (ts_collapsed.cpp:122-123) → fully-uninformative data returns binary trees + inflated n_topologies instead of the star (reproduced). Degenerate input; advisor: total_words==0 is a cliff not a gradient → dropped the finder's "muted realistic version" speculation. (3) RT11-02 → FILED T-332 (P3). ts_collapse_pool hangs on a non-binary edge matrix (TreeSearch:::-reachable via RcppExports.R:222, NOT the public API) → same internal-boundary class as T-323/T-328 (advisor: file for precedent-consistency, not just a note). (4) RT11-04 (TreePool tie float-equality) NOT FILED — finder self-ruled-out on the mainline path, EW-immune. RULED OUT this round: IW/profile-vs-EW flag consistency (same Fitch mechanics, concavity is a downstream reweight), constraint-canonicalization rooting. Seam status: STILL YIELDING (P1 found → next visit stays opus). NOT reached: multifurcating TreeLength (issue #259); dedup over/under-merge re-hunt (prior T-325 territory).

area: 12 (Red-team process meta-review)
reviewed_by: sonnet finder (ace5e16b)
date: 2026-07-03
tier: sonnet
yield: 1 (RT12-01, process defect, fixed inline — not a to-do row; process fixes are applied directly, not dispatched)
notes: First-ever review of the meta-review area. **RT12-01 (high-conf) FOUND AND FIXED INLINE:** the rotation-formula description hardcoded `(last_focus mod 10) + 1` at log.md:5, but the table has grown to 13 rows — areas 11/12/13 were mathematically unreachable by normal rotation (every visit to them in the log was an explicit manual override; the area-9 and area-13 2026-07-02 entries even disagreed with each other about the "natural next" area, a direct symptom of the broken pointer). Fixed: log.md:5 now reads `(last_focus mod N) + 1` where N = current row count, with a note to recompute N when a row is added. **Other findings applied directly (documentation-only, no code-behaviour claim needing verification):** added `ts_tree.cpp/.h`+`ts_pool.cpp/.h` to area 2 (TreePool/tree-op seams were unowned despite being central to T-325/T-330's MPT-set mechanics); added `ts_ls.h/.cpp` to area 5 (LeastSquares kernel had a filed P1 — 06-15 CRAN sweep, RSS=0 on degenerate `dist` — but no owning area to catch a regression); renamed area 10 to "Alternative scoring kernels: Profile/IW/HSJ/XFORM" and added `ts_hsj.cpp/.h`+`ts_sankoff.cpp/.h` (exactly the kernels T-330 proved are a blind spot, previously unnamed anywhere). **Not applied (judgment calls, logged for awareness only):** RT12-02 (area 3 ratchet — 3 consecutive dry opus rounds vs. a "keeps recurring" rationale note; recommend one more round before downtiering, not an immediate action); RT12-03 (area 4 RNG — dry-round count is misleading, off-rotation T-309 already proves it's live; rationale already correctly overrides, no action needed); RT12-04 (area 6 R↔C++ — 2 consecutive dry sonnet rounds, area's own rationale already says "escalate if a sonnet pass goes dry" — standard skill logic already escalates next visit automatically, no doc change needed); RT12-06 (area 13 may be too broad for one round's budget — recommend future rounds split key-questions into sub-passes, e.g. "verify-before-capture audit" vs "impossible-constraint audit," rather than merging/splitting the area itself); RT12-07 (area 11's written rationale predates T-330/331/332 and is mildly stale — cosmetic only, left as-is). **Not yet checked (flagged for a future area-12 round):** `R/*.R` file coverage (e.g. `R/CharacterHierarchy.R`, HSJ weight-zeroing, directly implicated in T-330, is unnamed in any area) — this round only diffed `src/*` against area file lists. Seam status: still yielding (first-ever pass, found a structural rotation bug) — next visit stays sonnet.

area: 13 (Constrained search correctness)
reviewed_by: opus finder a5cbecdd + orchestrator source-trace verify (opus) + advisor high-sev cross-check
date: 2026-07-03
tier: opus (rotation: (12 mod 13)+1 = 13; continuing the seam the 2026-07-02 round left open)
yield: 1 (T-13-A → MERGED into existing T-324, augmented — not a new ID; verify-before-capture audit otherwise clean)
notes: Started exactly where the 2026-07-02 round deferred — verify-before-capture on EVERY impose_constraint() caller. **AUDIT RESULT (the core deliverable):** fuse (ts_driven.cpp:1042-1058) OK — maps constraint nodes, imposes, re-verifies, gates the pool add on `fused_ok` (orchestrator-confirmed by read). parallel-fuse (ts_parallel.cpp:84-94) OK — same map+re-check pattern (orchestrator-confirmed). sector (ts_sector.cpp:1367-1372, 1536-1541) OK — revert-on-violation (restore_clade+continue), never captures a violating tree (finder-reported, not independently re-traced). **Wagner build-retry NOT OK → the finding.** (1) **T-13-A is NOT a new ID — it re-discovers + DEEPENS existing T-324** (2026-06-16, area 9), so per "avoid re-reporting" I AUGMENTED T-324 rather than file T-333. Novel contributions beyond T-324's original `AdditionTree()`-scoped, warning-parity framing: (a) the `MaximizeParsimony()` per-replicate pool capture at ts_driven.cpp:929 is **ungated** — asymmetric to the fuse gate 100 lines below; (b) **confirmed NO downstream filter** (ts_rcpp.cpp post-:1390 + MaximizeParsimony.R post-search = collapse-protection only, :1028-1034) → a violating start reaches the user unflagged; (c) a violating start is NOT repaired by constrained TBR (regraft_violates_constraint returns true for all moves once constraint_node[s]<0, ts_constraint.cpp:354-360 → freeze), only conditionally by nni_perturb's impose_constraint (nni_perturb_per>0 + heuristic success); (d) **fix caveat (advisor-caught mis-patch trap):** the :929 gate must use `violates_constraint_posthoc`, NOT the fuse-style `constraint_node[s]<0` check, because a posthoc-only violation (all cn>=0 but fails full-Fitch) would slip a constraint_node gate. **Severity: P3-proven (missing gate + missing warning), escalates to P2 IFF reachability confirmed** — i.e. a satisfiable USER constraint whose violation survives all 100 independent reshuffles (has_posthoc=true only on user constraints; ts_rcpp.cpp:1390→build_constraint). The retry loop's existence proves pass-construction/fail-posthoc trees exist; open bit = 100-reshuffle persistence → **Hamilton hard-but-satisfiable-constraint probe recommended** (heavy compute, not local). Also updated T-324's STALE line numbers (767-780/731-737/554 → 745-754/784-797/571-576) to HEAD 4b833e7f. (2) **T-13-B RULED DOWN, not filed** — the parallel-fuse `violates_constraint_posthoc` short-circuit (returns false when !has_posthoc, ts_constraint.cpp:388) only bites an `auto_cd` constraint (build_constraint_from_bitsets, has_posthoc=false), and auto_cd (consensus-tightening heuristic) engages ONLY when there is no user constraint (use_auto_constraint = consensus_constrain && (!cd||!cd->active), ts_driven.cpp:724, and consensus_constrain defaults false). A tree "violating" auto_cd is search guidance escaping consensus, not a user-facing correctness bug. Advisor concurred. (3) **HIGH-SEV SIGNAL — escalates area 13 next round; SHARPENED this session from a vague lead to a precise, resolvable question (do NOT record CLEAN KILL).** The impose_one_pass stale-best_node concern is the machinery of ALREADY-FIXED T-327 (6b60f235: `reanchor_best_node` before every move + snapshot-validate-revert gated on `build_postorder().postorder.size()==n_internal`, ts_constraint.cpp:702-743). A THIRD dedicated opus finder this session emit-stalled TWICE (channel = binding constraint, per [[dispatch-gotchas]]), so the orchestrator traced it directly (opus) + advisor cross-check. **Refined question:** can the `postorder.size()!=n_internal` revert-guard (build_postorder, ts_tree.cpp:75-112 — DFS over left/right with NO visited-set, cap `preorder.size()>n_internal` at :103) be SLIPPED by a stale-`M` topology_spr corruption netting to EXACTLY n_internal (duplicate +k offset by orphan −k)? **Advisor-corrected invariant (my FIRST trace had a flawed mechanism — same trap as [[redteam-reverify-flawed-refutes]], caught before locking; I wrongly called Case A a double-parent/over-count when `above==parent(below)` re-parents `below` → floating cycle/under-count):** topology_spr preserves the child-slot BIJECTION (every non-root node in-degree 1) in NON-DEGENERATE position → any corruption is a floating rho-component → UNDER-count → caught. The bijection can ONLY break for DEGENERATE slot-collision graft targets (above/below coinciding with nx/nz/ns); one case checked (graft onto edge (nx,ns) → ns in-degree 2, nx orphaned — a duplicate+orphan, but ns's 2nd parent IS the orphan so still under-counts → caught). **Net-zero slip (a duplicate DOUBLY-reachable-from-root + a compensating orphan) remains UNPROVEN.** VERDICT: guard robust in the bijection-preserving regime; residual P1 risk confined to the degenerate slot-collision regime; NOT retired. **NEXT VISIT: NOT another finder — a BOUNDED EXHAUSTIVE HARNESS (mcmc-diagnostician / heavy-test under dev/red-team/heavy-tests/): for n_tip 4–8, enumerate trees × the exact (clip,above,below) triples impose_one_pass can emit, apply topology_spr, assert FULL validity (in-degree-1 + root-reachable + acyclic) whenever `try_move` returns true. Accept-on-invalid = confirming P1 repro; exhaustive small-n silence = strong kill. Directly tests whether `postorder.size()==n_internal ⟺ validity`, which IS the whole question — stop hand-enumerating (error-prone, already mis-traced once).** (4) NOT reached this round: Q3 (clip-gating FALSE-NEGATIVE — does regraft_violates_constraint ever ALLOW a violating regraft) and Q4 (laminar/nested-split consistency across TBR/Wagner/sector paths) — both untouched, open. Finder + orchestrator both traced against HEAD 4b833e7f (area-13 source unchanged since c74ee6e6; no in-flight edits). Finder's final emit stalled mid-stream (198 tokens) — recovered via SendMessage compact re-emit (76k tokens). Seam status: STILL YIELDING (deepened T-324 + open high-sev signal, now precisely characterized) → next area-13 visit = a BOUNDED VALIDITY HARNESS on the topology_spr / build_postorder-guard equivalence (see (3)), NOT another finder; the finder-shaped questions Q3 (clip-gating false-negative) and Q4 (laminar consistency) remain for a later opus finder round.

last_focus: 13
