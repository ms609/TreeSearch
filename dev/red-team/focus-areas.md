# Red-team focus areas — TreeSearch

Rotation table for the `/red-team` skill. Built once, edited rarely. Each `/red-team`
invocation reviews **one** area (the next in rotation, see `last_focus:` at the bottom of
`log.md`) at its earned tier, then records the round in `log.md`. Verified non-trivial
findings are filed as GitHub issues in `agent-issues/TreeSearch` (labelled `red-team`,
`sev:high|med|low`, `area:N`). Durable lessons (bug patterns, fragile areas) live in
`../expertise/red-team.md`.

## start_tier

`start_tier` is the tier a **never-visited or freshly-rotated** area starts at. Unlike the
skill's default (everything `sonnet`), these tiers **encode measured maturity** from the
ported round history (`log.md`): areas whose seams have only ever yielded subtle,
opus-class bugs start higher; immature seams that still bleed cheap bugs start at `sonnet`.
The rotation still adjusts per recorded yield — a dry round escalates one tier, a yielding
round re-visits at the same tier with a fresh agent, a high-severity signal escalates
immediately. Treat these as the starting point, not a ceiling.

**A tier is a rung, not a model.** A dry verdict is scoped to the *version* that produced it,
and a version bump at the same rung is a cheaper step than a rung bump — so `opus-4.8 dry` goes
to **Opus 5**, not to fable. The alias→version mapping lives in the model-version legend at the
top of `log.md`; seams that a version bump has made re-eligible are queued in
`escalation-backlog.md`.

| # | Area | Files | start_tier | Key questions |
|---|------|-------|-----------|---------------|
| 1 | **Fitch scoring correctness** | `src/ts_fitch.h/.cpp`, `src/ts_fitch_na.h`, `src/ts_fitch_na_incr.h`, `src/ts_fitch_na_dirty.h`, `src/ts_simd.h` (added 2026-08-04 — the bit-parallel SIMD portability layer every Fitch combine call goes through; UNMEASURED, no inherited maturity) | **opus** | Does incremental / dirty-set scoring match full `score_tree()`? Bounded variants bail correctly? NA three-pass edge cases? Write a targeted test if you find a gap. |
| 2 | **Search topology invariants** | `src/ts_tbr.cpp`, `src/ts_drift.cpp`, `src/ts_search.cpp`, `src/ts_tree.cpp/.h`, `src/ts_pool.cpp/.h`, `src/ts_tabu.h` (added 2026-08-04 — the tabu-list hash buffer that directly implements this row's own "symmetry-breaking hash collisions" question; UNMEASURED) | **opus** | After every rejected move, is topology fully restored? Undo stack correct? No stale `postorder`? Constraint metadata re-synced on *all* reject paths (incl. tabu)? Symmetry-breaking hash collisions? `TreePool` (dedup key, capacity eviction) consistent with the topology invariants above? |
| 3 | **Ratchet & perturbation** | `src/ts_ratchet.cpp`, `src/ts_sector.cpp`, `src/ts_fuse.cpp`, `src/ts_prune_reinsert.cpp` | **opus** | `active_mask`/`upweight_mask`/`flat_blocks` fully restored after perturbation? Sectorial reinsertion reverts on worse score? `build_reduced_dataset` copies all needed fields? Fuse handles tied scores? |
| 4 | **Parallelism & RNG** | `src/ts_parallel.cpp`, `src/ts_rng.h/.cpp`, `src/ts_driven.cpp`, `src/ts_resample.cpp`, `src/ts_heartbeat.cpp/.h`, `src/build_postorder.h`, `src/ts_strategy.h`, `R/Resample.R`, `R/Jackknife.R` (added 2026-08-04 — `ts_heartbeat`/`build_postorder.h` are exactly this row's bug class: main-thread-only R-API + `set.seed()`-reproducibility RNG state; `ts_strategy.h` is the bandit consumed by `ts_driven.cpp`; `R/Resample.R`/`R/Jackknife.R` are the R-level entries to the parallel resample path already named in T-336/T-337/T-398 but never an owned file. ALL UNMEASURED, no inherited maturity) | **opus** | Thread-local RNG set before any search call? **No R API (incl. `unif_rand`/`Get/PutRNGstate`) from worker threads** — note the resample path. Pool mutex correct? Atomic stop-flag races? Seeds drawn from R RNG before spawn? |
| 5 | **Data pipeline & simplification** | `src/ts_data.h/.cpp`, `src/ts_simplify.h/.cpp`, `src/ts_ls.h/.cpp`, `R/tree_length.R`, `R/IWScore.R`, `R/LeastSquares.R`, `R/PrepareData.R`, `R/data.R`, `R/data_manipulation.R`, `R/fractional-weights.R`, `R/length_range.R` (added 2026-08-04 — the R-layer half of `TreeLength`/`MinimumLength`/`CharacterLength`/IW scoring flagged unowned by the 2026-07-03 area-12 round and by `escalation-backlog.md` item 5(a); `R/tree_length.R` is directly implicated by open issue #16/T-400, sev:high. ALL UNMEASURED, no inherited maturity) | **opus** | `build_dataset` handles edge cases (all-ambiguous, single-state, zero-weight, `n_states==32` UBSAN)? `build_reduced_dataset` copies all fields? XPIWE `obs==0` division? Least-squares distance fitting (`ts_ls.cpp`) — degenerate `dist` (NA/Inf) handled (cf. filed P1: `LeastSquaresFit`/`LeastSquaresTree` RSS=0 garbage)? |
| 6 | **R ↔ C++ interface** | `src/ts_rcpp.cpp`, `src/TreeSearch-init.c`, `src/RcppExports.cpp` (added 2026-08-04 — the generated third leg of the same boundary triangle as `R/RcppExports.R`; a signature drift here is exactly what `compile-attrs.R` exists to catch), `R/RcppExports.R`, `R/MaximizeParsimony.R`, `R/SearchControl.R`, `R/ts-driven-compat.R` (added 2026-08-04 — the flat-argument-to-grouped-list compatibility shim T-398/#14 traced through; UNMEASURED) | **sonnet** | Arg counts match? Concavity sentinel translated? Edge-matrix conventions? Return value attributes/types set (frozen-API `logical` vs `integer`)? Parameter validation in R layer? |
| 7 | **Shiny module wiring** | `inst/Parsimony/server.R`, `inst/Parsimony/global.R`, `inst/Parsimony/ui.R`, `inst/Parsimony/server/mod_*.R`, `inst/Parsimony/server/app_state.R`, `inst/Parsimony/server/logging.R`, `inst/Parsimony/tests/testthat.R`, `inst/Parsimony/tests/testthat/*.R` (added 2026-08-04 — the app's own test suite, in neither this row nor area 8's `tests/testthat/` glob; includes two quarantined `_problems/test-app-smoke-*.R` worth reading for *why* they're quarantined. UNMEASURED) | **sonnet** | Forward-ref callbacks resolve? Cross-module `updateXxxInput` namespaces correct? Re-entrancy / double-launch guards? Stale dataset-hash on async tasks? `onStop` cleanup (cancel signal + temp files)? Orphaned observers? |
| 8 | **Test suite health** | `tests/testthat/*.R` (broadened 2026-08-04 — the literal `test-ts-*.R` glob excluded ~44 of ~110 files, e.g. `test-CustomSearch.R`, `test-tree_length.R`, `test-MaddisonSlatkin.R`, `test-Concordance.R`; issue #4/T-363 already treated a non-`ts`-prefixed file as in-area, so the row's own precedent was already broader than its text. The non-`ts` files carry NO inherited maturity — treat as never-reviewed) | **sonnet** | Tier guards correct? Vacuous (always-pass) assertions? Missing `TreeSearch:::` prefixes? `set.seed()` before `sample()`? Edge-case coverage gaps (3-tip, single-char, all-NA)? Enduring regression for incremental-rescore? |
| 9 | **Wagner & addition trees** | `src/ts_wagner.h/.cpp`, `R/AdditionTree.R`, `R/PolEscapa.R` | **opus** | NA-incremental scoring staleness acceptable? Constraint mapping (LCA-based) correct? Retry loop fires? 3-taxon base case handles all orderings? R-layer index/`sequence` validation (OOB-write guard)? |
| 10 | **Alternative scoring kernels: Profile/IW/HSJ/XFORM** | `src/ts_fitch.cpp` (IW/profile paths), `src/ts_data.cpp` (precompute), `src/ts_hsj.cpp/.h`, `src/ts_sankoff.cpp/.h`, **`R/recode_hierarchy.R`**, **`R/CharacterHierarchy.R`** (added 2026-08-03 — see the rationale note; treat both as UNMEASURED), plus the criterion's **consumers** where a non-Fitch objective meets Fitch-only machinery: `src/ts_tbr.cpp` (candidate scan / accept / `try_root_edge_moves`), `src/ts_rcpp.cpp` (`unpack_hsj`, `unpack_xform`) | **opus** | `e/(k+e)` delta correct? Profile `info_amounts` lookup + capping matches? `concavity = 1.0` sentinel activates weighted path? `precompute_profile_delta` includes `precomputed_steps` offset? Clipped-subtree homoplasy in screening? HSJ/XFORM (`ds.hierarchy_blocks`/`ds.sankoff_*`) scoring correctness in its own right (not just collapse-flag blindness, cf. T-330 area 11) — does anything else outside collapse assume `ds.blocks[]` is exhaustive? |
| 11 | **Zero-length-branch collapse (MPT set)** | `src/ts_collapsed.cpp/.h`, `src/ts_splits.cpp` (`compute_collapsed_splits`), `src/ts_rcpp.cpp` (`ts_collapse_flags_batch`), `src/ts_tbr.cpp` (enum `add_collapsed` sites), `R/MaximizeParsimony.R` (collapse block) | **opus** | DEFAULT-ON since 2026-06-24, so every `MaximizeParsimony` call exercises it. Does `compute_collapsed_flags_aggressive` flag the *correct* min-length-0 branches under **IW / profile / NA**, not just EW (verified)? Is it really rooting-invariant, or does tip-rooting+`RenumberTips(labs)` alignment break on constraint trees / user start trees / `RenumberTips` permutations (cf. [[na-validation-alignment-gotcha]])? Can the dedup key `write.tree(SortTree(unroot(t)))` over-merge (two distinct collapsed topologies → same key) or under-merge across rootings? `result$scores == best_score` float-equality safe under IW/profile? Degenerate inputs: star tree, single MPT, 3–4 tips, all-resolved (must be exact no-op), fully-unresolved? Does collapse ever produce a tree that violates an active `constraint`? |
| 12 | **Red-team process meta-review** | `dev/red-team/focus-areas.md`, `dev/red-team/log.md`, the `red-team` issue list in `agent-issues/TreeSearch`, `dev/red-team/README.md` | **sonnet** | Are any areas too broad — spanning multiple distinct seams such that a finder concentrating on one file family misses another? Are any too narrow — a single-feature scope that would be better merged into a neighbour? Do any areas overlap (same source files audited under two different area headings)? Has any area gone persistently dry (≥ 3 consecutive rounds with zero confirmed findings) — should it be retired, merged, or downtiered? Are there new code seams (recently merged features, new source files) not covered by any existing area? Are tier assignments calibrated to actual yield recorded in `log.md` — any area that keeps surprising at its current tier and should escalate, or one that has been consistently empty and should drop? Propose concrete restructuring actions (split, merge, retire, add, re-tier) with rationale tied to `log.md` yield history. |
| 13 | **Constrained search correctness** | `src/ts_constraint.h/.cpp`, `src/ts_nni_perturb.cpp`, constraint integration points in `src/ts_driven.cpp` (fuse), `src/ts_parallel.cpp` (parallel-fuse), `src/ts_wagner.cpp`/`src/ts_sector.cpp` (posthoc retry), `src/ts_tbr.cpp` (`regraft_violates_constraint`) | **opus** | Does every `impose_constraint()` caller verify-before-capture, not just trust an improved score (T-213 gap, fixed d9a4f827: `nni_perturb_search` was the one caller that didn't re-check `constraint_node[]` after repair — fuse/parallel-fuse already did)? Any other heuristic-repair or posthoc-retry caller (Wagner build retry, sector) that skips discard-on-failure? Is `impose_one_pass`'s `best_node` reference stale after its own move-out loop's `topology_spr()` calls relocate a node — traced mechanism, produced one `std::bad_alloc` crash under experimental code, did NOT reproduce in 600 stress-test seeds against shipped code; needs a targeted adversarial tree construction, not more random seeds, to confirm either way. Is `map_constraint_nodes`/DFS-timestamp resync correct on every topology-mutation path, including reject paths (cross-check vs area 2's tabu-reject question)? Are nested/overlapping constraint splits handled consistently across TBR clip-gating, Wagner retry, and sector/fuse posthoc paths? |
| 14 | **Statistics & support metrics** | `src/MaddisonSlatkin.cpp`, `src/expected_mi.cpp`, `src/ts_mc_fitch.cpp`, `src/quartet_concordance.cpp`, `R/Concordance.R`, `R/ParsSim.R`, `R/pp_info_extra_step.r`, `R/WideSample.R`, `R/Consistency.R`, `R/TaxonInfluence.R`, `R/ScoreSpectrum.R`, `R/RandomTreeScore.R`, `R/WhenFirstHit.R`, `R/QuartetResolution.R`, `R/PresentContra.R`, `R/ClusterStrings.R` (last two added 2026-08-05 — owned by no other row, and both reviewed by the first-ever round) | **sonnet** | Is the recursive Maddison–Slatkin DP correct at its recursion boundaries, and does its cache key everything the recurrence depends on? Does the factorial-cache log-space arithmetic under/overflow at realistic tip counts, and are log-space sums accumulated stably? When does the exact DP hand off to the Monte Carlo fallback, and is the fallback's estimator unbiased — or silently substituted without the caller being able to tell? Are concordance-factor statistics well-defined on polytomies, on single-taxon splits, and on characters with missing data? Do the R wrappers validate tip-label correspondence, or index by position (cf. the [[na-validation-alignment-gotcha]] class)? Is any of this reachable from `MaximizeParsimony()`'s default output path the way #16/T-400 was — and does it return a silently wrong number rather than erroring? (carried from the 2026-08-05 round, where this question produced three of the four `sev:high` findings) |
| 15 | **Legacy pure-R search API** | `R/CustomSearch.R` (`TreeSearch()`), `R/Ratchet.R`, `R/NNI.R`, `R/SPR.R`, `R/TBR.R`, `R/SuccessiveApproximations.R`, `R/tree_rearrangement.R`, `R/morphy-deprecated.R`, `R/Bootstrap.R` | **sonnet** | Is `EdgeListScore()` — the default `TreeScorer` for `TreeSearch()`/`Ratchet()`/`Jackknife()`, and one of the four entry points #16 confirms vulnerable — reachable with the out-of-bounds inputs #16 describes? Do the pure-R rearrangement samplers (`NNI`/`SPR`/`TBR`) generate only valid topologies, and do they cover the neighbourhood they claim? Does `SuccessiveApproximations` reweight consistently with the C++ IW kernel, or has it drifted? Do `Bootstrap`/`Jackknife` resample characters with the weights the user supplied? Does anything here still route through removed MorphyLib paths (`morphy-deprecated.R`)? |

### Maturity / tier rationale (one line each)

- **1 Fitch correctness — opus.** Crown jewel; T-300 (systematic delta=−3) and T-306 were
  opus-class subtle bugs. Dry at **opus-4.8** (2026-07-24) ⇒ next visit is **opus (Opus 5)**,
  fresh-angle; **fable** is the escalation only once Opus 5 *also* runs dry (version bump before
  rung bump — see the model-version legend at the top of `log.md`).
- **2 Topology invariants — opus.** Deep state-restore subtleties; T-235 (SPR stale state),
  T-316 (P1 stale constraint metadata after tabu rejection).
- **3 Ratchet & perturbation — opus.** Mature, but the `build_reduced_dataset` /
  mask-sync class (T-275, T-303) keeps recurring.
- **4 Parallelism & RNG — opus.** Looked mature until T-309 (P1: R RNG API on worker
  threads via the parallel `Resample()` path, 2026-06-15). Concurrency = P1-capable.
- **5 Data pipeline — opus.** DAT-001 (`1u<<32` UBSAN at `n_states==32`), DAT-002 (XPIWE
  `obs==0` division). Edge inputs still bite.
- **6 R ↔ C++ interface — sonnet.** Mostly mechanical arg-count / sentinel audits; cheap to
  run. T-310 (frozen-API `pruneReinsertNni` type) shows it still occasionally yields —
  escalate if a sonnet pass goes dry.
- **7 Shiny wiring — sonnet.** **Immature seam:** a Sonnet pass found 5 bugs on 2026-06-16
  (T-309…T-313). Keep mining cheap until it runs dry. *(Scope row corrected 2026-07-27: it
  named a phantom `server/events.R` — no such file — while omitting `global.R` (430 lines),
  `ui.R`, `app_state.R` and `logging.R`, ~850 real lines that no round had ever owned. The
  2026-06-16 round read only `server.R` + `mod_search.R`, so `mod_consensus.R` (1449),
  `mod_treespace.R` (776), `mod_data.R` (660), `mod_clustering.R`, `mod_downloads.R` and
  `mod_references.R` remain unreviewed at any tier.)*
- **8 Test suite health — sonnet.** Reliably yields inline fixes (`set.seed`, vacuous
  asserts) and test-gap notes (T-304).
- **9 Wagner & addition — opus.** Kernel code; WGN-01 (P1 OOB write via `AdditionTree(sequence=)`).
- **10 Alternative scoring kernels — opus. THE HIGHEST-YIELDING AREA ON RECORD: 4 P1s in one
  round** (2026-07-28, Opus 5 — T-373…T-376, all in HSJ/XFORM; T-376 shows the HSJ score is not
  even a function of the data, and is wrong on a shipped dataset through the public API).
  **Read the cautionary tale before trusting any tier state here.** The area looked mature —
  "numerical delta algebra, subtle conservative bugs" — and its Profile/IW half genuinely is:
  every static residual was re-derived clean in that same round. But `ts_hsj.*`/`ts_sankoff.*`
  were added to the scope row on 2026-07-03 and **no finder had read them at any tier**, so
  both of the area's dry rounds (2026-05-19, 2026-05-26 — themselves pre-tier and
  version-unrecorded, hence *zero* version-scoped dry verdicts) were dry about a **different
  half of the area than the row described**. Generalisable lesson: **a scope row that grows
  does not inherit the dry verdicts earned before it grew** — when widening a row, reset its
  measured maturity for the new files. Stays opus; next visit resumes on the HSJ/XFORM half
  (`ts_sankoff.*` got materially less attention than `ts_hsj.*`), and carries T-374's open
  design question — whether a rooted objective is *intended* for XFORM, and whether Hopkins &
  St John's dissimilarity is defined on a single MPR or on marginal MPR sets, which is a
  **literature** question, not a finder question.
  **UPDATED 2026-08-03 (second consecutive high-yield round: 8 filed + T-377 reopened).** The
  HSJ literature question is now **settled from the paper** (rooting-invariance is required) and
  the HSJ half of T-374 is **fixed and merged** (`34eea581`); the XFORM design question is with
  the user. Three things changed about how this area should be attacked:
  (a) **The live seam has moved OUT of the kernels and INTO the consumers.** `ts_sankoff.*` was
  re-read in full and came back with a **reasoned clean on its most likely defect** (cost-matrix
  orientation on the asymmetric matrix — no transposition, argument in `log.md`), while the
  round's headline (T-392, P1) is in `src/ts_tbr.cpp`: `try_root_edge_moves` adopts a strictly
  worse HSJ/XFORM score with no comparison and no restore, on the **shipped default** path
  (`.AutoRung` → `sprint` → `tabuSize = 0` for every dataset ≤30 taxa). Attack the machinery
  that *consumes* these criteria, not just the kernels that compute them.
  (b) **Two R files produced three findings while sitting in NO area's scope row** — T-393/T-394
  (`R/recode_hierarchy.R` derives a secondary's state space from token strings) and T-395
  (`R/CharacterHierarchy.R` rejects every nested hierarchy, so the documented example is
  unusable while `R CMD check` stays green). Both are now in the row above, and per this area's
  own 2026-07-28 lesson they carry **no inherited maturity**: a grown scope row does not inherit
  dry verdicts, so treat them as never-reviewed. The 2026-07-03 area-12 round had already
  flagged `R/CharacterHierarchy.R` as unowned — that note sat unactioned for a month and cost
  three findings' delay.
  (c) **The next visit should NOT be a finder.** T-392 + the reopened T-377 together say
  HSJ/XFORM search hill-climbs on the Fitch residue with the hierarchy term as a mere accept
  filter, and the reach cost is now measured at 25 tips — but the **fix shape is undecided** and
  the candidate gate costs 2.7–2.9× wall. What is owed is a **wall-matched A/B of the fix shapes
  on a ≥25-tip hierarchy-heavy matrix, Hamilton-class**, mirroring the area-13 precedent where
  the right next step was a bounded harness rather than another review.
- **11 Zero-length collapse — opus, NEVER REVIEWED.** New + default-on; implementation hit
  three subtle traps in one sitting (conservative flags rooting-sensitive; aggressive flags
  need tip-rooting; tip-data alignment via `RenumberTips`). Cross-mode correctness
  (IW/profile/NA) and dedup canonicalization are the under-verified seams. Start opus.
- **12 Red-team meta-review — sonnet.** *(Corrected 2026-08-04 — this line read "NEVER
  REVIEWED" despite a 2026-07-03 round having already run; see `log.md` for both rounds.)*
  Pure document review; no code to trace. Cheap per round; findings are restructuring
  proposals (split/merge/retire/add/re-tier) that improve every future round. Escalate to
  opus only if evaluating a proposed split requires reading source files to assess scope
  boundaries. Still yielding at sonnet on both rounds to date — the file-coverage diff alone
  (2026-08-04) surfaced ~37 of 44 `R/*.R` files and ~12 `src/*` files owned by no row.
- **13 Constrained search — opus.** *(Rationale corrected 2026-08-04 — the text below had
  read "NEVER REVIEWED" since 2026-07-02 through two rounds that reviewed it.)* Split out
  2026-07-02 after fixing T-213 (`nni_perturb_search` missing verify-before-capture on
  `impose_constraint()`'s heuristic repair — cf. [[impose-constraint-verify-gap]]). The
  2026-07-02 round traced but did NOT confirm a second failure mode (stale `best_node` across
  `topology_spr()` relocation inside `impose_one_pass`); the 2026-07-03 round's audit of every
  `impose_constraint()` caller sharpened that into a precise, still-open question (can the
  `postorder.size()==n_internal` revert-guard be slipped by a net-zero corruption?) and closed
  with **"NEXT VISIT: NOT another finder — a BOUNDED EXHAUSTIVE HARNESS"** on the
  `topology_spr`/`build_postorder`-guard equivalence, mirrored in `escalation-backlog.md`'s
  "Not in this backlog (deliberately)" section. **That verdict is still standing and rotation
  reaches area 13 next** (after this area-12 round) — whoever dispatches that round should read
  it before defaulting to a finder. Two findings were ALSO filed into this area's scope from a
  2026-07-28→2026-08-04 area-11 round without an area-13 finder ever running: **#18/T-402
  (sev:high)** — a `constraint` silently ignored when `tree=` supplies a violating start — and
  **#19/T-403 (sev:med)** — collapse's "enforced splits protected" promise has no access to
  `consZero` and can return constraint-violating trees under default `collapse=TRUE`
  (`escalation-backlog.md` item 7). Deliberately **not** `needs-escalation`-labelled: that flag
  means only "dispatch at `opus`+", which this row already is, and a label hit makes step 3 skip
  reading the backlog row that holds the actual ask (item 7 explains this at length). Whoever
  takes area 13 next must decide explicitly: harness first, or #18/#19 first — both are live,
  and the harness plan predates the two findings.
- **14 Statistics & support metrics — MEASURED 2026-08-05, still yielding heavily.** Added
  2026-08-05 from #42's scope-coverage diff: 5,553 lines across 14 files that were owned by no
  area and therefore never reviewed at any tier. **The gap has already cost a finding** — the
  arm64 `probe_slot()` hang in `src/MaddisonSlatkin.cpp` (fixed, PR #272,
  cf. [[maddisonslatkin-arm64-profile-hang]]) was found incidentally, not by rotation. The code
  is numerically dense — recursive DP, factorial caches, log-space arithmetic, Monte Carlo
  fallbacks — the profile the tier doctrine normally reserves for `opus`, and #42 recommended
  `opus` on that basis. **Deliberately starting at `sonnet` anyway** (maintainer decision,
  2026-08-05): density is a prediction about where bugs *hide*, not evidence that cheap sweeps
  are exhausted, and this area has no measured yield at all. **Overtaken by events:** the
  first-ever review had already run at `opus` on 2026-08-05, before this row merged, and returned
  **36 findings, 4 sev:high — the highest yield on record for this rotation** (see `log.md`).
  `start_tier` is left at `sonnet` as decided, but it is now inert: the seam is measured and
  yielding, so the routing rules keep the next visit at **opus** with a fresh agent.
  **Next visit starts here** (the round's own leads, and the reason it stays opus): the
  **array-dimension-drop pattern** — four independent instances in one round (`ConcordanceTable`,
  `ClusteringConcordance`, `Consistency`, `ClusterStrings`, all missing `drop = FALSE`), so treat
  it as a class and sweep for it rather than re-finding instances; and the **not-yet-examined
  `R/PresentContra.R` forest/reference-tip-mismatch angle** — read but never exercised against a
  forest whose trees have tips absent from the reference (it calls `KeepTip` first, which *should*
  be safe, but that is unproven). Its own test convention
  (`test-MaddisonSlatkin.R`, `test-Concordance.R`, `test-ParsSim.R`, `test-Consistency.R`,
  `test-ScoreSpectrum.R`, `test-QuartetResolution.R`, `test-TaxonInfluence.R`,
  `test-WideSample.R`, `test-pp-*.R`) is a useful first read.
- **15 Legacy pure-R search API — sonnet, UNMEASURED / no inherited maturity.** Added
  2026-08-05 from #42's scope-coverage diff: 2,183 lines across 9 files backing the
  still-shipped pre-C++-engine search functions, owned by no area. **Higher urgency than its
  size suggests:** #16 (`sev:high`) names `EdgeListScore()` as *"the default `TreeScorer` for
  `TreeSearch()`/`Ratchet()`/`Jackknife()`"* and one of four confirmed-vulnerable entry points,
  so this family is a second, wholly unreviewed exposure surface for an already-confirmed bug —
  take that question first. #42 offered "review once as frozen legacy, then deprioritise";
  **the maintainer chose a full rotation area instead (2026-08-05): keep revisiting until the
  seam stops yielding.** Legacy is not the same as clean, and this code is still shipped and
  still the documented entry point for users who have not moved to the C++ engine. Treat "it
  isn't growing" as a reason the seam should *exhaust* quickly, not as a reason to stop early.
