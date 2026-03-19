# TreeSearch Completed Tasks Archive

Tasks moved here from `to-do.md` on completion. Newest first.

---

## 2026-03-19

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-142 | Shiny: Add TreeSearch logo to app header | A | Inline SVG (magnifier + 3-tip tree) added to `inst/Parsimony/ui.R` line 14; flex div wraps icon + h1. No external asset file needed. |
| S-COORD | Coordination review round 8 | A | T-144 fixed → CRAN unblocked. ~9835 pass/0 fail. T-141/T-140/T-097 confirmed complete. Standing tasks P1. |
| T-144 | Fix 15 CRAN Tier 1 test failures (PrepareDataProfile regression) | A | Added binary-reduction warning to PrepareDataProfile; fixed empty-phyDat return (avoid dataset[0] crash in new TreeTools); updated test-data_manipulation.R (17→17 pass) and test-Concordance.R (67→67 pass) to match new behavior. 0 regressions. |
| T-129 | Shiny: Evaluate progressive search result display | A | Briefing in `.positai/briefing-progressive-results.md`. Rec: progress file polling via existing C++ progress_callback infrastructure (~50 lines). Do not stream partial trees (misleading). Filed T-141. |
| S-PROF | Standing: Performance profiling (round 3) | A | EW no regression post T-115–T-124. HSJ ~0.6× EW (Fitch screening effective). XFORM ~1.7× EW (Sankoff overhead in Ratchet/Drift). Hierarchical resampling faster per-rep but no inter-replicate parallelism (1.1× vs 2.5× for Brazeau). No new optimization tasks. |
| S-RED | Red-team focus 7: Shiny module wiring | B | Reviewed all 7 modules + server.R. Forward-ref callbacks correct, cross-module updates correct (1 fragile parent_session usage documented), no orphaned observers, isolate() patterns correct, new progress observer correctly gated. No bugs found. |
| T-143 | Shiny: Per-replicate search progress display | B | Progress file polling: R callback writes `rep max_rep score hits target` to TREESEARCH_PROGRESS_FILE on each replicate. Shiny polls every 500ms via invalidateLater. Fixed C++ callback to fire regardless of verbosity. Updated test-ts-progress.R expectations. 1681 ts-* + 42 mod-search pass. |
| S-RED | Red-team focus 5: Data pipeline & simplification | B | Fixed `build_reduced_dataset` missing `inapp_state` copy (latent HSJ bug). Added `n_states > MAX_STATES` guard in `build_dataset`. Verified simplification, constraint indexing, EW offset interaction. 10 new NA-sector tests. 1679 ts-* pass. |
| T-130 | Shiny: Search cancellation button | B | Already fully implemented: hidden Stop button shown during search, file-based cancel signal polled by C++ (serial + parallel), cleanup on completion. 3 dedicated testServer tests, 38/38 mod-search pass. |
| T-147 | Shiny: "Tips to show" starts at 0, dimension error | F | Same root cause as T-146: `UpdateKeepNTipsRange()` reactive never consumed on init. Added `observe(UpdateKeepNTipsRange())`. 92/92 module tests pass. |
| T-146 | Shiny: "Root on" selectize not populated on launch | F | `UpdateOutgroupInput()` was a reactive never consumed on init. Added `observe(UpdateOutgroupInput())` in mod_consensus.R. 92/92 module tests pass. |
| T-145 | Shiny: Defer profile prep to search start | F | Removed auto-trigger on mode change; profile prep only runs when user clicks Search in profile mode. Error notification → silent LogMsg. On success: "Profile scores ready — click Search to start." 92/92 module tests pass. |
| T-097 | Verify: Sun 2018 continued-search indicator | F | Verified via SearchLog integration test (4/4 pass): EW→IW consecutive search works, `searchCount` increments. Sun2018 loads/scores correctly (54t, 225c, EW=802). Core fix in T-090, additional fix in T-138. |
| T-140 | Shiny: Profile prep error on dataset selection | F | Resolved by T-138 (lazy async prep). Added re-invocation guard: if user changes dataset while profile prep is running, cancel in-flight task via signal file and retry after 200ms. 88/88 module tests pass. |
| T-138 | Shiny: "Continue search" doesn't start new search | F | Root cause: `bindCache` on `scores()` prevented re-evaluation after dataset switch. Fix: replaced with plain `reactive()`, added `r$searchInProgress` flag, async profile prep with cancel, cancel button UI, starting-tree error handling. 5 new testServer tests. |
| S-RED | Red-team review (focus 6) | E | Verified T-137/139/141 fixes. Scores correct (serial+parallel). hits≤reps confirmed. 1676 ts-* pass, 88 module pass. 42 pre-existing failures from human's profile parsimony refactor (tracked by T-144). No new bugs. |
| T-141 | Shiny: Sun 2018 blank plot area despite trees | E | When loading a dataset whose file has no trees, stale trees from the previous dataset persisted with incompatible tips → blank plot. Fix: clear old trees when they don't match the new dataset. 88 module tests pass. |
| T-139 | Shiny: Hit/rep count forgets previous searches | E | Bug: `parallel_driven_search()` captured `hits_to_best` AFTER MPT enumeration, inflating count. Fixed to capture before enumeration (matching serial path). Added clamp in `SearchConfidenceText()` as safety net. 1669 ts-* + 275 parallel pass. |
| T-137 | Shiny: Stop button leaves UI stuck on "Stopping…" | E | C++: threaded cancel callback into `ratchet_search()`/`drift_search()` for faster response. Shiny: cancel observer removes notification immediately; result observer uses `searchInProgress` flag for robust cleanup. 1669 ts-* + 88 module tests pass. |
| T-136 | Shiny: Use `WideSample()` for tree thinning | F | Replaced `seq.int()` with `WideSample()` maximin selection in `mod_data.R`. Updated Distribution snapshots. 119/120 Shiny pass. |
| T-124 | Inapplicable: Hierarchical resampling in `Resample()` | F | `hierarchy`, `inapplicable`, `hsj_alpha` params added. Resamples at unit level (free chars + hierarchy blocks). `.HierarchicalResampleWeights()` + `.ResampleHierarchy()`. 72 tests, 1669 ts-* pass. |
| T-123 | Inapplicable: `TreeLength()` HSJ + xform extension | D | Added `hierarchy`, `inapplicable`, `hsj_alpha` to all TreeLength methods. HSJ via `ts_hsj_score()`, xform via Fitch + `ts_sankoff_test()`. 24 new tests, 73/73 tree_length pass. |
| T-125 | Inapplicable: Documentation & vignette | F | `vignettes/inapplicable.Rmd`: Brazeau/HSJ/xform approaches, hierarchy specification, worked example, comparison table. Updated MaximizeParsimony docs. |
| T-122 | Sankoff: Tests against Goloboff et al. examples | F | 80 assertions in test-ts-xform.R: gain/loss asymmetry, secondary variation penalty, HSJ vs xform cross-validation, Fitch+Sankoff mixed scoring, 8-tip search, deterministic seeds, gain cost scaling. |
| T-121 | Sankoff: Integration with search pipeline | F | `ScoringMode::XFORM`, `score_tree()` dispatch (Fitch+Sankoff), `ts_driven_search` bridge, `MaximizeParsimony()` xform wiring. End-to-end xform search works. |
| T-120 | Sankoff: C++ optimization engine | F | Fixed incorrect asymmetric cost test expectation (score 2, not 3). 24/24 Sankoff pass. |
| T-119 | Sankoff: R-level recoding function | F | `recode_hierarchy()`: primary+secondaries → Sankoff char, asymmetric cost, Hamming distance. 49 tests. |
| T-134 | HSJ: Fix secondary dissimilarity d always 0 | F | Added Fitch uppass to `fitch_label_char()`. 37/37 HSJ pass, 1509 ts-* pass. |
| T-133 | ParsSim: per-taxon/per-character missing rates | A | List/matrix inputs; Bernoulli per-cell sampling. 128/128 ParsSim pass. |
| T-132 | ParsSim: missing data (`?`) support — flat rate | A | `missing` param (0–1) injects `?` post-hoc. 97/97 ParsSim pass. |
| T-131 | Fix: 39 stale IW reference values in test-ts-iw.R | A | Recomputed after T-113 NA bit-stripping. 86/86 IW pass. |
| T-128 | Shiny: Rename "Mode" label to "Step weighting" | ? | Trivial label change. |
| T-120 | Sankoff: C++ optimization engine | F | `ts_sankoff.h/.cpp`, 24 unit tests. |
| T-118 | HSJ: End-to-end tests against paper examples | C | 123 assertions in test-ts-hsj.R. Found d=0 bug → T-134. |
| T-117 | HSJ: Wire into search pipeline + remove placeholder | D | `ScoringMode::HSJ`, `score_tree()` dispatch, driven search bridge. |
| T-116 | HSJ: Rcpp bridge + R marshalling | D | `ts_hsj_score()`, `build_tip_labels()`, `hierarchy_to_blocks()`. |
| T-115 | HSJ: Token mapping & DataSet partitioning | D | `inapp_state` in DataSet, `absent_state` in HierarchyBlock, `partition_weights()`. |
| T-114 | ParsSim: rootState vector handling + docs | D | Validate length/range, per-character scalar indexing. 80 ParsSim pass. |
| T-113 | Fix: T-097 NA ambiguity code fix incomplete | B | Bit-stripping in `build_dataset()`. 12/12 NA-ambig pass. |
| T-112 | Port missing MaddisonSlatkin tests from branch | D | All branch tests already present. 37/37 pass. |
| T-111 | ParsSim: extended test suite + edge cases | ? | 9 extended tests. 66/66 ParsSim pass. |
| T-110 | ParsSim: profile character selection | D | `concavity = "profile"` mode. 67 ParsSim pass. |
| T-109 | ParsSim: core implementation (EW + IW) | F | `R/ParsSim.R` with 7 internal helpers. 35 tests. |
| T-108 | Fix: pattern_freq multiplicative blowup in IW ratchet | B (S-RED) | `+= 1` instead of `*= 2`. 1404 ts-* pass. |
| T-105 | Docs + CRAN prep for multi-state profile | D | Updated vignette, roxygen, NEWS.md. |
| T-104 | Integration tests for multi-state profile parsimony | D | 6 test blocks in test-ts-profile.R. |
| T-103 | Generalize PrepareDataProfile() for multi-state | B | Removed binary decomposition, dynamic contrast matrix. |
| T-102 | Generalize StepInformation() for multi-state | B | 3–5 state dispatches to MaddisonSlatkin. 64 new tests. |
| T-101 | Port MaddisonSlatkin.cpp to main branch | D | 1401-line C++, 24 test assertions. |
| T-100 | Fix: `{-,X}` ambiguity tokens treated as applicable | B (S-RED) | Stale postorder in drift_phase. Also: bit-stripping in ts_data.cpp. |
| T-099 | Confidence text: hover tooltip | B | Title tooltip on summary + modal. 28/28 mod-search pass. |
| T-098 | Fix overconfident search confidence probability | B | `exp(-K)` conservative bound, Laplace smoothing. 76/76 mod pass. |
| T-097 | Fix: {inapplicable,state} ambiguity tokens scored as applicable | F | `{-,S}` tips collapsed to pure NA. 7 new tests. 1404 ts-* pass. |
| T-096 | Fix: upweight_mask missing in EW bounded/cached indirect scoring | C (S-RED) | 6 sites fixed. 1397 ts-* pass. |
| T-095 | Assess MorphyLib decoupling for custom search | ? | Thin wrappers around `ts_fitch_score()`. Design in morphy-migration.md. |
| T-094 | Un-deprecate custom search functions | B | Removed `.Deprecated()` from Ratchet, Jackknife, etc. |
| T-093 | Remove testthat-problems.rds artifact | B | Deleted + .gitignore. |
| T-092 | Fix: hits_to_best inflated by MPT enumeration | B | Capture hits_to_best before plateau walk. 1397 ts-* pass. |
| T-091 | Shiny: Terminology cleanup — "runs" not "replicates" | ? | Already resolved by T-088/T-089/T-090. |
| T-090 | Shiny: Search-in-progress indicator fix + results display refresh | B, C | Moved notification before tree selection; always call `DisplayTreeScores()`. 94/94 shinytest2 pass. |
| T-089 | Shiny: Confidence text rewrite + adaptive slider note | B | "Probability that a better tree exists: ~X%". 92/92 shinytest2 pass. |
| T-088 | Shiny: Search config modal cleanup | B | Renamed labels, removed helpText, restructured layout. 92/92 shinytest2 pass. |

## 2026-03-18

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-087 | Doc: `.AutoStrategy()` docs vs code mismatch | A | Fixed roxygen + AGENTS.md. |
| T-086 | Replace `suppressWarnings` with `expect_warning` in test-ts-profile.R | D | |
| T-085 | R CMD check refresh: 0E/0W/3N (benign) | C | Added test artifacts to `.Rbuildignore`. |
| T-084 | NEWS.md: Shiny/EasyTrees improvements section | F | |
| T-082 | Shiny: nThreads slider | — | Already implemented (triage). |
| T-083 | Shiny: Search config modal restructured | — | Already implemented. |
| T-081 | Shiny: Timeout slider verified done | F | Dataset-adaptive scaling. |
| T-080 | Shiny: renamed stop label | C | |
| T-079 | Shiny: Fix `maxProjDim` crash | A | |
| T-078 | Shiny: `PlotCharacter()` crash on multifurcating tree | C | |
| T-077 | Shiny: post-search confidence display | D | |
| T-076 | Wrap PrepareDataProfile warnings in test-ts-profile.R | C | |
| T-075 | Bench: ns/candidate cost linear in n_blocks | A | |
| T-072 | Shiny: result observer accumulates trees at matching score | C | |
| T-071 | Shiny: first-search trees unchanged (fixed via T-072) | C | |
| T-070 | Fix ts_wagner.cpp -Wcomment | C | |
| T-065 | Shiny mod: testServer() coverage pass | C | 92 assertions total. |
| T-064 | Shiny mod: Dissolve events.R + final integration | C | 75/75 shinytest2 pass. |
| T-063 | Shiny mod: mod_consensus | C | 1327-line module. 75/75 shinytest2 pass. |
| T-062 | Shiny mod: mod_data | B | Absorbs data.R + trees.R. 70/70 shinytest2 pass. |
| T-061 | Shiny mod: mod_search | A | 51/51 shinytest2 pass. |
| T-060 | Shiny mod: mod_clustering | B | 57/57 shinytest2 pass. |
| T-059 | Shiny mod: mod_treespace | A | |
| T-058 | Shiny mod: mod_downloads | E | 11-assertion testServer. |
| T-057 | Shiny mod: mod_references | ? | 4-assertion testServer. |
| T-056 | Shiny mod: AppState typed store | B | 34+ fields, 51/51 shinytest2 pass. |
| T-055 | Shiny mod: Three-file skeleton + source extraction | B | Split `app.R` → 3 files + 11 server/ files. |
| T-054 | Shiny mod: Capture shinytest2 baselines | B | FAIL 0 / WARN 0 / PASS 29. |

## 2026-03-17

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-052 | Shiny: seed=TRUE, cache keys, library() fixes | A | |
| T-051 | Shiny: Removed redundant TreeDist install.packages | A | |
| T-050 | Shiny: `[1:3]` → `head(, 3)` in Excel preview | A | |
| T-049 | Shiny: Dead `r$newTrees` — added missing write | A | |
| T-048 | Shiny: onStop cleanup — tempdir() prefix + Excel cleanup | A | |
| T-047 | Shiny: `which.min(NULL)` crash guard | A | |
| T-045 | Shinylive plan for EasyTrees | D | Plan written, sub-tasks TBD. |
| T-044 | Documentation: custom.Rmd vignette framing | A | |
| T-043 | Documentation: SuccessiveApproximations() parameter names | D | |
| T-042 | Documentation: DESCRIPTION, README, cross-references | D | |
| T-041 | Fix Shiny app (EasyTrees) for new MaximizeParsimony params | D | |
| T-040 | All-ambiguous phyDat guard in TreeLength/MaximizeParsimony | D | |
| T-039 | Fully-resolving constraint crash: column-major indexing bug | A+D | |
| T-038 | Final CRAN prep: version bump 1.9.0, R CMD check | D | |
| T-037 | Post-T-025 cleanup: sprFirst default, regression benchmark | C | |
| T-036 | R-level test coverage for new MaximizeParsimony features | E | |
| T-035 | Deprecate legacy search functions | C | |
| T-034 | Migrate RandomTreeScore to C++ engine | A | |
| T-033 | Write NEWS.md for next release | D | |
| T-032 | Disable CSS by default (cssRounds=0) | C | |
| T-031 | Named strategy presets in R interface | C | |
| T-030 | Tier 1 MorphyLib migration: TreeLength + CharacterLength → C++ | C | |
| T-029 | Default parameter tuning: driftCycles 6→2, ratchetCycles 10→5 | C | |
| T-028 | Pass `min_steps` for IW scoring in MaximizeParsimony() | C | |
| T-027 | progressCallback SIGSEGV (secondary symptom of T-025) | E | |
| T-026 | Performance regression benchmark script | C | |
| T-025 | PreallocUndo buffer overflow (P0 crash) | E | |
| T-024 | Parallel resample (jackknife/bootstrap) | E | |
| T-023 | Expose `maxSeconds` in MaximizeParsimony R interface | B | |
| T-022 | R CMD check preparation | D | |
| T-021 | Fix pre-existing test failures in legacy test files | A | |
| T-020 | Fix RcppExports/init.c mismatch for Wagner bridges | A+B | |
| T-019 | Migrate AdditionTree() to C++ Wagner engine | E | |
| T-018 | Pass user-supplied starting tree to C++ engine | B | |
| T-017 | Add test coverage for ambiguous-token simplification | B | |
| T-016 | Add `precomputed_steps` offset in `precompute_profile_delta` | B | |
| T-015 | Copy scoring_mode + simplification fields in build_reduced_dataset | B | |
| T-014 | Fix `compute_fixed_steps` for all-ambiguous characters | B | |
| T-013 | Fix `is_uninformative` for ambiguous tokens | B | |
| T-012 | SPR→TBR escalation in driven search | E | |
| T-011 | Documentation refresh (Morphy, MaximizeParsimony, vignette, README) | B | |
| T-010 | MorphyLib deprecation plan | B | |
| T-009 | Audit rearrange.cpp dead code | B | |
| T-008 | Fix test-ts-simd.R exit code 127 crash (fuseInterval=0) | B | |
| T-007 | Wagner NA-incremental scoring | E | |
| T-006 | Audit and clean R-level TODOs | C | |
| T-005 | Phase 6E: Adaptive strategy selection | A | |
| T-004 | Phase 6D: Benchmarking framework | D | |
| T-003 | Phase 6C: Define strategy space | B | |
| T-002 | Phase 6B: Curate benchmark dataset suite | B | |
| T-001 | Phase 6A: Per-phase timing instrumentation | A | |
