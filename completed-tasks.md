# TreeSearch Completed Tasks Archive

Tasks moved here from `to-do.md` on completion. Newest first.

---

## 2026-03-25

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-213 | Implement impose_constraint() for post-hoc topology repair | D (cleanup) | Already on cpp-search (a666918ed, PR #223). 88 tests pass. Formal closeout during S-COORD. |
| T-220 | [Shiny] Crash: searchExtendedIw not found when clicking Continue | D | Variable used in LogCode() before assignment. Moved snapshot above LogCode(). Direct fix on cpp-search. |
| T-219 | [Shiny] Dataset dropdown hover state visible | D | Selectize default hover (#f5f5f5) near-invisible on white. Added explicit hover CSS (#dde6ed). |
| T-211 | Stale `final_` in temper candidate scoring | C | Analyzed: conservative-only. Stale `final_` after clip/evaluate/restore biases Boltzmann screening but not verified acceptance (`temper_full_rescore` gates all accepted moves). Fix would require per-candidate full rescore or save/restore of all `final_` arrays — cost exceeds negligible SA benefit. Closed as not worth fixing. |

## 2026-03-24

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-190 | Adaptive starting-tree strategy mixing (bandit) | A | Thompson sampling over 4 fresh-start arms (3 Wagner variants + random tree). 33 unit tests, vignette section, strategy diagnostics attribute. Benchmarked on 75–88 tip datasets (neutral to slight benefit). Squash-merged via PR #214. |
| T-177 | Bug fix: mid-TBR/SPR timeout | G | `check_timeout` callback threaded through `tbr_search`, `spr_search`, `nni_search`. 282 targeted tests pass. |
| T-205 | Fix flaky test-pp-random-tree.R on Windows | G | MWC RNG in build_postorder.h uses static global state not seeded by set.seed(). Widened binomial bounds (stringency 0.005→1e-6) and increased nTrees (6000→12000, 12000→24000) across all tests. False-positive rate: ~0.0002% per run (was ~1%). GHA pass: run 23501977394. |
| T-203 | Simulated annealing for large trees | G | Linear cooling schedule (T_start→T_end over N phases) using stochastic TBR + Boltzmann acceptance. `ts_temper.h/.cpp` (Layer 1: stochastic_tbr_phase ported from T-198; Layer 3: anneal_search). Wired into driven pipeline between drift and final TBR polish. `SearchControl(annealPhases, annealTStart, annealTEnd, annealMovesPerPhase)`. `large` preset: drift disabled, 5 annealing phases T=20→0. 19 new tests; all pass. Merged to `cpp-search` (conflict with `enumTimeFraction` resolved). |
| T-197 | Fix `concavity = 0` NaN in `precompute_iw_delta` | D | C++ guard for e==0 avoids 0/0 NaN. R entry points already validate concavity>0. Added 8 new validation tests (MaximizeParsimony, SuccessiveApproximations, TreeLength, AdditionTree). 169 tests pass. |
| T-195 | GHA benchmark workflow | D | `agent-benchmark.yml` + `bench_regression.R` CLI args (`--datasets`, `--budget`, `--output`, `--threads`, `--lib`). 14 datasets with max_score/ref_time_s. CSV artifact upload. Commit `7a80e67a`. |
| T-202 | Fix MPT enumeration skipped on timeout | B | Two-phase timeout: main loop exits at `budget*(1-enumTimeFraction)`, reserving remainder for plateau walk. New `SearchControl(enumTimeFraction=0.1)`. PR #217 merged. |
| T-179 | Large-tree strategy preset (>=120 tips) | G | Tuned via systematic benchmarking on mbank_X30754 (180t, 418p). Key: NNI-perturb too expensive at 5.5s/cycle; ratchet 12, drift 4, no NNI-perturb, outerCycles=1, single biased Wagner, tbrMaxHits=1. 60s: median 1255 vs thorough 1259; 120s: tied at 1250 but 2 reps vs 0-1; 30s: 1276 vs 1283. Commit `fab1e52c`. |
| S-RED | Red-team focus 10: Profile & IW scoring | B | Filed T-196 (P2): `extract_divided_steps` NA+IW bug — four static copies read `local_cost` for NA blocks instead of three-pass corrected steps, mispricing IW candidate screening. Filed T-197 (P3): `concavity = 0` → NaN, no validation. Verified: profile `concavity = 1.0` sentinel correct; `precompute_profile_delta` precomputed_steps offset correct; all indirect IW variants structurally correct. |
| S-RED | Red-team focus 11: T-190/T-202/XPIWE merge review | F | Filed T-208 (P2): `random_topology_tree()` ignores constraints — bandit RANDOM_TREE arm can produce constraint-violating starting tree when `adaptiveStart=true` (thorough preset). TBR blocks all constraint-relevant moves, tree returned to user unvalidated. Verified: XPIWE scoring correctness (eff_k, phi, sector copy, resampling). T-202 two-phase timeout correct (both serial and parallel). T-190 StrategyTracker correct (Thompson sampling, decay, RNG safety). |
| T-194 | Stratified sample selection + dedup profiling | B | Profiled 35 multi-file projects; flagged 24 near-duplicates (≥95% char identity) as `dedup_drop`. Post-dedup: 659 usable (535 train, 124 val). Selected fixed 25-matrix training sample via max-min distance (`MBANK_FIXED_SAMPLE`). Dedup integrated into `build_mbank_catalogue.R`. Documented in `strategies.md` and `AGENTS.md`. |
| T-191 | MorphoBank matrix catalogue | B | Scanned 801 .nex files from neotrans/inst/matrices/. 797 parsed OK (4 failures). 683 usable after ntax≥20 filter: 554 training, 129 validation (project%5==0). Includes multi-matrix projects and 7 syab files. Output: `dev/benchmarks/mbank_catalogue.csv`. Script: `dev/benchmarks/build_mbank_catalogue.R`. |
| T-192 | External dataset loading functions | B | Added to bench_datasets.R: `load_mbank_catalogue()`, `load_mbank_datasets()`, `load_mbank_sample()` (stratified by tier), `load_mbank_split()`. Path auto-resolved to `neotrans/inst/matrices/`. |
| T-193 | MorphoBank benchmark runner integration | B | Added to bench_framework.R: `benchmark_mbank_sample()` (routine ~25 matrix training sample), `benchmark_mbank_sweep()` (full split), `benchmark_mbank_validation()` (one-way-door validation with warning). All results tagged with `source` column. End-to-end tested: 4 matrices x default x 5s. |

## 2026-03-23

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-177 | Bug fix: mid-TBR/SPR timeout callback | G | Verified `check_timeout` callback threaded through `tbr_search`, `spr_search`, and `nni_search` — all poll periodically (every `n_tip` clips) and bail mid-pass. Driven pipeline, ratchet, drift, NNI-perturb all pass callback. 282 targeted tests pass (0 fail). Closing as complete. |
| — | Ratchet perturbation tuning: 4%→25%, moves 20→5, cycles 5→10 | Human+AI | Systematic sweep across 14 datasets. 9 improved, 4 unchanged, 1 marginal at 10s (resolves at 20s). Commit `f1ae7edb`. |
| — | Drift→ratchet reallocation: driftCycles 4→2, ratchetCycles 10→12 | Human+AI | Drift ~0 per-replicate improvement; ratchet is strictly better use of budget. Commit `7ae01181`. |
| — | Large-tree profiling: 180-taxon dataset analysis | Human+AI | Discovered NNI essential at >100 tips, timeout bug in TBR, strategy presets not calibrated for large trees. Filed T-177 through T-183. |
| T-185 | Inspect IQ-TREE for parsimony search acceleration ideas | G | Reviewed `iqtree.cpp`/`iqtree.h` source. Top idea: stochastic NNI-perturbation (complement to ratchet). Also: diverse starting trees, adaptive perturbation scaling, perturbation-count stopping. Batch NNI not worthwhile (see `.positai/expertise/batch-nni.md`). |
| T-186 | Stochastic NNI-perturbation as escape mechanism | G | New `ts_nni_perturb.h/cpp`: random compatible NNI swaps on ~50% of branches + TBR re-optimization. Integrated between ratchet and drift in driven pipeline. `SearchControl(nniPerturbCycles, nniPerturbFraction)`. `thorough` preset: 5 cycles. 28 new test assertions; 1792 ts-* pass. |
| T-178 | NNI warmup in driven pipeline | G | NNI always-on (`nni_first = true` default). Each Wagner start NNI-optimized before selection (best of NNI-local optima). SPR auto-skipped when NNI active (NNI→TBR empirically optimal). Constraint guard: NNI warmup disabled when constraints active (nni_search lacks constraint support). All presets updated: `nniFirst = TRUE, sprFirst = FALSE`. 1846 ts-* pass. |
| T-156 | XPIWE C++ core | G | Already implemented in feature/xpiwe commit c7a41712. Verified: ScoringMode::XPIWE enum, eff_k[]/phi[] per-pattern vectors in DataSet, build_dataset() computes adjusted concavity, compute_iw()/precompute_iw_delta() use per-pattern eff_k[p]. Branch builds clean, 1677 ts-* pass. |
| T-157 | XPIWE Rcpp bridge | G | Already implemented in c7a41712. xpiwe bool + xpiwe_r + xpiwe_max_f + obs_count params in make_dataset, ts_fitch_score, ts_driven_search, ts_resample_search, ts_successive_approx, ts_parallel_resample. TreeSearch-init.c updated. |
| T-158 | XPIWE R API | G | Already implemented in c7a41712. extended_iw param in MaximizeParsimony(), TreeLength() (all S3 methods), Resample(), SuccessiveApproximations(). Silently ignored when EW/profile. SearchControl() correctly omitted (scoring property, not search control). |
| T-159 | XPIWE Tests | G | Already implemented in c7a41712. 18 XPIWE-specific tests in test-ts-xpiwe.R: formula unit tests, 8-taxon missing-data scenario, TNT validation gallery (Vinther2008, Sano2011, Sansom2010 stored reference k-values). All pass. |
| T-160 | XPIWE Docs + NEWS | G | Rd docs and NEWS already in c7a41712. Added vignette paragraph to profile-scores.Rmd explaining XPIWE formula (eff_k, phi, extrapolation factor). Added 'cdot' to WORDLIST. spell_check_package() clean. Commit ea602512. |
| T-161 | XPIWE Shiny GUI | G | Added "Implied (extended)" as default step weighting in Shiny search config modal. Both "Implied (extended)" and "Implied" share concavity slider. `extendedIw()` reactive threaded through scores(), searchTask, StartSearch(). Updated shinytest2 snapshots. 42 search module tests + 13 Distribution tests pass. Commit 6da9a861. |
| T-162 | XPIWE Shiny citation | G | Added Goloboff 2014 citation to global.R and references panel (Tree Search section). Always shown since XPIWE is default. Commit a553a325. |
| T-184 | maxTime → maxSeconds alias | G | Already implemented in commit fafd5d0e. Intercepts maxTime before Morphy detection, maps to maxSeconds with .Deprecated() warning, removes from dots. Removed maxTime from .morphyParams list. Verified working. |
| T-163 | Search confidence composite diagnostic | G | Replaced exp(-K) with tighter binomial bound (1-K/R)^R, falling back to exp(-K) when K==R. Added optional nTopologies/lastImprovedRep params (wired by T-164). Ruggedness warning when K/R < 0.3 and R >= 5. Limited independence flag when nTopologies==1. 58 search module tests pass. Commit 2d2115cb. |
| T-164 | Wire pool stats to Shiny search confidence | G | Added `count_at_best()` to TreePool. Initialized new DrivenResult fields in parallel path. Wired to Shiny: nTopologies=length(allTrees), lastImprovedRep from search attrs, reset on weighting/concavity/dataset change. 58 module tests pass. Commit 16c02dc7. |
| T-181 | Add 180-taxon dataset to benchmark suite | G | Added mbank_X30754 (180t, 425c, 11 states, 40% missing, 20.5% inapp) as large-tree benchmark tier. `LARGE_BENCHMARK_NAMES`, `load_large_benchmark_datasets()`, `benchmark_large()`. Commit adec48b6. |
| T-180 | Warm-start benchmark infrastructure | G | `bench_warmstart.R`: `compute_warmstart_tree()` (sprint→TBR optimum), `warmstart_run()` (single rep from warm start), `warmstart_benchmark()` (grid), `warmstart_summary()`. Isolates ratchet/drift escape from initial descent. Verified: Vinther2008 sprint→80, warm-start→79. Commit 13a019e3. |

---

## 2026-03-20

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-148 | Red-team: ParsSim log-space convolution (S-RED focus 8) | B | Fixed `.LogCumSumExp` NaN bug: `-Inf - (-Inf) = NaN` in IEEE 754 when both accumulator and new value are `-Inf`. Guard added; 7 new assertions pass. |
| T-165 | Shiny: reset run stats on concavity/weighting change | C | Added stat-reset to `observeEvent(input$concavity)` and `observeEvent(input$implied.weights)` in `mod_search.R`. Trees kept; hits/reps/bestScore cleared. 2 new tests (50 total pass). |
| T-163 | C++ cancel checks in MPT enumeration + final TBR polish | C | Root cause: a "continue" search fills the pool quickly, leading to MPT enumeration with no cancel checkpoint. Added `check_cancel()` between seeds in the MPT while-loop and `check_timeout()` after final TBR polish in `run_single_replicate`. 152 driven-search tests pass. |
| T-166 | Shiny crash: `tipLabels()` subscript out-of-bounds on empty `r$trees` | A | `mod_data.R#107`: added `if (!length(r$trees)) return(character(0L))` guard. Fixes both `DatasetMatchesTrees` (#126) and `UpdateActiveTrees` (#203) stack traces. Resolves T-169 as downstream consequence. |
| T-167 | Shiny: "cannot open file" warning spam in progress poll | A | `mod_search.R`: added `if (!file.exists(pf)) return()` after `invalidateLater(500)`. Observer now silently skips until C++ creates the file. |
| T-168 | Shiny: `LEFT == RIGHT` recycling warning from `ape::read.nexus` | A | Upstream `ape` bug (bracket-index comparison). Suppressed with `suppressWarnings()` at both `read.nexus` call sites in `mod_data.R`. |
| T-169 | Shiny: `updateSliderInput()` value outside `[min, max]` when trees → 0 | A | `mod_data.R#231`: guarded entire tree-range/count slider update block with `if (nTrees > 0L)`; controls hidden via `parentHide` anyway when no trees present. |
| T-170 | Shiny: "no data loaded" in Configure modal; search terminates after second data load | A | Root cause was T-166 crash propagating through `DatasetMatchesTrees()` inside `tryCatch` error handler; fixed as consequence of T-166. |
| T-171 | Shiny: `LengthAdded` matrix-recycling and unknown-scoring warnings | A | `mod_consensus.R`: `PolEscVal` reactive now guards with `setequal(tipLabels(), names(r$dataset))`; returns `NULL` when taxa sets differ (tree has superset taxa vs dataset). |
| T-172 | Shiny: run counter accumulates across concavity changes | A | Root cause: Configure modal created `sliderInput(value=1L)` on each open, sending `input$concavity=1` to server and firing reset observer spuriously. Fixed by initialising all modal inputs from current `input$*` values (`mod_search.R`). |
| T-173 | Shiny: stop-when-N probability always ~37% for any dataset when `targetHits=1` | A | `global.R`: `SearchConfidenceText()` now appends "— increase 'Stop when N runs hit best' for a tighter estimate" when K=R and R≤5, guiding user toward more informative settings. |
| T-164 | Shiny UX: confusing run count when reducing max_runs during continued search | A | `global.R`: `SearchConfidenceText()` now accepts `nSearches`; says "total runs across N searches" when N > 1. `mod_search.R`: tooltip updated to explain per-search vs cumulative distinction; `helpText` added below maxReplicates slider in Configure modal. |
| S-RED | Red-team focus 9: Wagner & addition trees | C | BUG FIXED: boundary-edge false positive in `wagner_edge_violates_constraint` (ts_wagner.cpp) and `regraft_violates_constraint` (ts_constraint.cpp). Both rejected the edge directly above the constraint clade for MUST_OUTSIDE elements (`is_ancestor_or_equal(cn,below)` true when `below==cn`). Fix: `&& below != cn` guard. Search quality improvement, no correctness impact. +2 Wagner tests (43 total). 152 driven-search + 18 constraint tests pass. |
| T-175 | Shiny: search progress indicators vanish after 3-4 seconds (worker startup) | A | `mod_search.R`: both `searchTask` and `profilePrepTask` result observers used `validation = function(e) req(FALSE)` to handle "task running" state. In Shiny 1.8+, `ExtendedTask$result()` throws class `c("shiny.silent.error", "shiny.output.progress", ...)` — NOT "validation". The wrong handler name caused `error = function(e) NULL` to catch every "still running" signal, triggering premature cleanup (notification removed, `r$searchInProgress = FALSE`, `DisplayTreeScores()` called) ~3–4 s into each search (when the worker process starts and status flips to "running"). Fix: rename both `validation =` to `shiny.silent.error =`. |
| T-176 | Shiny: misleading error + wrong dataset when data file uploaded to tree loader | A | `mod_data.R`: after all tree-load attempts fail, try `ReadTntAsPhyDat` then `ReadAsPhyDat` on the file. If either succeeds, set `r$dataset`, `r$chars`, `r$charNotes` directly and notify "No trees found — loaded N taxa and M characters as dataset"; `observeEvent(r$dataset)` clears incompatible trees. If neither succeeds, keep existing "Trees not in a recognized format" message. Fixes (1) misleading error wording and (2) search continuing on previously-loaded dataset. |
| T-174 | Shiny: spurious "Inferring tip labels from dataset" + console warning spam on tree file load | A | `mod_data.R`: (1) `readLines(tmpFile)` in NA/NaN retry branch wrapped in `suppressWarnings` — suppresses benign EOF warning leaking to console. (2) `withCallingHandlers` inside `ReadTntTree` tryCatch muffles "incomplete final line" EOF warnings before they reach the outer `warning` handler (which is for genuine TNT tip-label warnings only). Also added `error = function(e) NULL` handlers to both inner tryCatches (was bare `NULL`). |

---

## 2026-03-19

| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-151 | Shiny: dataset observer clears bundled trees (blank plot / 0 trees) | B | `observeEvent(r$dataset,…)` unconditionally cleared `r$allTrees` after `UpdateData()` loaded trees from same .nex; fix: guard with `HaveData() && !DatasetMatchesTrees()`. All 31 bundled datasets affected. 2 new regression tests; 11 mod-data + 1681 ts-* pass. |
| — | Fix inapplicable.Rmd vignette CSL URL (404) | A | Remote `raw.githubusercontent.com` CSL URL returned 404; changed to local `../inst/apa-old-doi-prefix.csl` (matching all other vignettes). Also un-claimed stale Agent C GUI issue in `issues.md`. |
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

## 2026-03-20
| ID | Description | Agent | Notes |
|----|-------------|-------|-------|
| T-149 | Profile MaddisonSlatkin → VTune hotspot analysis | A | logB/logPVec cache overhead dominates (53% DLL); skill updated |
| T-152 | LSEAccumulator long double→double + expl→exp | A | 26% avg speedup across k=3..5 boundary cases |
| T-153 | std::isfinite→NEG_INF compare in MaddisonSlatkin.cpp | A | bundled with T-152 |
| S-PROF | Threshold comment update + fpl_d/fpl_u batching in LogPVec | A | ~15% further speedup; k=3:27/k=4:19/k=5:13 new thresholds |

| T-154 | OAFlatMap for logB_cache/logPVec_cache in SolverT | A | Probe-layer OA map eliminates node ptr-chase; deque backing for logPVec stable refs; fixes latent LogB() dangling-ref UB; bd883459. Speedup vs post-fpl-batch baseline: k=3 n=27 ~4.4×, k=4 n=19 ~3.5×, k=5 n=13 ~12×. No formal A/B build (Windows AV makes double-build ~60 min); derived from S-PROF sweep timings. |

## 2026-03-22

| ID | Description | Agent | Notes |
| — | MaddisonSlatkin: test freeze fix (align pp-multistate with split_count gate) | A | Commit 4b95c37a |
| — | MaddisonSlatkin: Carter O(1) closed-form for k=2 + LSE lookup tables + buffer reuse | A | Commit 93851932. k=2 >100× speedup; k=3 ~5-19% |
| — | MaddisonSlatkin: threshold recalibration with correct bitmask encoding | A | Commit b8666825. OLD thresholds used wrong state encoding. k=3: 136→75, k=4: 100→50, k=5: 100→35 |
| — | MaddisonSlatkin: wall-clock time budget safety valve | A | Commit 536f7ff9. 2s chrono budget; every LogB/LogPVec call checks clock (~20ns overhead). Blowup returns NA + warning. 137 tests pass |
|----|-------------|-------|-------|
| — | Collapsed-region regraft merging + collapsed pool dedup | — | Goloboff (1996) approach: skip zero-length edges as clip candidates, skip interior collapsed regraft positions (boundary-only evaluation), diversity-aware pool eviction on ties, collapsed-topology pool dedup. 0% skip rate on standard morph datasets but improves pool diversity. 35b5ad99 |
| — | Strip dead CollapsedRegions code from hot path | — | region_id/n_regions never read by consumers; all callers now use compute_collapsed_flags() directly. a16373a7 |
| — | Collapsed dedup in parallel search path | — | ThreadSafePool::add_collapsed() + worker thread + fuse_round updated. 3648beb9 |
| — | Cross-replicate consensus constraint tightening | — | Opt-in consensusConstrain=TRUE: after ≥5 reps, lock pool strict consensus as topological constraints. Clears on new best score. build_constraint_from_bitsets(), extract_consensus_splits(). 09f69915 |
| — | Strategy preset tuning | — | default: wagnerStarts=3, sprFirst=TRUE, adaptiveLevel=TRUE. thorough: sprFirst=TRUE. Bundled in 09f69915 |

## 2026-03-23
| T-188 | Biased Wagner addition — API integration | Human+AI | `wagnerBias` (0=random, 1=Goloboff, 2=entropy) + `wagnerBiasTemp` in `SearchControl()` and `ts_driven_search`. First Wagner start uses biased order; remaining starts use random for diversity. Goloboff 2014 §3.3. Benchmarked: 80% Wagner→TBR gap reduction at 174t; marginal on ≤88t. |
| T-189 | Outer search cycle loop | Human+AI | `outerCycles` param in `SearchControl()`. Wraps [XSS+RSS+CSS → Ratchet → NNI-perturb → Drift → TBR] in outer loop, distributing cycles evenly. Matches TNT xmult pattern (Goloboff 1999 §2.3). `thorough` preset defaults to `outerCycles=2`. Backward-compatible: default=1. |
| T-184 | `maxTime` → `maxSeconds` alias | Human+AI | Already implemented in b8e56e2b. `maxTime` intercepted before `.morphyParams` check, mapped to `maxSeconds` with deprecation warning, routed to C++ engine. Verified: timings attribute present, Morphy() not called. |

## 2026-03-24 (evening)
| T-209 | NNI perturbation constraint guard | E | Gate on `(!cd || !cd->active)`. PR #220 merged. |
| T-206 | Outer cycle reset cap / maxOuterResets | E+A | Cap resets to maxOuterResets (default 3). PR #218 merged. |
| T-210 | SA best-tree tracking across phases | C | `anneal_search()` saves/restores best tree at phase boundaries. Commit `e204d0a0` on feature/pt-eval. |
| T-182 | Adaptive ratchet perturbation probability | G+E+A | hitRate-based tapering. PR #221 created. |

## 2026-03-25
| T-208 | random_topology_tree ignores constraints | G+A | WAGNER_RANDOM fallback. Cherry-picked to cpp-search (24427c9a). PR #219 closed. |
| T-211 | Stale final_ in temper candidate scoring | C | Closed: conservative-only impact, not worth fixing. |

## 2026-03-25 (morning)
| T-215 | cli progress bar `::` resolution fix | A | `pb_env` parent `baseenv()` → `environment()` in `MaximizeParsimony()`. Commit 908860d25. |
| T-216 | Shiny app `"brazeau"` → `"bgs"` | A | 8 occurrences in `mod_search.R` (comparisons, defaults, selectInput value). Commit 908860d25. |
| T-217 | `tree = NULL` in `MaximizeParsimony()` Morphy path | A | Added `!is.null(tree)` guard on Morphy delegation (line 485). Main path already correct. Commit 908860d25. |
| T-218 | Simplification transforms corrupt NA scoring | A | Commit `a48bfc4ad` let `?` tokens through Transforms 2/3 (not NA-safe). Fix: revert transform bypass, add conservative constant-char removal inside Phase 1 for `?`-only chars. Supersedes D's simple revert PR #224. Commit `08054102f`. |
| T-221 | [Shiny] Crash loop in cluster consensus concordance | B | `LabelConcordance()` guard `!is.null()` → `inherits(, "phylo")`. Commit `bc5313c22`. |
| T-222 | [Shiny] "Align tips" does nothing in Characters on trees | B | `Display` callback always set edge.length=1; now NULL when tipsRight checked. Commit `b23580823`. |
TASKEOF 2>&1
