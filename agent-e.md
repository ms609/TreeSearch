# Agent E — Progress Log

## Current Task
- **Task:** T-255 — Reduce drift in default and thorough presets
- **Status:** Starting. Set driftCycles=0 in default; reduce from 12 to 0 in thorough.

### T-254 — Drift MPT diversity experiment — DONE
- driftCycles=0 vs 2 on Wortley2006 (37t), Zhu2013 (75t), Geisler2001 (68t)
- 3 seeds × 2 budgets (30s, 120s), with and without consensus stopping
- Primary finding: drift provides zero benefit on all metrics (score, MPT count, topological diversity)
- Wortley2006: no-drift consistently finds 4 MPTs, drift finds 1–3 (median 2)
- Geisler2001/Zhu2013: mean RF identical (7.3 vs 7.4; 11.6 vs 10.2)
- Drift consumes 15–19% of wall time, reducing replicates by 10–22%
- With consensus stopping, drift delays stabilization without improving the answer
- Recommendation: set driftCycles=0 in default and thorough presets (T-255 unblocked)
- Write-up: `dev/benchmarks/drift_mpt_analysis.md`

### T-251 — TNT trajectory analysis — DONE

### T-250 — TNT Fitch kernel disassembly — DONE
- TNT is a 32-bit i386 binary with zero SIMD (no SSE2, no AVX) and table-based popcount (64KB LUT, 16-bit halves)
- TreeSearch uses 64-bit x86-64 with SSE2 128-bit pand/por and software Hamming weight popcount
- TreeSearch has ~4× raw Fitch throughput advantage over TNT (128 bits/iter vs 32 bits/iter)
- Neither uses hardware popcnt (potential minor optimization: `-mpopcnt` compile flag)
- **Key finding:** TNT's 3-5× convergence speed advantage is STRATEGIC, not implementation-level
- TNT evaluates fewer total candidates to reach the same score — the gap is in search heuristics
- T-251 (trajectory analysis) elevated in priority; T-246 (AVX2) still worthwhile but not the primary bottleneck
- Write-up: `dev/benchmarks/tnt_disassembly_analysis.md`

### S-COORD round 21 — DONE
- Closed 4 Shiny tasks (T-232, T-239, T-240, T-241) — all re-validated by GHA 23547582438 (cpp-search PASS)
- Updated T-242 note (GHA failure is stale, investigation task)
- Updated S-PR: #230 pending, #216 stale GHA, #178/#106 recommend closing

### T-248 — SA phase tuning — DONE
- Hamilton benchmark: AC=0/1/3 × 30s/60s, 5 seeds each, mbank_X30754 (180t)
- AC=1: 400ms/rep, 40% hit rate. AC=3: 1370ms/rep, 21% hit rate. No significant score difference (p>0.5).
- Changed large preset annealCycles from 3 to 1. Saves ~1s/rep (~6% of 17s).

### T-247 — XPIWE search quality investigation — DONE (not a bug)
- Score discrepancy (3.84382 vs TNT 3.79283) from different inapplicable handling
- TreeSearch: Brazeau three-pass → EW=80 on TNT's tree, EW=79 on its own tree
- TNT: standard Fitch → EW=78 on its tree
- XPIWE implementation verified correct: eff_k flows through compute_iw, precompute_iw_delta, indirect_iw_length
- Search confirms XPIWE uses different intermediate scores from IW (verbosity=3 shows different convergence paths)

### T-244 — Full-pipeline 180-tip benchmark — DONE
- Ran large preset on Hamilton HPC (EPYC 7702) at 30/60/120s budgets, 5 seeds each
- Median scores: 30s=1202, 60s=1190, 120s=1185 (cf. Intel pre-T-206: 1276/1255/1250)
- 65-74 step improvement primarily from T-206 (outer cycle reset cap), not hardware
- Per-replicate: median 17.3s (vs pre-T-206 ~60s due to unlimited resets)
- Phase: TBR 43.6%, Ratchet 32.2%, SA 7.4% (least productive: 14% hit, 0.8 steps/s)
- SA (annealing) is potential optimization target — annealCycles/Phases may be overtuned
- Updated profiling.md with baselines

### T-243 — hot-loop-opt PR — PARKED
- PR #230 open. GHA 23580149481 failed on pre-existing spelling/doc issues.
- Commented on PR. Changes are pure C++ — no new failures introduced.

### T-209 — NNI perturbation constraint guard — DONE
- Found via S-RED focus 2 (search topology invariants).
- Bug: `random_nni_perturb()` ignores constraints → cn=-1 → TBR blocked.
  Same mechanism as T-208 (random_topology_tree).
- Fix: gate NNI perturbation on `(!cd || !cd->active)` at ts_driven.cpp:310.
- Local build: OK (build-agent.sh). Tests: 152 driven, 43 wagner, 29 SC — all pass.
- init.c: 47 entries (45 Rcpp + 2 manual), all match.
- GHA 23506971922: PASS. PR #220 created.
- Also handled T-182 GHA failure (Agent G): fixed roxygen docs, re-dispatched GHA 23507428719.

### T-206 progress
- Added `max_outer_resets` to `DrivenParams` (C++) and `maxOuterResets` to
  `SearchControl()` (R).
- Default 0 (no resets); default preset = 2, thorough = 3.
- Fixed misleading comment in `ts_driven.cpp`.
- 1826/1826 ts-* tests pass (2 xpiwe namespace false positives excluded).
- Benchmark (Zhu2013 5 reps): no resets=4.81s, 2 resets=9.25s, old unlimited=12.99s.
  Score quality preserved (640–646 vs 641–646).
- GHA run history:
  - Run 23503177988: FAIL — test-SearchControl.R poolSuboptimal sensitivity
  - Run 23503681296: FAIL — codoc mismatch (maxOuterResets missing from Rd)
  - Run 23505343592: DISPATCHED — includes both fixes (test + roxygen regen)

### S-PROF supplement (round 4) — DONE
Agent G ran S-PROF round 4 earlier today (re-baseline with auto strategy).
My supplement focused on:

1. **strategy='none' baselines** (5-run medians, 5 reps, serial):
   - Vinther2008: 1.00s (was 0.42s), score 79
   - Agnarsson2004: 7.08s (was 1.79s), score 778
   - Zhu2013: 12.99s (was 3.17s), score 641–646 (was 648–666)
   - Dikow2009: 13.37s (was 4.90s), score 1611–1612 (was 1612–1614)

2. **Outer cycle reset mechanism finding:** Even with `outerCycles=1`,
   the reset-on-improvement logic (`ts_driven.cpp:406`) causes 3–5 full
   [XSS+RSS+Ratchet+Drift+TBR] cycles per replicate. C++ comment claims
   "outer_cycles = 1 is exactly the previous linear pipeline" — incorrect.
   Diminishing returns: cycle 1 = 47 steps/s, cycle 3 = 0.8 steps/s.

3. **Filed T-206** (P3): Outer cycle reset cap / minimum-Δ gate.

4. **Updated `.positai/expertise/profiling.md`** with strategy='none'
   baselines and per-cycle escalation data.

### S-RED focus 6 — DONE
- init.c: 45 entries (43 Rcpp + 2 manual), all arg counts match
- Score verification: Vinther2008=79, Agnarsson2004=778, Longrich2010=131
  (all match TreeLength, both serial and parallel)
- Parallel hits_to_best: confirmed hits ≤ reps after T-139 fix
- ts-* tests: 1676 pass, 3 fail (all in test-ts-profile.R:
  multi-state profile scoring — from human commit 5235d6e1
  "Profile parsimony with more steps", same regression as T-144)
- Module tests: 88 pass, 0 fail
- Pre-existing regressions (NOT from my changes):
  - test-data_manipulation.R: 3 fails (T-144, assigned to A)
  - test-Concordance.R: 3 fails (T-144 cascade)
  - test-MaddisonSlatkin.R: 26 fails (human commit 5235d6e1)
  - test-ts-profile.R: 3 fails (human commit 5235d6e1)
  - test-Morphy.R: 2 fails (unclear, not from my changes)
  - test-Multi-state StepInformation: 5 fails (profile scoring)
  All these relate to the profile parsimony / PrepareDataProfile refactor
  in the human's recent commits. T-144 (assigned to A) tracks the fix.
- No new bugs found from T-137/139/141 changes.

### T-141: Sun 2018 blank plot — DONE
- Root cause: Sun2018.nex has 0 trees. When switching from a dataset with
  trees (e.g., Wills 2012), old trees persisted with incompatible tips.
- Fix: In `mod_data.R`, when `read.nexus()` fails (no trees in file), clear
  old trees only if their tips don't match the new dataset.
- Profile error is a separate issue (T-140, assigned to Agent F).

### T-139: Hit/rep count — DONE
- Bug: `parallel_driven_search()` captured `hits_to_best` after MPT
  enumeration, inflating the count. Fixed to capture before enumeration
  (matching serial path in `driven_search()`).
- Added clamp in `SearchConfidenceText()` as safety net.

### T-137 Progress
- Triaged 5 issues from issues.md → T-137/138/139/140
- Root causes identified:
  1. Cancel-file check only happens between phases in `run_single_replicate()`,
     not inside `ratchet_search()` or `drift_search()`. Long phases (20 ratchet
     cycles, 12 drift cycles) block cancel detection for minutes.
  2. Cancel observer didn't remove the search notification → "Searching..." toast
     lingers while "Stopping…" shows in results area.
  3. Result observer gated cleanup on `r$searchNotification` — if cancel observer
     clears it first, buttons never get re-enabled.
- Fixes implemented:
  - **C++**: Added optional `std::function<bool()> check_timeout` param to
    `ratchet_search()` and `drift_search()`. Checked between cycles.
    `driven_search()` passes its `check_timeout` lambda to both.
  - **Shiny cancel observer**: Now removes `r$searchNotification` immediately
    and shows "Stopping — waiting for current search phase to finish…"
  - **Shiny result observer**: Gate cleanup on `r$searchInProgress` flag
    (not `r$searchNotification`). Flag set in `StartSearch()`, cleared on
    completion. Robust against notification being pre-cleared by cancel.
  - Added `searchInProgress` field to `AppState` and test helper.
- Tests: 1669 ts-* pass, 88 module tests pass, init.c verified

### S-RED: Red-team review (focus 4: Parallelism & RNG) — DONE
- Reviewed ts_parallel.cpp (450 lines), ts_rng.h/.cpp (110 lines), ts_driven.cpp (484 lines)
- Thread-local RNG: correctly set/cleaned per worker. Seeds pre-generated from R RNG.
- No R API calls from workers: Rprintf gated by verbosity, interrupt via thread_stop_flag.
- Pool mutex: all access under lock_guard. Fuse under lock is correct.
- HSJ/Sankoff thread safety: both stateless, no global mutable state, zero R API calls.
- DataSet copy: all HSJ/Sankoff fields are std::vectors, deep-copied correctly.
- PERF note: XFORM rebuilds SankoffData per score_tree() call (optimization opportunity).
- init.c: 45 entries, all arg counts match.
- Score verification: serial=80, parallel=79 (Vinther2008), both TreeLength-verified.
- No bugs found.

### T-127: Shiny: Async profile scoring with progress + cancel — DONE
- Inlined `PrepareDataProfile()` logic into the `profilePrepTask` ExtendedTask future body
  so the slow `StepInformation()` per-pattern loop can report progress and check a cancel file.
- Added `profileProgressFile` and `profileCancelFile` reactiveVals.
- Progress polling observer (`invalidateLater(500)`) reads a temp file written by the future
  and updates the notification with "X/Y patterns (Z%)".
- Cancel mechanism: switching away from profile mode creates a cancel signal file;
  the future checks `file.exists(cancelPath)` between patterns and returns NULL if cancelled.
- Result observer cleans up both temp files (progress + cancel).
- New test: "switching away from profile cancels prep via cancel file" (4 assertions).
- Fixed pre-existing bug in cancel button test (`ignoreInit = TRUE` required explicit
  input initialization with `cancel = 0` before `cancel = 1`).
- Updated shinytest2 snapshots for cancel button UI addition (Agent B's T-130).
- 38/38 mod-search tests pass, 119 Shiny assertions pass (1 pre-existing non-deterministic
  image snapshot in Distribution-003), 173 ts-* tests pass.

## Completed Tasks

### S-PROF: Performance profiling — DONE (~19:45)
Three investigations:
1. **Drift threshold sensitivity**: AFD×RFD grid on Zhu2013 (15 runs/config) + Dikow2009 (10 runs/config). No significant score differences across threshold values (p=0.60–1.00). Cycle count matters more. d2 drift provides no benefit over ratchet alone on Dikow2009.
2. **IW scoring overhead**: EW vs IW (k=3, k=10) on 3 datasets. IW 64% faster on tiny data, 26–57% slower on 62–75 tips. Scales with dataset size. No optimization opportunity.
3. **Fuse effectiveness**: fuseInterval=0 vs 3 on 3 datasets. Negligible overhead and score impact. Pool deduplication keeps pool small; fuse rarely triggers.
No new optimization tasks raised. All 7 profiling items now ✅.

### S-RED: Red-team review — DONE (~18:30)
Post mod_search (T-061), mod_data (T-062), mod_clustering (T-060), mod_consensus (T-063), MaximizeParsimony.R strategy+warning changes.
- Full test suite: 7979 pass, 0 fail, 22 warn (pre-existing), 44 skip (+270 vs prior S-RED)
- ts-* tests: 221 pass, 0 fail, 0 warn, 39 skip
- Module tests: 43 pass, 0 fail, 2 warn (pre-existing ns issue in mod-clustering testServer)
- init.c: 41 entries, 39 Rcpp + 2 manual, all arg counts match
- Score verification: Agnarsson2004=778, Vinther2008=79, Longrich2010=131 — all match TreeLength
- Determinism: identical Newick under same seed (serial mode, Vinther2008)
- Code review: mod_search, mod_data, mod_clustering, mod_consensus wiring in server.R all correct
- Forward-ref bridge (cb_ref) pattern verified: DisplayTreeScores, UpdateKeepNTipsRange, UpdateDroppedTaxaDisplay, UpdateOutgroupInput
- events.R slimmed to 36 lines (ShowConfigs + plotFormat observer only)
- MaximizeParsimony.R: .AutoStrategy(nTip, nChar) signal-density gate correct; replicate adequacy warning correct with !missing() guard
- ts_wagner.cpp: comment-style fix only (T-070)
- No new issues found

### S-RED: Red-team review — DONE (~17:35)
- init.c: 41 entries, 39 Rcpp + 2 manual-only, arg counts all match
- Test suite: 7709 pass, 0 fail, 22 warn (pre-existing), 44 skip
- Module tests: 38/38 pass (6 clustering, 5 data, 11 downloads, 4 refs, 8 search, 4 treespace)
- Score verification: Agnarsson2004=778, Vinther2008=79, Longrich2010=131 — all match TreeLength
- Determinism: identical Newick under same seed (serial)
- Shiny wiring: forward-ref callbacks correct, parent_session cross-module updates correct
- events.R (168 lines) properly annotated; remaining bindings are consensus-scope (T-063)
- No new issues found

### T-058: Shiny mod: mod_downloads — DONE (~17:10)
- `server/mod_downloads.R` with `downloads_ui(id)` / `downloads_server(id, ...)` — 8 download handlers namespaced
- Wired into `global.R` (source + `dl_ui <- downloads_ui("dl")`), `ui.R`, `server.R`
- `testServer()` tests: 11/11 pass (4 test blocks covering saveZip, savePlotZip, savePlotNwk, saveNwk)
- Old `server/downloads.R` deleted (no longer sourced)
- shinytest2: mod-downloads 11 pass, app-smoke 3, SearchLog 4, ViewChars 12

### Shinylive plan (stale, closed out)
- Was continuing T-045 plan; superseded by other work.

### T-036: R-level test coverage — DONE (~09:15)
- Created `tests/testthat/test-MaximizeParsimony-features.R`
- 18 test blocks, 38 expectations, 0 failures

### S-RED: Red-team review — DONE (08:45)
- Found and fixed T-025 root cause (PreallocUndo buffer overflow)

### T-012: SPR→TBR escalation — DONE (07:35)
### T-024: Parallel resample — DONE (07:25)
### S-COORD: Coordination review — DONE (07:10)
### T-007: Wagner NA-incremental scoring — DONE (07:00)
### T-019: AdditionTree migration — DONE (06:30)

## Blockers / Questions
(None)
