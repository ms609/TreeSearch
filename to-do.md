# TreeSearch Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (d1)` (or `d2`,
  `d3`, … — ephemeral dispatcher IDs issued by `dispatch.sh`).
- When a task is being developed in a **git worktree**, set its status to
  `WORKTREE (name)` where *name* is the worktree directory (e.g.
  `WORKTREE (TS-CID-cons)`). This distinguishes human/long-running worktree
  work from agent assignments and prevents double-claiming.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md` (see workflow in AGENTS.md).
- Tasks awaiting GHA results: `PARKED (d1, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (d1)`.
  S-COORD cleans these up after merge.
- The `Notes` column may include a bracketed model/effort hint, e.g.
  `[m:haiku e:low]`, `[m:sonnet e:medium]`, `[m:opus e:high]`. The
  dispatcher's ranker honours these hints and they override its default choice.
- Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
  completed, reset it to OPEN. Their effective priority is dynamic:
  - ≥6 OPEN specific tasks → standing tasks are P3
  - 3–5 OPEN specific tasks → standing tasks are P2
  - <3 OPEN specific tasks → standing tasks are P1

---

## Active Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-150 | P2 | PR #213 (F) | — | **CID-optimal consensus tree search** | PR #213. Vignette fix (TreeTools::Consensus) commit f8bfee49. GHA 23650002703. |
| T-204 | P2 | PR #216 (F) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | GHA 23649607006 PASSED. Ready for merge. |


### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-303 | P2 | PR #247 (t303) | — | **Sector heuristic degrades silently on HSJ/XFORM datasets.** `build_reduced_dataset` in `src/ts_sector.cpp:421-444` does not copy `hierarchy_blocks`, `tip_labels`, `n_orig_chars`, `hsj_alpha`, or `sankoff_*` fields. `rd.data.scoring_mode` IS copied, so internal `score_tree(rd.subtree, rd.data)` dispatches `hsj_score()`/Sankoff with empty hierarchy/Sankoff data, degrading to Fitch-only. Final acceptance scores correct (use full ds) — only the sector's internal accept/reject heuristic is wrong → missed improvements, possible accept-then-revert churn. Same class as T-275 guard. **Fix (PR #247):** rss/xss already guarded on cpp-search (e5ff2942, approach a); css_search unaffected (scores full ds, no reduced dataset) — documented + sectorial HSJ regression test added. Approach (b) intractable: HTU pseudo-tip has no valid HSJ tip_labels / Sankoff tip_costs. | Found by /red-team area 5 (2026-05-26). PROFILE+IW are fine. [m:sonnet e:medium] |
| T-304 | P2 | PR #248 (t304) | — | **T-300 dirty-set rescore has no enduring regression test.** EW+NA dirty-set rescores wired into `tbr_search` SPR accept path (`src/ts_tbr.cpp:1138-1180`). The `DEBUG_RESCORE`/`DEBUG_NA_RESCORE`/`DEBUG_NNI_RESCORE` cross-checks that validated them were fully removed in 5b210fdd, 44a4ebeb, 2be8228d. Previous incremental attempt regressed with systematic delta=-3 (b7303ee5 revert). Need a Tier-2 test driving many SPR accepts (small n, weak signal, many maxHits) asserting `result$score == ts_score(result_tree, ds)` across EW/IW/NA/NA-IW. | Found by /red-team area 8 (2026-05-26). Pattern: see `test-ts-spr-state-restore.R`. [m:sonnet e:medium] |
| T-306 | P3 | PR #249 (t306) | — | **HSJ/XFORM SPR/NNI accept-paths omit hierarchy DP contribution from `best_score`.** In `tbr_search` SPR accept (`src/ts_tbr.cpp:1146-1180`) and `nni_search` accept (`src/ts_search.cpp:79-95`), `best_score` is updated as Fitch-only delta (EW: `best_score + delta`; IW: `compute_weighted_score`). Neither calls `hsj_score()` nor adds Sankoff. For HSJ/XFORM modes (`use_iw = false` since concavity is HUGE_VAL), `best_score` therefore tracks Fitch+ew_offset only, not the topology-dependent `hsj_total`/Sankoff added by `score_tree`. Pre-T-300 the SPR path called `full_rescore` and was internally correct (but candidate evaluation in Phase 1 was already Fitch-only — a deeper structural issue: candidates aren't compared on full HSJ score, so accept/reject decisions never see hierarchy DP). User-visible scores remain correct because `run_single_replicate` always recomputes via `score_tree(tree, ds)` before pushing to pool (`ts_driven.cpp:181,247,259,...,595`). Search-quality regression only — missed/wrong accepts. Fix: gate dirty-set + delta path behind `ds.scoring_mode` being `EW`/`IW`/`PROFILE`/`XPIWE`, falling back to `full_rescore` for HSJ/XFORM. Even better: include `hsj_score`/Sankoff delta in candidate evaluation (broader fix). | Found by /red-team area 1 (2026-05-26). Empirical test on 15-tip HSJ dataset showed no user-visible score mismatch (final score recomputed via `score_tree`); search-quality impact is silent. Related: T-303 (sector path same family). [m:opus e:high] |
| T-322 | P3 | OPEN | — | **Wagner NA+IW regression test is tautological (omits `min_steps`).** `tests/testthat/test-ts-wagner.R:223-242` — the test "Wagner on NA + IW matches fitch_score" calls `TreeSearch:::ts_random_wagner_tree(...)` and `TreeSearch:::ts_fitch_score(...)` both with `concavity = k` but **omits `min_steps`** (defaults to `integer(0)`). The implied-weight homoplasy `h = steps − min_steps` is thus computed as `h = steps − 0` on both sides, so the cross-check (Wagner incremental score == independent Fitch rescore of the same tree) passes while validating a *non-production* formula. The real NA+IW path (`R/MaximizeParsimony.R:834`) always passes `min_steps = as.integer(MinimumLength(ds, compress = TRUE))`; Vinther2008 carries inapplicable characters so `MinimumLength` is non-zero and the tested formula genuinely differs from production. A regression in NA+IW `min_steps` handling would pass this test undetected. Fix: add `min_steps = as.integer(MinimumLength(pd, compress = TRUE))` to **both** calls (signature accepts it, RcppExports.R:147), re-run to confirm the cross-check still holds (same `min_steps` both sides → still valid, now exercises production scoring). Filed not fixed inline because it changes test numerics and must be verified by a test-run — and may itself surface a latent wagner NA+IW bug. | Found by /red-team area 8 (2026-06-16). Verified REAL (sonnet). Cross-links area 9 (Wagner). [m:sonnet e:low] |


### Shiny App

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-309 | P2 | OPEN | — | **EasyTrees: stale profile dataset scores wrong trees.** `inst/Parsimony/server/mod_search.R:440` — on `profilePrepTask` completion the code runs `profileDataHash(r$dataHash)`, stamping the *current* dataset hash at completion time, not the hash of the dataset that was actually prepared. Load dataset H2 while profile prep runs on H1 → completion sets `profileDataHash=hash(H2)` while `profileDataset=preparedFrom(H1)`; `StartSearch()` guard (`:640`, `identical(r$dataHash, profileDataHash())`) then skips re-prep, and `scores()` (`:475`, no hash check) scores the H2 search trees against the H1-derived profile dataset → researcher sees profile scores from the wrong dataset. `observeEvent(r$dataset)` (`:1128`) resets search stats but never clears `profileDataset()`/`profileDataHash()`. Fix: stamp `profileDataHash()` with the hash of the prepared dataset (snapshot at invoke time) and clear `profileDataset(NULL)`/`profileDataHash(NULL)` on data change. | Found by /red-team area 7 (2026-06-16). Verified REAL (opus). Data-integrity (publishable wrong numbers) but needs a mid-prep data swap. [m:sonnet e:medium] |
| T-310 | P2 | OPEN | — | **EasyTrees double-launch: no `searchInProgress` guard in `StartSearch()`.** `inst/Parsimony/server/mod_search.R:632` lacks a re-entrancy guard. `shinyjs::disable("go")` is an async browser round-trip, so a fast double-click fires `observeEvent(input$go, StartSearch())` twice. Verified vs shiny 1.13.0 `ExtendedTask` source: `invoke()` while running *queues* the second call. The 2nd `StartSearch()` overwrites the single `cancelFile()`/`progressFile()` reactiveVals and `r$searchNotification` (leaks the 1st notification — no `removeNotification` at `:719`) and re-enables Go mid-flight; when task 1 settles the result observer may delete task 2's signal files, or hit `searchTask$result()`→`req(FALSE)` and silently drop task 1's trees. Fix: `if (isTRUE(r$searchInProgress)) return(invisible())` at the top of `StartSearch()`. | Found by /red-team area 7 (2026-06-16). Verified REAL (opus). One-line fix. [m:sonnet e:low] |
| T-311 | P3 | OPEN | — | **EasyTrees: session disconnect never cancels the running search worker.** `inst/Parsimony/server.R:187` `onStop` cleans only file caches + cmd log; it never writes the `cancelFile()` signal the `future::future()` worker polls (`mod_search.R:710`). A user who disconnects mid-search leaves the worker consuming a core until it finishes its replicates or hits the timeout (up to ~60 min for "thorough"). Fix: write the active cancel signal in `onStop` (or expose a module `cancel()` for `server.R` to call). | Found by /red-team area 7 (2026-06-16). Verified REAL (haiku). [m:haiku e:low] |
| T-312 | P3 | OPEN | — | **EasyTrees: search temp files (`ts_*`) leak on session end.** `inst/Parsimony/server.R:192-194` — `onStop`'s `unlink(... pattern="^(data\|tree\|excel)File-")` does not match the temp files `mod_search.R` creates: `ts_cancel_*`, `ts_progress_*`, `ts_profile_prog_*`, `ts_profile_cancel_*`. The worker `on.exit` clears some on the normal path, but on error/interrupt/disconnect they accumulate in `tempdir()` (the documented "Issue 6" tempdir growth in `.positai/expertise/shiny-app.md`). Fix: add `unlink(list.files(tempdir(), pattern="^ts_(cancel\|progress\|profile_prog\|profile_cancel)_", full.names=TRUE))` to `onStop`. | Found by /red-team area 7 (2026-06-16). Verified REAL (haiku). [m:haiku e:low] |
| T-313 | P3 | OPEN | — | **EasyTrees: topology dedup includes branch lengths → inflated tree pool.** `inst/Parsimony/server/mod_search.R:1063-1066` — the "topology string" dedup uses `write.tree(ape::ladderize(t))`, but `write.tree()` serialises branch lengths when present. After `combined <- c(r$allTrees, newTrees)` mixes user-loaded trees (which may carry BLs) with parsimony trees (no BLs), topologically identical trees with different BLs are not deduplicated, inflating the pool and the displayed tree count. Fix: strip branch lengths before serialising (drop `$edge.length`, or use a topology-only key). | Found by /red-team area 7 (2026-06-16). Verified REAL (haiku). [m:haiku e:low] |


### Alternative Homologies (Goloboff 2026) — `feature/alt-homology` / `TS-AltHom`

Ref: Goloboff (2026) *Cladistics* doi:10.1111/cla.70033.
Plan: `dev/plans/2026-03-27-1415-implement-goloboff-2026-alternative-homologies-with-step-matrix-recoding.md`

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-280 | P3 | OPEN | — | **AltHom Phase 1: `AlternativeHomology` S3 class & core recoding (MVP).** Create `R/AlternativeHomology.R` (constructor, validation, print), `R/recode_alt_homology.R` (correspondence enumeration, morphotype states, cost matrix, tip assignment). Wire into `TreeLength()` for scoring on a fixed tree. Reproduce paper's Definition 1 cost matrix + Table 1 as tests. | WORKTREE (TS-AltHom). Invertible, no external constraints, two part-types only. |
| T-281 | P3 | OPEN | T-280 | **AltHom Phase 2: Constraints & options.** Non-invertible (`>`), adjacent (`>>`), restricted homology (`!`), configurable part transformation costs, adjacent-loss merging (`<`). Reproduce Definitions 2–3 and their cost matrices. | WORKTREE (TS-AltHom). |
| T-282 | P3 | OPEN | T-280 | **AltHom Phase 3: Wire into `MaximizeParsimony()` search pipeline.** Accept `AlternativeHomology` in `hierarchy` param, prepare xformArgs, end-to-end search. Also wire `Resample()` and `SuccessiveApproximations()`. | WORKTREE (TS-AltHom). |
| T-283 | P3 | OPEN | T-280 | **AltHom Phase 4: External inapplicability.** An external character can make individual characters, parts, or entire part sets inapplicable. Expand state enumeration for externally-disabled states. | WORKTREE (TS-AltHom). |
| T-284 | P3 | OPEN | T-280 | **AltHom Phase 5: Combination pruning.** Implement `xlinks&` (pairwise compatibility), `xlinks!` (observed-state-only), `xlinks@` (uninformative-state restriction) to reduce supercharacter state count. Verify same optimal trees as unpruned. | WORKTREE (TS-AltHom). |
| T-285 | P3 | OPEN | T-280 | **AltHom Phase 6: Implied weighting support.** Compute combined minimum steps across all valid alignments (not sum of per-char minima). Required for correct IW homoplasy counts. | WORKTREE (TS-AltHom). |
| T-286 | P3 | OPEN | T-280 | **AltHom Phase 7: Mixed `AlternativeHomology` + `CharacterHierarchy`.** Support datasets with both simple hierarchy blocks and alternative homology blocks in one analysis. | WORKTREE (TS-AltHom). |
| T-287 | P3 | OPEN | T-284 | **AltHom Phase 8: Static alignment fallback.** For datasets where supercharacter exceeds practical state limit, generate alternative static datasets (one per alignment) and search each. | WORKTREE (TS-AltHom). |
| T-288 | P3 | OPEN | T-282 | **AltHom Phase 9: Documentation & vignette.** `vignettes/alternative-homologies.Rmd`, roxygen docs for all new exports, `inst/REFERENCES.bib` entry. | WORKTREE (TS-AltHom). |

### Deferred / Future Directions

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-290 | — | DEFERRED | — | **GPU-accelerated batch tree scoring.** Evaluate many TBR/SPR candidate rearrangements in a single GPU kernel launch (parallelism across *trees*, not within one tree). For a 180-leaf tree the TBR neighborhood is O(n³) ≈ millions of candidates — enough to saturate GPU hardware. Main challenges: (1) per-candidate work is tiny for Fitch+bitwise (~50 word ops), so GPU arithmetic intensity is very low; (2) tree data structures need flat-array redesign for coalesced GPU memory access; (3) for morphological data sizes (≤500 chars, k ≤ 10) CPU OpenMP parallelism across candidates likely captures most of the win with far less effort. GPU becomes more compelling for Sankoff/implied-weights scoring (O(k²) per node per char) or phylogenomic-scale data (10k+ chars). A hybrid design (CPU manages search logic, GPU batch-scores candidates) is more practical than porting the full search engine to CUDA. **References:** Santander-Jiménez et al. (2020) *J Supercomput* 76:9827 (GPU Fitch parsimony, Kepler→Turing); Santander-Jiménez & Vega-Rodríguez (2025) *Future Gen Comput Syst* (OpenMP/OpenACC/SYCL multi-platform parsimony scoring); Ayres et al. (2019) *Syst Biol* 68:1052 (BEAGLE 3 — GPU likelihood, architectural lessons). | Research: MkPrime `.agent-d.md` 2026-03-29. |
| T-291 | — | DEFERRED | — | **GPU-parallel independent search replicates.** Run 100+ search replicates simultaneously on GPU SMs (one replicate per SM; modern GPUs have 60–128 SMs). Shared read-only character matrix fits in GPU L2 cache. Main obstacle: tree search has highly irregular, data-dependent control flow (rearrangement selection, acceptance decisions, ratchet perturbation) which causes warp divergence and poor GPU utilization. Branch-and-bound in sectorial search has the same problem. CPU multicore parallelism (8–16 cores via `future`/`parallel::mclapply`, or 100+ via HPC SLURM array jobs) is far simpler and more efficient per-replicate. GPU replicates only become attractive if per-replicate arithmetic is heavy enough to dominate over control flow overhead (e.g., large Sankoff matrices). **References:** same as T-290. | Research: MkPrime `.agent-d.md` 2026-03-29. |

### TNT Comparison & Strategy Learning

### Strategy Tuning


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-298 | P3 | PR #242 (d8) | — | **Profile and optimize `quartet_concordance.cpp` matrix allocation** | GHA 25777319791. Resize-hoist committed, benchmarked, PR open. |
| T-300 | P3 | PARKED (d7, other ) | — | **Lazy `apply_tbr_move` rescore in `tbr_search`.** After the `score_fresh` flag was wired (companion to T-187/PSF work), the trailing `full_rescore` at function exit is now skipped when states are coherent. The remaining redundancy is the `full_rescore(tree, ds)` call at `ts_tbr.cpp:1134`, run after **every** successful `apply_tbr_move` to obtain the authoritative score for the acceptance check. Each call is O(n_node × total_words). Since the move is local (clip + reroot + regraft), the indirect-evaluation pre-check at `ts_tbr.cpp:767-772` already shows that `fitch_incremental_downpass/uppass` from the join point gives the correct score in O(affected_subtree_depth × total_words) instead. **Plan:** make `apply_tbr_move` push touched nodes onto the prealloc_undo stack, return the join node, and replace line 1134's `full_rescore` with `fitch_incremental_downpass` from that node. Estimated savings: O(n_char) per accepted move × ~10–100 accepted moves per replicate. **Risk:** medium — `apply_tbr_move` is the hot correctness-critical path; need careful unit tests covering NA/non-NA, IW/EW, constrained/unconstrained, equal-accept paths. Validate by comparing scores against current unconditional rescore on a battery of datasets before committing. |






### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-28 focus 31 by G. ts_prune_reinsert.h/.cpp (583 lines): G-006 found + now fixed. Next: ts_search.cpp (NNI/SPR, 421 lines) and ts_nni_perturb.h/.cpp (unreviewed). [m:opus e:high] |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-05-12 round 7 by d6. T-260 hotspot audit: std::fill (9.1%) fixed T-261 ✓; StateSnapshot per-candidate save (14.6%) mitigated by opt #7 (once-per-pass) ✓; full_rescore line 1137 (~28%) → T-300 in progress (d7). No new tasks; re-profile with VTune after T-300 lands. |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-05-12 round 47 by d5. u.118 triaged → T-301 (progress ticker multi-thread). PR #210 (cpp-search→main): 2 CI failures are infra (ASAN vignettes: missing pkgdown; Windows: code-coverage only — R CMD check passes). PR #213 (T-150): CONFLICTING, no recent CI — needs human to resolve merge conflicts. PR #216 (T-204): agent-check 23649607006 PASSED; full R-CMD-check had failures Mar 2026 on ASAN/Windows/ubuntu-old — needs re-trig or human review. Active: d1(T-294), d2(T-298), d3(T-299), d4(S-RED parked). T-280–288 WORKTREE awaiting. S-PROF/S-PR OPEN. [m:haiku e:low] |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-05-12 round 48 by d5. PR #210 (cpp-search→main): MERGEABLE, fresh checks 2026-05-12 confirm 2 infra-only failures (Windows=code-coverage step, ASAN=vignettes infra) — all R-CMD-check PASS; ready for human to un-DRAFT and merge. PR #216 (T-204, feature/native-search→cpp-search): CONFLICTING — cpp-search has ~10 new commits since last merge; needs rebase then re-trig. PR #213 (T-150, feature/cid-consensus→cpp-search): CONFLICTING, no CI — needs rebase onto cpp-search. [m:sonnet e:medium] |
