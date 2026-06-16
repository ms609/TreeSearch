# Red-team findings — TreeSearch (OPEN)

Verified, non-trivial red-team findings that are **not yet resolved**. Trivial issues are
fixed inline during the round (and noted in `log.md`), never filed here. A finding is
filed here only **after** verification (see the `red-team-verifier` step in `/red-team`).

**Relationship to `to-do.md`:** the root `to-do.md` is the *dispatcher's* operational
queue (all task types). This file is the *red-team* view of the same open findings, so the
next rotation round can see at a glance what's already filed and avoid re-reporting it. When
a finding lands (PR merged), **remove the row** — resolved history lives in git, in
`to-do.md`/`completed-tasks.md`, and in the round entry in `log.md`. Do not keep a resolved
"trophy" table here.

Severity key: **P1** = wrong user-visible result / crash / desk-reject · **P2** = wrong on
edge input / frozen-shape inconsistency / search-quality · **P3** = robustness / polish.

| ID | Sev | Area | Title | File:line | Detail |
|----|-----|------|-------|-----------|--------|
| T-309 | P2 | 7 (Shiny) | EasyTrees: stale profile dataset scores wrong trees | `inst/Parsimony/server/mod_search.R:440` | On `profilePrepTask` completion the code runs `profileDataHash(r$dataHash)`, stamping the *current* dataset hash at completion time, not the hash of the dataset actually prepared. Load H2 while prep runs on H1 → `profileDataHash=hash(H2)` but `profileDataset=preparedFrom(H1)`; the `StartSearch()` guard (`:640`) then skips re-prep and `scores()` (`:475`, no hash check) scores H2 trees against the H1 profile dataset → wrong profile scores shown. `observeEvent(r$dataset)` (`:1128`) never clears `profileDataset/profileDataHash`. Fix: stamp the hash of the prepared dataset (snapshot at invoke), clear both on data change. Verified REAL (opus). Needs a mid-prep data swap to trigger. |
| T-310 | P2 | 7 (Shiny) | EasyTrees double-launch: no `searchInProgress` guard | `inst/Parsimony/server/mod_search.R:632` | `StartSearch()` lacks a re-entrancy guard; `shinyjs::disable("go")` is an async round-trip, so a fast double-click fires `observeEvent(input$go, StartSearch())` twice. Verified vs shiny 1.13.0 `ExtendedTask`: `invoke()` while running *queues* the 2nd call, which overwrites the single `cancelFile()`/`progressFile()` reactiveVals and `r$searchNotification` (leaks 1st notification — no `removeNotification` at `:719`) and re-enables Go mid-flight; on settle, task 1's trees may be silently dropped. Fix: `if (isTRUE(r$searchInProgress)) return(invisible())` at the top. One-line. Verified REAL (opus). |
| T-311 | P3 | 7 (Shiny) | EasyTrees: session disconnect never cancels the running worker | `inst/Parsimony/server.R:187` | `onStop` cleans file caches + cmd log but never writes the `cancelFile()` signal the `future::future()` worker polls (`mod_search.R:710`). A mid-search disconnect leaves the worker burning a core until it finishes its replicates or hits the timeout (~60 min for "thorough"). Fix: write the active cancel signal in `onStop` (or expose a module `cancel()`). Verified REAL (haiku). |
| T-312 | P3 | 7 (Shiny) | EasyTrees: search temp files (`ts_*`) leak on session end | `inst/Parsimony/server.R:192-194` | `onStop`'s `unlink(pattern="^(data\|tree\|excel)File-")` doesn't match the temp files `mod_search.R` creates: `ts_cancel_*`, `ts_progress_*`, `ts_profile_prog_*`, `ts_profile_cancel_*`. The worker `on.exit` clears some on the normal path; on error/interrupt/disconnect they accumulate in `tempdir()` (the documented "Issue 6" tempdir growth in `../expertise/shiny-app.md`). Fix: add a `ts_(cancel\|progress\|profile_prog\|profile_cancel)_` unlink to `onStop`. Verified REAL (haiku). |
| T-313 | P3 | 7 (Shiny) | EasyTrees: topology dedup includes branch lengths → inflated pool | `inst/Parsimony/server/mod_search.R:1063-1066` | The "topology string" dedup uses `write.tree(ape::ladderize(t))`, but `write.tree()` serialises branch lengths when present. After `combined <- c(r$allTrees, newTrees)` mixes user-loaded trees (which may carry BLs) with parsimony trees (no BLs), topologically identical trees with different BLs aren't deduplicated → inflated pool and displayed tree count. Fix: strip branch lengths before serialising (drop `$edge.length`, or use a topology-only key). Verified REAL (haiku). |
| T-322 | P3 | 8 (Tests) | Wagner NA+IW regression test is tautological (omits `min_steps`) | `tests/testthat/test-ts-wagner.R:223-242` | The test "Wagner on NA + IW matches fitch_score" calls `ts_random_wagner_tree(...)` and `ts_fitch_score(...)` both with `concavity = k` but **omits `min_steps`** (defaults to `integer(0)`), so the implied-weight homoplasy `h = steps − min_steps` is computed as `h = steps − 0` on *both* sides. The cross-check (Wagner incremental score == independent Fitch rescore of the same tree) therefore passes, but validates a non-production formula: the real NA+IW path (`R/MaximizeParsimony.R:834`) always passes `min_steps = as.integer(MinimumLength(ds, compress = TRUE))`, and Vinther2008 has inapplicable characters so `MinimumLength` is non-zero — a regression in NA+IW `min_steps` handling would slip through. Fix: pass `min_steps = as.integer(MinimumLength(pd, compress = TRUE))` to *both* calls (the fn accepts it, RcppExports.R:147) and re-run; the cross-check stays valid (same `min_steps` both sides) while now exercising production scoring. Filed not fixed inline: changes test numerics, needs a test-run to confirm (may itself surface a latent wagner NA+IW bug). Verified REAL (sonnet). |

<!--
Filing template (one row per verified finding):
| T-NNN | P1/P2/P3 | <area #> | **Title.** | `path:line` | Detail + fix + verifier verdict. |
-->
