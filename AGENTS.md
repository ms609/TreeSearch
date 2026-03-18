# TreeSearch Multi-Agent Development Notes

## Build isolation — per-agent library directories

Each agent **must** build and test to its own private library:

```bash
R CMD INSTALL --library=.agent-X .
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='ts-')"
```

**Never** install to the default library. On Windows, a loaded DLL locks
the file and blocks other agents.

**Never** use `devtools::load_all()` or `pkgbuild::compile_dll()` — these
target a shared temp location and will conflict.

## Build failure recovery

### Debug `.o` contamination

`roxygen2::roxygenise()` (default mode) calls `pkgbuild::compile_dll(debug=TRUE)`,
which leaves debug `.o` files in `src/`. Subsequent `R CMD INSTALL` reuses them,
producing a DLL that crashes at runtime (exit code 127/139).

**Fix:** `rm -f src/*.o src/*.dll` then rebuild.

**Prevention:** Never use bare `roxygen2::roxygenise()`. To regenerate docs:
```bash
Rscript -e ".libPaths(c('.agent-X', .libPaths())); roxygen2::roxygenise(load_code = roxygen2::load_installed)"
```

### DLL lock

If `R CMD INSTALL` fails with "Access is denied", another R process has the
DLL loaded. Kill it or wait, then retry.

### `TreeSearch-init.c` arg count mismatch

After `Rcpp::compileAttributes()`, **always** run `Rscript check_init.R` to
verify arg counts match between `RcppExports.cpp` and `TreeSearch-init.c`.

### Quick recovery

```bash
rm -f src/*.o src/*.dll
R CMD INSTALL --library=.agent-X .
Rscript check_init.R
```

## CPU limits — max 2 cores per agent

Use `nThreads = 2L` at most in tests/benchmarks. Never `nThreads = 0L`
(auto-detect). Use `-j2` at most for make.

## Shared files — coordination rules

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` are modified by every agent.
**Append only** — add new entries at the end. Do not reformat or reorder.

### `concavity` sentinel

Rcpp can't translate `R_PosInf`. All Rcpp-exported functions use
`concavity = -1.0` as the C++ default (sentinel for "equal weights / Inf").
Conversion `if (concavity < 0) concavity = HUGE_VAL` happens at three
gateway points in `ts_rcpp.cpp`: `make_dataset()`, `ts_resample_search()`,
`ts_successive_approx()`.

### `src/Makevars.win`

**Never leave a `src/Makevars.win` in place.** Debug/PGO/UBSan flags cause
crashes or miscompilation. Delete after any profiling session.

## Multi-agent workflow protocol

### Assignment

On `/assign X`:

1. Read `agent-X.md`. If a task is already in-progress, resume it.
2. Otherwise, check `issues.md` **before** `to-do.md`:
   a. If `issues.md` contains any unclaimed issues (blocks whose first line
      does **not** start with `CLAIMED`), **claim the bottom-most unclaimed
      issue** by prepending `CLAIMED (X):` to its first line.
   b. Triage the claimed issue: determine what needs doing, then add one or
      more discrete tasks to `to-do.md` (assign appropriate IDs and
      priorities — issues may be P0). Begin work on the first task.
   c. Once the `to-do.md` tasks are created, delete the entire issue block
      (including its `---` separator) from `issues.md`.
   d. **While `issues.md` still has unclaimed issues, triaging them takes
      priority over picking up existing `to-do.md` tasks** (an issue may
      contain a P0).
3. If `issues.md` is empty or all issues are already claimed, claim the next
   OPEN task from `to-do.md` as before.

Set `CONVERSATIONSUMMARY` to `Agent X: <task description>`.

> **Concurrency guard:** Only the bottom-most *unclaimed* issue may be
> claimed. Because agents always target the bottom and mark it `CLAIMED (X)`
> immediately, two agents will never parse the same issue. If an agent sees
> the bottom issue is already `CLAIMED`, it moves up to the next unclaimed
> one.

### During work

- Update `agent-X.md` after every significant step (crash-recovery record).
- All work uses `.agent-X/` as library directory.
- **All builds, tests, and benchmarks in bash subprocesses** — never in the
  RStudio R session.

### On task completion

1. Move task to Completed in `to-do.md`.
2. Set `agent-X.md` to IDLE.
3. Append a brief entry to this file documenting what changed.
4. Update `coordination.md` if strategic objectives are affected.
5. Take next task.

### Standing tasks

| ID | Type | Expertise file |
|----|------|---------------|
| S-RED | Red-team review | `.positai/expertise/red-team.md` |
| S-PROF | Performance profiling | `.positai/expertise/profiling.md` |
| S-COORD | Coordination review | `.positai/expertise/coordination.md` |

Priority: P3 when ≥6 OPEN tasks, P2 when 3–5, P1 when <3.

### Key files

| File | Purpose |
|------|---------|
| `issues.md` | Human-entered issues (agents triage → `to-do.md`) |
| `to-do.md` | Task queue |
| `coordination.md` | Strategic plan |
| `agent-X.md` | Agent progress log |
| `AGENTS.md` | Conventions + architecture reference |
| `.positai/expertise/*.md` | Standing task methodology |

## Test file conventions

All `tests/testthat/test-ts-*.R` files must use `TreeSearch:::` to call
internal C++ bridge functions. Define short local wrappers for readability.

Shared helpers are in `tests/testthat/helper-ts.R` (`make_ts_data()`,
`ts_score()`, `validate_result()`, `skip_extended()`).

**Never use `%in%` on Splits objects in test files** — S3 dispatch fails
in the cloned namespace created by `test_check()`. Use `as.logical()`
matrix comparison instead.

### Test tiering

Every new `test-ts-*.R` file must be assigned to one of three tiers.
See `tests/testing-strategy.md` for the full rationale.

| Tier | Guard | When it runs | Use for |
|------|-------|-------------|---------|
| 1 — CRAN | none | always (CRAN + CI + local) | Fast (< ~2 s) API and data-structure unit tests |
| 2 — CI | `skip_on_cran()` at **file level** (first executable line) | CI + local | C++ engine correctness, scoring, search algorithms |
| 3 — Extended | `skip_extended()` at **file level** | `TREESEARCH_EXTENDED_TESTS=true` only | Stress tests, benchmarks, timing measurements |

**Default for new `test-ts-*` files: Tier 2.** Add `skip_on_cran()` as the
very first executable line (before any helpers or `test_that()` calls):

```r
# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
```

Use Tier 3 only for tests that take > ~10 s or are sensitive to machine load.

## R source file ordering

`DESCRIPTION` has an explicit `Collate:` field. When adding a new `.R` file,
**update the Collate field** — otherwise R sources alphabetically, which can
break if one file's top-level code depends on a later file.

## Architecture reference

### R-level API

| Function | Engine | Purpose |
|----------|--------|---------|
| `MaximizeParsimony()` | C++ driven search | Primary search (EW, IW, profile, constraints) |
| `Morphy()` | R-loop + MorphyLib | Legacy search (custom stopping, per-iteration callbacks) |
| `MaximizeParsimony2()` | — | Deprecated alias for `MaximizeParsimony()` |
| `Resample()` | C++ | Jackknife/bootstrap resampling |
| `SuccessiveApproximations()` | C++ | Successive approximations weighting |
| `TreeLength()` | C++ `ts_fitch_score` | Score one or more trees |
| `FastCharacterLength()` | C++ `ts_char_steps` | Per-character step counts |
| `AdditionTree()` | C++ `ts_wagner_tree` | Wagner tree construction |
| `RandomTreeScore()` | C++ (phyDat) or MorphyLib (morphyPtr) | Score a random tree |
| `TaxonInfluence()` | C++ via `MaximizeParsimony()` | Per-taxon search |
| `SearchControl()` | — | Expert parameter constructor for `MaximizeParsimony()` |

`MaximizeParsimony()` has a backward-compatibility shim: passing old
Morphy-style parameters (`ratchIter`, `tbrIter`, etc.) triggers a deprecation
warning and delegates to `Morphy()`. Scheduled for removal in 2028.

### Driven search pipeline per replicate

1. Random Wagner tree → optional SPR → TBR to local optimum
2. XSS sectorial search (if tree large enough)
3. RSS random sectorial search
4. CSS constrained sectorial search
5. Ratchet perturbation to escape local optima
6. Drift search (accept suboptimal moves)
7. Final TBR polish
8. Add to pool
9. Fuse against pool (every `fuse_interval` replicates)

Post-search: TBR plateau enumeration from all pool seeds to find MPTs.

### Strategy presets (auto-selected by `NTip` and signal density)

| Preset | Condition | Key settings |
|--------|-----------|-------------|
| sprint | ≤30 tips | 3 ratchet, 0 drift, XSS only |
| default | 31–74 tips; or ≥75 tips with ≥5 chars/taxon | 5 ratchet, 2 drift, XSS+RSS |
| thorough | ≥75 tips with <5 chars/taxon | 20 ratchet (adaptive), 12 drift, XSS+RSS+CSS |

Signal-density gate: datasets with many characters per taxon converge
more easily, so stay on "default" even at larger sizes.

### C++ module map

| Module | Header/Source | Purpose |
|--------|--------------|---------|
| Fitch scoring | `ts_fitch.h/.cpp` | Downpass, uppass, incremental, indirect |
| NA scoring | `ts_fitch_na.h` | Three-pass inapplicable algorithm (Brazeau et al. 2019) |
| NA incremental | `ts_fitch_na_incr.h` | Incremental NA-aware scoring for TBR/drift |
| SIMD | `ts_simd.h` | SSE2/NEON portability layer for bit-parallel ops |
| Data | `ts_data.h/.cpp` | `DataSet`, `CharBlock`, `build_dataset`, simplification |
| Tree | `ts_tree.h/.cpp` | `TreeState`, topology manipulation, `PreallocUndo` |
| Constraint | `ts_constraint.h/.cpp` | Topological constraint enforcement |
| TBR | `ts_tbr.h/.cpp` | TBR search (with sector_mask for CSS) |
| SPR/NNI | `ts_search.h/.cpp` | SPR and NNI search (standalone, not in driven pipeline) |
| Ratchet | `ts_ratchet.h/.cpp` | Perturbation (zero/upweight/mixed, adaptive) |
| Drift | `ts_drift.h/.cpp` | Accept suboptimal moves within AFD/RFD limits |
| Wagner | `ts_wagner.h/.cpp` | Greedy addition tree (incremental scoring, NA-aware) |
| Sectorial | `ts_sector.h/.cpp` | RSS, XSS, CSS |
| Fuse | `ts_fuse.h/.cpp` | Tree fusing (in-place exchange) |
| Pool | `ts_pool.h/.cpp` | Dedup via split hashing, score-based eviction |
| Splits | `ts_splits.h/.cpp` | Bipartition computation and comparison |
| Driven | `ts_driven.h/.cpp` | Multi-replicate orchestrator |
| Resample | `ts_resample.h/.cpp` | Jackknife, bootstrap, successive approximations |
| Parallel | `ts_parallel.h/.cpp` | `std::thread` inter-replicate parallelism |
| RNG | `ts_rng.h/.cpp` | Thread-safe RNG (`thread_local` dispatch) |
| Simplify | `ts_simplify.h/.cpp` | Character compression and uninformativeness checks |
| Rcpp bridge | `ts_rcpp.cpp` | All Rcpp-exported functions |

### Scoring modes

`ScoringMode` enum in `ts_data.h`: `EW`, `IW`, `PROFILE`.
- **EW**: standard Fitch parsimony
- **IW**: implied weights via `e/(k+e)` where `e = steps - min_steps`
- **PROFILE**: lookup in `info_amounts` table (structurally identical to IW pipeline)

Profile mode sets `ds.concavity = 1.0` (finite sentinel) so existing
`isfinite()` checks activate the weighted pipeline without code duplication.

### Parallelism design

- `std::thread` (not OpenMP) to avoid R memory allocator conflicts
- Per-thread: `DataSet` copy, `ConstraintData` copy, `std::mt19937` RNG
- Shared: `ThreadSafePool` (mutex-guarded), atomic stop flag
- Main thread: pre-generates seeds from R's RNG, polls
  `R_CheckUserInterrupt()` and timeout every 200ms
- Worker threads make no R API calls — `ts_rng.h` provides `thread_local`
  dispatch (null → R API for serial; set → thread-local for parallel)

### NA scoring notes

- `.h` file changes (`ts_fitch_na.h`, `ts_fitch_na_incr.h`) may require
  `touch src/ts_fitch.cpp` before rebuild if the build system doesn't track
  header dependencies.
- Wagner NA incremental scoring has minor `subtree_actives` staleness
  (pessimistic but correct — `score_tree()` gives authoritative final result).
- TBR uses incremental NA heuristic for candidate screening + full three-pass
  verification on acceptance.

### Constraint enforcement

- `build_constraint()` reads R split matrix with **column-major** indexing:
  `split_matrix[s + n_splits * t]`.
- Wagner uses LCA-based constraint mapping (`wagner_map_constraint_nodes`)
  since splits aren't fully present during incremental construction.
- Wagner has a posthoc retry loop (up to 100 random addition orders) as a
  safety net for edge cases.

## Exported Rcpp functions

All registered in `ts_rcpp.cpp` and `TreeSearch-init.c`. Run
`Rscript check_init.R` to verify consistency.

| Function | Module | Purpose |
|----------|--------|---------|
| `ts_fitch_score` | ts_fitch | Score a tree |
| `ts_char_steps` | ts_rcpp | Per-pattern step counts (with simplification offsets) |
| `ts_na_debug_char` | ts_fitch_na | Per-node debug for a single pattern |
| `ts_na_char_steps` | ts_fitch_na | Per-pattern step counts (raw, no offsets) |
| `ts_debug_clip` | ts_fitch | Debug SPR clip/regraft |
| `ts_test_indirect` | ts_fitch | Debug indirect length |
| `ts_nni_search` | ts_search | NNI hill-climbing |
| `ts_spr_search` | ts_search | SPR hill-climbing |
| `ts_tbr_search` | ts_tbr | TBR with plateau exploration |
| `ts_ratchet_search` | ts_ratchet | Ratchet perturbation |
| `ts_drift_search` | ts_drift | Drift search |
| `ts_wagner_tree` | ts_wagner | Wagner tree (specified addition order) |
| `ts_random_wagner_tree` | ts_wagner | Wagner tree (random order) |
| `ts_compute_splits` | ts_splits | Bipartition splits from edge matrix |
| `ts_trees_equal` | ts_splits | Compare two trees |
| `ts_pool_test` | ts_pool | Pool deduplication test |
| `ts_tree_fuse` | ts_fuse | Fuse two trees |
| `ts_sector_diag` | ts_sector | Sectorial search diagnostics |
| `ts_rss_search` | ts_sector | Random Sectorial Search |
| `ts_xss_search` | ts_sector | Exclusive Sectorial Search |
| `ts_driven_search` | ts_driven | Full driven search |
| `ts_resample_search` | ts_resample | One jackknife/bootstrap replicate |
| `ts_successive_approx` | ts_resample | Successive approximations |
| `ts_parallel_resample` | ts_parallel | Batch resample with parallelism |
| `ts_bench_tbr_phases` | ts_rcpp | TBR phase timing diagnostic |

## MorphyLib deprecation status

Migration plan in `inst/deprecation/morphy-migration.md`.

**Already migrated to C++:** `MaximizeParsimony`, `AdditionTree`, `Resample`,
`SuccessiveApproximations`, `TreeLength`, `CharacterLength`,
`FastCharacterLength`, `RandomTreeScore`, `TaxonInfluence`.

**Still using MorphyLib:** Legacy search functions (`Ratchet`, `Jackknife`,
`MorphyBootstrap`, `CustomSearch`), R-level tree rearrangement functions.
These are candidates for deprecation rather than migration.

## Shiny app (`inst/Parsimony/`)

Decomposed from monolithic `app.R` into three-file Shiny convention:
- `global.R` — library calls, constants, helpers, colours, citations
- `ui.R` — `fluidPage(...)` definition
- `server.R` — server function shell, `reactiveValues()`, `source()` calls
- `server/*.R` — source files loaded with `source(local = TRUE)`
- `server/mod_*.R` — Shiny modules (`NS()`/`moduleServer()`)

**Completed modules:**
- `mod_references.R` — references panel (no state)
- `mod_treespace.R` — tree space visualization + plot settings.
  Returns 16 reactives (`distances`, `mapping`, `dims`, `saveDetails`,
  `TreespacePlot`, etc.) that `server.R` assigns to local scope for use
  by still-source'd files (clustering.R, consensus.R, downloads.R).
  Receives `clusterings`, `silThreshold`, `scores`, `concavity`,
  `distMeth`, `plotFormat` as reactive args + logging functions via
  `log_fns` list.

**Important:** Server source files are in `server/` NOT `R/`. Shiny 1.5+
auto-sources all `.R` files in an app's `R/` directory at startup (before
any session exists), which crashes on references to `output`/`input`/`session`.

Test suite: `NOT_CRAN=true` required for shinytest2 (4 test files, 33 assertions).
Module tests: `test-mod-references.R` (4 assertions), `test-mod-treespace.R` (4 assertions).

## Version and CRAN status

- **Version**: 2.0.0 (major bump for new `MaximizeParsimony()` API)
- **R CMD check**: 0 ERRORs, 0 WARNINGs, 1 NOTE (R 4.5.2 internal bug)
- **Test suite**: ~960 pass, 0 fail, ~18 skip

## Key design decisions (reference)

1. **PreallocUndo** (`ts_tree.h`): Pre-allocated flat buffers for TBR/drift
   undo stack. Uses `grow()` to dynamically expand when capacity exceeded
   (NA uppass saves both internal nodes and tips). Initial capacity `3 * n_node`.

2. **TBR symmetry breaking** (`ts_tbr.cpp`): FNV-1a hash deduplication of
   `virtual_prelim` vectors to skip redundant rerooting evaluations.

3. **Bounded indirect scoring**: All search modules use `_bounded` variants
   that bail out when accumulated score exceeds best candidate.

4. **Profile parsimony**: Reuses IW indirect pipeline unchanged; only delta
   precomputation differs. `ds.concavity = 1.0` sentinel activates weighted
   path. Max 2 informative states per character; inapplicable → ambiguous.

5. **MPT enumeration**: Post-search TBR plateau walk from all pool seeds.
   `tbr_search()` accepts optional `TreePool* collect_pool` parameter.

6. **All-ambiguous phyDat guard**: `TreeLength()` and `MaximizeParsimony()`
   check for `levels = NULL` / 0-column contrast matrix before calling C++.

## Benchmarks and profiling

Benchmark scripts in `inst/benchmarks/`. Key files:
- `bench_regression.R` — CI regression test (score quality + timing bounds)
- `bench_framework.R` — Dataset × strategy × replicate grid
- `strategies.md` — Strategy space documentation

Profiling baselines in `.positai/expertise/profiling.md`. Current phase
distribution (d2_r5 defaults, EW): TBR 11–33%, sectorial 13–26%,
ratchet 25–40%, drift 20–28%. Per-candidate indirect scoring is at
memory-throughput limit (~23 ns at 75 tips).
