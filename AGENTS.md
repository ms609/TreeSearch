# TreeSearch Multi-Agent Development Notes

## Current phase: bug-fixing / pre-release (as of 2026-03-20)

The project is in a **bug-fixing and stabilisation phase** with the goal of
shipping the package. Agents should:

- Monitor `to-do.md` as usual for task selection.
- **Prioritise bug fixes, test failures, documentation issues, and R CMD check
  problems** over new functionality.
- **Do not implement new features on `cpp-search` or `main`.**  
  Feature work is allowed only on dedicated `feature/<name>` branches, and
  only when the task is explicitly labelled as a feature and has been approved
  for active development.
- When in doubt, prefer a conservative fix (minimal diff, no API changes) over
  an ambitious refactor.

This phase ends when a clean `R CMD check` (0 errors, 0 warnings) is confirmed
and the maintainer signals readiness to tag a release.

---

## Build isolation — tarball workflow (mandatory)

Multiple agents share the same `src/` directory. In-place `R CMD INSTALL .`
compiles `.o` files and links the DLL directly in `src/`, causing races:
concurrent builds corrupt each other's object files, and a running R process
locks the DLL on Windows, blocking everyone else.

**Always build via tarball** so compilation happens in an isolated temp
directory:

```bash
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-X "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
```

Key points:
- `rm -f src/*.o src/*.dll` **must** precede every build — stale artifacts slow traversal and corrupt DLLs.
- Build into an agent-specific `$TMPBUILD` outside the source tree — avoids tarball collision when multiple agents build concurrently (R CMD build has no `--output=` flag; the tarball always lands in the working directory).
- `--no-resave-data` skips unnecessary `.rda` re-saving (not needed for dev installs).

To run tests:
```bash
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='ts-')"
```

**Never** use `R CMD INSTALL --library=.agent-X .` (in-place build).

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
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-X "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
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

### `src/TreeSearch-win.def`

**Keep this file.** It explicitly exports `R_init_TreeSearch` for Windows
DLL builds. Without it, the default `nm | sed` pipeline generates a
`tmp.def` that truncates long C++ mangled symbols, causing linker failures
or corrupt DLLs (especially under `pkgbuild::compile_dll(debug=TRUE)`).

## Branch structure

```
main              ← stable, taggable; receives only reviewed bug fixes
  └─ cpp-search   ← integration branch; all feature work merges here
       ├─ feature/cid-consensus
       ├─ feature/hsj-polish
       └─ feature/<name>   (one per major feature)
```

### Rules

- **`main`**: bug fixes and release tags only. No experiments.
- **`cpp-search`**: integration target. Bug-fix agents (S-RED, S-PROF, S-COORD,
  and ad-hoc fixes) work directly here. Feature branches merge here when ready.
- **`feature/*`**: branch from `cpp-search`; contain **code changes only**.
  Each feature branch is owned by a single agent at a time.

### Coordination files live on `cpp-search` only

`to-do.md`, `issues.md`, `agent-X.md`, `completed-tasks.md`, `coordination.md`,
and `AGENTS.md` are **never committed on feature branches**. When an agent
working on a feature branch needs to log progress or claim a task, they commit
those changes directly to `cpp-search` (coordination-only commit), keeping the
feature branch clean.

To read coordination files while on a feature branch without switching:
```bash
git show cpp-search:to-do.md
git show cpp-search:agent-X.md
```

To update a coordination file from a feature branch:
```bash
git checkout cpp-search -- agent-X.md   # pull latest into working tree
# edit, then:
git stash                               # stash any code changes first
git add agent-X.md && git commit -m "chore: agent X progress note"
git push origin cpp-search
git stash pop                           # restore code work
```

### Shared files at merge time

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` use the existing append-only
convention — merge conflicts resolve cleanly by keeping both appended blocks.
`DESCRIPTION` (Collate field) and `NAMESPACE` require a manual merge pass;
this is expected and should be done carefully at feature-merge time.

### Feature branch lifecycle

1. `git checkout cpp-search && git checkout -b feature/<name>`
2. Claim task on `cpp-search`'s `to-do.md` (coordination commit).
3. Do all code work on `feature/<name>`.
4. When complete: `git checkout cpp-search && git merge feature/<name>`.
5. Resolve any Collate/NAMESPACE conflicts, rebuild, run tests.
6. Log completion in `completed-tasks.md` on `cpp-search`.
7. Delete feature branch.

---

## Multi-agent workflow protocol

### Worktree tasks

Tasks with status `WORKTREE (name)` are actively developed in a dedicated git
worktree (e.g. `C:/Users/pjjg18/GitHub/TS-CID-cons`). **Do not claim or
modify these tasks.** They are reserved for the human developer working in
that worktree. To mark a task as in-flight on a worktree, set its status to
`WORKTREE (name)` where *name* matches the worktree directory basename.

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

1. **Delete** the task row from `to-do.md`. If the task was the last open
   row in a section/group, delete the section header too.
2. **Append** a summary row to `completed-tasks.md` under the current date
   heading (create a new `## YYYY-MM-DD` heading if needed):
   `| T-nnn | Short description | X | Brief notes |`
3. Set `agent-X.md` to IDLE.
4. Append a brief entry to this file documenting what changed.
5. Update `coordination.md` if strategic objectives are affected.
6. Take next task.

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
| `to-do.md` | Task queue (active/open tasks only) |
| `completed-tasks.md` | Archive of completed tasks |
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

## Documentation checks (mandatory)

After any change to a function signature or roxygen block, run:

```r
devtools::check_man()
```

After writing or updating documentation prose, also run:

```r
spelling::spell_check_package()
```

Both should be clean before committing. These are fast and catch issues
(`check_man` catches Rd parse errors, cross-ref failures, `\usage` mismatches;
`spell_check_package` catches typos in `@description`/`@details`/`@param` text).

References are added using Rdpack's \insertCite{}, with
\insertAllCited{} in the references section.

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
| `ParsSim()` | Pure R | Simulate datasets under parsimony (EW/IW/profile) |

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
| sprint | ≤30 tips | 3 ratchet, 0 drift, XSS only, consensus-stop 3 |
| default | 31–64 tips; or ≥65 tips with <100 char patterns | 5 ratchet, 2 drift, XSS+RSS, consensus-stop 3, Wagner×3, SPR-first, adaptive level |
| thorough | ≥65 tips with ≥100 char patterns | 20 ratchet (adaptive), 12 drift, XSS+RSS+CSS, consensus-stop 3, Wagner×3, SPR-first |

All presets set `consensusStableReps = 3`: search stops early if the strict
consensus of best-score pool trees is unchanged for 3 consecutive replicates.

`default` and `thorough` presets set `sprFirst = TRUE` (SPR warmup before
TBR) and `wagnerStarts = 3` (pick best of 3 random Wagner trees). `default`
also enables `adaptiveLevel = TRUE` (scale ratchet/drift by hit rate);
`thorough` omits it because high base cycle counts (20 ratchet, 12 drift)
already cover hard landscapes, and 1.5× scaling on top causes excessive
per-replicate time without commensurate score improvement.

Signal-density gate: datasets with few character patterns (<100) have flat
parsimony landscapes where intensive search adds no benefit.

### Adaptive sectorial search

XSS and CSS use **adaptive early-exit**: after each round of sector searches
+ global TBR polish, if the overall best score did not improve, remaining
rounds are skipped. This avoids wasting ~7% of replicate time on datasets
where sectorial search is unproductive (e.g. Dikow2009). On productive
datasets (e.g. Zhu2013), the early exit never fires.

### Conflict-guided RSS

RSS uses **conflict-guided sector selection**: before each replicate's RSS
phase, `driven_search()` computes a `SplitFrequencyTable` from the pool's
best-score trees. Within `rss_search()`, each internal node's "conflict
score" is `1 − (fraction of pool trees containing that split)`.
Max-descendant conflict is propagated upward, and eligible sector roots
are sampled via `std::discrete_distribution` with weight `1 + 3 × conflict`.
Falls back to uniform selection when the pool has <2 best-score trees or
when conflict variation is negligible.

### Consensus-stability stopping

After each replicate, if `consensus_stable_reps > 0` (default 3 in all
presets), the pool's strict consensus hash is compared to the previous
replicate's. If unchanged for `consensus_stable_reps` consecutive
replicates, the search terminates early. `compute_consensus_hash()` uses
XOR of per-split FNV-1a hashes for O(pool × splits) cost.

### Adaptive search level

When `adaptive_level = true`, ratchet and drift cycle counts are scaled
each replicate based on the cumulative hit rate:
- hit_rate > 0.7 → 0.5× (easy landscape)
- hit_rate > 0.4 → 0.75×
- hit_rate < 0.15 → 1.5× (hard landscape)
- else → 1.0×

### TBR zero-length clip skipping + regraft merging (collapsed flags)

`compute_collapsed_flags()` (`ts_collapsed.h/.cpp`) identifies edges where
clipping provably cannot improve score. Checks 5 conditions: (1) zero
standard-block cost at parent, (2) zero NA-block cost at parent, (3) prelim
preservation (`prelim[sibling] == prelim[parent]`), (4) down2 preservation
(NA), (5) subtree_actives preservation (NA). Works for EW, IW, Profile,
and NA-aware scoring. Integrated into TBR, SPR, and drift search.
Disabled during MPT enumeration (equal-score topologies may exist).
Recomputed after every accepted move.

**Regraft merging** (Goloboff 1996): within a collapsed region (connected
set of nodes linked by zero-length edges), all regraft positions yield the
same full score. Only boundary edges (entering the region) are evaluated;
interior collapsed edges are skipped via `if (collapsed[below]) continue`.
TBR, SPR, and drift all use this. The `CollapsedRegions` struct exists in
the header but callers use `compute_collapsed_flags()` directly (the
`region_id` field is unused — only the boolean flag array matters).

**Collapsed-topology pool dedup**: `compute_collapsed_splits()` in
`ts_splits.cpp` produces the split set excluding collapsed edges. Two
binary trees differing only in zero-length resolutions produce the same
collapsed split set → treated as duplicates by `TreePool::add_collapsed()`.
Both serial (`driven_search`) and parallel (`ThreadSafePool`) paths use
collapsed dedup.

**Benchmark results** (2026-03-22, 4 standard datasets, 3 seeds each):
Skip rate = 0% on all datasets (Vinther2008 23t, Agnarsson2004 62t,
Zhu2013 75t, Dikow2009 88t). Near-optimal trees in these morphological
datasets have negligible zero-length edges. Overhead from flag computation
is negligible. Score equivalence confirmed (enabled vs disabled produce
identical best scores). Benefit expected on sparse/synthetic data.

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
| Sectorial | `ts_sector.h/.cpp` | RSS (conflict-guided), XSS, CSS; from-above HTU |
| Fuse | `ts_fuse.h/.cpp` | Tree fusing (in-place exchange) |
| Pool | `ts_pool.h/.cpp` | Dedup, eviction, consensus hash, split frequency table |
| Splits | `ts_splits.h/.cpp` | Bipartition computation, comparison, `hash_single_split()` |
| Driven | `ts_driven.h/.cpp` | Multi-replicate orchestrator |
| Resample | `ts_resample.h/.cpp` | Jackknife, bootstrap, successive approximations |
| Parallel | `ts_parallel.h/.cpp` | `std::thread` inter-replicate parallelism |
| RNG | `ts_rng.h/.cpp` | Thread-safe RNG (`thread_local` dispatch) |
| Simplify | `ts_simplify.h/.cpp` | Character compression and uninformativeness checks |
| Collapsed | `ts_collapsed.h/.cpp` | Zero-length edge detection for clip skipping |
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

### Scoring notes

- `.h` file changes (`ts_fitch_na.h`, `ts_fitch_na_incr.h`) may require
  `touch src/ts_fitch.cpp` before rebuild if the build system doesn't track
  header dependencies.
- Incremental scoring is a **screening heuristic** for candidate selection;
  `full_rescore()` / `score_tree()` is always authoritative.
- See `.positai/expertise/fitch-scoring.md` for detailed invariants:
  uppass correctness proof, NA staleness analysis, `upweight_mask` audit.

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

Fully modularized from monolithic `app.R` into Shiny modules:
- `global.R` — library calls, constants, helpers, colours, citations, module UI instantiation
- `ui.R` — `fluidPage(...)` definition using module UI elements
- `server.R` — `AppState()` + module wiring + `ShowConfigs` observer + `onStop()`
- `server/app_state.R` — `AppState()` typed `reactiveValues()` constructor
- `server/logging.R` — session logging infrastructure
- `server/mod_*.R` — 7 Shiny modules (`NS()`/`moduleServer()`)

**All server logic now lives in modules.** The old `events.R` has been
dissolved; its `ShowConfigs` function and `plotFormat` observer are inlined
in `server.R` (they operate on top-level DOM elements).

**Modules:**
- `mod_references.R` — references panel (no state)
- `mod_downloads.R` — all 8 download handlers
- `mod_data.R` — data loading + tree management (9 returned reactives).
  Uses `cb_ref` forward-reference env for circular deps with consensus module.
- `mod_clustering.R` — clustering analysis + tree distances (5 returned reactives)
- `mod_search.R` — search engine, scoring, weighting.
  Owns ExtendedTask, search config modal, result accumulation.
- `mod_treespace.R` — tree space visualization + plot settings (14 returned reactives)
- `mod_consensus.R` — consensus plotting, character mapping, stability/rogue analysis,
  concordance, cluster consensus, main plot dispatch, plot logging (1327 lines).
  Returns `MainPlot`, `RCode`, `UpdateKeepNTipsRange`,
  `UpdateDroppedTaxaDisplay`, `UpdateOutgroupInput`.

**Important:** Server source files are in `server/` NOT `R/`. Shiny 1.5+
auto-sources all `.R` files in an app's `R/` directory at startup (before
any session exists), which crashes on references to `output`/`input`/`session`.

Test suite: `NOT_CRAN=true` required for shinytest2 integration tests.
Run from `inst/Parsimony/`:
```bash
NOT_CRAN=true Rscript -e "testthat::test_dir('tests/testthat')"
```
`setup.R` loads `library(shinytest2)` for `AppDriver` availability.

**Important:** Integration tests trigger `pkgbuild::compile_dll(debug=TRUE)`
via `load_all()`. `src/TreeSearch-win.def` prevents linker failures from
corrupted auto-generated `tmp.def` on Windows.

Module tests: `test-mod-references.R` (4), `test-mod-data.R` (9),
`test-mod-clustering.R` (12), `test-mod-treespace.R` (5),
`test-mod-downloads.R` (11), `test-mod-search.R` (28),
`test-mod-consensus.R` (9).
Integration tests: `test-app-smoke.R` (3), `test-Distribution.R` (13),
`test-SearchLog.R` (4), `test-ViewChars.R` (12). Total: 110 assertions.

## Version and CRAN status

- **Version**: 2.0.0 (major bump for new `MaximizeParsimony()` API)
- **R CMD check**: 0 ERRORs, 0 WARNINGs, 1 NOTE (R 4.5.2 internal bug)
- **Test suite**: ~9200 R-level + 1859 ts-* + 128 ParsSim + 37 MaddisonSlatkin + 49 recode-hierarchy pass

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

7. **From-above HTU for sectorial search** (`ts_sector.cpp`):
   `compute_from_above_for_sector()` computes `from_above[sector_root]` —
   the Fitch state-set the rest of the tree sends *down* to the sector
   boundary, excluding the sector's own contribution. Used instead of
   `final_[parent]` in `build_reduced_dataset()`. O(depth × total_words).

8. **Split frequency table** (`ts_pool.h/.cpp`): `SplitFrequencyTable` maps
   per-split FNV-1a hash → occurrence count across best-score pool trees.
   Used by conflict-guided RSS to weight sector selection. The same FNV-1a
   hash (`hash_single_split()` in `ts_splits.h`) is used by consensus
   hashing and split frequency counting — must stay consistent.

9. **Consensus-stability hash** (`ts_pool.cpp`): XOR of FNV-1a hashes of
   splits present in ALL best-score trees. Updated after each replicate.
   Hash collision false-matches are conservative (over-count stability).

10. **Diversity-aware pool eviction** (`ts_pool.cpp`): When the pool is full
    and a new tree ties the worst score, the entry most similar to the new
    tree (most shared splits, counted via per-split FNV-1a hash set
    membership) is evicted. This maintains topological diversity in the pool,
    improving fusing effectiveness. Falls back to arbitrary worst entry when
    the new tree is strictly better.

11. **Cross-replicate consensus constraint tightening** (`ts_driven.cpp`):
    When `consensus_constrain = true` and no user constraint is supplied,
    after ≥5 replicates, unanimous pool splits are extracted and enforced
    as topological constraints via `build_constraint_from_bitsets()`. The
    TBR/SPR search then avoids breaking established consensus clades.
    Constraints are cleared and rebuilt whenever the best score changes.
    Sector/fuse operations do not enforce auto-constraints (no posthoc
    DataSet is built).

## Alternative inapplicable-handling algorithms (in progress)

Plan: `.positai/plans/2026-03-19-0643-alternative-inapplicable-handling-algorithms.md`

Adding HSJ (Hopkins & St. John 2021) and step-matrix/x-transformation
(Goloboff et al. 2021) scoring as alternatives to the existing Brazeau
et al. (2019) three-pass algorithm. Both require an explicit character
hierarchy specification.

### New files

| File | Purpose | Status |
|------|---------|--------|
| `R/CharacterHierarchy.R` | `CharacterHierarchy` S3 class, `validate_hierarchy()`, `hierarchy_from_names()`, `hierarchy_chars()`, `hierarchy_controlling()`, `non_hierarchy_weights()` | Complete, 34 tests passing |
| `tests/testthat/test-CharacterHierarchy.R` | Unit tests for hierarchy specification + weight partitioning | Complete |
| `src/ts_hsj.h` | `HierarchyBlock` struct (with `absent_state`), `hsj_score()` declaration, `partition_weights()` | Complete |
| `src/ts_hsj.cpp` | `partition_weights()`, `fitch_label_char()` (with uppass), `score_hierarchy_block()`, `hsj_score()` | Complete (full-rescore only; not wired to search pipeline) |
| `src/ts_sankoff.h` | `SankoffChar`, `SankoffData` structs, `sankoff_score()`, `sankoff_score_char()`, `sankoff_uppass()` | Complete |
| `src/ts_sankoff.cpp` | Sankoff downpass, uppass, root forcing | Complete |
| `R/recode_hierarchy.R` | `recode_hierarchy()`: x-transformation recoding (Goloboff et al. 2021) | Complete, 49 tests |
| `tests/testthat/test-recode-hierarchy.R` | Unit tests for recode_hierarchy() | Complete |
| `inst/REFERENCES.bib` | Added `Goloboff2021b` entry | Complete |

### Modified files

| File | Change |
|------|--------|
| `DESCRIPTION` | Added `CharacterHierarchy.R` to Collate field |
| `R/MaximizeParsimony.R` | Added `hierarchy`, `inapplicable`, `hsj_alpha` params with validation; non-brazeau methods currently `stop()` with "not yet implemented" |
| `src/ts_data.h` | Added `inapp_state` field to `DataSet` (for HSJ) |
| `src/ts_data.cpp` | Populate `inapp_state` in `build_dataset()` |

### Design decisions

- `hierarchy` is a **separate argument** to `MaximizeParsimony()` (not a phyDat attribute)
- `inapplicable` and `hsj_alpha` are **top-level args** alongside `concavity`
- Default `hsj_alpha = 1.0`
- IW + hierarchy and Profile + hierarchy: **deferred**
- Constraint interaction: **ignored** for now
- Resampling: **hierarchical** — resample top-level chars; when a controlling primary is sampled, also resample within its block; recurse for nested hierarchies

### Resampling with hierarchy (T-124)

`Resample()` now accepts `hierarchy`, `inapplicable`, and `hsj_alpha`
parameters. When `inapplicable != "brazeau"`, resampling is hierarchy-aware:

- **Resampling units**: each non-hierarchy character = 1 unit; each
  top-level hierarchy block (primary + all dependents) = 1 atomic unit.
- **Jackknife**: retain `proportion` of units without replacement.
- **Bootstrap**: sample `n_units` units with replacement (blocks can be
  duplicated).
- Per replicate: `.HierarchicalResampleWeights()` computes pattern weights
  for non-hierarchy chars and per-block sample counts. `.ResampleHierarchy()`
  calls `ts_driven_search` per replicate with filtered HSJ blocks or xform
  chars.
- **No C++ changes**: reuses existing `ts_driven_search` HSJ/xform infrastructure.
- **Parallelism**: serial R loop over replicates (C++ inter-search parallelism
  via `nThreads` still available within each replicate). Adding C++-level
  inter-replicate parallelism is a future optimization.

### Remaining work (Phase 1c–f)

1. ~~Pass `absent_token`, `n_tokens` from R to C++~~ **Done** (T-115): `absent_state` in HierarchyBlock, `inapp_state` in DataSet.
2. ~~Partition original characters into hierarchy vs non-hierarchy sets~~ **Done** (T-115): `partition_weights()` (C++) and `non_hierarchy_weights()` (R).
3. ~~Implement `hsj_score()` core algorithm in `ts_hsj.cpp`~~ **Done** (T-116): `fitch_label_char()` + `score_hierarchy_block()` + `hsj_score()`.
4. ~~Add Rcpp bridge function for HSJ scoring in `ts_rcpp.cpp`~~ **Done** (T-116): `ts_hsj_score()` registered in init.c.
5. ~~R-side marshalling~~ **Done** (T-116): `build_tip_labels()`, `hierarchy_to_blocks()`.
6. ~~Remove placeholder `stop()`~~ **Done** (T-117): HSJ wired into `score_tree()` dispatch, `ts_driven_search()` bridge, `MaximizeParsimony()`. End-to-end test against paper examples (T-118)

### Key algorithm notes (HSJ)

- Paper's Algorithm 1 initializes `a(l) = p(l) = 0` for all leaves. This is
  incorrect for enforcing observed leaf states. Correct initialization:
  leaf with primary absent → `a(l) = 0, p(l) = INF`; primary present →
  `a(l) = INF, p(l) = 0`. Verified against hand-computed example.
- `score_hierarchy_block()` operates per hierarchy block. Non-hierarchy
  characters use standard Fitch. Total = Fitch(non-hierarchy) + Σ HSJ(blocks).
- Secondary character labels at internal nodes from Fitch first-pass
  (inapplicable treated as a separate state).
- HSJ is full-rescore only (no incremental variant). Performance mitigation:
  candidate screening via Fitch, full HSJ only for promising candidates.

### Phase 2 (step-matrix) — Complete (end-to-end functional)

Sankoff engine (`ts_sankoff.h/.cpp`) implements downpass, uppass, root forcing.
R-level `recode_hierarchy()` combines primary + secondaries into composite
step-matrix character with asymmetric costs (gain:loss = n+1:1). Multistate
secondaries supported (state count = ∏k_i + 1). Nested hierarchies deferred.
Integration complete: `ScoringMode::XFORM` in `score_tree()` dispatches
Fitch(non-hierarchy) + Sankoff(recoded). `MaximizeParsimony()` accepts
`inapplicable = "xform"`. End-to-end search verified.

## Search optimization roadmap

Plan: `.positai/plans/2026-03-21-search-optimizations.md`

Ranked by priority:
1. ~~Consensus-guided sector targeting~~ — **Done**: RSS weighted by
   pool split conflict scores
2. ~~Diverse pool maintenance~~ — **Done**: evict most-similar entry on ties
3. ~~Cross-replicate constraint tightening~~ — **Done**: opt-in via
   `consensusConstrain = TRUE`
4. ~~Collapsed-tree clip skipping~~ — **Done**: zero-length edges skipped
   in TBR, SPR, and drift. Benchmark shows 0% skip rate on standard
   morphological datasets (Vinther2008, Agnarsson2004, Zhu2013, Dikow2009)
   because near-optimal trees have few zero-length edges. Negligible
   overhead. Benefit expected primarily on sparse/synthetic data.
   Full polytomy search remains post-2.0.0.
5. ~~Collapsed-region regraft merging + pool dedup~~ — **Done**: within
   collapsed regions (connected zero-length edges), only the boundary
   regraft position is evaluated (Goloboff 1996). Collapsed-topology
   pool dedup treats trees differing only in zero-length resolutions as
   duplicates. Parallel path also uses collapsed dedup. Diversity-aware
   pool eviction selects most-similar entry on ties.
6. ~~Strategy preset tuning~~ — **Done**: `default` preset now uses
   `wagnerStarts=3`, `sprFirst=TRUE`, `adaptiveLevel=TRUE`; `thorough`
   preset uses `sprFirst=TRUE`.

## Benchmarks and profiling

Benchmark scripts in `inst/benchmarks/`. Key files:
- `bench_regression.R` — CI regression test (score quality + timing bounds)
- `bench_framework.R` — Dataset × strategy × replicate grid
- `strategies.md` — Strategy space documentation

Profiling baselines in `.positai/expertise/profiling.md`. Phase distribution
(default strategy, EW, 2026-03-21):

| Phase | Zhu2013 (75t) | Dikow2009 (64t) | Vinther2008 (23t) |
|-------|:---:|:---:|:---:|
| Ratchet | 27% | 36% | 42% |
| Drift | 21% | 21% | 22% |
| TBR | 38% | 29% | 17% |
| XSS | 8% | 9% | 12% |
| RSS | 4% | 3% | — |
| Final TBR | 2% | 2% | 6% |

Per-candidate indirect scoring is at memory-throughput limit (~23 ns at
75 tips).
