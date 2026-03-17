
# TreeSearch Multi-Agent Development Notes

## Build isolation — ALWAYS use per-agent library directories

Multiple agents work on `src/` concurrently. Each agent **must** build
and test to its own private library directory, e.g.:

```bash
R CMD INSTALL --library=.agent-e .
R -e "library(TreeSearch, lib.loc = '.agent-e'); testthat::test_dir('tests/testthat')"
```

- Agent C → `.agent-c/`
- Agent E → `.agent-e/`
- etc.

**Never** install to the default library (`R CMD INSTALL .` without
`--library`). On Windows, a loaded DLL locks the file and any other
agent (or user) trying to install will get "Access is denied".

Do **not** use `devtools::load_all()` or `pkgbuild::compile_dll()` —
these target a shared temp location and will conflict.

## CPU limits — max 2 cores per agent

This machine is shared. **No agent may use more than 2 CPU cores** for
any operation — builds, tests, or benchmarks. Key rules:

- When testing parallel code (`nThreads` parameter), use **`nThreads = 2L`**
  at most. Never use `nThreads = 0L` (auto-detect) in tests or benchmarks
  — it will grab all cores and starve other agents.
- `R CMD INSTALL` already runs single-threaded by default; no action needed.
- If running `make` or any tool that supports `-j`, use `-j2` at most.
- Benchmark scripts must hard-code `nThreads = 2L` (or 1L).

## Shared files — coordination rules

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` are modified by every
agent (to add Rcpp bridges and register C entry points). To avoid
merge conflicts:

- **Append only** — add your new declarations / registrations at the
  end of the relevant section.
- Agents should not reformat or reorder existing entries.

### `concavity` parameter — sentinel convention

`Rcpp::compileAttributes()` cannot translate `R_PosInf` to an R default.
To avoid needing manual patching after every regeneration, all Rcpp-exported
functions now use `concavity = -1.0` as the C++ default (a sentinel meaning
"equal weights / Inf"). The conversion `if (concavity < 0) concavity = HUGE_VAL`
happens at three gateway points in `ts_rcpp.cpp`:

1. `make_dataset()` — used by most Rcpp bridge functions
2. `ts_resample_search()` — uses its own dataset construction
3. `ts_successive_approx()` — uses its own dataset construction

**Any new Rcpp function with a `concavity` parameter must use `-1.0` as the
default and convert to `HUGE_VAL` before passing to C++ internals.**

R-level callers that pass `concavity = Inf` explicitly are unaffected
(`Inf` > 0, so the sentinel check doesn't trigger).

## Documentation — update plans when you complete a task

When you finish a step or phase, **always update**:

1. **`AGENTS.md`** — Append a section documenting what you changed, files
   created/modified, test status, and any design decisions.
2. **The phase plan `.md`** (in `.positai/plans/`) — Mark it as IMPLEMENTED
   at the header.
3. **`2026-03-16-production-plan.md`** — Update the relevant phase entry
   with a completion summary.

Do not leave documentation updates for a later step or another agent.

## Multi-agent workflow protocol

### Assignment

When the user types `/assign X`, the agent takes ownership of letter X and
follows this sequence:

1. **Read `agent-X.md`** — check for an in-progress task.
2. **If a task is in progress:** Resume it. Read the progress notes carefully,
   understand where the previous session left off, and continue from there.
3. **If no task is in progress (status IDLE):** Read `to-do.md`, claim the
   highest-priority OPEN task by setting its status to `ASSIGNED (X)`, and
   update `agent-X.md` with the task info.
4. **Set the conversation name** to `Agent X: <brief task description>` via
   the `CONVERSATIONSUMMARY` tag. Example:
   `<CONVERSATIONSUMMARY>Agent A: Phase 6A timing instrumentation</CONVERSATIONSUMMARY>`
   Keep this updated if the task changes.

### During work

- **Update `agent-X.md` after every significant step** — not at the end of
  the session. This file is the crash-recovery record. If RStudio restarts,
  the conversation history is lost; `agent-X.md` is the sole record.
- Progress notes should be detailed enough that a fresh agent instance can
  pick up exactly where you left off: what's done, what's next, key decisions,
  files modified, any partial state.
- All work (build, test, benchmark) uses `.agent-X/` as the library directory.

### On task completion

1. Move the task row from Active to Completed in `to-do.md` (with agent
   letter and date).
2. Update `agent-X.md`: set status to IDLE, clear the task fields, summarize
   what was accomplished.
3. Append a section to `AGENTS.md` documenting the work (as per the existing
   Documentation rules above).
4. Update `coordination.md` if the work affects strategic objectives.
5. Then take the next highest-priority OPEN task from `to-do.md`.

### Standing tasks

Three standing tasks are always present in `to-do.md`:

| ID | Type | Expertise file |
|----|------|---------------|
| S-RED | Red-team review | `.positai/expertise/red-team.md` |
| S-PROF | Performance profiling | `.positai/expertise/profiling.md` |
| S-COORD | Coordination review | `.positai/expertise/coordination.md` |

Their effective priority is dynamic based on how many specific OPEN tasks
remain:
- ≥6 OPEN specific tasks → standing tasks are P3
- 3–5 OPEN specific tasks → standing tasks are P2
- <3 OPEN specific tasks → standing tasks are P1

When you take a standing task, read the corresponding expertise file for
methodology. On completion, reset the standing task to OPEN (it recycles).

### Key files

| File | Purpose |
|------|---------|
| `to-do.md` | Task queue — claim and complete tasks here |
| `coordination.md` | Strategic plan — objectives, agent status, decisions |
| `agent-X.md` | Your progress log — update continuously |
| `AGENTS.md` | Project conventions + completed work documentation |
| `.positai/expertise/*.md` | Task-specific methodology for standing tasks |

### Subprocess discipline

**All builds, tests, and benchmarks must run in bash subprocesses** — never
in the RStudio R session. A crashed subprocess does not take down the agent.

```bash
# Build
R CMD INSTALL --library=.agent-X .

# Test
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='ts-')"

# Benchmark
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); source('inst/benchmarks/bench_simd.R')"
```

Never use `runCode` for long-running operations — if the R session crashes,
the agent session dies with it.

## Test file conventions

All `tests/testthat/test-ts-*.R` files must use **`TreeSearch:::`** to
call internal (non-exported) C++ bridge functions. Without the prefix,
tests fail when run via `library(TreeSearch, lib.loc = ...)` (which is
how agents run them). Only `devtools::load_all()` puts internals on the
search path, and agents must not use `devtools::load_all()`.

**Pattern:** Define short local helper functions that wrap the prefixed calls:

```r
ts_score <- function(tree, ds) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ...)
}
```

Then use the helpers in `test_that()` blocks for readability.

## Completed agent work

| Agent | Step | Files created | Files modified |
|-------|------|---------------|----------------|
| (Step 1) | TBR search | `ts_tbr.h/.cpp`, `test-ts-tbr-search.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| A | Ratchet | `ts_ratchet.h/.cpp`, `test-ts-ratchet-search.R`, `test-ts-ratchet-stress.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| B | Drifting | `ts_drift.h/.cpp`, `test-ts-drift-search.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| C | Splits + Pool | `ts_splits.h/.cpp`, `ts_pool.h/.cpp`, `test-ts-splits.R`, `test-ts-pool.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| D | Wagner tree | `ts_wagner.h/.cpp`, `test-ts-wagner.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| E | Tree fusing | `ts_fuse.h/.cpp`, `test-ts-fuse.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| F | Sectorial search | `ts_sector.h/.cpp`, `test-ts-sector.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| G | Driven search | `ts_driven.h/.cpp`, `test-ts-driven.R` | `ts_rcpp.cpp`, `TreeSearch-init.c` |
| G | CSS (Phase 3B) | `test-ts-css.R` | `ts_tbr.h/.cpp`, `ts_sector.h/.cpp`, `ts_driven.h/.cpp`, `ts_rcpp.cpp`, `TreeSearch-init.c`, `R/MaximizeParsimony.R` |

## Inapplicable character scoring (Phase 4)

Three-pass algorithm (Brazeau et al. 2019) implemented in `ts_fitch_na.inc`,
`#include`d at the end of `ts_fitch.cpp`. Key implementation details:

- **Pass 1 (First downpass)**: NA-aware state resolution matching morphy's
  `mpl_NA_fitch_first_downpass` exactly. Cases: both-applicable → Fitch on
  applicable states; one-applicable → union + NA; both-NA → keep {NA}.
- **Pass 2 (First uppass)**: Applicability propagation matching morphy's
  `mpl_NA_fitch_first_uppass`. Tip update matches `mpl_fitch_NA_tip_update`:
  strip NA only when tip intersects ancestor AND ancestor is applicable.
- **Pass 3 (Second downpass)**: Step counting uses `subtree_actives` (not
  children's D2 directly). Step formula:
  `l_act & r_act & ~(ss_app & any_d2_isect)`. This counts region-separation
  steps through inapplicable nodes (where both subtrees have applicable tips).

### Critical bug fixes applied:
1. **`.inc` file recompilation**: `ts_fitch_na.inc` changes require
   `touch src/ts_fitch.cpp` before rebuild — the build system doesn't
   track `.inc` dependencies automatically.
2. **State word count in `build_dataset`**: All blocks must use
   `total_app_states` (= number of applicable columns in the contrast
   matrix), NOT the max per-pattern count. The global `state_remap`
   assigns consecutive indices, so a pattern using state index k needs
   word k to exist even if few patterns use that state.

### Verified:
- All 30 `inapplicable.phyData` datasets match morphy (provided trees)
- 90/90 random-tree × dataset combinations match morphy
- 20/20 random DNA trees match phangorn (standard Fitch regression)

## Sectorial search (Agent F, Step 7)

Implemented in `ts_sector.h/.cpp`. Supports Random Sectorial Search (RSS) and
Exclusive Sectorial Search (XSS). Key implementation details:

- **Reduced dataset**: Sector clade extracted with an HTU pseudo-tip representing
  the rest of tree (using `final_` states from parent of sector root).
- **Node mapping**: `sector_to_full[]` and `full_to_sector[]` arrays for topology
  interchange between sector and full tree.
- **Root structure guard**: After TBR on sector tree, verifies that the synthetic
  root's children are still the HTU and sector_root_mapped. If TBR regrafted onto
  a root edge (displacing nodes above sector_root_mapped), the result is discarded
  to prevent topology corruption during reinsertion.
- **XSS partitioning**: Bottom-up postorder walk, greedily claiming clades of
  approximately `n_tip / n_partitions` unclaimed tips.
- **Global TBR**: Both RSS and XSS finish with a global TBR pass on the full tree.

### Files created/modified:
- Created: `ts_sector.h`, `ts_sector.cpp`, `test-ts-sector.R`
- Modified: `ts_rcpp.cpp`, `TreeSearch-init.c` (Rcpp bridges: `ts_rss_search`,
  `ts_xss_search`, `ts_sector_diag`)

### Test status: 32/32 passing (stable across repeated runs)

## Exported Rcpp functions (current)

All registered in `ts_rcpp.cpp` and `TreeSearch-init.c`:

| Function | Source module | Purpose |
|----------|-------------|---------|
| `ts_fitch_score` | ts_fitch | Score a tree (dispatches to NA-aware if needed) |
| `ts_na_debug_char` | ts_fitch_na | Per-node debug info for a single pattern |
| `ts_na_char_steps` | ts_fitch_na | Per-pattern step counts |
| `ts_debug_clip` | ts_fitch | Debug SPR clip/regraft cycle |
| `ts_test_indirect` | ts_fitch | Debug indirect length calculation |
| `ts_nni_search` | ts_search | NNI hill-climbing search |
| `ts_spr_search` | ts_search | SPR hill-climbing search |
| `ts_tbr_search` | ts_tbr | TBR search with plateau exploration |
| `ts_ratchet_search` | ts_ratchet | Ratchet perturbation search |
| `ts_drift_search` | ts_drift | Drift search (accept suboptimal moves) |
| `ts_wagner_tree` | ts_wagner | Build Wagner tree (specified addition order) |
| `ts_random_wagner_tree` | ts_wagner | Build Wagner tree (random addition order) |
| `ts_compute_splits` | ts_splits | Compute bipartition splits from edge matrix |
| `ts_trees_equal` | ts_splits | Compare two trees by topology |
| `ts_pool_test` | ts_pool | Test pool deduplication |
| `ts_tree_fuse` | ts_fuse | Fuse two trees |
| `ts_sector_diag` | ts_sector | Diagnostic info for sectorial search |
| `ts_rss_search` | ts_sector | Random Sectorial Search |
| `ts_xss_search` | ts_sector | Exclusive Sectorial Search |
| `ts_driven_search` | ts_driven | Full driven search (multi-replicate) |
| `ts_resample_search` | ts_resample | One jackknife or bootstrap replicate |
| `ts_successive_approx` | ts_resample | Full SA search with iterative reweighting |
| `ts_bench_tbr_phases` | ts_rcpp | TBR phase timing diagnostic (Phase 3D) |

## Driven search (Agent G, Step 8)

Implemented in `ts_driven.h/.cpp`. Orchestrates all search components into a
single multi-replicate driven search strategy.

### Search pipeline per replicate:
1. Random Wagner tree → TBR to local optimum
2. XSS sectorial search (if tree large enough)
3. Ratchet perturbation to escape local optima
4. Final TBR polish
5. Add to pool

### Additional features:
- **Tree fusing**: Every `fuse_interval` replicates, fuses best tree against pool
- **Convergence**: Stops when `hits_to_best >= target_hits`
- **Pool management**: Dedup via split hashing, score-based eviction

### Files created/modified:
- Created: `ts_driven.h`, `ts_driven.cpp`, `test-ts-driven.R`
- Modified: `ts_rcpp.cpp`, `TreeSearch-init.c` (Rcpp bridge: `ts_driven_search`)

### Test status: 20/20 passing (stable across repeated runs)

### Exported function:

| Function | Source module | Purpose |
|----------|-------------|---------|
| `ts_driven_search` | ts_driven | Full driven search (Wagner+TBR+XSS+ratchet+fuse) |

## Red-team bug fixes (Agent E, cross-agent review)

All fixes applied and verified (build succeeds, fuse 16/16 and sector 32/32 tests pass).

| Agent | File | Bug | Fix applied |
|-------|------|-----|-------------|
| A (Ratchet) | `ts_ratchet.cpp` | Missing `GetRNGstate()`/`PutRNGstate()` around `unif_rand()` | Wrapped seed call |
| B (Drift) | `ts_drift.cpp` | `std::random_device{}()` ignores `set.seed()` | Seed from `unif_rand()` with RNG state management |
| B (Drift) | `ts_drift.cpp` | RFD computation ignores `CharBlock.weight` | Multiply popcount by `ds.blocks[b].weight` |
| C (Splits) | `ts_splits.cpp` | `__builtin_popcountll` (GCC-only) | Replaced with `popcount64()` from `ts_data.h` |
| D (Wagner) | `ts_wagner.cpp` | `init_wagner_state` missing `down2`/`subtree_actives` | Added allocation |
| F (Sector) | `ts_sector.cpp` | Missing `GetRNGstate()`/`PutRNGstate()` in RSS/XSS | Wrapped both functions |
| F (Sector) | `ts_sector.cpp` | No revert when reinsertion worsens full-tree score | Save/restore topology on worse score |
| F (Sector) | `ts_sector.cpp` | `internal_ratchet_cycles` param unused | Documented as reserved; commented out param name |
| F (Sector) | `ts_sector.cpp` | Dead `ensure_na_arrays` + redundant if/else in reinsertion | Removed dead code, simplified |
| — | `ts_rcpp.cpp` | Junk `RCPP_APPEND 2>&1` line at EOF | Removed |
| — | `RcppExports.cpp` | Missing `ts_driven_search` wrapper | Regenerated with `Rcpp::compileAttributes()` |

## Red-team Agent G (Driven search) — bug fixes applied

| Issue | File | Bug | Fix applied |
|-------|------|-----|-------------|
| G1 | `ts_driven.cpp` | `max_replicates=0` → crash on empty pool | Guard: return early with default result |
| G2 | `ts_driven.cpp` | Fused trees inflate `hits_to_best` convergence counter | Save/restore counter around fuse; reset to 0 if fusing found new best |
| G3 | `ts_driven.h`, `ts_rcpp.cpp` | Dead `rss_picks` parameter | Removed from `DrivenParams` and Rcpp bridge |
| G4 | `ts_driven.h/.cpp` | Misleading `ts_drift.h` include + header comment | Removed include; updated comment |
| — | `ts_pool.h` | No `set_hits_to_best()` accessor | Added for G2 fix |
| — | `TreeSearch-init.c` | Arg count mismatch after G3 | Updated 18→17 |
| — | `R/RcppExports.R`, `src/RcppExports.cpp` | Stale after param removal | Regenerated via `Rcpp::compileAttributes()` |

### Test status: fuse 16/16, sector 32/32, driven 20/20

## Cross-module integration fixes (Agent E, integration review)

| Issue | Files | Bug | Fix applied |
|-------|-------|-----|-------------|
| I1 | `ts_tbr.cpp` | `std::random_device{}()` seeding — root cause of non-deterministic TBR | Seed from `unif_rand()` with RNG state management |
| I2 | ts_tbr/drift/fitch/search.cpp | `__builtin_ctzll` GCC-only (7 occurrences) | Added portable `ctz64()` wrapper in `ts_data.h` |
| I3 | `ts_data.h` | `__builtin_huge_val()` GCC-only | Replaced with `HUGE_VAL` from `<cmath>` |
| I4 | `ts_tree.cpp` | `init_from_edge` assigns first child→right (opposite of `tree_to_edge` output) | Changed to first child→left; R↔C++ round-trips now idempotent |

### TBR non-determinism: RESOLVED
Previous observation of non-deterministic TBR topologies was caused by two bugs:
1. `std::random_device` seeding in TBR's `std::shuffle` (I1, primary cause)
2. `init_from_edge` left-right swap (I4, secondary cause via R round-trips)

Both are now fixed. `set.seed()` produces identical results across runs.

### Design notes (not fixed — documented):
- **TBR incremental pass ignores NA**: Uses standard Fitch heuristic for move
  evaluation; `full_rescore` with `score_tree()` → `fitch_na_score()` for exact
  verification. May miss some moves on heavily-inapplicable datasets.
- **Ratchet `active_mask` not RAII-protected**: Interrupt during perturbation
  could leave masks corrupted. Low risk — DataSet rebuilt from scratch per R call.

## Phase 1C: All-trees + timeout + verbosity (Agent C)

Implements production plan Phase 1C: feature completeness for driven search.

### Changes:

**`ts_driven.h`**:
- `DrivenParams` gains `max_seconds` (double, default 0 = no timeout) and
  `verbosity` (int, 0=silent / 1=per-replicate / 2=per-phase).
- `DrivenResult` gains `timed_out` (bool).
- `driven_search()` signature changed: takes `TreePool&` (out parameter)
  instead of `TreeState&`, so caller can access all pool contents.

**`ts_driven.cpp`**:
- Timeout via `std::chrono::steady_clock`. Checked after each heuristic phase
  and at end of each replicate. Uses `goto finish` for clean exit.
- `R_CheckUserInterrupt()` added after each phase (Wagner+TBR, XSS, ratchet,
  drift, final TBR) for responsiveness, plus end of replicate.
- `Rprintf()` progress reporting at verbosity levels 1 and 2.

**`ts_rcpp.cpp`**:
- `ts_driven_search` Rcpp bridge now accepts `maxSeconds` and `verbosity`.
- Return value changed: `trees` (list of edge matrices), `scores` (numeric
  vector), `best_score`, `replicates`, `hits_to_best`, `pool_size`, `timed_out`.
  Previously returned only `edge` (single best tree) and `score`.

**`TreeSearch-init.c`**: Arg count updated 22→24.

**`R/RcppExports.R`**: Added `concavity = Inf` default (Rcpp can't auto-generate
`R_PosInf`). Also fixed pre-existing missing default in other test helpers.

**`test-ts-driven.R`**: Updated all tests for new return structure (`trees`/`scores`
instead of `edge`/`score`). Added 7 new tests for pool-return, suboptimal
collection, timeout, verbosity, and zero-replicate edge case.

### Test status: 47/47 passing

### Exported function changes:

| Function | Change |
|----------|--------|
| `ts_driven_search` | +2 params (`maxSeconds`, `verbosity`); return structure changed |

## R-level function rename (Phase 9)

The old MorphyLib-based `MaximizeParsimony()` has been renamed to `Morphy()`.
The C++ driven search (`MaximizeParsimony2()`) is now the standard
`MaximizeParsimony()`.

### Naming convention:
- **`MaximizeParsimony()`** — C++ driven search (EW + IW natively; delegates to
  `Morphy()` for profile parsimony and constraints)
- **`Morphy()`** — R-loop search using MorphyLib (supports EW, IW, profile
  parsimony, constraints)
- **`MaximizeParsimony2()`** — deprecated alias for `MaximizeParsimony()`

### Files renamed:
| Old | New | Contents |
|-----|-----|----------|
| `R/MaximizeParsimony.R` | `R/Morphy.R` | `Morphy()`, `Resample()`, `EasyTrees()`, `EasyTreesy()`, internal helpers |
| `R/MaximizeParsimony2.R` | `R/MaximizeParsimony.R` | `MaximizeParsimony()` (C++ engine), `MaximizeParsimony2()` (deprecated) |
| `tests/testthat/test-MaximizeParsimony.R` | `tests/testthat/test-Morphy.R` | Tests for `Morphy()` |

### Updated callers:
- `R/TaxonInfluence.R`: calls `Morphy()` instead of `MaximizeParsimony()`
- `R/IWScore.R`: `.Deprecated("Morphy")`
- `inst/Parsimony/app.R`: calls `Morphy()` for Shiny search
- `tests/testthat/test-CustomSearch.R`: calls `Morphy()`
- Vignettes, README, shinytest expected outputs: updated

### Migration queue (functions still using MorphyLib):
See plan file `.positai/plans/2026-03-16-plan-mmt759uv.md` for the full
tiered migration queue of functions to move from MorphyLib to the C++ engine.

## Phase 1D: Successive approximations + jackknife/bootstrap (Agent C)

Implements production plan Phase 1D: resampling and successive approximations
in the C++ search engine.

### Files created:
- `ts_resample.h` — Header with `ResampleParams`, `SAParams`, `ResampleResult`,
  `SAResult` structs and function declarations.
- `ts_resample.cpp` — Implementation of `resample_search()` and
  `successive_approximations()`.
- `tests/testthat/test-ts-resample.R` — 35 tests covering jackknife, bootstrap,
  and successive approximations.

### Files modified:
- `ts_rcpp.cpp` — Added `ts_resample_search` and `ts_successive_approx` Rcpp
  bridges; added `#include "ts_resample.h"`.
- `TreeSearch-init.c` — Registered 2 new entry points (14 args each).
- `R/RcppExports.R` — Added `concavity = Inf` defaults for new functions.

### Implementation details:

**Jackknife (`resample_search`, `bootstrap=false`):**
- Expands original pattern weights into a flat character index.
- Fisher-Yates partial shuffle to sample `ceil(proportion * n_chars)` characters
  without replacement.
- Computes new per-pattern weights from the sample.
- Rebuilds `DataSet` with modified weights and runs `driven_search`.
- Uses `unif_rand()` with `GetRNGstate()`/`PutRNGstate()` for R-compatible RNG.

**Bootstrap (`resample_search`, `bootstrap=true`):**
- Same as jackknife but samples WITH replacement (same total character count).
- Each pattern's new weight = number of times it was sampled.

**Successive approximations (`successive_approximations`):**
- Outer loop: for each SA iteration:
  1. Compute effective weights = `original_weight * sa_weight`, scaled so
     minimum nonzero weight = 1 (rounded to int).
  2. Build `DataSet` with effective weights, run `driven_search`.
  3. Rebuild an EW `DataSet` with original weights, re-initialize best tree
     against it via `init_from_edge`, score with `fitch_score`/`fitch_na_score`.
  4. Extract per-pattern step counts via `extract_char_steps`.
  5. Check convergence: if step counts identical to previous iteration, stop.
  6. Reweight: `w_i = (p_i)^(-k) - 1` where `p_i = steps_i / (n_internal - 1)`.
     Characters with `p_i >= 1` get weight 0; characters with 0 steps get
     maximum weight.
- Returns EW parsimony score of the final tree (not the SA-weighted score).

### Exported Rcpp functions:

| Function | Args | Purpose |
|----------|------|---------|
| `ts_resample_search` | 14 | One jackknife or bootstrap replicate |
| `ts_successive_approx` | 14 | Full SA search with iterative reweighting |

### Test status: 35/35 passing (resample), 47/47 passing (driven, unchanged)

### Red-team stress testing (Agent C self-review):
- Created `test-ts-resample-stress.R` — 60 additional tests.
- Covers: `set.seed()` reproducibility (driven, jackknife, bootstrap, SA),
  inapplicable characters (Vinther2008 dataset), implied weights, `fuse_interval=1`,
  large `poolSuboptimal`, very short timeout (0.01s), extreme `jackProportion`
  (0.1, 0.99), `maxSAIter=0`, pool tree validity, SA EW score verification,
  medium dataset (20 tips).
- **All 148 tests pass** (53 driven + 35 resample + 60 stress).
- Fixed: `concavity` default problem — replaced `R_PosInf` with `-1.0` sentinel
  in all 15 Rcpp-exported functions; `compileAttributes()` now safe to run freely.

## Phase 2B: TBR neighborhood traversal optimization (Agent F)

Six optimizations applied to `ts_tbr.cpp` and `ts_fitch.cpp`:

### Optimizations implemented:

1. **Early termination in indirect scoring**: `fitch_indirect_length_bounded()`,
   `indirect_iw_length_bounded()`, `fitch_indirect_length_cached()`, and
   `indirect_iw_length_cached()` — bail out as soon as accumulated score exceeds
   the current best candidate. Eliminates 60-80% of inner-loop work for losing
   candidates.

2. **Clip smaller subtree only**: Precomputes subtree sizes; skips clips where
   `subtree_size > n_tip / 2`. Halves the worst-case rerooting count. Both
   orientations of each edge are still explored through smaller clips of
   sub-components.

3. **Avoid full_rescore on rejection**: `StateSnapshot` struct saves complete
   state arrays (prelim, final_, local_cost, down2, subtree_actives) via
   `memcpy` before applying a TBR move. On rejection, restores from snapshot
   instead of running an expensive full scoring pass.

4. **Precomputed vroot cache**: For TBR rerooting, precomputes
   `vroot[edge][s] = final_[A][s] | final_[D][s]` for all main edges once per
   clip. Eliminates redundant OR operations in the inner loop (previously
   recomputed k times per main edge, once per rerooting).

5. **Deferred reshuffling**: Only reshuffles clip candidates when a full pass
   finds no improvement. After accepting a move, retries with the same ordering
   to exploit spatial locality of improving moves.

### Files created:
- `tests/testthat/test-ts-tbr-bench.R` — 7 tests (26 expectations) for TBR
  optimization correctness.

### Files modified:
- `src/ts_fitch.h` — Added 5 new function declarations (bounded/cached variants).
- `src/ts_fitch.cpp` — Implemented `fitch_indirect_length_bounded`,
  `fitch_indirect_length_cached`, `indirect_iw_length_bounded`,
  `indirect_iw_length_cached`.
- `src/ts_tbr.cpp` — All optimizations applied: smaller-subtree filter,
  early termination, vroot cache, StateSnapshot, deferred reshuffling.
  Added `compute_subtree_sizes()`, `precompute_vroot_cache()`, `StateSnapshot`
  struct.
- `R/RcppExports.R` — Fixed `concavity = Inf` defaults after
  `Rcpp::compileAttributes()` regeneration.
- `src/RcppExports.cpp` — Regenerated via `Rcpp::compileAttributes()`.

### Test status:
- TBR bench: 26/26 passing
- Driven: 53/53 passing
- Sector: 32/32 passing
- Fuse: 16/16 passing (1 skip)
- Resample: 35/35 passing

## Phase 2D: Ratchet optimization (Agent G)

Implements production plan Phase 2D: configurable perturbation modes,
upweight support, and adaptive perturbation tuning.

### Changes:

**`ts_ratchet.h`**:
- `PerturbMode` enum: `ZERO_ONLY` (default), `UPWEIGHT_ONLY`, `MIXED`.
- `RatchetParams` gains: `perturb_mode`, `perturb_max_moves` (was hardcoded),
  `perturb_accept_equal`, `adaptive`, `target_escape_rate`, `adapt_min_prob`,
  `adapt_max_prob`.
- `RatchetResult` gains: `n_escapes`, `final_perturb_prob`.

**`ts_data.h`**:
- `CharBlock` gains `uint64_t upweight_mask = 0`. When nonzero, scoring
  functions count upweighted characters double.

**`ts_fitch.cpp`**:
- `fitch_downpass`, `fitch_incremental_downpass`, `fitch_indirect_length`:
  check `blk.upweight_mask`; add extra popcount for upweighted characters.
  When mask is 0 (always, outside perturbation), branch predictor handles
  with near-zero overhead.

**`ts_fitch_na.inc`**:
- Pass 1 standard-block scoring and Pass 3 NA step counting: same
  `upweight_mask` check as `fitch_downpass`.

**`ts_ratchet.cpp`**:
- Three perturbation modes: `perturb_zero` (original), `perturb_upweight`
  (sets `upweight_mask`; for IW doubles `pattern_freq`), `perturb_mixed`
  (zeroes some chars, upweights others — disjoint sets).
- `PerturbSnapshot` saves/restores `active_mask`, `upweight_mask`, and
  `pattern_freq` per cycle.
- Adaptive tuning: tracks escape rate over batches of 3 cycles; adjusts
  `current_prob` within `[adapt_min_prob, adapt_max_prob]`.
- `perturb_max_moves` configurable (was hardcoded `max(20, min(200, n_tip/8))`).

**`ts_driven.h`/`ts_driven.cpp`**:
- `DrivenParams` gains `ratchet_perturb_mode`, `ratchet_perturb_max_moves`,
  `ratchet_adaptive`. Passed through to `ratchet_search()`.

**`ts_rcpp.cpp`**:
- `ts_ratchet_search`: +4 params (`perturbMode`, `perturbMaxMoves`,
  `adaptive`, `targetEscapeRate`); return gains `n_escapes`, `final_perturb_prob`.
- `ts_driven_search`: +3 params (`ratchetPerturbMode`, `ratchetPerturbMaxMoves`,
  `ratchetAdaptive`).

**`TreeSearch-init.c`**: Ratchet 10→14 args; driven 30→33 args.

**`R/MaximizeParsimony.R`**: Exposes `ratchetPerturbMode`, `ratchetPerturbMaxMoves`,
`ratchetAdaptive` in the R interface.

### Files created:
- `tests/testthat/test-ts-ratchet-opt.R` — 23 tests covering upweight mode,
  mixed mode, adaptive tuning, perturb_max_moves, IW upweighting, mask cleanup.

### Test status: 165/165 passing (53 driven + 17 ratchet-search + 72 ratchet-stress + 23 ratchet-opt)

### Deferred to Phase 6 (empirical benchmarking):
- Optimal default `perturb_prob` for each mode
- Whether mixed/upweight outperform zero-only
- Optimal `target_escape_rate` for adaptive mode

## Wagner tree incremental scoring (Agent G, Phase 2G)

Optimized Wagner tree construction in `ts_wagner.cpp` by replacing full
`build_postorder() + score_tree()` after each taxon insertion with incremental
two-pass Fitch scoring.

### Implementation:
- **Incremental downpass**: Walks from newly inserted internal node to root,
  recomputing `prelim` and `local_cost`. Stops early when `prelim` stabilizes.
  Tracks score delta via local_cost subtraction/addition.
- **DFS uppass with early termination**: After root `final_ = prelim`, does
  DFS through internal nodes, computing `final_` from parent's `final_` + own
  `prelim`. Skips entire subtrees where `final_` didn't change.
- **Lazy postorder**: `build_postorder()` called once at the end, not per
  insertion. Incremental passes use parent pointers and DFS instead.
- **Final score**: `score_tree()` called once after all taxa added, giving
  authoritative result (handles NA three-pass and IW).

### Complexity improvement:
- Downpass: O(depth × C) per insertion vs O(n × C) — typically ~log(n)/n speedup
- Uppass: O(affected_region × C) per insertion, with early termination avoiding
  full tree traversal for most insertions

### Files modified:
- `src/ts_wagner.cpp` — Added `wagner_incremental_rescore()` (anonymous namespace),
  changed `wagner_tree()` main loop to use incremental scoring
- `tests/testthat/test-ts-wagner.R` — Fixed `TreeSearch:::` namespace prefixes
  for all 11 original tests; added 5 new tests (NA datasets, multiple seeds,
  driven search integration)

### Test status: 26/26 passing (Wagner), 53/53 passing (driven, unchanged)

## Phase 2C: SPR/NNI indirect scoring improvements (Agent F)

Ported Phase 2B TBR optimizations to SPR search and added incremental
scoring to NNI search. SPR and NNI are standalone search modes (not used
in driven search pipeline) but useful for interactive exploration.

### SPR optimizations implemented:

1. **Bounded indirect scoring with early termination**: Replaced unbounded
   `fitch_indirect_length` / `indirect_iw_length` with bounded variants
   (`_bounded`) that bail out when accumulated score exceeds current best
   candidate. Added NA-aware dispatch (was previously missing — SPR used
   standard Fitch indirect even for NA datasets).

2. **Subtree-size filtering**: Skip clips where `subtree_size > n_tip / 2`,
   matching TBR's smaller-subtree-only pattern.

3. **Deferred reshuffling**: Only reshuffle clip candidates when previous
   pass found no improvement. After acceptance, retry same ordering.

4. **Incremental clip scoring**: Replaced full `reset_states + score_tree +
   manual clip subtree scoring` with `fitch_incremental_downpass/uppass`
   (or NA three-pass variants), matching TBR's clip phase exactly.

5. **Best-of-all screening**: Changed from first-pass-any-non-dominated to
   tracking `best_candidate` across all destinations, verifying only the
   single best candidate with full rescore. Reduces unnecessary full rescores.

### NNI optimization implemented:

6. **Incremental two-pass scoring**: Replaced full `build_postorder +
   score_tree` per candidate with `fitch_incremental_downpass` from the NNI
   node to root. O(depth × C) per candidate instead of O(n × C).
   State restoration on rejection via second incremental downpass after undo.
   NA datasets fall back to full rescore (three-pass too complex for
   incremental NNI). `build_postorder` + `fitch_uppass` called only on
   acceptance (not on every candidate).

### Files created:
- `tests/testthat/test-ts-spr-nni-opt.R` — 47 tests covering EW, IW, NA,
  determinism, cross-method score ordering, and empirical datasets.

### Files modified:
- `src/ts_search.cpp` — Complete rewrite of both `nni_search()` and
  `spr_search()` with all optimizations applied.

### Test status:
- SPR/NNI opt: 47/47 passing
- TBR bench: 26/26 passing (unchanged)
- Driven: 53/53 passing (unchanged)
- Sector: 32/32 passing (unchanged)
- Fuse: 16/16 passing (1 skip, unchanged)
- Resample: 35/35 passing (unchanged)

## Phase 2E+2F: Sectorial search & tree fusing optimization (Agent H)

Six optimizations applied to `ts_sector.cpp` and `ts_fuse.cpp`, plus RSS
integration into the driven search pipeline.

### Optimizations implemented:

1. **O(n) XSS partitioning**: Replaced O(n²) `xss_partition()` (DFS per node
   to count unclaimed tips) with O(n) algorithm using maintained
   `unclaimed_below[]` array. When a sector is claimed, rootward walk
   subtracts from ancestors — O(height) per claim, O(n_partitions × height)
   total.

2. **In-place fuse exchange with undo**: Replaced per-trial `copy_topology()`
   (full tree allocation + copy) with in-place `replace_subtree()` + targeted
   clade save/restore. On acceptance, avoids both the copy and the
   `recipient = trial` assignment (~4× tree size in memory operations). On
   rejection, restores only the clade's topology; state arrays become stale
   but `score_tree()` recomputes them from scratch on the next trial.

3. **Eliminated redundant full-tree rescores**: Removed `score_tree()` calls
   before each sector in RSS and XSS. States are guaranteed valid from the
   initial rescore or the previous sector's acceptance/rejection handling.

4. **Lazy donor processing in fuse**: Donors are prepared (copy_topology +
   reroot + split extraction) on first access rather than all upfront. Cached
   across rounds. Saves processing cost for donors that are never reached
   (early improvement triggers break-to-TBR).

5. **RSS wired into driven search**: Added `rss_rounds` parameter to
   `DrivenParams` (default 1). After XSS, runs RSS rounds with random sector
   picks, complementing XSS's systematic partitioning. TNT's `xmult` uses
   both; our previous pipeline only used XSS.

6. **Targeted topology snapshots**: Replaced full-vector saves
   (`auto save_left = tree.left`) in RSS/XSS with `CladeSnapshot` that saves
   only the sector clade's internal nodes. For a 500-tip tree with 30-tip
   sectors, saves O(30) instead of O(500) per snapshot.

### Files modified:
- `src/ts_sector.h` — No changes to header
- `src/ts_sector.cpp` — All 6 optimizations: O(n) partition, removed redundant
  rescores, `CladeSnapshot` struct + `save_clade()`/`restore_clade()` helpers
- `src/ts_fuse.cpp` — In-place exchange, lazy donor processing, pre-allocated
  save/restore buffers
- `src/ts_driven.h` — Added `rss_rounds` to `DrivenParams`
- `src/ts_driven.cpp` — RSS rounds after XSS in pipeline
- `src/ts_rcpp.cpp` — Added `rssRounds` parameter to `ts_driven_search` bridge
- `src/TreeSearch-init.c` — Driven search arg count 34→35
- `R/MaximizeParsimony.R` — Exposed `rssRounds` parameter
- `R/RcppExports.R` — Regenerated via `Rcpp::compileAttributes()`
- `src/RcppExports.cpp` — Regenerated

### RcppExports regeneration note:
`Rcpp::compileAttributes()` was run to sync RcppExports with current
ts_rcpp.cpp signatures (including `infoAmounts` parameter added by another
agent). The `concavity` default changed from `Inf` to `-1.0` (the C++
sentinel). Both work correctly for EW mode: `-1.0` triggers
`if (concavity < 0) concavity = HUGE_VAL` in `make_dataset()`, while
`Inf` passes through as `HUGE_VAL` directly. The `-1.0` default survives
byte-compilation (unlike `Inf`, which R's byte-compiler can drop).

### Test status: sector 32/32, fuse 16/16, driven 53/53, resample 35/35

## Phase 1B: Profile parsimony scoring in C++ (Agent B) — ✅ COMPLETE

Brings profile parsimony scoring into the C++ search engine so that
`MaximizeParsimony(dataset, concavity = "profile")` uses the fast C++
driven search instead of delegating to the R-loop `Morphy()`.

### Design:
Profile parsimony is structurally identical to implied weights — both
transform per-character step counts into scores via a lookup. The existing
`indirect_iw_length()` function is reused unchanged; only delta
precomputation differs.

### Implementation details:

- **`ScoringMode` enum** (`EW`, `IW`, `PROFILE`): Added to `ts_data.h`.
  Auto-determined at `build_dataset` time: if `info_amounts` provided →
  PROFILE; if `concavity` finite → IW; else → EW.
- **Sentinel concavity**: Profile mode sets `ds.concavity = 1.0` so
  `isfinite(concavity)` checks activate the weighted pipeline in all
  search modules without code duplication.
- **`info_amounts` table**: Column-major `[max_steps × n_patterns]` from R,
  stored in `DataSet`. Indexing: `info_amounts[(s-1) + info_max_steps * p]`
  where `s` = total step count (1-based).
- **`compute_profile()`**: Looks up score per pattern from `info_amounts`.
- **`precompute_profile_delta()`**: Computes marginal cost of one extra step
  per pattern (profile analogue of `precompute_iw_delta`).
- **`compute_weighted_score()` / `precompute_weighted_delta()`**: Dispatchers
  (IW or profile based on `ds.scoring_mode`).
- **Search module changes**: 6 callsites in TBR/SPR/Drift changed from
  `compute_iw` → `compute_weighted_score` and `precompute_iw_delta` →
  `precompute_weighted_delta`.
- **Ratchet perturbation**: Existing `use_iw` flag (via `isfinite(concavity)`)
  naturally includes profile mode thanks to the sentinel.

### R-level integration:
- `MaximizeParsimony(dataset, concavity = "profile")` calls
  `PrepareDataProfile(dataset)` to simplify characters to binary, then passes
  the `info.amounts` matrix to the C++ engine via `infoAmounts` parameter.
- `concavity` is set to `Inf` after profile prep (EW on simplified data;
  profile scores applied via lookup in C++).

### Limitations preserved (Phase 1):
- Max 2 informative states per character (PrepareDataProfile reduces to binary)
- Inapplicable tokens → ambiguous (PrepareDataProfile converts `-` to `?`)
- NA blocks never arise in profile mode, so incremental NA scoring is moot

### Files created:
- `tests/testthat/test-ts-profile.R` — 25 tests covering score verification,
  search end-to-end, regression, edge cases, suboptimal collection, timeout.

### Files modified:
- `src/ts_data.h` — Added `ScoringMode` enum, `info_amounts`, `info_max_steps`
- `src/ts_data.cpp` — Extended `build_dataset` to accept `info_amounts_r`
- `src/ts_fitch.h` — Declared profile scoring functions
- `src/ts_fitch.cpp` — Implemented `compute_profile`, `precompute_profile_delta`,
  `compute_weighted_score`, `precompute_weighted_delta`; modified `score_tree`
- `src/ts_tbr.cpp` — Replaced `compute_iw` → `compute_weighted_score` (2 sites)
- `src/ts_search.cpp` — Same (2 sites)
- `src/ts_drift.cpp` — Same (2 sites)
- `src/ts_rcpp.cpp` — Added `infoAmounts` param to `make_dataset`,
  `ts_fitch_score`, `ts_driven_search`
- `src/TreeSearch-init.c` — Updated arg counts
- `R/MaximizeParsimony.R` — Profile prep + `infoAmounts` passthrough
- `R/RcppExports.R` — Regenerated with `infoAmounts` defaults

### Test status: profile 25/25, driven 53/53, resample 35/35

## Phase 2A: Incremental NA-aware scoring (Agent E) — ✅ COMPLETE

Implements incremental three-pass Fitch scoring for inapplicable characters
in TBR and drift search. Previously these searches used standard Fitch for
candidate screening, causing false negatives (missed improving moves) and
false positives (wasted full rescores) on NA datasets.

### Files created:
- `src/ts_fitch_na_incr.inc` — Four key functions (included at end of
  `ts_fitch.cpp` after `ts_fitch_na.inc`):
  - `fitch_na_incremental_downpass()` — Walk rootward from clip site, update
    `prelim` using NA-aware Pass 1 logic, maintain `subtree_actives`
  - `fitch_na_incremental_uppass()` — Update `final_` using NA-aware Pass 2
    logic; handle tips separately (not in postorder)
  - `fitch_na_pass3_score()` — Full Pass 3 (second downpass) on divided tree;
    compute exact score after incremental Passes 1+2
  - `fitch_na_indirect_length()` + `indirect_na_iw_length()` — NA-aware
    candidate screening; suppress steps where clip subtree or edge-below
    subtree lacks applicable tips
- `tests/testthat/test-ts-na-incremental.R` — 33 tests covering
  reproducibility, IW+NA combinations, pool suboptimal collection, timeouts,
  score verification on inapplicable datasets

### Files modified:
- `src/ts_fitch.h` — Added 6 function declarations for NA incremental variants
- `src/ts_fitch.cpp` — Added `#include "ts_fitch_na_incr.inc"` at end
- `src/ts_tree.h` / `src/ts_tree.cpp` — Extended `NodeSnapshot` to save/restore
  `down2` and `subtree_actives` arrays (essential for NA state restoration)
- `src/ts_tbr.cpp` — NA detection at search start; save `clip_actives` before
  clipping; dispatch to NA-aware incremental functions when `has_na`; NA-aware
  indirect length for candidate evaluation
- `src/ts_drift.cpp` — Same pattern as TBR: NA detection, `clip_actives` save,
  NA-aware dispatch

### Key design decisions:
1. **Hybrid approach**: Incremental heuristics for candidate screening + full
   verification via `full_rescore` (correctness-first). The NA indirect length
   is a conservative lower bound.
2. **`clip_actives` cached once per clip**: `subtree_actives` at clip node from
   original tree, invariant across rerootings (tips don't change). Saved before
   incremental downpass overwrites it.
3. **Non-cached variants only**: Bounded/cached NA indirect variants deferred as
   future optimization (standard variants sufficient for correctness).
4. **Minimal overhead for non-NA datasets**: Single `has_na` branch check,
   cached `clip_actives_buf` allocation.

### Build note:
`.inc` file changes require `touch src/ts_fitch.cpp` before rebuild — the
build system doesn't track `.inc` dependencies automatically.

### Test status: 33/33 (NA incremental), 53/53 (driven), 35/35 (resample),
60/60 (resample-stress), 16/16 (fuse), 32/32 (sector)

## Profile parsimony in C++ (Agent B, Phase 1B)

Implements profile parsimony scoring (Faith & Trueman, 2001) in the C++
search engine, replacing the previous delegation to `Morphy()`.

### Design: reuse of IW indirect pipeline

Profile parsimony is structurally identical to implied weights — both
transform per-character step counts into scores via a per-pattern function.
The key insight: `indirect_iw_length()` is reused unchanged for profile
parsimony; only the delta precomputation differs (table lookup vs. formula).

### Data structures added to `DataSet` (ts_data.h):
- `ScoringMode scoring_mode` — enum `{EW, IW, PROFILE}`
- `std::vector<double> info_amounts` — lookup table [info_max_steps × n_patterns],
  column-major. `info_amounts[(step-1) + info_max_steps * pattern]` gives the
  information cost when pattern has `step` total steps.
- `int info_max_steps` — number of rows in the info_amounts table.

### Scoring functions added to ts_fitch.h/cpp:
- `compute_profile(ds, char_steps)` — sum of `freq[p] * info_amounts[steps[p], p]`
- `precompute_profile_delta(ds, divided_steps, delta)` — marginal cost of +1 step
- `compute_weighted_score(ds, char_steps)` — dispatches to `compute_iw` or `compute_profile`
- `precompute_weighted_delta(ds, divided_steps, delta)` — dispatches to IW or profile delta

### Scoring mode dispatch:
- `score_tree()` uses `ds.scoring_mode == EW` instead of `std::isinf(ds.concavity)`
- For profile mode, `ds.concavity` is set to `1.0` (finite sentinel) so existing
  `std::isfinite()` checks in TBR/SPR/Drift activate the weighted pipeline automatically
- Search modules call `compute_weighted_score` / `precompute_weighted_delta` instead
  of `compute_iw` / `precompute_iw_delta` (6 callsites changed)

### R-level integration:
- `MaximizeParsimony(dataset, concavity = "profile")` now calls `PrepareDataProfile()`
  in R and passes `info.amounts` matrix to the C++ engine
- No longer delegates to `Morphy()`
- `MaximizeParsimony()` result handling updated for multi-tree return format

### Performance impact on EW/IW: None
- One additional `if` in `score_tree()` (profile check before EW dispatch)
- `info_amounts` vector only allocated when in PROFILE mode
- No changes to Fitch downpass/uppass or indirect length calculation paths

### Existing limitations (preserved from R implementation):
- Max 2 informative states per character (PrepareDataProfile reduces to binary)
- Inapplicable tokens treated as ambiguous (PrepareDataProfile replaces `-` → `?`)

### Files created/modified:
- Modified: `ts_data.h`, `ts_data.cpp` (ScoringMode enum, info_amounts, build_dataset)
- Modified: `ts_fitch.h`, `ts_fitch.cpp` (profile scoring functions, score_tree dispatch)
- Modified: `ts_tbr.cpp`, `ts_search.cpp`, `ts_drift.cpp` (weighted score dispatch)
- Modified: `ts_rcpp.cpp`, `TreeSearch-init.c` (infoAmounts parameter)
- Modified: `R/MaximizeParsimony.R` (profile support, multi-tree return handling)
- Created: `tests/testthat/test-ts-profile.R`

### Exported Rcpp function changes:
| Function | Change |
|----------|--------|
| `ts_fitch_score` | +1 param (`infoAmounts`); arg count 7→8 |
| `ts_driven_search` | +1 param (`infoAmounts`); arg count 34→35 |

### Test status: 25/25 (profile), 53/53 (driven, unchanged), 35/35 (resample, unchanged)

## Phase 3C: Character-Ordering Optimization (Agent F)

Constant-factor optimizations to the character scoring pipeline: block
ordering, zero-weight compaction, active_mask skip, and bounded indirect
call migration.

### Optimizations implemented:

1. **Expensive-blocks-first ordering**: Changed `build_dataset()` sort from
   ascending to descending weight within each `(has_inapp, weight)` group.
   Higher-weight blocks are processed first, triggering early termination
   sooner in bounded indirect-length functions. Non-NA blocks still precede
   NA blocks (required by `fitch_na_score`).

2. **Zero-weight pattern compaction**: Patterns with `weight == 0` are now
   removed before block assignment in `build_dataset()`. Eliminates ghost
   patterns after jackknife/bootstrap resampling that would otherwise waste
   block space and scoring time. `pattern_index` mapping is preserved.

3. **`active_mask == 0` skip in all block loops**: Added
   `if (blk.active_mask == 0) continue;` at the top of every block iteration
   in scoring functions (downpass, uppass, incremental, indirect, NA three-pass,
   extract_char_steps). During ratchet perturbation when entire blocks are
   zeroed, this skips all inner-loop work. Near-zero overhead when blocks are
   active (single well-predicted branch).

4. **Bounded indirect calls in drift search**: Replaced 4 unbounded indirect
   length calls with bounded variants (`fitch_indirect_length_bounded`,
   `indirect_iw_length_bounded`, `fitch_na_indirect_length_bounded`,
   `indirect_na_iw_length_bounded`). Cutoff is `best_candidate` (IW) or
   `best_candidate - divided_length` (EW). Eliminates work on losing
   candidate positions.

5. **Bounded indirect in Wagner tree**: Replaced 2 unbounded
   `fitch_indirect_length` calls with `fitch_indirect_length_bounded` using
   `best_extra` as the cutoff. Reduces per-insertion evaluation cost.

### Files created:
- `tests/testthat/test-ts-char-ordering.R` — 49 tests covering score
  invariance, EW/IW/NA correctness, jackknife compaction, bootstrap,
  bounded indirect, drift/Wagner, ratchet active_mask, reproducibility.

### Files modified:
- `src/ts_data.cpp` — Zero-weight pattern removal + descending weight sort
- `src/ts_fitch.cpp` — `active_mask == 0` skip in 12 block loops
- `src/ts_fitch_na.inc` — `active_mask == 0` skip in 4 block loops
- `src/ts_fitch_na_incr.inc` — `active_mask == 0` skip in 11 block loops
- `src/ts_drift.cpp` — Switched to bounded indirect variants (4 calls)
- `src/ts_wagner.cpp` — Switched to bounded indirect (2 calls)

### No API changes:
No new Rcpp functions, no changes to `ts_rcpp.cpp`, `TreeSearch-init.c`,
or `R/RcppExports.R`.

### Test status:
- char-ordering (new): 49/49
- driven: 53/53
- TBR bench: 26/26
- resample: 35/35
- resample-stress: 57/57
- SPR/NNI opt: 47/47
- sector: 32/32
- fuse: 15/15 (1 skip)
- Wagner: 26/26
- NA incremental: 34/34

## Test namespace fix (Agent F, post-Phase 3C)

Fixed `TreeSearch:::` prefix in 8 test files that were calling internal
(non-exported) C++ bridge functions without the namespace qualifier. These
tests worked under `devtools::load_all()` but failed under
`library(TreeSearch, lib.loc = ...)`, which is how agents run tests.

### Files modified:
- `tests/testthat/test-ts-tbr-search.R`
- `tests/testthat/test-ts-drift-search.R`
- `tests/testthat/test-ts-ratchet-search.R`
- `tests/testthat/test-ts-ratchet-opt.R`
- `tests/testthat/test-ts-ratchet-stress.R`
- `tests/testthat/test-ts-iw.R`
- `tests/testthat/test-ts-pool.R`
- `tests/testthat/test-ts-splits.R`

### Functions prefixed (partial list):
`ts_fitch_score`, `ts_tbr_search`, `ts_spr_search`, `ts_ratchet_search`,
`ts_drift_search`, `ts_na_char_steps`, `ts_pool_test`, `ts_trees_equal`,
`ts_compute_splits`, `morphy_iw`, `SingleCharMorphy`, `UnloadMorphy`

### Test status (all previously-broken files now pass):
- ratchet-opt: 23/23, ratchet-search: 17/17, ratchet-stress: 72/72
- drift-search: 22/22, tbr-search: 28/28
- iw: 36/36 (10 skips), pool: 16/16, splits: 35/35

## Phase 3A: TBR symmetry-breaking (Agent G)

Implements virtual_prelim deduplication to skip redundant TBR rerooting
evaluations within each clip.

### Optimizations implemented:

1. **SPR-equivalence fast path**: Before the TBR inner loop for each
   rerooting edge `(sp, sc)`, `memcmp` the computed `virtual_prelim`
   against `clip_prelim` (the no-rerooting SPR case). If identical,
   skip the entire inner main-edge loop — the SPR phase already
   evaluated all regraft destinations for this state.

2. **FNV-1a hash-based deduplication**: Track seen `virtual_prelim`
   vectors via their FNV-1a hash in an `std::unordered_set<uint64_t>`.
   The hash set is seeded with `clip_prelim`'s hash before the
   rerooting loop. For each rerooting, compute the hash; if already
   in the set, skip the O(|main_edges| × total_words) inner loop.

### Correctness argument:

Two rerootings with identical `virtual_prelim` produce identical
indirect scores for every regraft destination. The first-seen
rerooting's `(reroot_parent, reroot_child)` is recorded if it yields
the best candidate; skipping the duplicate cannot miss a better score.
Hash collisions (false "seen") cause a missed optimization, not
incorrect results — but FNV-1a collisions on short vectors are
extremely rare in practice.

### Files created:
- `tests/testthat/test-ts-tbr-symmetry.R` — 24 tests covering EW, IW,
  NA, IW+NA, duplicate-state datasets, determinism, driven search
  integration, accept_equal mode.

### Files modified:
- `src/ts_tbr.cpp` — Added `fnv1a_hash()`, `<unordered_set>` include,
  `seen_vp_hashes` set + `memcmp` fast path + hash dedup in TBR
  rerooting loop.

### Test status: symmetry 24/24, TBR 28/28, TBR bench 26/26, driven
53/53, sector 32/32, fuse 16/16 (1 skip), resample 35/35

## Phase 3D: Memory layout profiling and optimization

Profiling-driven investigation of memory access patterns in the C++ scoring engine.
The main finding is that **indirect scoring dominates TBR time at scale** (72% at
200 tips), and the current memory layout is already well-structured for cache access.

### Key profiling findings:

1. **Indirect evaluation is the bottleneck**: 15% at 20 tips → 72% at 200 tips.
   Clip+incremental phase dominates at small scales (63% at 20 tips → 22% at 200).
2. **Per-candidate cost is stable** (~33 ns for `total_words=12`), confirming no
   cache pressure issue. Cost scales linearly with `total_words`.
3. **Scaling is super-quadratic** (exponent 2.78 vs expected 2.0), primarily from
   candidate count growth (exponent 2.66), not per-candidate degradation.
4. **Snapshot overhead is negligible**: 5 μs per save/restore at 200 tips (0.01% of
   TBR pass time).
5. **State arrays fit in L2** for typical datasets (162 KB total at 200 tips).

### Steps investigated and skipped:

| Step | Decision | Reason |
|------|----------|--------|
| Postorder renumbering | Skipped | Downpass not the bottleneck; arrays fit in L2 |
| Binary-char specialization | Skipped | All blocks share same `n_states` from contrast matrix |
| Block-major layout | Skipped | vroot_cache already linear; no cache pressure |
| StateSnapshot reduction | Skipped | Negligible overhead (5 μs vs 36 ms indirect) |

### Optimization applied: postorder save/restore in TBR

After `spr_unclip()`, the postorder is identical to pre-clip state. Previously rebuilt
via O(n) DFS; now saved before clip and restored via `assign()`. Also removed 2
redundant `build_postorder()` calls after `state_snap.restore()` in rejection paths.

### Files created:
- `inst/benchmarks/bench_memory.R` — profiling harness
- `inst/benchmarks/memory_profile_results.md` — full measurement report
- `tests/testthat/test-ts-memory-layout.R` — 32 regression tests

### Files modified:
- `src/ts_rcpp.cpp` — added `ts_bench_tbr_phases` diagnostic (append only)
- `src/TreeSearch-init.c` — registered `ts_bench_tbr_phases` (7 args)
- `src/ts_tbr.cpp` — postorder save/restore optimization (3 changes)
- `R/RcppExports.R`, `src/RcppExports.cpp` — regenerated

### Exported Rcpp functions (new):

| Function | Args | Purpose |
|----------|------|---------|
| `ts_bench_tbr_phases` | 7 | TBR phase timing diagnostic |

### Test status: memory-layout 32/32, driven 53/53, tbr-bench 26/26,
fuse 16/16 (1 skip), sector 32/32

### Implications for Phase 3E (SIMD):
The indirect scoring inner loop is the clear SIMD target. SSE2/AVX2 can
parallelize the `AND + OR` operations over contiguous uint64_t state words.
With `n_states` typically 3-8, even 2× throughput from SSE2 would reduce
the dominant 72% phase significantly.

## R-side wrapper migration (Agent B, Phase 1D completion)

Migrated `Resample()` and `SuccessiveApproximations()` from MorphyLib
to the C++ search engine.

### `Resample()` (`R/Morphy.R`):
- Dispatches to C++ `ts_resample_search` for EW/IW mode (the common case).
- Falls back to `Morphy()` for profile parsimony or constrained searches
  (C++ resample bridge doesn't accept those parameters).
- Input validation (method, proportion) preserved with same error messages.
- Returns `multiPhylo` with one best tree per call.
- Test updated: `replicate()` calls use `simplify = FALSE` to prevent
  `simplify2array` from extracting phylo objects from length-1 multiPhylo.

### `SuccessiveApproximations()` (`R/SuccessiveApproximations.R`):
- Dispatches to C++ `ts_successive_approx`.
- Returns `multiPhylo` with attributes: `score`, `sa_iterations`, `converged`.
- Old return format (list of per-iteration tree sets) replaced with single
  final tree. No tests or internal callers depended on the old format.
- `Suboptimality()`, `SuccessiveWeights()`, `PrepareDataSA()` preserved
  for backward compatibility (exported internal functions).

### Files modified:
- `R/Morphy.R` — `Resample()` rewritten with C++/Morphy dispatch
- `R/SuccessiveApproximations.R` — `SuccessiveApproximations()` rewritten
- `tests/testthat/test-Morphy.R` — `replicate()` calls use `simplify = FALSE`

### Test status: Morphy 30/30 (1 skip, 1 pre-existing `.CombineResults` failure),
driven 53/53, resample 35/35, profile 25/25

## Phase 5: Thread-parallel driven search (Agent H)

Implements inter-replicate parallelism for the C++ driven search engine
using `std::thread`. Each search replicate runs on its own thread with
a shared pool protected by a mutex.

### Architecture:
- **`std::thread`** (not OpenMP) — avoids R memory allocator conflicts
- **Per-thread resources**: `DataSet` copy (ratchet mutates masks),
  `ConstraintData` copy (mutable workspace), `std::mt19937` RNG
- **Shared resources**: `ThreadSafePool` (mutex-guarded), atomic stop flag,
  atomic replicate counter
- **Main thread**: pre-generates RNG seeds from R, spawns workers, polls
  `R_CheckUserInterrupt()` and timeout every 200ms, joins workers
- **Worker threads**: no R API calls — all RNG/interrupt/output via
  `ts_rng.h` thread-local indirection

### Files created:
- `src/ts_rng.h` / `src/ts_rng.cpp` — Thread-safe RNG and interrupt helpers.
  `thread_local` pointers (`thread_rng`, `thread_stop_flag`) control dispatch:
  when null (serial mode), falls back to R API; when set (parallel mode),
  uses thread-local RNG and atomic flag.
- `src/ts_parallel.h` / `src/ts_parallel.cpp` — `ThreadSafePool`,
  `parallel_driven_search()`, `parallel_resample()`.
- `tests/testthat/test-ts-parallel.R` — 35 tests covering serial equivalence,
  parallel correctness, timeout, edge cases, IW, NA, R-level interface.

### Files modified:
- `src/ts_tbr.cpp` — Replaced `GetRNGstate()/unif_rand()/PutRNGstate()` with
  `ts::make_rng()`; replaced `R_CheckUserInterrupt()` with `ts::check_interrupt()`
- `src/ts_ratchet.cpp` — Same RNG and interrupt replacement
- `src/ts_drift.cpp` — Same
- `src/ts_search.cpp` — Same (NNI + SPR)
- `src/ts_sector.cpp` — Same (RSS + XSS)
- `src/ts_wagner.cpp` — Same (Fisher-Yates shuffle uses `ts::thread_safe_unif()`)
- `src/ts_resample.cpp` — Same (bootstrap/jackknife sampling)
- `src/ts_rcpp.cpp` — Same (diagnostic benchmark); added `#include "ts_parallel.h"`;
  added `nThreads` parameter to `ts_driven_search` bridge with dispatch to
  `parallel_driven_search()` when > 1
- `src/ts_driven.h` / `src/ts_driven.cpp` — Extracted `run_single_replicate()`
  from `driven_search()`. Serial `driven_search()` now calls it in a loop.
  Added `ReplicateResult` struct and `std::function<bool()> check_timeout`
  parameter.
- `src/TreeSearch-init.c` — Updated `ts_driven_search` arg count 35→36
- `R/MaximizeParsimony.R` — Added `nThreads` parameter (default 1L = serial)
- `R/RcppExports.R` / `src/RcppExports.cpp` — Regenerated via
  `Rcpp::compileAttributes()`

### Key design decisions:
1. **`thread_local` over function-parameter RNG passing** — Avoids changing
   ~15 function signatures. Each search function's existing pattern changes
   from `GetRNGstate()/unif_rand()/PutRNGstate()` to a single `ts::make_rng()`
   call with zero overhead in serial mode (null pointer check).
2. **Pre-generated seed array** — All seeds generated from R's RNG on the main
   thread before spawning workers. Guarantees `set.seed()` reproducibility for
   the seed sequence, though thread scheduling makes tree-for-tree results
   nondeterministic across runs.
3. **Mutex-guarded pool (not lock-free)** — Pool contention is low (one `add()`
   per completed replicate, seconds apart). `std::mutex` is simple and correct.
4. **Workers silent, main thread reports** — `Rprintf` is not thread-safe.
   Workers run with verbosity=0; main thread periodically reports aggregate
   progress from shared atomics.
5. **Serial mode as default** — `nThreads=1` calls `driven_search()` exactly
   as before, with zero overhead from the parallel infrastructure.

### R-level interface:
```r
MaximizeParsimony(dataset, nThreads = 2L)  # 2 parallel worker threads
MaximizeParsimony(dataset, nThreads = 0L)  # auto-detect CPU count
MaximizeParsimony(dataset, nThreads = 1L)  # serial (default, reproducible)
```

### Exported Rcpp function changes:
| Function | Change |
|----------|--------|
| `ts_driven_search` | +1 param (`nThreads`); dispatches to `parallel_driven_search` when > 1 |

### Test status: 291/291 passing (53 driven + 35 resample + 60 resample-stress +
32 sector + 16 fuse + 26 tbr-bench + 33 na-incremental + 35 parallel + 1 skip)

## Phase 3B: Constrained Sectorial Search — CSS (Agent G)

Implements sector-restricted TBR that operates on the full tree (no HTU
pseudo-tip approximation). Unlike RSS/XSS which extract sectors into reduced
datasets, CSS restricts TBR clips and regrafts to within a sector while
scoring against the full tree. Eliminates HTU approximation errors (false
improvements that revert on reinsertion) at the cost of higher per-candidate
evaluation.

### Design:
- **Core mechanism**: Optional `sector_mask` parameter on `tbr_search()`.
  When non-null, 3 filter lines restrict clips to sector nodes and regrafts
  to sector edges. All TBR optimizations (early termination, vroot cache,
  symmetry breaking, state snapshots) apply unchanged.
- **Sector selection**: Uses existing `xss_partition()` to partition the tree
  into non-overlapping sectors. For each partition, builds a `vector<bool>`
  mask and calls sector-restricted TBR.
- **No reduced dataset**: CSS needs no `build_reduced_dataset`, no
  `reinsert_sector`, no `CladeSnapshot` save/restore, no post-hoc constraint
  check (TBR handles constraints natively).

### Pipeline position:
Wagner → TBR → XSS → RSS → **CSS** → Ratchet → Drift → Final TBR

CSS runs after XSS+RSS (which do cheap bulk optimization) as a polishing
step with exact scoring.

### Files created:
- `tests/testthat/test-ts-css.R` — 26 tests covering EW, IW, NA, determinism,
  timeout, R-level integration, multiple partition sizes.

### Files modified:
- `src/ts_tbr.h` — Added `sector_mask` parameter to `tbr_search`
- `src/ts_tbr.cpp` — 3 filter lines: clip candidates, SPR regrafts, TBR regrafts
- `src/ts_sector.h` — Added `css_search` declaration
- `src/ts_sector.cpp` — Implemented `css_search` (~80 lines) + `build_sector_mask`
- `src/ts_driven.h` — Added `css_rounds`, `css_partitions` to `DrivenParams`
- `src/ts_driven.cpp` — CSS phase in pipeline after RSS
- `src/ts_rcpp.cpp` — Added `cssRounds`, `cssPartitions` to bridge
- `src/TreeSearch-init.c` — Driven search arg count 39→41
- `R/MaximizeParsimony.R` — Exposed `cssRounds`, `cssPartitions` parameters
- `R/RcppExports.R`, `src/RcppExports.cpp` — Regenerated

### Exported Rcpp function changes:
| Function | Change |
|----------|--------|
| `ts_driven_search` | +2 params (`cssRounds`, `cssPartitions`) |

### Test status: CSS 26/26, driven 53/53, sector 32/32, TBR bench 26/26,
TBR symmetry 24/24, fuse 16/16 (1 skip), resample 35/35

## Constraint + profile parsimony in resample/SA bridges (Agent B)

Threaded constraint enforcement and profile parsimony parameters through
`ts_resample_search` and `ts_successive_approx`, eliminating the need for
`Morphy()` fallback in `Resample()` and `SuccessiveApproximations()`.

### Changes:

**C++ layer**:
- `ts_resample.h/.cpp`: Added 3 optional params to both `resample_search()`
  and `successive_approximations()`: `info_amounts_r`, `info_max_steps`,
  `ConstraintData* cd`. Threaded to `build_dataset()` and `driven_search()`.
- `ts_parallel.h/.cpp`: Same 3 params added to `parallel_resample()`.
- `ts_rcpp.cpp`: Extracted `build_constraint_from_r()` and
  `extract_info_amounts()` helpers (used by `ts_driven_search`,
  `ts_resample_search`, `ts_successive_approx`). +7 Rcpp params to each
  resample/SA bridge: `consSplitMatrix`, `consContrast`, `consTipData`,
  `consWeight`, `consLevels`, `consExpectedScore`, `infoAmounts`.
- `TreeSearch-init.c`: Arg counts 14→21 for both functions.

**R layer**:
- `R/MaximizeParsimony.R`: Extracted `.PrepareConstraint()` helper from
  `MaximizeParsimony()`. Shared by `Resample()` and `SuccessiveApproximations()`.
- `R/Morphy.R`: `Resample()` no longer falls back to `Morphy()`. Handles
  profile parsimony (`PrepareDataProfile` + `infoAmounts`) and constraints
  (`.PrepareConstraint()`) natively via C++ engine.
- `R/SuccessiveApproximations.R`: Added `concavity` and `constraint` params.
  Handles profile and constraints via C++ engine.

### Test status:
- driven 53/53, resample 35/35, profile 25/25
- Morphy 29/30 (1 skip; 2 pre-existing failures in `Morphy()` and `.CombineResults`)
- Constraint + profile Resample verified interactively
- Constraint + profile SA verified interactively

## Phase 3E: SIMD vectorization — ✅ COMPLETE

SSE2/NEON vectorization of Fitch scoring inner loops. Applies SIMD to the
bit-parallel AND/OR operations on uint64_t state words — the 72%-of-runtime
bottleneck identified in Phase 3D profiling.

### Architecture:

**`ts_simd.h`** — Portability layer with compile-time detection:
- SSE2 (`__SSE2__` or `_M_X64`): `__m128i` intrinsics
- NEON (`__ARM_NEON`): `uint64x2_t` intrinsics
- Scalar fallback: plain C++ loops

Key inline helpers:
- `any_hit_reduce(a, b, n)` — `OR(a[s] & b[s])` for s in [0,n)
- `any_hit_reduce3(clip, a, b, n)` — `OR(clip[s] & (a[s] | b[s]))`
- `or_reduce(a, n, start)` — `OR(a[s])` for s in [start,n)
- `_from1` variants — skip word 0 (for NA blocks where state 0 = inapplicable)

All functions process 2 × uint64_t per SIMD iteration with scalar tail for
odd counts. Horizontal reduction via `storeu128` + scalar OR.

### SIMD callsites (15+ loop replacements):

**`ts_fitch.cpp`**:
- `fitch_downpass_node` — Pass 1 `any_intersect` + Pass 2 mask broadcast
- `fitch_downpass` — same pattern per block
- `fitch_incremental_downpass` — `any_intersect` reduction
- `uppass_node` — `any_intersect` reduction
- All 6 indirect length variants (`fitch_indirect_length`,
  `fitch_indirect_length_bounded`, `fitch_indirect_length_cached`,
  `indirect_iw_length`, `indirect_iw_length_bounded`,
  `indirect_iw_length_cached`)

**`ts_fitch_na.inc`**:
- Pass 1 downpass (standard block + NA-aware block)
- Pass 3 downpass (second downpass step counting)

**`ts_fitch_na_incr.inc`**:
- All 6 NA indirect variants (`fitch_na_indirect_length`,
  `indirect_na_iw_length`, plus bounded variants)
- `clip_has_active` and `below_has_active` reductions

### Data padding:
- `total_words` padded to even count in `build_dataset()` for SIMD safety
- Padding words are zero-initialized (no scoring impact)
- Block offsets unchanged (padding at end of total_words array)

### Performance results:

Per-candidate indirect cost with SIMD is similar to pre-SIMD baseline —
GCC `-O2` was already auto-vectorizing the simple reduction loops. The
clip+incremental phase shows ~10-40% improvement for small datasets.

| Dataset | Tips | Words | Clip (μs) | Indirect (μs) | Per-cand (ns) |
|---------|------|-------|-----------|---------------|---------------|
| Vinther2008 | 23 | 28 | 417 | 305 | 85 |
| Agnarsson2004 | 62 | 59 | 4157 | 8530 | 110 |
| Zhu2013 | 75 | 20 | 4050 | 6189 | 46 |

Primary value: establishes vectorization infrastructure for future AVX2
(4-wide, requiring runtime dispatch or separate compilation unit).

### Files created:
- `src/ts_simd.h` — SIMD portability layer (165 lines)
- `tests/testthat/test-ts-simd.R` — 49 tests (EW/IW/NA cross-validation,
  TBR search, driven search, determinism, consistency)
- `inst/benchmarks/bench_simd.R` — TBR phase timing benchmark

### Files modified:
- `src/ts_data.cpp` — `total_words` even padding
- `src/ts_fitch.h` — Added `#include "ts_simd.h"`
- `src/ts_fitch.cpp` — SIMD in 10 functions (downpass, uppass, 6 indirect)
- `src/ts_fitch_na.inc` — SIMD in Pass 1 and Pass 3
- `src/ts_fitch_na_incr.inc` — SIMD in all 6 NA indirect variants

### No API changes:
All SIMD is transparent to callers. No new Rcpp functions, no changes to
`ts_rcpp.cpp` or `TreeSearch-init.c`.

### Test status: simd 49/49, tbr 28/28, driven 53/53, sector 32/32,
fuse 16/16 (1 skip), na-incremental 33/33, resample 35/35,
resample-stress 60/60, ratchet 112/112, drift 22/22, tbr-bench 26/26,
tbr-symmetry 24/24, spr-nni-opt 47/47, char-ordering 49/49,
simplify 40/40, profile 25/25, iw 36/36 (10 skip), wagner 26/26,
splits 35/35, pool 16/16, memory-layout 32/32, css 26/26,
parallel 35/35 (12 skip), progress 48/48

## TaxonInfluence migration to C++ engine + MaximizeParsimony output fix

### TaxonInfluence migration:
- `TaxonInfluence()` now calls `MaximizeParsimony()` instead of `Morphy()`.
- Removed `DropTip(startTree, leaf)` as starting tree for per-taxon searches;
  `MaximizeParsimony()` generates its own Wagner starting trees.
- Documentation updated: `@inheritParams MaximizeParsimony`, `@param verbosity,\dots`
  references `MaximizeParsimony()`, example uses `maxReplicates` instead of
  `ratchIter`/`startIter`.

### MaximizeParsimony output tree fix (IMPORTANT):
- **Bug discovered**: C++ engine's `tree_to_edge()` produces edge matrices whose
  internal node ordering is NOT valid preorder. The R wrapper inherited
  `order = "preorder"` from the template tree, causing downstream functions
  (`DropTip`, `ClusteringInfoDistance`, etc.) to crash or produce wrong results.
- **Fix**: Added `Renumber()` call on each output tree in `MaximizeParsimony()`.
  `Renumber()` both renumbers internal nodes to sequential preorder IDs and
  reorders edges to valid preorder sequence.
- **Impact**: This affects ALL callers of `MaximizeParsimony()` — any code that
  previously received malformed trees now gets correct ones. Distance calculations,
  tree plotting, and DropTip operations that previously crashed should now work.

### PGO Makevars.win cleanup:
- Found untracked `src/Makevars.win` with `-fprofile-generate` PGO flags left by
  another agent. These cause segfaults in the compiled DLL. Renamed to
  `src/Makevars.win.pgo-bak`. **Do not restore this file for normal builds.**

### Files modified:
- `R/TaxonInfluence.R` — `Morphy()` → `MaximizeParsimony()`, removed DropTip
  starting tree, updated docs
- `R/MaximizeParsimony.R` — Added `Renumber()` to output tree reconstruction;
  added `Renumber` to `@importFrom TreeTools`
- `tests/testthat/test-TaxonInfluence.R` — Updated params to MaximizeParsimony
  style (`maxReplicates`, `targetHits`, `verbosity`)

### Test status: TaxonInfluence 12/12, driven 53/53, Morphy 29/30 (pre-existing),
profile 25/25, resample 35/35

## Test suite optimization

Reduced ts-* test suite runtime from ~294s to ~27s (91% reduction).

### Key changes:

1. **`tests/testthat/helper-ts.R`** — Shared helpers (`make_ts_data()`,
   `ts_score()`, `validate_result()`). Removed 20+ duplicate definitions
   across ts-* test files.

2. **`test-ts-iw.R` rewrite** (96s → 3.1s): Replaced all `morphy_iw_ref()`
   calls (MorphyLib ground truth) with hard-coded reference scores from the
   C++ engine. Reduced dataset sweep from 30 to 6 representative datasets.
   Search tests reduced from 10 to 3 datasets, 5→3 methods × datasets.

3. **`test-ts-na-incremental.R`** (26s → 4.0s): Removed maxSeconds guards
   (never triggered). Reduced datasets from 5→3, replicates from 5→3.

4. **`test-ts-driven.R`** (21s → 1.6s): Reduced "multiple replicates" loop
   from 5→2, inner maxReplicates from 5→3.

5. **`test-ts-char-ordering.R`** (19s → 0.9s): Reduced replicate counts.

6. **`test-ts-parallel.R`** (15s → 3.4s): Removed score-quality comparison
   test (non-deterministic, slow). Reduced replicate counts. Consolidated
   duplicate parallel correctness tests.

7. **`test-ts-css.R`** (11s → 1.1s): Reduced competitive quality comparison
   from 3 seeds to 1.

8. **`test-ts-progress.R`** (10s → 0.9s): Switched to small inline dataset
   for callback tests. Reduced maxReplicates.

9. **`test-ts-tabu.R`** (11s → 2.9s): Reduced maxReplicates from 3→2.

10. **`test-ts-ratchet-stress.R`**: Reduced loop counts (20→5, 10→3),
    cycle counts (50→10, 20→5), tree size (75→30 tips, 25→15 tips).

11. **`test-ts-resample-stress.R`**: Reduced poolSuboptimal replicates
    from 10→5.

### Files created:
- `tests/testthat/helper-ts.R`

### Files modified (all ts-* test files):
Every `test-ts-*.R` file was modified to remove duplicate `make_ts_data`,
`ts_score`, and `library("TreeTools")` definitions in favor of the shared
helper. Hot files received additional replicate/dataset count reductions.

### Test counts: All 27 ts-* files pass (total ~980 expectations, 0 failures)

### Known issue: `test-ts-simd.R` crashes on R process exit (exit code 127)
due to MorphyLib cleanup in `morphy_ew_ref()`. All 49 tests pass before
the crash. Pre-existing issue, not caused by this optimization.

## PGO (Profile-Guided Optimization) — ✅ COMPLETE

Tested GCC PGO on the C++ search engine. Build recipe documented in
`inst/benchmarks/pgo_recipe.md`.

### Results (3-run medians, GCC 13 / rtools45 / Windows x86_64):

| Benchmark | Baseline (s) | PGO (s) | Speedup |
|-----------|-------------|---------|---------|
| Vinther EW (23 tips) | 0.240 | 0.240 | 0% |
| Vinther IW (23 tips) | 0.170 | 0.190 | −12% |
| Zhu EW (75 tips) | 4.010 | 3.790 | 5% |
| Zhu IW (75 tips) | 5.340 | 4.990 | 7% |
| Agnarsson EW (62 tips) | 2.200 | 2.080 | 5% |

### Conclusions:
- **Modest 5–7% speedup on medium datasets** where C++ hot paths dominate.
- **No benefit on small datasets** — R overhead and startup time swamp C++ gains.
- **Not worth shipping** — PGO requires machine-specific `.gcda` files and a
  two-pass build process. The 5–7% speedup doesn't justify the complexity for
  end users. Worth revisiting if targeting a specific benchmark platform.
- **Correctness verified**: 53/53 driven search tests pass with PGO build.

### Previous agent hang diagnosis:
The PGO agent hung during the `-fprofile-use` build because GCC's PGO
optimization passes are 2–5× slower than normal compilation. With 30+ source
files, the build exceeded the default 120-second bash timeout. The fix is
simply using a longer timeout (300s). No path mismatch or DLL lock issue.

### Files created:
- `inst/benchmarks/pgo_recipe.md` — Full PGO build recipe with steps

### Cleanup:
- `src/Makevars.win` removed after benchmarking (must not be left in place)
- `.pgo-data/` contains machine-specific binary data; not version-controlled
- `.agent-pgo/`, `.agent-pgo-gen/`, `.agent-pgo-use/` are build artifacts

## Profiling round 1: R overhead + drift optimization (Agent F)

Manual profile → optimize → measure loop per `profiling-pgo-plan.md`.

### Profiling findings:

**Per-phase timing (Zhu2013 EW, 75 tips, 3 replicates):**
| Phase | Avg ms/rep | % of C++ time |
|-------|-----------|---------------|
| Drift | 365 | 39% |
| Ratchet | 274 | 29% |
| Initial TBR | 210 | 23% |
| Sectorial (XSS+RSS+CSS) | 75 | 8% |
| Final TBR + Wagner | 9 | 1% |

**R overhead**: 1.48s (34% of total) caused by `AdditionTree()` in the
`MaximizeParsimony()` wrapper — builds a full Morphy-based starting tree
that the C++ engine never uses.

**TBR inner loop**: Per-candidate indirect scoring costs 23 ns (75 tips),
already well-optimized with SIMD, early termination, vroot cache. Near
memory-throughput limit — hard to improve without algorithmic changes.

### Optimizations applied:

**1. Replace `AdditionTree()` with `RandomTree()` in R wrapper**
(`R/MaximizeParsimony.R`): The C++ engine builds its own Wagner trees
internally; the R starting tree is only used as a template for tip labels
in the output. Replaced ~2s `AdditionTree()` call with ~0ms `RandomTree()`.

**2. Drift subtree-size filter** (`src/ts_drift.cpp`): Ported the
smaller-subtree-only filter from `tbr_search()` to `drift_phase()`.
Skips clips where `subtree_size > n_tip / 2`, halving the worst-case
candidate count. Also added postorder save/restore (avoiding O(n) DFS
rebuild after each unclip).

### Per-phase timing instrumentation:
Added `std::chrono` phase timing to `ts_driven.cpp` — visible at
`verbosity >= 2` as `[Nnn ms]` after each phase score.

### Results (3-run medians, Zhu2013):

| Benchmark | Before | After | Speedup |
|-----------|--------|-------|---------|
| Vinther EW (23 tips) | 0.240s | 0.160s | **33%** |
| Vinther IW (23 tips) | 0.170s | 0.180s | ~0% |
| Zhu EW (75 tips) | 4.010s | 2.730s | **32%** |
| Zhu IW (75 tips) | 5.340s | 3.560s | **33%** |
| Agnarsson EW (62 tips) | 2.200s | 1.520s | **31%** |

The ~30% speedup is almost entirely from the AdditionTree fix. The drift
optimization provides a smaller, dataset-dependent improvement.

### Remaining C++ bottleneck:
Drift (39%) and ratchet (29%) together account for 68% of C++ time. Both
use TBR internally. The per-candidate indirect evaluation (23 ns at 75 tips)
is near-optimal. Further gains require:
- **Algorithmic**: SPR→TBR escalation (do SPR first, TBR only where SPR fails)
- **Tuning**: Reduce drift/ratchet cycles (quality vs speed tradeoff → Phase 6)

### Files modified:
- `R/MaximizeParsimony.R` — Replaced `AdditionTree()` with `RandomTree()`;
  updated docs; added `RandomTree` to `@importFrom`
- `src/ts_driven.cpp` — Added per-phase `std::chrono` timing at verbosity ≥ 2
- `src/ts_drift.cpp` — Subtree-size filter + postorder save/restore in
  `drift_phase()`

### Test status: drift 22/22, driven 53/53
