
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

## Shared files — coordination rules

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` are modified by every
agent (to add Rcpp bridges and register C entry points). To avoid
merge conflicts:

- **Append only** — add your new declarations / registrations at the
  end of the relevant section.
- Agents should not reformat or reorder existing entries.

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

## Upcoming agents

| Agent | Step | New files | Notes |
|-------|------|-----------|-------|
| G | Driven search (Step 8) | `ts_driven.h/.cpp`, `test-ts-driven.R` | Orchestrates everything: Wagner start → TBR → ratchet/drift/sector/fuse loop |
