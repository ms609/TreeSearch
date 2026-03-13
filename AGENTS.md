
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

## Upcoming agents

| Agent | Step | New files | Notes |
|-------|------|-----------|-------|
| E | Tree fusing (Step 5) | `ts_fuse.h/.cpp`, `test-ts-fuse.R` | Depends on C's splits/pool |
| F | Sectorial search (Step 7) | `ts_sector.h/.cpp`, `test-ts-sector.R` | Most complex; recursive search |
| G | Driven search (Step 8) | `ts_driven.h/.cpp`, `test-ts-driven.R` | Orchestrates everything |
