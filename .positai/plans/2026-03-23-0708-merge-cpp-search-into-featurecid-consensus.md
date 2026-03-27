# Plan: Merge cpp-search into feature/cid-consensus

## Context

`feature/cid-consensus` is 6 commits ahead and 16 commits behind `cpp-search`.
The branch has uncommitted work (now stashed as `stash@{0}`) adding ~400 lines
of CID wiring (Rcpp bridge, drift two-baseline, sector CID, ratchet CID).

Trial merge shows **1 conflict** (AGENTS.md only). All source files auto-merge
cleanly. The stash pop will likely conflict on `ts_drift.cpp` and
`ts_sector.cpp` (both modified by cpp-search's collapsed-edge features and by
the stash's CID wiring).

## Steps

### 1. Merge cpp-search

```bash
git merge cpp-search
```

### 2. Resolve AGENTS.md (only conflict)

Two conflict regions:

- **C++ module map** (~line 489): CID branch adds a CID row, cpp-search adds
  a Collapsed row. **Keep both rows.**

- **Benchmarks & profiling** (~line 871): CID branch has CID scoring
  optimization notes; cpp-search has updated profiling phase table + ratchet
  tuning validation. **Keep both sections** (CID optimizations section
  followed by ratchet tuning section, or vice versa).

Commit the merge.

### 3. Pop stash

```bash
git stash pop
```

Expected conflicts: `ts_drift.cpp`, `ts_sector.cpp` (cpp-search added
collapsed-edge skipping; stash adds CID two-baseline tracking). These are
independent concerns operating on different code paths — resolution is to
keep both.

Other stashed files (`ts_rcpp.cpp`, `ts_ratchet.cpp`, `ts_fitch.cpp`,
`ts_data.h`, `RcppExports.*`, `TreeSearch-init.c`, NAMESPACE, man pages)
should pop cleanly since cpp-search didn't touch them.

### 4. Resolve stash conflicts (if any)

For `ts_drift.cpp`: cpp-search added `#include "ts_collapsed.h"` and
collapsed-edge logic; stash added `#include "ts_cid.h"` and CID
two-baseline tracking. Keep both — they're additive.

For `ts_sector.cpp`: cpp-search rewrote sector internals (conflict-guided
RSS, collapsed dedup); stash added CID-mode dispatch. Integrate carefully.

### 5. Build and test

```bash
rm -f src/*.o src/*.dll
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-P "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
```

Then run tests:
```bash
Rscript -e "library(TreeSearch, lib.loc='.agent-P'); testthat::test_dir('tests/testthat')"
```

### 6. Commit the stash work

Stage and commit the popped stash changes (CID wiring) on top of the merge.

## Risk assessment

- **Low risk**: AGENTS.md resolution is purely editorial.
- **Medium risk**: stash pop conflicts on ts_drift.cpp and ts_sector.cpp
  require careful integration of two independent features (collapsed edges +
  CID scoring). Both are additive — no logical conflicts, just textual ones.
- **Build risk**: The merged code introduces new `#include` dependencies
  (ts_collapsed.h) into files that also get CID includes. Compilation will
  surface any header issues immediately.
