# Plan: Enable Sectorial Search for CID Consensus

**Date:** 2026-03-20
**Context:** CIDConsensus currently disables all sectorial search (XSS, RSS, CSS
rounds hardcoded to 0). This plan re-enables them.

---

## Analysis: Why sectorial search was disabled

Two bugs prevent sectorial search from working in CID mode:

### Bug 1: Reduced dataset inherits CID scoring mode

`build_reduced_dataset()` (line 270 of `ts_sector.cpp`) copies `scoring_mode`
from the full dataset:
```cpp
rd.data.scoring_mode = ds.scoring_mode;  // copies CID
```

When `score_tree(rd.subtree, rd.data)` is later called on the sector (line 528),
it dispatches to `cid_score(tree, *ds.cid_data)` — but the reduced dataset's
`cid_data` is null (never set). This would crash.

### Bug 2: `final_` states not populated in CID mode

`build_reduced_dataset()` reads `tree.final_[...]` to construct the HTU
pseudo-tip's state (line 298). But in CID mode, `score_tree()` dispatches to
`cid_score()` which computes splits from topology — it never runs Fitch, so
`final_` is never populated (or is zeroed by `reset_states()`). The HTU would
get all-zero states, making sector search useless.

---

## Approach: MRP-Fitch proxy + CID verification

The CID DataSet already contains MRP (Matrix Representation with Parsimony)
binary characters for TBR screening. These encode each input tree split as a
binary Fitch character. Sectorial search can use these same MRP characters for
within-sector optimization, with full CID verification after reinsertion.

This extends the dual-layer pattern already used by TBR:
- **Fast layer**: Fitch on MRP characters (additive, supports HTU decomposition)
- **Accurate layer**: Full CID scoring on the complete tree

### Why MRP is a good proxy for CID in sectors

1. Each MRP character represents one split from one input tree. Fitch
   parsimony on MRP measures compatibility with input tree structure — closely
   related to CID.
2. The HTU `final_` states from Fitch uppass encode the optimal MRP
   assignment at the sector boundary — the principled decomposition that
   sectorial search was designed for.
3. Full-tree CID verification after reinsertion catches any proxy errors.
   We never accept a move that worsens actual CID.

### CSS needs zero changes

CSS uses sector_mask TBR on the **full** tree (no reduced dataset, no HTU).
TBR already works in CID mode (MRP incremental screening + CID verification).
The sector_mask just restricts clip/regraft sites. Enabling CSS is purely an
R-wrapper change.

---

## Changes

### 1. `src/ts_fitch.cpp` — maintain Fitch state arrays in CID mode

In `score_tree()`, run `fitch_score()` before `cid_score()` so that
`prelim`/`final_` are always valid after a CID-mode `score_tree()` call:

```cpp
double score_tree(TreeState& tree, const DataSet& ds) {
  if (ds.scoring_mode == ScoringMode::CID) {
    // Run Fitch on MRP characters to maintain prelim/final_ state arrays.
    // Needed by build_reduced_dataset() for HTU construction, and by
    // TBR incremental scoring as a valid baseline.
    fitch_score(tree, ds);
    return cid_score(tree, *ds.cid_data);
  }
  // ... existing dispatches unchanged
}
```

**Cost**: One extra Fitch pass per `score_tree()` call in CID mode.
For typical sizes (50 tips, ~5000 MRP chars), this is O(250K) operations
(microseconds). CID scoring itself is O(n_trees × n_splits² × LAP) —
milliseconds. Overhead < 1%.

**Correctness**: `fitch_score()` calls `fitch_downpass()` (traverses
`tree.postorder`) then `fitch_uppass()`. All existing callers ensure
postorder is valid before `score_tree()` (TBR calls `build_postorder()`
at line 843 before `full_rescore()`; sector code calls `build_postorder()`
at lines 312, 545, 663).

### 2. `src/ts_sector.cpp` — override scoring mode for CID reduced datasets

After the existing scoring_mode copy (line 270), add:

```cpp
// CID mode: sector-internal search uses Fitch on MRP characters.
// Full CID is verified on the complete tree after reinsertion.
if (rd.data.scoring_mode == ScoringMode::CID) {
  rd.data.scoring_mode = ScoringMode::EW;
  rd.data.cid_data = nullptr;
}
```

This makes `score_tree(rd.subtree, rd.data)` use Fitch parsimony within
sectors, which is correct for MRP characters and supports the HTU
decomposition.

### 3. `R/CIDConsensus.R` — enable sector parameters

Replace hardcoded disabled sectors with SearchControl-driven defaults.
In `.CIDDrivenSearch()`:

```r
xssRounds = .NullOr(ctrl[["xssRounds"]], 3L),     # was 0L
xssPartitions = .NullOr(ctrl[["xssPartitions"]], 4L),
rssRounds = .NullOr(ctrl[["rssRounds"]], 3L),      # was 0L
cssRounds = .NullOr(ctrl[["cssRounds"]], 2L),       # was 0L
cssPartitions = .NullOr(ctrl[["cssPartitions"]], 4L),
```

Also in `.TopologySearch()` (the light re-optimization used by rogue
dropping), enable sectors when tree is large enough:

```r
xssRounds = if (nTip >= 12L) 1L else 0L,
rssRounds = 0L,   # keep light for re-optimization
cssRounds = 0L,
```

### 4. Tests (in `tests/testthat/test-ts-cid.R` or new file)

**Tier 2** (skip on CRAN).

1. **Fitch state validity**: After `ts_cid_consensus()` on a small case,
   verify the returned tree's MRP Fitch score is consistent (call
   `ts_fitch_score()` on the MRP data and confirm it's a reasonable
   integer).

2. **Sector smoke test**: Run `ts_cid_consensus()` with
   `xssRounds = 2, rssRounds = 2, cssRounds = 1` on a 20-tip, 30-tree
   case. Verify it returns a valid tree and score ≤ the score from
   sectors-disabled run.

3. **Regression**: Compare CIDConsensus scores with sectors enabled vs
   disabled on the existing test cases. Sectors-enabled should be ≤.

4. **Small tree guard**: Verify that on a 6-tip case (below
   `2 * sectorMinSize`), sectors are skipped gracefully.

### 5. Update comments

Remove the "Sectorial search disabled" comments in:
- `R/CIDConsensus.R` (lines 262-263)
- Update `AGENTS.md` architecture notes

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| MRP proxy doesn't correlate well with CID for sector-level moves | Full CID verification after reinsertion ensures correctness. If sectors rarely help, they waste time but don't harm quality. Empirical benchmarking will determine optimal defaults. |
| Extra Fitch in `score_tree()` adds overhead during TBR | Overhead is < 1% of CID scoring. Measurable with `ts_bench_tbr_phases`. |
| `fitch_score()` in `score_tree()` assumes valid postorder | All existing callers already ensure this (verified by code audit). Add a debug assertion if desired. |
| Sector search produces poor results and wastes wall-clock time | Default to modest sector rounds (3 XSS, 3 RSS, 2 CSS) rather than parsimony defaults. Can be tuned via SearchControl. |

---

## Implementation order

1. Fix `score_tree()` in `ts_fitch.cpp` (change #1)
2. Fix `build_reduced_dataset()` in `ts_sector.cpp` (change #2)
3. Build and run existing tests (verify no regression)
4. Update R wrapper defaults (change #3)
5. Add new tests (change #4)
6. Benchmark: compare CID scores with sectors enabled vs disabled
7. Update comments and AGENTS.md

---

## Results

### Implementation notes

Change #1 was revised during implementation.  Running `fitch_score()` inside
every `score_tree()` call triggers a pre-existing memory issue (heap
fragmentation from accumulated CID evaluations via LAP + split temporaries)
sooner.  Instead, Fitch state preparation (`prepare_cid_states()`) is called
only in `rss_search()` and `xss_search()` before sector extraction and after
each reinsertion.  This avoids extra overhead in non-sector code paths.

Sector defaults are gated on `nTip >= 20` because CID scoring is expensive
per evaluation and small trees converge quickly with TBR alone.

### Benchmark results (30-tip, 50 input trees, 5 replicates)

| Seed | No sectors | With sectors | Time (ns) | Time (s) |
|------|-----------|-------------|-----------|----------|
| 1001 | 16.29 | 2.99 | 0.61 | 0.42 |
| 2002 | 18.14 | 3.08 | 0.64 | 0.50 |
| 3003 | 16.04 | 3.05 | 0.61 | 0.40 |
| 4004 | 16.01 | 3.91 | 0.68 | 0.50 |
| 5005 | 23.30 | 3.07 | 0.73 | 0.45 |

**~5× CID score improvement** with sectors enabled.  Convergence is also
faster (0.40–0.50s vs 0.61–0.73s) because XSS finds good sector topologies
early, breaking TBR local optima.

### Test status

- 274 C++ engine tests pass (28 CID + 42 sector + 152 driven + 52 TBR)
- CIDConsensus R-level tests: 121/125 pass; last 4 trigger pre-existing
  `std::bad_alloc` from accumulated CID searches (see known issue in AGENTS.md)
