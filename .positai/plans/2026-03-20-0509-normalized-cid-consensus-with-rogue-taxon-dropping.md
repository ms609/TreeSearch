# Plan: Normalized CID Consensus with Rogue Taxon Dropping

**Date:** 2026-03-20
**Status:** Empirical investigation complete; implementation in progress

---

## Resolved Design Decisions

### 1. Scoring formulation

Add a `normalize` parameter (default `FALSE`):

| Mode | Objective | Score formula (minimized) | Notes |
|------|-----------|---------------------------|-------|
| `normalize = FALSE` | Maximize total shared info | `-mean(MCI_i)` | Default; maximizes bits of input information captured |
| `normalize = TRUE` | Maximize fraction captured | `1 - mean(MCI_i / CE_i)` | Required for meaningful rogue dropping |

**Fast path**: both modes use precomputed input splits and CEs. The
per-candidate loop computes `MCI_i` via `MutualClusteringInfoSplits()`
and either sums directly or divides by `CE_i`. Overhead: one extra
division per input tree when normalized.

**Rationale** (from empirical testing on 50-tip/100-tree primates data):
- Raw CID (CE_cand + CE_i - 2*MCI) is NOT suitable as a scoring
  objective for rogue dropping: 41/50 tips "improve" when dropped,
  because fewer tips = less possible distance.
- `-mean(MCI)` maximizes information content directly. Dropping a tip
  always loses some MCI, so rogue dropping is dubious without
  normalization.
- `1 - mean(MCI/CE_i)` correctly identifies ~29/50 tips as beneficial
  to drop. Rogue taxon 'r' in the toy example goes from 0.30 to 0.00.
- `-mean(MCI)` and `1 - mean(MCI/CE_i)` have rank correlation rho=0.06
  on random candidate topologies for a fixed tip set, so they genuinely
  optimize for different things. The user should choose.

### 2. neverDrop semantics

| Value | Meaning |
|-------|---------|
| `TRUE` | No dropping; current behaviour |
| `FALSE` (default) | All taxa droppable |
| `character(...)` / `integer(...)` | Listed tips protected |

If `neverDrop != TRUE` and `normalize = FALSE`, emit a warning that raw
MCI rogue dropping may not be meaningful.

### 3. Pre-screening

Tested two approaches on 50-tip / 100-tree primates data:

| Method | Time | Spearman rho with marginal NID |
|--------|------|-------------------------------|
| **Marginal NID** (drop each tip, rescore) | 9.3 s | 1.0 (by definition) |
| **TipInstability** (Rogue package) | 0.04 s | 0.36 |

TipInstability is a poor proxy (rho=0.36). Use marginal NID as the
primary pre-screen. For very large datasets (>200 tips), consider
TipInstability to narrow to top-k candidates, then marginal NID.

### 4. Architecture: outer loop (not inner-loop move proposal)

Rogue dropping is Phase 3, running after topology search + collapse/resolve.
Each accepted drop triggers a lighter re-optimization (half ratchet iters).
Trees with <8 tips skip Ratchet and use only collapse/resolve.

This avoids EdgeSwapper interface changes and is sufficient: each drop
changes the landscape enough to warrant re-optimizing topology.

---

## Empirical Results (primates 50-tip / 100-tree)

### Greedy sequential dropping (normalized, from MR consensus)

| Step | Taxon | Score (1-MCI/CE) | Tips |
|------|-------|-----------------|------|
| 0 | (start) | 0.4356 | 50 |
| 1 | Cercopithecus_mitis_stuhlmanni | 0.4193 | 49 |
| 2 | Cercopithecus_doggetti | 0.4049 | 48 |
| 3 | Cercopithecus_erythrotis_camerunensis | 0.3905 | 47 |
| 4 | Cercopithecus_roloway | 0.3783 | 46 |
| 5 | Cercopithecus_albogularis_francescae | 0.3666 | 45 |

All top rogues are Cercopithecus species (biologically plausible: closely
related, branching order uncertain across bootstrap replicates).

### Formulation agreement

All three normalization schemes (raw CID, NVI, MCI/CE) agree on rogue
ranking with high Spearman correlation (0.91-0.97 between formulations).

### CIDConsensus timing

17s for 3 ratchet iterations on 50-tip x 100-tree data (raw CID mode).
Score improved from 20.60 (MR) to 11.57 (after collapse/resolve).

---

## Implementation

### Modified: `R/CIDConsensus.R`

**New parameters:**
- `normalize = FALSE` in `CIDConsensus()`
- `neverDrop = FALSE` in `CIDConsensus()`

**New internal functions:**
- `.TopologySearch()` -- factored out of `CIDConsensus()` for reuse
- `.RogueRefine()` -- Phase 3 outer loop
- `.PrescreenMarginalNID()` -- marginal NID pre-screen
- `.ScoreTree()` -- score a phylo against cidData
- `.PruneTrees()` -- drop tips from multiPhylo
- `.BestInsertion()` -- find best edge to graft a restored tip
- `.InsertTipAtEdge()` -- graft a tip onto a specific edge

**Modified internal functions:**
- `.MakeCIDData()` -- gains `normalize` parameter
- `.CIDScorer()` -- supports both `-mean(MCI)` and `1-mean(MCI/CE_i)`
- `.CIDScoreFast()` -- same
- `.CIDBootstrap()` -- simplified (no `meanInputCE` to track)

### Safety guards

- Trees with <8 tips skip Ratchet re-optimization (only collapse/resolve)
- Re-optimization uses `ratchIter %/% 2` to avoid excessive compute
- Minimum 4 tips to allow dropping
- Warning when `neverDrop != TRUE` and `normalize = FALSE`

---

## Test data

Primates bootstrap trees (189 tips, 100 trees) from Takazawa et al. (2026)
PhyloCRISP repository, subsampled to 50 tips.

Saved at: `inst/test-data/takazawa/primates_bootstrap.nw`
Reference: `inst/test-data/takazawa/primates_reference.nw`

---

## Remaining work

1. Write the implementation (code written but R session hung before testing)
2. Test on toy rogue example (both normalize modes)
3. Test on primates data (profile, compare normalize modes)
4. Add tests to `tests/testthat/test-CIDConsensus.R`
5. Update DESCRIPTION (Rogue to Suggests if TipInstability used)
6. Stretch: clade dropping (cherries with unstable position)
