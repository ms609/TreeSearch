# Agent F — Progress Log

## Current State

- **Status:** IDLE
- **Date:** 2026-03-28 ~17:45 GMT

### T-245 — COMPLETE (PR #238)

TBR 4-wide candidate batch + flat-variant switch. Branch: feature/tbr-batch. Commit 038e00a8.

**Changes:**
- `ts_fitch.h/cpp`: added `fitch_indirect_cached_flat_x4()` (EW) and
  `fitch_na_indirect_cached_flat_x4()` (NA). Each processes 4 vroot_cache
  rows simultaneously per block, with shared early-exit when all 4 exceed
  cutoff (bitwise-AND test to avoid branch overhead).
- `ts_tbr.cpp`: `use_flat` flag computed once per `tbr_search` call
  (weight==1, no upweight_mask). SPR loop switched to flat variants.
  TBR rerooting inner loop: `use_flat && !use_iw` → batch-of-4 while loop;
  IW and ratchet-upweight paths fall back to existing scalar loop unchanged.

**Local tests:** 28 test-ts-tbr-search + 23 constraint-small PASS (NOT_CRAN=true).
**GHA:** 23690208221 PASS.
**PR #238** open → cpp-search.
**Hamilton benchmark:** pending (feature/tbr-batch vs cpp-search, mbank_X30754 +
syab07205_206t, 60s/120s, 10 seeds, EW).

### Also completed this session

- **spelling.R covr fix** (f1e9c4c5 on cpp-search): `withCallingHandlers` to
  muffle the "Failed to find package source directory" warning that was
  causing PR #210 R-CMD-check to fail in the Code coverage step.
- **S-RED focus 29**: reviewed ts_tbr.cpp (T-245 batch loop) + ts_driven.cpp
  (T-289 Stage 4 preset). No correctness bugs. Fixed stale PR comment in
  MaximizeParsimony.R (commit f1e9c4c5).
