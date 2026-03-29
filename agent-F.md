# Agent F — Progress Log

## Current State

- **Status:** IDLE
- **Date:** 2026-03-29 ~14:00 GMT

### F-030 — COMPLETE (PR #239)

TBR clip-ordering: Phase 2 full propagation + documentation fix.
Branch: feature/weighted-clip-order. Commits 5a060b92 + ca8f4f08.

**Root cause of Phase 1 null result:**
`clip_order` was only wired to ~10% of TBR calls (Wagner warmup + final
TBR polish). The dominant phases — ratchet (76% of replicate time) and
sectorial search (XSS/RSS/CSS) — always used RANDOM, making all ordering
variants empirically inert.

**Changes:**
- `ts_ratchet.h`: added `clip_order` field to `RatchetParams`
- `ts_sector.h`: added `clip_order` field to `SectorParams`
- `ts_driven.cpp`: propagated from `SearchControl` to all construction
  sites for both structs
- `ts_ratchet.cpp`: both initial TBR and perturbed-phase TBR now use
  `clip_order`
- `ts_sector.cpp`: added `clip_order` param to `search_sector()` helper;
  propagated to all 6 internal TBR call sites
- `R/SearchControl.R`: added `@param clipOrder` documentation
- `man/SearchControl.Rd`: corrected stale `\usage` block + added
  `\item{clipOrder}`

**Benchmark (5 seeds × 30s):**
- Zhu2013 75t thorough: +12% replicate throughput (TIPS_FIRST)
- Dikow2009 88t thorough: +8% replicate throughput (TIPS_FIRST)
- Agnarsson2004 62t default: −2% (TIPS_FIRST slightly harmful at low
  enrichment; neutral at this scale)

**No preset changes.** `clipOrder = 0L` (RANDOM) remains default.
Users can opt in with `SearchControl(clipOrder = 2L)`.

**GHA:** 23708949592 PASS.
**PR #239** open → cpp-search.

### Previous: T-245 — COMPLETE (PR #238, merged 2026-03-28)

TBR 4-wide candidate batch + flat-variant switch.

**Hamilton benchmark:** pending (feature/tbr-batch vs cpp-search,
mbank_X30754 + syab07205_206t, 60s/120s, 10 seeds, EW).
