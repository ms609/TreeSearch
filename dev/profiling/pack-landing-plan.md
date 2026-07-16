# TS_PACK_LOCAL landing plan (default-ON) — 2026-07-16

## EXECUTION STATUS (2026-07-16, in-session)
- **Step 1 DONE** — `CharBlock::plane_state[MAX_STATES]` (ts_data.h), populated in build_dataset
  (ts_data.cpp, inverts `loc[]` through `state_remap`), consumed at ts_rcpp.cpp:494. Audit for
  other plane→label emitters: CLEAN (only `state_str` in `ts_na_debug_char`; all other consumers
  use the loop index for index-agnostic bitmask access — verified).
- **Step 2 DONE** — default-ON flip at ts_data.cpp (`TS_PACK_LOCAL=0` forces off); harness OFF-arms
  switched from `unsetenv` to `setenv("0")` in `reeval/pack_{gate,reach,mpt_gate}.R`.
- **Step 3 (test gate) — key finding + results:**
  - `ts_na_debug_char` (the only C++ state-LABEL emitter) is DEBUG-ONLY, package-INTERNAL (not in
    NAMESPACE), and called by NO test or R business logic. The user-facing reconstruction
    `PlotCharacter` is PURE R (Brazeau2019 passes on its own state matrix) — packing does NOT touch
    it. So packed reconstruction changes nothing users see.
  - Packed reconstruction ambiguous SETS legitimately differ from unpacked (`0/1` vs `0/1/2/3`):
    packing restricts `?` to the char's local alphabet, so states never observed for a char (2,3)
    are dropped — which is CORRECT (score + candidates_evaluated byte-identical proves those states
    are never in an MPR). NOT the plane-relabel bug. So "reconstruction byte-identical" was the
    WRONG oracle. The relabel fix is validated by EXACT-label assertion instead (below).
  - NEW permanent test `tests/testthat/test-ts-pack-local.R`: (1) score-exact TreeLength ON vs OFF
    on Sansom2010; (2) plane_state relabel via a hand-built non-contiguous-alphabet fixture (block
    packs 3→2 planes; a "2"-tip must print "2" not the plane index "1"). 6/6 assertions PASS. This
    is the ONLY gate that exercises plane_state[] (scorer is index-agnostic).
  - `pack_mpt_gate.R` Sansom2010 EW/IW/NA: BYTE-IDENTICAL ON vs OFF (score+ntree+cand, 3 seeds each;
    NA reproduced the 66-tree MPT set). Score-exactness holds with the flip.
  - Full testthat suite default-ON: RUNNING (flip-safety gate) — commit only after it passes clean.
- **Step 4** — commit after suite passes; then collect 17895591/17887591 and apply the rule.

---


**Decision (user-approved):** make the `state_str` reconstruction fix, then **immediately flip
`TS_PACK_LOCAL` default-ON**. Let the in-flight 5432 run confirm post-hoc; if 5432 shows a
*systematic* reach loss, revisit (revert default / keep opt-in). Rationale + evidence:
`kernel-1p9ns-findings.md` "RULING" section. This file is the executable checklist (written
pre-compaction so it survives).

## Why default-ON is justified now
Score-exact per tree (PROVEN: `TreeLength` byte-identical OFF vs ON on every dataset incl. 5432).
The *only* thing packing perturbs is which equal-cost tie the Wagner **start** breaks — an unbiased
reshuffle (mean-score Δ = 0.00 across 60 small-corpus runs + 6 char-rich runs). Reach 18/18 keys
equal so far; speed 1.10× (small) → 1.27×/1.49× (char-rich 716c/2954c), climbing with nChar. Only
`project510` (decision-critical, most char-rich) and `project4085` of the big-tip set have landed;
the rest are confirmatory. Packing can NEVER find a *higher* score — so the only possible veto is a
reach *shortfall*, and there's no mechanism for a systematic one.

## Step 1 — fix `state_str` (ts_rcpp.cpp:488) — plane→global relabel
**Bug:** in `ts_states_at_pattern`, applicable plane `st` is labelled `st - has_inapp`. That equals
the global applicable index ONLY when planes use the global alphabet. Under packing, applicable
plane `p` (= `st - has_inapp`) = **the p-th set bit of the block's local alphabet mask**
(`block_app_mask[b]`), i.e. a different global state → mislabel.

**Fix (mode-agnostic, precomputed map):**
1. `src/ts_data.h:74` `struct CharBlock` (after `pattern_index[]`, ~line 87): add
   `int plane_state[MAX_STATES];  // applicable plane p -> global-applicable label to display`.
2. `src/ts_data.cpp` build_dataset, right after `blk.n_states = ...` (line 305) and with
   `blk_app`/`block_app_mask` + `state_remap` in scope (state_remap is built at line 346 — either
   move plane_state population below it, or recompute the skip-inapp index inline). For each
   applicable plane `p` in `0..(local/max_app-1)`:
   - packed (`pack_local`): global state `s` = the p-th set bit of `blk_app`; `plane_state[p] =
     state_remap[s]` (= its global applicable index, the label the unpacked build shows).
   - unpacked: `plane_state[p] = p` (identical to current `st - has_inapp` behaviour).
3. `src/ts_rcpp.cpp:494` change the applicable-state label line to
   `s += std::to_string(blk.plane_state[st - (blk.has_inapplicable ? 1 : 0)]);`
   (inapp plane 0 still prints "-"). Now correct in BOTH modes with no runtime-flag logic in the
   consumer.

**AUDIT before flipping:** grep for OTHER plane→global label sites that assume plane==global
(ancestral-state reconstruction, MPT state output, any `std::to_string` of a plane index in
ts_rcpp.cpp / reconstruct paths). `state_str` is the known one (per memory) — confirm no others.

## Step 2 — flip default-ON (ts_data.cpp:131)
- Current: `const bool pack_local = std::getenv("TS_PACK_LOCAL") != nullptr;`
- New (default ON, explicit OFF via `=0`):
  `const char* pl = std::getenv("TS_PACK_LOCAL");`
  `const bool pack_local = (pl == nullptr) || !(pl[0]=='0' && pl[1]=='\0');`
- **Consequence:** A/B harnesses that used *unset = OFF* must switch the OFF arm to
  `Sys.setenv(TS_PACK_LOCAL="0")`: `dev/profiling/reeval/pack_{gate,reach,mpt_gate}.R`,
  `.agent-pack/*`. (`TS_PACK_SORT` gate at :136 unaffected.)
- **DO NOT rebuild `/nobackup/pjjg18/packbuild/lib`** until the in-flight jobs finish — they use the
  OLD default-OFF lib with unset=OFF/=1=ON semantics; rebuilding mid-flight corrupts their OFF arm.

## Step 3 — test gate (before commit)
- Clean build (stale-ABI gotcha: `CCACHE_DISABLE=1 R CMD INSTALL --preclean`, or fast build then a
  clean validate). Watch for `~/.R/Makevars.win` zeroing flags.
- `dev/profiling/reeval/pack_mpt_gate.R`: **Sansom2010 EW/IW/NA byte-identical** default-ON vs
  forced-OFF (the score-exact-per-tree property; this is the correctness anchor).
- Full `testthat` — reconstruction/state-inspector tests now hit the packed path by DEFAULT; the
  Step-1 fix must make them pass. Run a multi-block dataset through the `ts_states_at_pattern` R
  wrapper packed vs `TS_PACK_LOCAL=0` → labels identical.
- Sanity: a couple of `MaximizeParsimony` scores match a known floor.

## Step 4 — commit + close
- Commit on branch `claude/mission-b-kernel-goloboff-6f5c0e` (never `git add -A` — concurrent-session
  hazard; stage `src/ts_data.cpp src/ts_data.h src/ts_rcpp.cpp` + the harness OFF-arm edits + docs).
- Task #7 → completed. Update `mission-b-kernel-refuted.md` (lever LANDED default-ON).

## In-flight jobs to collect post-compaction (Hamilton)
- **17895591** `packtailA` (`shared`/48h, skip-if-CSV guard): 10 tail keys ×3 seed ×2 mode. DONE
  already: project4085 (reach 4051 3/3 Δ0, 1.27×), project510 (18345 3/3 Δ0, 1.49×). Pull
  `/nobackup/pjjg18/packbuild/out/pt_*.csv`; aggregate via scratchpad `tail_agg.R`.
- **17887591** `packtailB` (`long`, project5432 482t, ~8.5h+ elapsed): 6 tasks → `pt_project5432_*`.
  THE genuine reach test (only reach-limited key). Apply pre-committed rule: net reach-neutral →
  default-ON stands; systematic loss → revisit.
- Ruling doc + memory already record the 12-cell + char-rich results and the 3× walltime lesson
  ([[hamilton-partition-choice]]).
