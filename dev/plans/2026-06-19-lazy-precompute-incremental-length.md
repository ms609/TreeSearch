# Lazy / incremental-length scoring — the per-candidate PRECOMPUTE lever (task #46, "Big & Hard")

**This file is the complete, self-contained spec for an UNATTENDED autonomous
session.** It launches fresh at ~21:32 on 2026-06-19 with no memory of the
conversation that wrote it. Read this top-to-bottom first, then
`dev/profiling/findings.md` (row **T-P5k**), then the
`MEMORY.md` index in the project memory dir. Everything you need is here or in
those two places.

**▶ SPAWNED-TASK NOTE — READ FIRST (you were launched as a background task chip).**
You are running in YOUR OWN fresh git worktree, spun off the `cpp-search` branch.
Everywhere this plan says "the l3b worktree," substitute **your own spawned
worktree** as your WORKING / build location and commit target. The
`claude/l3b-edgeset-timing` worktree
(`C:\Users\pjjg18\GitHub\worktrees\TreeSearch\l3b-edgeset-timing`) is only a
READ-ONLY reference for the footprint-instrumentation template (its
`src/ts_tbr.cpp`, the `#ifdef TS_EDGESET_FOOTPRINT` block ~line 1568); your M1
`-DTS_EDGESET_CONSUMED` measurement is NEW code you add in your own worktree.
This plan file, `findings.md`, and `MEMORY.md` live in the MAIN checkout
(`C:\Users\pjjg18\GitHub\TreeSearch\...`) — read them, and APPEND your PROGRESS
LOG to THIS file, by absolute path. **Before starting, read the PROGRESS LOG at
the bottom: if a prior run already completed milestones, resume from the last
checkpoint or, if M5 is reached, stop and report "already complete" (do not
redo).** Commit only to your own spawned worktree's branch; the
no-push / never-touch-cpp-search-or-main / never-merge rules in §3 apply
UNCHANGED.

---

## 0. Mission & where this sits

**Mission:** close TreeSearch's per-iteration wall-clock gap to TNT 1.6 on
**equal-weights (EW) Fitch** parsimony. The *quality* gap is closed (gapB = 0);
only throughput remains. NA/IW is OUT OF SCOPE (another agent owns it). Fitch is
triggered in drivers by replacing `-` with `?` (`m[m == "-"] <- "?"` →
`has_na = FALSE` → the standard Fitch path).

**The finding that motivates this work (T-P5k, VTune, post-L1, Zanol2014 full-EW):**
The per-candidate scoring cluster is ≈ ½ of EW CPU, and it splits:

| Half | What | Share of cluster |
|---|---|---|
| **PRECOMPUTE** | `compute_insertion_edge_sets` builds the FULL ~210-block directional view (`edge_set_buf`) for **every** candidate edge — two `combine` operator() lambdas + an O(N) up-pass | **~56%** |
| **CONSUMPTION** | `fitch_indirect_length_cached` scores each candidate, but **bails after ~2.5 of ~210 blocks** (95–99% of scorings bail early) | **~44%** |

**The lever (this task):** We eagerly materialise a ~210-block directional view
per candidate, but the bail-fast scorer reads only ~2.5 blocks of it. **~98% of
the eager `combine` work builds blocks no scorer ever reads.** That eager build is
exactly what TNT's incremental-length walk avoids (Goloboff quick-TBR). The ~2.5×
per-candidate throughput gap maps onto this avoidable work.

**Honest ceiling (do not over-sell):** the per-candidate cluster is only ~½ of EW
CPU (the rest = ratchet reweight ~60% of mission wall, sectorial ~30%, hashing, R
glue). So even a large win here yields **~1.2–1.5× end-to-end wall**, not 2.5×.
That is still a real, banked win worth pursuing — but record the realistic figure.

**The CONSUMPTION half (~44%) is being handled SEPARATELY and interactively** (a
scalar-inline micro-bench of `fitch_indirect_length_cached`). **Do not work on the
consumption/scorer path. Your job is the PRECOMPUTE half only.**

---

## 1. The approach — exact lazy-per-block precompute (measure-gated)

There are three distinct sub-levers. **You are building (a). (b) and (c) are NOT
for unattended execution** — if (a)'s gate fails, scope them in prose and STOP.

- **(a) Lazy-per-block EXACT precompute [THIS TASK].** Reproduce `edge_set_buf`
  **bit-identically**, but only compute the directional-view blocks that some
  surviving candidate actually consumes before bailing. Exact by construction —
  **no admissible bound needed, no scoring/accept change, trajectory stays
  bit-identical.** The win is skipping the combine work for blocks in no
  candidate's consumed set.
- **(b) Incremental-length / TNT quick-TBR.** Restructure the candidate loop to
  slide the regraft point along tree edges with O(1) deltas. Exact, bigger
  restructure, high integration risk. **Out of scope for unattended.**
- **(c) Bound-then-verify.** A cheap admissible Fitch-insertion bound screens
  candidates; exact `edge_set` only for survivors. Needs a derived+oracle-checked
  admissible bound. **Out of scope for unattended** (the bound is the risky front
  step the supervisor must check).

### Why (a) might or might not win — the GATE

The scorer bails per-candidate at ~2.5 blocks, but **different candidates bail at
different blocks**. What matters is the **UNION of consumed blocks across ALL
candidates in a clip loop**:

- If the union is a **small** fraction of the 210 blocks (the same few low-index
  blocks dominate every candidate's bail) → lazy-per-block skips most combine work
  → **GO**.
- If the union approaches **all** blocks (candidates bail on diverse, scattered
  blocks) → lazy-per-block saves ~nothing → **DEAD**; record and stop.

This is the same measure-first discipline that killed L3b (cross-clip reuse) by
direct measurement. **Reuse the existing L3b footprint instrumentation plumbing**
(`-DTS_EDGESET_FOOTPRINT`, `ts_edgeset_footprint_report()` Rcpp export, the
`l3b_regime` gate) — it is already wired in this worktree.

---

## 2. Milestones — each ends with a DURABLE checkpoint written to disk

A headless session can die on token exhaustion at any point (see MEMORY:
`dispatch-gotchas`). **After every milestone, append your findings + state to the
bottom of THIS file** (a `## PROGRESS LOG` section) so a relaunch — or the
supervising user — can pick up exactly where you stopped. Checkpoint *before* you
think you need to.

### M0 — Orient (cheap, ~5 min)
Read T-P5k in findings.md, the L3b plan (`dev/plans/2026-06-19-lever3-incremental-edgeset.md`),
this file, and MEMORY entries: `profiling`, `tbr-rooted-vs-unrooted`,
`fast-iteration`, `feedback-no-local-heavy-compute`, `concurrent-session-git-hazard`.
Confirm the worktree is `claude/l3b-edgeset-timing`
(`C:\Users\pjjg18\GitHub\worktrees\TreeSearch\l3b-edgeset-timing`). Confirm
`src/ts_fitch.cpp:479 compute_insertion_edge_sets` and the per-candidate loop in
`src/ts_tbr.cpp` (the `fitch_indirect_length_cached` call site, ~line 1528) are as
described. Write M0 checkpoint.

### M1 — THE GATE: measure the consumed-block union (build nothing else first)
Instrument the existing eager path (behind a new `-DTS_EDGESET_CONSUMED` flag,
mirroring the `-DTS_EDGESET_FOOTPRINT` plumbing) to record, **per clip loop**:
- For each candidate scored, which block index `b` it was on when it bailed (i.e.
  the highest block index `fitch_indirect_length_cached` touched before
  `extra_steps >= cutoff` or loop end). Instrument INSIDE that loop in
  `src/ts_fitch.cpp:452-474` under the flag, or pass back the bail-block via a
  thread-local counter (NOT `thread_local` in the hot path of the loaded DLL on a
  threaded run — but M1 runs `nThreads=1L`, so a plain file-scope static is fine).
- Accumulate the **union** of consumed block indices across all candidates in the
  clip loop, and `n_blocks_total`. Report median `|union| / n_blocks` per dataset.

Driver `dev/profiling/drivers/consumed_union.R`: full default-preset
`MaximizeParsimony(phy, nThreads=1L, verbosity=0L)` on Wortley2006 / Zanol2014 /
Zhu2013 (Fitch via `-`→`?`), `set.seed` fixed. Per-agent install (NEVER
`load_all`). Durable CSV to `dev/profiling/consumed_union.csv`.

**Decision (write the verdict to the PROGRESS LOG before acting on it):**

| median |union| / n_blocks | meaning | action |
|---|---|---|
| **small** (< ~0.35) | few blocks dominate all bails | **GO — build M2 (Scheme a).** |
| **≈ 1** | every block consumed by some candidate | **DEAD.** Record AT-LIMIT; lazy-per-block cannot win on this dataset class. Scope (b)/(c) in prose, STOP. |
| **mid** (0.35–0.7) | partial | GO but flag the realisable win is proportionally smaller; build and measure honestly. |

### M2 — Build the lazy-per-block prototype behind a flag (iff M1 = GO)
In `src/ts_fitch.cpp` add `compute_insertion_edge_sets_lazy(...)` (or a `lazy`
bool param) that, given the set of block indices to materialise (the consumed
union, or computed-on-demand), runs the up-pass + combine **only for those
blocks**, leaving other blocks of `edge_set_buf` untouched. Gate the call site in
`src/ts_tbr.cpp` behind `l3b_regime` + a new env/param flag
(`TS_LAZY_PRECOMPUTE`). The non-lazy path stays the default and untouched.

Two viable designs — pick the simpler that M1 supports:
- **Two-pass within a clip:** pass 1 = cheap scoring with whatever blocks are
  already materialised, collecting which blocks each candidate needs; pass 2 =
  materialise the union, then score exactly. (Risk: double-touch.)
- **On-demand single-pass:** the scorer requests block `b` of candidate `D`; a
  per-clip "is block b's full-tree up-pass done?" bitmap triggers a one-time O(N)
  up-pass+combine for block `b` across all nodes the first time any candidate needs
  it. Subsequent candidates reuse it. (This is the clean version — amortises each
  consumed block's up-pass once per clip.)

The on-demand single-pass is preferred: it makes the saving exactly "skip the
up-pass+combine for blocks in no candidate's consumed set," which is what M1
measures.

### M3 — Oracle: bit-identical correctness (MANDATORY before any wall claim)
Build an oracle config (`-DTS_LAZY_ORACLE`, asserts on): run BOTH the lazy path and
the from-scratch `compute_insertion_edge_sets` into separate buffers; `memcmp` the
consumed blocks every clip; `assert` bit-identical on Wortley2006 / Zanol2014 /
Zhu2013 × ≥3 seeds. Then the full kernel test suite on a per-agent install:
`R CMD INSTALL --library=.agent-lazy .` then
`testthat::test_dir('tests/testthat')` — must pass unchanged. **No wall claim from
the oracle/asserts build.**

### M4 — Wall A/B (local iterate-tier; deterministic; bit-identical trajectory)
Separate asserts-OFF wall build. Verify the from-scratch combine is actually
skipped (else no saving). Direct C-vs-A timing on the 3 datasets, single-thread,
fixed seed, **≥3 medians** (median of medians). Because Scheme (a) is exact and
order-preserving, **`candidates_evaluated` and final `score` must be bit-identical**
between lazy and baseline — assert this; any drift = bug, STOP and record.

Wall measurement is a single-threaded ~5s-bare run × a few medians — this is
iterate-tier and **stays local** (it is NOT CPU-greedy parallel compute; see
MEMORY `feedback-no-local-heavy-compute`). **Do NOT depend on Hamilton** (VPN may
be down, user unavailable). Record that final authoritative wall should be
Hamilton-confirmed later by the supervisor.

### M5 — Verdict
- **WIN:** record measured wall delta (honestly vs the ~1.2–1.5× ceiling), leave
  the prototype on the `claude/l3b-edgeset-timing` branch behind its flag,
  committed (NOT pushed). Write a findings.md row (T-P5l) and a clear
  "READY FOR SUPERVISOR REVIEW" note in the PROGRESS LOG with the exact build/run
  commands to reproduce.
- **DEAD:** record AT-LIMIT in findings.md + PROGRESS LOG with the measured numbers
  that killed it. Keep the instrumentation (gitignored Makevars; flags compile
  out). STOP.

---

## 3. AUTONOMY CONTRACT (you cannot ask the user — these are standing rules)

**You MAY, without asking:**
- Build per-agent libraries, run oracle + the ~30s targeted test suite locally.
- Run single-threaded deterministic wall A/B locally (iterate-tier).
- Write/append freely to: THIS file's PROGRESS LOG, `dev/profiling/findings.md`,
  `dev/profiling/*.csv`, new `dev/profiling/drivers/*.R`, new instrumentation in
  `src/*` **in the worktree**, the gitignored worktree `src/Makevars.win`.
- `git commit` to the **worktree branch `claude/l3b-edgeset-timing`** (staging ONLY
  named files you changed).

**You MUST NOT, under any circumstance:**
- `git push` anything, anywhere.
- Touch, commit to, merge into, or check out `cpp-search` or `main`. The prototype
  stays on the worktree branch for human review.
- `git add -A` / `git add .` / `git checkout -- .` / any broad git op — cpp-search
  is shared by parallel sessions (MEMORY `concurrent-session-git-hazard`). Stage
  **named files only.**
- Commit `src/Makevars.win` (gitignored, worktree-only — MEMORY/profile skill).
- Run CPU-greedy or long parallel local compute (MEMORY
  `feedback-no-local-heavy-compute`). No Hamilton dependency.
- Use `thread_local` for hot-path scratch in the loaded DLL (MinGW emutls —
  MEMORY `parallel-na-crash`). M1/M4 run `nThreads=1L` anyway; file-scope statics
  under the instrumentation flag are fine for the single-thread measurement.
- `devtools::load_all()` / `compile_dll()` for anything perf-measured — per-agent
  install only (MEMORY `fast-iteration`, r-conventions skill).

**Deliverable:** a measured verdict + (if WIN) a flagged prototype on the branch.
**Never a merge, never "ship it."** The supervisor decides merge.

**If blocked** (build won't compile, oracle won't pass, gate is ambiguous): record
the exact blocker + what you tried in the PROGRESS LOG and STOP cleanly. A clean,
well-documented stop is a success; a forced, untrustworthy result is a failure.

---

## 4. Conventions (project standing rules)

camelCase vars / BigCamelCase funcs / dot-private; UK 'ize'; `return x;` in C++
(the `# Return:` idiom is **R-only**); `@importFrom` not `@import`; base over
tidyverse; TreeTools over ape; `|>` not magrittr. Build with the symboled/release
flags from the `/profile` skill if profiling; otherwise a normal per-agent
`R CMD INSTALL`. Fast build: see MEMORY `fast-iteration` (ccache + incremental).

## 5. Key files
- `src/ts_fitch.cpp:452` — `fitch_indirect_length_cached` (the bail-fast scorer; M1
  instruments its bail-block).
- `src/ts_fitch.cpp:479` — `compute_insertion_edge_sets` (the eager precompute; M2
  adds the lazy variant beside it).
- `src/ts_tbr.cpp` ~1504 / ~1528 — per-clip `compute_insertion_edge_sets` call +
  the `fitch_indirect_length_cached` candidate loop (the call site to gate).
- `src/ts_rcpp.cpp` ~1655 + `src/TreeSearch-init.c` — mirror
  `ts_edgeset_timing_report` / `ts_edgeset_footprint_report` for any new report
  export.
- `dev/profiling/findings.md` — T-P5k (read), T-P5l (you write).
- `dev/plans/2026-06-19-lever3-incremental-edgeset.md` — the L3b plan + its DEAD
  verdict (the measurement template to copy).

---

## PROGRESS LOG
<!-- Append a dated checkpoint after EACH milestone. Newest at the bottom.
     Format: ### <ISO timestamp> — M<n> — <one-line state> then the details. -->

### 2026-06-19T22:?? — M0 — Oriented; OWN worktree created; strong analytic prior the gate will read DEAD

**Environment reality (IMPORTANT for any relaunch):** This session was launched
*in the MAIN checkout* (`C:/Users/pjjg18/GitHub/TreeSearch`, branch `cpp-search`,
carrying the parent session's uncommitted WIP) — NOT in an auto-created spawned
worktree as the launcher text implied. Per the AUTONOMY CONTRACT (never commit to
cpp-search), I created my OWN worktree:
- **Worktree:** `C:/Users/pjjg18/GitHub/worktrees/TreeSearch/lazy-precompute`
- **Branch:** `claude/lazy-precompute-m46` (off cpp-search HEAD da0f203f, clean)
- This is my WORKING/build/commit target. The PROGRESS LOG (this file) and the
  reference docs live in the MAIN checkout and I only READ them there + APPEND
  here (this file is untracked `??` in main, so appending never touches
  cpp-search's tracked state; I never `git add` it).
- Read-only reference for instrumentation plumbing:
  `C:/Users/pjjg18/GitHub/worktrees/TreeSearch/l3b-edgeset-timing` (commit
  00d73d6a) — `src/ts_edgeset_footprint.h` + `ts_edgeset_footprint_report()` in
  `src/ts_rcpp.cpp:1689` + registration in `src/TreeSearch-init.c:53`.

**Code confirmed (cpp-search da0f203f, src/ matches HEAD — not modified in main):**
- `src/ts_fitch.cpp:452` `fitch_indirect_length_cached(clip_prelim, vroot, ds,
  cutoff)`: loops `b` over `ds.n_blocks`, `continue` on `blk.active_mask==0`,
  accumulates `extra_steps`, and **`if (extra_steps >= cutoff) return;`** at
  ts_fitch.cpp:470. ⇒ each candidate consumes a **prefix** `[0, bail_block]` of
  the character blocks; the per-clip consumed UNION = `[0, max_bail_block]`.
- `src/ts_fitch.cpp:479` `compute_insertion_edge_sets`: eager — two `combine`
  lambdas build `edge_set[D]` for ALL non-root in-tree nodes D over ALL `nb`
  blocks (the up-pass + the edge-set pass). Called once per clip at
  `src/ts_tbr.cpp:1504` (gated `use_directional = !has_na`).
- **THE CRUX (kills design (a) a priori, pending measurement):** the SPR
  candidate loop at `src/ts_tbr.cpp:1537` initialises `best_candidate = HUGE_VAL`
  (ts_tbr.cpp:1529), so the **FIRST SPR candidate every clip is scored with
  `cutoff = INT_MAX`** (ts_tbr.cpp:1578-1580) → it never bails → it walks ALL
  blocks. The TBR reroot loop (ts_tbr.cpp:1819) and all later SPR candidates are
  bounded (best_candidate already set). So per clip there is ALWAYS ≥1 full-scan
  candidate ⇒ the raw consumed-union ≈ 1.0 by construction, compounded by the
  ~1% genuinely fully-scored candidates (T-P5h4). The on-demand single-pass lazy
  design (M2) builds block b for all edges on first demand, so ONE candidate that
  demands the last block forces all blocks built ⇒ **lazy-per-block saves ≈ 0 on
  the raw union.**

**Plan adjustment (still within (a)/M1):** M1 will measure and report BOTH:
1. **all-union** (raw, incl. the INT_MAX first candidate) — expected ≈ 1.0 ⇒
   confirms design (a) DEAD *as specified*.
2. **bounded-union** = max bail over candidates with `cutoff < INT_MAX` only —
   the informative "what the union would be IF the cutoff were seeded with the
   incumbent" number. This isolates whether the consumed set is genuinely
   concentrated (a smarter scheme could win) vs intrinsically broad.
3. Per-candidate bail-block distribution (histograms, split bounded/unbounded),
   fully-scanned fraction, and a **cost-weighted** union (weight each clip by
   `n_node`, the precompute-work proxy) so the realisable-saving estimate is honest.
This gives a defensible GO/DEAD verdict either way and correctly identifies that
the real enabling lever (seed the cutoff) is a *bound* change (sub-lever (c)
territory), not lazy precompute (a) — to be scoped in prose if (a) is DEAD.

**Guardrails reloaded:** per-agent `R CMD INSTALL` only; nThreads=1L; file-scope
statics (NOT thread_local) for instrumentation; stage named files only; no push;
no Hamilton dependency; never touch cpp-search/main. (MEMORY: fast-iteration,
feedback-no-local-heavy-compute, concurrent-session-git-hazard, profiling.)

**Next:** advisor sanity-check on the "first-candidate-unbounded ⇒ union≈1 ⇒ (a)
DEAD" reasoning + the bounded-union metric, THEN build the `-DTS_EDGESET_CONSUMED`
instrumentation (M1).

### 2026-06-19T22:?? — M1 — GATE = **DEAD** (and the spec's premise was factually wrong)

**Advisor (pre-build):** crux holds; gate on the RAW union (the spec's literal
metric); bounded-union is a diagnostic only; do NOT build M2 even if bounded-union
is small (a bit-identical-trajectory lazy prototype still eats the full
first-candidate scan); minimise build surface (avoid RcppExports/init.c codegen);
run Zanol first; document the "seed best_candidate" spin-off but don't implement.

**Built (worktree `claude/lazy-precompute-m46`, off cpp-search da0f203f):**
`-DTS_EDGESET_CONSUMED` instrumentation — NO Rcpp/codegen surface:
- NEW `src/ts_edgeset_consumed.h` — `EdgeSetConsumed` accumulator (per-clip union
  + per-candidate bail histograms), inline methods, env-var CSV writer.
- `src/ts_fitch.cpp` — global `g_edgeset_consumed` + `record()` hook in
  `fitch_indirect_length_cached` (captures bail block index + bounded/full flags).
- `src/ts_tbr.cpp` — `new_clip()` hook at the `compute_insertion_edge_sets` call
  site (ts_tbr.cpp:1504), brackets each clip's SPR+TBR scoring; counts n_active.
- `src/ts_rcpp.cpp` — RAII flusher in `ts_driven_search` → writes a key,value CSV
  named by env `TS_CONSUMED_OUT` on any return path.
- `src/Makevars.win` (gitignored) — `PKG_CPPFLAGS = -DTS_EDGESET_CONSUMED`.
- Driver `dev/profiling/drivers/consumed_union.R`; results
  `dev/profiling/consumed_union.csv` (+ per-dataset raw CSVs). Built per-agent
  `R CMD INSTALL --library=../_lib_lazy_consumed`; flag confirmed live on the
  compile line; nThreads=1L; default preset; maxReplicates=2; seed=1.

**THE FINDING THAT KILLS IT (and corrects T-P5k/T-P5h4): `n_blocks` ≈ 2–4, NOT
~210.** The "~210 blocks" in T-P5k/T-P5h4 is a FACTUAL ERROR — 210 was the
*character* count (Zanol). A `CharBlock` packs up to 64 characters
(`ts_data.h:73-78`, `n_chars` "1..64"), and the `FlatBlock` comment
(`ts_data.h:91,98`) literally reads *"For 4 blocks this fits in a single 64-byte
cache line."* The scorer loop `for (b < ds.n_blocks)` runs ~2–4 iterations, not
210. Measured `mean_nblocks`: Wortley **2**, Zanol **4**, Zhu **4**
(`max_union_all` matches: 2/4/4 — two independent confirmations).

| dataset | tips | n_blocks | **med per-clip union (GATE)** | mean union | cost-wt saving | bounded-only union | per-cand mean bail frac | n_clips | n_calls |
|---|---|---|---|---|---|---|---|---|---|
| Wortley2006 | 37 | 2 | **0.975** | 1.000 | **0** | 0.9999 | 0.755 | 11 812 | 0.76 M |
| Zanol2014 | 74 | 4 | **0.975** | 1.000 | **0** | 0.9999 | 0.653 | 49 680 | 28.7 M |
| Zhu2013 | 75 | 4 | **0.975** | 0.9993 | **0** | 0.9993 | 0.398 | 54 054 | 8.5 M |

**VERDICT: DEAD.** Gate threshold for GO was median(|union|/n_blocks) < 0.35;
measured **0.975** on all three (≫ threshold). The on-demand single-pass
lazy-per-block precompute (M2 design) can only skip combine work for blocks in NO
candidate's consumed set — but with only 2–4 blocks total and ~150–580
candidates/clip each consuming a prefix `[0, bail]`, the per-clip UNION of demanded
block *indices* saturates ~all blocks every clip ⇒ `saving_all_costwt = 0`. The
premise "~98% of eager combine builds blocks no scorer reads" is FALSE: there are
~4 blocks and ~all are read per clip. **M2/M3/M4 SKIPPED** (no prototype to build).

**The "seed the cutoff" idea is also killed — by a seed-INVARIANT argument**
(corrected; the earlier "bounded-union is a large upper bound on the seeded union
⇒ seeding fails" reasoning was a non-sequitur — a large upper bound says nothing
about the bounded quantity). The correct argument: the BEST regraft of a clip has
the lowest total extra-steps, so its partial sums never reach ANY
incumbent-derived cutoff → it is always fully scanned to the LAST block → the
union includes the last block → it saturates under *any* seed. Empirically
corroborated by `frac_bounded_full` = **1–8 %** of calls (Wortley 8.0 %, Zanol
1.0 %, Zhu 3.8 %) ⇒ ~several fully-scanned candidates per clip even with bounds
active. (Seeding `best_candidate` may still cut per-candidate WORK independently —
a separate, trajectory-changing lever — but it does NOT enable block-laziness.)
The bounded-only union numbers (0.9993–0.9999) are reported as corroboration, not
as the proof.

**NOT killed (scoped for supervisor, NOT built — autonomy contract):** a *different*
lever the data supports — **per-candidate lazy/fused edge-set**: fuse the edge-set
combine into the bail-fast scorer so `edge_set[D]` block b is materialised
on-demand during scoring and blocks beyond `bail(D)` are never combined (keeps the
full up-pass; no path re-walk; bit-identical). This is governed by the
PER-CANDIDATE bail (1 − mean_bail_frac unread), NOT the union: Wortley 24%, Zanol
35%, **Zhu 60%** of the edge-set-combine pass skippable, + eliminates the
`edge_set_buf` store. BUT: the edge-set pass is only ~10% of EW (T-P5k), the
combine is bandwidth/cache-bound (L2 AT-LIMIT) and the buffer is cache-resident at
37–75t (L3a AT-LIMIT washed for exactly this reason), so the realisable wall win is
uncertain (~2–6% EW, may wash). It is a moderate restructure (scorer signature +
the TBR `vroot_cache` path) ⇒ out of unattended scope; measurement-decidable, flag
for human. Per-candidate bail numbers above bound it.

**Next:** commit deliverables to the worktree branch; advisor reconcile (premise
was wrong + is per-candidate-lazy worth flagging as promising vs likely-washes);
write M5 verdict + findings row T-P5l.

### 2026-06-19T22:?? — M5 — **VERDICT: DEAD (measured). Lazy-per-block precompute saves zero. No prototype. STOP.**

Deliverables committed to worktree branch `claude/lazy-precompute-m46` @ `7c762305`
(NOT pushed, NOT merged): instrumentation (measurement-only, behind gitignored
`-DTS_EDGESET_CONSUMED`), driver `consumed_union.R`, results CSVs. Reproduce:
`R CMD INSTALL --library=../_lib_lazy_consumed --no-docs .` (with the worktree's
`src/Makevars.win`) then `TS_LIB=../_lib_lazy_consumed TS_DS=<ds> TS_REPS=2
Rscript dev/profiling/drivers/consumed_union.R`.

**HEADLINE (count-invariant, seed-invariant): the per-clip consumed-block union
fraction is ≈ 1.0** (`mean_union_all`: Wortley 1.000, Zanol 1.000, Zhu 0.9993).
Lazy-per-block precompute can only skip combine work for blocks in NO candidate's
consumed set; at union ≈ 1.0 that set is empty ⇒ `saving = 0`. **This verdict does
NOT depend on the block count** — whether n_blocks is 4 or 210, a saturated union
means zero saving. (The median 0.975 is the same quantity, just histogram-quantized
to "median clip lands in the top bin"; the fraction is the cleaner statement.)

**Premise correction (PROVEN, not inferred):** T-P5k/T-P5h4's "~210 blocks" is a
factual error — n_blocks ≈ 2–4. A `CharBlock` packs ≤64 characters
(`ts_data.h:73-78`); `FlatBlock` comment says "For 4 blocks ... single cache
line" (`ts_data.h:91,98`). Character/pattern counts (measured): Wortley **105**
→ ceil/64 = 2 blocks; Zanol **210** patterns → 4; Zhu **253** → 4; Giles **236**
→ 4. The 210 (Zanol) and 236 (Giles) match T-P5h4's two numbers EXACTLY ⇒ those
were pattern counts mislabeled "blocks". Self-corroborating: bail ≈ 2.5 blocks AND
union ≈ all blocks jointly force n_blocks ≈ 4 — so "2.5 of 210 blocks" is
internally impossible (210 real blocks couldn't saturate at a 2.5 bail).

**(a) lazy-per-block precompute [THIS TASK] — DEAD.** Gate threshold GO < 0.35;
measured union ≈ 1.0. M2/M3/M4 correctly skipped (no prototype). The on-demand
single-pass builds block b for all nodes the first time any candidate demands it;
with ~150–580 candidates/clip each consuming a prefix [0, bail] and only 2–4
blocks, the union of demanded indices spans all blocks every clip.

**Spec-required scoping of the other sub-levers on gate-fail:**
- **(b) incremental-length / cross-clip full-view maintenance — already DEAD**
  (L3b, T-P5j): per-clip changed-value footprint 41–68 %, Euler descend-delta
  1.14–1.24× the footprint ⇒ no cross-clip locality on this dataset class. Not
  reopened here.
- **(c) bound-then-verify / lazy-exact — UNTESTED, blocked on a prerequisite:**
  needs a cheap *admissible* (never-screens-an-improver) lower bound on Fitch
  insertion cost. None is established (the old union-of-finals bound OVERcounts —
  the very bug the directional fix cured). Out of unattended scope (the bound is
  the risky front step a supervisor must derive + oracle-check).
- **(d) per-candidate lazy/fused edge-set [NEW, the one the data actually
  supports, NOT built] — measurement-decidable, lean SKEPTICAL.** Distinct from
  (a): governed by the PER-CANDIDATE bail, not the union. Fuse the edge-set
  combine into the bail-fast scorer so `edge_set[D]` block b is materialised
  on-demand during scoring; blocks beyond `bail(D)` are never combined (keeps the
  full up-pass; no path re-walk; bit-identical). Skippable fraction of the
  edge-set-combine pass = 1 − mean_bail_frac: Wortley 24 %, Zanol 35 %, **Zhu
  60 %**; also eliminates the `edge_set_buf` store. BUT the edge-set pass is only
  ~10 % of EW (T-P5k), the combine is bandwidth-bound (**L2 AT-LIMIT**), and this
  exact buffer was found **cache-resident at 37–75t — L3a fused the two sweeps and
  saved NOTHING (T-P5i)**; design (d) adds per-block branching to a currently tight
  scalar loop, so it could wash or regress. Moderate restructure (scorer signature
  + the TBR `vroot_cache` path). Flag for human; do not assume a win.

**Mission context unchanged:** per-candidate PRECOMPUTE is NOT reducible via
block-laziness on this dataset class. The CONSUMPTION-half micro-bench (scalar vs
SIMD, handled separately/interactively) and recipe composition + the sectorial
~30 % (#39/#40) remain where the wall actually lives (T-P5c: ratchet ~60 %,
sectorial ~30 %, all TBR < 8 %).

---

#### findings.md row T-P5l  ·  **TRANSCRIBE TO MAIN findings.md** (do NOT let an unattended session edit the parent's uncommitted findings.md — concurrent-git-hazard)

```
| T-P5l | **P0** | **lazy-precompute lever (#46 (a)) DEAD by direct measurement; T-P5k premise ("~210 blocks") was a factual error** | T-P5k | [Algorithm, DECISIVE] Per-clip consumed-block UNION gate (worktree claude/lazy-precompute-m46, -DTS_EDGESET_CONSUMED, driver consumed_union.R, Wortley/Zanol/Zhu ×reps2) | **n_blocks ≈ 2–4, NOT ~210** — a CharBlock packs ≤64 chars (ts_data.h:73-78; "For 4 blocks ... single cache line"). T-P5k/T-P5h4's "210"/"236" were the Zanol/Giles PATTERN counts (measured 210/236 exactly), mislabeled "blocks"; "2.5 of 210 blocks" is internally impossible. **Gate = per-clip union fraction; GO needed <0.35, measured ≈1.0** (mean_union_all Wortley 1.000 / Zanol 1.000 / Zhu 0.9993; cost-wt saving = 0) on every dataset ⇒ the on-demand lazy-per-block precompute (skip combine for blocks no candidate reads) saves ZERO: with 2–4 blocks and ~150–580 cand/clip each consuming a prefix [0,bail], the union of demanded indices saturates all blocks. **Verdict count-invariant** (saturated union ⇒ 0 saving regardless of block count). Seed-INVARIANT too: the best regraft of each clip is always fully scanned to the last block (frac_bounded_full 1–8%) ⇒ seeding best_candidate cannot shrink the union. M2/M3/M4 skipped (no prototype). **(b) already DEAD (T-P5j cross-clip no-locality); (c) bound-then-verify still blocked on an unestablished admissible Fitch-insertion bound; (d) NEW per-candidate lazy/fused edge-set (skippable 1−bail = Wortley 24%/Zanol 35%/Zhu 60% of the ~10%-of-EW edge-set pass + buffer-store elimination) is measurement-decidable but lean-skeptical — L2 bandwidth-bound + L3a found this buffer cache-resident at 37–75t (fuse saved nothing) + adds per-block branching; moderate restructure, flag for human.** Instrumentation MEASUREMENT-ONLY (gitignored Makevars), nothing to merge; committed 7c762305 (NOT pushed). |
```

**READY FOR SUPERVISOR REVIEW.** Nothing to merge (measurement-only). The branch
`claude/lazy-precompute-m46` holds the instrumentation + driver + results for
reproduction/audit. Recommend: (1) transcribe T-P5l to findings.md; (2) correct
the "~210 blocks" wording in T-P5k/T-P5h4/T-P5j; (3) if any per-candidate lever is
pursued, it is (d) — measure it, do not assume a win.
