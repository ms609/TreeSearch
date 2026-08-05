# S7 — sizing the plain-EW TBR fast path (#7) + pinning the TNT target

**Date:** 2026-07-13/14 · **Base tip:** cpp-search `d8c59998` (agent worktree, not pushed) ·
**Mission link:** the ~24× per-move gap to TNT on project5432 (482t morphological).

Two linked measurements:
- **Task A** — how much of TS's per-candidate cost is *strippable general-engine orchestration*
  that a plain-EW search pays for nothing (the #7 lever, hypothesised ~10×)?
- **Task B** — disassemble the **64-bit** Hamilton TNT (the one that hits ~1.9 ns) to resolve,
  not infer, how it gets the speed.

---

## VERDICT (headline)

- **Task A — #7 REVISED DOWN, hard, from "~10× big lever" to ≈ 5% of whole-search per-move cost.**
  A byte-identical strip of every dead-in-plain-EW per-candidate branch (category (ii):
  `sector_mask`, `constrained`, `use_collapsed`, `has_na`/`use_iw` dispatch, `b2_ceiling`) recovers
  only **~0.5–1.1 ns/candidate on the SPR scan (~4–6% of that path)**. Adding the per-block
  bookkeeping strip (category (iii): the flat scorer drops CharBlock deref + weight-multiply +
  `active_mask` check — byte-identical on unit-weight data) roughly triples the SPR saving to
  **~17% of the SPR scan ≈ 4.7% whole-search** (Zhu2013). The SPR scan is only ~19–23% of
  candidates, and the dominant reroot path is *already* specialised (branch-free x4 kernel when
  `all_weight_one`; a skip-predicate-hoisted scalar loop otherwise). Because branch/bookkeeping cost
  is **fixed** (data-/size-independent) while the scorer's gather grows with n_tips, this fraction
  *shrinks* at 5432 scale. **Not a per-move lever (~1.05×, not 10×).** The residual per-candidate
  cost is the scorer's gather + SIMD reduce + popcount — representation/memory, the irreducible
  category (i)+(iii)-core.
- **Task B — the task's premise is CORRECTED; the mechanism is narrowed, not fully pinned.** The
  64-bit TNT is **SSE2-vectorised** (auto-vectorised GCC), not scalar as the 32-bit T-250 disasm
  was; it uses **no hardware `popcnt`, no 64 KB LUT, no Hamming-weight SWAR**, and packs characters
  into **2-bit / 4-bit fields** (up to 64 / 32 chars per 128-bit vector). But it is *stripped*, and
  the **per-reconnection incremental scorer itself was not located** — the loop read is the
  full-tree downpass, not the per-candidate kernel. So "TNT-64 uses SIMD" is a fact, and it is
  explicitly **not** the source of the 1.9 ns; the residual advantage is consistent only with
  **touching far fewer words per candidate** (incremental reoptimisation + score cutoff), which a
  stripped-binary static read cannot prove and which prior TS work ([[tbr-per-move-kernel-gap]],
  L3b/quick-TBR) found infeasible on morphological data (no cross-clip locality).

Net: **#7 is not the lever.** The per-move gap is neither strippable orchestration (measured ~1–4%)
nor per-word arithmetic (T-250: TS SSE2 ~4×/word faster). It is the constant-factor cost of TS's
per-candidate re-score-with-gather vs TNT's incremental few-word update. This *sharpens* the prior
"per-candidate path at-limit" conclusion and re-points the only real headroom at a
representation/incrementalisation rewrite — which prior work already judged infeasible here.

---

## TASK A — sizing the strippable overhead

### Per-candidate work, categorised (`src/ts_tbr.cpp`, `src/ts_fitch.cpp`)

Plain-EW per-move work = an **SPR scan** (scalar, one candidate at a time, ~19–23% of candidates)
+ a **reroot loop** (~77–81%), around a per-clip precompute. Feature flags (`has_na`, `use_iw`,
`use_flat`, `use_collapsed`, `sector_mask`, `constrained`, `b2_ceiling`, `revert_check`, …) are all
hoisted to loop-invariant `const bool` at the function top (each `getenv` read once); the historic
per-clip `getenv("TS_REVERT_CHECK")` (13–19% Windows wall) is already merged (beb52138/846cc773) —
a **per-clip**, platform-specific cost, not per-candidate orchestration.

**SPR scan — per candidate:** identity skip (i); then the (ii) dead checks
`sector_mask && …` (null), `constrained && …` (false), `use_collapsed && collapsed[below]`
(non-empty vector but **all flags 0** for morphological data), `if(has_na)…else if(use_iw)…else`
(false/false); then the scorer `fitch_indirect_length_cached(clip_prelim, edge_set_buf[below*tw])`
(i)+(iii); `++n_evaluated`; `if(b2_ceiling){…}` (false); accept `if(candidate<best_candidate)`.

**Reroot loop — dataset-dependent:**
- `all_weight_one` (e.g. Zhu2013): `use_flat=true` → the T-245 **x4 batch kernel**
  `fitch_indirect_cached_flat_x4` — four data-independent `any_hit_reduce`+`popcount` with one
  combined-cutoff bail. Skip predicate hoisted once/clip into `kept_ei`. **Branch-free; nothing
  to strip.**
- weighted (e.g. Zanol2014, Dikow2009): `use_flat=false` → a **scalar** `kept_ei` loop that still
  carries a per-candidate `has_na`/`use_iw` dispatch (~2 dead branches), but the skip predicate is
  hoisted. (This is why the flat scorer can't be reused for weighted data — see below.)

### Experiment — controlled A/B, trajectory-identical by construction

Runtime toggle `TS_EW_STRIP` (uncommitted probe in `ts_tbr.cpp`, default OFF) selects a specialised
SPR loop keeping only the identity skip + **the same** `fitch_indirect_length_cached` scorer +
accept, removing all (ii) branches. Because the removed branches never fire in plain EW, the
candidate sequence and every strict-`<` accept are identical → **score + `candidates_evaluated`
byte-identical**. The untouched reroot path is the **control**: its BASE↔STRIP delta = residual bias.

**Positive control (the essential check — this experiment initially mis-fired).** The gate first
required `!use_collapsed`, but `compute_collapsed_flags` does `collapsed.assign(n_node,0)`, so
`use_collapsed==true` even for morphological data (all flags 0) → the fast loop **never ran** and an
early "0% effect" was a **BASE-vs-BASE artifact**. Fixed the gate to `collapsed_all_zero` (checked
once; the vector isn't mutated in the clip loop) + added a fire-counter (`sprfast=N/clips`) and a
deliberately-wrong mode (`TS_EW_STRIP_WRONG`, `extra+=1`). On Dikow2009:

| mode | score | tbrMoves | sprfast |
|---|---|---|---|
| BASE  | 1608 | 232 | **0**/5854 |
| STRIP | 1608 | 232 | **5854**/5854 |
| WRONG | 1609 | 231 | 4655/4655 |

STRIP fires on 100% of clips **and** matches BASE (path live + output used correctly); WRONG diverges
(output genuinely drives the search). Only after this did the timing A/B become valid.

- **Panel:** Zhu2013 (75t, `all_weight_one`), Zanol2014 (74t, weighted), Dikow2009 (88t, weighted),
  each recoded `-`→`?` (`to_fitch`) so `has_na=false / use_iw=false` (verified levels `0,1,2,3`,
  `has_dash=FALSE`). Panel is the **conservative** place to cap #7: branch cost is fixed, gather cost
  grows with n_tips, so branch-fraction(panel 74–88t) ≥ branch-fraction(5432 482t).
- **Instrument:** in-tree `TS_IW_TIMING` chrono (per-candidate ns for SPR and REROOT separately).
- **Method:** in-process interleaving (BASE/STRIP ratchet calls back-to-back, same seed → identical
  trajectory), single-thread, 40 pairs/dataset. The box was under heavy variable load (per-rep swings
  ±20–40% on *both* paths), so wall medians are unusable; **load only ever adds time**, so
  **minimum-of-runs** is the load-robust true-cost estimator and the untouched REROOT min-delta
  corrects residual bias. (Byte-identity re-verified: BASE==STRIP score+tbrMoves for every dataset.)

### Results — byte-identical dead-branch strip (min-of-runs, ns/candidate)

| Dataset (n_tip, weights) | SPR base→strip | REROOT (control) | **SPR-specific** | SPR abs. saving | SPR share |
|---|---|---|---|---|---|
| Zhu2013 (75, unit) | 14.99→13.93 (−7.0%) | 9.06→9.00 (−0.6%) | **−6.5%** | **1.06 ns** | 18.5% |
| Zanol2014 (74, wtd) | 18.42→17.71 (−3.9%) | 15.06→15.13 (+0.4%) | **−4.3%** | **0.71 ns** | 23.4% |
| Dikow2009 (88, wtd) | 25.72→25.20 (−2.0%) | 20.68→21.20 (+2.5%) | **−4.6%** | **0.52 ns** | 18.7% |

Control (untouched reroot) stays within ±2.5%; SPR-specific = SPR delta − control delta. The strip
recovers **~0.5–1.1 ns/candidate (~4–6%) on the SPR scan**. ~0.5–1 ns ≈ ~2–4 cycles — consistent
with removing ~6 predicted-not-taken branches + the b2 block, i.e. work that mostly overlaps the
scorer's memory latency. (The absolute ns varies ~2× across datasets and is comparable to the
residual control drift — that is the estimator's floor, not signal; the hedged ranges are the honest
resolution.)

### Results — (iii) per-block bookkeeping strip (flat scorer), byte-identical on unit-weight data

Additionally swapping the SPR scorer to `fitch_indirect_cached_flat` drops the per-block CharBlock
deref, weight-multiply and `active_mask` check. For `all_weight_one` data this is **exactly**
`fitch_indirect_length_cached` (byte-identical — verified on Zhu2013: single trajectory 626/221 for
all 40 pairs), so it is a legitimate category-(iii) strip there:

| Dataset | maximal (ii)+(iii) SPR base→strip | REROOT ctrl | **SPR-specific** | byte-identical? |
|---|---|---|---|---|
| Zhu2013 (unit) | 15.85→13.13 (−17.2%) | 9.55→9.57 (+0.1%) | **−17.4%** | **yes** (626/221) |
| Dikow2009 (wtd) | 27.07→23.36 (−13.7%) | 22.01→22.25 (+1.1%) | −14.8% | **no** — weight-blind (307↔299) |

So on unit-weight data the (iii) bookkeeping is worth ~11% of SPR on top of the (ii) branches
(−6.5% → −17.4%). On **weighted** data the current flat kernel is **weight-blind** → wrong lengths,
so it reaches the same score via a different trajectory (Dikow deterministically 307 vs 299): NOT
safe as-is; a weight-aware flat block would be needed, and the reroot path is scalar there too so
both paths carry the same (iii) headroom.

### #7 ceiling

- **SPR path (~19–23%):** (ii) branches ~0.5–1.1 ns (~4–6%); (ii)+(iii) with flat scorer ~17% of
  SPR (Zhu, unit-weight, abs 2.72 ns), byte-identical.
- **Reroot path (~77–81%):** `all_weight_one` → branch-free x4, **already (ii)+(iii)-optimal**;
  weighted → scalar loop with the same (ii)+(iii) headroom (unmeasured, bounded ≤ SPR).
- **Whole-search ceiling:** Zhu2013 (unit-weight, reroot already x4) = **~4.7%** (SPR 2.72 ns ×
  18.5% / 10.7 ns all-in). Call it **~5% best case, byte-identical**, ≈ a **1.05× speedup**. At
  5432 (482t) SPR candidates cost ~46 ns and the gather dominates further (prior perf: ~16 L1 + ~4
  LLC misses, ~211 cyc/candidate), so the *fixed* branch+bookkeeping saving is a *smaller* fraction —
  the panel cap holds a fortiori.

**#7 (branch + per-block-bookkeeping strip) is refuted as a mission lever (~5%, not ~10×).**
Concordant with the independent `tbr-microlever-sweep` (the same (ii)-class hoists measured +0.00%;
the only real win was a per-clip getenv) and the `abapprox` A/B (the shared reduce+gather is the
floor). The recoverable overhead is fixed and ≤~5%; the per-candidate floor is irreducible scorer
work (gather + SIMD reduce + popcount). The real headroom, if any, is a representation/
incrementalisation rewrite (Task B) — which prior work already judged infeasible on morphological
data.

---

## TASK B — the 64-bit TNT per-move kernel (static disassembly on Hamilton)

**Binary:** `/nobackup/pjjg18/TreeSearch/tnt/TNT-bin/tnt` — ELF **64-bit** LSB PIE, x86-64, dynamic,
**stripped** (no symbols), 4,373,832 bytes. **Disassembly (durable on Hamilton):**
`/nobackup/pjjg18/tnt_disasm_64.txt` (AT&T, 807,780 lines) + `…_intel.txt`. Hot loops located by
density-scanning `and`/`or`/`movzx` clusters at short backward jumps (no symbols to anchor on).

### Grep-proven, high confidence

- **(a) Step counting: NEITHER hardware popcnt NOR LUT NOR Hamming-SWAR.** `popcnt` count = **0**;
  no 64 KB byte-LUT in the hot regions; the SWAR masks `0x33333333`, `0x0f0f0f0f` and the
  `imul …,0x01010101` finaliser are **all absent**. Only `0x55555555`/`0x77777777`/`0x88888888`
  appear — field-*select* tricks (detect empty 2-bit / nonzero 4-bit fields), not a step counter.
  **The length-accumulation routine was not located; its counting mechanism is undetermined.**
- **(b) SIMD, not scalar — SSE2 baseline only.** The Fitch state-set combine is genuinely
  SSE2-vectorised (`movdqa/movdqu`, `pand/pandn/por/pxor`, `pslld/psrld`, 128-bit stores),
  auto-vectorised (`cmp count,0x5; seta` + pointer-alias `setae` + scalar 32-bit remainder = the GCC
  auto-vec signature). **No AVX/AVX2/AVX-512** (`ymm`=`zmm`=0), **no SSE4**
  (`pshufb`/`ptest`/`pmovmskb`/`blendv`=0), **no BMI** (`lzcnt`/`tzcnt`/`pdep`/`pext`/`andn`=0).
  **⇒ corrects the task's "prior 32-bit was scalar" premise for the 64-bit build.**
- **(c) State representation: WIDE PACKED BLOCKS**, two widths coexisting, each with its own kernel:
  **2-bit fields** (mask `0x55555555`, `psrld/pslld 1`) ≈ 64 chars / 128-bit vector, and **4-bit
  fields** (`0x77777777`/`0x88888888`, `psrld 3`) ≈ 32 chars / 128-bit. This is the *opposite* of
  "compact per-character scalar with early bail." The located loop has **no per-character
  early-exit** — only a word-count bound (`cmp count,idx; ja`).
- **(e) Throughput of the located loop:** ~22 instructions per 128-bit iteration, i.e.
  ~0.34 instr/char (2-bit) or ~0.7 instr/char (4-bit), fully branchless.

### What the disasm does NOT resolve (honest limit)

The densest loop read is the **full-tree preliminary-state downpass**, **not** the per-reconnection
candidate scorer (with its length accumulator + cutoff early-bail), which lives in a distinct routine
not locatable in a stripped binary. So **(d) incremental-slide vs full combine** and the
per-reconnection instruction count / whether it is incremental remain **open**.

### Reconciliation (do not over-read)

TNT-64 uses SSE2 yet is ~24–30× faster per move than TS (also SSE2/AVX2, and per T-250 ~4×/word
faster arithmetically). Both hold only if **TNT's advantage is touching far fewer words per
candidate** — incremental reoptimisation + a score cutoff — **not** per-word SIMD throughput. That is
exactly the "quick-TBR / incremental views" family prior TS work found **infeasible on morphological
data** (no cross-clip locality; L3b footprint 41–68%; the irreducible up-pass TNT also pays). So the
64-bit SIMD finding is genuine and corrects a premise, but it is **not** the explanation for 1.9 ns
and does **not** reopen a buildable TS lever. Fully pinning the loop shape needs a **symbolised/debug
TNT** or **dynamic profiling** (`perf record`/Callgrind on a live search) — static grep on a stripped
binary has reached its limit. Given prior work already judged the copy-target infeasible here,
stopping static analysis is a decision, not an exhausted method.

---

## Harness / provenance

- Probe: `src/ts_tbr.cpp` runtime toggles `TS_EW_STRIP` (default OFF; byte-identical branch-strip)
  + `TS_EW_STRIP_WRONG` (positive-control corruption) + `sprfast=N/clips` in the `TS_IW_TIMING` line
  + `collapsed_all_zero` gate. Uncommitted; not a production path.
- Scripts (session scratchpad): `s7_ab.R` (multi-process), `s7_inproc.R` (in-process interleaved),
  `s7_parse.R` / `s7_parse_paired.R` / `s7_min.R` (aggregate / paired / load-robust min-of-runs).
- Per-agent lib `.agent-s7` (full install + `dev/build-fast.R` hot-swap). Single-thread throughout.
- TNT disasm durable on Hamilton: `/nobackup/pjjg18/tnt_disasm_64.txt` (+ `_intel`).
- Builds on: `tbr-per-move-kernel-gap`, `tbr-microlever-sweep.md`, `getenv-ucrt-cost`,
  `sectorial-throughput-at-limit`.
