# Component-isolation profiling program (2026-06-19)

## STATUS (updated 2026-06-21)

Progress against the per-component "two gates" (AT-LIMIT VTune + shared-start TNT
race) and the standing ordering:

| Component | Gate 1 — AT-LIMIT | Gate 2 — TNT race | Recipe lever | Verdict |
|---|---|---|---|---|
| **scoring (Fitch EW)** | ✅ AT-LIMIT (T-P5l: AVX2 reduce optimal at n_states=9) | n/a (TNT exposes no scoring loop) | — | **CLOSED** |
| **TBR (keystone)** | ✅ kernel at-limit (T-P5l) · precompute lazy/incremental-VIEW dead by measurement (M46/T-P5j) · scaffolding below-floor (T-P5m) · re-survey de-opaqued, no opaque-bucket prize (T-P5o) | ✅ done (T-P5h3/h4): **quality gapB=0, efficiency≈1** (the "2–4× candidates" was a counting artifact); residual = per-candidate **throughput 1.3–2.3×** on heavy multistate only | getenv hoist **banked, ~20–26 % MISSION wall** (T-P5n/T-S6d) | **CLOSED — incl. MIDDLE-LEVEL algorithm** (T-P5p, 21-agent audit): TS **already implements** quick-TBR's incremental-length method — the directional up-pass IS the slide at the one-combine-per-node floor; `ts_rate` flat-in-N proves the t² asymptotic; residual ~2× = accepted constant factor (mechanism unpinned — reduce + combine each at-limit). **Closed because the cross-cutting kernel (≈½ EW CPU, 96% of sectorial wall is `tbr_search`) is at-limit — NOT because the phase is <8%.** T-P5n/T-P5o "contradiction" = labeling mismatch → incremental-length is **dead-by-solid-argument**. (c) bound-then-verify now **SETTLED dead-by-proof-plus-magnitude (T-P5q, #51)**; only (d) fused edge-set remains (refutable, ~2-6% wash, low-priority, flagged-for-human). |
| **sectorial (xss/rss/css)** | ✅ Round 6 (T-S6a–d): ~96 % of isolated sectorial wall = the inner `tbr_search`; sectorial-specific scaffolding ≤2 %; byte-identical micro-levers ~2.8 % banked (T-S6c) | ✅ **probe-closed (T-S6e), branch `sect-profile-da0f203f` FULLY MERGED at `00967d77`**: the efficiency axis was probed (suppress-trailing-TBR; without-replacement picks) → **AT-LIMIT for safe/behaviour-neutral wins**; the one real lever (consolidate the 3× sequential trailing TBRs in xss→rss→css) is a recipe redesign **handed to #40**. CAVEAT: probe-verdict, **not** a literal TNT-`sectsch` head-to-head race (accepted — other agent's domain, at-limit-by-inheritance). | sector-resolve at parity (#24 float-HTU gated out at plateau) | **BOTH gates closed (sectorial agent); residual lever → #40** |
| **ratchet** | ~AT-LIMIT by inheritance (it is reweight + `tbr_search`; throughput rides the now-closed TBR kernel) | ✅ **DONE (2026-06-21, job `17533025`)**: TS `ts_ratchet_search` vs TNT `ratchet=iter 30` from a shared Wagner start, seeds 1–5 — **cycle-quality PARITY** (same score @ fixed iters: Zanol 1262=1262, Zhu 625=625, Giles 670=670) ⇒ **TNT does NOT reach the optimum in fewer reweight cycles**; wall ~1.8–2.6× = at-limit throughput, no ratchet-specific lever. (Examined-candidate efficiency unmeasured — `RatchetResult` lacks the counter; score+wall are valid.) | `ratchetCycles` 12→6 **banked**, ~20–38 % wall, no quality loss (T-P5d) | **BOTH gates done; recipe lever banked; ratchet = at-limit + cycle-parity** |
| **fuse / drift** | ✅ AT-LIMIT-by-inheritance (#52: `tree_fuse`/`drift_search` both wrap the closed `tbr_search` kernel; drift = the TBR kernel duplicated + an accept rule; no getenv/hidden alloc) | ⏳ low-priority Hamilton-confirmatory (throughput inherited from the closed kernel; fuse race intrinsically awkward — needs a diverse POOL not a single shared start) | fuse >64-tip reroot **crash fix LANDED** (ac8e808a) | **gate-1 done; gate-2 confirmatory/low-priority** |
| **connective tissue (phase 0)** | ✅ AT-LIMIT (read 2026-06-20): in production (default `verbosity=1L`) the per-phase `score_tree` prints are `verbosity>=2`-gated (OFF); only un-gated full rescores are `score_before_cycle`+`score_after_cycle` for the convergence/reset check = **2/outer-cycle**, ~µs each over ~1–few cycles/replicate ⇒ **~0.001% of wall** (one is redundant — `score_before_cycle`≡prior cycle's `score_after_cycle` — but sub-floor, not worth the convergence-logic risk). R.dll 12% already T-P5o'd as amortizable/startup-inflated. Step-switching: each phase owns its state; only orchestrator rebuild = intra-fuse `build_postorder+reset_states` (preset-only, 1/cycle). | — | — | **CLOSED — no addressable production fat** |

### Mission KPI re-measure (2026-06-21) — REFRAMES "the gap" (see dev/profiling/kpi-2026-06-21.md)

Fresh Hamilton run on post-fix cpp-search `5ee3ba3c` (getenv hoist + sector levers
+ ratchet 12→6; freshness-asserted). Two robust conclusions + one corrected
overreach:

1. **QUALITY CLOSED, BANKED (budget-independent).** TS reaches the optimum on
   every dataset/seed; TNT's fast configs miss by +1; on Zanol (ns=9) **TS is the
   only reliably-1261 config (3/3)** — TS is the *more reliable* engine on hard
   data. This is the solid half of parity.
2. **The wall gap is NOT algorithmic.** Candidate-efficiency (COUNT-based,
   throughput-independent; `headtohead_phase0.csv`) is `cand_ratio` ≈ 1.2–1.9×
   near-parity; per-candidate throughput ≈ 2× at-limit. The KPI's eye-catching
   8–110× is a **default-budget mismatch** (TS `default` = heavy search; TNT
   `xmult` default = light), not inefficiency.
3. **Composition #40 is a HYPOTHESIS, not an order-of-magnitude prize** (advisor
   correction to my first write-up): the ratio is biggest where wall is cheapest
   (Wortley/Giles, seconds); on Zanol — where wall actually hurts — the
   thoroughness is **load-bearing** for the reliable optimum. Proven head-room is
   only "thorough→default" (same score, ~2× wall = pure waste); whether there is
   more *below* default without losing reliability is exactly #40's open question.
   Opening diagnostic dispatched: fresh converge-mode h2h (gapB=0 + current
   `cand_ratio`, job `17533024`) + the queued ratchet probe (job `17533025`,
   coarse — units/work-per-iter confounded, order-of-magnitude only).

### Structural clarifications (answering the supervising questions, 2026-06-20)

- **"Thin Sectorial" — is sectorial done, or is there a "fat sectorial" to follow?**
  "Thin Sectorial" is just shorthand for the *lean isolation pass* on the ONE
  sectorial component (the same treatment TBR got) — there is **no separate "fat
  sectorial" component** coming. Sectorial = the single component covering all
  three TNT varieties (**XSS / RSS / CSS**). Gate-1 (AT-LIMIT VTune) is done; what
  remains is gate-2 (the TNT `sectsch` race) + the efficiency loose-ends, both
  owned by the sectorial agent. The *heavier* sectorial questions that surfaced
  (consolidating the 3 sequential trailing TBRs in xss→rss→css, sector-size tuning
  — T-S6e) are **RECIPE composition (#40)**, a separate axis, not a second
  sectorial element.
- **Are fuse / rss / etc. covered, or do they have their own slots?**
  - **RSS** (Random Sectorial Search) is **not** a separate component — it is one of
    the three sectorial varieties, covered under "sectorial" (Round 6 instrumented
    `rss_search` directly).
  - **FUSE** (tree fusing) and **DRIFT** each have **their own slot** (component 4,
    "fuse / drift — later") and are **not yet isolated/raced**. Drift had QUALITY
    work (#25 TNT-faithful drift for +1 datasets); fuse has a pending >64-tip
    reroot-crash fix to port ([[fuse-reroot-segfault]]). Both still owe the two
    isolation gates.

### Status: ALL COMPONENTS CLOSED — program complete (2026-06-21)

Every component is through **both** gates and measured at-limit:
- **scoring** ✅ · **TBR** (kernel + precompute + scorer + middle-level algorithm,
  T-P5p) ✅ · **sectorial** ✅ (probe-closed + branch merged `00967d77`) ·
  **ratchet** ✅ (cycle-parity isolated race, 2026-06-21, job `17533025`) ·
  **fuse/drift** ✅ (at-limit-by-inheritance; gate-2 low-priority confirmatory) ·
  **connective tissue** ✅.
- The 2026-06-21 mission KPI (above) confirms the synthesis end-to-end: quality
  ≥ TNT, candidate-efficiency ~1.5× near-parity, throughput ~2× at-limit ⇒ the wall
  gap is **budget/composition**, NOT per-component throughput.

Residual TBR thread = only sub-lever (d) per-candidate fused edge-set (bit-identical,
~2-6% predicted wash, flagged-for-human — not a blocker). lever-c SETTLED dead
(T-P5q/#51). Data-class reopens recorded in T-P5p/T-P5q (large-N/molecular → revive
incremental-VIEW; binary/DNA ns≤4 → scalar scorer / S2 split shift; S1 + S3 lemma
data-independent).

### GATE BEFORE COMPOSITION #40 — fresh-eyes component re-audit (2026-06-21, user-ordered)

Before tuning the recipe: **"what did we MISS in the individual components?"** An
adversarial completeness pass — independent auditors per component, tasked to
*break* the at-limit verdicts (find an untested lever, a wrong assumption, an
uncovered data-class, a measurement blind-spot), NOT to re-confirm them. The getenv
hoist (~20-26% mission wall, VTune-invisible) is the precedent: the biggest win of
the program was something the standard measurement *missed*. Survivors that are
genuinely new feed back into the relevant component; composition #40 begins only once
this pass is dry.

**RESULTS (2026-06-21, workflow `wf_24dc492a`, 27 agents, 8 component lenses):
18 candidates → 3 survived adversarial verification → 15 killed as rediscoveries/
refuted. The core kernel/TBR THROUGHPUT verdicts STAND — no second getenv-class
hidden hotspot.** Confirmed solid-at-limit, nothing new: scoring-kernel, TBR
precompute (incremental-length=quick-TBR already done), ratchet (12→6 banked),
starting-trees/Wagner. The 3 survivors (all MODEST — none a confirmed multi-x win):
- **#55 (rank 1, fuse, HIGH conf) — getenv-class-in-KIND:** all fuse-VALUE evidence
  on the >64t mission datasets predates the 2026-06-20 reroot fix (`ac8e808a`); the
  recipe's "fuse is free / +1 intraFuse regression" rests on pre-fix runs whose
  multi-round path was skipped/truncated/corrupting. Fuse is **unmeasured on correct
  code**. (Default `poolSuboptimal=0` ⇒ pool size 1 ⇒ inter-replicate fuse SKIPPED
  entirely.) Re-measure dispatched (Hamilton job `17533029`): count productive
  `Fuse improved` events with `poolSuboptimal=5`+`intraFuse` on Zanol/Zhu/Giles.
  Binary → either a #40 simplification (drop wasted fuse) or a recovered quality
  lever on the hardest datasets. **Direct #40 input.**
- **#56 (rank 2, sectorial, MED conf) — NEW throughput lever:** `build_reduced_dataset`
  (ts_sector.cpp:431-440) copies the full block structure; `active_mask` is GLOBAL ⇒
  constant-within-sector-but-globally-informative columns scanned at every inner-sector
  node (~96% of sectorial wall). Offline (reproduced): ~40-60% fewer SIMD blocks on
  Zhu/Zanol/Giles (weakest on Zanol). Column analog of the row-only
  `sectorCollapseTarget`. GATED by early-abandonment ⇒ realizable only if front-packing
  cuts blocks-reached-before-bail; correctness needs the HTU pseudo-tip state. Net
  realistic low-single-digit to ~10%.
- **#57 (rank 3, tbr-scaffold, LOW conf, likely sub-floor):** x4 reroot batch scores
  every member to the deepest-bailing member's depth; gross ceiling ~1-2% EW,
  ILP-confounded. Cheap wasted-block counter as a kill-gate.

Net: the at-limit picture holds; the addressable wall stays in orchestration (#40),
and the strongest survivor (fuse) is itself a #40 input. Tasks #55-57; #40 blocked-by
#55,#56.

### AUDIT FOLLOW-UPS — RESOLVED (2026-06-21, all measured on Hamilton)

- **#55 fuse → DROP (dead weight).** Probe `17533029`: fuse FIRES on the >64t mission
  class (7 attempts/run, 60-70 exchanges/run) but **0 improvements** across
  Zanol/Zhu/Giles × pool/intra × 2 seeds. Not pool-collapse — genuinely useless. #40
  input: drop fuse / raise `fuseInterval` for this class (it is already off-by-default
  since `poolSuboptimal=0` ⇒ pool size 1).
- **#56 sectorial column reduction → SHIPPED opt-in (`830b8cc3`, `TS_SECT_COLREDUCE`,
  off by default).** `reduce_sector_columns_ew` drops constant-within-{sector tips+HTU}
  chars (0 Fitch steps ⇒ scores exact) + repacks survivors into fewer n_states-grouped
  blocks. Adversarial review (`wf_3727ea63`) caught a CRITICAL stale-`rd.subtree`-stride
  OOB that the in-process-toggle A/B had FALSE-PASSED (flag read once at static init ⇒
  both arms ran OFF); fixed. Re-validated (`17533059`): dScore=0 9/9, valgrind clean,
  review-verified invariance+bit-arithmetic. **Saving (rss-isolated): Giles 17%, Zhu 9%,
  Zanol ~0%** (uniform ns=9 = least reduction = the load-bearing case). `dCand≠0` on
  mixed-n_states (block reorder shifts bail timing ⇒ equally-optimal different path;
  Zanol uniform = byte-identical) ⇒ **OPT-IN, never a default flip.** Before default-on
  for any class: run a sector-score ORACLE (reduced vs full, same topology, mixed state)
  — an accept-gated search can't discriminate a masked packing bug. **#56 = a #40
  ingredient (enable per-class where it helps; never Zanol), not standalone.**
- **#57 x4 reroot waste → SETTLED: x4-optimal, force-scalar REJECTED.** Counter probe
  (`17533033`) measured X4_WASTE frac=0.137 = ~1.9% EW gross ceiling. The force-scalar
  A/B (`17533065`, runtime flag `TS_REROOT_SCALAR`, separate processes) settles the sign:
  GATE PASSED (dScore=0 **and** dCand=0 9/9 = byte-identical score+candidates), wall
  speedup x4/scalar = Giles 0.939, Zhu 0.945 (scalar **5-6% slower**), Zanol 1.001 (dead
  heat); overall 0.946. ⇒ the ~1.9% ceiling is **not realizable** — the x4 ILP (4
  independent `any_hit_reduce` chains) more than covers it; forfeiting it loses 5-6% on
  mixed-state and breaks even on ns=9. Flag reverted (measurement-only). Closed.

### CROSS-CUTTING LEVER characterised (2026-06-21): `clipOrder=2L` = per-class, Zanol-only safe

The switches reference flagged tips-first clip ordering as an untested cheap throughput
trial. Now measured (`17533071`@20-rep + `17541277`@40-rep, 3 seeds, EW): `clipOrder=2L`
is ~1.25× faster / ~26% fewer candidates overall, but it **biases the search trajectory**
(not byte-identical) and is a **per-class TRADEOFF, not a global win**:
- **Zanol (uniform ns=9): CLEAN win** — 3/3 reach 1261, consistently ~1.5× faster.
- **Zhu (mixed): quality tradeoff** — loses +1 on 1 seed *even at 40 reps* (doubling the
  budget did NOT recover it ⇒ a genuine trajectory effect, not under-budget).
- **Giles (mixed): wall unstable** (one seed examined 60% more candidates).
⇒ #40 may enable `clipOrder=2L` **on Zanol-type data only**; it complements
`TS_SECT_COLREDUCE` (clipOrder helps the uniform-ns case col-reduce can't, and hurts the
mixed-state case col-reduce helps). Default stays `0L`. Recorded in the switches doc §3a.

**Audit follow-ups closed. #40 composition is the next deliberate, supervised move
(gated: recompose-from-scratch on any step-cost change ⇒ all pieces finished first).**

## Why this reframe

The previous round optimised "the expensive phase of the current recipe"
(ratchet, ~60%). That produced a real recipe win — `ratchetCycles` 12→6 is
~20–38% wall with no quality loss (findings T-P5d) — **but recipe tuning only
reshuffles component *proportions*; it cannot address the core belief that TNT
is faster *per iteration*.** The framing decomposition left throughput as a
~1.4–2.3× same-machine residual (32-bit lower bound) but never localised it to a
component.

A winning search combines scoring + TBR + sectorial + ratchet (+ fuse/drift) in
proportions that vary by dataset. We can only responsibly *compose*
"proven-at-limit" components once each has independently been (a) profiled to its
own performance ceiling and (b) raced head-to-head against TNT's equivalent **in
isolation**.

## STEP 0 (BLOCKS the build) — size the prize on 64-bit first

Advisor course-correction: do NOT build the shared-start harness until the
per-iteration gap is pinned on the hardware that counts. "TNT faster per
iteration" = rearrangements/second; `framing.R` already has both rates on 32-bit
(thr 1.36–2.30). The only missing number is **64-bit TNT rate on Hamilton**, and
it needs no new harness — the existing `bench_tnt_headtohead.R` budget mode gives
it. It *sizes the whole program*:
- 64-bit TNT ≈ 1.3× our rate ⇒ gap is efficiency/quality; swing there, the
  component build is NOT worth it.
- 64-bit TNT ≈ 4–9× ⇒ throughput is the prize; the component build is justified.

Step-0 deliverables (Hamilton, 64-bit TNT, EW-fitch, ≥3 seeds, several sizes):
1. **Score at equal wall** (budget mode) — unit-free bottom line: does TNT beat
   TS at equal 64-bit wall, and by how much?
2. **rearr/sec for both** on the same node — the throughput ratio (mind the unit
   caveat below).
3. **Confirm gapB=0 at FULL budget** vs the *same* well-configured TNT — the
   premise that the residual is throughput, not quality. (Ablation's 2/3 Zanol
   misses were weak-budget; re-test at full budget.)

Only if Step 0 shows a large 64-bit throughput gap do we proceed to build.

## Two questions per component

1. **AT-LIMIT?** VTune the component's hot path in isolation. At the
   AVX2 / compiler / memory-bandwidth ceiling, or is there a real optimisation?
2. **vs TNT per-iteration?** Feed the *same starting tree* to both engines; run N
   iterations of ONLY that component; compare:
   - **score reached** — bitness-independent ⇒ correctness/quality of the
     component's neighbourhood (does our TBR reach the same local optimum TNT's
     does, from the same start?).
   - **count examined** — rearrangements / candidates; bitness-independent ⇒
     efficiency (how many moves to get there?).
   - **wall** — 32-bit local TNT = directional + LOWER bound only; 64-bit
     Hamilton = authoritative ⇒ THE "faster per iteration" test.

## Components, isolation entries, metrics

| component | TS isolated entry | TNT isolated invocation | bitness-free metric |
|---|---|---|---|
| scoring (Fitch EW) | `bench_score_micro.R` / `std::chrono` | `length;` (no loop exposed — hard to race fairly) | ns/score; prior AT-LIMIT T-S3b/c |
| **TBR** (keystone) | `ts_tbr_diagnostics(tree=…)` | `tread <start>; bbreak=tbr;` (NO xmult) | score@opt, rearr-to-opt |
| sectorial | ✅ `rss_search` instrumented directly (no export needed) | `sectsch` settings, no ratchet/fuse (see [[tnt-sectorial-recipe]]) | score, #sectors |
| ratchet | ✅ `ts_ratchet_search` (exported, RcppExports.R:135) | ✅ `ratchet=iter N;` via **STDIN pipe** (runfile-arg → curses, fails headless) | score@iters, wall (examined-count N/A — `RatchetResult` lacks it) |
| fuse / drift | later | `tfuse` / drift flags | later |

## Shared-start plumbing (NEW — does not exist in bench_tnt_headtohead.R)

- Build ONE start tree in TS (e.g. Wagner via the existing builder), per
  (dataset, seed).
- TS side: feed via `MaximizeParsimony(tree=…)` or the component diagnostic's
  tree argument (`ts_tbr_diagnostics` already takes a tree).
- TNT side: `tread "(newick);"` then the single-component command. **Verify TNT
  tree-read format + taxon-index↔name mapping matches `WriteTntCharacters`
  ordering** (off-by-one taxon maps would silently invalidate the race).
- Reuse the existing TNT plumbing: `WriteTntCharacters`, alphabetic `.run`
  filename, `iconv(…sub="")`, regex on "Best score:" / "Total rearrangements
  examined:" (bench_tnt_headtohead.R:56–87).

## Ordering

1. **TBR keystone** — shared-start race vs `bbreak=tbr`. Most decisive for
   "faster per iteration", and TBR underlies sectorial + ratchet, so its
   per-candidate cost propagates everywhere. Do FIRST.
2. **Scoring** — largely settled AT-LIMIT (Round 3); confirm the cross-program
   angle only if a fair isolation is feasible.
3. **Sectorial ‖ ratchet** — composition overhead + candidate-selection, on top
   of whatever TBR turns out to cost.
4. **Compose** dataset-size-tailored recipes from the proven-at-limit elements
   (step (x)).

## Caveats / gates

- Local TNT is **32-bit** ⇒ wall directional + lower-bound only; authoritative
  wall race = **Hamilton 64-bit**. Counts + scores are bitness-independent and
  valid locally.
- Commensurability: TNT "Total rearrangements examined" ≈ our `++n_evaluated`
  (confirmed prior, headtohead_phase0). Pin the unit for ratchet "iterations"
  and sectsch "sectors" before racing those.
- **EW-fitch only** — NA/inapplicable path is owned by another agent.
- A finding that a component is AT-LIMIT only counts with a micro-bench; a
  cross-program score/count parity only counts with the *same* start tree
  verified fed to both.
- **Separate throughput from acceptance policy.** TS TBR and TNT `bbreak` differ
  in first-vs-best-improvement, clip order, and accept-equal ⇒ they reach
  *different* local optima from the same start for reasons that are neither bugs
  nor throughput. Report throughput (rearr/sec or rearr-to-fixed-target) and
  quality (*which* optimum) on separate axes; pin TNT `bbreak` settings
  explicitly, never inherit xmult defaults.
- Stochastic components (ratchet, sectorial) → race as seed *distributions*
  (≥3 seeds), not point comparisons.
- Don't race scoring cross-program — TNT exposes no scoring loop. If the TBR
  race shows equal rearrangements but slower TS wall, scoring throughput is
  implicated for free.
