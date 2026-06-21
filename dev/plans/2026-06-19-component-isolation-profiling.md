# Component-isolation profiling program (2026-06-19)

## STATUS (updated 2026-06-20)

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

### Next task (this program)
TBR is now CLOSED on **all** gates including the middle-level algorithm (T-P5p
verdict landed: TS already runs quick-TBR's incremental-length method; no buildable
algorithmic lever on this data class). **Critically, "redirect from TBR" does NOT
mean escaping the TBR kernel** — it is cross-cutting (≈½ EW CPU; 96 % of sectorial
wall is `tbr_search`) and now proven at-limit, so what is left to optimise is the
**ORCHESTRATION around it** (recipe / candidate-selection / the T-S6e redundant
trailing TBRs), i.e. **composition #40** — not more per-component throughput (which
is why sectorial also reads at-limit). Sectorial is owned by the sectorial agent;
ratchet is TBR-inherited and paired with it under #39. So the **cleanest unowned
element-isolation slot is fuse / drift (component 4)** — currently untouched (fuse
also has a pending >64-tip reroot-crash fix to port, [[fuse-reroot-segfault]]) —
ahead of **recipe composition (#40)**. **lever-c (bound-then-verify) is now SETTLED
dead-by-proof-plus-magnitude (T-P5q / #51, proof artifact
`dev/red-team/proofs/lever-c-bound-then-verify.md`)** — so there is no longer any
open TBR-side *algorithmic* item; the only residual TBR thread is sub-lever (d)
per-candidate fused edge-set (bit-identical, refutable-not-proven, ~2-6% predicted
wash — flagged for human, not a blocker). Data-class reopen conditions recorded in
T-P5p/T-P5q (large-N/molecular → revive incremental-VIEW; binary/DNA ns≤4 → S2 split
shifts / scalar scorer; S1 + the S3 lemma are data-independent).

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
| sectorial | xss/rss/css R entry — **TBD (may need export)** | `sectsch` settings, no ratchet/fuse (see [[tnt-sectorial-recipe]]) | score, #sectors |
| ratchet | `ratchet_search` R entry — **TBD (may need export)** | ratchet-only command — **syntax TBD** | score, #cycles |
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
