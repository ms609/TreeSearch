# Chasing TNT's 1.9 ns/candidate: how far the kernel can go (2026-07-14)

> **CORRECTION (2026-07-15, advisor-caught).** The char-packing dismissal below rested on a METRIC
> BUG: the per-char state count was inflated by counting the all-states `?` token as all global
> states. The TRUE per-char alphabet (states a char can actually take; `?` = the char's own states,
> which is Fitch-EXACT) is **mean ~2.6–3.4 states** vs a global 9–10 that bit-slicing pads to. A
> per-char/per-block LOCAL alphabet cuts words touched per candidate by a **corrected median 0.40,
> up to 0.62 headroom** (5432: 0.47) — a REAL, EXACT packing lever, corpus-wide. So packing is NOT
> refuted. Its *speedup* is regime-dependent: ~1.2× compute-bound (fixed per-block overhead
> dominates at low planes) to ~2× cache-bound (words = cache misses; the EPYC/large-buffer regime).
> Stacks with incremental (1.36×) → ~1.6–2.7× total. Details: this file's tail + `mission-b-goloboff-gates.md`
> Gate-2 correction. The "25× is artifact" conclusion stands; the "~1.4× ceiling" does NOT.


Directive: get TS's per-candidate TBR kernel as close to TNT's ~1.9 ns as the algorithm allows,
profound rewrite permitted. This is the empirical result of that chase (baselines
`kernel-speed-baseline.md`; design workflow `kernel-reslice-plan.md`; bail-depth probe below).

## The decisive new measurement — bail-depth distribution (real weighted data, near-optimal)

Env-gated `TS_AUDIT_PROBE` histogram of blocks-scanned-per-candidate before the cutoff bail, in
the scalar per-candidate scorer (`fitch_indirect_length_cached`), on the actual path real
morphological (weighted) data takes:

| dataset | nBlocks | **mean bail** | % of tree scanned | char-granular-bail ceiling |
|---|---|---|---|---|
| project4327 (823c) | 13 | 3.1 blk | 24% | ~16% |
| project2668 (1227c) | 20 | 7.4 blk | 37% | ~7% |
| project970 (1844c) | 29 | **14.5 blk** | **50%** | **~3%** |

**Candidates DEEP-scan** — 24–50% of all blocks before busting cutoff, not the ~1 block the
"shallow bail" model assumed. This is because big-subtree clips have large cutoffs
(`best − divided_length`), so a candidate must accumulate many step-adders before it can be
rejected. It directly explains the ∝-total_words baseline.

## What this refutes (the levers that would approach 1.9 ns), with reasons

1. **Char-packed early-bail — REFUTED for this corpus.** Char-granular bail saves only the tail of
   the *last* block: **3–16%** (the mean-bail reciprocal above), because the candidate scans
   14/7/3 full blocks first regardless. AND the density penalty is adverse: 5–6-state data rounds
   up to an 8-bit field (1.3–1.6× inflation), so char-packed touches *more* words than bit-slice on
   the deep scan. Net ≈ neutral-to-negative. (Wins only for binary / power-of-2-state corpora.)
2. **Block-major layout — likely NET-NEGATIVE for deep scans.** A depth-first deep scan reads one
   node's blocks 0..14 *contiguously* (node-major); block-major scatters those into 14 separate
   block-planes. Only a *breadth-first* block-plane traversal streams well — a real restructure
   whose ceiling is still the compute floor (below).
3. **Union-construct screen — DEFEATED by incongruence.** It is a *sound* lower bound (proven,
   `union-construct-lower-bound.md`, with true finals) and it IS the mechanism behind TNT's low
   ns/examined (bulk-reject far reinsertions of big clips). BUT on incongruent data the union set
   (all descendant states) is large → the bound is loose → it rarely exceeds cutoff → it rarely
   rejects → pure overhead (Goloboff 1996 explicit: incongruent data "did not save any time — or
   even made the program slower"). **All 13 large corpus datasets measured are incongruent
   (CI 0.11–0.35).**

## Why 1.9 ns is not a like-for-like target on incongruent data (the core finding)

- The per-candidate cost of the exact indirect-length method is **O(blocks-until-cutoff)** — a
  compute+cache floor. For project970 the deep scan is ~14 blocks × 6 words ≈ 84 words of SIMD
  reduce+popcount ≈ **~35 ns irreducible compute floor**. No representation change avoids counting
  to the cutoff on a candidate you must fully score.
- **TNT uses the SAME indirect method** (Goloboff 1996). Its ~1.9 ns is a *per-"rearrangement
  examined"* average whose denominator includes the huge cheaply-rejected neighborhood (length-check
  early-outs on far reinsertions + union-construct bulk-rejects). It is **not** a per-*fully-scored*
  kernel number. On incongruent data the screen can't bulk-reject, so **even TNT must deep-scan the
  near-optimal candidates** — its per-fully-scored cost there is tens of ns too, not 1.9.
- Corroboration: the head-to-head shows TS reaches the optimum *faster* than TNT on 5/6 large
  datasets — impossible if TS's per-*useful*-candidate throughput were 25× worse. The 25× is a
  per-examined counting artifact, not a per-useful-candidate deficit.

## LINCHPIN RESOLVED — 5432 measured directly (2026-07-15)

Pulled `project5432.nex` (Hamilton `lgsweep/matrices/`) and measured the ACTUAL target:

- **482 tips, 189 chars, 9 states, 3 blocks, floor 1943, CI = 0.159** → **highly incongruent**
  (min-steps 309 / 1943). Same regime as all 13 others. **⇒ the union-construct screen is DEFEATED
  on 5432** (loose bound, rarely rejects — Goloboff's incongruent-data case). The "if congruent,
  build the screen" branch is closed: it is not congruent.
- **Near-optimal per-candidate = 14.3 ns LOCAL** (SPR 22, REROOT 13.5), bail at **block 2 of 3**
  for 95% of candidates (~18 words, cache-bound gather from the 482-node buffer). NOT 47 ns.
- **The "~47 ns" that motivated the mission is a confound**, already flagged in the harness itself
  (`micro_tbr_k.R` header, prior advisor 2026-07-12): 47–60 ns is the **far-from-optimal** regime
  (huge cutoff, scan never bails); near-optimal is ~14–20 ns. 5432 has only 3 blocks, so it
  physically cannot deep-scan — its cost is the cache-cold gather, not block count.

## The "25×" is largely a measurement artifact (harness-confirmed)

Read the fingerprint harness (`reeval/fp_run.sh` + `micro_tbr_k.R`):
- **TNT** rate = `rearrangements EXAMINED / wall` → 1.4–1.9 ns is **per-examined**, and "examined"
  includes the whole neighborhood cheaply rejected by the length-check + union screen.
- **TS** rate = `wall / candidates_evaluated` → 47 ns is **per-FULLY-SCORED**.
  These denominators are **not commensurable** (KPI already flagged this). TNT counts a
  union-screen bulk-reject of ~17 destinations as 17 "examined" at ~1 node's cost; TS scores each.
- Corroboration: the head-to-head shows TS reaching the optimum **faster** than TNT on 5/6 large
  datasets — impossible if TS's per-*useful*-candidate throughput were 25× worse.

**Net: the 25× conflates (a) far-from-optimal vs near-optimal regime (~47 vs ~14–20 ns), (b)
per-examined vs per-fully-scored denominators, and (c) a screen that incongruent 5432 defeats for
TS AND TNT alike.** The real near-optimal per-fully-scored cost on 5432 is ~14–20 ns, scanning
~2 cache-bound blocks — near the floor for the exact indirect method on incongruent data.

## (Superseded) the conditional path

The one thing that flips this: **is the mission's actual 5432 (482t×189c) congruent or
incongruent?** If congruent, the union-construct screen (sound, proven) works on it and IS the
profound rewrite that reaches TNT-class ns/examined — worth building. If incongruent (as ALL 13
comparable datasets are, and as a hard basin-capture dataset likely is), 1.9 ns/fully-scored is
unreachable and the exact ceiling is ~1.4–1.5×.

Measuring 5432's CI needs the Hamilton neotrans matrix; **Hamilton SSH is currently
`Permission denied (publickey)` = a dropped university VPN session** (worked earlier this session,
no key change). ACTION: restore the VPN, then `ci_fix.R`-style CI on 5432.

## Packing PROTOTYPE — BUILT + MEASURED (2026-07-15, `TS_PACK_LOCAL`)

Implemented per-block-local-alphabet packing in `build_dataset` (`ts_data.cpp`, default-OFF):
each block carries only its chars' union alphabet (not global n_states); `?` → the char's own
alphabet (Fitch-exact); scorer/combine unchanged (index-agnostic, fewer planes). Split the
alphabet-clustering reorder behind `TS_PACK_SORT`. Gate `dev/profiling/reeval/pack_gate.R`.

- **EXACTNESS: PASS** — `TreeLength` of fixed random trees byte-identical pack OFF vs ON
  (40/40 on 5432, 40/40 project970, 30/30 project510). The representation change is score-exact.
- **TRAJECTORY-IDENTICAL control (cost-order-preserved):** a full `ts_tbr_diagnostics` descent OFF
  vs ON is byte-identical in `n_candidates_evaluated` (5432: 297,721,463; project510: 47,892,620)
  AND converged length — so the ~1.2× is a clean *same-candidate-sequence* per-candidate speedup,
  and collapse did not diverge in the descent.
- **SCOPE CAVEAT (must foreground): scorer-exact, standard-Fitch only.** The gate recodes `-`→`?`,
  so only the pure-applicable path ran. The `has_inapplicable` NA/BGS branch (NA at plane 0,
  applicable at `loc[s]+1`) and IW/HSJ/XFORM have NOT executed. This is a **scorer-exact
  PROTOTYPE, not a deployable exact win.** Deploying needs the per-block plane→global map wired
  through collapse / MPT / reconstruction / output + a full `MaximizeParsimony` byte-identity gate
  on EW + IW + NA. All numbers are LOCAL (EPYC would firm the packing figure as it did the baseline).
- **SPEED — PACK with cost-order PRESERVED (no alphabet reorder), near-optimal, EXACT, scales with
  n_chars:**
  - 5432 (189c, buffer ~207 KB): **1.08×**
  - project970 (1844c): **1.21×**
  - project510 (2954c, buffer ~1.4 MB): **1.23×** (SPR 1.15, REROOT 1.25)
- **The alphabet-clustering reorder is COUNTERPRODUCTIVE — do NOT do it:** 5432 drops 1.08→1.05×,
  project510 to ~0.95×. Clustering by alphabet unlocks more union-headroom but reorders away from
  cost-ordering → deeper early-bail → the lost bail exceeds the extra word-saving. The earlier
  "wash" reading was entirely this reorder; the right design keeps cost-ordering and takes the
  partial per-block headroom that cost-ordered blocks already have (most 64-char blocks union to
  < the global 9–12 states because morphological chars are ~2–3 state).
- **Why it scales with n_chars:** more chars ⇒ more blocks ⇒ the padding saved (global − local
  planes per block) accumulates, and the larger edge_set buffer makes the shrunk gather matter more.
  So packing pays most exactly where the kernel is slowest (char-rich large matrices).

**Net corrected ceiling:** packing 1.08–1.23× (exact, cost-order-preserved, scaling with n_chars)
× incremental precompute 1.36× ≈ **~1.47× (5432) to ~1.65× (char-rich)** per-candidate, EXACT. The
user's instinct that the kernel has real structural slack was RIGHT — the global-alphabet padding
is genuine waste (my "refuted" was a metric bug), and cost-order-preserved packing monetises a
meaningful, exact fraction of it. It does NOT reach 1.9 ns (that remains a per-examined /
regime / incongruent-screen artifact), but it is a real multi-percent-to-20%+ per-candidate win.

## Behavioral validation (the RIGHT bar — reach + MPT-spread + speed, not byte-identity)

Byte-identity is not required (user, 2026-07-15): the goal is fast-to-best-score, good MPT spread,
reach preserved corpus-wide. Wagner is NOT broken — score-exact per tree means it never
miscomputes; under packing it breaks equal-cost insertion **ties** differently (plane-dependent),
taking a different *valid* path. Worst-case dataset (Dikow2009, the biggest byte-identity
divergence), 8 seeds, clean probe-free build:

| metric | OFF | ON |
|---|---|---|
| reach (min score) | 1606 | **1606** (8/8 seeds) |
| mean score | 1606.00 | **1606.00** (Δ 0) |
| MPT spread (mean trees) | 100 | **100** |
| total wall | 1401 s | **926 s (0.66× → ~1.5× faster)** |

Reach + spread preserved; ~1.5× faster (per-candidate ~1.2× × a favorable trajectory realignment
on this dataset). Corpus-wide aggregate in progress (Hamilton array `17874774`: 23 matrices ntax
30–482 / states 2–16 + char-rich + 5432, OFF vs ON × 5 seeds) to confirm the aggregate is
reach-neutral + faster (a few unfavorable realignments tolerated if the net is a win).

## RULING (2026-07-16) — the lever is confirmed; the default-flip has one async gate left

Split the question (per advisor): **is the lever real/exact/sound? YES — rule now. Is it safe to
flip default-ON unconditionally? Not quite — one confirmatory measurement, in flight.**

**Corpus campaign `17874774` completed 12 of 23 cells** (all ≤114 tips; the 11 keys >114 tips —
970, 4085, 510, 2668, 4327, 3806, syab07204/05, 4550, 2183, **5432** — hit the 8h walltime and
produced no output, because a cell ran 5 seeds × 2 modes = 10 full `MaximizeParsimony` searches).
The 12 completed (5 seeds each, OFF vs ON):

| key | nTip | nChar | reach OFF→ON | mean-score Δ | MPT OFF→ON (mean) | speedup |
|---|---|---|---|---|---|---|
| project3210 | 37 | 70 | 335 = 335 | 0 | 49.4 → 49.4 | 1.10× |
| project4363 | 36 | 76 | 234 = 234 | 0 | 1.0 → 1.0 | 1.13× |
| project2798 | 76 | 92 | 439 = 439 | 0 | 21.4 → 21.4 | 1.02× |
| project4044 | 30 | 93 | 177 = 177 | 0 | 8.0 → 8.0 | 1.00× |
| project4461 | 44 | 95 | 336 = 336 | 0 | 47.4 → 47.4 | 1.01× |
| project2769 | 102 | 219 | 2039 = 2039 | 0 | 68.8 → 67.0 | 1.06× |
| project2648 | 95 | 272 | 1242 = 1242 | 0 | 62.6 → 66.2 | 1.08× |
| project3656 | 61 | 339 | 1216 = 1216 | 0 | 14.0 → 14.0 | 1.11× |
| project4271 | 105 | 410 | 1727 = 1727 | 0 | 41.8 → 52.2 | 1.19× |
| project198 | 83 | 412 | 1281 = 1281 | 0 | 45.8 → **34.0** | 1.12× |
| project691 | 103 | 446 | 2164 = 2164 | 0 | 76.6 → 89.0 | 1.20× |
| project2099 | 114 | 555 | 2972 = 2972 | 0 | 17.8 → 18.2 | **1.33×** |

- **Reach: 12/12 equal, and *mean* score identical to the decimal on every key (Δ = 0.00 across
  60 runs).** Mean-equality (not just min-reach) is strong evidence ON's Wagner-start draw is
  statistically indistinguishable from OFF's — the perturbation is a different valid draw, not a
  biased-worse one (reseeding already redraws every replicate). No reach regression anywhere measured.
- **Speed: per-key median 1.10× / mean 1.11×; pooled (wall-weighted) 1.20×**, climbing monotonically
  with nChar (1.00× @70c → 1.33× @555c). Consistent with the per-candidate proof (1.2–1.23× on the
  char-rich keys) — the whole-search wall converts the per-candidate win at ~60–100% efficiency.
- **MPT spread: net-neutral *reshuffle*, not loss.** project4271 +25%, project691 +16%, project2648
  +6% UP; project198 −26% DOWN (the lone >10% drop). Several keys sit at the pool cap (691 = 100
  both), so raw ntree is a weak completeness proxy. Within the stated "a couple of matrices"
  tolerance; logged, not a gate.

**Confirmatory tail re-run submitted (2026-07-16)** — the 8h timeout ate exactly the >114-tip
regime that discriminates (harder landscapes, char-rich, budget-sensitive). One
`MaximizeParsimony` **per array task** (≈10× headroom vs the dead 10-run cells), 3 seeds × 2 modes:
- `packtailA` (job **17887590**, `shared`/8h, 60 tasks): the 10 non-5432 tail keys.
- `packtailB` (job **17887591**, 36h, 6 tasks): project5432 (482t) isolated. [Submitted to `long`
  in error — a 36h walltime belongs on `shared` (72h max); `long` is >3-day only. Left running since
  it had already claimed nodes; script corrected to `-p shared`. A single 5432 run is ~8–12h.]
- Decision-critical subset carrying the ruling: **5432 + 970 + 510 + 2668**.

**Tail round-1 result (2026-07-16): 2 char-rich keys landed and CONFIRM the ruling; 46/60 tasks
timed out at 8h (walltime too short again) → resubmitted at 48h.** The two that finished are the
strongest evidence yet, because they are the char-rich big-tip regime the small-corpus set could
not reach:

| key | nTip | nChar | reach OFF→ON (3 seeds) | MPT OFF→ON | pooled speedup |
|---|---|---|---|---|---|
| project4085 | 164 | 716 | 4051 = 4051 (3/3, Δ 0) | 52/41/49 → 39/74/74 (reshuffle) | **1.27×** |
| project510 | 188 | 2954 | 18345 = 18345 (3/3, Δ 0) | 100/100/100 → 100/100/100 (cap) | **1.49×** |

- **Reach: 6/6 paired runs equal, Δ = 0** — extends the reach-neutral finding into the char-rich
  regime. **project510 is the most char-rich key in the corpus and one of the 4 decision-critical
  keys; it now reads reach-neutral + 1.49× faster.**
- **Speed is *larger* here than the small end** (pooled 1.20× on ≤114t → 1.27× / 1.49× on 716c /
  2954c), confirming the "climbs with nChar" law and that whole-search wall converts (and compounds)
  the per-candidate 1.2–1.23×.
- **Weak corroboration ON is faster:** project970 s1 and project4550 s2 finished their ON run but
  their OFF run timed out — i.e. ON crossed the 8h line where OFF did not (no reach comparison, but
  directionally consistent).

**PROCESS LESSON (3rd timeout — recorded).** `MaximizeParsimony` wall is driven by search
*difficulty* (replicate/ratchet count until the stopping rule), NOT predictable from tips×chars:
project970 (1844c/157t) timed out at 8h while project510 (2954c/188t) finished in ~4h. Tuning
walltime tight repeatedly backfires. **Fix:** resubmit `packtailA` on `shared` at **`-t 48:00:00`**
(shared max 72h, so generous provisioning is free) with an **idempotent skip-if-CSV-exists guard**
in the cell (job **17895591**) — the 14 finished cells skip instantly, the 46 timed-out cells
re-run. Provision full-search benchmark tasks generously by default; see
[[hamilton-partition-choice]].

**Pre-committed decision rule:** tail reach net-neutral (a miss or two either way tolerated as
seed noise) → **flip `TS_PACK_LOCAL` default-ON**; a *systematic* tail reach loss → **keep opt-in**
(land the code behind the flag regardless — the lever is confirmed). The `ts_rcpp.cpp:488`
`state_str` plane→global relabel must be wired (or gated to un-packed) before packed
reconstruction is exposed, independent of the default choice.

## Phase-2 deployment gate — the byte-identity view (packing is score-exact, tie-realigning)

Full-pipeline `MaximizeParsimony` C-vs-A (pack OFF vs ON, fixed seed, nThreads=1;
`dev/profiling/reeval/pack_mpt_gate.R`):
- **Score-exact per tree, confirmed everywhere:** `TreeLength` byte-identical 30–40/40 on every
  dataset; a single `ts_tbr_diagnostics` descent is trajectory-identical (n_evaluated + length
  byte-identical: Dikow/Conrad/Zhu/5432/510).
- **Sansom2010 (23t): EW + IW + NA all BYTE-IDENTICAL** through the full pipeline (incl. a 66-tree
  MPT enumeration and the never-before-run NA/BGS packing branch). So packing's data path is
  correct across modes.
- **Dikow2009 (88t): EW + IW DIFFER** — same best score with enough search (6 reps → 1606 both),
  but the MPT set and candidate trajectory diverge (ntree 100 vs 52), and on a *single* replicate
  ON lands on a worse local optimum (1608 vs 1606).
- **Localized:** the divergence is **within a single replicate** (maxReplicates=1 diverges), **not
  collapse** (`collapse=FALSE` identical), **not ratchet** (ratchet0 diverges), and **not the
  descent** (fixed-start descent is byte-identical). ⇒ the **Wagner greedy start** builds a
  different tree under the per-block-local planes (a plane-dependent choice/tie-break in the
  addition scorer). Score-invariant per tree, but it changes *which* trees the stochastic search
  lands on — budget-sensitive (worse/fewer trees on limited effort on hard data).

**Consequence:** packing is a **score-exact-per-tree, search-trajectory-perturbing** change — the
same CLASS as the deployed cost-ordering lever ([[char-cost-ordering]]), NOT a transparent speedup.
Landing it needs ONE of:
- **(A) byte-identity fix** — make the Wagner start's plane-dependent choice canonical (global-state
  / pattern-index order), so packing becomes a pure speedup. Localized to the Wagner addition path.
- **(B) reach/completeness validation** — treat it as a cost-ordering-class lever and confirm
  reach + MPT-completeness neutrality across many seeds/corpus on Hamilton (the bar cost-ordering
  met). The single-replicate worse-score on Dikow is a yellow flag that must be cleared here.

The known-safe R-facing consumer that WOULD mislabel under packing is the per-pattern state
inspector `ts_rcpp.cpp:488` (`state_str` maps plane→global label) — gate it to the un-packed build
or wire the per-block plane→global map before exposing packed reconstruction.

## Achievable now (exact, banked)

- **Incremental edge-set (lever-6): ~1.36× raw TBR** (built, `claude/lever6-edgeset`), attacks the
  ~17 ns precompute third; independent of scan depth. The single biggest exact per-candidate win.
- **+ block-major/breadth-first streaming reroot**: exact, targets the deep-scan gather cache cost;
  unmeasured, ceiling ~1.5× (compute-floor-bound), needs EPYC A/B.
- Stacked exact ceiling ≈ **1.4–1.5×** (51→~34 ns), bounded by the ~1.5× whole-edge-set-machinery
  limit and the ~35 ns compute floor of the deep reduce.

## Provenance

Probe: `src/ts_fitch.cpp` `TS_AUDIT_PROBE` bail-depth histogram (record_bail + dump_bail_hist);
built `.agent-probe` with `PKG_CPPFLAGS=-DTS_AUDIT_PROBE`. Driver `.agent-probe/kern_probe.R`;
raw histograms `.agent-probe/project*.err`. Baselines `dev/profiling/reeval/kern_real.R`.
