# Drift EW scorer: does the union-of-finals over/under-count matter?

**Companion to `exactness-gate.md`.** That note proved that the deployed
union-of-finals insertion cost (`fitch_indirect_length[_bounded]`, reading the
engine's simplified `tree.final_`) is **not** the exact EW-Fitch insertion cost,
and the main TBR/SPR path was migrated to the exact directional edge set
`E(D)=combine(prelim(D),up(D))` in **2b299e4b**. The **drift** phase
(`src/ts_drift.cpp`) and the standalone **`spr_search`** (`src/ts_search.cpp`)
were never migrated: they still score pure-EW candidates with
`fitch_indirect_length_bounded` and select `argmin`. This note answers whether
that matters, and lands the exact scorer opt-in behind a flag.

## TL;DR

1. **No inexact score is ever leaked.** Every drift accept path re-scores the
   applied tree with `drift_full_rescore` (an exact `reset_states`+`score_tree`),
   so reported/best scores stay exact. Same for `spr_search` (`full_rescore`
   verifies before accept). The over/under-count only distorts *which* candidate
   is selected and *whether the AFD gate admits it* — a heuristic-effectiveness
   bug, not a correctness bug.
2. **But the distortion is large, and biased in the harmful direction.** On the
   deployed path the union error is **two-sided** (both over- and under-counts,
   up to tens of steps) — not the strict over-count the exactness gate measured
   with *fresh* finals. Selecting `argmin(divided_length + u_cost)` is subject to
   the **optimizer's curse**: it preferentially picks candidates whose `u_cost`
   is *below* their true `e_cost`, so the selected move's scan error is biased
   negative. A truly-expensive move can look cheap, pass the AFD gate, and be
   applied outside drift's intended drift envelope.
3. **The exact directional scan is exact here too** (validated:
   `applied_mismatch == 0` on every applied move, reroots included).
4. **Search-quality effect (30-seed Hamilton matched-wall gate, extended to
   nc1024): a real frontier CROSSOVER → land opt-in.** Exact wins the
   `auto`/wall-race regime decisively on all three (~2.7–3× faster, no
   time-to-optimum regression). Union wins the `thorough`/saturated regime on
   2/3 (Zhu2013 clearly, Zanol2014 narrowly) at *genuinely* matched wall — this
   survived the censoring check: extending exact to nc1024 shows it does not
   plateau (it keeps descending, e.g. Zhu 624.80→624.40) but descends *slower*
   than union, which stays ahead at matched wall throughout the tail (Zhu floor
   624.0 vs exact ~624.4). `TS_DRIFT_EXACT` therefore lands **opt-in**, not as the
   default. (`spr_search`/`TS_SPR_EXACT` is the opposite — a large unambiguous
   reach win, no crossover.)

## Instrumentation (`TS_DRIFT_SCANCHK`, cf. `TS_IW_SCANCHK`)

`drift_search` reads two env flags **once** at entry (never per candidate;
`getenv` is ~2.4µs on ucrt):

- `TS_DRIFT_EXACT=1` — route the pure-EW drift scorer to
  `fitch_indirect_length_cached(prelim, &edge_set_buf[below*tw], …)` with a
  per-clip `compute_insertion_edge_sets`, mirroring `tbr_search`
  (`ts_tbr.cpp:1749`+`:1835`). **Opt-in**; union-of-finals stays the default so
  committed behaviour is unchanged until the gate clears it. NA/IW drift paths
  are untouched (gate `!has_na && !use_iw`).
- `TS_DRIFT_SCANCHK=1` — per-candidate two-sided scan-error census plus per-clip
  decision metrics, printed as one `DRIFT-SCANCHK` line per `drift_search`.

`spr_search` has the analogous opt-in flag `TS_SPR_EXACT=1` (no census).

Reproduce: `dev/profiling/drift-exactness-census.R` (output:
`drift-exactness-census.out`). Datasets recoded `-`→`?` (`m[m=="-"]<-"?"`) so
they are pure Fitch/EW — verified `hasGap=FALSE`, and the census's `cands>0`
confirms the EW path ran (a stray inapplicable would route to NA and score 0
candidates). Start = random tree → `ts_tbr_search` to a local optimum → 12
drift cycles, `afdLimit=3`.

## Results (`drift-exactness-census.out`)

Per-candidate error `= union − exact`; `select_flip` = union's chosen move is
exact-suboptimal; `env_violation` = union admits a move whose TRUE cost exceeds
`afdLimit` (drift leaves its envelope); `applied under` = applied moves that are
truly WORSE than the scan claimed (optimizer's-curse signature).

| dataset (EW) | path | over % (max) | under % (max) | select_flip | env_violation | applied under/checked | applied maxdiff |
|---|---|---|---|---|---|---|---|
| Zhu2013 (75t) | union | 45.2% (69) | **49.0% (100)** | 80.4% | **62.3%** | 127/138 | 83 |
| Zhu2013 | exact | 43.0% (39) | 49.7% (86) | 78.4% | 68.3% | **0/360** | **0** |
| Zanol2014 (74t) | union | 49.9% (30) | 44.2% (52) | 78.6% | 52.1% | 128/135 | 42 |
| Zanol2014 | exact | 47.9% (27) | 45.5% (54) | 72.9% | 54.9% | **0/286** | **0** |
| Dikow2009 (88t) | union | 66.2% (15) | 8.6% (16) | 31.6% | 7.6% | 129/170 | 7 |
| Dikow2009 | exact | 65.0% (16) | 10.6% (15) | 30.7% | 11.7% | **0/344** | **0** |

- **Exact path validates the wiring**: `applied_mismatch == 0`, `maxdiff == 0`
  on all three (the exact scan's predicted `best_candidate` equals the
  post-apply `full_rescore` for every applied move — SPR and reroot).
- **Union is two-sided**: e.g. Zhu2013 under-counts 49% of candidates, by up to
  **100 steps**. (The over/under % on the *exact* path measure what union WOULD
  have done along that trajectory — same conclusion, different walk.)
- **Optimizer's curse is explicit**: of the union path's mis-scored *applied*
  moves, 127/128 (Zhu), 128/132 (Zanol), 129/170 (Dikow) are UNDER-counts — the
  selected moves are overwhelmingly the ones the scan was optimistic about.
- **Envelope violations**: on Zhu2013/Zanol2014 the union path admits a move
  whose true cost exceeds `afdLimit` on the *majority* of clips (62%, 52%).

Single-seed final drift scores (NOT a gate; fixed cycles, not matched wall):
Zhu2013 625≡625, Zanol2014 union **1264** vs exact 1266, Dikow2009 union 1609 vs
exact **1606**. Mixed → consistent with "diversifier: expect wash". Note the
interleaved-TBR workload differs hugely (union `tbr_moves` 331/251/45 vs exact
23/26/26): union drift lands in worse places that the exact TBR phase then has
to climb out of.

## Why the exactness gate said "strict over-count" but drift shows under-counts

The gate (`exactness-gate.md`, P2) recomputed **fresh** finals on the clipped
tree, isolating the *simplified-final formula* error (a bounded over-count ≤5).
Deployed drift reads `tree.final_` maintained by `fitch_incremental_uppass`,
which refreshes only the clip-path nodes; off-path finals carry stale content,
so the union scan errs in both directions by much larger margins. **The precise
per-candidate error source (stale finals vs. the union formula on large
multistate trees) is not separately verified here and is not needed:** the
harm follows from selection bias regardless of source, and the exact scan
removes it by construction (`applied_mismatch == 0`).

## Fix scope

- `ts_drift.cpp` — pure-EW SPR (`:505`→`score_ew`) and reroot (`:563`→`score_ew`)
  loops, `edge_set_buf`/`up`/`pre` hoisted to function scope, per-clip
  `compute_insertion_edge_sets`. EW-only; NA/IW unchanged.
- `ts_search.cpp` `spr_search` — pure-EW candidate scan (`:365`), same shape,
  behind `TS_SPR_EXACT`.
- Verified diagnostic-only (unchanged): `ts_rcpp.cpp:794` (returns a
  `match` field), `ts_rcpp.cpp:2660` (timing harness, return discarded).
  `ts_temper.cpp:308` is a legitimate approximate ranker with exact re-score.

## Local directional pilot (`drift-exactness-gate-bench.R`, `-pilot.out`)

8 seeds x nCycles {8,24,64}, same per-seed local-opt start, identical RNG stream
per path (`set.seed` before each drift run) so the ONLY difference is the scorer.
`nCycles` is the replicate-like policy knob; wall is the eval metric
(`policy-in-replicates-not-seconds`). Mean score / mean wall(s):

| dataset | nCyc | exact score | union score | exact wall | union wall |
|---|---|---|---|---|---|
| Zhu2013 | 8 | **625.50** | 626.75 | **0.079** | 0.220 |
| Zhu2013 | 24 | **625.13** | 625.25 | **0.224** | 0.668 |
| Zhu2013 | 64 | 625.13 | **624.38** | **0.555** | 1.698 |
| Zanol2014 | 8 | **1263.13** | 1263.88 | **0.099** | 0.273 |
| Zanol2014 | 24 | **1262.38** | 1262.63 | **0.325** | 0.784 |
| Zanol2014 | 64 | 1262.25 | **1261.88** | **0.823** | 2.244 |
| Dikow2009 | 8 | **1607.38** | 1608.38 | **0.201** | 0.210 |
| Dikow2009 | 24 | **1606.75** | 1607.50 | **0.599** | 0.663 |
| Dikow2009 | 64 | **1606.25** | 1606.63 | 1.615 | 1.668 |

### Matched-WALL frontier (score achievable by wall ≤ W, mean over 8 seeds)

The fair comparison is score-at-matched-**wall**, not matched-cycles (exact is
faster per cycle, so matched-cycles understates it at low wall and overstates it
at high wall). Frontier from the same CSV:

| dataset | W=0.3s | W=0.6s | W=1.0s | W=1.6s | W=2.5s |
|---|---|---|---|---|---|
| Zhu2013 union | 626.71 | 626.75 | 625.25 | **624.62** | **624.38** |
| Zhu2013 exact | **625.12** | **625.12** | **625.12** | 625.12 | 625.12 |
| Zanol2014 union | 1263.71 | 1263.88 | 1262.62 | 1262.62 | **1261.88** |
| Zanol2014 exact | **1263.12** | **1262.38** | **1262.25** | **1262.25** | 1262.25 |
| Dikow2009 union | 1608.38 | 1608.25 | 1607.50 | 1606.88 | 1606.62 |
| Dikow2009 exact | **1607.43** | **1607.25** | **1606.75** | **1606.50** | **1606.25** |

**There is a frontier CROSSOVER, not a uniform win** (advisor-flagged; my initial
matched-cycles read missed it):

1. **Fast/wall-race regime — exact wins on all three.** Exact reaches a
   given score at less wall (it amortizes `compute_insertion_edge_sets`/clip
   over cheap cached scans; the ~2.9× whole-`drift_search` wall gap on Zhu/Zanol
   also reflects union drifting to worse places so the interleaved TBR phases do
   more cleanup — 331 vs 23 `tbr_moves` — i.e. faster *and* shallower are one
   coin, not two wins). This is the `auto` objective: **no time-to-optimum
   regression.**
2. **Saturated/thorough regime — union wins 2/3.** Exact *plateaus* (Zhu 625.12
   from nc≥24; Zanol 1262.25), while union keeps *climbing* (Zhu → 624.38; Zanol
   → 1261.88) at ~3× the wall. This is structural across the whole sweep, not
   nc-64 noise: union's productive wrongness admits falsely-cheap moves →
   over-diversifies → escapes to a deeper optimum the interleaved exact-TBR then
   captures. Exact's cleaner hill-climb converges faster but shallower. Dikow2009
   shows no crossover — exact dominates every budget.

This is precisely the `auto` vs `thorough` split (`auto-vs-thorough-objective`):
exact is the better `auto` scorer; union has a small saturated-reach edge on 2/3
for `thorough`.

### Authoritative Hamilton gate — 30 seeds, nCycles {8..256} (CONFIRMS the crossover)

Run on Hamilton (build job 17870086 → array 17870087, one task/dataset;
`dev/benchmarks/hamilton_drift_exactness_{build,array}.sh`; raw CSVs
`drift-exactness-gate-hamilton-*.csv`; reanalyse with
`drift-exactness-frontier.R`). Matched-**wall** frontier, mean score over 30
seeds:

| dataset | W=0.2s | W=0.4s | W=0.8s | W=1.6s | W=3.2s | W=6.4s | winner path |
|---|---|---|---|---|---|---|---|
| Zhu2013 union | 626.6 | 626.1 | 625.3 | **624.6** | **624.4** | **624.2** | exact<0.8s, union≥1.6s |
| Zhu2013 exact | **626.2** | **625.6** | **625.2** | 624.9 | 624.8 | 624.8 | (plateaus 624.8) |
| Zanol2014 union | 1264.0 | 1263.5 | 1262.7 | 1262.0 | **1261.7** | **1261.5** | exact≤1.6s, union≥3.2s |
| Zanol2014 exact | **1263.3** | **1262.6** | **1262.3** | **1261.9** | 1261.9 | 1261.9 | (plateaus 1261.9) |
| Dikow2009 union | 1608.1 | 1608.1 | 1607.6 | 1607.2 | 1607.1 | 1606.8 | **exact all budgets** |
| Dikow2009 exact | **1608.0** | **1607.8** | **1607.2** | **1606.9** | **1606.7** | **1606.6** | (no crossover) |

- **Exact wins the `auto`/wall-race regime on all three** — and by a large wall
  margin (at nc=256: Zanol exact 3.08s vs union 8.28s; Zhu 2.03s vs 6.27s;
  ~2.7–3× faster per matched cycle). No time-to-optimum regression anywhere.
  This is solid and unaffected by the caveat below.
- **The nc≤256 saturated tail was censoring-confounded — extended to nc1024 to
  settle it.** Both scorers stop at nc256, but exact runs ~3× faster, so its
  nc256 point (~2s) was being compared against union's ~6–8s point. So exact+union
  were re-run at nc {512,1024} on Zhu2013/Zanol2014 (30 seeds;
  `drift-exactness-gate-hamilton-*-ext.csv`; combined reanalysis
  `drift-exactness-frontier-combined.R`). Result: exact does **not** plateau (Zhu
  nc256→512→1024 = 624.80→624.73→624.40; Zanol 1261.87→1261.73→1261.50) — it keeps
  descending, just *slower* than union. At genuinely matched wall (both with real
  points in 3–12s) **union stays ahead throughout the tail**: Zhu W=6.4s union
  624.17 vs exact 624.73, W=12s 624.17 vs 624.40, union floor 624.0 vs exact ~624.4;
  Zanol union ahead by a narrower ~0.13. So the crossover is **real, not an
  artifact** — union genuinely wins the thorough regime on both (clearly on Zhu).
- **best-of-30 is equal on all three** (Zhu 624, Zanol 1261, Dikow 1606): both
  scorers *can* reach the same optimum — a reliability/wall difference, not
  reachability.

### spr_search (`TS_SPR_EXACT`) — validated, and a large win (no crossover)

`spr_search` is a pure hill-climb with no compensating exact-TBR phase. Smoked on
Zhu2013 (5 seeds, run to convergence): the exact-path returned score equals an
independent `ts_fitch_score` of the returned tree on all runs (wiring valid), and
**exact converges far deeper than union** (e.g. 639/631/626/634/630 vs
810/682/717/690/698). The union over-count inflates `best_candidate` past the
`dominated` gate (`ts_search.cpp:379`), so improving SPR moves are skipped
without a full re-score → premature convergence. With nothing to clean it up,
the fix is an unambiguous reach gain for `spr_search` — a stronger case than
drift for eventually flipping its default (its own gate still recommended).

## Gate — CLOSED. Verdict: LAND OPT-IN, do not flip the drift default

Settled by the 30-seed Hamilton gate + the nc1024 censoring check:
- **Exact wins the `auto`/wall-race regime decisively on all three** (no
  time-to-optimum regression, ~2.7–3× faster) — justifies landing `TS_DRIFT_EXACT`
  and, as a follow-up, having the `auto`/fast presets enable it.
- **Union wins the `thorough`/saturated regime on 2/3** (Zhu2013 clearly,
  Zanol2014 narrowly; Dikow no crossover) at genuinely matched wall — confirmed
  real, not a censoring artifact.
Per the rule (`auto-vs-thorough-objective`: thorough disqualifier = saturated-reach
regression), a scorer that regresses saturated reach on 2/3 must not be the global
default. `TS_DRIFT_EXACT` lands **opt-in**; union's productive-wrongness
diversification stays the default for `thorough`.

`spr_search` is the opposite case — a large unambiguous reach win, no crossover
(the union over-count trips the `dominated` gate and prematurely converges, with
no exact-TBR phase to recover). Flipping its default is well-motivated but
warrants its own matched-wall gate; kept opt-in (`TS_SPR_EXACT`) here.

The correctness/severity finding (no score leak; two-sided error + optimizer's
curse) is independent of and unaffected by the gate outcome.
