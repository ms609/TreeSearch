# Wagner insertion-cost formula is wrong → +30% starting trees (2026-06-18)

## TL;DR
TreeSearch's RAS Wagner trees score **~+30% over the optimum** on Zanol2014
(mean 1664; optimum ~1261) where TNT's no-swap RAS Wagner scores **~+3%**
(1283–1325). Proven cause: the candidate-edge insertion-cost function
`fitch_indirect_length` (src/ts_fitch.cpp) computes the edge "passing set" as the
**union of the two endpoints' FINAL states**, `Y = final(A) | final(D)`. The union
is a superset of the true edge set → it **undercounts** insertion cost → too many
positions look free → greedy stepwise addition degenerates toward arbitrary
(first-found) placement. The error is largest when state sets are most ambiguous
(early Wagner steps), which is why Wagner is hit so hard.

This affects EVERY search start (RANDOM, GOLOBOFF=the production default, ENTROPY).

## Evidence (all on canonical Zanol2014, EW Fitch, dev/benchmarks/t0/)
- `bench_wagner.R` (K=8): TSrand mean **1656** (sd 66) vs TNT no-swap RAS
  **1300.9** (sd 13); KS D=1.0, p=1.6e-4. Diversity also differs (TSrand
  meanPairwise CID 0.77 vs TNT 0.40 — but that's a *symptom* of near-random
  placement, not a virtue).
- `diag_wagner_verify.R`: kernel's own score == TreeLength(reconstruction) for
  every seed (MATCH) → NOT a reconstruction artifact. The tested public
  `AdditionTree()` path reproduces it (~1515–1667). RandomTree ref = 2295.
- `diag_wagner_exact.R` (**decisive**): an EXACT-insertion RAS Wagner (try every
  edge, full TreeLength, true argmin) reaches **1295–1309** — TNT parity — on the
  SAME addition orders where the fast formula gives 1644–1678 (+342…+370).
  ⇒ the algorithm is fine; the fast cost formula is the bug.
- `diag_wagner_bias_scores.R`: RANDOM 1664 / GOLOBOFF(default) 1661 / ENTROPY 1479
  — bias changes only the order, so all inherit the bug.

## The formula
`fitch_indirect_length(clip_prelim, A, D)` (ts_fitch.cpp:380):
```
Y = final(A) | final(D)            // per character, OR of state words
needs_step = ~any_hit(clip_prelim, Y) & active
extra = popcount(needs_step) * weight
```
Comment claims union is "exact for non-additive (Goloboff 1996); intersection
would overcount." Empirically the opposite: we UNDERcount.

## Proposed fix (to validate against the exact-insertion oracle BEFORE shipping)
Classical Fitch result: the set of states on edge (A,D) in MPRs is
`(final(A) ∩ final(D))` if that intersection is non-empty, else
`(final(A) ∪ final(D))` — i.e. **intersect-else-union**, computed per character.
Adding a tip with downpass set T costs `[ T ∩ E == ∅ ]` per character.

So the candidate fix is to replace the pure union with a per-character
intersect-else-union of the two endpoint finals — the same combine logic the
downpass already uses. Likely localized, but `fitch_indirect_length` has several
siblings that must all change consistently:
- `fitch_indirect_length` / `_bounded` / `_cached` (ts_fitch.cpp)
- `fitch_indirect_bounded_flat` (+ any flat/EW specialisation)
- the NA variant in `ts_fitch_na_incr.h`
- wherever a precomputed `vroot` edge set is built for `_cached` (TBR) — it must
  use intersect-else-union too.

VALIDATION GATE: a fast-formula Wagner must reach the exact-insertion oracle
(~1300 on Zanol) before the fix is accepted. If intersect-else-union of *finals*
doesn't get there, fall back to maintaining a directional uppass view and
combining `prelim[D]` with the incoming view at D.

## Scope / risk
`fitch_indirect_length*` is shared by Wagner, TBR, sector, prune-reinsert, drift,
temper. A correct (tighter) cost estimate should only HELP candidate ranking, but:
- must re-run the full testthat suite,
- must re-benchmark a short EW search (Wagner+TBR multistart) to confirm
  end-to-end improvement and no regression,
- the TBR chip (task: "Compare TBR ensemble: TNT vs TreeSearch") should re-test
  after the fix — the SAME formula drives TBR reinsertion scoring, so TBR may be
  partially degraded too (less than Wagner, since full-tree final sets are less
  ambiguous).

## CORRECTED fix + VALIDATION (2026-06-18, later)
The "intersect-else-union of FINALS" guess above is WRONG (finals are contaminated
by D's own subtree; it gave +150 over oracle). The proven-correct edge set is the
**directional** message combine:
```
down[D] = prelim[D]
up[D]   = combine(up[parent], prelim[sibling])   // root degree-2: up[child]=prelim[other child]
E(A,D)  = (down[D] ∩ up[D]) if non-empty else (down[D] ∪ up[D])   // per character
cost    = #chars where  T & E == 0
```
Reference kernel `ts_wagner_tree_dir` (ts_wagner.cpp) + per-edge probe
`ts_reinsert_scan` (ts_rcpp.cpp) added behind the worktree build.

VALIDATION (strict gate met):
- Per-edge (ts_reinsert_scan, clip+rescore truth): directional == actual **71/71**
  edges (74-tip tree) and **9/9** (12-tip tree); union matches only 4–6/71.
- End-to-end: directional RAS Wagner mean **1308** (Zanol) / **659** (Zhu) ==
  brute oracle (~1300/657) and TNT band (1283–1325); buggy union = 1631/1189.
- Same-order vs brute: diffs ±27, BOTH directions = pure tie-break noise (greedy
  tie-break sensitivity measured at ±15). Tie-break is NOT the lever.
- Speed: directional kernel only **1.7×** slower than buggy (2.0 vs 1.2 ms/tree,
  n=74) — full down+up recompute per step is cheap enough for production.

Two separate defects clarified:
- Cost formula (union→directional): the bug in PRODUCTION wagner_tree.
- Pendant-edge scan: was a bug in my REFERENCE kernel only (used postorder =
  internal nodes only). Production's DFS already scans tip edges → not affected.

NUANCE: on a RESOLVED tree the union ARGMIN is already correct (picks the optimal
edge; only magnitudes wrong). The bug bites during CONSTRUCTION (ambiguous partial
trees). ⇒ TBR/sector (resolved trees) likely far less affected; decide separately
whether their fitch_indirect_length/vroot need the directional set.

## PRODUCTION PORT — DONE (2026-06-18)
Shared helper `compute_insertion_edge_sets` (ts_fitch.cpp/.h) builds the exact
per-node edge set E[D]=combine(prelim[D],up[D]); callers score with
`fitch_indirect_length_cached(clip_prelim, &E[child], ds, cutoff)` (that helper
already computes [T ∩ vroot == 0], so passing E as the vroot gives the exact
directional cost through existing code).

- `wagner_tree`: candidate DFS scan now uses E[below] (constraint filter +
  incremental rescore unchanged).  AdditionTree / random_wagner / biased_wagner
  all fixed: RANDOM 1664→1310, GOLOBOFF(default) 1661→1306, ENTROPY 1479→1304.
- `tbr_search`: EW-only path (`ew_directional = !has_na && !use_iw`) — the SPR
  scan and the rerooting `vroot_cache` now use E[]; NA (three-pass) and
  implied-weights keep union-of-finals (their cached scorers require it).
  vroot_cache[ei] = E[main_edges[ei].second].
- Test `tests/testthat/test-wagner-quality.R`: mean of 8 RAS addition trees
  within 8% of the MPT (Zanol 1261, Zhu 624).  Fixed ~+4–6%; bug was +30%.

VERIFICATION:
- Full testthat: **0 failed expectations**.  6 file-level errors are pre-existing
  (test-CharacterHierarchy.R / test-LeastSquares.R call bare unexported `.fns`,
  invisible under test_dir on an installed pkg — unrelated to this change).  EW
  score checks pass (Vinther2008 TBR/XSS/Ratchet = 79).
- END-TO-END payoff (pure Wagner+TBR multistart, ratchet/sectorial OFF, Zanol):
  | starts | fixed | buggy era | target |
  |--------|-------|-----------|--------|
  | 1      | 1267  | 1315      | 1261   |
  | 5      | 1264  | 1306      | 1261   |
  | 20     | 1264  | 1287      | 1261   |
  One fixed start now beats twenty buggy starts; +3 over the optimum vs +26.

## TODO before merge to cpp-search
- Remove validation scaffolding: debug exports `ts_wagner_tree_dir`,
  `ts_reinsert_scan`, the reference `directional_wagner_tree`, and the
  `TS_WAGNER_UNION` env diagnostic (+ their init.c registrations).
- Optional perf: `compute_insertion_edge_sets` allocates its up[] scratch per
  call; reuse a buffer if a /profile pass flags TBR overhead.
- IW/NA insertion cost still uses union-of-finals (out of scope; separate task if
  it matters for those objectives).

## Status
SHIPPED in worktree (fix + vroot + test, 0 regressions, end-to-end gap +26→+3 on
the core engine). Not yet on cpp-search. Unrooted TBR handled separately by the
chip ([[tbr-rooted-vs-unrooted]]).
