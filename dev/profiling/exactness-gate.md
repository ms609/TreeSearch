# Exactness gate: is fresh rooted union-of-finals an EXACT Fitch insertion cost?

**Question (gates a perf refactor).** TS's TBR reconnection scoring pays an
expensive *directional edge-set* precompute (`compute_insertion_edge_sets`,
`src/ts_fitch.cpp:521`) and then, per candidate edge (A,D), counts characters
where the clipped subtree's prelim set is disjoint from `edge_set[D]`. A cheaper
candidate is **union-of-finals**: cost at edge (A,D) = #chars where
`clipRootStates ∩ (final(A) | final(D)) = ∅` (`fitch_indirect_length`,
`src/ts_fitch.cpp:397`). Prior work abandoned union-of-finals after it
mis-counted, but possibly only because finals were stale or the framework was
unrooted. **Is a CLEAN, freshly-recomputed ROOTED union-of-finals insertion cost
EXACT, or only a bound?**

## Verdict (three-way)

Tested empirically against full-rescore ground truth (`TreeSearch::TreeLength`)
over **6066 reconnection candidates** spanning binary/multistate, missing-data,
clip sizes 1–10, and multiple rootings. Fresh Fitch passes recomputed on the
clipped tree throughout — never incremental/stale.

| Predictor | What it is | Result | Verdict |
|---|---|---|---|
| **P1** `edge_set[D] = combine(prelim[D], up[D])` | engine's **exact** path, per-child directional edge set | signed error `added − pred = 0` for **6066 / 6066** | **EXACT** (rooting-invariant) |
| **P3** union of **true** edge sets `E(A) \| E(D)` | union-of-finals using the *true* directional finals | signed error ∈ **{0,+1,+2,+3,+4}**, never < 0 | **SOUND LOWER BOUND** (under-counts; never exact) |
| **P2** union of the engine's **stored** `final_` | the **deployed** `fitch_indirect_length` reads `tree.final_` | signed error ∈ **{−5,−4,−3,−2,−1,0}**, never > 0 | **UNSAFE — strict OVER-count** (upper bound) |

**Bottom line: a cheap "root like TNT → exact union-of-finals" refactor is a
MIRAGE.** Union-of-finals is *never* exact for EW-Fitch insertion cost — not even
on binary unambiguous data, and not under any rooting. At best it is a sound
*lower* bound, and only when computed over the *true* directional edge sets
`combine(prelim, up)` — which cost exactly the same directional up-pass the
refactor was trying to avoid. The only exact insertion cost is the per-child
edge set `E(D)`, which the engine already computes.

Worse, the union variant that is actually **deployed** today
(`fitch_indirect_length`, which reads the incrementally-maintained
`tree.final_`) is **not even a sound bound: it OVER-counts** the true insertion
cost (by up to +5 here), because `tree.final_` is a *simplified* final set.

## Why the deployed union over-counts (mechanism, not staleness)

The stored `final_` (see `uppass_node`, `src/ts_fitch.cpp:41-70`) uses a
**simplified** rule:

```
final(root)  = prelim(root)
final(node)  = final(parent) ∩ prelim(node)   if non-empty
             = prelim(node)                    otherwise      // NEVER unions
```

Because it never performs Fitch's *union* step, `final_(node) ⊆ prelim(node)`
always — it is generally a **strict subset** of the true MPR/edge-set final.
Union of two too-small sets is too small ⇒ more characters look "disjoint" ⇒
over-count. Concrete minimal case found in the sweep (binary, fresh finals):

```
char: clipRootStates X = {1},  prelim(D) = {0}
      true edge set  E(D) = combine(prelim(D), up(D)) = {0,1}   (up-pass unions in state 1)
        -> X ∩ E(D) = {1} != empty  ->  TRUE added step = 0
      deployed simplified finals: final_(A) = {0}, final_(D) = {0}, union = {0}
        -> X ∩ {0} = empty  ->  DEPLOYED union predicts step = 1   (OVER-COUNT)
```

Finals here are freshly recomputed, so this is a **formula limitation of the
simplified up-pass**, not a staleness artifact. The theoretical "union
under-counts" argument (a union of finals is a superset of the edge set, so it
only *adds* available states) is correct **only for true finals** (that is what
P3 confirms); it does **not** hold for the simplified `final_` the engine stores
and the union screen reads.

### Reconciliation with the code comment
`fitch_indirect_length`'s comment claims the union "UNDER-counts the true
insertion cost (never over-counts)." That is **false for the deployed code**: it
holds for a union of *true* edge sets (P3), but `tree.final_` is the simplified
final, and the deployed union over-counts (P2). Any caller that treats
`fitch_indirect_length` as an *admissible lower bound* (e.g. to prune/skip a
candidate whose bound ≥ current best) is **unsound** — it can discard the true
optimum. It is only safe in its documented use as an *approximate ranker*
followed by exact re-scoring, and even then it can mis-rank (demote a genuinely
good reconnection).

## Equality-on-binary (the target-dataset question)

project5432 is mostly binary, so the key sub-question was: does union == exact on
unambiguous binary characters (where one might hope finals are singletons)?
**No.** Internal Fitch sets are frequently non-singleton even with unambiguous
binary tips, and both union variants fail on binary:

- Binary regimes (2896 candidates): **P1 exact 100%**; P2 (deployed) exact only
  65.2% (**over-counts 34.8%**); P3 exact 72.4% (**under-counts 27.6%**).

## Effect of conditions

Exactness of **P1 is unconditional** (100% across every cell). The union
predictors are never exact and degrade with alphabet size:

| regime | candidates | P1 exact | P2 (deployed) exact / over | P3 exact / under |
|---|---|---|---|---|
| binary, no missing | 1410 | 100% | 68.5% / over ≤3 | 70.8% / under ≤4 |
| binary, 15% missing | 1486 | 100% | 62.0% / over ≤4 | 73.9% / under ≤3 |
| 4-state, no missing | 1564 | 100% | 51.8% / over ≤5 | 60.6% / under ≤4 |
| 4-state, 20% missing | 1606 | 100% | 42.8% / over ≤5 | 67.1% / under ≤3 |

- **Rooting:** `E(D)` is a property of the edge (bipartition), so P1 is
  rooting-invariant — verified exact under 3–4 distinct rootings of the same
  clipped tree across 8 further trees. The union predictors are rooting-*dependent*
  (magnitudes shift), but P2 still over-counts and P3 still under-counts under
  every rooting tried.
- **Clip size / leaf-vs-clade:** single-leaf (SPR) and multi-tip clade (TBR)
  reinsertions both included; no effect on the verdict.
- **Multifurcations:** out of scope — the engine reconnects onto strictly binary
  trees during search; the exact-vs-bound question is about binary base trees.

## Method & reproducibility

Pure R, independent of the C++ engine (guarantees genuinely fresh finals and an
independent oracle). Files (this dir): `exactness-gate.R` (self-contained
harness). Pipeline:

1. **Ground truth = full re-score.** For rooted tree T, clip clade S at node s →
   base T′. For each candidate edge (parent(D), D) in T′, build the reconnected
   tree by explicit edge-matrix surgery (`regraft`) and score with
   `TreeSearch::TreeLength`. `added(e) = L(reconnected) − L(T′) − L_within(S)`.
2. **Predictions from FRESH passes on T′.** `prelim` = my Fitch down-pass;
   `E = combine(prelim, up)` with `up[D] = combine(up[parent], prelim[sibling])`
   (root child ⇒ `prelim[other child]`) — this exactly reproduces
   `compute_insertion_edge_sets`. `Fe` = the engine's simplified `final_` rule.
   `X` = fresh prelim at S's root.
3. Compare `added` to P1 `|X ∩ E(D)|=∅`, P2 `|X ∩ (Fe(A)|Fe(D))|=∅`,
   P3 `|X ∩ (E(A)|E(D))|=∅`.

**Validation of the pipeline itself (so the ground truth is trustworthy):**
- My Fitch down-pass length == `TreeSearch::TreeLength` == `phangorn::fitch`
  (after fixing a tip-label/node-index alignment: `ape::rtree` permutes tip
  labels — the oracle matches by label, so tip masks must be aligned to
  `tree$tip.label`, not to data row order).
- `regraft` well-formed and both scorers (TreeSearch + phangorn) agree on
  **3902/3902** reconnected trees; regrafting at the original clip location
  recovers L(T) on **307/307** checks.
- P1 == full-rescore on all 6066 candidates independently cross-validates the
  edge-set formula, the regraft surgery, and the decomposition simultaneously.

**Caveat on "true MPR final ≠ edge set."** A node's MPR set (states optimal at
the node) is *not* equal to `combine(prelim, up)` in general (the naive "one
extra step" reasoning fails for node states — forcing a cherry's ancestor to the
wrong state can cost +2). But the object relevant to *insertion at an edge* is
the edge set `combine(prelim(D), up(D))`, and inserting a subtree adds ≤ 1
step/char; P1's 100% exactness confirms the edge set (not the node MPR set) is
the right and exact quantity.

## Implication for the refactor

- **Do NOT** replace the exact per-child edge-set scoring with a union-of-finals
  screen expecting exactness — it is a bound at best, under every rooting,
  including on binary data.
- **Do NOT** rely on the deployed `fitch_indirect_length` as an admissible lower
  bound for pruning; it over-counts and can reject the optimum. If a sound lower
  bound is wanted, it must union the *true* edge sets `combine(prelim,up)` (P3) —
  but that requires the same directional up-pass as the exact path, so it buys no
  precompute savings.
- The exact, cheap-per-candidate quantity is the per-child edge set `E(D)`, which
  the engine already precomputes. The genuine perf lever is making that
  precompute cheaper, not swapping in union-of-finals.
