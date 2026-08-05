# T-374b — XFORM rooting policy: decision document

**Date:** 2026-07-29
**Finding:** T-374 (P1, red-team area 10), XFORM half
**Scope:** investigation and recommendation only. No `src/` change was made.
**Branch:** `claude/t374b-xform-rooting-policy` (worktree from local `cpp-search`; rebased onto `45a3a038` before landing)
**Evidence script:** [`dev/red-team/heavy-tests/xform-rooting-oracle.R`](../red-team/heavy-tests/xform-rooting-oracle.R) — pure R, no build, independent re-implementation of the Sankoff DP and the x-transformation cost matrix.

---

## Recommendation, up front

**Option 3: document XFORM as an approximation whose reported score is a rooting-dependent upper bound, canonicalise the rooting at the single place the user-visible discrepancy is produced, and do not pin `sankoff_forced_root`.**

Concretely, and in priority order:

1. **Fix the reported-score/`TreeLength()` discrepancy, which is the actual P1 user-visible defect**, by making the report path and `TreeLength()` agree on one rooting — not by making the objective rooted. This is the lead evidence (Q-E): the quantity the search already optimises is a **valid upper bound** on the well-defined unrooted objective (min over rootings), it is **tight for 87–98% of rootings**, mean overstatement 0.02–0.17 steps, and the worst case is bounded by `nSec` per hierarchy block (Q-B). So agreement at the boundary is cheap and provably close, and nothing about the search needs to change to get it.
2. **Do not set `forced_root_state = 0` as a standalone change.** This was the attractive cheap fix, and it is **incoherent on its own** — not merely suboptimal. Pinning the root *state* makes the root meaningful while leaving the root *position* arbitrary and, at four sites, actively moving (Q3). Measured: the pinned-state criterion's value varies across root positions on **84–117/120** random 9-tip topologies, with spread up to 5 — so the one-liner does not remove root-sensitivity from a rerooting pipeline, it relocates and enlarges it. **Read that number correctly** (Q-C): `forced_root_state = 0` defines a *different, explicitly rooted* criterion, perfectly well-defined at any fixed rooting, so this is not evidence that a rooted criterion is wrong. It is evidence that pinning the state is only meaningful **together with** pinning the position — i.e. it is Option 2, not a one-line fix.
3. **Do not "rethink TBR", and do not run the rooting-pinned A/B.** See Q3 and §"The A/B" — TBR's fragment reroot is not in conflict with a pinned root, the whole-tree rerooting sites are a short enumerable list, and the already-filed **T-377** is a first-order defect upstream that makes any second-order rooting measurement uninterpretable (see the "Blocker" section).
4. **Document XFORM as rooting-sensitive** in `?MaximizeParsimony` and `?RecodeHierarchy`, stating the bound.

The rest of this document answers questions 1–5 and records the blocker (which turns out to be the already-filed T-377).

---

## Blocker found during this investigation: the Sankoff term is absent from XFORM's candidate screen

This is upstream of, and larger than, the rooting question. It is **not** part of T-374 — **it is the already-filed T-377**, "TBR's candidate scan is hierarchy-blind, so moves that improve the HSJ DP or Sankoff term at a cost in Fitch steps are never *proposed*" (P2, verified, `ts_tbr.cpp:2206`). I re-derived it independently from code below before finding that row, and briefly mis-filed it as a new T-383; that duplicate is deleted and T-377 now carries the three refinements this section adds. **T-377 is being worked on concurrently** (worktree `claude/t377-hierarchy-screen`), so coordinate rather than duplicate.

What this section adds beyond T-377's existing row: (a) the convergence sweep's blindness is **`has_na`-conditional**, and `has_na` is data-dependent under XFORM; (b) the `dominated` test is systematically **over-permissive**, a wall-clock cost T-377 does not record; (c) T-377 **gates** T-374's XFORM A/B.

In `src/ts_tbr.cpp`:

- `incremental_ok` (`:1503-1508`) deliberately excludes `ScoringMode::XFORM`, so the **accept gate** `actual = full_rescore(tree, ds)` (`:2717`) is authoritative and includes the Sankoff term. Accepted scores are therefore *correct*.
- But `best_candidate` — which drives both the `dominated` test (`:2640-2642`) and the **argmin over candidates within a clip** — comes from the Fitch/NA indirect scan and contains **no Sankoff term at all**. `score_tree()` (`src/ts_fitch.cpp:1584-1629`) shows the XFORM total is `fitch_score_ew(...) + sankoff_score(...)`; only the first summand reaches the scan.

Consequences: (i) the incumbent `best_score` includes the Sankoff term while the candidate estimate does not, so `dominated` is systematically over-permissive — nearly every candidate looks improving and is applied-then-rescored, which is also a wall-clock cost; (ii) the candidate actually *selected* within each clip is the Fitch argmin, so **XFORM's search is Fitch-guided with a Sankoff correctness filter bolted onto the accept path.**

**The `has_na` question, which sets the severity.** `ew_directional = !has_na && !use_iw` (`:1739`). For XFORM `use_iw` is false (implied weighting is rejected outright, `R/MaximizeParsimony.R:1197-1199`), so the branch turns entirely on `has_na`:

- `has_na` is **data-dependent**, not fixed by the mode. It is true iff some retained pattern has a genuine `-` (`src/ts_data.cpp:152` ← `ts_simplify.cpp:164` `has_genuine_inapp` → `:309` `blk.has_inapplicable` → `ts_tbr.cpp:1489-1492`). Under XFORM, hierarchy characters are zero-weighted and erased, but **non-hierarchy characters keep their gaps** — the `.GapsAsMissing()` recode fires only for `inapplicable = "missing"` (`R/MaximizeParsimony.R:1179-1182`). So a dataset whose gaps are all inside declared hierarchy blocks gives `has_na = FALSE`; a dataset with any undeclared gap-bearing character gives `has_na = TRUE`.
- `has_na = FALSE` → `ew_directional = TRUE` → the convergence sweep uses `try_root_edge_moves`'s fast additive Fitch path (`:700-704`), so the Sankoff term appears **nowhere** in move evaluation, only in the accept gate.
- `has_na = TRUE` → `ew_directional = FALSE` → the sweep routes to `try_root_edge_moves_rescore` (`:704`, `:615`, `:653`) and `exact_verify_sweep` (`:1093`), both of which call `full_rescore` and therefore **do** see the Sankoff term at convergence.

So the honest statement is: **the Sankoff term is never in the inner-loop candidate screen; whether it enters the convergence sweep depends on whether the user's non-hierarchy characters happen to contain gaps.** That is a worse shape than either branch alone, because search behaviour changes qualitatively with an incidental property of the data.

**Why this dominates the rooting decision.** The rooting drift is bounded by `nSec` per block and is an *overstatement of a valid bound* (Q-E below). The screen omission is unbounded in its effect on which topologies are visited. Pinning the rooting would be polishing a mechanism that is blind to the term whose rooting-dependence is the finding.

---

## Q1. Is XFORM's objective intended to be rooted?

**The maths is intrinsically asymmetric; the user's tree is not rooted; and the source paper does not resolve the difference. Answer: asymmetric by design, but not *specified* as rooted, and the engine treats it as unrooted throughout.** These are three separate claims and the code conflates them.

**The matrix is asymmetric by construction.** `R/recode_hierarchy.R:101-118` builds, per hierarchy block with `nSec` secondary characters:

| from → to | cost |
|---|---|
| absent → present | `nSec + 1` (gain) |
| present → absent | `1` (loss) |
| present → present | Hamming distance over secondaries |

Asymmetric whenever `nSec ≥ 1`; **symmetric when `nSec = 0`** (gain = loss = 1). That asymmetry is the whole point of the recoding — it implements the x-transformation of Goloboff, De Laet, Ríos-Tamayo & Szumik (2021), *Cladistics* 37: 596–629 (`inst/REFERENCES.bib:337-346`), in which the first gain of the controlling primary character must pay for the secondaries it brings into existence. The gain/loss ratio is the mechanism, not an incidental detail.

**Asymmetric step matrices make length root-dependent.** This is standard: an asymmetric cost matrix means the minimum-cost labelling total depends on how edges are oriented, hence on where the root sits; PAUP\* therefore requires a rooted tree (and optionally a specified ancestral state) for asymmetric step matrices. The oracle reproduces this directly.

**But "the maths is asymmetric" is not "the user's tree is rooted".** The engine's rooting is an *artefact*, not a hypothesis:

- `MaximizeParsimony` supplies `TreeTools::RandomTree(nTip, root = TRUE)` when no start tree is given (`R/MaximizeParsimony.R:1219`) — the root position is whatever the RNG produced.
- `TreeLength()` scores at whatever rooting the user's `phylo` happens to carry (`R/tree_length.R:198-204`, `:347-354`); it does not canonicalise.
- Wagner addition, sector search, fusing and `ts_collapse_pool` all move the root freely, on the stated premise that length is root-invariant (`ts_tbr.cpp:123-124`, `:2999`, and the already-recorded `ts_collapse_pool` comment).
- `forced_root_state` is `-1` for every block (`R/recode_hierarchy.R:176`), i.e. "min over root states" — a deliberate attempt to *neutralise* the rooting rather than to honour it.

**And the paper does not licence a rooted reading of an arbitrary rooting.** The x-transformation is offered as an *approximation* to inapplicable-aware parsimony. Nothing in it identifies the engine's incidental root with a biological root. Pinning "absent at the root" would be an extra assumption — often defensible for a neomorphic primary character, wrong for a transformational one — and the code has no way to know which the user has.

**Verdict:** the objective is *asymmetric*, therefore mathematically rooted; the *pipeline* is unrooted and has no root to honour. The coherent unrooted reading of an asymmetric criterion is **min over rootings**, and the current score is an upper bound on it (Q-E). That is the reading to document, not a root to invent.

## Q2. If rooted: what rooting, where pinned, what breaks?

For completeness, since the answer is "don't":

- `forced_root_state` is already plumbed end to end and needs no kernel work: `R/recode_hierarchy.R:176` → `src/ts_rcpp.cpp:1815` (`ds.sankoff_forced_root = fr_vec`) → `src/ts_data.h:185` → `src/ts_fitch.cpp:1600` → `src/ts_sankoff.cpp:79-82`. Setting it to `0L` in R is a one-line change.
- `-1` currently means "take the min over root states" (`ts_sankoff.cpp:84-87`), i.e. free the root's own label.
- **What breaks is that this is only half of the problem.** `forced_root_state` pins the root *state*; the rooting-dependence comes from the root *position*. Pinning the state while the position stays arbitrary makes the root meaningful exactly where the pipeline treats it as an artefact — so the criterion becomes well-defined only *relative to a root the user never chose*, and the oracle measures how much that unchosen choice is worth: up to 5 steps, on 84–117/120 topologies (Q-C). A pinned root state is therefore coherent **only** in combination with a pinned root position, which requires the user to supply a rooted tree, requires every whole-tree reroot site to be disabled (Q3), and requires the root position itself to become a searched parameter or an outgroup-fixed one. That is a feature, not a bug fix.

## Q3. What does pinning cost TBR?

**Much less than the finding's framing assumes — TBR is not the obstacle.** The task row says pinning "would additionally require rethinking TBR, since it reroots the clipped subtree". That conflates two different rerootings:

- **TBR's fragment reroot is a genuine unrooted-topology change**, not a rooting change. Reattaching a clipped fragment by a different internal edge yields a different *unrooted* tree. It is fully compatible with a pinned whole-tree root, and needs no change. (`apply_tbr_move` step 2, `ts_tbr.cpp:437-514`; `reroot_fragment`, `:544`.)
- **Whole-tree rerooting** is the part that moves the objective, and in the **default** path it does not happen: the `!phys_reroot` branch (`:2967-2984`) enumerates root-edge moves via `try_root_edge_moves` / `exact_verify_sweep` instead of physically rerooting. The physical sweep at `:2997` is legacy, behind `TS_PHYS_REROOT`.

The complete list of sites that move the whole-tree rooting is short and enumerable:

| Site | What it does |
|---|---|
| `src/ts_tbr.cpp:2997` | legacy physical-reroot sweep, `TS_PHYS_REROOT` only |
| `src/ts_fuse.cpp:20` (`reroot_at_tip0`), called at `:344`, `:381`, `:398` | reroots recipient and donors before fusing |
| `src/ts_rcpp.cpp:2088` | `ts_collapse_pool`'s tip-0 canonicalisation — the T-374 trigger site |
| `src/ts_sector.cpp:1101` | reroots at the sector HTU |

So the honest cost of pinning is "audit and gate four sites", not "rethink TBR". That is what makes Option 2 below *possible* rather than absurd — but see the caveats there; the reason to decline is Q-C and the blocker, not TBR.

Three code comments assert root-invariance and are false for HSJ/XFORM: `ts_tbr.cpp:123-124`, `ts_tbr.cpp:2999` (`// root-invariant; refreshes states`), and the already-recorded `ts_collapse_pool` comment in `ts_rcpp.cpp`. **Annotate-on-fix, not new findings.** (`ts_tbr.cpp:802-818` is *correctly* root-dependence-aware and needs no change.)

## Q4. If unrooted: is min-over-root-states enough?

**No.** `ts_sankoff.cpp:74-87` minimises over the root's own *state* at a fixed root *position*. That removes one degree of freedom out of two; the position remains, and the position is what varies.

What *would* restore invariance is **min over rootings as well as root states**: for an unrooted topology, minimise the total over all `2n − 3` placements of the root and all root states. This is a well-defined unrooted objective, and the oracle shows the current single-rooting score is always an **upper bound** on it (Q-E) — an overstatement, never an understatement. That matters: it means the existing kernel is a *sound but sometimes loose* evaluator of the correct unrooted objective, and it is tight for 87–98% of rootings.

For an asymmetric matrix that satisfies the triangle inequality this is the standard construction (equivalently: attach a hypothetical ancestor of unspecified state, and let the DP place it). Its cost is a factor `2n − 3` on the Sankoff term, which is why it belongs on the **report path only**, not in the search loop.

## Q5. Is there a cheaper acceptable outcome?

Yes, and it is the recommendation. The **actual P1 defect** in T-374 is not "the objective is rooting-dependent" — it is "`MaximizeParsimony` reports a best score that `TreeLength()` of its own returned trees does not reproduce", and "4 of 6 trees in one MPT set do not share a score under a common rooting". Both are *reporting* defects, and both are fixed by making one rooting authoritative at the boundary:

- Have the report path and `TreeLength()` agree — either both canonicalise (tip-0 rooting, as `ts_collapse_pool` already does) or, better, both evaluate the min-over-rootings objective for the returned pool only (`|pool| × (2n − 3)` Sankoff evaluations, negligible against a search).
- The min-over-rootings variant additionally makes the MPT set internally consistent by construction, since the score no longer depends on the representation.
- Document the bound: for a single block the search-time score can exceed the min-over-rootings objective by at most `nSec` (measured, below), and does so for 2–13% of rootings.
- Do **not** error out on unrooted input. XFORM has no root to demand, and requiring one would ask users to supply a hypothesis the criterion does not actually use.

---

## Evidence: `dev/red-team/heavy-tests/xform-rooting-oracle.R`

Pure R; independent re-implementation of the cost matrix (`recode_hierarchy.R:101-118`) and the DP (`ts_sankoff.cpp`), so agreement with the kernel is informative rather than circular. 120 random unrooted 9-tip topologies per scenario, scored under every one of the `2n − 3` edge rootings. Exits 0 when every stated prediction holds; it currently does.

| Scenario | dependent, `forced_root = -1` | max spread (bound `nSec`) | root-*position*-sensitive, `forced_root = 0` † | max spread † | % rootings at min | mean overstatement |
|---|---|---|---|---|---|---|
| `nSec=1` (2 levels) | 16/120 | 1 (≤ 1) | 108/120 | 3 | 96% | 0.042 |
| `nSec=2` (2×2) | 10/120 | 2 (≤ 2) | 117/120 | 4 | 98% | 0.022 |
| `nSec=3` (2×2×2) | 17/120 | 2 (≤ 3) | 115/120 | 5 | 96% | 0.042 |
| `nSec=2`, 20% ambiguous (`-1`) | 34/120 | 2 (≤ 2) | 97/120 | 4 | 87% | 0.167 |
| `nSec=2`, 20% present-unknown (`-2`) | 8/120 | 2 (≤ 2) | 113/120 | 4 | 98% | 0.018 |
| `nSec=0` control (symmetric) | **0/120** | 0 | 84/120 ‡ | 2 ‡ | 100% | 0 |

† **The two arms are not commensurable as "better/worse".** `forced_root = 0` defines a *different, explicitly rooted* criterion; its spread measures the cost of leaving that criterion's root **position** arbitrary, not a degradation of the free-root objective. The script verifies the pinned criterion is perfectly reproducible at a fixed canonical rooting, i.e. it is well-defined — the spread is entirely the unchosen-root cost.

‡ **The `nSec = 0` row is an asymmetry control for the FREE arm only** (0/120, spread 0 — see Q-A). Its pinned cell is *expected* to be non-zero and carries no information about asymmetry: constraining the root's label on a symmetric matrix necessarily makes the root position matter. Do not cite it as evidence that pinning is harmful.

**Q-A — confirmed, and the mechanism is localised.** The `nSec = 0` row is the control: with gain = loss = 1 the matrix is symmetric and the score is invariant across all 120 topologies. Every asymmetric row is dependent on a subset of topologies. So the dependence is caused by the gain/loss asymmetry, exactly as expected, and nothing else in the DP contributes.

**Q-B — the source of the dependence, with one corrected prediction.** The asymmetric part is a *gradient*: with `f(absent) = 0`, `f(present) = nSec/2`, we have `c(i,j) = s(i,j) + f(j) − f(i)` for symmetric `s = (c + cᵀ)/2`. Summing over edges oriented parent→child,

```
total  =  Σ_edges s(u,v)  −  Σ_{internal} f(u)  −  f(root)  +  Σ_tips f
```

`Σ_tips f` is constant and `s` carries no orientation, so the root node — charged `f` twice, once inside `Σ_internal` and once as `−f(root)` — is the entire source of the dependence. **An earlier, tighter prediction of `nSec/2` was falsified by this script** (`nSec = 1` gives spread 1; `nSec = 2` gives spread 2): it overlooked that rooting an unrooted tree *subdivides* an edge, so the edge set, and hence `Σ_edges s`, is not itself rooting-invariant. The retained and tested bound is **`nSec`** — `2·f(present)` — and it held in every scenario. Per-block, so a dataset with `k` blocks has a worst case of `Σ nSec_b`.

**Q-C — the cheap fix is incoherent on its own, and the numbers must be read carefully.** `forced_root_state = 0` ("absent at the root") does not restore rooting-invariance. Pinning the root *state* does nothing about the root *position*, and it removes the min-over-states freedom that was previously absorbing part of the positional difference: the pinned criterion's value then varies across root positions on 84–117/120 topologies, with spread up to 5.

**What that does and does not establish.** It does *not* establish that a rooted criterion is wrong, and the earlier draft of this document over-claimed exactly that. `forced_root_state = 0` defines a **different, explicitly rooted** criterion, and the script confirms it is perfectly reproducible at a fixed canonical rooting — it is well-defined. Its spread across rootings is therefore the price of leaving *its* root position arbitrary, and it is not commensurable with the free-root arm's spread as "worse". Nor is the `nSec = 0` pinned cell (84/120) evidence about asymmetry: constraining the root's label on a *symmetric* matrix must make the root position matter, so a non-zero value there is a tautology, not a contradiction. That row is a control for the free arm only.

**What it does establish, and it is enough to reject Option 1.** Shipping the one-liner into *today's* pipeline — which leaves `forced_root_state`'s root position unchosen, and actively moves it at the four sites in Q3 — would make the score depend on a root the user never selected, by up to 5 steps. The change is only meaningful **in combination with** pinning the position, i.e. it is Option 2 in disguise, not a one-line fix. Option 1 is rejected as incoherent, not as harmful.

**Consequently the recommendation rests on Q-E and the T-377 blocker, not on this section.** Q-E shows the current score is a tight, sound upper bound on a well-defined unrooted objective; the blocker shows the search is not even using the Sankoff term to choose moves. Those two carry the decision independently of how Q-C is read.

**Q-D — ambiguity is the aggravating factor, and does not break the bound.** Fully ambiguous tips (`-1`) roughly triple the dependence rate (34/120 vs 10/120 at `nSec = 2`) and raise the mean overstatement almost eightfold (0.167 vs 0.022), because a free tip lets the optimal labelling shift with the orientation. Present-but-unknown tips (`-2`) do not (8/120) — they still exclude "absent", which is where the asymmetry lives. The `nSec` bound survived both.

**Q-E — the current score is a loose upper bound, not a wrong number.** In every scenario an arbitrary rooting either attains the min-over-rootings objective or exceeds it; 87–98% of rootings attain it, and the mean overstatement is 0.02–0.17 steps. This is what makes Option 3 sound: the search is optimising a valid upper bound on a well-defined unrooted objective, and canonicalising the report path costs at most `nSec` per block.

**Reconciling with T-374's 165/300.** T-374 measured 55% of random topologies rooting-dependent under XFORM; this oracle sees 7–28% per single block. Not a contradiction: T-374 used real data with more tips and multiple blocks, and a dataset is dependent if *any* block is. Per-block rates of 10–30% compose to roughly half at three or four blocks.

---

## Options considered

**Option 1 — pin `forced_root_state = 0` only.** *Rejected as incoherent.* Q-C: it makes the root meaningful while leaving its position unchosen and, at four sites, moving — so the score comes to depend on a root the user never selected, by up to 5 steps on 84–117/120 topologies. Pinning the state is only meaningful together with pinning the position, which is Option 2. Do not ship it as a standalone one-liner.

**Option 2 — full rooted criterion: require a rooted input tree, pin state and position, gate the four whole-tree reroot sites.** *Rejected on cost/benefit, not feasibility.* Q3 shows the mechanical cost is moderate. But it (i) asks users for a root the criterion does not need, (ii) makes the root position an unsearched nuisance parameter whose choice changes the answer by up to `Σ nSec_b`, (iii) shrinks the reachable move set at the four gated sites — including fusing, whose whole value is topological diversity, and `ts_collapse_pool`, whose canonicalisation the collapse logic depends on, and (iv) delivers no user-visible benefit over Option 3, because Option 3 already removes the reported discrepancy. Reconsider only if a user presents a hierarchy where the root state is genuinely known *and* the tree genuinely rooted.

**Option 3 — recommended: unrooted reading, canonicalised reporting, documented bound.** Keep `forced_root_state = -1`. Keep the search as-is (an upper-bound evaluator). Make `MaximizeParsimony`'s reported score and `TreeLength()` agree on one rooting; prefer min-over-rootings on the returned pool, which additionally makes the MPT set self-consistent. Document XFORM as rooting-sensitive with the `nSec`-per-block bound. Fix the three false root-invariance comments as part of whatever touches those files.

**Option 4 — implement min-over-rootings inside the search.** Correct but `(2n−3)×` on the Sankoff term. Not warranted while the Sankoff term is absent from the candidate screen (blocker). Revisit only after that is fixed and only if measurement shows the loose bound misleads the search.

---

## The A/B measurement: recommend **against**, for now

The task asks whether a matched A/B against a rooting-pinned variant is worth running to settle whether TBR's accept/reject comparisons are incoherent across moves. **No — do not run it, and do not treat this as an open question.**

Three reasons:

1. **The only cheap pinned variant available is Option 1, and it is not the hypothesis you want to test.** An A/B against `forced_root_state = 0` would compare the status quo against a *different criterion* whose root position is itself unchosen and moving (Q-C) — so a score difference would not distinguish "pinning helps" from "the arm happened to draw a favourable root". Building a variant that actually tests the pinning hypothesis means pinning the position too, i.e. Option 2 — the kernel and reroot-site work this task explicitly forbids.
2. **The effect being measured is second-order behind a first-order defect.** The Sankoff term is not in the candidate screen at all (T-377). Any measured difference would be confounded with, and probably swamped by, Fitch-only candidate selection.
3. **The accept path is not actually incoherent in the way feared.** `actual = full_rescore` is authoritative at every accept (`ts_tbr.cpp:2717`, `:2761`, `:2787`), and the default path does not physically reroot mid-search (Q3), so accepted scores are consistent at a *stable* rooting within a TBR pass. Incoherence enters across the sites in the Q3 table — fusing, sector search, `ts_collapse_pool` — not within TBR's accept loop. That is a narrower and cheaper thing to reason about than a search-wide A/B.

If T-377 is fixed and a measurement is still wanted, the informative one is different: **the Sankoff-in-screen A/B** (Fitch-only screen vs Sankoff-aware screen), which is first-order, and which should be Hamilton-class. No local heavy compute either way.

---

## Follow-ups this document generates

- **Refinements folded into the existing T-377 row** (not a new finding — see the blocker section): the convergence sweep's blindness is `has_na`-conditional and `has_na` is data-dependent under XFORM; the `dominated` test is systematically over-permissive; and T-377 gates T-374's XFORM A/B. `ts_tbr.cpp:1503-1508`, `:1739`, `:2640-2642`, `:2717`. Coordinate with the concurrent `claude/t377-hierarchy-screen` worktree.
- **T-374 row:** point to this document for the XFORM half.
- **Annotate-on-fix:** false root-invariance comments at `ts_tbr.cpp:123-124` and `ts_tbr.cpp:2999`, alongside the already-recorded `ts_collapse_pool` one.
- **Docs:** `?MaximizeParsimony` and `?RecodeHierarchy` should state that XFORM scores are rooting-sensitive, bounded by `nSec` per block, and that the reported score is an upper bound on the min-over-rootings objective.
