# Cheap directional incremental scoring for NA — feasibility analysis

**Question (project lead, 2026-06-19):** design a cheap directional incremental
scoring approach for the inapplicable (NA) path — the EW/IW kind that makes
candidate evaluation O(1) after O(n) preprocessing — or prove why it can't be
done.

**Bottom line.** The *exact* O(1) "EW additive" directional scan is provably
**not available** for NA, for a concrete, code-confirmed reason (below). A richer
fixed-size-message DP that recovers O(1) is **not fundamentally impossible** (NA
parsimony is a linear-time tree DP), but it is **unsolved and research-grade**:
the one directional NA message that exists today (`fitch_na_indirect_length`) is
deliberately *approximate*, precisely because the exact context is non-local. The
**practical** cost lever is therefore incremental *exact* rescore — O(affected),
not O(1) — building on machinery that already exists (`fitch_na_dirty_*`). This is
the recommended direction for task #18.

---

## 1. Why EW/IW directional is O(1) per candidate

For equal/implied weights the per-edge cost is **2-local**: clipping edge (u,w)
and rejoining, the length is

    len(A) + len(B) + join(prelim_A, prelim_B)

where `len(A)`, `len(B)` are constants and `join` depends only on the two
preliminary state-sets at the cut. Two facts make the all-candidates scan cheap:

1. **2-locality** — the cost the candidate adds is a function of just two
   fixed-size sets at the broken edge.
2. **Exact, one-pass from-above** — `compute_from_above` yields, at every node,
   the prelim set of "everything outside" that node, so the presented set at
   *any* rerooting is `join(from_above[sc], prelim[sc])` in O(1).

`try_root_edge_moves` (EW) uses exactly this: `base_split = best_score - rootjoin`
is constant, and each candidate costs `base_split + fitch_join(stateL, stateR)`.

## 2. Why the EW additive trick is provably unavailable for NA

Two independent obstructions, both visible in the code:

**(a) The NA per-node step is NOT 2-local.** Pass 3 of the Brazeau three-pass
counts a step at a node from (`ts_fitch_na_incr.h:359`):

    needs_step = l_act & r_act & ~(ss_app & any_isect)

where `ss_app` comes from `final_` (the Pass-2 *uppass*, i.e. whole-tree context)
and `l_act/r_act` from `subtree_actives`. The count depends on the global
applicability resolution (which tree regions are "applicable"), not just two
prelim sets at the cut. The code states this directly: the additive split
"does NOT hold for IW (concave) or NA (3-pass)" (`ts_tbr.cpp`, `try_root_edge_moves`
header), and routes NA to apply+`full_rescore` (`try_root_edge_moves_rescore`).
This is *why* the root-edge completeness fix (2026-06-19) had to use rescore.

**(b) Structural reason — the region term is not a 2-local edge cost.** The NA
score includes the number of *applicable regions* (maximal connected applicable
components). This is a global connectivity functional. It is not a sum of
symmetric per-edge costs: a single applicable region with 3 inapplicable
neighbours (3 boundary edges, 1 region) and three separate applicable leaves
(3 boundary edges, 3 regions) have identical boundary-edge multisets but
different region counts. So no fixed *symmetric* cost matrix on {states, ⊥}
reproduces the term, and the EW "constant base + 2-local join" form cannot
represent a quantity that changes with global region structure as the halves are
rerooted.

## 3. Could a richer fixed-size-message DP recover exact O(1)?

Not impossible in principle. NA parsimony is a **linear-time tree DP** (down1 /
up1 / down2), and any linear tree DP is, by construction, a fixed-size-message
scheme — so an "all-rerootings" directional form is conceivable, with a message
richer than a single set:

- message = cost table indexed by (boundary state, region-status bit), i.e.
  "min internal cost of this fragment given its boundary resolves to state s and
  its region is open/closed across the cut";
- the region term *can* be charged locally in a **directed/rooted** formulation
  (charge +1 at each applicable node whose parent is inapplicable = one charge
  per region's top), which is fixed-size;
- combine at a join would re-optimise the boundary in O(table²).

**But three things make this research-grade, not a quick win:**

1. **The existing directional NA message is already approximate.**
   `fitch_na_indirect_length` *is* the attempt at an O(blocks)/candidate NA
   message. It approximates the candidate's context by reusing the **base tree's**
   `final_` sets at the attachment nodes (`tree.final_[a_base]`,
   `tree.final_[d_base]`) — exact for Fitch's from-above, but only a proxy for the
   NA uppass, which is not reconstructable as a simple set. That a *sound* version
   was not used (and the expensive exact sweep was built instead) is direct
   evidence the exact O(1) message is non-trivial.
2. **Exact from-above is the crux.** EW gets exact context from one
   `compute_from_above` pass. NA would need an exact *oriented* from-above message
   for the down1/up1 algebra AND a region-aware combine; deriving and proving
   these correct (especially that (down, up) messages are a sufficient statistic
   for re-optimising across an arbitrary new join, including region merge/split)
   is the open problem.
3. **Validation burden.** Any candidate design must match `exact_verify_sweep`
   (now complete, post the 2026-06-19 root-edge fix) bit-for-bit on a corpus, plus
   the oracle 0/N. The history here (T-300 "unresolved −3", the approximate
   indirect) shows how easy it is to get subtly wrong.

**Status: not disproven, but unsolved and high-risk. Do not attempt before the
cheaper lever below is exhausted.**

## 4. Recommended lever for #18: incremental EXACT rescore (O(affected), not O(1))

The sweep currently does `apply_tbr_move + full_rescore` (full Passes 1+2+3,
O(n)) for *every* candidate. Two sub-levers, both reusing existing machinery:

- **Pruning via the approximate scan.** Run `fitch_na_indirect_length` (cheap)
  first and `full_rescore` only candidates it cannot rule out. **Precondition:**
  confirm the approximation's *direction* — it is only a safe filter if it never
  *over*-estimates the true improvement (i.e. it must not hide a real improver).
  If it can under-estimate the resulting length, it can be used as an admissible
  lower bound; if it can over-estimate, it cannot prune safely. Verify before use.
- **Localised Pass-3 delta.** `fitch_na_dirty_*` already does incremental Passes
  1+2 over the union of affected rootward paths (built for the SPR accept path).
  The remaining O(n) is the Pass-3 step recount (full postorder). Extend the dirty
  machinery to accumulate only the Pass-3 *delta* over affected nodes. **Caveat:**
  the Pass-2 uppass propagates context from the changed path *down into off-path
  subtrees*, so the affected set for an exact NA recount can exceed the path; bound
  it carefully (this is exactly where the region term bites).

Either reduces per-candidate cost from O(n) toward O(depth)/O(affected) — a large
win on the 74/88-tip datasets — without the directional-message risk.

## 5. Validation harness (already in place)

- `dev/benchmarks/tbr_oracle_na.R` (real data) and
  `dev/benchmarks/tbr_oracle_na_small.R` (fast, high-N) — completeness 0/N.
- `tests/testthat/test-ts-na-complete.R` — pins the Zanol2014 start-#14 optimum.
- For any cost change: assert per-tree score equality vs the current
  `full_rescore` path on a corpus (the cost change must be score-transparent),
  then re-run the oracles.

## 6. Recommendation

1. Do **not** build the exact O(1) directional NA scan now — provably can't reuse
   the EW additive form, and the sound richer-message version is research-grade.
2. Pursue **incremental exact rescore** for #18: pruning first (cheapest, verify
   approximation direction), then localised Pass-3 delta on top of
   `fitch_na_dirty_*`.
3. Keep `exact_verify_sweep` (now complete) as the exact ground-truth oracle that
   any faster path is validated against.
