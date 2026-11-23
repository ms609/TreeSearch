# Full Polytomy Search for TreeSearch C++ Engine

**Status:** IN PROGRESS  
**Target branch:** `feature/polytomy-search` (from `cpp-search`)  
**Target worktree:** `../TS-Polytomy`

## Motivation

The TNT benchmark (2026-03-20, `TS-TNT-bench` worktree) shows TreeSearch
falls 1–14 steps behind TNT on datasets with ≥50 taxa. The TNT
outperformance analysis identifies **tree collapsing during search** as the
single biggest remaining algorithmic gap:

> "Searches that collapse branches with minimum possible length produce more
> effective searches than criteria which collapse fewer branches, both in
> terms of time needed to complete searches, and ability to find shortest
> trees." — Goloboff (2023), Cladistics 39: 229–238

The existing `ts_collapsed.h/.cpp` (clip-skipping) was a partial step
toward this, but benchmarks showed 0% skip rate on standard morphological
data because near-optimal **binary** trees have few zero-length edges. The
key insight is that collapsing those edges into polytomies *changes the
search topology space*, making TBR/SPR more efficient by eliminating
distinctions that carry no phylogenetic signal.

### Key literature

| Reference | Key contribution |
|-----------|-----------------|
| Goloboff (1996), "Methods for faster parsimony analysis", Cladistics 12: 199–220 | §"Collapsing The Trees": partial reoptimization, shortest-path shortcut, asymmetric reachability |
| Goloboff & Farris (2001), "Methods for quick consensus estimation", Cladistics 17: S26–S34 | TBR-collapsing rule: collapse all nodes between source/dest for equal-length rearrangements |
| Goloboff (2023), "Searches, implied weights, and tree collapsing", Cladistics 39: 229–238 | Empirical comparison of collapsing criteria during search; "minimum possible length 0" recommended |
| Day et al. (1985) / TreeDist | O(n·k) strict consensus via compatible-splits method; available in TreeDist |

### Detailed literature notes (from PDF review 2026-03-22)

**Goloboff 1996 — §"Collapsing The Trees" (pp. 213–218)**

1. *Shortest-path test (approximate)*: If no node in the path between the
   clipped subtree's original position and the destination is "supported"
   (has character-state change), the rearranged tree collapses to the same
   polytomy as the original. The tree can be discarded without full
   reoptimization. This is the core shortcut that our collapsed-region
   skipping approximates.

2. *Asymmetric reachability*: The shortest-path shortcut creates directed
   connectivity — swapping on tree A may find B, but swapping on B may not
   find A. Goloboff gives an explicit example (5 taxa, `x000 a100 b011
   c111 d111`) where the dichotomous tree A can reach the trichotomous
   tree B, but no resolution of B can reach A because the movement would
   cross only unsupported nodes. He argues this is acceptable: "heuristic
   searches cannot guarantee finding all of the optimal trees, or even any
   of them—with or without shortcuts."

3. *Efficient collapsing via final states*: For characters where the final
   state sets don't change after rearrangement (checked by comparing basal
   node of clipped subtree against ancestor/descendant of destination
   branch), only 10–20% of characters need reoptimization for collapsing.
   This maps to our incremental scoring infrastructure.

4. *Union construct method*: A further optimization that evaluates
   destinations "en masse" by computing union state sets for subtrees and
   rejecting entire branches when the union construct produces suboptimal
   length. Achieved 50% time reduction on congruent datasets (168 taxa),
   but no gain on incongruent data.

**Goloboff & Farris 2001 — "Methods for quick consensus estimation"**

1. *TBR-collapsing rule*: "when a rearrangement produces a tree of the
   same length as the one being swapped, collapsing all of the nodes
   between source and destination (and new root, in the case of TBR)."
   This is equivalent to saving all equal-length trees and computing their
   strict consensus, but uses no extra memory and less time.

2. *SPR vs TBR collapsing*: TBR-based collapsing eliminates more
   spurious groups than SPR-based, with minimal loss of correct groups.
   On Zilla (500 taxa): SPR collapsing gives 79.6% true nodes recovered
   with 0.63% error rate; TBR gives 79.0% true nodes with 0.48% error
   rate. Net effect: TBR collapsing is more reliable.

3. *RFD (Relative Fit Difference)*: Extends collapsing to suboptimal
   trees by measuring `(F-C)/F` where F = favorable fit, C = contradictory
   fit. Nodes with RFD below a threshold Q are collapsed. When calculating
   rearrangement length, as soon as length increase X > D/(1-Q), the
   rearrangement can be abandoned. For Q=0.10, tree collapsing takes only
   5% additional time. This could be a future extension (post-2.0.0).

4. *Pool benefit*: Collapsing trees during swapping means different
   dichotomous trees that differ only in "minor" rearrangements collapse
   to the same polytomy. The pool then stores more topologically diverse
   trees, improving search effectiveness. This directly validates our
   Phase 5 (collapsed-topology pool dedup).

**Goloboff & Morales 2023 — TNT version 1.6**

1. *Consensus stabilization*: TNT's driven search can stop when the
   strict consensus is stable after N hits — analogous to TreeSearch's
   `consensusStableReps`. TNT's parallel mode has a coordinator that
   centralizes consensus calculation.

2. *Parallel architecture*: "Builders" create trees via Wagner+TBR+
   sectorial/ratchet/drift, pass them to a "fuser" task. Similar to
   TreeSearch's `ThreadSafePool` pattern but using PVM processes rather
   than threads.

3. *Fast consensus*: The user notes that Day et al. (1985) O(n·k) strict
   consensus is available via the TreeDist package. This could replace or
   supplement the XOR-hash consensus approximation in `ts_pool.cpp` for
   more accurate stability detection. Not needed for the polytomy search
   itself, but relevant for improving consensus-stability stopping.

### What TNT does

TNT collapses zero-length branches **during search** by default (`collapse
3;` = TBR-rule). After each TBR rearrangement is accepted, zero-length
edges are contracted into polytomies. TBR then operates on the collapsed
(non-binary) tree, which has fewer edges to clip and regraft through. The
key benefits are:

1. **Fewer TBR candidates**: a polytomous tree with k collapsed edges has
   ~2k fewer clip candidates and ~2k fewer regraft positions per clip.
2. **Pool deduplication**: collapsed trees that differ only in unsupported
   resolution are identical, preventing the pool from filling with
   trivially different trees.
3. **Better convergence**: the search explores "real" topological
   differences rather than wasting effort on unsupported resolutions.

---

## Design decision: Approach B (collapsed-edge set, binary internals)

After reviewing the codebase, **Approach A** (replacing `left[]`/`right[]`
with multi-child representation) would require rewriting every module — TBR,
SPR, NNI, Fitch scoring, NA scoring, incremental scoring, undo stacks,
Wagner construction, constraint checking, sectorial search, fusing, splits.
This is estimated at 10+ weeks and carries extreme regression risk.

**Approach B** is both faster to implement and closer to what TNT actually
does. TNT stores trees as binary internally but maintains a set of
"collapsed" edges that modify candidate enumeration and pool comparison.
The binary topology is always available for scoring; collapsed edges just
indicate which resolutions are unsupported.

### Core idea

Maintain a `std::vector<uint8_t> collapsed` flag array alongside the
existing binary `TreeState`. After each accepted TBR/SPR move + full
rescore:

1. **Recompute collapsed flags** (already implemented in `ts_collapsed.cpp`)
2. **Skip collapsed clips** in TBR/SPR/drift candidate enumeration
3. **Skip collapsed regraft distinctions**: when regrafting into a region
   of consecutive collapsed edges, all positions within that region
   produce the same score — evaluate only one representative position
4. **Pool comparison uses collapsed form**: two binary trees that collapse
   to the same polytomy are treated as duplicates

### Why this works without changing TreeState

- Scoring uses the binary tree (exact Fitch downpass/uppass, unchanged)
- Topology manipulation uses binary operations (SPR clip/regraft, unchanged)
- Only candidate **enumeration** changes (skip/merge collapsed regions)
- Pool comparison adds a collapsed-topology hash alongside the existing
  binary split hash

The binary tree is always there as a "refinement" of the collapsed tree.
When a move is accepted that resolves a polytomy (puts signal on a
previously zero-length edge), the collapsed flag simply clears.

---

## Implementation plan

### Phase 1: Collapsed-region identification (extend existing code)

**Files:** `src/ts_collapsed.h`, `src/ts_collapsed.cpp`

The existing `compute_collapsed_flags()` already identifies edges where
clipping cannot improve score. Extend this to also identify **collapsed
regions** — maximal connected subsets of collapsed edges forming a
polytomy:

```cpp
struct CollapsedRegion {
  int representative;    // one node in the region (for regraft targeting)
  int n_edges;           // number of collapsed edges in this region
  std::vector<int> nodes; // all nodes with collapsed[node] == 1 in region
};

struct CollapsedInfo {
  std::vector<uint8_t> collapsed;       // per-node flag (existing)
  std::vector<int> region_id;           // per-node: which region (-1 if not collapsed)
  std::vector<CollapsedRegion> regions; // the collapsed regions
  int n_collapsed = 0;                  // total collapsed edges
};

void compute_collapsed_info(
    const TreeState& tree,
    const DataSet& ds,
    CollapsedInfo& info);
```

This is a simple post-processing step after the existing flag computation:
BFS/DFS from each collapsed node, grouping connected collapsed edges.

**Estimated effort:** 1–2 days

### Phase 2: TBR clip skipping (already partially done)

**Files:** `src/ts_tbr.cpp`

The current code already skips collapsed clips when `!collect_pool`. Verify
this is working correctly and add a **diagnostic counter** (`n_collapsed_skipped`)
to the TBR return value for benchmarking.

No code change needed beyond the diagnostic counter — Phase 1's extended
flags subsume the existing implementation.

**Estimated effort:** 0.5 days

### Phase 3: TBR regraft region merging (the main win)

**Files:** `src/ts_tbr.cpp`

This is the key new optimization. When evaluating regraft positions for a
non-collapsed clip:

**Current behavior:** enumerate all main-tree edges as regraft candidates,
evaluate each independently.

**New behavior:** for each collapsed region, evaluate only **one
representative regraft position** within the region. All positions within a
collapsed region produce identical scores (because the intermediate nodes
have zero cost and identical state sets — exactly the conditions verified
by `compute_collapsed_flags()`).

Implementation in the TBR regraft loop:
```cpp
for (auto& [above, below] : main_edges) {
  // Skip redundant positions within collapsed regions
  if (collapsed_info.collapsed[below] &&
      collapsed_info.region_id[below] == last_evaluated_region) {
    continue;  // same region, same score — skip
  }
  last_evaluated_region = collapsed_info.region_id[below];

  // ... evaluate regraft as before ...
}
```

**Correctness argument:** Within a collapsed region, all edges have:
- Zero local cost at parent (condition 1–2 of collapsed flags)
- `prelim[sibling] == prelim[parent]` (condition 3)
- `down2[sibling] == down2[parent]` (condition 4, NA)
- `subtree_actives[sibling] == subtree_actives[parent]` (condition 5, NA)

Therefore the `final_` states used by `fitch_indirect_length()` at any
edge within the region produce the same `vroot` value, giving identical
scores for all regraft positions in the region.

**Important subtlety:** The best regraft position's `(above, below)` pair
matters for the actual topology after the move. When a collapsed-region
regraft is chosen, we regraft at the representative position. The resulting
tree will have a different binary resolution of the polytomy, but the same
score and the same collapsed topology. This is equivalent to TNT's behavior.

**Estimated effort:** 3–5 days (careful correctness verification needed)

### Phase 4: SPR and drift integration

**Files:** `src/ts_search.cpp`, `src/ts_drift.cpp`

Apply the same clip-skipping (already in Phase 2) and regraft-merging
(Phase 3 pattern) to SPR search and drift search.

For drift: suboptimal-acceptance moves should still skip collapsed clips
(a collapsed clip cannot improve OR change the score, so accepting it
is always a no-op). Regraft merging applies identically.

**Estimated effort:** 2–3 days

### Phase 5: Pool deduplication using collapsed form

**Files:** `src/ts_pool.h`, `src/ts_pool.cpp`, `src/ts_splits.h`

Currently pool deduplication uses binary split hashes. Two trees that
differ only in unsupported resolution have different split hashes but
should be considered duplicates.

**Add collapsed-topology hashing:**
1. After computing collapsed flags, identify the "collapsed splits" —
   the splits that remain after contracting all collapsed edges.
2. Hash only the non-collapsed splits for pool dedup.
3. Use this hash as the primary dedup key; fall back to binary hash
   for trees with no collapsed edges (fully resolved).

Implementation:
```cpp
uint64_t compute_collapsed_hash(
    const TreeState& tree,
    const CollapsedInfo& info,
    int n_tip);
```

This is a filtered version of the existing `compute_splits()` +
`hash_single_split()` pipeline — just skip splits corresponding to
collapsed edges.

**Estimated effort:** 2–3 days

### Phase 6: Ratchet interaction

**Files:** `src/ts_ratchet.cpp`

During ratchet perturbation, character weights change, which means
collapsed flags must be recomputed after perturbation. The ratchet already
calls `tbr_search()` which recomputes flags after each accepted move, so
this should work automatically.

**One subtlety:** After ratchet perturbation (upweighting/zeroing chars),
some previously collapsed edges may become non-collapsed (the perturbed
weights create artificial signal). This is correct behavior — the
perturbation should explore the full binary space.

After ratchet un-perturbation (restoring original weights), the full
rescore will re-establish correct collapsed flags.

**Estimated effort:** 1 day (verification + edge case testing)

### Phase 7: Sectorial search interaction

**Files:** `src/ts_sector.cpp`

For sectorial search, collapsed flags should be computed on the full tree
and passed to the sector TBR. Within a sector:
- Clip candidates that are collapsed in the full tree remain collapsed
- Regraft merging applies within the sector

Collapsed flags for the **reduced dataset** (sector subproblem) should be
recomputed from the sector's own scoring, not inherited from the full tree.

**Estimated effort:** 2–3 days

### Phase 8: Wagner tree collapsing

**Files:** `src/ts_wagner.cpp`

After Wagner tree construction, compute collapsed flags before the first
TBR pass. Wagner trees typically have many zero-length edges (the greedy
construction often creates unsupported resolutions), so this is where
collapsed-region merging may have the biggest per-tree impact.

**Estimated effort:** 0.5 days

### Phase 9: Testing

**Files:** `tests/testthat/test-ts-polytomy-search.R` (Tier 2)

1. **Region identification:** hand-built trees with known collapsed
   regions; verify region count and membership.
2. **Regraft merging correctness:** verify that evaluating all positions
   vs. one-per-region gives identical best scores.
3. **Pool collapsed-hash dedup:** two trees differing only in zero-length
   resolution are treated as duplicates.
4. **Score equivalence:** driven search with collapsed optimization
   produces same or better scores than without.
5. **IW/Profile mode compatibility.**
6. **NA dataset compatibility.**
7. **Ratchet interaction:** collapsed flags correctly update after
   perturbation and un-perturbation.
8. **End-to-end regression:** run existing benchmark datasets, verify
   no score degradation.

**Estimated effort:** 3–4 days

### Phase 10: Benchmarking

Re-run the TNT benchmark comparison with collapsed search enabled:
- Same 14 datasets, EW Fitch, 10s and 30s timeout
- Compare scores, timing, and replicates completed
- Focus on the 5 datasets where TreeSearch fell behind

Also measure:
- Collapsed edge percentage per dataset (at optimum)
- Regraft candidates skipped per TBR pass
- Pool duplicate reduction

**Estimated effort:** 1–2 days

---

## Risk assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| Regraft merging incorrectly skips a productive position | HIGH | Formal correctness proof + extensive unit tests; conservative fallback to evaluate all if collapsed count is low |
| Collapsed flags stale after ratchet perturbation | MEDIUM | Flags always recomputed after full_rescore; verify in ratchet tests |
| Pool collapsed-hash collisions (different topologies hash same) | LOW | Conservative direction (over-dedup); hash collision = treat as duplicate = miss one tree, not wrong scores |
| Negligible benefit on dense morphological data | MEDIUM | TNT benchmarks show the benefit is real; if our data shows otherwise, document and stop |
| Interaction with MPT enumeration | HIGH | Collapsed optimizations MUST be disabled during `collect_pool` (equal-score exploration); already guarded in existing code |

---

## Estimated total effort

| Phase | Days | Cumulative |
|-------|------|------------|
| 1. Collapsed regions | 1–2 | 1–2 |
| 2. TBR clip (existing) | 0.5 | 1.5–2.5 |
| 3. TBR regraft merging | 3–5 | 4.5–7.5 |
| 4. SPR + drift | 2–3 | 6.5–10.5 |
| 5. Pool dedup | 2–3 | 8.5–13.5 |
| 6. Ratchet | 1 | 9.5–14.5 |
| 7. Sectorial | 2–3 | 11.5–17.5 |
| 8. Wagner | 0.5 | 12–18 |
| 9. Testing | 3–4 | 15–22 |
| 10. Benchmarking | 1–2 | 16–24 |

**Total: 16–24 agent-days.** Substantially less than the 9–13 weeks
estimated for Approach A (full polytomy representation).

---

## Literature review — COMPLETE (2026-03-22)

All three papers reviewed from PDF. Key algorithmic details extracted
in the "Detailed literature notes" section above. The Goloboff (2023)
paper on collapsing criteria was not available in PDF but its core
recommendation ("minimum possible length 0" during search) is documented
in the AGENTS.md architecture reference.

---

## Success criteria

1. **Score parity or improvement** on all 14 TNT benchmark datasets
   (no regressions)
2. **Measurable collapsed-edge skip rate** (>0%) on at least the harder
   datasets (Wortley2006, Eklund2004, Zanol2014, Zhu2013, Giles2015)
3. **All existing tests pass** (1859 ts-* tests + full R-level suite)
4. **New test file** with ≥15 assertions covering all phases

---

## References

- Goloboff, P. A. (1996). Methods for faster parsimony analysis. Cladistics, 12, 199–220.
- Goloboff, P. A. & Farris, J. S. (2001). Methods for quick consensus estimation. Cladistics, 17, S26–S34.
- Goloboff, P. A. (2023). Searches, implied weights, and tree collapsing. Cladistics, 39, 229–238.
- Goloboff, P. A. & Catalano, S. A. (2016). TNT version 1.5. Cladistics, 32, 221–238.
