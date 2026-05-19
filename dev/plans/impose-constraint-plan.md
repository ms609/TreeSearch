# Plan: C++ `impose_constraint()` for post-hoc topology repair

## Motivation

Several operations in the search pipeline can produce constraint-violating
trees but currently lack a way to *repair* violations cheaply:

| Operation | Current handling | Cost |
|-----------|----------------|------|
| **NNI perturbation** | Disabled entirely when constraints active (T-209) | Loses primary topology-space escape in `thorough` preset |
| **Fuse** | Posthoc check then revert (discard fused tree) | Wastes fuse work; reduces pool diversity under constraints |
| **Drift** | Move rejection (same as TBR) | Narrower exploration under constraints |

A C++ `impose_constraint(TreeState&, const ConstraintData&)` function would
take an existing tree with minor violations and minimally rearrange it to
satisfy all constraint splits. This is the "repair" complement to the
existing "prevention" approach (`regraft_violates_constraint`).

**Not in scope:** TBR/SPR candidate screening. Move rejection remains the
right approach there (O(1) per candidate vs. O(n) for fixup + rescore).

## Existing infrastructure to reuse

- `map_constraint_nodes(tree, cd)` — identifies which splits are satisfied
  (`constraint_node[s] >= 0`) and which are violated (`== -1`). Internally
  computes per-node subtree tip bitmasks via postorder traversal, but this
  `node_tips` array is a **local variable** (not stored on ConstraintData).
  **Refactoring needed:** extract the tip-bitmask computation into a shared
  helper `compute_node_tips(tree, n_words)` that both `map_constraint_nodes`
  and `impose_constraint` can call.
- `spr_clip(node)` / `spr_regraft(above, below)` — existing SPR primitives
  that detach a subtree and reattach it at a new edge. Work on tips and
  internal nodes alike. When clip_node is a tip, clip detaches the tip and
  frees its parent node; regraft reuses that parent as the new internal node.
  **No new topology-manipulation primitives are needed.**
- `split_tips[s * n_words .. (s+1)*n_words - 1]` — target tip set per split
  (canonicalized: tip 0 always "outside")
- `update_constraint(tree, cd)` — combined remap + DFS timestamp refresh
- TreeTools `ImposeConstraint()` — R reference implementation (polytomy
  backbone + resolution). Our approach is different: minimal surgical repair
  rather than full rebuild, preserving perturbation diversity.

### Thread safety

Each worker thread in `ts_parallel.cpp` makes a local copy of
`ConstraintData` (lines 95–99). `impose_constraint` mutates `cd`
(via `map_constraint_nodes` and `update_constraint`), which is safe
since each thread operates on its own copy.

### State array lifecycle

`impose_constraint` only modifies topology (parent/left/right) and
constraint metadata. It does **not** touch Fitch state arrays (prelim,
final_, local_cost, etc.). The caller is responsible for calling
`tree.reset_states(ds)` + `score_tree(tree, ds)` after repair, which
rebuilds all state arrays from scratch. This means intermediate state
array inconsistency during the SPR moves is harmless.

## Algorithm

### High-level

For each violated split, find the internal node whose subtree is closest to
the target tip set. Identify **subtrees** of misplaced taxa and SPR-move them
to the correct side of the tree. Process splits from smallest to largest (this
is provably safe for compatible constraints; see correctness note below).

### Detailed steps

```
impose_constraint(TreeState& tree, ConstraintData& cd):

  1. map_constraint_nodes(tree, cd)
     -> constraint_node[s] == -1 for violated splits
     -> node_tips[] bitmask array computed as side effect

  2. If no violations, return (common case — free)

  3. Collect violated splits; sort by popcount ascending (smallest first)

  4. For each violated split S (in ascending size order):

     a. Rebuild node_tips bitmasks via postorder traversal
        (reuse the same buffer; needed because previous split's
        moves changed the topology)

     b. Find the "best candidate node" N — the internal node that
        MINIMIZES |symmetric_difference(subtree(N), target(S))|:
          cost(N) = popcount(node_tips[N] XOR split_tips[S])
        Iterate postorder; keep track of minimum.

     c. Compute misplaced tip sets via bitmask:
          move_out_mask = node_tips[N] AND NOT split_tips[S]
            (tips in N's subtree that shouldn't be)
          move_in_mask  = split_tips[S] AND NOT node_tips[N]
            (tips outside N's subtree that should be inside)

     d. Find maximal misplaced subtrees (not individual tips):
        For each direction (move_out, move_in):
          In postorder within the relevant tree region, find nodes
          whose subtrees are entirely contained in the misplaced set:
            (node_tips[v] & ~move_xxx_mask) == 0
          Keep only maximal ones (parent's subtree is NOT entirely
          contained). These are the subtrees to clip.

     e. For each misplaced subtree root M:
          - Skip if M is a direct child of the tree root (see edge
            cases below)
          - Pick a random target edge:
            * move_out: DFS from root, collect edges NOT in N's
              subtree; pick one uniformly at random
            * move_in:  DFS from N, collect edges within N's
              subtree; pick one uniformly at random
          - tree.spr_clip(M)
          - tree.spr_regraft(target_above, target_below)
          (Each clip-regraft pair is a self-contained SPR.
          The single clip_state slot is overwritten each time,
          which is fine since we never undo these moves.)

     f. tree.build_postorder()
        (makes tree valid for next split's bitmask computation)

  5. update_constraint(tree, cd)
     (remaps constraint nodes + refreshes DFS timestamps for
     subsequent TBR/SPR/temper/drift move screening)

  6. Full rescore (caller's responsibility)
```

### Why minimum symmetric difference, not maximum overlap

The candidate selection criterion is `popcount(node_tips XOR split_tips)`,
which counts the total number of misplaced tips (both directions).
"Maximum overlap" (`popcount(node_tips AND split_tips)`) can prefer nodes
with many extra tips that require more moves:

| Node | Subtree | Target | Overlap | Sym. diff |
|------|---------|--------|---------|-----------|
| X | {A,B,C} | {A,B,C,D} | 3 | 1 move |
| Y | {A,B,C,D,E,F} | {A,B,C,D} | 4 | 2 moves |

Max-overlap picks Y; min-symmetric-diff picks X (correct).

### Why subtrees, not individual tips

`random_nni_perturb()` swaps subtrees at each edge. After perturbation,
misplaced items are typically contiguous subtrees in the current tree.
With `fraction = 0.5` on a 100-tip tree and a 20-tip constraint clade,
5–10 misplaced tips might share 2–3 subtree roots. Moving subtrees
instead of tips halves the number of SPR operations, and each SPR has
the same cost regardless of subtree size.

Finding maximal subtrees is cheap with the bitmask infrastructure:
```cpp
// For move_out direction within N's subtree:
for (int node : postorder_within_N) {
  uint64_t* nt = &node_tips[node * n_words];
  bool all_misplaced = true;
  for (int w = 0; w < n_words; ++w) {
    if (nt[w] & ~move_out_mask[w]) { all_misplaced = false; break; }
  }
  if (all_misplaced) {
    // Check parent isn't also all-misplaced (maximality)
    // If maximal: add to clip list
  }
}
```

### Correctness of smallest-first ordering

For compatible constraints (required by tree construction), splits are either
nested or disjoint:
- **Nested (S1 ⊂ S2):** Fixing S1 first moves tips within S2's boundary.
  When S2 is processed, S1's tips are already correctly placed, so S2's
  repair only touches non-S1 tips — no interaction.
- **Disjoint:** Fixing S1 moves tips that are not in S2's target set
  and vice versa — no interaction.

Therefore smallest-first is **provably correct for compatible constraints**;
fixing one split cannot violate another. No re-checking of previously
fixed splits is needed.

## Edge cases

### Root-adjacent clips

`spr_clip()` has an awkward path when `clip_node` is a direct child of the
root (`parent[clip_node] == root == n_tip`). The existing code (lines
304–320 of ts_tree.cpp) handles this but the comments acknowledge it's
unusual. If `impose_constraint` ever needs to move a root child, the tree
is likely so heavily violated that repair is the wrong strategy.

**Guard:** Skip subtrees whose parent is the root. If any remain after
this filter, fall back to `random_constrained_tree()` (full rebuild).

### Best candidate IS the root

The root's subtree tip mask has all bits set. Since constraint splits are
canonicalized with tip 0 outside, the root will always have a large
symmetric difference with any split. So the root is never the best
candidate. No special handling needed.

### Bail-out for heavy violations

If the total number of subtree moves exceeds a threshold (e.g., n_tip / 4),
the tree is so disrupted that surgical repair offers little advantage over
building fresh. In this case, fall back to `random_constrained_tree()`.
This also guards against pathological NNI perturbation scenarios.

## Integration points

### 1. NNI perturbation (highest value)

In `nni_perturb_search()` (ts_nni_perturb.cpp), between
`random_nni_perturb()` (line 82) and `tree.reset_states(ds)` (line 92):

```cpp
int n_swaps = random_nni_perturb(tree, params.perturb_fraction);
// NEW: repair constraint violations from blind NNI perturbation
if (n_swaps > 0 && cd && cd->active) {
  impose_constraint(tree, *cd);
}
tree.reset_states(ds);
score_tree(tree, ds);
// TBR to local optimum (existing code, constraint-aware)
```

Remove the `(!cd || !cd->active)` gate in `run_single_replicate()`
(ts_driven.cpp:322).

### 2. Fuse (medium value)

In `driven_search()` (ts_driven.cpp:719–727), replace the discard with
repair:

```cpp
bool fuse_ok = true;
if (cd && cd->active) {
  fuse_ok = !violates_constraint_posthoc(fused, *cd);
  if (!fuse_ok) {
    impose_constraint(fused, *cd);
    fused_score = score_tree(fused, ds);
    fuse_ok = true;  // repaired
  }
}
if (fuse_ok) {
  std::vector<uint8_t> fused_collapsed;
  compute_collapsed_flags(fused, ds, fused_collapsed);
  pool.add_collapsed(fused, fused_score, fused_collapsed);
}
```

### 3. Sector search (low value — optional)

In `xss_search()` and `rss_search()` (ts_sector.cpp:746, 876), the
posthoc violation check currently reverts the sector to its previous
topology. An alternative: repair + keep. However, sector violations arise
from a local rebuild that changed the sector–rest-of-tree relationship,
so repair might undo the sector's improvements. **Defer to benchmarking
before committing to this integration.**

### 4. Random tree (not needed)

`random_constrained_tree()` already handles this via the polytomy
approach, which is better for building from scratch (constructive,
properly random). `impose_constraint` is not needed here.

## Testing strategy

1. **Unit test:** Build a tree that violates a known constraint.
   Call `impose_constraint`. Verify all splits satisfied.
   Verify score matches `score_tree()` of the result.

2. **Subtree grouping test:** Build a tree where two tips from the same
   subtree are on the wrong side of a constraint. Verify that
   `impose_constraint` clips the shared subtree once rather than
   processing tips individually (check via move count or topology).

3. **Round-trip test:** Start from a valid constrained tree. Apply
   `random_nni_perturb`. Call `impose_constraint`. Verify constraints
   satisfied. Verify the topology is not identical to the original
   (perturbation diversity is preserved).

4. **Multiple constraints test:** Tree violates two nested constraints.
   Verify both are repaired in a single `impose_constraint` call.

5. **Integration test:** Run `MaximizeParsimony` with constraints +
   `nniPerturbCycles > 0`. Verify output trees satisfy constraints.
   (Currently this combination is impossible because NNI perturb is
   disabled.)

6. **Fuse test:** Run constrained search, verify fuse doesn't discard
   all exchanges (check that pool diversity is maintained).

7. **Bail-out test:** Heavily scramble a constrained tree (fraction ~1.0).
   Verify `impose_constraint` falls back to `random_constrained_tree()`
   and still produces a valid tree.

## Prerequisite refactoring

Extract the node tip-bitmask computation from `map_constraint_nodes()`
(ts_constraint.cpp:130–151) into a standalone helper:

```cpp
// Compute per-node subtree tip bitmasks via postorder traversal.
// Returns array of size n_node * n_words.
std::vector<uint64_t> compute_node_tips(const TreeState& tree, int n_words);
```

Then `map_constraint_nodes()` becomes:
```cpp
void map_constraint_nodes(const TreeState& tree, ConstraintData& cd) {
  if (!cd.active) return;
  auto node_tips = compute_node_tips(tree, cd.n_words);
  // ... search for exact matches (existing code) ...
}
```

This is a pure extraction — no behaviour change, no new tests needed.

## Estimated scope

| Component | Lines | Complexity | Notes |
|-----------|-------|------------|-------|
| `compute_node_tips` helper | ~25 | Low | Extract from `map_constraint_nodes` |
| `impose_constraint` | ~100 | Medium | Min-sym-diff selection, subtree grouping, SPR clip/regraft |
| NNI perturbation integration | ~10 | Low | Remove gate, add call |
| Fuse integration | ~10 | Low | Replace revert with repair |
| Tests | ~100 | Low | 7 test cases |
| **Total** | **~245** | | |

No new topology-manipulation primitives needed (reuses `spr_clip`/`spr_regraft`).

## Complexity

- **Per violated split:** O(n × n_words) for node_tips rebuild +
  O(n_internal × n_words) for candidate search + O(k) SPR operations
  where k = number of maximal misplaced subtrees (typically 1–3)
- **Total:** O(v × n × n_words) where v = number of violated splits
  (usually 0–2 after NNI perturbation)
- **Common case (no violations):** O(n_internal × n_splits × n_words)
  for `map_constraint_nodes()`, then immediate return.  This is the same
  cost as `update_constraint()` which already runs after every accepted
  TBR/SPR move, so impose_constraint adds no new overhead in the
  no-violation case.

## Risks

- **SPR clip/regraft could corrupt tree state** if edge cases aren't
  handled (root-adjacent clips, invalid postorder between moves).
  Mitigated by: (a) guarding against root-child clips, (b) rebuilding
  postorder after each split's moves, (c) assertions on small trees.
- **Repair might undo most of the perturbation** if constraints are
  very tight (many splits, large clade coverage). In that case the
  NNI perturbation + fixup is little better than the current "skip
  entirely" approach. Mitigated by bail-out threshold + monitoring
  via benchmarking.
- **Heavy violation = too many moves.** If total moves exceed n_tip/4,
  we fall back to `random_constrained_tree()` rather than performing
  a large number of SPRs that would effectively rebuild the tree.
