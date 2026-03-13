# Step 1: TBR Implementation Spec

## Status: COMPLETE — all 29 tests passing

Files created: src/ts_tbr.h, src/ts_tbr.cpp, tests/testthat/test-ts-tbr-search.R
Files modified: src/ts_rcpp.cpp, src/TreeSearch-init.c

### Key bugs fixed during implementation:
1. **Segfault from tip sibling**: `fitch_incremental_downpass` crashes if
   started at a tip node. Fixed by starting at nz (grandparent) instead of
   ns (sibling), since nz is always internal.
2. **Wrong divided_length**: Must subtract nx's local_cost from the original
   score, since nx is removed from the divided tree.
3. **Topology corruption on undo**: The original approach of unclip/re-clip/
   reroot/regraft was fragile. Replaced with save_topology/restore_topology
   snapshots for clean undo.

## What TBR is

TBR = SPR + subtree rerooting. After clipping a subtree, SPR tries
regrafting it at every edge in the main tree. TBR additionally tries
every possible rerooting of the clipped subtree before each regraft.

## Current codebase (read these files for reference)

- `src/ts_tree.h/cpp` — TreeState with parent/left/right arrays, SPR clip/regraft/unclip
- `src/ts_fitch.h/cpp` — Fitch scoring: full downpass/uppass, incremental, indirect length
- `src/ts_search.h/cpp` — NNI and SPR search loops (first-improvement hill-climbing)
- `src/ts_rcpp.cpp` — R bridge: ts_fitch_score, ts_nni_search, ts_spr_search
- `src/ts_data.h/cpp` — DataSet, CharBlock, bit-packed character data

### Key conventions
- Tips: 0..n_tip-1, Internal: n_tip..2*n_tip-2, Root: n_tip
- parent[root] = root
- left[]/right[] indexed by (node - n_tip)
- State arrays: node-major, prelim[node * total_words + offset]
- SPR clip saves state to clip_undo_stack, unclip restores

## Algorithm: TBR search with indirect calculation

### Per-clip workflow

```
for each clip_node (skip children of root):
  1. spr_clip(clip_node)           — detach subtree
  2. build_postorder()              — main tree only
  3. incremental_downpass(sibling)  — reoptimize main tree
  4. incremental_uppass(sibling)    — get final states in main tree
  5. divided_length = best_score + delta

  6. SPR candidates (no rerooting):
     clip_prelim = &prelim[clip_node * total_words]
     For each edge (above, below) in main tree:
       extra = fitch_indirect_length(clip_prelim, tree, ds, above, below)
       candidate = divided_length + extra
       If candidate improves: verify with full rescore, accept if confirmed

  7. TBR candidates (with rerooting, only if clip_node is internal):
     Compute from_above[] for all nodes in clipped subtree
     For each edge (P, C) in clipped subtree:
       virtual_prelim = fitch_join(from_above[C], prelim[C])
       For each edge (above, below) in main tree:
         extra = fitch_indirect_length(virtual_prelim, ..., above, below)
         ...same acceptance logic...

  8. If accepted: keep new topology, full rescore for consistency
     If not: spr_unclip() to restore
```

### from_above computation (key to TBR efficiency)

For each node C in the clipped subtree, from_above[C] represents the
Fitch preliminary state looking from C's parent toward the rest of the
subtree (excluding C's own subtree). This lets us compute virtual root
states for any rerooting in O(total_words) per edge.

```
Subtree root = clip_node (= Nm)
Nm has children L, R:
  from_above[L] = prelim[R]
  from_above[R] = prelim[L]

For any other node C with parent P, sibling S (preorder from Nm):
  from_above[C] = fitch_join(from_above[P], prelim[S])

Virtual root at edge (P, C):
  virtual_prelim = fitch_join(from_above[C], prelim[C])
```

The fitch_join operation: intersection if non-empty, else union (per char).
This is identical to fitch_downpass_node but we only need the state set,
not the step count.

### Verification strategy

The indirect calculation gives the exact tree length for standard Fitch EW.
However, for safety (and future compatibility with IW/inapplicables where
indirect may be approximate), we verify the best candidate with full_rescore
before accepting. This adds O(n) per accepted move — negligible since
acceptance is rare relative to evaluation.

### First-improvement strategy

Like SPR: shuffle clip candidates each pass, accept first improving move
(break to next pass). For TBR, we find the BEST candidate across all
(reroot, regraft) combinations for each clip, then verify that single best.
This avoids expensive full rescores on every candidate.

## TBRParams interface (flexible for ratchet/drifting)

```cpp
struct TBRParams {
    bool accept_equal = false;     // accept Δ=0 moves?
    int max_accepted_changes = 0;  // 0 = no limit (run to convergence)
    int max_hits = 1;              // equal-score hits before stopping pass
};

struct TBRResult {
    int best_score;
    int n_accepted;
    int n_evaluated;
    bool converged;
};
```

Start minimal — only the fields needed for standard search. Ratchet/drifting
fields (afd_limit, rfd_limit, block_weights, etc.) will be added by their
respective agents. The interface is designed so those additions are
backward-compatible (new fields with defaults).

## Files to create/modify

### NEW: src/ts_tbr.h
- TBRParams, TBRResult structs
- tbr_search() declaration
- Helper declarations: compute_from_above(), collect_subtree_edges()

### NEW: src/ts_tbr.cpp
- tbr_search() implementation
- compute_from_above(): preorder traversal of clipped subtree, computing
  from_above states using fitch_join
- collect_subtree_edges(): DFS from clip_node, collecting (parent, child) pairs
- fitch_join_state(): compute Fitch join of two state sets (reusable helper)
- Full search loop with first-improvement + indirect filtering + verification

### MODIFY: src/ts_rcpp.cpp
- Add ts_tbr_search() export function, same pattern as ts_spr_search

### NEW: tests/testthat/test-ts-tbr-search.R
Tests:
1. TBR finds same score as SPR on simple trees (SPR ⊂ TBR)
2. TBR finds better score than SPR on trees where SPR gets stuck
3. Score matches TreeLength() on the result tree
4. Tree topology is valid after search
5. Flexible params: accept_equal works, max_accepted_changes stops early
6. Various tree sizes: small (7 tips), medium (20), larger (50+)

## Implementation notes

- Store from_above in a temporary vector<uint64_t> sized n_node * total_words.
  Only the clipped subtree's entries are populated; others are unused.
- For collecting subtree edges: DFS from clip_node through left[]/right[].
  These arrays still point to valid children even after clip (clip only
  modifies the connection above clip_node's parent).
- The indirect calculation uses tree.final_[] which has been updated by
  the incremental uppass for the main tree. Clipped subtree nodes retain
  their original prelim[] from the pre-clip full downpass.
- After accepting a move: physically apply it (spr_regraft + optional
  subtree reroot), then full_rescore to get all states consistent.

## Subtree rerooting (physical, for applying accepted TBR moves)

When accepting a TBR move with reroot at edge (P, C) in the clipped subtree:
We need to physically reroot before regrafting. The reroot changes which
node of the subtree connects to Nx (the clip parent).

Implementation: reverse parent-child links along the path from clip_node
down to P. After reversal, the subtree is rooted such that P's children
are: C and the reversed chain back to clip_node's original root.

For the first implementation, a simpler approach: don't physically reroot.
Instead, after finding the best (reroot_edge, regraft_edge) via indirect
calc, apply the full TBR move by:
1. spr_unclip() to restore original tree
2. Apply the TBR as a single topology modification
3. Full rescore

Actually simplest: just apply spr_regraft as usual, then reroot the
subtree in-place. OR: since we verify with full_rescore anyway, we can
directly construct the new topology.

DECISION: For the first implementation, use a helper `tbr_apply_move()`
that takes (clip_node, reroot_edge, regraft_edge) and constructs the
new topology from scratch, then full_rescore. This is clean and correct.
The helper reverses parent-child links in the subtree along the path to
the reroot point, then does the standard regraft.
