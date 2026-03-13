# Plan: Phase 4 — Inapplicable Character Scoring

## Goal

Implement the Brazeau et al. (2019) three-pass inapplicable scoring algorithm
using bit-parallel operations, integrated into the existing SPR search.
Standard Fitch blocks (no inapplicable data) continue using the one-pass
algorithm; inapplicable blocks use the three-pass algorithm. Verify against
morphy on all `inapplicable.phyData` datasets.

## Background

Morphy's per-character implementation (in `src/fitch.c`) processes characters
one at a time with branches per character. Our bit-packed approach processes
64 characters per block using branchless mask operations. The plan already
contains verified mask algebra for all three passes.

### Morphy's state representation
- Bit 0 = inapplicable (NA); bits 1..k = applicable states
- `ISAPPLIC = ~NA = all bits except bit 0`
- Our representation is identical: state 0 = NA in `has_inapplicable` blocks

### Morphy's arrays per node (per character)
- `downpass1` → first downpass preliminary states (NA-aware)
- `uppass1` → first uppass result (applicability determination)
- `downpass2` → second downpass states (corrected, used for scoring)
- `subtree_actives` → applicable states present anywhere in the subtree below

## Implementation

### Step 1: Add new state arrays to TreeState

In `ts_tree.h`, add two new arrays to `TreeState`:

```cpp
// Second downpass states (only computed for inapplicable blocks)
std::vector<uint64_t> down2;

// Subtree actives: applicable states present in the subtree below.
// Only word indices for applicable states (s > 0) are meaningful.
std::vector<uint64_t> subtree_actives;
```

Both use the same layout as `prelim` — `node * total_words + offset`.
Allocate in `init_from_edge()` and clear in `reset_states()`.

For tips, `subtree_actives` is initialised from `tip_states` with the NA
word (state 0) zeroed out. This happens in `load_tip_states()` or
`reset_states()`.

### Step 2: Implement the three passes in ts_fitch.cpp

Add three new functions for inapplicable blocks. All operate on one block
at a time within the main scoring loop.

#### 2a. First downpass (NA-aware)

Replace standard Fitch for inapplicable blocks. No step counting.
From the plan's verified mask algebra:

```
I[s] = L[s] & R[s]
I_app = I[1] | ... | I[k-1]
L_app = L[1] | ... | L[k-1]
R_app = R[1] | ... | R[k-1]
both_app = L_app & R_app
case_keep = I_app | (I[0] & ~I_app & ~both_app)
case_strip = ~I[0] & ~I_app & both_app
result[s>0] = (I[s] & case_keep) | (union[s] & ~case_keep)
result[0]   = (I[0] & case_keep) | (union[0] & ~case_keep & ~case_strip)
```

Also compute `subtree_actives`:
```
subtree_actives[s>0] = left_actives[s] | right_actives[s]
subtree_actives[0]   = 0  // NA is never "active"
```

#### 2b. First uppass (applicability propagation)

Determines whether each node is applicable or inapplicable, propagating
from root downward. Root's uppass = its downpass1 states (same as standard).

From the plan's verified mask algebra (for internal nodes):
```
npre_has_NA  = N[0]      (= prelim has NA bit set)
npre_has_app = N[1] | ... | N[k-1]
anc_app      = A[1] | ... | A[k-1]
anc_is_NA    = A[0] & ~anc_app    (ancestor is purely NA)

case_passthrough = ~npre_has_NA
case_strip_NA    = npre_has_NA & npre_has_app & ~anc_is_NA
case_children    = npre_has_NA & ~npre_has_app & ~anc_is_NA & children_app
case_force_NA    = active_mask & ~case_passthrough & ~case_strip_NA & ~case_children

result[0]   = (N[0] & case_passthrough) | case_force_NA
result[s>0] = (N[s] & (case_passthrough | case_strip_NA))
            | ((L[s]|R[s]) & case_children)
```

Where `children_app = (L[1]|R[1]) | ... | (L[k-1]|R[k-1])`, computed from
the first downpass prelim states of the children.

For tips: if tip has NA and ancestor says applicable, strip NA. If ancestor
is NA, force NA. (Same logic as internal nodes but without children lookup.)

#### 2c. Second downpass (corrected scoring)

The scoring pass. Uses the first uppass result (`final_`) to determine
applicability, then counts steps only at applicable nodes. From verified
mask algebra, with reference to morphy's `mpl_NA_fitch_second_downpass`:

```
// ss_app: is this node applicable (from first uppass)?
ss_app = final_[1] | ... | final_[k-1]

// Children's second-downpass states
I_app = (L2[1] & R2[1]) | ... | (L2[k-1] & R2[k-1])
L2_app = L2[1] | ... | L2[k-1]
R2_app = R2[1] | ... | R2[k-1]
both_app = L2_app & R2_app

if node is applicable (ss_app):
  if applicable intersection exists (I_app):
    down2[s>0] = (L2[s] & R2[s]) & has_I_app   (strip NA from intersection)
    down2[0]   = 0
  else:
    down2[s>0] = (L2[s] | R2[s]) & applicable_mask
    down2[0]   = 0
    // Count step if: both children applicable, or both have active descendants
    needs_step = (both_app | (l_actives_any & r_actives_any)) & ~I_app & ss_app
else:
  // Node is inapplicable
  down2 = final_ states (i.e., NA)
  // Count step if both children have applicable descendants
  needs_step = l_actives_any & r_actives_any & ~ss_app
```

This is the most complex pass. I need to carefully verify the step-counting
conditions against morphy's code:

From morphy's `mpl_NA_fitch_second_downpass` (line 453-494):
```c
if (nifin[j] & ISAPPLIC) {        // node is applicable
  if ((temp = (left[j] & right[j]))) {
    if (temp & ISAPPLIC) {
      npre[j] = temp & ISAPPLIC;   // applicable intersection
    } else {
      npre[j] = temp;              // NA-only intersection (shouldn't happen if applicable?)
    }
  } else {                         // no intersection
    npre[j] = (left[j] | right[j]) & ISAPPLIC;
    if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
      steps += weights[i];         // both children have applicable states → step
    } else if (lacts[j] && racts[j]) {
      steps += weights[i];         // both children have active descendants → step
    }
  }
} else {                           // node is inapplicable
  npre[j] = nifin[j];             // keep NA
  if (lacts[j] && racts[j]) {
    steps += weights[i];           // both children have active descendants → step
  }
}
stacts[j] = (lacts[j] | racts[j]) & ISAPPLIC;
```

So the step-counting conditions are:
1. **Applicable node, no applicable intersection, both children applicable**: step
2. **Applicable node, no applicable intersection, both children have active descendants**: step
3. **Inapplicable node, both children have active descendants**: step

Translating to bit-parallel:
```
// Per-character masks
ss_app = final_[1] | ... | final_[k-1]      // node applicable?
I = L2[s] & R2[s] for each s
I_app = I[1] | ... | I[k-1]                  // applicable intersection?
L2_app = L2[1] | ... | L2[k-1]               // left child has applicable?
R2_app = R2[1] | ... | R2[k-1]               // right child has applicable?
l_act_any = left_actives[1] | ... | left_actives[k-1]  // left has active descendants?
r_act_any = right_actives[1] | ... | right_actives[k-1]

// Step conditions:
step_case1 = ss_app & ~I_app & L2_app & R2_app
step_case2 = ss_app & ~I_app & ~(L2_app & R2_app) & l_act_any & r_act_any
step_case3 = ~ss_app & l_act_any & r_act_any

needs_step = (step_case1 | step_case2 | step_case3) & active_mask
// Simplifies to:
needs_step = (~I_app & ((ss_app & (L2_app & R2_app | l_act_any & r_act_any))
           | (~ss_app & l_act_any & r_act_any))) & active_mask
// Which further simplifies (since L2_app implies l_act_any for the same subtree):
needs_step = (~I_app & l_act_any & r_act_any
           & (ss_app | active_mask)) & active_mask
// Wait, that's not right because ss_app & L2_app & R2_app is a subset of
// ss_app & l_act_any & r_act_any. Let me re-derive...
```

Actually, cases 1 and 2 together are equivalent to:
`ss_app & ~I_app & (L2_app & R2_app | l_act_any & r_act_any)`

But `L2_app & R2_app` implies `l_act_any & r_act_any` (if a child has
applicable states in its down2, it certainly has active descendants).
So cases 1+2 simplify to: `ss_app & ~I_app & l_act_any & r_act_any`.

Combined with case 3: `~ss_app & l_act_any & r_act_any`.

So the complete step mask is simply:
```
needs_step = ~I_app & l_act_any & r_act_any & active_mask
```

Wait — that ignores `ss_app` entirely. Let me check: if `ss_app` is true
and `I_app` is true, no step. If `ss_app` is true and `I_app` is false and
`l_act_any & r_act_any`, step. If `ss_app` is false and `l_act_any &
r_act_any`, step.

So: needs_step = `~I_app & l_act_any & r_act_any` regardless of `ss_app`.

But we also need `I_app` to be false. If the node is applicable and has
an applicable intersection, no step regardless. If the node is inapplicable,
`I_app` would be based on children's down2 states... which may be NA. Let me
think: if the node is inapplicable, then `down2 = NA`, and the children's
down2 states feed into the intersection. But actually, `I_app` is computed
from the children's states, not the node's. So `I_app` could be non-zero
even when the node is inapplicable.

Hmm, but morphy checks `I_app` (via `temp = left & right`) only in the
applicable branch. In the inapplicable branch, it just checks `lacts && racts`.
So the full logic is:

```
needs_step = (ss_app & ~I_app & l_act_any & r_act_any)
           | (~ss_app & l_act_any & r_act_any)
           = l_act_any & r_act_any & (~I_app | ~ss_app)
           = l_act_any & r_act_any & ~(I_app & ss_app)
```

This is: step if both children have active descendants, AND NOT (node is
applicable AND children have applicable intersection).

That's a clean formula for the bit-parallel version.

For the down2 state assignment:
```
if ss_app:
  if I_app:
    down2[s>0] = (I[s] & I_app_mask) for applicable states
    down2[0] = 0
  else:
    down2[s>0] = (L2[s] | R2[s]) & applicable (union of applicable)
    down2[0] = 0
else:
  down2 = final_ (which is NA)
```

And subtree_actives update: `stacts[s] = (l_acts[s] | r_acts[s]) & ISAPPLIC`
(same as in the first downpass; just re-propagate from the corrected states).

### Step 3: Integrate into fitch_score()

Modify `fitch_score()` (or create a new `fitch_na_score()`) to:

1. For standard blocks: standard downpass (as now)
2. For inapplicable blocks: first downpass (NA-aware, no step counting)
3. Standard uppass for standard blocks
4. First uppass for inapplicable blocks (applicability propagation)
5. For inapplicable blocks: second downpass (step counting)
6. Total score = standard block steps + inapplicable block steps

The traversal order: all blocks share the same postorder, so we can
interleave standard and NA blocks within a single traversal.

### Step 4: Integrate into SPR search

The SPR search currently does:
1. Clip → full two-pass on divided tree → indirect calc for screening → verify

For inapplicable blocks, the indirect calculation is not straightforward
(the three-pass structure means changes propagate globally). Strategy:

- **Standard blocks**: Continue using indirect calc (union-based) for
  screening, with verification rescore
- **Inapplicable blocks**: Always contribute to the verification rescore;
  don't attempt indirect screening for these blocks

In practice: the indirect calc screens based on standard blocks only. For
any candidate that passes the standard-block screen, do a full three-pass
rescore of the whole tree (both standard and inapplicable blocks) to get
the actual score. This is conservative but correct.

Simpler alternative (for now): since we already do verification rescores
for all accepted candidates, the indirect calc just needs to not reject
moves that are actually improvements. Since the indirect calc only
considers standard blocks, its estimate is a lower bound on the extra
steps (inapplicable blocks may add more). So using the standard-block
indirect calc as a screen is valid — it won't reject genuinely good moves.

Actually, this is wrong. The indirect calc computes `divided_length + extra`
where `divided_length` includes ALL blocks (standard + inapplicable). But
`extra` only accounts for standard blocks. So the candidate score would
underestimate the actual score for inapplicable chars, which means it might
accept candidates that are actually worse. This is the opposite of
conservative.

**Correct approach**: Compute indirect extra for standard blocks only. For
inapplicable blocks, add a conservative estimate (e.g., 0 — assume no extra
steps from inapplicable chars). Then the candidate score ≤ actual score,
making the screen conservative. But the divided_length already includes the
inapplicable block contributions, so:

```
candidate_score = divided_length + extra_standard_blocks + 0_for_inapp_blocks
```

This is ≤ actual score (since inapplicable extra ≥ 0). Conservative: may
let through candidates that are actually worse, but verification rescore
catches them. Better than screening them out incorrectly.

Wait, I just said "may let through candidates that are actually worse" — but
that's what verification is for. The current search already does verification
for all candidates that pass screening. So this approach is fine.

### Step 5: Testing

1. **Score test**: For each dataset in `inapplicable.phyData`, score the
   first reference tree with `ts_fitch_score()` and compare to
   `TreeLength(tree, dataset)` (morphy). Must match exactly.

2. **Random tree test**: Generate 10 random trees per dataset, score with
   both methods, compare.

3. **SPR search test** (if time): Run SPR search on inapplicable data,
   verify result tree scores match morphy.

## Files to modify

| File | Changes |
|------|---------|
| `src/ts_tree.h` | Add `down2`, `subtree_actives` arrays |
| `src/ts_tree.cpp` | Allocate new arrays, initialise in `reset_states()` |
| `src/ts_fitch.h` | Declare new NA functions |
| `src/ts_fitch.cpp` | Implement three-pass algorithm for NA blocks |
| `src/ts_search.cpp` | Adjust divided_length/extra computation for mixed blocks |
| `src/ts_rcpp.cpp` | No changes needed (existing `ts_fitch_score` will use the updated scoring) |

## Risks

1. **Mask algebra bugs**: The plan says verified in R, but the bit-packed
   implementation could have off-by-one or masking errors. Mitigation:
   test against morphy on diverse datasets.

2. **Performance of full three-pass**: Each inapplicable block requires
   3 passes instead of 1. But the bit-parallel approach still processes
   64 characters per operation, so it should be fast. The main concern is
   that the SPR search does more verification rescores (since indirect
   screening is less effective for inapplicable data).

3. **Edge cases**: All-NA characters, single-taxon subtrees with NA,
   polymorphic terminals with mixed applicable/inapplicable states. These
   are covered by the `inapplicable.phyData` test datasets.
