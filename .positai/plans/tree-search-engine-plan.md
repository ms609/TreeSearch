# Plan: High-Performance C++ Tree Search Engine for TreeSearch

## Motivation

Currently, `MaximizeParsimony()` coordinates tree search in R, with C code
handling only the Fitch scoring (via morphy). Tree rearrangements are generated
by Rcpp but the search loop, move acceptance, and ratchet logic all live in R.
Implied weighting requires one morphy object per character pattern, adding
overhead. We want to move the entire search loop into C++, keeping the existing
R infrastructure as a correctness reference.

## Goals

1. **All-C++ heuristic search**: TBR/SPR/NNI rearrangement + scoring + move
   acceptance in a single compiled loop, with no R callbacks per move.
2. **Inapplicable support**: Implement the Brazeau et al. (2019) algorithm
   in our own data structures, verified against morphy.
3. **Implied weights**: Per-character scoring with early termination.
4. **Partial rescoring**: Store per-node state sets so that after a
   rearrangement, only affected nodes need rescoring.
5. **Heuristic strategies**: Ratchet, sectorial search, and adaptive search
   parameters — all in C++, exposed to R for configuration.

---

## Data layout: bit-packed state sets

### Packing across characters

Pack up to 64 characters of the same number of states into a single block of
`uint64_t` words. Bit *i* of word *j* represents "character *i* can have
state *j*". The Fitch downpass for 64 characters at once:

```cpp
// For k-state characters packed into k words:
uint64_t any_intersect = 0;
for (int s = 0; s < k; ++s) {
    intersect[s] = left[s] & right[s];
    any_intersect |= intersect[s];
}
uint64_t needs_union = ~any_intersect;
// Mask to only count characters actually present in this block:
needs_union &= active_mask;
steps += popcount(needs_union);  // 64 characters in one popcount
for (int s = 0; s < k; ++s) {
    result[s] = (intersect[s] & any_intersect)
              | ((left[s] | right[s]) & needs_union);
}
```

### Why this wins even for local_reopt and IW

Morphy's per-character `local_reopt` loop processes characters one at a time
with a branch each. Its early termination (bail out when `steps > maxlen`)
saves work when the cutoff is hit quickly. But with bit packing, one block
operation processes 64 characters for roughly the cost of morphy's one. Even
if we process all characters before checking the cutoff, 4 block iterations
(for 200 characters) is cheaper than morphy processing 30 characters with
a branch per character. The 64x throughput advantage overwhelms any early
termination benefit for realistic dataset sizes (50–500 characters).

For `local_reopt` specifically:
```cpp
// Bit-packed local_reopt: 64 chars per block
for (int b = 0; b < n_blocks; ++b) {
    for (int s = 0; s < block_n_states[b]; ++s) {
        tgt[s] = tgt1_up[b * max_k + s] | tgt2_up[b * max_k + s];
    }
    uint64_t no_intersect = active_mask[b];
    for (int s = 0; s < block_n_states[b]; ++s) {
        no_intersect &= ~(src_down[b * max_k + s] & tgt[s]);
    }
    steps += block_weight[b] * popcount(no_intersect);
    if (steps > maxlen) return steps;  // block-level early exit
}
```

For EW, `block_weight` is uniform so `popcount` directly gives the count.
The block-level early exit is coarser than per-character but still effective,
and the cost per check is negligible compared to the 64x throughput.

### Per-character step counts for IW

IW requires per-character extra-step counts to compute `Σ e_i/(k + e_i)`.
We extract these from the bit-packed representation using a scatter-add
during the downpass:

```cpp
// At each internal node, after computing needs_union for block b:
uint64_t mask = needs_union;
while (mask) {
    int bit = ctz(mask);          // index of lowest set bit
    char_steps[block_start + bit] += 1;
    mask &= mask - 1;             // clear lowest set bit
}
```

This only iterates over characters that actually have a step at this node
(typically a minority), so the scatter-add is cheap relative to the bulk
Fitch operations. The `char_steps` array is maintained alongside the
bit-packed state sets and reset at the start of each full score.

For the IW score itself, evaluate characters in order of decreasing expected
contribution and bail out early if the partial score exceeds the target.

### Handling inapplicable data — verified mask algebra

The inapplicable bit occupies state 0 in each block. The Brazeau et al.
three-pass algorithm's per-character branching reduces to branchless mask
operations. This has been **verified in R** against morphy's per-character
implementation on all edge cases.

State 0 = inapplicable (NA). States 1..k-1 = applicable.

#### First downpass (NA-aware)

```cpp
// Intersection per state
I[s] = L[s] & R[s];                         // for each state s
I_app = I[1] | I[2] | ... | I[k-1];         // applicable intersection?
L_app = L[1] | ... | L[k-1];                // left has applicable?
R_app = R[1] | ... | R[k-1];                // right has applicable?
both_app = L_app & R_app;

// Key mask: which characters keep the intersection result?
// Keep if: has applicable intersection, OR intersection is {NA} + not both applicable
case_keep = I_app | (I[0] & ~I_app & ~both_app);
// Strip NA from union: empty intersection + both applicable
case_strip = ~I[0] & ~I_app & both_app;

// Results (64 characters in parallel):
result[s] = (I[s] & case_keep) | (union[s] & ~case_keep);   // s > 0
result[0] = (I[0] & case_keep) | (union[0] & ~case_keep & ~case_strip);
```

#### First uppass (applicability propagation)

```cpp
npre_has_NA  = N[0];
npre_has_app = N[1] | ... | N[k-1];
anc_app      = A[1] | ... | A[k-1];
anc_is_NA    = A[0] & ~anc_app;
children_app = (L[1]|R[1]) | ... | (L[k-1]|R[k-1]);

case_passthrough = ~npre_has_NA;
case_strip_NA    = npre_has_NA & npre_has_app & ~anc_is_NA;
case_children    = npre_has_NA & ~npre_has_app & ~anc_is_NA & children_app;
case_force_NA    = everything else;

result[0] = (N[0] & case_passthrough) | case_force_NA;
result[s] = (N[s] & (case_passthrough | case_strip_NA))
          | ((L[s]|R[s]) & case_children);                  // s > 0
```

#### Second downpass (corrected scoring)

```cpp
ss_app = SS[1] | ... | SS[k-1];     // uppass1 says node is applicable?
I_app  = (L[1]&R[1]) | ...;         // applicable intersection of children?
both_app = L_app & R_app;
union_has_NA = L[0] | R[0];

// Step counting mask: applicable node, no applicable intersection,
// both children applicable, union doesn't contain NA
needs_step = ss_app & ~I_app & both_app & ~union_has_NA;

result[s] = (applic_intersect[s] & has_I_app) | (union[s] & ~has_I_app);
result[0] = ~ss_app | (union[0] & ~has_I_app);  // inapplicable → NA
```

Each pass is ~4k + 10 operations per block, processing 64 characters.
Compare with morphy: ~5 operations × 64 characters × branch misprediction.

#### Partial rescoring with inapplicables

The three-pass structure means a local rearrangement can propagate globally:
the first downpass goes up, the uppass comes back down (potentially to all
nodes), the second downpass goes up again. Morphy's own `local_reopt` for
inapplicable characters (`mpl_fitch_NA_local_reopt`) doesn't attempt partial
rescoring — it marks all characters as needing update.

**Strategy**: For inapplicable blocks, always do a full three-pass rescore
when a move passes `local_reopt`. The `local_reopt` itself is still cheap
(one block-pass per candidate). For standard Fitch blocks only, partial
rescoring (rescore only affected path, stop when state sets stabilize) is
safe.

This is acceptable because `local_reopt` rejects the vast majority of
candidates before any full rescore is needed.

### Character weight handling

phyDat compresses duplicate character patterns and stores a `weight` vector
(pattern frequency). Typical values are 1–8. For bit-packed scoring, we want
uniform weight per block so that `popcount` directly gives the score.

**Decision: expand patterns by weight.** Duplicate each pattern `weight`
times so all characters have weight 1. Checked on all 30 `inapplicable.phyData`
datasets: expansion is minimal (e.g. 360 patterns → 363 expanded characters
for the largest dataset). Maximum expansion across all datasets: ~6 blocks
of 64. This eliminates weight complications from the entire scoring pipeline.

For IW, each expanded character maps back to its original pattern index
(stored in a lookup array) so that per-character step counts can be
aggregated back to per-pattern counts for the IW fit calculation.

### Loading from phyDat

phyDat stores:
- Per-taxon integer vectors indexing into `allLevels`
- `contrast` matrix: rows = allLevels tokens, columns = states (including `-`)
- `levels`: the state labels (e.g. `"-" "0" "1" "2" "3"`)
- `index`: mapping from actual characters to unique patterns
- `weight`: pattern frequencies

Conversion to bit-packed form:
1. For each taxon × pattern, look up the contrast row → get the set of
   possible states → pack into a uint64_t word per state.
2. Expand patterns by weight.
3. Partition into blocks: first all characters with inapplicable data
   (has `-` state), then all without. Within each partition, group by
   number of applicable states. Pack up to 64 characters per block.

### Block structure

```cpp
struct CharBlock {
    int n_chars;              // characters in this block (≤ 64)
    int n_states;             // including NA if has_inapplicable
    bool has_inapplicable;    // needs the NA-aware algorithm
    uint64_t active_mask;     // bits 0..n_chars-1 set
    // For IW: map back to original pattern index
    int pattern_index[64];    // which original pattern each char came from
};

struct DataSet {
    int n_tips;
    int n_blocks;
    std::vector<CharBlock> blocks;
    int total_words;          // sum of n_states across all blocks
    // Tip state data, flattened:
    // tip_states[tip * total_words + word_offset]
    std::vector<uint64_t> tip_states;
    // IW metadata (per original pattern):
    int n_patterns;
    std::vector<int> min_steps;     // minimum steps per pattern
    std::vector<int> pattern_freq;  // original weight (for score reporting)
};
```

---

## Tree representation: `TreeState`

### Topology

Flat arrays, as in `build_postorder.h`:
- `parent[0..2n-2]`, `left[0..n-2]`, `right[0..n-2]` (internal nodes
  indexed from n_tip)
- `n_tip` tips indexed `0..n_tip-1`, internal nodes `n_tip..2*n_tip-2`
- Root at index `n_tip`; `parent[root] = root`

### Per-node state sets

All state data in a flat allocation for cache locality:
```cpp
// downpass1[node * total_words + word_offset]
std::vector<uint64_t> downpass1;
std::vector<uint64_t> uppass1;
std::vector<uint64_t> downpass2;  // only for inapplicable partitions
// Temp buffers for tentative moves:
std::vector<uint64_t> temp_down1;
std::vector<uint64_t> temp_up1;
```

For 100 taxa, 200 characters in ~4 blocks of 2–4 states (~12 words per
node): each state array is ~200 nodes × 12 × 8 bytes ≈ 19 KB. All arrays
together fit comfortably in L1/L2 cache.

### Per-character step counters (for IW)

```cpp
// char_steps[char_idx] — extra steps for character, accumulated during scoring
std::vector<int16_t> char_steps;
```

Reset and populated during each full score. For partial rescores, we
subtract the old contribution of changed nodes and add the new.

### Topology deduplication

The Tromp ID (`tree_num_t` in TreeTools) supports up to ~51 leaves. For
larger trees, use **split-based hashing**:

- Each internal edge defines a bipartition (split) of the taxa.
- Represent each split as a bitset of `ceil(n_tip / 64)` uint64_t words.
- Hash the sorted set of splits for the tree.
- Two trees are identical iff they have the same split set.
- Confirm identity on hash collision by comparing sorted split vectors.

Splits can be maintained incrementally: when a subtree is pruned/regrafted,
only splits involving edges in the affected region change.

For small trees (n_tip ≤ 51), Tromp ID is available as an exact fast path
via `TreeTools::tree_number.h`.

---

## Scoring

### Standard Fitch (EW, no inapplicables)

Bit-packed downpass over all internal nodes. One popcount per block per node.
Total EW score = sum of popcount results × block weight.

### Inapplicable Fitch (EW, with inapplicables)

Three passes, reimplementing the Brazeau et al. algorithm using bit-packed
mask operations (as sketched above):
1. First downpass: NA-aware intersection with mask-based case selection
2. First uppass: applicability propagation
3. Second downpass: corrected scoring

### `local_reopt` — the key optimization for SPR/TBR

After pruning a subtree, for each candidate regraft edge (u,v):
```cpp
for (int b = 0; b < n_blocks; ++b) {
    // Quick test: does subtree root intersect with regraft point?
    uint64_t no_hit = active_mask[b];
    for (int s = 0; s < block_n_states[b]; ++s) {
        no_hit &= ~(src_down[offset + s]
                   & (tgt1_up[offset + s] | tgt2_up[offset + s]));
    }
    steps += popcount(no_hit);
    if (steps > maxlen) return steps;
}
```

### Implied weights

After a full downpass populates `char_steps[]`:
```cpp
double iw_score = 0;
for (int i : char_order) {  // sorted by decreasing expected contribution
    int e = char_steps[i] - min_steps[i];
    if (e > 0) iw_score += weight[i] * e / (k + e);
    if (iw_score > target) return INFINITY;
}
```

---

## Indirect calculation and incremental optimization

*Reference: Goloboff (1996), "Methods for faster parsimony analysis."*

The core idea behind efficient TBR/SPR is **indirect tree length calculation**
(Goloboff 1994, 1996; independently used in Hennig86). When a subtree rooted
at Nm is clipped:

1. Compute final (uppass) state sets for both parts of the divided tree.
2. For each candidate regraft edge (A, D), derive the state set of a virtual
   root Y = union of final states of A and D. For non-additive characters
   this is `SA ∪ SD`; for additive characters, the smallest closed interval
   containing both.
3. Compute the length increase from joining the clipped subtree's basal node
   X to Y using the standard downpass rule — checking only 3 nodes (X, A, D)
   per character.

This is exact (proved in Goloboff 1996). For TBR, rerooting the clipped
subtree changes only the basal node's state set, which can be derived from
final states without re-optimizing the subtree.

### Incremental two-pass optimization (Shortcut C)

The bottleneck in indirect calculation is computing final states for the
divided tree after each clip. A naive two-pass over all characters at every
clip is expensive. Goloboff (1996) describes three shortcuts of increasing
sophistication; the third (**incremental two-pass**, Shortcut C) is the target:

1. **Incremental downpass**: Starting from Nz (ancestor of the clip point),
   walk rootward recomputing preliminary states. Stop when the new preliminary
   set equals the old one — changes cannot propagate further. Track the
   length delta by subtracting old local costs and adding new ones.

2. **Incremental uppass**: Reassign final states for all nodes visited in
   step (1), then propagate downward through any node whose ancestor's final
   states changed (compared to the whole tree). Also reassign the sister of
   the clip point (Ns) if the final states of its new ancestor (Nz) differ
   from those of its old ancestor (Nx). For the clipped subtree, the root
   node's preliminary states become its final states; propagate downward
   similarly.

3. **Short-circuit**: If the preliminary states of Ns equal those of Nm,
   skip steps (1) and (2) entirely — no states change.

Performance: Goloboff reports 0.15–0.20 o/t/c (optimizations per tree per
character), with search times ~45% of full two-pass reoptimization. The
clipping phase drops to 5–15% of total TBR time, making the rearrangement
phase (destination scoring) the dominant cost.

**Implementation target**: Use Shortcut C for Phases 2–3. Requires storing,
per node per character: (a) final state sets for the whole tree (buffered for
restoration on unclip), (b) preliminary state sets (shared storage for whole
and clipped tree, with the changed nodes buffered in a single variable during
the incremental downpass), and (c) per-character local costs for the whole
tree (for length delta during incremental downpass).

### Union construct — bulk rejection of destinations

Once the clipping phase is fast, the rearrangement phase (checking each
candidate destination) dominates. The **union construct** (Goloboff 1996)
provides a top-down pruning strategy:

- Precompute a **union set** at each internal node: the union of all states
  present in any descendant.
- For a candidate node N in the main subtree, form a "union construct" =
  union set of N's descendant + final states of N's ancestor.
- Check the clipped subtree's basal states against this union construct,
  exactly as in indirect calculation.
- If this produces a suboptimal length, **all** descendants of N can be
  skipped — they can only be worse.

Goloboff reports 45–57% of internal nodes rejected this way, with ~50% time
savings on congruent data. Less effective for highly incongruent data or very
suboptimal trees (e.g. early Wagner tree swapping).

An extension: check a union construct of the entire clipped subtree against
a destination node's union construct. If suboptimal, skip that destination
for all rerootings of the clipped subtree (TBR only).

**Implementation target**: Phase 3 or later; profile first to confirm the
rearrangement phase is the bottleneck before investing in union construct
infrastructure. Requires per-node union sets (incrementally maintainable
via an incremental downpass analogous to Shortcut C).

### Polymorphic terminal handling

Polymorphic terminals (multiple observed states for a character) need special
treatment during indirect calculation. The final state set of a polymorphic
terminal is its observed set (unmodified by the uppass). However, using those
raw states for destination scoring can under-count steps — the states
responsible for minimum length in the main subtree may differ from those
shared with the clipped subtree (Goloboff 1996, Fig. 2).

**Correction**: Treat a polymorphic terminal as if it were an internal node
whose two children both have the terminal's observed state set, then apply
the uppass rule. This need only be done when the polymorphic terminal is the
descendant node of a destination branch *and* it shares no states with the
ancestral node's final set — a rare case.

For missing-data terminals (`?`): treat their state set as identical to
their ancestor's final set for indirect calculation purposes.

---

## Zero-length branch collapsing

*Reference: Goloboff (1996), §"Collapsing The Trees."*

When zero-length branches are to be collapsed (critical for search
effectiveness — avoids filling the tree buffer with trivially different
dichotomous resolutions), final states from the divided tree enable fast
collapsing:

1. **Shortest-path shortcut** (approximate): If no node is supported by
   character change along the shortest path between the original and
   destination positions of the clipped subtree, assume the rearranged tree
   collapses to the same polytomy as the original. Skip the tree without
   further work. This shortcut is not exact — it can miss rare cases where
   ambiguous optimizations differ — but in practice it rarely causes errors.
   Note: this creates **asymmetric reachability**: swapping on tree A may
   find tree B, but swapping on B may not find A (Goloboff 1996, §"Collapsing
   and Islands"). This is a known limitation of heuristic search and does not
   affect correctness of accepted trees.

2. **Partial reoptimization for collapsing**: For characters where the
   clipped subtree's basal state set is identical to one of the destination
   branch's node state sets, the state assignments are unchanged and no
   reoptimization is needed. Only characters where the basal set differs from
   *both* destination nodes need reoptimization — typically 10–20% of
   characters (fewer for more taxa). This makes collapsing up to 10× faster
   than full reoptimization.

3. With incremental two-pass (Shortcut C) providing exact final states, the
   collapsed trees are guaranteed correct, with no need for a post-hoc
   verification pass.

**Implementation target**: Phase 6 (tree pool), as collapsing is needed for
effective pool management and deduplication. Consider making it optional
(configurable) since not all search modes require collapsed trees.

---

## Rearrangement operations

All operations work on `TreeState` in-place with undo support.

### SPR (primary workhorse)

```
prune(subtree_root):
    save parent/sibling linkage + state sets of affected nodes
    detach subtree, connect sibling to grandparent
    incremental two-pass (Shortcut C) on main subtree → final states for
        divided tree + length of divided tree
    if length_divided == best_length and not saving multiple trees:
        skip rearrangement phase entirely (no improvement possible)

    for each candidate regraft edge (A, D):
        indirect length calculation: check 3 nodes (X, A, D) per character
            → bail out as soon as accumulated length exceeds best_length
        if equal or better:
            accept; update state sets + topology

undo_prune():
    restore saved linkage and state sets (buffered during incremental pass)
```

The incremental two-pass (see "Indirect calculation" section above) makes
the clipping phase fast (~5–15% of TBR time). The rearrangement phase
(destination scoring) dominates and benefits from early length-check bailout
and, optionally, the union construct for bulk destination rejection.

### TBR

Like SPR, but after bisection, also considers rerooting the pruned subtree
at each internal edge before regrafting. Each rerooting changes the subtree's
root state set, derivable from final states without re-optimizing the
subtree (Goloboff 1994, 1996). This means TBR examines more rearrangements
per (expensive) clipping, making TBR more efficient *per clip* than SPR.

The union construct extension (checking clipped-subtree union against
destination union) can skip entire destination subtrees for *all* rerootings
at once.

### NNI

Swap two subtrees adjacent to an internal edge. Implemented for completeness
but **not used as a search strategy**. See note below on TBR dominance.

### Why TBR dominates SPR and NNI

With indirect calculation (Goloboff 1994, 1996) and union construct bounding,
TBR evaluates individual rearrangements *faster* than SPR, and SPR faster
than NNI. This is because more reinsertions per (expensive) clipping amortizes
the clipping cost. For large structured datasets, TBR completes full swap
cycles faster than SPR in absolute wall-clock time (Goloboff 2015). With
bounding, the advantage grows with dataset size.

Consequence: **TBR is the only rearrangement used during search.** SPR and
NNI are implemented (SPR as an intermediate step to TBR; NNI for testing)
but not exposed as search-time options. "Variable neighbourhood" strategies
(NNI→SPR→TBR escalation) are counterproductive in an efficient
implementation and are not implemented.

---

## Search strategies

### Main search loop

```
SearchEngine::run(starting_tree, dataset, params):
    tree = starting_tree
    full_score(tree)

    // Phase 1: Initial descent
    local_search(tree, TBR, maxHits=params.startHits)

    // Phase 2: Ratchet iterations
    for i in 1..ratchIter:
        perturbed_weights = perturb(dataset)
        local_search(tree, TBR, perturbed_weights)  // TBR even for escape phase
        restore_weights()
        local_search(tree, TBR, original_weights)
        record_if_best(tree)
        check_interrupt()  // R_CheckUserInterrupt()

    // Phase 3: Final refinement
    local_search(tree, TBR, maxHits=params.finalHits)

    return best_trees
```

### Parsimony ratchet

- Resample characters (with replacement) or perturb weights
- **Cheap ratchet option**: During the escape phase, use standard Fitch (no
  inapplicable correction, no IW) for speed. The perturbed landscape doesn't
  need accurate scoring — it just needs to push us somewhere different.

### Sectorial search

- Select a random subtree of ~6–25% of taxa
- Fix the rest of the tree; intensive rearrangement within the sector
- Accept if sector-local score improves
- Useful for large trees where full TBR is too expensive

### Tree fusing

- Maintain a pool of good trees
- For pairs of trees, find shared clades; swap compatible subtrees
- Accept if improvement found

### Adaptive parameters

- Track improvement rate (improvements per rearrangement evaluated)
- **When to ratchet**: When improvement rate drops below threshold
- **Ratchet intensity**: Start mild (resample 25% of characters); increase
  if consecutive iterations fail to find new basins

---

## Implementation Phases

### Phase 0: Data structures and scoring engine ✓
- [x] `DataSet`: Load from phyDat, bit-pack into blocks by state count
- [x] `TreeState`: Flat topology arrays + contiguous per-node state storage
- [x] Standard Fitch downpass (bit-packed, EW, no inapplicables)
- [x] R interface: `fitch_score(tree, dataset)` for testing
- [x] **Test**: scores match `preorder_morphy()` on datasets without
      inapplicable characters

### Phase 1: NNI search loop ✓
- [x] NNI rearrangement on `TreeState` (apply + undo)
- [x] Hill-climbing search loop (first-improvement)
- [x] R interface: `nni_search(starting_tree, dataset)` → edge list
- [x] **Test**: 20/20 random trees exact match with phangorn

### Phase 2: SPR with indirect calculation (partially complete)
- [x] SPR prune/regraft/undo on `TreeState`
- [x] Full two-pass after clip (not yet incremental Shortcut C)
- [x] Indirect tree length screening for rearrangement phase
- [x] SPR hill-climbing search with first-improvement
- [x] **Test**: 20 random DNA trees + multi-state morphological data match phangorn
- [ ] **FIX**: Exact indirect calculation (use union, not Fitch downpass —
      see "Lessons learned" §3)
- [ ] Incremental two-pass optimization (Shortcut C) for clipping phase:
  - [ ] Incremental downpass (rootward propagation with early stop)
  - [ ] Incremental uppass (propagate final state changes downward)
  - [x] Buffer/restore infrastructure for unclipping
  - [x] Per-character local cost storage
- [ ] Polymorphic terminal correction (uppass-as-internal-node)
- [ ] `local_reopt`: bit-packed quick lower bound for candidate regraft points

### Phase 2b: Exact indirect calculation ← **NEXT PRIORITY**
- [ ] Fix virtual root formula: `Y = final(A) ∪ final(D)` (plain union),
      NOT `fitch(final(A), final(D))` (intersection-then-union)
- [ ] Remove full rescore verification — indirect calc should be exact
- [ ] **Test**: indirect length = actual regrafted length on all test cases
      (zero tolerance; any discrepancy is a bug)
- [ ] **Test**: SPR search finds same optima without verification rescores

### Phase 2c: Group-by-weight (replaces weight expansion)
- [ ] Partition characters by weight, then by has_inapp, then by n_states
- [ ] Store per-block weight (uniform within block)
- [ ] Score: `total += block_weight * popcount(needs_union)`
- [ ] Update IW metadata (pattern_index mapping unchanged)
- [ ] **Test**: scores match weight-expanded version on all datasets

### Phase 3: TBR (separate agent)
- [ ] TBR bisect/reconnect/undo (extends SPR with subtree rerooting)
- [ ] Subtree rerooting via final states (no re-optimization needed)
- [ ] Full TBR search loop
- [ ] **Test**: matches R-side TBR search results
- [ ] (Optional) Union construct for bulk destination rejection:
  - [ ] Per-node union sets (incrementally maintained)
  - [ ] Top-down pruning of destination subtrees
  - [ ] Profile to confirm rearrangement phase is bottleneck before investing

### Phase 4: Inapplicable characters ← **PRIORITY after 2b/2c**
- [ ] Bit-packed NA-aware first downpass with mask-based case selection
- [ ] First uppass (applicability propagation)
- [ ] Second downpass (corrected scoring)
- [ ] NA-aware indirect calculation (virtual root with union, adapted for
      three-pass final states)
- [ ] Integrate into SPR search loop (inapplicable blocks scored with
      three-pass; standard blocks with one-pass)
- [ ] **Test**: exact match with `morphy_length()` on all
      `inapplicable.phyData` datasets
- [ ] **Test**: SPR search on inapplicable data finds same optima as
      existing R-side search

### Phase 5: Implied weights
- [ ] Per-character step extraction via scatter-add from bit-packed downpass
- [ ] IW scoring with early termination
- [ ] Character ordering heuristic (most variable first)
- [ ] IW-aware `local_reopt` (convert to fit as we accumulate)
- [ ] **Test**: IW scores match existing `morphy_iw()`

### Phase 6: Ratchet and heuristics
- [ ] Ratchet: character resampling + weight perturbation
- [ ] Cheap ratchet: standard Fitch for escape phase
- [ ] Sectorial search
- [ ] Tree pool with split-based deduplication
- [ ] Zero-length branch collapsing:
  - [ ] Shortest-path shortcut (approximate, configurable)
  - [ ] Partial reoptimization using final states from divided tree
  - [ ] Collapsed-tree comparison before pool insertion
- [ ] Progress reporting back to R
- [ ] `R_CheckUserInterrupt()` integration

### Phase 7: Adaptive search and polish
- [ ] Improvement rate tracking
- [ ] Adaptive ratchet intensity
- [ ] `MaximizeParsimony2()` R wrapper
- [ ] Constraint support
- [ ] Profile parsimony support
- [ ] Benchmark suite vs. existing `MaximizeParsimony()`
- [ ] Documentation

---

## File plan

| File | Purpose |
|------|---------|
| `src/ts_data.h` | `CharBlock`, `DataSet` — bit-packed character storage |
| `src/ts_data.cpp` | Loading from R, bit-packing, character compression |
| `src/ts_tree.h` | `TreeState` — topology + state sets, rearrangement ops |
| `src/ts_tree.cpp` | NNI/SPR/TBR apply/undo, partial rescore |
| `src/ts_fitch.h` | Fitch scoring: standard, inapplicable, local_reopt |
| `src/ts_fitch.cpp` | Bit-packed Fitch implementations |
| `src/ts_score.h` | Scoring interface: EW, IW, profile |
| `src/ts_score.cpp` | IW scoring with early termination |
| `src/ts_search.h` | `SearchEngine`, `SearchParams` |
| `src/ts_indirect.h` | Indirect calculation, incremental two-pass, union construct |
| `src/ts_indirect.cpp` | Incremental downpass/uppass, polymorphic correction |
| `src/ts_collapse.h` | Zero-length branch collapsing with partial reoptimization |
| `src/ts_search.cpp` | Main search loop, ratchet, sectorial, adaptive |
| `src/ts_splits.h` | Split-based topology hashing and deduplication |
| `src/ts_rcpp.cpp` | Rcpp exports: R ↔ C++ bridge |
| `R/MaximizeParsimony2.R` | R wrapper for new search engine |

The `ts_` prefix keeps our files clearly separated from existing morphy
source.

---

## Testing strategy

Each phase is tested against the existing morphy-based implementation:

- **Score equivalence**: For every tree scored, verify our score matches
  `preorder_morphy()` / `morphy_iw()`. Run on all 30 `inapplicable.phyData`
  datasets plus random trees.
- **Search equivalence**: On small datasets (6–10 tips), verify that our
  search finds exactly the same set of optimal trees.
- **Split hash correctness**: Verify deduplication against `as.TreeNumber()`
  for small trees and against explicit split comparison for larger trees.
- **Partial rescore correctness**: After every partial rescore, periodically
  do a full rescore and assert they match (debug mode only).
- **Fuzz testing**: Generate random trees, random rearrangements, verify
  scores throughout.

---

## Resolved questions

1. **Inapplicable mask algebra**: All three passes (first downpass, first
   uppass, second downpass) reduce to branchless mask operations. **Verified
   in R** against morphy's per-character implementation on all edge cases.
   ~4k + 10 operations per block of 64 characters, per node, per pass.

2. **Character weights**: Expand patterns by frequency so all characters
   have weight 1. Expansion is minimal on real datasets (max ~6 blocks
   of 64 for the largest dataset). Eliminates weight complications from
   the entire scoring and `local_reopt` pipeline.

3. **Partial rescoring with inapplicables**: The three-pass structure means
   changes can propagate globally via the uppass. Strategy: always full
   rescore for inapplicable blocks after a move passes `local_reopt`;
   partial rescore for standard Fitch blocks only. Morphy takes the same
   approach (its `mpl_fitch_NA_local_reopt` marks all characters for update).

4. **phyDat conversion path**: Straightforward. Contrast matrix gives state
   sets per token; index + weight give pattern expansion; levels identify
   the inapplicable state (`"-"`). Typical datasets: 7–65% of patterns have
   inapplicable data; the rest go into standard Fitch blocks.

5. **Indirect calculation algorithm** (Goloboff 1996): Use incremental
   two-pass optimization (Shortcut C) for the clipping phase of SPR/TBR.
   This is exact, fast (0.15–0.20 o/t/c), and provides the final state sets
   needed for both destination scoring and zero-length branch collapsing.
   Gladstein's (1996) incremental downpass-only method is 2.5–5× slower
   because it lacks final states and requires more work per destination.

6. **Zero-length branch collapsing**: Use final states from the divided tree
   to collapse with only 10–20% of characters needing reoptimization.
   Shortest-path shortcut for fast rejection of duplicate collapsed
   topologies. Important for search effectiveness (Goloboff 1996): prevents
   the tree pool from filling with trivially different dichotomous
   resolutions. Note asymmetric reachability as a known limitation.

## Open questions

1. **Wagner (ordered) characters**: Don't fit the Fitch intersection/union
   model. Defer until after Phase 5. When added, they'll need their own
   block type with a different downpass kernel.

2. **Profile parsimony**: Needs per-character step counts mapped through a
   lookup table. Same extraction mechanism as IW. Defer to Phase 7.

3. **Memory layout tuning**: Node-major storage (all words for one node
   contiguous) is natural for tree traversals. The inner loops are stride-1
   accesses that compilers should auto-vectorize. Profile before considering
   alternatives.

4. **Compiler portability of `popcount`**: C++20 has `std::popcount`.
   Older standards need `__builtin_popcountll` (GCC/Clang) or equivalent.
   R packages must compile on all CRAN platforms. Use a wrapper with
   fallback.

5. **Rcpp vs raw `.Call`**: Rcpp adds convenience but also compile-time
   overhead and some runtime overhead from SEXP wrapping. For the main
   search entry point (called once, runs for minutes), Rcpp is fine.
   For any per-tree-scored functions exposed for testing, raw `.Call`
   may be leaner. Decision: use Rcpp unless profiling shows a reason not to.
