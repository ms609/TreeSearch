# Plan: Heuristic Search Tricks for the C++ Tree Search Engine

This document specifies the heuristic search strategies ("tricks") to be
implemented in C++ on top of the core engine described in
`tree-search-engine-plan.md`. These tricks are what make the difference
between a basic TBR search and an effective search on datasets with 50+
taxa and significant character conflict.

**Key references:**
- Goloboff (1999) — original algorithms: SS, DFT, TF
- Nixon (1999) — parsimony ratchet
- Goloboff & Pol (2007) — SS vs Rec-I-DCM3; we follow SS
- Goloboff, "From Observations to Optimal Phylogenetic Trees", Ch. 5

**Current R-side status:** `MaximizeParsimony()` implements a 4-phase
search (initial TBR → ratchet → final TBR) with character-resampling
ratchet. Sectorial search code exists but is commented out and
non-functional. No tree fusing or tree drifting.

---

## Prerequisites from the core engine

These tricks assume the core engine (Phases 0–5 of
`tree-search-engine-plan.md`) provides:

| Requirement | Used by |
|---|---|
| `TreeState` with flat topology arrays + per-node state sets | All |
| Full Fitch downpass (EW/IW, ± inapplicables) | All |
| `local_reopt` for quick move evaluation | SS, TBR in ratchet/drift |
| SPR/TBR apply + undo on `TreeState` | All |
| Partial rescore after rearrangement | SS, TF |
| Per-character step counts (`char_steps[]`) | Ratchet (reweighting), IW |
| Split-based topology hashing | Tree pool dedup |

Additionally, the tricks need two facilities not yet in the core plan:

1. **Downpass state extraction at arbitrary nodes** — sectorial search
   needs to read off the first-pass (downpass) state sets at the nodes
   bordering a sector to construct the reduced dataset's HTU terminals.

2. **Subtree score accounting** — tree fusing needs the subtree length
   (sum of steps within a clade) stored per node, updated incrementally.

---

## 1. Sectorial search (SS)

### 1.1 Rationale

The key insight (Goloboff 1999; Goloboff & Pol 2007) is that a large tree
can be viewed as composed of many sub-problems. RAS + TBR may solve each
sub-problem individually with probability p, but the probability of solving
all of them simultaneously is p^k — exponentially small. Sectorial search
solves them one at a time.

Unlike Rec-I-DCM3 (which selects terminal taxa and changes the search
landscape), SS uses **first-pass HTU state sets** at sector boundaries.
This means any improvement of N steps in the sector is exactly N steps
improvement in the full tree. No global TBR is needed to "fix up" after
a sector improvement (though periodic global TBR can find further
improvements enabled by sector changes).

### 1.2 Sector selection

Three modes, following TNT:

#### 1.2.1 Random sectorial search (RSS)

Pick a random internal node; the sector is the clade rooted at that node.
Accept if sector size is in `[min_sector_size, max_sector_size]`; reject
and re-pick otherwise. For a tree of T taxa, do ~2T/S to 3T/S random
picks to probabilistically cover the whole tree.

Default size range: `[min(6, T/10), max(T/4, 50)]`.

#### 1.2.2 Exclusive sectorial search (XSS)

Partition the tree into P non-overlapping sectors of approximately equal
size. This guarantees full coverage. The partition is computed by a
greedy algorithm:

```
xss_partition(tree, P):
    targets = T / P  (target sector size)
    Walk tree bottom-up; when a subtree reaches ~targets nodes,
    mark it as a sector boundary.
    Adjust to minimize size variance.
```

The user specifies a range for P (e.g. 4–8); each round picks a random
P in that range, so successive rounds use different sector boundaries.

#### 1.2.3 Constrained sectorial search (CSS)

Uses a consensus tree (of the current best tree and the tree from the
previous replicate) as a constraint. Polytomies in the consensus with
degree above a threshold (default 10) identify poorly-resolved regions.
Sectors are created around these polytomies.

This is lower priority for initial implementation. RSS and XSS are the
workhorses.

### 1.3 Reduced dataset construction

Given a sector (a clade rooted at node `v`):

1. **Reroot** the full tree so that `v`'s parent becomes a leaf-like
   terminal in the reduced dataset. Specifically: compute the downpass
   state sets for the node that represents "the rest of the tree" by
   rerooting at the branch connecting the sector to the rest. The
   state sets at the sister/ancestor nodes bordering the sector become
   "HTU terminals" in the reduced dataset.

2. **Terminals** of the reduced dataset:
   - All tips within the sector (real OTUs)
   - HTU nodes at the sector boundary (1–3 nodes, carrying downpass
     state sets from the full tree)

3. **Character filtering**: Exclude characters that are uninformative
   for the reduced terminal set. This is a quick bitwise check per
   block.

4. **Topology**: The initial topology for the sector is the current
   subtree topology, already available in `TreeState`.

```cpp
struct ReducedDataset {
    DataSet data;             // bit-packed, subset of characters
    TreeState subtree;        // initial topology for the sector
    int sector_root;          // node index in full tree
    int n_real_tips;          // OTUs (not HTUs)
    int n_htus;               // boundary HTU count
    // Mapping back to full tree:
    std::vector<int> tip_to_full;   // reduced tip → full tree node
    std::vector<int> char_to_full;  // reduced char → full char index
};
```

### 1.4 Searching the reduced dataset

The search strategy for a sector depends on its size:

| Sector size | Strategy |
|---|---|
| ≤ 50 tips | 3 × RAS + TBR; if lengths differ, 3 more |
| 51–200 tips | RAS + TBR + tree drifting (6 cycles) |
| > 200 tips | RAS + TBR + SS (recursive) + TF |

These thresholds are configurable. The key point: the same search engine
is used recursively on the reduced dataset.

### 1.5 Reinsertion

If the sector search finds a tree shorter than the original resolution:

1. Map the improved sector topology back to full tree node indices.
2. Replace the subtree in `TreeState` (update `parent[]`, `left[]`,
   `right[]`).
3. Partial rescore: only nodes along the path from sector root to tree
   root need state-set updates.
4. Update `char_steps[]` if using IW.

If the score is equal, optionally accept the alternative resolution
(important for exploring flat landscapes in morphological datasets).

### 1.6 Global TBR between rounds

After each complete round of XSS (all sectors processed), or after every
~2T/S RSS picks, run a round of global TBR. Sector improvements may
enable further improvements elsewhere in the tree.

### 1.7 Performance expectations

From Goloboff Ch. 5: RSS on top of RAS + TBR finds near-optimal trees
~60× faster (at cost of 3.3× more time per replicate). XSS finds
near-optimal trees ~190× faster (at cost of 21× more time per
replicate). The net gain is enormous for datasets above ~100 taxa.

---

## 2. Parsimony ratchet

### 2.1 Algorithm

The ratchet (Nixon 1999, with TNT modifications from Goloboff 1999)
alternates perturbation and search phases:

```
ratchet(tree, n_cycles):
    best = tree
    for i in 1..n_cycles:
        // Perturbation phase
        perturbed_weights = perturb_characters(dataset)
        tree = tbr_search(tree, perturbed_weights,
                          accept_equal=true,
                          stop_after=T/8 changes)

        // Search phase
        tree = tbr_search(tree, original_weights)
        if score(tree) < score(best):
            best = tree

    return best
```

### 2.2 TNT-style modifications (important)

1. **Perturbation phase accepts equal-score rearrangements.** This is
   crucial — it allows drifting through equally-optimal trees under the
   perturbed weights, producing larger topological changes than only
   accepting improvements.

2. **Low reweighting probability.** Because equal-score moves are
   accepted, only ~4% of characters need reweighting (not 50% as in
   original Nixon). Higher reweighting + equal-acceptance drifts too far.

3. **Stopping rule for perturbation phase.** Stop after accepting T/8
   rearrangements (min 20, max 200), or after 99% of possible clippings
   have been tried. Do not run TBR to completion.

4. **Cheap ratchet option (from core plan).** During the perturbation
   phase, use standard Fitch (no inapplicable correction, no IW) for
   speed. The perturbed landscape only needs to push the tree somewhere
   different, not be scored accurately.

### 2.3 Character perturbation

Two modes:

- **Reweighting** (default): For each character, with probability p
  (default 0.04), multiply weight by a random factor (e.g. uniform in
  [0, 2]) or set to 0. With bit-packed scoring, this is done per-block
  by modifying `block_weight[]`.

- **Resampling** (Nixon original): Resample characters with
  replacement. This is equivalent to reweighting with Poisson-
  distributed weights. The existing R `Ratchet()` uses this approach
  via `tabulate(sample(...))`.

### 2.4 Integration with IW

When using implied weights, the ratchet can optionally use equal weights
during the perturbation phase (`ratchEW` option in current R code). This
is faster and the perturbation landscape doesn't need to match the real
scoring function.

---

## 3. Tree drifting

### 3.1 Algorithm

Like the ratchet, drifting alternates perturbation and search phases.
But instead of reweighting characters, the perturbation phase accepts
**suboptimal rearrangements** with probability that decreases with
suboptimality.

```
drift(tree, n_cycles):
    best = tree
    for i in 1..n_cycles:
        // Perturbation phase: accept suboptimal moves
        tree = tbr_drift(tree,
                         afd_limit,    // max absolute fit difference
                         rfd_limit,    // max relative fit difference
                         max_changes=T/8)

        // Search phase: standard TBR
        tree = tbr_search(tree, accept_equal=false)
        if score(tree) < score(best):
            best = tree

    return best
```

### 3.2 Acceptance criterion

For a rearrangement that increases score by Δ = Y - X:

1. If Δ = 0: always accept (allows drifting through equal-score space).
2. If Δ > `afd_limit`: always reject.
3. Otherwise, compute the **relative fit difference (RFD)**:
   - C = sum of weighted step decreases across characters
   - F = sum of weighted step increases across characters
   - RFD = (F - C) / F
4. If RFD > `rfd_limit`: reject.
5. Otherwise, accept with probability that decreases with Δ and
   increases with how close RFD is to 0. Include a drift-distance
   penalty factor `f` that increases as the tree drifts further from
   the starting score.

### 3.3 RFD computation

RFD requires per-character score changes, not just total score. During
TBR evaluation, for each candidate rearrangement we already compute
per-block step counts. For blocks where steps changed:

```cpp
// During local_reopt or full rescore:
for each block b:
    old_steps_b = ...; new_steps_b = ...;
    delta_b = new_steps_b - old_steps_b;
    if (delta_b > 0) F += delta_b * block_weight[b];
    if (delta_b < 0) C += (-delta_b) * block_weight[b];
```

This is approximate (block-level, not character-level) but sufficient for
the acceptance decision and much cheaper than character-level accounting.

### 3.4 Efficiency shortcut

Before computing RFD, check: if the score increase Y from reinsertion
satisfies Y ≥ X / (1 - rfd_limit), then RFD must exceed the limit
regardless of C. This rejects most suboptimal moves cheaply, based
solely on the reinsertion cost.

### 3.5 Alternation with equal-score drifting

TNT alternates drift cycles (accepting suboptimal moves) with
equal-score drift cycles (accepting only Δ = 0 moves). The latter
is far more efficient than saving and swapping multiple trees: a few
accepted equal-score moves can reach trees 20+ TBR moves away from
the start almost instantly.

### 3.6 Advantages over ratchet

- No need to modify character weights (simpler internal state).
- Directly addresses character conflict via RFD.
- Equal-score drifting is very efficient for flat landscapes
  (common in morphological datasets).

Disadvantages:
- `afd_limit` needs tuning for datasets with many characters.
- For phylogenomic datasets, even one step is a large absolute
  difference, requiring large `afd_limit` (but keep `rfd_limit` low,
  ~0.1–0.2).

---

## 4. Tree fusing (TF)

### 4.1 Algorithm

Tree fusing combines the best parts of multiple suboptimal trees into a
single better tree.

```
tree_fuse(tree_pool):
    recipient = best tree in pool
    improved = true
    while improved:
        improved = false
        for each donor in pool (donor ≠ recipient):
            consensus_groups = strict_consensus(recipient, donor)
            // Try exchanges in bottom-up order:
            for each shared group g (smallest first):
                score_delta = exchange_score(recipient, donor, g)
                if score_delta < 0:  // improvement
                    apply_exchange(recipient, donor, g)
                    improved = true
                    skip all ancestor groups of g
            // Optionally accept equal-score exchanges too

        if improved:
            recipient = tbr_search(recipient)  // global TBR

    return recipient
```

### 4.2 Shared group identification

Two trees share a group (clade) when they both have a bipartition
splitting the same set of taxa. Finding shared groups:

1. Compute splits for both trees.
2. Intersect the split sets.
3. Shared splits define the nodes eligible for exchange.

With the split-based hashing infrastructure from the core plan, this is
efficient. For the consensus guide tree, only the shared splits matter.

### 4.3 Fast exchange scoring

For a candidate exchange of group `g` between donor and recipient:

1. The donor's subtree for `g` has known downpass state sets and subtree
   length (stored from its last full score).
2. Remove the recipient's subtree for `g`.
3. Insert the donor's subtree. Score the change by walking from the
   exchange point to the root, using the pre-computed downpass states
   of the donor subtree and the recipient's sister-group states.

This is O(depth) per exchange, not O(T×C). With downpass states cached,
most exchanges are evaluated in microseconds.

### 4.4 Bottom-up ordering

Try smaller shared groups first. If exchanging a small group improves
the tree, skip all ancestor groups that contain it (they may no longer
be comparable after the exchange). This maximizes the number of
independent improvements found per fusion round.

### 4.5 Polytomy skipping

Skip shared groups that form polytomies of degree ≤ 3 in the consensus.
These are trivially resolved by TBR and unlikely to represent genuine
local optima.

### 4.6 Equal-score exchanges

For morphological datasets with flat landscapes, accepting equal-score
exchanges is valuable. It moves the recipient tree to a different region
of the equally-optimal plateau, where subsequent TBR may find
improvements.

### 4.7 Integration with search

TF is most powerful when applied to trees from independent replicates
(different RAS starting points). Each replicate may have solved different
parts of the tree optimally. TF assembles the globally optimal
combination.

**Critical lesson from Goloboff & Pol (2007):** Roshan et al.'s TNT
runs discarded previous results with each new `xmult` call. The simple
fix — fusing results across replicates — dramatically improved
performance. Our search loop must always pool trees for TF.

---

## 5. Tree pool management

### 5.1 Pool structure

```cpp
struct TreePool {
    std::vector<TreeState> trees;
    std::vector<double> scores;
    std::vector<uint64_t> split_hashes;  // for dedup
    int max_size;          // e.g. 100
    double suboptimal;     // tolerance for keeping suboptimal trees
};
```

### 5.2 Deduplication

Before adding a tree to the pool, check its split hash against existing
entries. On hash collision, compare full split vectors. This prevents
the pool from filling with duplicates, which wastes TF effort.

For small trees (≤ 51 tips), Tromp ID from TreeTools is an exact fast
path.

### 5.3 Pool policy

- Keep all trees within `suboptimal` steps of the best.
- When pool is full, evict the worst tree.
- After TF produces an improvement, evict all trees worse than the new
  best minus `suboptimal`.

---

## 6. Combined search strategy ("driven search")

### 6.1 Main loop

Following the TNT `xmult` / `combosearch` model:

```
driven_search(dataset, params):
    pool = empty TreePool

    for rep in 1..max_replicates:
        // Build starting tree
        tree = wagner_addition(dataset, random_sequence)
        tree = tbr_search(tree)

        // Sectorial search
        tree = xss_search(tree, dataset,
                          min_partitions=params.min_P,
                          max_partitions=params.max_P)

        // Ratchet or drifting (user choice or both)
        tree = ratchet(tree, n_cycles=params.ratch_cycles)
        // or: tree = drift(tree, n_cycles=params.drift_cycles)

        pool.add(tree)

        // Periodic tree fusing
        if rep % params.fuse_interval == 0 or pool.size >= params.fuse_threshold:
            best = tree_fuse(pool)
            pool.add(best)

        // Convergence check
        if pool.hits_to_best >= params.target_hits:
            break

        // Consensus stability check (optional)
        if params.check_consensus and rep > params.min_reps:
            if consensus_stable(pool):
                break

        R_CheckUserInterrupt()

    return pool.best_trees()
```

### 6.2 The XSS + TF combo for very large trees

For datasets with 500+ taxa, the combosearch structure from Goloboff &
Pol (2007) is most effective:

```
combosearch(tree, dataset, params):
    for round in 1..:
        steps_saved = 0
        for xss_round in 1..3:
            P = random_int(params.min_P, params.max_P)
            sectors = xss_partition(tree, P)
            for sector in sectors:
                // Analyze sector with xmult-like search
                improved = search_sector(tree, sector, dataset)
                steps_saved += improved
            tree = tbr_search(tree)  // global TBR

        if steps_saved < params.min_step_decrease:
            consecutive_failures++
            if consecutive_failures >= 3:
                // Pool tree and restart or fuse
                pool.add(tree)
                tree = tree_fuse(pool)
                consecutive_failures = 0
        else:
            consecutive_failures = 0
```

### 6.3 Minimum step decrease (MSD)

The MSD parameter controls when to give up on the current XSS cycle and
try something different (fusing, new starting point). When three
consecutive XSS rounds save fewer than MSD steps, the strategy switches.

Recommended defaults (from Goloboff & Pol 2007):
- DS-11K (11,361 taxa): MSD = 10
- Typical morphological (50–500 taxa): MSD = 1–5

### 6.4 Adaptive parameter selection

Track improvement rate (improvements per rearrangement evaluated):

- **High rate** (early search): use TBR, large sectors, fewer ratchet
  cycles.
- **Low rate** (approaching optimum): use SS with smaller sectors,
  more ratchet/drift cycles, more frequent TF.
- **Near zero** (at or very near optimum): emphasize TF with
  equal-score exchanges, consensus stability checking.

---

## 7. Flat landscape strategies for morphological data

Morphological datasets often have many MPTs and flat landscapes where
equally-optimal trees span large regions of tree space. This creates
specific challenges (Goloboff Ch. 5, §5.7.3):

### 7.1 Equal-score acceptance everywhere

- Ratchet: already accepts equal-score moves in perturbation phase.
- Drifting: alternates suboptimal-acceptance with equal-score-only
  phases.
- SS: accept alternative sector resolutions of equal score.
- TF: accept equal-score exchanges between trees.

### 7.2 Increased perturbation length, decreased intensity

In flat landscapes, the perturbation phase should accept more changes
(larger `max_changes`) but with lower character-reweight probability
(ratchet) or lower `afd_limit` (drifting). This drifts further through
equal-score space without jumping to distant score basins.

### 7.3 Zero-length branch collapsing

During search, collapse zero-length branches (using dual representation:
store both binary and collapsed forms). When saving trees to the pool,
compare collapsed forms. This prevents the pool from filling with trees
that differ only in unsupported resolution, and speeds up TBR by
allowing earlier abandonment of rearrangement scoring (full buffer
effect).

### 7.4 Consensus-driven termination

For morphological data, the practical goal is often a stable strict
consensus rather than an exact set of MPTs. Drive the search toward
consensus stability: stop adding replicates when additional independent
hits to minimum length no longer reduce consensus resolution.

---

## 8. Implementation phases

These build on top of the core engine phases (0–5).

### Phase 6a: Ratchet
- [ ] Character reweighting (per-block weight modification)
- [ ] TNT-style perturbation phase (equal-score acceptance, T/8 stop)
- [ ] Cheap ratchet mode (standard Fitch during perturbation)
- [ ] Integration with existing ratchet R parameters
- [ ] **Test**: matches or beats R-side `Ratchet()` scores

### Phase 6b: Tree drifting
- [ ] Suboptimal-acceptance TBR with AFD/RFD criteria
- [ ] Block-level RFD computation
- [ ] Alternating drift/equal-score phases
- [ ] Drift-distance penalty factor
- [ ] **Test**: escapes local optima that ratchet alone cannot

### Phase 6c: Tree fusing
- [ ] Shared-group identification via split intersection
- [ ] Fast exchange scoring (root-ward walk with cached states)
- [ ] Bottom-up exchange ordering with ancestor skipping
- [ ] Equal-score exchange option
- [ ] Tree pool with split-based deduplication
- [ ] **Test**: fusing 5 independent RAS+TBR trees beats any single one

### Phase 6d: Sectorial search
- [ ] Reduced dataset construction with HTU state sets
- [ ] RSS: random sector selection within size bounds
- [ ] XSS: even partitioning algorithm
- [ ] Recursive search on reduced datasets
- [ ] Sector reinsertion with partial rescore
- [ ] Global TBR between rounds
- [ ] **Test**: SS(xmult) outperforms plain xmult on 100+ taxon datasets

### Phase 7: Combined search and R interface
- [ ] Driven search loop with convergence criteria
- [ ] Combosearch structure for very large trees
- [ ] MSD-based adaptive stopping
- [ ] Adaptive parameter selection
- [ ] Zero-length branch collapsing (dual representation)
- [ ] Consensus stability checking
- [ ] `MaximizeParsimony2()` R wrapper with all options
- [ ] Progress reporting via `R_CheckUserInterrupt()` + callbacks
- [ ] **Test**: benchmark suite vs existing `MaximizeParsimony()`

### Phase ordering rationale

6a (ratchet) first because it's simplest and the R code already has the
logic. 6b (drifting) next as it shares infrastructure with the ratchet.
6c (TF) before 6d (SS) because TF is simpler and dramatically improves
multi-replicate searches even without SS. SS is most complex (reduced
dataset construction, recursive search) and benefits from having TF
available for sector analysis.

---

## 9. File plan (additions to core engine)

| File | Purpose |
|---|---|
| `src/ts_ratchet.h/.cpp` | Ratchet: character perturbation + modified TBR |
| `src/ts_drift.h/.cpp` | Tree drifting: suboptimal acceptance + RFD |
| `src/ts_fuse.h/.cpp` | Tree fusing: shared groups, exchange scoring |
| `src/ts_sector.h/.cpp` | Sectorial search: sector selection, reduced datasets |
| `src/ts_pool.h/.cpp` | Tree pool: storage, dedup, eviction policy |
| `src/ts_driven.h/.cpp` | Driven/combined search: main loop, adaptive params |

---

## 10. Design decisions and non-goals

### Decisions

- **No Rec-I-DCM3.** Goloboff & Pol (2007) convincingly showed it's
  inferior to SS for the same search effort, and after a few iterations
  behaves as a perturbation method rather than true divide-and-conquer.
  The random dichotomization of the merged supertree destroys
  information, requiring expensive global TBR to recover.

- **SS with HTU state sets, not terminal selection.** This is the core
  advantage: sector improvements transfer exactly to the full tree.

- **Block-level RFD for drifting.** Character-level RFD would require
  unpacking bit-packed blocks. Block-level (64 characters at once) is
  a reasonable approximation and much faster.

- **Ratchet uses low reweight probability (4%) with equal-score
  acceptance.** Following TNT, not Nixon's original 50%.

### Non-goals (for now)

- **Tree hybridization** (Goloboff 2015): useful for very poorly
  structured data but complex. Defer.
- **Tree rebuilding** (Goloboff 2015): prune + reinsert random taxa.
  Simpler but lower priority than SS/TF.
- **Piñón fijo** (Goloboff 2015): Rec-I-DCM3-inspired perturbation.
  Low priority given SS superiority.
- **Parallel search**: multiple threads running independent replicates.
  The architecture supports this (independent `TreeState` per thread)
  but thread management is deferred.
