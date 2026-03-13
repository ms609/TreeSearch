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

### Phase 7: Starting trees and combined search
- [ ] C++ Wagner tree with indirect length calculation
- [ ] Random addition sequence support
- [ ] **Test**: C++ Wagner trees match R `AdditionTree()` scores
- [ ] Driven search loop with convergence criteria
- [ ] Combosearch structure for very large trees
- [ ] MSD-based adaptive stopping
- [ ] Adaptive parameter selection
- [ ] Zero-length branch collapsing (dual representation)
- [ ] Consensus stability checking
- [ ] `MaximizeParsimony2()` R wrapper with all options
- [ ] Progress reporting via `R_CheckUserInterrupt()` + callbacks
- [ ] **Test**: benchmark suite vs existing `MaximizeParsimony()`
- [ ] **Test**: Wagner vs NJ vs random starting trees comparison

### Phase 8: Parallelization
- [ ] Thread-safe tree pool (`ThreadSafePool`)
- [ ] Replicate-level parallelism with `std::thread` or `RcppThread`
- [ ] Atomic progress counter + main-thread interrupt checking
- [ ] `n_threads` parameter in R interface
- [ ] **Test**: verify deterministic results with fixed seeds per thread
- [ ] **Test**: scaling benchmarks (1, 2, 4, 8 threads)
- [ ] (Optional) Sector-level parallelism for XSS
- [ ] (Optional) Async TF background thread

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
| `src/ts_wagner.h/.cpp` | Wagner tree (greedy addition) with indirect scoring |
| `src/ts_parallel.h/.cpp` | Thread pool, thread-safe pool, parallel replicate runner |

---

## 10. Starting tree generation

### 10.1 Current implementation

`AdditionTree()` in R is already a Wagner tree (greedy addition): it
adds each taxon to the most parsimonious available position, with
randomized addition sequence by default. It supports IW scoring and
topology constraints. This is exactly what Goloboff Ch. 5 §5.3.1
describes as the standard approach for generating starting points.

### 10.2 C++ reimplementation

The C++ engine needs its own Wagner tree builder because the R version
uses `AddTipEverywhere()` + `TreeLength()` per taxon — O(T) candidate
trees scored from scratch at each step. A C++ version can use indirect
length calculation: at each step, evaluate all 2i−3 insertion positions
using the downpass state sets of the growing tree, scoring each insertion
in O(C) via `local_reopt`-style comparison rather than a full O(TC)
downpass. This reduces per-taxon cost from O(T×C) to O(C) per position
× O(T) positions = O(TC) total, but with much smaller constants.

```cpp
TreeState wagner_tree(const DataSet& data, std::vector<int> sequence) {
    // Start with 3-taxon pectinate tree
    TreeState tree = three_taxon_tree(data, sequence[0..2]);
    full_score(tree, data);

    for (int i = 3; i < n_tips; ++i) {
        int best_edge = -1;
        int best_delta = INT_MAX;
        // Try all insertion positions
        for (int e = 0; e < tree.n_edges(); ++e) {
            int delta = insertion_cost(tree, data, sequence[i], e);
            if (delta < best_delta) {
                best_delta = delta;
                best_edge = e;
            }
        }
        insert_tip(tree, sequence[i], best_edge);
        update_downpass(tree, data, best_edge);  // incremental
    }
    return tree;
}
```

### 10.3 Alternative starting trees to test

| Method | Description | When useful |
|---|---|---|
| **Wagner (RAS)** | Current approach; different random addition sequences | Default; well-understood |
| **NJ starting tree** | `TreeTools::NJTree()` already available | Fast; good for molecular-like data with branch-length signal |
| **Random tree** | Uniform random topology | Baseline; much further from optimal, wastes TBR time |
| **Informative addition** | Goloboff (2014a): add taxa in order of informativeness | May help in specific datasets; not widely tested |

Goloboff Ch. 5 §5.4.1 notes that random trees are generally worse
starting points than Wagner trees (further from optimal → more TBR time
wasted on initial descent). NJ trees could be useful as one starting
point in a pool but should not replace RAS Wagner trees as the primary
strategy.

**Testing plan:** Compare, for a suite of empirical datasets:
1. Time to reach best-known score from Wagner (RAS) vs NJ vs random
   starting trees, each followed by TBR.
2. Whether mixing starting tree types (e.g. 8 Wagner + 2 NJ) in the
   pool for TF improves convergence.

---

## 11. Parallelization

### 11.1 Where parallelism fits

The driven search (§6) is naturally parallel at the **replicate level**.
Each replicate (RAS + TBR + SS + ratchet/drift) is independent until
the TF step. This is the simplest and most effective parallelization
strategy.

```
parallel_driven_search(dataset, params, n_threads):
    pool = shared TreePool (thread-safe)

    parallel for rep in 1..max_replicates (n_threads workers):
        tree = wagner_tree(dataset, random_sequence)
        tree = tbr_search(tree)
        tree = xss_search(tree, dataset)
        tree = ratchet(tree)
        pool.add(tree)  // lock-free or mutex-protected

    // TF is sequential (reads from pool, modifies recipient)
    best = tree_fuse(pool)
    return best
```

### 11.2 Levels of parallelism

| Level | What's parallel | Complexity | Payoff |
|---|---|---|---|
| **Replicate-level** | Independent RAS+TBR+SS runs | Low (shared-nothing except pool) | High: linear speedup up to ~8–16 threads |
| **Sector-level** | XSS sectors analyzed concurrently | Medium (sectors share full-tree state) | Moderate: useful for very large trees where sectors are expensive |
| **Move-level** | TBR rearrangements evaluated in parallel | High (fine-grained, synchronization overhead) | Low: individual moves are cheap; overhead may dominate |
| **Block-level (SIMD)** | 64 characters via popcount, already implicit | Zero (compiler auto-vectorizes) | Already captured by bit-packing design |

**Replicate-level parallelism is the clear first target.** It requires
almost no synchronization: each thread owns its own `TreeState` and
`DataSet` copy (the dataset is read-only and can be shared). The only
shared mutable state is the tree pool, which needs a mutex for
`add()`/`evict()` but these are infrequent operations.

### 11.3 Implementation plan

#### Phase 1: Thread-safe pool + replicate parallelism

```cpp
struct ThreadSafePool {
    std::mutex mtx;
    TreePool pool;

    void add(TreeState&& tree, double score) {
        std::lock_guard<std::mutex> lock(mtx);
        pool.add(std::move(tree), score);
    }

    double best_score() const {
        std::lock_guard<std::mutex> lock(mtx);
        return pool.best_score();
    }
};
```

Each worker thread:
1. Copies `DataSet` (or shares read-only reference).
2. Allocates its own `TreeState` + temp buffers.
3. Runs the full replicate pipeline.
4. Pushes result to the shared pool.

TF runs periodically on a single thread (or a designated "fuser"
thread), pulling trees from the pool.

Use `std::thread` or a simple thread pool. R integration via
`RcppParallel` or `RcppThread` (the latter is simpler and avoids
TBB dependency). **Must not call any R API from worker threads** —
all R interaction happens on the main thread.

#### Phase 2: Sector-level parallelism (optional)

Within a single XSS round, non-overlapping sectors are independent.
They can be farmed out to a thread pool. Each sector analysis needs
its own `ReducedDataset` + `TreeState`, which are independent once
constructed. The only shared state is the full tree, which is read-only
during sector analysis (sector results are applied sequentially after
all sectors complete).

This is useful for very large trees (500+ taxa) where individual sector
analysis takes seconds to minutes. For typical morphological datasets
(50–200 taxa), replicate-level parallelism is sufficient.

#### Phase 3: Async TF (optional)

A background "fuser" thread continuously:
1. Waits for new trees in the pool.
2. Runs TF on the current pool.
3. If TF finds an improvement, broadcasts the new best score to worker
   threads so they can update their `stopAtScore` targets.

This overlaps TF with ongoing replicates, reducing idle time.

### 11.4 R integration considerations

- `R_CheckUserInterrupt()` must be called from the main thread only.
  Worker threads signal a shared atomic flag; main thread checks it
  periodically.
- Progress reporting: workers increment an atomic counter; main thread
  reads it for progress callbacks to R.
- `n_threads` parameter exposed in R, defaulting to
  `parallel::detectCores() - 1` or 1 for sequential mode.
- All morphy objects are per-thread (morphy is not thread-safe).
  In the new C++ engine, the bit-packed `DataSet` is read-only and
  safely shared.

### 11.5 Expected speedup

For replicate-level parallelism with k threads:
- Near-linear speedup (0.85–0.95 × k) for the search phase.
- TF is sequential but fast (typically <1% of total time).
- Pool synchronization overhead is negligible.
- Memory scales linearly: each thread needs ~100 KB for a 100-taxon
  dataset's `TreeState` + temp buffers.

On a typical 8-core desktop, expect 6–7× speedup for the driven search.
This compounds with the algorithmic improvements from SS/TF/ratchet.

---

## 12. Design decisions and non-goals

### Decisions

- **Wagner tree in C++.** Reimplement `AdditionTree()` with indirect
  length calculation for speed. Keep the R version as correctness
  reference.

- **Replicate-level parallelism first.** Simplest, highest payoff,
  minimal synchronization. Sector-level and async TF are optional
  extensions for very large datasets.

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

- **Tree hybridization** (Goloboff 2015): For poorly structured datasets
  where tree fusing fails (near-optimal trees share few groups). Takes
  two trees, identifies two evenly-sized partitions with similar taxon
  composition, exchanges halves (A1A2 + B1B2 → A1B2 + A2B1), reinserts
  missing taxa by random addition + TBR. More effective than standard GA
  crossover because even-sized partition exchange preserves relative
  distances. Combined with cyclic taxon pruning/reinsertion in TNT
  scripts, outperforms Sampars/Hydra on random datasets by 16–30×
  (Goloboff 2015). Defer: mainly useful for random/unstructured data,
  not our primary use case.
- **Piñón fijo (pfijo)** (Goloboff 2015): Cyclic perturbation method
  borrowing from Rec-I-DCM3. Creates reduced datasets from tree sectors,
  but instead of true HTU downpass states (which give exact evaluations),
  uses states of one descendant of the HTU with probability depending on
  constant K. Lower K = stronger perturbation (states further from root).
  After searching reduced datasets (4 rounds RAS+TBR + tree-fusing with
  original resolution), does full TBR + drifting. Repeat. Best for
  random/poorly structured datasets where true SS doesn't help. Defer:
  lower priority than SS/TF/ratchet/drifting for structured morphological
  data.
- **Tree rebuilding** (Goloboff 2015): prune + reinsert random taxa.
  Simpler but lower priority than SS/TF. Can be combined with
  hybridization for a powerful cyclic perturbation on unstructured data.
- **Parallel search**: multiple threads running independent replicates.
  The architecture supports this (independent `TreeState` per thread)
  but thread management is deferred.

## Future experiments

Ideas worth testing once the core tricks are working. These are not
commitments — they're hypotheses to evaluate empirically.

### Split-novelty bonus during drifting

**Idea:** During tree drifting phases, bias move acceptance toward
rearrangements that introduce splits not yet seen in the tree pool.

**Motivation:** Barker et al. 2025 (larch) found that scoring SPR moves
by `parsimony_change − number_of_novel_splits` outperformed pure
parsimony scoring for exploring the space of optimal trees. Their
implementation requires a full DAG of all optimal histories, but a
lightweight version may be feasible using our split-hash infrastructure.

**Why it might be cheap for us:** The tree pool already maintains a set
of split hashes for deduplication. During drifting, when we're already
accepting some suboptimal moves, we could check whether a proposed
rearrangement introduces splits absent from the pool's hash set. If so,
give the move a small acceptance bonus (e.g., treat it as if it were
`bonus` steps better than it actually is). This requires only hash
lookups against the existing pool structure — no DAG, no RF distance
computation.

**Concrete test:** Compare driven search (ratchet + drift + TF + SS)
with and without the novelty bonus on a handful of morphological
datasets (varying size and homoplasy). Measure: (a) number of distinct
optimal topologies found, (b) best score found, (c) time to best score.

**Risk:** On morphological data with implied weights or inapplicable
characters, the parsimony landscape may be structured differently from
the molecular DNA datasets where larch was tested. The bonus could
waste time exploring genuinely bad regions rather than crossing
meaningful valleys. The experiment should include a control where
drifting uses the standard RFD-based acceptance criterion.

**When to try:** After Phase 6 (combined driven search) is working and
benchmarked. The drifting infrastructure and pool split-hash tracking
must both be in place.

---

## Papers reviewed

Papers consulted during plan development, with notes on influence:

| Paper | Influence |
|---|---|
| Goloboff 1996, "Methods for faster parsimony analysis" | **Major** — indirect calculation, Shortcut C, union construct, zero-length collapsing. All added to core engine plan. |
| Goloboff & Pol 2007, "On divide-and-conquer strategies…" | **Major** — established SS with HTU state sets as superior to Rec-I-DCM3. Shaped sectorial search design. |
| Goloboff 2015, Chapter 5 (textbook) | **Major** — comprehensive reference for ratchet, drifting, tree fusing, sectorial search, combined strategies. |
| Nixon 1999, parsimony ratchet | **Moderate** — original ratchet paper; TNT's refinements (4% reweight, equal-score acceptance) supersede the original parameters. |
| Moilanen 1999, "Searching for Most Parsimonious Trees with Simulated Evolutionary Optimization" | **None** — GA + local search hybrid (PARSIGAL). Conceptually superseded by TNT's targeted strategies. Bit-packing idea already standard. |
| Andreatta & Ribeiro 2002, "Heuristics for the Phylogeny Problem" | **None** — Wagner construction variants, NNI/SPR local search, VND, GRASP, VNS. All standard or superseded by TNT-class strategies (ratchet, SS, TF). VND (escalate NNI→SPR on stagnation) is a weaker version of the intuition behind sectorial search. Small instances, not competitive. |
| Komusiewicz et al. 2023, "On the Complexity of Parameterized Local Search for Maximum Parsimony" | **None (theoretical validation)** — Shows searching k-NNI/SPR/TBR/ECR neighborhoods is W[1]-hard (brute-force essentially optimal). But k-sECR (contracted edges form a connected subtree) is FPT: k^O(k)·\|I\|^O(1). Paper explicitly identifies sECR = sectorial search (Goloboff 1999). Confirms our focus on SS is theoretically well-founded, but offers no new implementation strategies. |
| Barker et al. 2025, "larch: mapping the parsimony-optimal landscape of trees for directed exploration" (bioRxiv) | **None (interesting future idea)** — MADAG (mutation-annotated DAG) compactly stores large families of equally parsimonious trees for molecular/viral data. Key innovation: novelty-directed search — accept SPR moves introducing new splits even at slight parsimony cost; preferentially optimize topologically distant trees from the ensemble. Claims to outperform TNT via matOptimize, but tested only on molecular DNA data (standard Fitch), not morphological data with inapplicables/implied weights. MADAG structure exploits densely-sampled molecular regime (1–2 mutations per edge), not applicable to morphological data. Novelty-directed valley-crossing is conceptually interesting as a future refinement of our flat landscape / drifting strategies, but ratchet + drifting + TF already cover local optima escape. Could revisit if split-hash pool tracking makes novelty scoring cheap. |
| Charleston 2001, "Hitch-Hiking: A Parallel Heuristic Search Strategy" | **None** — Parallel search where "driver" solutions share accepted perturbations with "hitcher" solutions. NNI re-coded as "Short Path Shift" (3-leaf identifier) for cross-tree transferability. Tested on 10-taxon simulated data with Great Deluge; marginal improvement over independent runs. A less targeted predecessor of the idea behind tree fusing (cross-pollination between parallel searches). Our replicate-level parallelism + TF already covers this more effectively. |
| Goloboff 2015, "Computer science and parsimony: a reappraisal" | **Moderate** — Demonstrates TNT outperforms all CS parsimony methods (Hydra, Sampars, GA+PR+LS) by 1–3 orders of magnitude. Key technical insight: with indirect calculation + union construct, TBR completes full swap cycles faster than SPR in absolute wall-clock time, making "variable neighbourhood" strategies (NNI→SPR escalation) counterproductive. Led to removing SPR/NNI as search strategies from core engine plan (TBR only). Also provides concrete details on pfijo (degraded sectorial search as cyclic perturbation) and tree hybridization (even-partition exchange between trees) for poorly structured datasets — enriched non-goals entries. |
| Gladstein 1997, "Efficient Incremental Character Optimization" | **None** — Formalizes incremental downpass evaluation using attribute evaluation theory (Hudson 1991). After topology change, walk from modified node to root, recomputing preliminary state sets per character, stopping when state set is unchanged. Averages < 3 node visits per character per tree. This is the "incremental downpass only" approach already noted in the core engine plan as inferior to Goloboff 1996's Shortcut C (incremental two-pass: 0.15–0.20 o/t/c, 2.5–5× faster). Clean formalization but no new algorithmic content. |
| Vazquez-Ortiz 2016 thesis + Vazquez-Ortiz & Rodriguez-Tello 2011 (CIB) + Richer et al. 2012 (SAMPARS) | **None** — SA and ILS metaheuristics with NNI/SPR neighborhood combinations (switching by temperature, iteration count, or probability). Also: bottom-up path-relinking between trees (weaker than tree fusing — no exact score evaluation of exchanges), parsimony score prediction via regression (academic, not relevant to search), GPU Fitch evaluation (full downpass parallelized per character — irrelevant with indirect calculation). These are the specific CS methods Goloboff 2015 showed TNT outperforms by 30–40× on the same benchmark datasets. |
| Ribeiro & Vianna 2005, "A GRASP/VND heuristic for the phylogeny problem" | **None** — GRASP + VND with SPR and 2-SPR (two successive SPR moves). Uses indirect calculation for reconnection but no TBR, no union construct. VND escalation (SPR→2-SPR) is a weaker version of what TBR achieves natively. One of the CS papers Goloboff 2015 showed to be orders of magnitude slower than TNT. |
| Richer 2008, "Three new techniques to improve phylogenetic reconstruction with MP" (Hydra tech report) | **None** — Three techniques: (1) Progressive Neighborhood — start SPR/TBR, shrink to NNI; counterproductive with efficient TBR. (2) DiBIP crossover — average topological distance matrices of two parents, rebuild via UPGMA; less targeted than tree fusing or hybridization. (3) SSE2 vectorization of full Fitch downpass; our bit-packing (64 chars per uint64_t) achieves the same effect, and with indirect calculation the bottleneck moves to the 3-node comparison anyway. Hydra explicitly lacks indirect calculation, explaining its speed disadvantage vs TNT despite SIMD optimization. |
| Goëffon, Richer & Hao 2005, "Local Search for the Maximum Parsimony Problem" (ICNC proceedings) | **None** — Short (6pp) conference paper introducing an array-based tree representation and a new subtree swapping neighborhood. Tested on random benchmarks and 8 real sequences. Superseded by the authors' own 2008 journal paper (below). |
| Goëffon, Richer & Hao 2008, "Progressive Tree Neighborhood Applied to the Maximum Parsimony Problem" (IEEE/ACM TCBB) | **None (interesting observation, not actionable)** — Full journal version of the 2005 conference paper. Key idea: Parametric Progressive Neighborhood (PPN) — start with full SPR, linearly shrink regraft distance *d* to NNI over *M* iterations. Empirical justification (Section 6): best-improvement descent shows improving neighbors get closer as search progresses, so large moves are only useful early on. PPN beats fixed NNI and fixed SPR in simple first-improvement descent; finds Zilla best-known score (16,218). However: (1) only tested with simple descent, not within ratchet/SS/drifting framework; (2) Goloboff 2015 showed TBR with indirect calculation is faster in wall-clock time than SPR, making SPR→NNI shrinkage counterproductive when TBR is available; (3) the Richer 2008 Hydra tech report already lists PPN as technique (1) and it was part of the system TNT outperforms by 30–40×. The underlying observation (improving moves get more local as search converges) is real but already exploited by our phased strategy (broad TBR → ratchet perturbation → final TBR refinement). |
| Blazewicz et al. 2011, "Adaptive memory programming: local search parallel algorithms for phylogenetic tree construction" (Ann. Oper. Res.) | **None** — Parallel tabu search / adaptive memory programming applied to MP with standard NNI/SPR/TBR neighborhoods. Claims superlinear speedup from parallelization. Tested on simulated data and hepatitis C sequences. Only 1 citation; no novel rearrangement operators or scoring tricks. Cites both Goëffon et al. 2005 and 2008. Low relevance: parallelization wrapper over standard neighborhoods, no ideas beyond what ratchet + drifting + TF already cover. |
| Hoang et al. 2018, "MPBoot: fast phylogenetic maximum parsimony tree inference and bootstrap approximation" (BMC Evol. Biol.) | **None** — Adapts IQ-TREE's UFBoot approach to parsimony. Search uses SPR (radius 3 or 6) + ratchet (50% site duplication) + candidate set of 5 trees (GA-inspired). Main contribution is fast bootstrap via REPS (resampling parsimony score — reweight site-pattern scores instead of re-searching each bootstrap MSA). Finds scores better than fast-TNT but comparable to intensive-TNT. Search algorithm is a strict subset of what TreeSearch plans: SPR not TBR, no sectorial search, no tree fusing, no drifting. The REPS trick is bootstrap-specific and not relevant to tree search. |
| Goloboff & Morales 2023, "TNT version 1.6, with a graphical interface for MacOS and Linux, including new routines in parallel" (Cladistics) | **None (software release)** — GUI update (GTK3 for Linux/Mac), PVM-based parallel execution, 64-bit binaries. No new search heuristics or scoring algorithms beyond what was in TNT 1.5 and Goloboff 2015. |
| Togkousidis et al. 2023, "Adaptive RAxML-NG: accelerating phylogenetic inference under maximum likelihood using dataset difficulty" (MBE) | **None (ML-specific)** — Predicts dataset "difficulty" via ML model (Pythia), adjusts SPR radius and number of search iterations accordingly. For easy datasets, reduces search effort; for hard ones, increases it. Not transferable: difficulty metrics are ML-specific (based on likelihood convergence patterns), and TNT/TreeSearch already adapts effort based on hitting minimum length independently. |
| Togkousidis, Stamatakis & Gascuel 2025, "Accelerating maximum likelihood phylogenetic inference via early stopping to evade (over-)optimization" (Syst. Biol.) | **None (ML-specific)** — Observes that ML tree search can over-optimize: continuing to search finds marginally better likelihood but no better topology. Suggests stopping early. Not transferable: parsimony scores are discrete (integer steps), so there is no analogous over-optimization — a tree with fewer steps is genuinely better. |
| Azouri et al. 2021, "Harnessing machine learning to guide phylogenetic-tree search algorithms" (Nature Comms.) | **None (ML-specific)** — Trains ML model to predict which SPR moves will improve likelihood, evaluates only top 5%. Not transferable: parsimony evaluation is already cheap (indirect calculation + union construct), so ML prediction overhead would exceed the cost of just evaluating all moves directly. |
| Nguyen et al. 2015, "IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies" (MBE) | **None (ML-specific)** — Search strategy: NNI hill-climbing + stochastic NNI perturbation + candidate set of top trees. Conceptually equivalent to ratchet (perturbation + hill-climbing) but with a smaller neighborhood (NNI vs TBR) and less targeted perturbation (random NNI vs character reweighting). TreeSearch's planned TBR + ratchet + drifting is strictly more powerful. |
