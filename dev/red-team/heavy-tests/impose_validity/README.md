# impose_one_pass / T-327 guard — exhaustive validity harness

**Area 13 (red-team).** Investigates one question about the constraint-repair
guard added by the T-327 fix (commit `6b60f235`):

> Is `build_postorder().postorder.size() == n_internal` equivalent to "the tree
> is fully structurally valid," or can a *corrupt* tree slip that equality — and
> if so, does it reach the advertised `std::bad_alloc` / a wrong answer from the
> public `MaximizeParsimony(constraint = …)`?

## The mechanism

`impose_one_pass()` (`src/ts_constraint.cpp`) repairs constraint violations with
`topology_spr()` SPR moves. Each move is guarded (`src/ts_constraint.cpp:735`):

```cpp
topology_spr(tree, M, above, below);
tree.build_postorder();                                  // src/ts_tree.cpp:75
if ((int)tree.postorder.size() != tree.n_internal) {     // ACCEPT iff ==
  /* restore parent/left/right; */ return false;         // else REVERT
}
```

`build_postorder` is a root-anchored DFS over `left[]`/`right[]` with **no
visited-set** and a cap `if (preorder.size() > n_internal) break;`. So an
*over*-count is always caught (size becomes `n_internal+1`). The only way an
invalid tree can score exactly `n_internal` is a **net-zero corruption**: one
node reached twice (`+1`) compensated by one node never reached (an orphan,
`-1`).

The guard author's comment claims it "guarantees no corrupt tree ever reaches a
DFS helper." This harness tests that completeness claim exhaustively.

## Method (`driver.cpp`)

Standalone C++ (no R / Rcpp / Morphy / SIMD). It compiles the **real**
`src/ts_tree.cpp` (`build_postorder`, `init_from_edge`) and `#include`s the
functions-under-test **extracted verbatim** from `git show
HEAD:src/ts_constraint.cpp` at build time (`extract_funcs.sh` →
`extracted_spr.gen.inc`): `topology_spr`, `collect_edges_in_subtree`,
`collect_edges_outside_subtree`, `compute_node_tips`, `find_maximal_subtrees`.
Verbatim extraction means the code under test **cannot drift** from what the
release ships — if a signature changes, extraction yields nothing and the build
fails loudly.

For every rooted binary tree on `n_tip = 4..8` tips (enumerated exhaustively by
sequential leaf addition; counts checked against `(2n-3)!!` = 15/105/945/10395/135135),
the harness applies `topology_spr` and compares two verdicts:

* **guard**  = `build_postorder().size() == n_internal`  (what the guard sees)
* **truth**  = `full_validity(tree)` — an *independent* check:
  `parent[root]==root`; slots in range and `left!=right`; `parent[]`
  bidirectionally consistent with `left`/`right`; in-degree exactly 1 for every
  non-root node; a root DFS visits every node exactly once (acyclic + no orphan).

`full_validity` deliberately checks `parent[]` too, because `topology_spr` reads
`parent[clip]`/`parent[nx]` at the top of the *next* loop iteration
(`src/ts_constraint.cpp:432,439`) even though the guard and the DFS helpers read
only `left`/`right`.

### Passes

| pass | what it enumerates | role |
|------|--------------------|------|
| **A**  | geometric superset: every non-root clip × every real edge + degenerate `nx/nz/ns` targets | shows the guard is *not* a full structural validator; most triples are unreachable (`collect_edges_*` can't emit them) |
| **C1** | FAITHFUL first-move emission: exact `best_node`=argmin-symdiff, masks, verbatim `find_maximal_subtrees` → move roots, verbatim `collect_edges` targets (all, not RNG) | the exact set of moves emitted on the *first* move of a split |
| **D**  | FAITHFUL model of the **whole** `impose_constraint` (≤ `n_splits+1` passes) for every single-tip-split constraint | the **decisive** reachability + wrong-answer test |
| C2   | crude arbitrary-`M` stale-root superset | off by default (superseded by D) |
| A-perm | Pass A under all internal-id relabellings (`n≤5`) | relabel-invariance belt-and-suspenders |

**Pass D** reproduces the entire `impose_constraint` loop (`:768-773`). Each
`impose_one_pass` computes its move-root list **once** (`:665/669`) and does not
refresh it between moves (so later moves carry a *stale* `M`), while `best_node`
**is** re-anchored before each move (`:721`) and targets are the verbatim
`collect_edges` of the *current* tree. The code picks one **random** target per
move (`:726`); D explores every RNG-reachable outcome — recurse on each
**accepted** target (continuing *past* accepted-but-corrupt trees, exactly as the
real loop does), and take the **skip** branch only when ≥1 target is **rejected**
by the guard. At end-of-list it recurses into the next pass; the pass loop's
own `moves==0`/bail terminations are modelled. The safety cap `n/4+2` (`:675`)
is applied.

Every accepted-invalid tree is classified:
* `invalidity_type` — **type-1** (`left/right` itself corrupt: double-ref +
  orphan) vs **type-2** (`left/right` a valid arborescence, only `parent[]`
  desynced). This distinguishes "the guard admits structural corruption" from
  "the guard fully validates structure; only `parent[]` is stale."
* `has_root_reachable_cycle` (3-colour DFS) — the property that would let the
  **uncapped** downstream traversals (`collect_edges_*`, `compute_dfs_timestamps`)
  loop without bound → `std::bad_alloc`.
* **`g_returnable_corrupt`** — the decisive wrong-answer count: a *final* tree
  (at an `impose_constraint` return point) that is corrupt **and** on which the
  split still maps, so the callers' verify-and-discard (`ts_driven.cpp:1056-1059`,
  `ts_nni_perturb.cpp:119-122`, `ts_parallel.cpp:92-95`) KEEPS it → scored/returned.

## Results

```
n_tip=4..8, all (2n-3)!! trees (counts verified 15/105/945/10395/135135):
Pass A  (geometric superset)                 : ~1.33M   guard is NOT a full validator
Pass C1 (faithful first move)                : 1890     first-move slips at n=8 (0 for n<=7)
Pass D  (full impose_constraint model, n<=6) : 44,496 reachable corrupt trees (9.98M probes)
  of which the split maps ANYWHERE           : 0
  RETURNABLE (final corrupt AND survives verify): 0     => no malformed tree ever returned
invalidity type (all reachable C1+D)         : 100% TYPE-1 (genuine left/right corruption)
root-reachable left/right cycles             : 0        => root-anchored DFS helpers can't loop
parent[] cycles                              : 0        => parent-ascending consumers can't loop
                                                          (tbr_search/reroot_at_tip, run pre-verify)
```

One witness hand-verified end-to-end (n=8, `clip=14 above=8 below=14`): the
`topology_spr` **root-child case** grafts `clip` onto its own parent edge, which
absorbs the sibling into the root, re-points `right[root]` and **orphans** node
13 while making `left[12]==right[12]==14` (a **double-reference**). Net-zero →
`postorder.size()==n_internal` → guard accepts. Confirmed structural (type-1),
acyclic (no `bad_alloc`).

## Verdict

**The hypothesised P1 (`std::bad_alloc`) is refuted; a milder incompleteness is real.**

1. **No crash (airtight, two vectors).** (a) *left/right DFS helpers*:
   `build_postorder`'s cap is *unconditional*, so any move producing a
   **root-reachable** cycle is rejected (`size != n_internal`); every DFS helper
   (`collect_edges_*`, `compute_node_tips`, `compute_dfs_timestamps`) is
   root/postorder-anchored, so none visits a non-root-reachable structure —
   **0** root-reachable left/right cycles across ~1.33M accepted-invalid trees.
   (b) *parent-ascending consumers*: `tbr_search`→`reroot_at_tip`
   (`ts_tbr.cpp:106-107`, an unguarded `while(cur!=root) path.push_back(cur);
   cur=parent[cur]`) runs on the post-impose tree at `ts_nni_perturb.cpp:105`
   **before** verify-and-discard, so a `parent[]` cycle there would `bad_alloc`.
   Empirically: **0** `parent[]` cycles across all accepted-invalid trees. So no
   `std::bad_alloc` from `MaximizeParsimony(constraint=…)` on either vector.
2. **The guard admits genuine TYPE-1 structural corruption** (net-zero double-ref
   + orphan) on reachable moves — first-move at `n≥8` (C1), and within/across
   passes at `n≤6` (D). So the guard-comment "no corrupt tree ever reaches a DFS
   helper" is inaccurate: corrupt trees *do* reach `compute_node_tips` /
   `map_constraint_nodes` — they just don't crash (root-anchored, arrays sized
   `n_node`).
3. **No wrong answer.** Modelling the full `impose_constraint` to its return
   points (single-split, `n≤6`), **`g_returnable_corrupt = 0`**: every corrupt
   tree that reaches a return point fails the callers' `map_constraint_nodes`
   verify-and-discard (all three call sites gate on it). Defense-in-depth
   (unconditional cap + caller re-verify) closes the gap.

### Residual (documented)

The wrong-answer closure (3) is exhaustive for **single-split** constraints at
**`n≤6`** within the probe budget. Not exhaustively closed: multi-split
constraints (a different split could map on a tree corrupted while repairing
another) and `n≥7` continuation of the full model. The **no-crash** result (1) is
per-move and has no such caveat. Any residual, if real, is a
correctness/robustness issue (a malformed tree returned), **not** the
`std::bad_alloc` P1.

## Running it

```sh
bash build_and_run.sh          # n_tip = 4..8 (C1 all n; D n<=6)
bash build_and_run.sh 4 5      # narrower range
```

`d_max` (in `main`, faithful sequential sim depth-gate) and `c2_max` (crude
superset, off) are cost knobs. `n_tip<=7` is local-safe seconds-to-minutes;
`n_tip=8` Pass D branching is gated off by default.
