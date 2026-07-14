#include "ts_tbr.h"
#include "ts_fitch.h"
#include "ts_collapsed.h"
#include "ts_rng.h"
#include "ts_tabu.h"
#include "ts_splits.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>
#include <unordered_set>
#include <climits>
#include <cmath>
#include <cstring>
#include <chrono>

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

namespace ts {

// --- Fast hash for virtual_prelim deduplication (Phase 3A) ---
// Word-at-a-time multiply-xor hash (faster than byte-by-byte FNV-1a).

static uint64_t fast_hash(const uint64_t* data, int n_words) {
  uint64_t hash = 14695981039346656037ULL;
  for (int i = 0; i < n_words; ++i) {
    hash ^= data[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

// Reusable open-addressed hash set for per-clip virtual_prelim dedup (Phase 3A).
// S-PROF round 3 / Tier 2: replaces a per-clip std::unordered_set<uint64_t>,
// which heap-allocated a bucket array plus one node per insert on every clip of
// the TBR hot loop. Declared once as a plain local before the clip loop — NOT
// static thread_local: insert() runs per reroot candidate, and MinGW resolves
// thread_local in a loaded DLL via emutls (a function call per access), which
// cancelled the win in an earlier thread_local variant. A plain local in
// tbr_search is already per-thread-safe (each thread owns its call frame) and
// has zero TLS cost. Generation stamping makes reset() O(1): a slot counts as
// occupied only while its stamp equals the current generation, so per-clip
// clearing just bumps the generation. Dedup semantics are identical to
// unordered_set<uint64_t> — insert() returns true iff this exact 64-bit key was
// not already present this generation (distinct keys never merge; identical keys
// always do; collisions resolve by linear probing). fast_hash is FNV-1a (weak
// low-bit avalanche), so probe from Fibonacci-mixed high bits, not key & mask.
struct VpHashSet {
  std::vector<uint64_t> keys;
  std::vector<uint32_t> stamp;
  uint32_t cur = 0;
  int shift = 64;
  size_t mask = 0;

  void reset(size_t expected) {
    size_t cap = 16;
    while (cap < (expected + 1) * 2) cap <<= 1;   // load factor < 0.5
    if (keys.size() < cap) {
      keys.assign(cap, 0);
      stamp.assign(cap, 0);
      cur = 0;                                     // all stamps 0 after grow
    }
    mask = keys.size() - 1;
    shift = 64 - popcount64(mask);                 // keep top log2(cap) bits
    if (++cur == 0) {                              // generation wrapped (2^32)
      std::fill(stamp.begin(), stamp.end(), 0);
      cur = 1;
    }
  }

  // True if newly inserted; false if `key` was already present this generation.
  bool insert(uint64_t key) {
    size_t i = (key * 0x9E3779B97F4A7C15ULL) >> shift;
    while (stamp[i] == cur) {
      if (keys[i] == key) return false;
      i = (i + 1) & mask;
    }
    keys[i] = key;
    stamp[i] = cur;
    return true;
  }
};

// --- Helpers (file-local) ---

static double full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

// Re-root the tree so tip `t` is a direct child of the root pseudo-node n_tip.
// Parsimony length is root-invariant, so this only changes the representation
// (which edges are clippable and where the root edge sits) — it lets the search
// reach moves the current rooting hides.  Rebuilds postorder; does NOT refresh
// Fitch state arrays, so the caller must full_rescore() afterwards.
// Generalises reroot_at_tip0() in ts_fuse.cpp to an arbitrary tip.
// Declared in ts_tbr.h (used by the output-collapse kernel in ts_rcpp.cpp).
void reroot_at_tip(TreeState& tree, int t) {
  const int n_tip = tree.n_tip;
  const int root = n_tip;
  if (tree.parent[t] == root) return;            // already a child of root

  std::vector<int> path;                          // t's parent .. child-of-root
  int cur = tree.parent[t];
  while (cur != root) { path.push_back(cur); cur = tree.parent[cur]; }
  std::reverse(path.begin(), path.end());

  const int path_len = static_cast<int>(path.size());
  const int root_ni = 0;
  const int root_other = (tree.left[root_ni] == path[0])
                             ? tree.right[root_ni] : tree.left[root_ni];
  for (int i = 0; i < path_len; ++i) {
    const int node = path[i];
    const int ni = node - n_tip;
    const int toward = (i + 1 < path_len) ? path[i + 1] : t;
    const int replacement = (i == 0) ? root_other : path[i - 1];
    if (tree.left[ni] == toward) tree.left[ni] = replacement;
    else                          tree.right[ni] = replacement;
    tree.parent[replacement] = node;
  }
  const int last_path = path[path_len - 1];
  tree.left[root_ni] = t;
  tree.right[root_ni] = last_path;
  tree.parent[t] = root;
  tree.parent[last_path] = root;
  tree.build_postorder();
}

// Collect (parent, child) edge pairs reachable from root of main tree.
static void collect_main_edges(
    const TreeState& tree,
    std::vector<std::pair<int,int>>& edges)
{
  edges.clear();
  // Reusable per-thread DFS stack (Tier 1): avoids a heap alloc per clip.
  std::vector<int> stack;
  stack.clear();
  stack.push_back(tree.n_tip);

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) continue;

    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    edges.push_back({node, lc});
    edges.push_back({node, rc});
    stack.push_back(lc);
    stack.push_back(rc);
  }
}

// Collect (parent, child) edge pairs within a subtree rooted at subtree_root.
static void collect_subtree_edges(
    const TreeState& tree, int subtree_root,
    std::vector<std::pair<int,int>>& edges)
{
  edges.clear();
  if (subtree_root < tree.n_tip) return;

  // Reusable per-thread DFS stack (Tier 1): avoids a heap alloc per clip.
  std::vector<int> stack;
  stack.clear();
  stack.push_back(subtree_root);

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) continue;

    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    edges.push_back({node, lc});
    edges.push_back({node, rc});
    stack.push_back(lc);
    stack.push_back(rc);
  }
}

// Compute fitch_join of two state sets.
static void fitch_join_states(
    const uint64_t* state_a,
    const uint64_t* state_b,
    uint64_t* out,
    const DataSet& ds)
{
  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int offset = ds.block_word_offset[b];

    uint64_t any_isect = 0;
    for (int s = 0; s < blk.n_states; ++s) {
      any_isect |= (state_a[offset + s] & state_b[offset + s]);
    }
    uint64_t no_isect = ~any_isect & blk.active_mask;

    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t isect = state_a[offset + s] & state_b[offset + s];
      uint64_t uni   = state_a[offset + s] | state_b[offset + s];
      out[offset + s] = (isect & any_isect) | (uni & no_isect);
    }
  }
}

// Compute from_above states for all nodes in a clipped subtree.
static void compute_from_above(
    const TreeState& tree, const DataSet& ds,
    int subtree_root,
    std::vector<uint64_t>& from_above)
{
  int tw = tree.total_words;

  // Reusable per-thread scratch (Tier 1): avoids two heap allocs per clip.
  std::vector<int> preorder;
  preorder.clear();
  {
    std::vector<int> stack;
    stack.clear();
    stack.push_back(subtree_root);
    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();
      if (node < tree.n_tip) continue;
      preorder.push_back(node);
      int ni = node - tree.n_tip;
      stack.push_back(tree.right[ni]);
      stack.push_back(tree.left[ni]);
    }
  }

  if (preorder.empty()) return;

  // Root's children: from_above = sibling's prelim
  int root = preorder[0];
  int ri = root - tree.n_tip;
  int lc = tree.left[ri];
  int rc = tree.right[ri];

  for (int w = 0; w < tw; ++w) {
    from_above[static_cast<size_t>(lc) * tw + w] =
        tree.prelim[static_cast<size_t>(rc) * tw + w];
    from_above[static_cast<size_t>(rc) * tw + w] =
        tree.prelim[static_cast<size_t>(lc) * tw + w];
  }

  // Remaining internal nodes
  for (size_t i = 1; i < preorder.size(); ++i) {
    int node = preorder[i];
    int ni = node - tree.n_tip;
    int lc_i = tree.left[ni];
    int rc_i = tree.right[ni];

    fitch_join_states(
        &from_above[static_cast<size_t>(node) * tw],
        &tree.prelim[static_cast<size_t>(rc_i) * tw],
        &from_above[static_cast<size_t>(lc_i) * tw], ds);
    fitch_join_states(
        &from_above[static_cast<size_t>(node) * tw],
        &tree.prelim[static_cast<size_t>(lc_i) * tw],
        &from_above[static_cast<size_t>(rc_i) * tw], ds);
  }
}

// --- Topology snapshot for safe undo ---

struct TopoSnapshot {
  std::vector<int> parent;
  std::vector<int> left;
  std::vector<int> right;
};

static void save_topology(const TreeState& tree, TopoSnapshot& snap) {
  snap.parent = tree.parent;
  snap.left = tree.left;
  snap.right = tree.right;
}

static void restore_topology(TreeState& tree, const TopoSnapshot& snap) {
  tree.parent = snap.parent;
  tree.left = snap.left;
  tree.right = snap.right;
}

// --- Full state snapshot for undo without rescore ---

struct StateSnapshot {
  std::vector<uint64_t> prelim;
  std::vector<uint64_t> final_;
  std::vector<uint64_t> local_cost;
  std::vector<uint64_t> down2;
  std::vector<uint64_t> subtree_actives;
  std::vector<int> postorder;
  bool has_na_arrays;

  void allocate(const TreeState& tree, bool has_na) {
    size_t state_sz = static_cast<size_t>(tree.n_node) * tree.total_words;
    size_t cost_sz = static_cast<size_t>(tree.n_node) * tree.n_blocks;
    prelim.resize(state_sz);
    final_.resize(state_sz);
    local_cost.resize(cost_sz);
    has_na_arrays = has_na;
    if (has_na) {
      down2.resize(state_sz);
      subtree_actives.resize(state_sz);
    }
    postorder.resize(tree.postorder.size());
  }

  void save(const TreeState& tree) {
    size_t state_bytes = prelim.size() * sizeof(uint64_t);
    size_t cost_bytes = local_cost.size() * sizeof(uint64_t);
    std::memcpy(prelim.data(), tree.prelim.data(), state_bytes);
    std::memcpy(final_.data(), tree.final_.data(), state_bytes);
    std::memcpy(local_cost.data(), tree.local_cost.data(), cost_bytes);
    if (has_na_arrays) {
      std::memcpy(down2.data(), tree.down2.data(), state_bytes);
      std::memcpy(subtree_actives.data(), tree.subtree_actives.data(),
                   state_bytes);
    }
    std::memcpy(postorder.data(), tree.postorder.data(),
                 tree.postorder.size() * sizeof(int));
  }

  void restore(TreeState& tree) const {
    size_t state_bytes = prelim.size() * sizeof(uint64_t);
    size_t cost_bytes = local_cost.size() * sizeof(uint64_t);
    std::memcpy(tree.prelim.data(), prelim.data(), state_bytes);
    std::memcpy(tree.final_.data(), final_.data(), state_bytes);
    std::memcpy(tree.local_cost.data(), local_cost.data(), cost_bytes);
    if (has_na_arrays) {
      std::memcpy(tree.down2.data(), down2.data(), state_bytes);
      std::memcpy(tree.subtree_actives.data(), subtree_actives.data(),
                   state_bytes);
    }
    // Restore postorder size AND data (clip may have shrunk the vector)
    tree.postorder.resize(postorder.size());
    std::memcpy(tree.postorder.data(), postorder.data(),
                 postorder.size() * sizeof(int));
  }
};

// --- Topology validation (debug, catches bugs before they crash R) ---
// [[maybe_unused]] because the sole call site (in tbr_search) is gated under
// NDEBUG (Tier 3a) — in release builds this function is never referenced.

[[maybe_unused]] static bool validate_topology(const TreeState& tree) {
  int root = tree.n_tip;

  // DFS from root with visited check to detect cycles
  std::vector<bool> visited(tree.n_node, false);
  std::vector<int> stack;
  stack.push_back(root);
  int n_visited = 0;

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    if (node < 0 || node >= tree.n_node) return false;
    if (visited[node]) return false;  // cycle!
    visited[node] = true;
    ++n_visited;

    if (node >= tree.n_tip) {
      int ni = node - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }

  return (n_visited == tree.n_node);
}

// Apply a TBR move directly to the tree topology.
//
// Detaches the subtree at clip_node, optionally reroots it at
// (reroot_parent, reroot_child), and regrafts at (above, below).
//
// The caller must save a topology snapshot first and full_rescore after.
static bool apply_tbr_move(
    TreeState& tree,
    int clip_node,
    int reroot_parent, int reroot_child,  // -1 if SPR (no reroot)
    int above, int below)
{
  int nx = tree.parent[clip_node];  // clip parent (spare node)
  int nz = tree.parent[nx];        // grandparent
  int nxi = nx - tree.n_tip;
  int ns;                           // sibling of clip_node
  ns = (tree.left[nxi] == clip_node) ? tree.right[nxi] : tree.left[nxi];

  // Step 1: Detach — connect sibling to grandparent, freeing nx
  tree.parent[ns] = nz;
  if (nz >= tree.n_tip) {
    int nzi = nz - tree.n_tip;
    if (tree.left[nzi] == nx) {
      tree.left[nzi] = ns;
    } else {
      tree.right[nzi] = ns;
    }
  }

  // Step 2: Reroot the clipped subtree if needed
  int new_subtree_root = clip_node;
  if (reroot_parent >= 0 && reroot_parent != clip_node) {
    // Find path from clip_node DOWN to reroot_parent via DFS
    std::vector<int> path;
    {
      // Build parent map within the subtree via DFS
      std::vector<int> dfs_stack;
      std::vector<int> sub_parent(tree.n_node, -1);
      dfs_stack.push_back(clip_node);

      while (!dfs_stack.empty()) {
        int node = dfs_stack.back();
        dfs_stack.pop_back();
        if (node == reroot_parent) break;
        if (node < tree.n_tip) continue;
        int ni = node - tree.n_tip;
        int lc = tree.left[ni];
        int rc = tree.right[ni];
        sub_parent[lc] = node;
        sub_parent[rc] = node;
        dfs_stack.push_back(lc);
        dfs_stack.push_back(rc);
      }

      // Reconstruct path: reroot_parent -> ... -> clip_node
      int cur = reroot_parent;
      while (cur != clip_node && cur >= 0) {
        path.push_back(cur);
        cur = sub_parent[cur];
      }
      if (cur < 0) return false;  // reroot_parent not in subtree!
      path.push_back(clip_node);
      std::reverse(path.begin(), path.end());
      // path = [clip_node, ..., reroot_parent]
    }

    if (path.size() < 2) return false;

    // Reverse parent-child links along the path.
    // After this, reroot_parent becomes the subtree root.
    for (size_t j = 0; j + 1 < path.size(); ++j) {
      int A = path[j];       // current parent
      int B = path[j + 1];   // current child (will become parent)

      int ai = A - tree.n_tip;
      int bi = B - tree.n_tip;

      // Find B's child NOT on the path (B_off_path)
      int B_off_path;
      if (j + 2 < path.size()) {
        int next_on_path = path[j + 2];
        B_off_path = (tree.left[bi] == next_on_path)
                     ? tree.right[bi] : tree.left[bi];
      } else {
        // B is reroot_parent; off-path child is non-reroot_child
        B_off_path = (tree.left[bi] == reroot_child)
                     ? tree.right[bi] : tree.left[bi];
      }

      // In A: replace child B with B_off_path
      if (tree.left[ai] == B) {
        tree.left[ai] = B_off_path;
      } else {
        tree.right[ai] = B_off_path;
      }
      tree.parent[B_off_path] = A;

      // In B: replace B_off_path with A
      if (tree.left[bi] == B_off_path) {
        tree.left[bi] = A;
      } else {
        tree.right[bi] = A;
      }
      tree.parent[A] = B;
    }

    new_subtree_root = reroot_parent;
  }

  // Step 3: Regraft — insert nx between (above, below)
  if (above >= tree.n_tip) {
    int ai = above - tree.n_tip;
    if (tree.left[ai] == below) {
      tree.left[ai] = nx;
    } else {
      tree.right[ai] = nx;
    }
  }
  tree.parent[nx] = above;

  tree.left[nxi] = new_subtree_root;
  tree.right[nxi] = below;
  tree.parent[new_subtree_root] = nx;
  tree.parent[below] = nx;

  return true;
}

// Reroot the subtree rooted at `frag_root` so `reroot_parent` becomes its new
// top, by reversing parent/child links along the path frag_root..reroot_parent.
// `reroot_child` selects reroot_parent's off-path child. Touches ONLY
// fragment-internal links; the CALLER must set the new top's parent pointer.
// Returns the new top (reroot_parent) or -1 on failure.  Faithful copy of
// apply_tbr_move() Step 2 — keep the (reroot_parent,reroot_child)/from_above
// pairing identical so the root-edge scan score equals the applied score.
// (TODO: de-duplicate with apply_tbr_move once both are settled.)
static int reroot_fragment(TreeState& tree, int frag_root,
                           int reroot_parent, int reroot_child) {
  if (reroot_parent < 0 || reroot_parent == frag_root) return frag_root;

  std::vector<int> path;
  {
    std::vector<int> dfs_stack;
    std::vector<int> sub_parent(tree.n_node, -1);
    dfs_stack.push_back(frag_root);
    while (!dfs_stack.empty()) {
      int node = dfs_stack.back();
      dfs_stack.pop_back();
      if (node == reroot_parent) break;
      if (node < tree.n_tip) continue;
      int ni = node - tree.n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];
      sub_parent[lc] = node;
      sub_parent[rc] = node;
      dfs_stack.push_back(lc);
      dfs_stack.push_back(rc);
    }
    int cur = reroot_parent;
    while (cur != frag_root && cur >= 0) {
      path.push_back(cur);
      cur = sub_parent[cur];
    }
    if (cur < 0) return -1;          // reroot_parent not in subtree
    path.push_back(frag_root);
    std::reverse(path.begin(), path.end());
  }
  if (path.size() < 2) return -1;

  for (size_t j = 0; j + 1 < path.size(); ++j) {
    int A = path[j];
    int B = path[j + 1];
    int ai = A - tree.n_tip;
    int bi = B - tree.n_tip;
    int B_off_path;
    if (j + 2 < path.size()) {
      int next_on_path = path[j + 2];
      B_off_path = (tree.left[bi] == next_on_path)
                   ? tree.right[bi] : tree.left[bi];
    } else {
      B_off_path = (tree.left[bi] == reroot_child)
                   ? tree.right[bi] : tree.left[bi];
    }
    if (tree.left[ai] == B) tree.left[ai] = B_off_path;
    else tree.right[ai] = B_off_path;
    tree.parent[B_off_path] = A;
    if (tree.left[bi] == B_off_path) tree.left[bi] = A;
    else tree.right[bi] = A;
    tree.parent[A] = B;
  }
  return reroot_parent;
}

// IW/NA variant of the root-edge enumeration: the EW additive split
// (base_split + Fitch join) is invalid for implied weights (concave) and
// inapplicables (3-pass), so here we score each candidate by APPLYING it and
// full-rescoring (exact for any scorer).  Same (L-rooting x R-rooting)
// enumeration as the EW path; only the scoring/acceptance differs.  Lazy
// (convergence only).  Snapshot the clean topology once, then
// restore->apply->rescore per candidate so each starts clean.
static bool try_root_edge_moves_rescore(TreeState& tree, const DataSet& ds,
                                        double& best_score) {
  const int n_tip = tree.n_tip;
  const int cL = tree.left[0];
  const int cR = tree.right[0];
  const double eps = std::isfinite(ds.concavity) ? 1e-10 : 0.0;

  best_score = full_rescore(tree, ds);

  // Rerooting metadata for each half: identity {-1,-1} plus each internal edge
  // (sp,sc) with sp != croot (same distinct-rootings set as the EW path).
  std::vector<std::pair<int,int>> metaL, metaR, edges;
  auto build_meta = [&](int croot, std::vector<std::pair<int,int>>& meta) {
    meta.clear();
    meta.push_back({-1, -1});
    if (croot < n_tip) return;
    collect_subtree_edges(tree, croot, edges);
    for (auto& [sp, sc] : edges) if (sp != croot) meta.push_back({sp, sc});
  };
  build_meta(cL, metaL);
  build_meta(cR, metaR);
  const int nL = static_cast<int>(metaL.size());
  const int nR = static_cast<int>(metaR.size());

  TopoSnapshot snap;
  save_topology(tree, snap);

  auto rejoin = [&](int li, int ri) -> bool {
    int topL = (metaL[li].first < 0) ? cL
               : reroot_fragment(tree, cL, metaL[li].first, metaL[li].second);
    int topR = (metaR[ri].first < 0) ? cR
               : reroot_fragment(tree, cR, metaR[ri].first, metaR[ri].second);
    if (topL < 0 || topR < 0) return false;
    tree.left[0] = topL; tree.right[0] = topR;
    tree.parent[topL] = n_tip; tree.parent[topR] = n_tip;
    tree.build_postorder();
    return true;
  };

  double best_cand = best_score;
  int bestLi = -1, bestRi = -1;
  for (int li = 0; li < nL; ++li) {
    for (int ri = 0; ri < nR; ++ri) {
      if (li == 0 && ri == 0) continue;            // identity = current tree
      if (rejoin(li, ri)) {
        double cand = full_rescore(tree, ds);
        if (cand < best_cand - eps) { best_cand = cand; bestLi = li; bestRi = ri; }
      }
      restore_topology(tree, snap);                // back to clean state
    }
  }

  // After the scan loop the topology is the restored clean one, but its
  // postorder is stale (the last candidate's rejoin overwrote it and
  // restore_topology does NOT restore postorder).  score_tree walks the
  // postorder, so every full_rescore below MUST rebuild it first or best_score
  // silently drifts from the returned tree (the bug that made IW look
  // incomplete).  The improving branch rebuilds it inside rejoin().
  if (bestLi < 0) {
    tree.build_postorder();
    best_score = full_rescore(tree, ds);
    return false;
  }

  if (!rejoin(bestLi, bestRi)) {              // defensive: restore clean state
    restore_topology(tree, snap);
    tree.build_postorder();
    best_score = full_rescore(tree, ds);
    return false;
  }
  best_score = full_rescore(tree, ds);        // rejoin() already rebuilt postorder
  return true;
}

// Direct in-pass enumeration of the ONE unrooted edge the root pseudo-node
// n_tip sits on (cL—cR), which the rooted clip loop structurally skips
// (parent==n_tip guard).  This is the whole residual that the physical-reroot
// sweep used to cover at O(n_tip) full rescores — here it costs one edge.
//
// PURE-EW ONLY: Fitch length is root-invariant, so each half-fragment's
// internal length is constant under rerooting, and a root-edge TBR move's
// length decomposes exactly as
//     base_split + fitch_join(stateL, stateR),
// base_split = best_score - join(prelim[cL], prelim[cR])  (constant),
// stateL/stateR = each half rerooted at the chosen connection edge.  This
// additive split does NOT hold for IW (concave) or NA (3-pass), so the caller
// gates this on ew_directional and routes IW/NA to the physical-reroot fallback.
//
// Reuses compute_from_above + fitch_indirect_length_cached exactly as the
// normal rerooting scan.  Applies the best reconnection (if it beats the
// current root join) via reroot_fragment on each half, rejoining at n_tip, and
// returns true; else leaves the tree unchanged and returns false.
static bool try_root_edge_moves(TreeState& tree, const DataSet& ds,
                                double& best_score, bool ew_directional) {
  // IW/NA cannot use the additive base_split decomposition below; route them to
  // the apply+rescore variant (same enumeration, exact scoring).
  if (!ew_directional) return try_root_edge_moves_rescore(tree, ds, best_score);

  const int n_tip = tree.n_tip;
  const int tw = tree.total_words;
  const int cL = tree.left[0];      // root internal index = n_tip - n_tip = 0
  const int cR = tree.right[0];

  // Refresh states so prelim[] is current for cL/cR and every fragment node.
  best_score = full_rescore(tree, ds);

  const uint64_t* pL = &tree.prelim[static_cast<size_t>(cL) * tw];
  const uint64_t* pR = &tree.prelim[static_cast<size_t>(cR) * tw];
  const int rootjoin = fitch_indirect_length_cached(pL, pR, ds, INT_MAX);
  const double base_split = best_score - rootjoin;

  std::vector<uint64_t> from_above;
  std::vector<uint64_t> rowsL, rowsR;
  std::vector<std::pair<int,int>> metaL, metaR;
  std::vector<std::pair<int,int>> edges;
  from_above.assign(static_cast<size_t>(tree.n_node) * tw, 0ULL);

  // Build the rerooting state-sets for one half.  Row 0 = identity (the half's
  // current root prelim); rows 1.. = join(from_above[sc], prelim[sc]) for each
  // internal edge (sp,sc) with sp != croot.  The sp==croot edges are skipped:
  // from_above[child] = sibling prelim there, so their join reproduces the
  // identity state-set (this is the normal scan's sp==clip_node guard, and it
  // yields exactly the 2k-3 distinct rootings of a k-tip half).
  auto build_half = [&](int croot, std::vector<uint64_t>& rows,
                        std::vector<std::pair<int,int>>& meta) {
    rows.clear();
    meta.clear();
    const uint64_t* prc = &tree.prelim[static_cast<size_t>(croot) * tw];
    rows.insert(rows.end(), prc, prc + tw);          // row 0: identity
    meta.push_back({-1, -1});
    if (croot < n_tip) return;                        // single-tip half
    compute_from_above(tree, ds, croot, from_above);
    collect_subtree_edges(tree, croot, edges);
    for (auto& [sp, sc] : edges) {
      if (sp == croot) continue;
      size_t base = rows.size();
      rows.resize(base + tw);
      fitch_join_states(&from_above[static_cast<size_t>(sc) * tw],
                        &tree.prelim[static_cast<size_t>(sc) * tw],
                        &rows[base], ds);
      meta.push_back({sp, sc});
    }
  };
  build_half(cL, rowsL, metaL);
  build_half(cR, rowsR, metaR);

  const int nL = static_cast<int>(metaL.size());
  const int nR = static_cast<int>(metaR.size());

  int best_join = rootjoin;
  int bestLi = -1, bestRi = -1;
  for (int li = 0; li < nL; ++li) {
    const uint64_t* sL = &rowsL[static_cast<size_t>(li) * tw];
    for (int ri = 0; ri < nR; ++ri) {
      if (li == 0 && ri == 0) continue;               // identity = current tree
      const uint64_t* sR = &rowsR[static_cast<size_t>(ri) * tw];
      int j = fitch_indirect_length_cached(sL, sR, ds, best_join);
      if (j < best_join) { best_join = j; bestLi = li; bestRi = ri; }
    }
  }

  if (bestLi < 0 || best_join >= rootjoin) return false;   // no improvement

  // Apply: reroot each half to its chosen connection edge, rejoin at n_tip.
  int topL = cL, topR = cR;
  if (metaL[bestLi].first >= 0)
    topL = reroot_fragment(tree, cL, metaL[bestLi].first, metaL[bestLi].second);
  if (metaR[bestRi].first >= 0)
    topR = reroot_fragment(tree, cR, metaR[bestRi].first, metaR[bestRi].second);
  if (topL < 0 || topR < 0) return false;                  // defensive

  tree.left[0] = topL;
  tree.right[0] = topR;
  tree.parent[topL] = n_tip;
  tree.parent[topR] = n_tip;
  tree.build_postorder();

  const double actual = full_rescore(tree, ds);

  // Degree-2 tripwire (env TS_TBR_ASSERT): the applied length MUST equal the
  // scan's prediction.  Any off-by-one in the reversal or the n_tip rewire
  // trips here on the first accepted move.
  if (std::getenv("TS_TBR_ASSERT")) {
    double predicted = base_split + best_join;
    if (std::fabs(actual - predicted) > 0.5) {
      REprintf("ROOT-EDGE MISMATCH: actual=%.1f predicted=%.1f (base=%.1f join=%d)\n",
               actual, predicted, base_split, best_join);
    }
  }

  best_score = actual;
  return true;
}

// FNV-1a hash over each internal node's (min,max) child pair.  Sorting the pair
// canonicalizes only the LEFT/RIGHT child order WITHIN a node (the same rooted
// tree hashes identically however its two children happen to be stored); it does
// NOT make the hash root-independent.  Different rootings renumber the internal
// nodes and flip parent/child directions along the reroot path, so they hash
// differently (verified: two rootings of one unrooted tree differ at 55/73
// nodes).  That root-DEPENDENCE is REQUIRED, not a limitation: exact_verify_sweep
// is itself root-dependent — it skips root-child clips, so each rooting has a
// different neighbourhood and a separately-valid optimum verdict (the residual
// completeness gap, task #19) — so each rooting MUST be cached separately.  Safe
// because the NA path never physically reroots mid-search (the legacy reroot
// sweep is TS_PHYS_REROOT-only) and restores preserve the exact rooting, so a
// converged optimum is always re-verified at the rooting it was cached under.
// Do NOT "canonicalize" this to a root-independent hash unless/until
// exact_verify_sweep is made root-complete (task #19), or cross-rooting hits
// would suppress improvers one rooting finds and another misses.
static inline uint64_t tree_topo_hash(const TreeState& tree) {
  uint64_t h = 14695981039346656037ULL;
  for (int i = 0; i < tree.n_tip - 1; ++i) {
    int a = tree.left[i], b = tree.right[i];
    if (a > b) { int t = a; a = b; b = t; }
    h ^= (uint64_t)a * 2246822519ULL;
    h *= 1099511628211ULL;
    h ^= (uint64_t)b * 2654435761ULL;
    h *= 1099511628211ULL;
  }
  return h;
}

// Dataset fingerprint: mix n_tips, n_blocks, and every tip_states word so that
// any dataset change (adding a char, changing a tip state) produces a new key.
// The cache is cleared whenever the fingerprint changes.
static inline uint64_t ds_fingerprint(const DataSet& ds) {
  uint64_t h = (uint64_t)ds.n_tips * 2654435761ULL
             ^ (uint64_t)ds.n_blocks * 2246822519ULL;
  for (uint64_t v : ds.tip_states) { h ^= v; h *= 1099511628211ULL; }
  return h;
}

// Weighting-regime fingerprint.  The ratchet mutates the per-block active_mask
// and upweight_mask and the IW pattern_freq IN PLACE on the live DataSet (see
// save_perturb_state / restore_perturb_state in ts_ratchet.cpp — these three
// fields ARE the "scoring state that varies mid-search"), and NA scoring reads
// all three (ts_fitch_na*.h).  exact_verify therefore scores a topology
// differently in the perturbed vs base regime, so the regime must be part of
// the cache key.  Mixing the same three fields here covers the regime by
// construction.
static inline uint64_t weight_fingerprint(const DataSet& ds) {
  uint64_t h = 14695981039346656037ULL;
  for (const auto& blk : ds.blocks) {
    h ^= blk.active_mask;   h *= 1099511628211ULL;
    h ^= blk.upweight_mask; h *= 1099511628211ULL;
  }
  for (int f : ds.pattern_freq) { h ^= (uint64_t)(uint32_t)f; h *= 1099511628211ULL; }
  return h;
}

// Shared cache key for exact_verify_sweep's optimum memoization (declared in
// ts_tbr.h).  Both the cache below AND the regression probe
// (ts_ev_cache_key_probe -> test-ts-na-evcache.R) call THIS one function, so
// dropping any term here — topology, dataset, or weighting regime — is caught
// by a deterministic test rather than only by a silent search-quality drop.
uint64_t exact_verify_cache_key(const TreeState& tree, const DataSet& ds) {
  return tree_topo_hash(tree) ^ ds_fingerprint(ds) ^ weight_fingerprint(ds);
}

// Exact full-neighbourhood TBR verification, for scorers whose indirect scan is
// only APPROXIMATE (inapplicables / NA).  Under Brazeau's three-pass the
// divided+reconnect decomposition is not exact (the clipped subtree's internal
// count is attachment-dependent — down2 reads whole-tree uppass context), so the
// fast inner clip loop can declare convergence while improving moves remain
// (confirmed by dev/benchmarks/tbr_oracle_na.R: both the direct scan AND the
// physical-reroot path leave improving NA neighbours).  At convergence this
// sweeps every NON-root edge's TBR neighbourhood (the 2n-4 edges not incident on
// the display root), scoring each candidate EXACTLY via apply_tbr_move +
// full_rescore, applies the first strict improver found (first-improvement: the
// cheap approximate loop re-climbs between calls), and returns true.  The clip
// loop structurally skips root-child clips, so the ONE root edge (cL-cR) is
// enumerated separately at the optimum exit via try_root_edge_moves_rescore;
// together they cover all 2n-3 unrooted edges, so a false return is a genuine
// unrooted-TBR optimum.  The clip enumeration is clip_node x {identity + fragment
// rerootings} x divided-tree regraft edges, with the regraft edges built directly
// from the unclipped tree (every (parent[c], c) with c outside the clipped
// subtree, plus the merged (nz, ns) edge that replaces nz-nx-ns), so
// apply_tbr_move re-clips from the original tree exactly as the scan's accept
// path does.  EW/IW never use this (their scan is exact); the caller gates it on
// has_na, so the default and EW/IW paths are byte-identical.
static bool exact_verify_sweep(TreeState& tree, const DataSet& ds,
                               double& best_score) {
  const double eps = std::isfinite(ds.concavity) ? 1e-9 : 0.5;
  // Plain locals (not thread_local): MinGW emutls thread_local teardown across
  // std::thread spawn/exit corrupted the heap on the parallel path.  Each worker
  // owns its call frame, so plain locals are per-thread-safe; the per-clip
  // (re)allocation is in the noise (measured <=1.6% on 88-tip data).
  TopoSnapshot snap;
  std::vector<std::pair<int,int>> sub_edges;
  std::vector<char> in_sub;
  std::vector<int> dfs, marked;

  // TS_NA_INCR_AUDIT: validate the incremental dirty-rescore of each candidate
  // against full_rescore (the prospective fast path for exact_verify). Decisions
  // still use full_rescore — this only cross-checks the incremental score, so it
  // is SAFE while developing. The incremental path uses the 3-seed dirty passes
  // (nz, nx, clip_node) + Pass3 + (IW: extract+compute_iw | EW: +ew_offset),
  // mirroring the SPR accept path but covering reroot via the clip_node seed.
  const bool na_incr_audit = std::getenv("TS_NA_INCR_AUDIT") != nullptr;
  // Incremental dirty rescore is the PRODUCTION DEFAULT for exact_verify
  // candidates (validated byte-identical to legacy full_rescore: per-candidate
  // audit + 180/180 roster climbs + 40/40 mission cells; ~25-30% native-NA
  // mission wall). Kill-switch TS_NA_NOINCR restores the legacy full_rescore
  // path. Disabled in audit mode (which uses full_rescore for decisions).
  // getenv is read per exact_verify call (per-convergence, not per-candidate).
  const bool na_incr = std::getenv("TS_NA_NOINCR") == nullptr && !na_incr_audit;
  // Mirror score_tree/fitch_score_ew's dispatch exactly (NOT isfinite(concavity)):
  // weighted modes (IW/XPIWE/PROFILE) extract per-pattern steps + compute_weighted;
  // EW adds ew_offset. HSJ/XFORM have extra DP/Sankoff scoring score_tree wraps
  // around fitch_score_ew, so the incremental path can't reproduce them -> skip.
  const bool ev_weighted = (ds.scoring_mode == ScoringMode::IW ||
                            ds.scoring_mode == ScoringMode::XPIWE ||
                            ds.scoring_mode == ScoringMode::PROFILE);
  const bool ev_incr_ok = ev_weighted || ds.scoring_mode == ScoringMode::EW;
  const double na_incr_tol = ev_weighted ? 1e-6 : 0.5;
  std::vector<int> incr_char_steps;

  tree.build_postorder();
  best_score = full_rescore(tree, ds);   // sync to the current (converged) tree

  // Topology cache: exact_verify is a pure function of (topology, dataset,
  // weighting regime).  A FALSE result for topology T means every unrooted-TBR
  // neighbour is no better UNDER THE CURRENT WEIGHTS, and no state changes
  // between calls could make that untrue.  The ratchet mutates the weighting in
  // place (active_mask / upweight_mask / pattern_freq; ts_ratchet.cpp) and runs
  // NA TBR under both perturbed and base weights within one cycle, so the
  // weighting MUST be in the key — otherwise a base-regime "optimal" verdict
  // would be reused during a perturbed pass (or vice-versa), silently skipping
  // the very improving moves the ratchet exists to find.  Key on
  // hash(child-pairs) XOR dataset-fingerprint XOR weight-fingerprint; only the
  // dataset-fingerprint is the clear-trigger (a true dataset switch), so
  // base-regime entries survive across perturbation excursions and are reused.
  // Memoization lives on the (per-worker) DataSet, NOT a function-local
  // thread_local — MinGW emutls thread_local teardown across std::thread
  // spawn/exit corrupted the heap.  ds_local has the same per-thread,
  // cross-replicate lifetime, so the cache's persistence is unchanged.  See
  // the evs_false_cache / evs_last_fp comment in ts_data.h.
  std::unordered_set<uint64_t>& evs_false_cache = ds.evs_false_cache;
  const uint64_t fp = ds_fingerprint(ds);   // clear-trigger: a true dataset switch
  if (fp != ds.evs_last_fp) { evs_false_cache.clear(); ds.evs_last_fp = fp; }
  const uint64_t cache_key = exact_verify_cache_key(tree, ds);

  // TS_EV_AUDIT (dev/bench only): distrust cache hits.  On a hit, run the full
  // sweep anyway and abort if it finds an improver the cached FALSE claimed
  // absent — the live tripwire for weighting-regime contamination of the key.
  // Off by default: a hit returns FALSE immediately and the default path is
  // unchanged.  getenv is read only on a hit (convergence-frequency, never in
  // the inner loop) and not cached, so the env var toggles reliably.  The
  // deterministic guard is test-ts-na-evcache.R.
  const bool cache_hit = evs_false_cache.count(cache_key) != 0;
  if (cache_hit && !std::getenv("TS_EV_AUDIT")) return false;

  save_topology(tree, snap);
  in_sub.assign(tree.n_node, 0);

  for (int clip_node = 0; clip_node < tree.n_node; ++clip_node) {
    if (clip_node == tree.n_tip) continue;            // display root
    int nx = tree.parent[clip_node];
    if (nx < 0 || nx == tree.n_tip) continue;         // root child (scan parity)
    int nz = tree.parent[nx];
    int nxi = nx - tree.n_tip;
    int ns = (tree.left[nxi] == clip_node) ? tree.right[nxi] : tree.left[nxi];

    // Mark the clipped subtree (invalid regraft targets), tracking set nodes
    // for O(subtree) reset.
    dfs.clear(); dfs.push_back(clip_node); marked.clear();
    while (!dfs.empty()) {
      int nd = dfs.back(); dfs.pop_back();
      in_sub[nd] = 1; marked.push_back(nd);
      if (nd >= tree.n_tip) {
        int ni = nd - tree.n_tip;
        dfs.push_back(tree.left[ni]);
        dfs.push_back(tree.right[ni]);
      }
    }

    collect_subtree_edges(tree, clip_node, sub_edges);   // fragment rerootings
    const int n_reroot = 1 + static_cast<int>(sub_edges.size());

    bool found = false;
    for (int ri = 0; ri < n_reroot && !found; ++ri) {
      int rp = (ri == 0) ? -1 : sub_edges[ri - 1].first;   // -1 => SPR (no reroot)
      int rc = (ri == 0) ? -1 : sub_edges[ri - 1].second;
      for (int below = 0; below < tree.n_node; ++below) {
        if (below == tree.n_tip || in_sub[below] || below == nx) continue;
        int above = tree.parent[below];
        if (below == ns) above = nz;          // merged edge after the clip
        if (above < 0 || above == nx) continue;
        if (ri == 0 && below == ns) continue; // identity (no reroot, original spot)

        if (!apply_tbr_move(tree, clip_node, rp, rc, above, below)) {
          restore_topology(tree, snap);
          continue;
        }
        tree.build_postorder();
        double s;
        if (na_incr_audit && ev_incr_ok) {
          // INCREMENTAL AUDIT: recompute this candidate via the 3-seed dirty
          // rescore from a clean BASE state, and compare to full_rescore.
          // (The legacy loop maintains only base topology, not base state, so we
          // reconstruct base state here.) full_rescore at the end overwrites the
          // dirty modifications, so no save/restore-undo is needed in audit mode;
          // we only clear the undo buffer to stop it growing across candidates.
          restore_topology(tree, snap);
          tree.build_postorder();
          full_rescore(tree, ds);                                 // base state
          apply_tbr_move(tree, clip_node, rp, rc, above, below);  // re-apply (known ok)
          tree.build_postorder();
          if (tree.prealloc_undo) tree.prealloc_undo->clear();
          else tree.clip_undo_stack.clear();
          int third = (clip_node >= tree.n_tip) ? clip_node : -1;
          fitch_na_dirty_downpass(tree, ds, nz, nx, third);
          fitch_na_dirty_uppass(tree, ds, nz, nx, third);
          double s_incr;
          if (ev_weighted) {
            incr_char_steps.resize(ds.n_patterns);
            fitch_na_pass3_score(tree, ds, &incr_char_steps);   // fused (= fast path)
            s_incr = compute_weighted_score(ds, incr_char_steps);
          } else {
            s_incr = static_cast<double>(fitch_na_pass3_score(tree, ds)) + ds.ew_offset;
          }
          if (tree.prealloc_undo) tree.prealloc_undo->clear();
          else tree.clip_undo_stack.clear();
          s = full_rescore(tree, ds);                             // authoritative
          if (std::fabs(s_incr - s) > na_incr_tol) {
            Rcpp::stop("TS_NA_INCR_AUDIT mismatch clip=%d ri=%d below=%d rp=%d rc=%d "
                       "incr=%.6f full=%.6f", clip_node, ri, below, rp, rc, s_incr, s);
          }
        } else if (na_incr && ev_incr_ok) {
          // FAST incremental path: 3-seed dirty rescore is the decision score.
          // State stays valid on accept; on reject it is undone below (restore_
          // prealloc_undo / restore_saved_states) before restore_topology.
          if (tree.prealloc_undo) tree.prealloc_undo->clear();
          else tree.clip_undo_stack.clear();
          int third = (clip_node >= tree.n_tip) ? clip_node : -1;
          fitch_na_dirty_downpass(tree, ds, nz, nx, third);
          fitch_na_dirty_uppass(tree, ds, nz, nx, third);
          if (ev_weighted) {
            // Fused: Pass3 fills char_steps in its own walk; skip extract.
            incr_char_steps.resize(ds.n_patterns);
            fitch_na_pass3_score(tree, ds, &incr_char_steps);
            s = compute_weighted_score(ds, incr_char_steps);
          } else {
            s = static_cast<double>(fitch_na_pass3_score(tree, ds)) + ds.ew_offset;
          }
        } else {
          s = full_rescore(tree, ds);
        }
        if (s < best_score - eps) {           // keep this improver applied
          if (cache_hit) {                    // TS_EV_AUDIT: cache lied
            Rcpp::stop("TS_EV_AUDIT: exact_verify cache returned FALSE (optimum) "
                       "for a topology with an improving neighbour (%.4f < %.4f) "
                       "- weighting-regime contamination in the cache key.",
                       s, best_score);
          }
          best_score = s;
          found = true;
          break;
        }
        // Reject: fast path must undo the dirty-rescore state before topology,
        // so the base state stays valid for the next candidate's dirty rescore.
        if (na_incr && ev_incr_ok) {
          if (tree.prealloc_undo) tree.restore_prealloc_undo();
          else tree.restore_saved_states();
        }
        restore_topology(tree, snap);
      }
    }
    for (int nd : marked) in_sub[nd] = 0;
    if (found) return true;
  }

  restore_topology(tree, snap);               // clip loop found no improver
  tree.build_postorder();
  best_score = full_rescore(tree, ds);

  // The clip loop skips root-child clips (the nx==n_tip guard above), so the ONE
  // unrooted edge n_tip sits on (cL-cR) is never enumerated — one blind edge with
  // a large neighbourhood, the residual that let poor NA starts converge with a
  // root-edge improver still present (tbr_oracle_na.R: was 2/20 on Zanol2014).
  // Enumerate it exactly here (apply + full_rescore, the same path IW uses at
  // convergence) before declaring an optimum, so the verdict — and the memoized
  // FALSE — means a TRUE unrooted-TBR optimum over all 2n-3 edges, not just the
  // 2n-4 non-root ones.
  if (try_root_edge_moves_rescore(tree, ds, best_score)) return true;

  evs_false_cache.insert(cache_key);          // memoize: true unrooted-TBR optimum
  return false;
}

// --- Edge length computation ---

// --- Subtree size computation ---

// Compute the number of tips in the subtree below each node.
// Result indexed by node id. Tips have size 1.
static void compute_subtree_sizes(const TreeState& tree,
                                  std::vector<int>& sizes) {
  sizes.assign(tree.n_node, 0);
  for (int i = 0; i < tree.n_tip; ++i) sizes[i] = 1;
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    sizes[node] = sizes[tree.left[ni]] + sizes[tree.right[ni]];
  }
}

// Add the clipped subtree's INTERNAL Fitch steps (per pattern, standard blocks)
// to char_steps.  The IW indirect scan needs base_iw = weighted(divided_tree +
// clip_internal): the EW divided_length (best_score + delta - nx_cost) keeps the
// clipped subtree's internal length, but extract_char_steps over the *clipped*
// postorder omits it (spr_clip removed those nodes), so base_iw systematically
// UNDER-counts by the clipped subtree's internal IW cost (small for tip/2-tip
// clips, large for the L872 clip-both large fragment) — corrupting cross-clip
// move ranking and leaving improving unrooted-TBR moves unreached.  X's internal
// per-character step count is root-invariant and its nodes retain valid
// local_cost after spr_clip (the clip only detaches clip_node; the incremental
// downpass from nz never descends into X), so a plain DFS sum is exact.
static void add_clip_internal_steps(const TreeState& tree, const DataSet& ds,
                                    int clip_node,
                                    std::vector<int>& char_steps) {
  std::vector<int> stack;
  stack.clear();
  stack.push_back(clip_node);
  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();
    if (node < tree.n_tip) continue;            // tips carry no internal step
    int ni = node - tree.n_tip;
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable || blk.active_mask == 0) continue;
      uint64_t mask =
          tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      while (mask) {
        int c = ctz64(mask);
        char_steps[blk.pattern_index[c]] += 1;
        mask &= mask - 1;
      }
    }
    stack.push_back(tree.left[ni]);
    stack.push_back(tree.right[ni]);
  }
}

// NOTE (NA clip-internal): the same base-score omission exists for inapplicable
// data — `fitch_na_pass3_score` over the clipped postorder drops the clipped
// subtree's internal NA homoplasy, so the NA indirect scan under-counts (scan
// vs full_rescore mispredicts up to ~200 on Vinther2008 via TS_IW_SCANCHK). A
// Pass-3 analog of add_clip_internal_steps was prototyped and DID collapse the
// error (NA+IW 6–14 → ±0.5), but the residual is MIXED-SIGN: NA's Pass-3 reads
// down2 (whole-tree uppass context), so the clipped subtree's internal count is
// attachment-DEPENDENT and the indirect decomposition cannot be made exact the
// way Fitch/IW can. The NA oracle (dev/benchmarks/tbr_oracle_na.R) confirms the
// direct scan stays incomplete even with the prototype, while the physical
// reroot path is 0-improving. So NA reaches true unrooted-TBR optima via the
// physical-reroot path (has_na dispatch in tbr_search), and the NA scan is left
// at production baseline. Making the rooted-NA scan more accurate (the prototype)
// is a SEPARATE change that needs its own search-quality validation.

// --- Precompute vroot cache for main edges ---

// For each main edge (A, D), compute vroot[s] = final_[A][s] | final_[D][s].
// vroot_cache layout: vroot_cache[edge_idx * total_words + s]
static void precompute_vroot_cache(
    const TreeState& tree,
    const std::vector<std::pair<int,int>>& main_edges,
    std::vector<uint64_t>& vroot_cache) {
  int tw = tree.total_words;
  size_t n_edges = main_edges.size();
  vroot_cache.resize(n_edges * tw);

  for (size_t ei = 0; ei < n_edges; ++ei) {
    int a = main_edges[ei].first;
    int d = main_edges[ei].second;
    size_t a_base = static_cast<size_t>(a) * tw;
    size_t d_base = static_cast<size_t>(d) * tw;
    size_t out_base = ei * tw;

    for (int s = 0; s < tw; ++s) {
      vroot_cache[out_base + s] = tree.final_[a_base + s]
                                | tree.final_[d_base + s];
    }
  }
}


// --- Clip ordering ---

// Apply the selected clip ordering strategy to clip_candidates.
// subtree_sizes must be up-to-date. n_tip is the number of tips.
static void order_clips(
    std::vector<int>& clips,
    const std::vector<int>& subtree_sizes,
    int n_tip,
    ClipOrder order,
    std::mt19937& rng)
{
  switch (order) {
    case ClipOrder::RANDOM:
      std::shuffle(clips.begin(), clips.end(), rng);
      break;

    case ClipOrder::INV_WEIGHT: {
      // Weighted random: w = 1/(1+s). Full Fisher-Yates with weighted draw.
      int n = static_cast<int>(clips.size());
      std::vector<double> w(n);
      for (int i = 0; i < n; ++i) {
        w[i] = 1.0 / (1.0 + subtree_sizes[clips[i]]);
      }
      for (int i = 0; i < n - 1; ++i) {
        double total = 0.0;
        for (int j = i; j < n; ++j) total += w[j];
        if (total <= 0.0) break;
        std::uniform_real_distribution<double> dist(0.0, total);
        double r = dist(rng);
        double cumul = 0.0;
        int pick = i;
        for (int j = i; j < n; ++j) {
          cumul += w[j];
          if (cumul >= r) { pick = j; break; }
        }
        std::swap(clips[i], clips[pick]);
        std::swap(w[i], w[pick]);
      }
      break;
    }

    case ClipOrder::TIPS_FIRST: {
      // Partition: tips first, then internal. Shuffle within each group.
      auto mid = std::partition(clips.begin(), clips.end(),
          [n_tip](int node) { return node < n_tip; });
      std::shuffle(clips.begin(), mid, rng);
      std::shuffle(mid, clips.end(), rng);
      break;
    }

    case ClipOrder::BUCKET: {
      // Three buckets: tips (s=1), small (2 <= s <= sqrt(n)), large (s > sqrt(n))
      int sqrt_n = static_cast<int>(std::sqrt(static_cast<double>(n_tip)));
      if (sqrt_n < 2) sqrt_n = 2;

      auto large_start = std::partition(clips.begin(), clips.end(),
          [&subtree_sizes, sqrt_n](int node) {
            return subtree_sizes[node] <= sqrt_n;
          });
      auto small_start = std::partition(clips.begin(), large_start,
          [n_tip](int node) { return node < n_tip; });

      // clips = [tips | small internal | large]
      std::shuffle(clips.begin(), small_start, rng);
      std::shuffle(small_start, large_start, rng);
      std::shuffle(large_start, clips.end(), rng);
      break;
    }

    case ClipOrder::ANTI_TIP: {
      // Non-tip clips (shuffled) first, tip clips (shuffled) last.
      // Hypothesis: tips are under-productive; deprioritise them.
      // Inverse of TIPS_FIRST.
      auto tip_start = std::partition(clips.begin(), clips.end(),
          [n_tip](int node) { return node >= n_tip; }); // non-tips first
      std::shuffle(clips.begin(), tip_start, rng);
      std::shuffle(tip_start, clips.end(), rng);
      break;
    }

    case ClipOrder::LARGE_FIRST: {
      // Large (>√n) clips first, then small (2..√n), then tips; random within.
      // Hypothesis: large clips are enriched relative to their clip-fraction.
      int sqrt_n = static_cast<int>(std::sqrt(static_cast<double>(n_tip)));
      if (sqrt_n < 2) sqrt_n = 2;

      // Partition: [large | rest]
      auto rest_start = std::partition(clips.begin(), clips.end(),
          [&subtree_sizes, sqrt_n](int node) {
            return subtree_sizes[node] > sqrt_n;
          });
      // Partition rest: [large | small-internal | tips]
      auto tip_start = std::partition(rest_start, clips.end(),
          [n_tip](int node) { return node >= n_tip; }); // non-tip non-large = small

      // clips = [large | small-internal | tips]
      std::shuffle(clips.begin(), rest_start, rng);
      std::shuffle(rest_start, tip_start, rng);
      std::shuffle(tip_start, clips.end(), rng);
      break;
    }
  }
}

// --- Lever #7: monomorphized plain-EW SPR candidate scan ---
//
// The general SPR loop in tbr_search re-tests, for EVERY candidate, checks that
// are constant for the whole search: the criterion flavour (has_na / use_iw /
// all_weight_one) and the search-mode gates (sector_mask / constrained /
// collapsed / b2_ceiling). This template lifts that dispatch out of the hot
// loop: it is instantiated ONCE per weight-class at dispatch (a runtime switch
// in tbr_search selects the instantiation), so the dead-in-plain-EW branches
// compile away and the compiler picks a specialised kernel.
//
// The template is called ONLY in the plain-EW regime (no NA, no IW, no sector /
// constraint / collapsed / b2), so keeping only the identity skip + scorer +
// strict-< accept is trajectory-identical to the general loop: the removed
// branches never fire there, so the candidate sequence, every accept, and the
// n_evaluated count are BYTE-IDENTICAL. Only the per-candidate instruction
// count differs.
//
// UseFlat is the payoff: the weight-blind flat kernel (fitch_indirect_cached_
// flat) drops the CharBlock deref + weight-multiply + active_mask==0 skip. It
// is UNSAFE as a global toggle (weight-blind), but SAFE as a template
// instantiation because dispatch selects UseFlat=true ONLY for all_weight_one
// data, where it computes exactly the same extra_steps as the general
// fitch_indirect_length_cached (verified byte-identical, dev/profiling/
// s7-fastpath-sizing.md). UseFlat=false (weighted / ratchet upweight) keeps the
// general scorer, so that path is byte-identical too.
//
// Wrong is a positive-control ONLY: the <..,true> instantiations are emitted
// solely under -DTS_EW_MONO_WRONG (see the dispatch), so a production binary
// has NO code path — and no env var — that can corrupt the scorer. Under that
// build flag it corrupts the output so a live specialised path provably changes
// the result, distinguishing "path fired + output used" from a no-op that
// passes the gate falsely (the false-0% trap documented in s7-fastpath-sizing).
template<bool UseFlat, bool Wrong>
static inline void spr_scan_plain_ew(
    const std::vector<std::pair<int,int>>& main_edges,
    int nz, int ns,
    const uint64_t* clip_prelim,
    const std::vector<uint64_t>& edge_set_buf,
    int tw, const DataSet& ds,
    double divided_length,
    double& best_candidate,
    int& best_above, int& best_below,
    int& best_reroot_parent, int& best_reroot_child,
    int& n_evaluated, int& cutoff) {
  for (auto& [above, below] : main_edges) {
    if (above == nz && below == ns) continue;
    const uint64_t* vroot = &edge_set_buf[static_cast<size_t>(below) * tw];
    int extra;
    if constexpr (UseFlat) {
      // Weight-blind flat kernel — safe ONLY because dispatch selects UseFlat
      // for all_weight_one data, where it returns exactly the same extra_steps
      // as the general scorer below (verified byte-identical: the accept
      // decision is identical because a bailed candidate, extra >= cutoff, is
      // always rejected, and a non-bailed candidate returns the true total).
      extra = fitch_indirect_cached_flat(clip_prelim, vroot, ds, cutoff);
    } else {
      extra = fitch_indirect_length_cached(clip_prelim, vroot, ds, cutoff);
    }
    // Positive control only (never in prod): a DATA-DEPENDENT, always->=1
    // corruption. Unlike a uniform +1 (which preserves intra-SPR ordering and
    // so only flips SPR-vs-reroot ties), this perturbs each candidate by a
    // different amount, so it changes the SPR argmin whenever an SPR move can
    // win — a faithful stand-in for "the flat kernel returns wrong values".
    if constexpr (Wrong) extra += static_cast<int>(below & 7) + 1;
    double candidate = divided_length + extra;
    ++n_evaluated;
    if (candidate < best_candidate) {
      best_candidate = candidate;
      best_above = above;
      best_below = below;
      best_reroot_parent = -1;
      best_reroot_child = -1;
      cutoff = static_cast<int>(best_candidate - divided_length + 1);
    }
  }
}

// Monomorphized plain-IW SPR candidate scan (lever #7, IW/XPIWE analogue of
// spr_scan_plain_ew). Same dead-branch strip (no NA / no EW / no sector /
// constraint / collapsed / b2), dispatched once at the scan site for a pure-IW
// (has_na=false, use_iw=true) search. There is NO flat-kernel bonus here: the
// implied-weights scorer must apply the per-pattern concavity weights (iw_delta)
// on every candidate, so — unlike EW — there is no weight-blind variant to swap
// in. The win is only the ~5%-of-SPR dead-branch strip, on a heavier
// double-accumulator scorer, so it is smaller than the EW-flat win.
//
// Byte-identical to the general loop's IW branch: the scorer
// (indirect_iw_length_cached, bounding on best_candidate) and the accept block
// are reproduced verbatim — including the `cutoff` update, which is dead on the
// IW path (the IW reroot bounds on best_candidate, never cutoff) but is kept
// line-for-line so byte-identity does not rest on a "cutoff is IW-dead"
// argument. Wrong is the compile-gated positive control (see spr_scan_plain_ew).
template<bool Wrong>
static inline void spr_scan_plain_iw(
    const std::vector<std::pair<int,int>>& main_edges,
    int nz, int ns,
    const uint64_t* clip_prelim,
    const std::vector<uint64_t>& edge_set_buf,
    int tw, const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double divided_length,
    double& best_candidate,
    int& best_above, int& best_below,
    int& best_reroot_parent, int& best_reroot_child,
    int& n_evaluated, int& cutoff) {
  for (auto& [above, below] : main_edges) {
    if (above == nz && below == ns) continue;
    double candidate = indirect_iw_length_cached(
        clip_prelim, &edge_set_buf[static_cast<size_t>(below) * tw],
        ds, base_iw, iw_delta, best_candidate);
    if constexpr (Wrong) candidate += static_cast<double>((below & 7) + 1);
    ++n_evaluated;
    if (candidate < best_candidate) {
      best_candidate = candidate;
      best_above = above;
      best_below = below;
      best_reroot_parent = -1;
      best_reroot_child = -1;
      cutoff = static_cast<int>(best_candidate - divided_length + 1);
    }
  }
}

// --- Main TBR search ---

TBRResult tbr_search(TreeState& tree, const DataSet& ds,
                     const TBRParams& params,
                     ConstraintData* cd,
                     const std::vector<bool>* sector_mask,
                     TreePool* collect_pool,
                     std::function<bool()> check_timeout) {
  double best_score = full_rescore(tree, ds);
  // Tracks whether `best_score` is the authoritative score of the current
  // (tree, state arrays). Each accepted move and each state_snap.restore
  // re-establishes the invariant; apply_tbr_move temporarily breaks it
  // until the following full_rescore + best_score update completes.
  // The trailing full_rescore at function exit is gated on this flag,
  // skipping a redundant O(n_node x n_char) pass when states are coherent.
  bool score_fresh = true;

  // No informative characters: all trees have the same score.
  if (ds.total_words == 0) {
    return {best_score, 0, 0, 0, true};
  }

  // Initialize constraint mapping if active
  bool constrained = cd && cd->active;
  if (constrained) {
    update_constraint(tree, *cd);
  }
  int n_accepted = 0;
  int n_evaluated = 0;
  int n_zero_skipped = 0;
  int hits = 1;
  const bool use_iw = std::isfinite(ds.concavity);
  // Floating-point tolerance for score equality
  const double eps = use_iw ? 1e-10 : 0.0;

  // DIAGNOSTIC (env TS_IW_TIMING): isolate IW-specific cost under REAL bail
  // (production cutoffs), EW vs IW on the same data. Read ONCE (no per-clip
  // getenv). t_clip_ns = the bail-INDEPENDENT per-clip IW precompute
  // (extract_char_steps + compute_iw base + iw_delta; no EW analog); t_scan_ns
  // = the bail-DEPENDENT per-candidate scan (EW x4-flat/popcount vs IW
  // scalar-gather), normalised by n_evaluated. No-op unless the var is set.
  const bool iw_timing = std::getenv("TS_IW_TIMING") != nullptr;
  long long t_clip_ns = 0, t_spr_ns = 0, t_rer_ns = 0;
  long long n_clips_t = 0, n_spr_t = 0, n_rer_t = 0;
  // Pure-IW / NA-IW reroot 4-wide batching (T-245 ILP ported to IW). On by
  // default (fires on IW+XPIWE, byte-identical to scalar per the dirty-rescore
  // guard + the 3725746a validation); kill-switch TS_IW_NOX4 reverts to the
  // scalar reroot path (for A/B). Read once (no per-clip getenv, per the
  // getenv-cost lesson).
  const bool iw_x4 = std::getenv("TS_IW_NOX4") == nullptr;
  // Pure-IW extract_char_steps dirty-region: derive per-clip divided_steps
  // incrementally (full_char_steps + cs_delta - nx) instead of the O(n_node)
  // walk. On by default; kill-switch TS_IW_NODIRTY reverts to extract+add.
  const bool iw_dirty = std::getenv("TS_IW_NODIRTY") == nullptr;
  // DIAGNOSTIC: per-clip compare dirty divided_steps vs extract+add baseline.
  const bool iw_dirtychk = std::getenv("TS_IW_DIRTYCHK") != nullptr;

  // Check if any block has inapplicable characters (for state snapshot)
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  // T-306: the SPR accept-path computes `actual` as an incremental delta on
  // top of best_score — a Fitch-only EW delta (best_score + delta) or an
  // IW/profile rescore from per-pattern step counts.  For HSJ/XFORM scoring,
  // score_tree() additionally adds a topology-dependent hierarchy-DP (HSJ) or
  // Sankoff (XFORM) term that neither delta captures, so best_score would
  // drift from the authoritative score and corrupt accept/reject decisions.
  // Restrict the incremental fast path to scoring modes whose total equals the
  // Fitch/IW result; HSJ and XFORM fall back to full_rescore (the same scoring
  // -mode classification used by the T-275/T-303 guards).
  const bool incremental_ok =
      ds.scoring_mode == ScoringMode::EW ||
      ds.scoring_mode == ScoringMode::IW ||
      ds.scoring_mode == ScoringMode::XPIWE ||
      ds.scoring_mode == ScoringMode::PROFILE;

  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  // Tabu list: prevent cycling during plateau exploration
  TabuList tabu(params.tabu_size);
  if (tabu.active()) {
    tabu.insert(hash_tree(tree));
  }

  // Candidate clip nodes: all non-root nodes
  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    clip_candidates.push_back(node);
  }

  // Collapsed flags: edges that provably cannot yield an improvement
  // (clip skipping + regraft merging).  Disabled during MPT enumeration
  // (collect_pool) where equal-score topologies are collected.
  //
  // TS_COLLAPSE_AGGRESSIVE (default OFF, prototype): use the TNT `collapse 3`
  // criterion (collapse min-length-0 INTERNAL branches) instead of the exact
  // score-identity criterion.  This is a HEURISTIC neighbourhood reduction
  // (Goloboff asymmetric reachability) — scoring stays exact, but some improving
  // moves may be skipped.  Aimed at large/molecular/congruent datasets where
  // zero-length branches are common (the morphological roster has ~0 of them).
  // See dev/plans + b2_collapsed_kernel_validate.R.  Scoped to tbr_search's
  // neighbourhood reduction only; pool dedup keeps the exact flags.
  const bool collapse_aggr = std::getenv("TS_COLLAPSE_AGGRESSIVE") != nullptr;
  std::vector<uint8_t> collapsed;
  if (!collect_pool) {
    if (collapse_aggr) compute_collapsed_flags_aggressive(tree, ds, collapsed);
    else               compute_collapsed_flags(tree, ds, collapsed);
  }
  // Reusable scratch for MPT-enumeration pool dedup (collect_pool set).  The
  // enumeration deduplicates on COLLAPSED topology to match the main-loop pool
  // (add_collapsed in ts_driven.cpp / ts_parallel.cpp); a plain add() would
  // over-collect zero-length-resolution variants of one collapsed topology,
  // inflating n_topologies and the returned multiPhylo on soft-polytomy
  // matrices (see memory mpt-enumeration-dedup-asymmetry).  Conservative flags
  // (not the aggressive criterion) to match the main loop exactly; on data with
  // no collapsible branches this is a byte-for-byte no-op (all-zero flags ->
  // compute_collapsed_splits == compute_splits).
  //
  // TS_ENUM_RESOLVED (default OFF): kill-switch reverting enumeration to the
  // pre-fix plain add() (resolved-topology dedup).  The `collect_pool &&`
  // short-circuit keeps the getenv off the main-search hot path entirely
  // (collect_pool is nullptr everywhere except MPT enumeration — verified at
  // every tbr_search call site).
  std::vector<uint8_t> collect_collapsed;
  const bool enum_resolved =
      collect_pool && std::getenv("TS_ENUM_RESOLVED") != nullptr;
  // Hoisted per-call-invariant flags. collapsed's empty-ness never changes
  // within a tbr_search call (the reroot-loop recompute is guarded on
  // non-empty and only refreshes contents); the diagnostic env vars are
  // constant for the process. Avoids a per-candidate vector::empty() and a
  // per-clip / per-accept / per-reroot getenv scan. Byte-identical (gates the
  // same checks). getenv is ~2.4 us/call on Windows/ucrt (linear env-block
  // scan); its cost hides in VTune's ucrtbase self-time, so a per-clip read is
  // a real, profiler-invisible cost (see findings.md T-P5n).
  const bool use_collapsed = !collapsed.empty();
  // Lever #7 dispatch gate: the collapsed vector is ALWAYS sized n_node
  // (compute_collapsed_flags does assign(n_node, 0)), so use_collapsed is true
  // even for morphological EW where every flag is 0 and the per-candidate
  // collapsed[below] check never skips. The monomorphized plain-EW scan can
  // strip that check only when every flag is CURRENTLY 0. Gating on
  // collapsed_all_zero (not !use_collapsed) is what avoids the false-BASE-vs-
  // BASE trap that a first A/B hit (dev/profiling/s7-fastpath-sizing.md).
  //
  // collapsed IS mutated across clips: compute_collapsed_flags is re-run after
  // every accepted move (post-accept, and after each root-edge re-descent),
  // and a new tree can introduce zero-length branches => set flags. So this
  // bool must be REFRESHED at every recompute, not cached once — else the gate
  // goes stale, the template wrongly skips the collapsed check, and it
  // evaluates edges the general loop skips (a byte-identity break the
  // MaximizeParsimony gate caught on Zhu2013). Cost: an O(n_node) scan per
  // accept — negligible beside the candidate scan.
  bool collapsed_all_zero = true;
  auto refresh_collapsed_all_zero = [&collapsed, &collapsed_all_zero]() {
    collapsed_all_zero = true;
    for (size_t ci = 0; ci < collapsed.size(); ++ci)
      if (collapsed[ci]) { collapsed_all_zero = false; break; }
  };
  refresh_collapsed_all_zero();
  const bool revert_check = std::getenv("TS_REVERT_CHECK") != nullptr;
  const bool iw_scanchk = std::getenv("TS_IW_SCANCHK") != nullptr;
  // TS_PHYS_REROOT selects the legacy physical-reroot reference path; it is read
  // once per outer reroot-loop iteration below (>=1/call), so hoist it too.
  const bool phys_reroot = std::getenv("TS_PHYS_REROOT") != nullptr;
  // B2 ceiling probe (measurement only, env-gated => production byte-identical):
  // count SPR-regraft candidates that pass collapse Condition 1 (zero parent
  // cost) but are NOT conservatively collapsed (Condition 3 fails) — the set an
  // aggressive TNT-style collapse-on-equal-length could skip but the exact merge
  // cannot.  Work-weighted against evaluated candidates = upper bound on B2's
  // realizable regraft saving.  See dev/benchmarks/b2_collapsed_density.R.
  const bool b2_ceiling = std::getenv("TS_B2_CEILING") != nullptr;
  long long n_b2_eval = 0, n_b2_aggr = 0;

  // Lever #6 / L3b (dev/plans/2026-07-14-lever6-incremental-edgeset-land.md):
  // maintain the directional insertion edge set INCREMENTALLY per clip
  // (patch-from-full "Scheme 1") instead of the from-scratch O(n_node)
  // compute_insertion_edge_sets recompute.  Env-gated OFF by default; setting
  // TS_L3B_INCREMENTAL turns on the incremental path (config C).  With the flag
  // OFF the search is byte-identical to production (config A) — so C-vs-A is a
  // pure wall-clock comparison on an identical trajectory (the ship gate).  The
  // per-clip oracle (equality vs the from-scratch result) fires when
  // TS_L3B_ORACLE is set OR under an asserts-on (NDEBUG-off) build.
  const bool l3b_incremental = std::getenv("TS_L3B_INCREMENTAL") != nullptr;
#ifdef NDEBUG
  const bool l3b_oracle = std::getenv("TS_L3B_ORACLE") != nullptr;
#else
  const bool l3b_oracle = true;
#endif
  long long l3b_patch_clips = 0, l3b_changed_sum = 0, l3b_edges_sum = 0;

  // Lever #7 (dev/plans/2026-07-14-lever7-scorer-monomorphization.md): the
  // monomorphized plain-EW SPR scan is the DEFAULT path; TS_EW_MONO_OFF reverts
  // to the general per-candidate-dispatch loop (baseline for the byte-identity
  // gate + A/B, and a production kill-switch). Exact/byte-identical, so it is
  // safe as the default.
  const bool ew_mono = std::getenv("TS_EW_MONO_OFF") == nullptr;
#ifdef TS_EW_MONO_WRONG
  // Positive-control build only (compile-time gated so there is NO env-var
  // path to a corrupted scorer in a production binary — the same discipline as
  // the removed TS_EW_MONO_ASSERT probe). Rebuild with -DTS_EW_MONO_WRONG in
  // PKG_CPPFLAGS, then TS_EW_MONO_WRONG=1 selects the corrupted instantiation
  // (must DIVERGE — proves the specialised path's output drives the search).
  const bool ew_mono_wrong = std::getenv("TS_EW_MONO_WRONG") != nullptr;
#endif
  // Positive-control fire counters, split by weight class: proving A fast path
  // fired is NOT enough (the UseFlat=false path calls the same scorer as
  // baseline, so it passes the gate trivially). To know the FLAT kernel — the
  // actual #7 payoff — executed, spr_mono_flat_fired must be nonzero on an
  // all_weight_one dataset. Reported in the TS_IW_TIMING line.
  long long spr_mono_flat_fired = 0, spr_mono_cached_fired = 0, spr_mono_iw_fired = 0;

  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;
  std::vector<int> kept_ei;  // per-clip non-skipped main-edge indices (reroot)

  // Temporary buffers
  std::vector<uint64_t> from_above(
      static_cast<size_t>(tree.n_node) * tree.total_words, 0);
  std::vector<uint64_t> virtual_prelim(tree.total_words);

  // IW buffers (allocated once, reused per clip)
  std::vector<int> divided_steps;
  std::vector<double> iw_delta;
  if (use_iw) {
    divided_steps.resize(ds.n_patterns, 0);
    iw_delta.resize(ds.n_patterns, 0.0);
  }
  // IW dirty-region buffers: F = full-tree char_steps (refreshed when the tree
  // changes, f_dirty), cs_delta = per-clip path delta (filled by the downpass),
  // nx_cs = clip-parent's per-pattern contribution. Pure-IW only.
  // The char_steps dirty-region + x4 reroot opts apply to BOTH implied-weights
  // modes: plain IW and extended XPIWE. Both score via compute_iw /
  // precompute_iw_delta from the SAME per-pattern char_steps (extract_char_steps
  // is weighting-agnostic; min_steps identical), differing ONLY in the per-pattern
  // eff_k[p]/phi[p] arrays — which those downstream functions already consume. So
  // the dirty-region derivation (divided_steps = F + cs_delta - nx) and the x4
  // gather are byte-identical across IW and XPIWE. PROFILE is excluded: it scores
  // via compute_profile (info_amounts table + precomputed_steps offset), a
  // genuinely different char_steps convention. NOTE: production MaximizeParsimony
  // defaults to XPIWE (xpiwe=TRUE + obs_count; ts_data.cpp), so gating on plain
  // IW alone left these opts DEAD on the production path for every dataset.
  const bool iw_family = (ds.scoring_mode == ScoringMode::IW ||
                          ds.scoring_mode == ScoringMode::XPIWE);
  const bool iw_dirty_active = iw_family && !has_na && iw_dirty;
  std::vector<int> full_char_steps, cs_delta, nx_cs;
  if (iw_dirty_active) {
    full_char_steps.resize(ds.n_patterns, 0);
    cs_delta.resize(ds.n_patterns, 0);
    nx_cs.resize(ds.n_patterns, 0);
  }

  // Subtree sizes for smaller-subtree filtering
  std::vector<int> subtree_sizes(tree.n_node, 0);

  // Vroot cache for TBR rerooting (precomputed per clip)
  std::vector<uint64_t> vroot_cache;

  // State snapshot for rejection without full_rescore (optimization #3)
  StateSnapshot state_snap;
  state_snap.allocate(tree, has_na);

  // Pre-allocated undo stack: eliminates ~50 heap allocs per clip.
  // save_node_state() writes to flat buffers instead of allocating vectors.
  TreeState::PreallocUndo fast_undo;
  // Capacity must cover downpass + uppass + tips: up to 3 * n_node saves per clip
  fast_undo.init(3 * tree.n_node, tree.total_words, tree.n_blocks, has_na);
  tree.prealloc_undo = &fast_undo;

  // Pre-allocated work buffer for build_postorder_prealloc
  std::vector<int> work_stack;
  work_stack.reserve(tree.n_node * 2);

  // Pre-allocate postorder save buffer for clip/unclip cycle
  std::vector<int> saved_postorder;
  saved_postorder.reserve(tree.postorder.size());

  // Pre-allocated clip_actives buffer (NA indirect evaluation)
  std::vector<uint64_t> clip_actives_buf;
  if (has_na) clip_actives_buf.resize(tree.total_words);

  // Pre-allocated below_actives cache (NA TBR rerooting)
  std::vector<uint64_t> below_actives_cache;

  // Per-clip rerooting dedup table (Tier 2). Declared once here, reset() per
  // clip — see VpHashSet above for why this is a plain local, not thread_local.
  VpHashSet seen_vp_hashes;

  // Use flat FlatBlock variants for indirect scoring when weight==1 and no
  // upweight_mask is active.  This is true for normal EW search; false during
  // ratchet phases that apply upweighting.  Checked once per tbr_search call.
  bool use_flat = ds.all_weight_one;
  if (use_flat) {
    for (int b_chk = 0; b_chk < ds.n_blocks; ++b_chk)
      if (ds.blocks[b_chk].upweight_mask) { use_flat = false; break; }
  }
  // Cache total_words as a local int to keep inner-loop expressions concise.
  const int tw = tree.total_words;

  // Pure-EW path scores candidate edges with the EXACT directional insertion
  // edge set (combine of the two directional Fitch messages) instead of the
  // union-of-finals (final[A] | final[D]) approximation, which overcounts and
  // can hide improving moves.  NA (three-pass) and implied-weights keep their
  // existing scorers; edge_set_buf is reused across clips.
  const bool ew_directional = !has_na && !use_iw;
  // Exact directional insertion edge sets apply to BOTH EW and IW (the join is
  // a Fitch operation on state sets; IW just weights the resulting per-char
  // steps).  Using them replaces the union-of-finals approximation
  // (final[a]|final[d]) that mis-counts and hides improving moves -- the same
  // bug the EW directional fix cured, now extended to IW.  NA keeps its own
  // 3-pass scorers.
  const bool use_directional = !has_na;
  std::vector<uint64_t> edge_set_buf;
  // Caller-owned scratch for compute_insertion_edge_sets, reused across clips
  // (size-ensured, non-zeroing) so the up-message buffer and preorder list are
  // not reallocated/zeroed every clip.
  std::vector<uint64_t> edge_set_up;
  std::vector<int> edge_set_pre;
  // L3b incremental path (lever #6): pristine intact-tree base recomputed once
  // per pass; the working edge_set_buf/edge_set_up are patched per clip and
  // restored (memcpy base -> buf over `l3b_changed`) between clips.  Oracle
  // scratch holds the from-scratch reference for the per-clip equality assert.
  std::vector<uint64_t> edge_set_base;
  std::vector<uint64_t> up_base;
  std::vector<int> l3b_changed;
  std::vector<int> l3b_worklist;
  std::vector<uint64_t> l3b_oracle_es, l3b_oracle_up;
  std::vector<int> l3b_oracle_pre;
  // Base-incremental-update (dev/profiling/l3b-land.md follow-on): after an
  // accepted SPR move, patch the base in place instead of the per-pass full
  // recompute.  base_current == "base + working buffers already reflect the
  // current tree" — set true after a successful incremental update, false on
  // first entry / after a reroot / after a TBR-reroot (non-SPR) accept, forcing
  // a full recompute at the next pass top.  Kill switch TS_L3B_NOBASEINCR.
  std::vector<uint64_t> l3b_base_tmp;
  bool l3b_base_current = false;
  const bool l3b_baseincr = std::getenv("TS_L3B_NOBASEINCR") == nullptr;

  TopoSnapshot snap;
  bool keep_going = true;
  bool need_shuffle = true;  // optimization #6: defer reshuffle
  // Poll timeout every n_tip clips to avoid overhead on small trees
  // while ensuring responsiveness on large ones.
  const int timeout_interval = std::max(tree.n_tip, 50);
  int clips_since_timeout_check = 0;
  bool timed_out = false;

  // Per-pass diagnostic counters
  int pass_index = 0;
  int pass_clips_tried = 0;
  int pass_candidates_evaluated = 0;
  int accepted_clip_size = 0;
  std::vector<TBRPassRecord> diag_records;

  // ===== Outer unrooted reroot loop (params.unrooted) =====
  // After the inner loop converges at one rooting, re-root at the next tip and
  // re-descend; stop when a full tip-sweep yields no strict improvement (=> a
  // true unrooted-TBR optimum).  Gated to the plain search — sector/constraint/
  // tabu/pool keep node-id-keyed state that re-rooting would invalidate.
  const bool do_reroot = params.unrooted && sector_mask == nullptr
      && cd == nullptr && params.tabu_size == 0 && collect_pool == nullptr
      && tree.n_tip >= 4;
  // `do_reroot` == "the unrooted-completeness mechanism is live here" (plain
  // search only).  It gates the three enumeration relaxations below AND the
  // root-edge branch; the per-scorer split (EW additive vs IW/NA apply+rescore)
  // happens inside try_root_edge_moves.  Sector/ratchet sub-searches keep the
  // default (smaller-side clip, nz/ns skip) so they pay no extra cost.
  int reroot_tip = 0;
  int reroot_clean = 0;            // consecutive reroots with no strict gain
  double reroot_prev = HUGE_VAL;
  bool first_descent = true;

  // L3b incremental edge-set maintenance is valid in the plain directional
  // regime: EW/IW (use_directional), no NA, no sector / constraint / tabu /
  // pool, no b2 probe.  The unrooted reroot loop (do_reroot) IS supported: it
  // still restores the tree to the pass-start base after every clip (so the
  // per-pass base is valid for all clips), re-roots only BETWEEN passes (base is
  // recomputed at each pass top), and its root-edge moves use their own buffers
  // — never edge_set_buf.  do_reroot only relaxes the smaller-side clip filter,
  // so it exercises larger clips with larger footprints (a wall effect the ship
  // test measures, not a correctness one).  Every other configuration falls back
  // to the from-scratch recompute.  Loop-invariant, so evaluated once.
  const bool l3b_active = l3b_incremental && use_directional
      && sector_mask == nullptr && cd == nullptr && params.tabu_size == 0
      && collect_pool == nullptr && !b2_ceiling
      && tree.n_tip >= 4;

  for (;;) {
   keep_going = true;
   // L3b: the tree is (re)established here — first entry, or re-rooted at the
   // end of the previous outer iteration (do_reroot).  Either invalidates the
   // incremental base, so force a full recompute at the first inner pass.
   l3b_base_current = false;
  while (keep_going && !timed_out) {
    keep_going = false;

    // Optimization #7: save state snapshot once per pass, not per candidate.
    // After a rejected move, state_snap.restore() returns the tree to exactly
    // the state saved here. The per-candidate save was redundant: consecutive
    // rejections all restore to the same state. Re-saving only happens when
    // the while loop restarts after an accepted move.
    save_topology(tree, snap);
    state_snap.save(tree);

    // IW dirty-region: refresh the full-tree char_steps cache F at the pass
    // start. clip-undo restores every clip in this pass to exactly this tree
    // (verified by TS_REVERT_CHECK), so F is valid for all clips in the pass;
    // per-clip derivation is F + cs_delta - nx_cs.  Refreshing here (not on
    // accept) aligns F with the snapshot the clips actually restore to.
    if (iw_dirty_active) extract_char_steps(tree, ds, full_char_steps);

    // Reset per-pass diagnostic counters
    pass_clips_tried = 0;
    pass_candidates_evaluated = 0;
    accepted_clip_size = 0;

    // Recompute subtree sizes (needed for smaller-subtree filtering
    // and for clip ordering strategies)
    compute_subtree_sizes(tree, subtree_sizes);

    // L3b: establish the pristine intact-tree base edge_set/up for this pass
    // (the tree is a full binary tree here, restored to this state after every
    // clip), then seed the working buffers from it.  Each clip patches the
    // working buffers and restores them (base -> buf over l3b_changed) before
    // the next clip, so they equal the base at every clip start.  When
    // l3b_base_current the base + buffers were already brought current by the
    // accept-time incremental update (base-incr follow-on), so the O(n_node)
    // full recompute + memcpy is skipped — this is what removes the per-pass
    // fixed cost that drowned the clip-patch win in accept-dense passes.
    if (l3b_active && !l3b_base_current) {
      compute_insertion_edge_sets(tree, ds, edge_set_base, up_base,
                                  edge_set_pre);
      const size_t nbuf = static_cast<size_t>(tree.n_node) * tw;
      if (edge_set_buf.size() < nbuf) edge_set_buf.resize(nbuf);
      if (edge_set_up.size() < nbuf) edge_set_up.resize(nbuf);
      std::memcpy(edge_set_buf.data(), edge_set_base.data(),
                  nbuf * sizeof(uint64_t));
      std::memcpy(edge_set_up.data(), up_base.data(), nbuf * sizeof(uint64_t));
      l3b_base_current = true;
    }

    // Optimization #6: only reorder when the previous pass found no
    // improvement. After an accepted move, retry with the same ordering
    // (the topology changed, so previously-failing clips may now succeed).
    if (need_shuffle) {
      order_clips(clip_candidates, subtree_sizes, tree.n_tip,
                  params.clip_order, rng);
    }
    need_shuffle = true;  // default: reorder next time (unless we accept)

    for (int clip_node : clip_candidates) {
      if (tree.parent[clip_node] == tree.n_tip) continue;

      // CSS: skip clips outside the sector
      if (sector_mask && !(*sector_mask)[clip_node]) continue;

      // Optimization #2: skip clips of the larger subtree (only clip the side
      // with fewer tips).  RELAXED under the opt-in complete-TBR path: the
      // Step-1 oracle factorial showed smaller-side-only is NOT complete —
      // clipping the larger side recovers a distinct ~5/100 of missed moves
      // (separate from the nz/ns and root-edge gaps; all three are needed for
      // 0/N).  Default keeps the filter (a perf win) and stays byte-identical.
      int clip_size = subtree_sizes[clip_node];
      if (!do_reroot && clip_size > tree.n_tip / 2) continue;

      // Skip collapsed edges: zero-length edge where clipping provably
      // cannot improve the score. Works for EW, IW, Profile, and NA.
      // Disabled during MPT enumeration (collapsed is empty).
      if (use_collapsed && collapsed[clip_node]) {
        ++n_zero_skipped;
        continue;
      }

      ++pass_clips_tried;
      int clip_evals_before = n_evaluated;

      // --- Phase 1: Clip + indirect evaluation ---

      // Save clip subtree's actives before clipping (needed for NA indirect)
      const uint64_t* clip_actives = nullptr;
      if (has_na) {
        size_t clip_sa_base =
            static_cast<size_t>(clip_node) * tree.total_words;
        std::memcpy(clip_actives_buf.data(),
                    &tree.subtree_actives[clip_sa_base],
                    tree.total_words * sizeof(uint64_t));
        clip_actives = clip_actives_buf.data();
      }

      // Save postorder before clip (restored after unclip instead of rebuild)
      saved_postorder.assign(tree.postorder.begin(), tree.postorder.end());

      fast_undo.clear();
      tree.spr_clip(clip_node);
      tree.build_postorder_prealloc(work_stack);

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;
      int nx = tree.clip_state.clip_parent;

      double divided_length;
      if (has_na) {
        // NA-aware incremental three-pass: correct prelim, final_,
        // subtree_actives, down2, and exact divided-tree score.
        fitch_na_incremental_downpass(tree, ds, nz);
        fitch_na_incremental_uppass(tree, ds, nz);
        divided_length = static_cast<double>(fitch_na_pass3_score(tree, ds));
      } else {
        int delta = fitch_incremental_downpass(tree, ds, nz,
            iw_dirty_active ? &cs_delta : nullptr);
        fitch_incremental_uppass(tree, ds, nz);

        if (iw_dirty_active) std::fill(nx_cs.begin(), nx_cs.end(), 0);
        int nx_cost = 0;
        for (int b = 0; b < ds.n_blocks; ++b) {
          uint64_t lc = tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b];
          int nu = popcount64(lc);
          if (ds.blocks[b].upweight_mask) nu += popcount64(lc & ds.blocks[b].upweight_mask);
          nx_cost += ds.blocks[b].weight * nu;
          // nx is spliced out by the divided tree; capture its per-pattern
          // contribution (raw local_cost bits). MUST mirror extract_char_steps,
          // which skips fully-inactive blocks (active_mask == 0): under a ratchet
          // perturbation that zeroes a block, F (the pass-top extract_char_steps
          // cache) omits that block's patterns, so subtracting nx_cs for them
          // drove divided_steps negative (DIRTYCHK: dirty=-1 ref=0 on
          // Dikow2009/Vinther2008 at ratchetCycles>=3). The candidate scan masks
          // needs_step by active_mask so the score was unaffected, but the
          // intermediate per-pattern step counts must stay consistent with F.
          // EW nx_cost above is intentionally left counting all blocks (its own
          // best_score/delta length convention is internally consistent).
          if (iw_dirty_active && ds.blocks[b].active_mask != 0) {
            const int* pidx = ds.blocks[b].pattern_index;
            uint64_t m = lc;
            while (m) { int c = ctz64(m); nx_cs[pidx[c]]++; m &= m - 1; }
          }
        }
        divided_length = best_score + delta - nx_cost;
      }

      // For weighted scoring (IW or profile): precompute base score and deltas
      double base_iw = 0.0;
      std::chrono::steady_clock::time_point _t_clip;
      if (iw_timing) _t_clip = std::chrono::steady_clock::now();
      if (use_iw) {
        if (iw_dirty_active) {
          // Dirty-region: divided_steps = (divided tree + clip internal) char
          // counts = F + cs_delta - nx_cs.  This is the per-pattern analog of
          // the EW divided_length = best_score + delta - nx_cost computed above
          // (which likewise keeps the attachment-invariant clip-internal steps),
          // so it equals extract_char_steps + add_clip_internal_steps but in
          // O(path) not O(n_node).  Pure-IW only (NA Pass-3 is attachment-dep).
          for (int p = 0; p < ds.n_patterns; ++p)
            divided_steps[p] = full_char_steps[p] + cs_delta[p] - nx_cs[p];
          if (iw_dirtychk) {
            // Oracle: per-clip compare against the extract+add baseline.
            std::vector<int> ref(ds.n_patterns, 0);
            extract_char_steps(tree, ds, ref);
            add_clip_internal_steps(tree, ds, clip_node, ref);
            for (int p = 0; p < ds.n_patterns; ++p)
              if (ref[p] != divided_steps[p]) {
                REprintf("DIRTYCHK clip=%d p=%d dirty=%d ref=%d\n",
                         clip_node, p, divided_steps[p], ref[p]);
                break;
              }
          }
        } else {
          std::fill(divided_steps.begin(), divided_steps.end(), 0);
          extract_char_steps(tree, ds, divided_steps);
          // extract_char_steps walks only the clipped (divided) postorder, so the
          // clipped subtree's internal steps are missing; add them back so base_iw
          // matches the EW divided_length convention (divided_tree + clip_internal)
          // and the indirect candidate base_iw + reconnect_delta is exact.
          // Pure-IW only.  The clipped subtree's internal Fitch length is
          // root- AND attachment-invariant, so adding it back makes the IW
          // indirect candidate (base_iw + reconnect_delta) EXACT.  NOT applied
          // under inapplicables: the NA Pass-3 internal step count is
          // attachment-DEPENDENT (down2 reads whole-tree uppass context), so it
          // cannot make the NA scan exact — it would be an unvalidated change to
          // production rooted-NA scoring.  NA instead reaches true unrooted
          // optima via the physical-reroot path (see the has_na dispatch below).
          if (!has_na) {
            add_clip_internal_steps(tree, ds, clip_node, divided_steps);
          }
        }
        base_iw = compute_weighted_score(ds, divided_steps);
        precompute_weighted_delta(ds, divided_steps, iw_delta);
      }
      if (iw_timing) {
        t_clip_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_clip).count();
        ++n_clips_t;
      }

      // Exact directional insertion edge sets for the pure-EW path; computed
      // once per clip from the current (clipped) main-tree downpass, then used
      // by both the SPR scan and the rerooting vroot cache below.
      if (use_directional) {
        if (l3b_active) {
          // L3b: patch the working buffers (== the per-pass intact base at
          // entry) to this divided tree, touching only the changed frontier.
          // l3b_changed records touched nodes for the restore after the scan.
          patch_insertion_edge_sets(tree, ds, nz, ns, up_base,
                                    edge_set_up, edge_set_buf,
                                    l3b_changed, l3b_worklist);
          if (l3b_oracle) {
            // The incremental result MUST equal a from-scratch recompute for
            // every node the divided tree reads (its preorder, minus root).
            // Abort on the first mismatch — trusting a wall number before this
            // passes is unsafe (the suppress-node transition is the risk).
            compute_insertion_edge_sets(tree, ds, l3b_oracle_es,
                                        l3b_oracle_up, l3b_oracle_pre);
            for (int D : l3b_oracle_pre) {
              if (D == tree.n_tip) continue;   // root has no edge_set
              size_t db = static_cast<size_t>(D) * tw;
              if (std::memcmp(&edge_set_buf[db], &l3b_oracle_es[db],
                              static_cast<size_t>(tw) * sizeof(uint64_t)) != 0) {
                Rf_error("L3B-ORACLE mismatch clip=%d nz=%d ns=%d node=%d "
                         "changed=%d (incremental edge_set != from-scratch)",
                         clip_node, nz, ns, D,
                         static_cast<int>(l3b_changed.size()));
              }
            }
            l3b_edges_sum +=
                static_cast<long long>(l3b_oracle_pre.size()) - 1;
          }
          ++l3b_patch_clips;
          l3b_changed_sum += static_cast<long long>(l3b_changed.size());
        } else {
          compute_insertion_edge_sets(tree, ds, edge_set_buf,
                                      edge_set_up, edge_set_pre);
        }
      }

      collect_main_edges(tree, main_edges);
      // Partial shuffle: seed the first few evaluation positions with edges
      // from across the tree so the bounded indirect scoring gets a tight
      // cutoff early.  Full O(n) shuffle has non-trivial overhead relative
      // to the per-candidate scoring cost; partial Fisher-Yates for a small
      // prefix keeps overhead negligible.
      {
        int ne = static_cast<int>(main_edges.size());
        int k = std::min(20, ne);
        for (int i = 0; i < k; ++i) {
          std::uniform_int_distribution<int> dist(i, ne - 1);
          std::swap(main_edges[i], main_edges[dist(rng)]);
        }
      }

      // Constraint: classify this clip against each constraint split
      if (constrained) {
        classify_clip_constraints(tree, clip_node, *cd);
      }

      // Find best (reroot, regraft) combination
      double best_candidate = HUGE_VAL;
      int best_above = -1, best_below = -1;
      int best_reroot_parent = -1, best_reroot_child = -1;

      // SPR candidates — with early termination (optimization #1)
      size_t clip_base = static_cast<size_t>(clip_node) * tree.total_words;
      const uint64_t* clip_prelim = &tree.prelim[clip_base];

      // EW/NA bail cutoff, maintained across this clip's SPR + reroot loops.
      // Recomputed ONLY inside an accept block (when best_candidate improves) —
      // divided_length is clip-constant, so this is byte-identical to the old
      // per-candidate `(best_candidate<HUGE_VAL)?(int)(best_candidate-
      // divided_length+1):INT_MAX`, just computed O(improvements) not
      // O(candidates).  Unused on the IW path (which bounds on best_candidate
      // directly).  best_candidate == HUGE_VAL here ⇒ INT_MAX matches the
      // old ternary's first-candidate branch.
      int cutoff = INT_MAX;

      std::chrono::steady_clock::time_point _t_scan;
      long long _e_spr = n_evaluated;
      if (iw_timing) _t_scan = std::chrono::steady_clock::now();
      // Lever #7 dispatch: in the "plain" regime (no NA, no sector / constraint
      // / collapsed / b2) EVERY per-candidate branch below except the identity
      // skip is provably dead, so route to the monomorphized template chosen
      // once by criterion flavour — EW (with the flat-kernel weight-class split)
      // or IW/XPIWE (strip-only, no flat variant). Byte-identical; ew_mono
      // default ON, TS_EW_MONO_OFF reverts to the general loop for BOTH.
      const bool mono_plain = ew_mono && !has_na
          && sector_mask == nullptr && !constrained && collapsed_all_zero
          && !b2_ceiling;
      const bool ew_mono_plain = mono_plain && !use_iw;
      const bool iw_mono_plain = mono_plain && use_iw;
      if (ew_mono_plain) {
        if (use_flat) {
          ++spr_mono_flat_fired;
#ifdef TS_EW_MONO_WRONG
          if (ew_mono_wrong)
            spr_scan_plain_ew<true, true>(main_edges, nz, ns, clip_prelim,
                edge_set_buf, tw, ds, divided_length, best_candidate,
                best_above, best_below, best_reroot_parent, best_reroot_child,
                n_evaluated, cutoff);
          else
#endif
            spr_scan_plain_ew<true, false>(main_edges, nz, ns, clip_prelim,
                edge_set_buf, tw, ds, divided_length, best_candidate,
                best_above, best_below, best_reroot_parent, best_reroot_child,
                n_evaluated, cutoff);
        } else {
          ++spr_mono_cached_fired;
#ifdef TS_EW_MONO_WRONG
          if (ew_mono_wrong)
            spr_scan_plain_ew<false, true>(main_edges, nz, ns, clip_prelim,
                edge_set_buf, tw, ds, divided_length, best_candidate,
                best_above, best_below, best_reroot_parent, best_reroot_child,
                n_evaluated, cutoff);
          else
#endif
            spr_scan_plain_ew<false, false>(main_edges, nz, ns, clip_prelim,
                edge_set_buf, tw, ds, divided_length, best_candidate,
                best_above, best_below, best_reroot_parent, best_reroot_child,
                n_evaluated, cutoff);
        }
      } else if (iw_mono_plain) {
        ++spr_mono_iw_fired;
#ifdef TS_EW_MONO_WRONG
        if (ew_mono_wrong)
          spr_scan_plain_iw<true>(main_edges, nz, ns, clip_prelim, edge_set_buf,
              tw, ds, base_iw, iw_delta, divided_length, best_candidate,
              best_above, best_below, best_reroot_parent, best_reroot_child,
              n_evaluated, cutoff);
        else
#endif
          spr_scan_plain_iw<false>(main_edges, nz, ns, clip_prelim, edge_set_buf,
              tw, ds, base_iw, iw_delta, divided_length, best_candidate,
              best_above, best_below, best_reroot_parent, best_reroot_child,
              n_evaluated, cutoff);
      } else {
      for (auto& [above, below] : main_edges) {
        if (above == nz && below == ns) continue;
        if (sector_mask && !(*sector_mask)[below]) continue;
        if (constrained && regraft_violates_constraint(below, *cd)) continue;

        // Collapsed-region regraft merging: skip interior collapsed edges.
        // If collapsed[below] == 1, the edge (above, below) is zero-length
        // and lies inside a collapsed region.  The boundary edge entering the
        // region (where collapsed[below] == 0 but the node is in the region)
        // is always evaluated, and it dominates interior positions because
        // its vroot includes states from outside the region.
        if (use_collapsed && collapsed[below]) {
          continue;
        }

        double candidate;
        if (has_na) {
          // NA-aware indirect with early termination
          if (use_iw) {
            candidate = indirect_na_iw_length_bounded(clip_prelim,
                clip_actives, tree, ds, above, below, base_iw, iw_delta,
                best_candidate);
          } else {
            int extra = use_flat
                ? fitch_na_indirect_bounded_flat(clip_prelim, clip_actives,
                      tree, ds, above, below, cutoff)
                : fitch_na_indirect_length_bounded(clip_prelim, clip_actives,
                      tree, ds, above, below, cutoff);
            candidate = divided_length + extra;
          }
        } else if (use_iw) {
          // Exact directional cost (mirrors the EW path): the edge set above
          // `below` is edge_set_buf[below], replacing the union-of-finals
          // approximation that hid improving IW moves.
          candidate = indirect_iw_length_cached(
              clip_prelim, &edge_set_buf[static_cast<size_t>(below) * tw],
              ds, base_iw, iw_delta, best_candidate);
        } else {
          // Exact directional cost: the edge set above `below` (= node_d) is
          // edge_set_buf[below], replacing the union-of-finals approximation.
          int extra = fitch_indirect_length_cached(
              clip_prelim, &edge_set_buf[static_cast<size_t>(below) * tw],
              ds, cutoff);
          candidate = divided_length + extra;
        }
        ++n_evaluated;
        if (b2_ceiling) {
          ++n_b2_eval;
          int pbelow = tree.parent[below];
          bool c1 = (pbelow != tree.n_tip);  // root-child edges aren't collapsible
          for (int b = 0; b < ds.n_blocks && c1; ++b) {
            if (ds.blocks[b].has_inapplicable) continue;
            if (tree.local_cost[static_cast<size_t>(pbelow) * ds.n_blocks + b])
              c1 = false;
          }
          if (c1) ++n_b2_aggr;  // Cond1 pass; reached here => not conservatively collapsed
        }
        if (candidate < best_candidate) {
          best_candidate = candidate;
          best_above = above;
          best_below = below;
          best_reroot_parent = -1;
          best_reroot_child = -1;
          cutoff = static_cast<int>(best_candidate - divided_length + 1);
        }
      }
      }  // end else (general SPR loop; ew_mono_plain template dispatch above)
      if (iw_timing) {
        t_spr_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_scan).count();
        n_spr_t += n_evaluated - _e_spr;
      }

      // TBR candidates (rerooting) — with vroot cache (optimization #4)
      if (clip_node >= tree.n_tip) {
        compute_from_above(tree, ds, clip_node, from_above);
        collect_subtree_edges(tree, clip_node, sub_edges);

        // Precompute vroot for all main edges (optimization #4).  EW and IW use
        // the exact directional edge set (vroot = edge_set_buf[below]); NA keeps
        // the union-of-finals form its cached scorer expects.
        int n_main = static_cast<int>(main_edges.size());
        if (use_directional) {
          vroot_cache.resize(static_cast<size_t>(n_main) * tw);
          for (int ei = 0; ei < n_main; ++ei) {
            int d = main_edges[ei].second;  // child endpoint (node_d)
            std::memcpy(&vroot_cache[static_cast<size_t>(ei) * tw],
                        &edge_set_buf[static_cast<size_t>(d) * tw],
                        static_cast<size_t>(tw) * sizeof(uint64_t));
          }
        } else {
          precompute_vroot_cache(tree, main_edges, vroot_cache);
        }

        // For NA: precompute per-edge below_actives (OR of applicable
        // subtree_actives words for node_d of each edge)
        if (has_na) {
          below_actives_cache.resize(
              static_cast<size_t>(n_main) * ds.n_blocks);
          for (int ei = 0; ei < n_main; ++ei) {
            int d = main_edges[ei].second;
            size_t d_base = static_cast<size_t>(d) * tree.total_words;
            for (int b_i = 0; b_i < ds.n_blocks; ++b_i) {
              if (!ds.blocks[b_i].has_inapplicable) {
                below_actives_cache[
                    static_cast<size_t>(ei) * ds.n_blocks + b_i] = 0;
                continue;
              }
              int off = ds.block_word_offset[b_i];
              int k = ds.blocks[b_i].n_states;
              uint64_t ba = 0;
              for (int s = 1; s < k; ++s) {
                ba |= tree.subtree_actives[d_base + off + s];
              }
              below_actives_cache[
                  static_cast<size_t>(ei) * ds.n_blocks + b_i] = ba;
            }
          }
        }

        // Phase 3A: Symmetry-breaking — deduplicate equivalent rerootings.
        // Seed with clip_prelim hash (SPR case already evaluated above).
        // Reset the pre-loop dedup table (O(1) generation bump, no alloc).
        seen_vp_hashes.reset(sub_edges.size());
        seen_vp_hashes.insert(fast_hash(clip_prelim, tree.total_words));

        // Precompute the non-skipped main-edge indices ONCE per clip.  The skip
        // predicate (nz/ns identity, sector_mask, constraint, collapsed) is
        // sub_edge-INVARIANT — it depends only on ei via main_edges[ei] — yet
        // both reroot inner loops below re-evaluate it for every (sub_edge, ei)
        // pair.  Filtering here removes that n_sub_edges x redundancy.  kept_ei
        // is in ascending ei order, so the candidate sequence (hence the
        // strict-< tie-break) is byte-identical to the per-candidate skips.
        kept_ei.clear();
        for (int ei = 0; ei < n_main; ++ei) {
          const int ab = main_edges[ei].first, bl = main_edges[ei].second;
          const bool skip = ((ab == nz && bl == ns) && !do_reroot)
              || (sector_mask && !(*sector_mask)[bl])
              || (constrained && regraft_violates_constraint(bl, *cd))
              || (use_collapsed && collapsed[bl]);
          if (!skip) kept_ei.push_back(ei);
        }

        std::chrono::steady_clock::time_point _t_scan2;
        long long _e_rer = n_evaluated;
        if (iw_timing) _t_scan2 = std::chrono::steady_clock::now();
        for (auto& [sp, sc] : sub_edges) {
          if (sp == clip_node) continue;

          fitch_join_states(
              &from_above[static_cast<size_t>(sc) * tree.total_words],
              &tree.prelim[static_cast<size_t>(sc) * tree.total_words],
              virtual_prelim.data(), ds);

          // Fast path: skip if virtual_prelim matches SPR case exactly
          if (std::memcmp(virtual_prelim.data(), clip_prelim,
                          tree.total_words * sizeof(uint64_t)) == 0) {
            continue;
          }

          // Hash-based dedup: skip if we've seen this virtual_prelim before
          uint64_t vp_hash = fast_hash(virtual_prelim.data(),
                                         tree.total_words);
          if (!seen_vp_hashes.insert(vp_hash)) {
            continue;  // Already evaluated an equivalent rerooting
          }

          if (use_flat && !use_iw) {
            // === EW flat 4-wide batch (T-245) ===
            // Collect up to 4 non-skipped candidates per iteration and
            // evaluate them simultaneously.  4 independent vroot_cache rows
            // are read in parallel, hiding L2 latency for large trees.
            // Batch 4 at a time from the pre-filtered kept_ei (skips already
            // applied once per clip above) — byte-identical b_ei sequence.
            size_t ki = 0;
            const size_t n_kept = kept_ei.size();
            while (ki < n_kept) {
              int b_ei[4];
              int b_n = 0;
              while (ki < n_kept && b_n < 4) {
                b_ei[b_n++] = kept_ei[ki++];
              }
              if (b_n == 0) break;

              // cutoff is maintained across the clip (recomputed only on
              // improvement); byte-identical to the old per-batch recompute.
              int cutoff_b = cutoff;
              // Initialise to cutoff_b so partial-batch trailing slots
              // never accidentally improve best_candidate.
              int scores[4] = {cutoff_b, cutoff_b, cutoff_b, cutoff_b};

              if (b_n == 4) {
                // Full 4-wide batch: all 4 vroot_cache rows in flight.
                if (!has_na) {
                  fitch_indirect_cached_flat_x4(
                      virtual_prelim.data(),
                      &vroot_cache[static_cast<size_t>(b_ei[0]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[1]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[2]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[3]) * tw],
                      ds, cutoff_b, scores);
                } else {
                  fitch_na_indirect_cached_flat_x4(
                      virtual_prelim.data(), clip_actives,
                      &vroot_cache[static_cast<size_t>(b_ei[0]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[1]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[2]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[3]) * tw],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[0]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[1]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[2]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[3]) * ds.n_blocks],
                      ds, cutoff_b, scores);
                }
              } else {
                // Scalar fallback for trailing partial batch (< 4 valid).
                for (int k = 0; k < b_n; ++k) {
                  scores[k] = has_na
                      ? fitch_na_indirect_cached_flat(
                            virtual_prelim.data(), clip_actives,
                            &vroot_cache[static_cast<size_t>(b_ei[k]) * tw],
                            &below_actives_cache[
                                static_cast<size_t>(b_ei[k]) * ds.n_blocks],
                            ds, cutoff_b)
                      : fitch_indirect_cached_flat(
                            virtual_prelim.data(),
                            &vroot_cache[static_cast<size_t>(b_ei[k]) * tw],
                            ds, cutoff_b);
                }
              }

              n_evaluated += b_n;
              for (int k = 0; k < b_n; ++k) {
                double candidate = divided_length + scores[k];
                if (candidate < best_candidate) {
                  best_candidate = candidate;
                  best_above = main_edges[b_ei[k]].first;
                  best_below = main_edges[b_ei[k]].second;
                  best_reroot_parent = sp;
                  best_reroot_child = sc;
                  cutoff = static_cast<int>(best_candidate - divided_length + 1);
                }
              }
            }
          } else if (iw_family && iw_x4) {
            // === IW/XPIWE 4-wide batch (pure-IW/XPIWE and IW+NA; T-245 ILP) ===
            // Same candidate-selection + bookkeeping as the EW flat-x4 above,
            // but scores via the IW double-accumulator kernels.  has_na selects
            // the NA-aware kernel (indirect_na_iw_cached_flat_x4 full batch /
            // indirect_na_iw_length_cached scalar tail), which ANDs in
            // clip_has_active + each candidate's below_actives row; the no-NA
            // path is byte-identical to before this branch handled NA.  When
            // iw_x4 is disabled (TS_IW_NOX4) this whole branch is skipped and
            // both NA and no-NA IW fall through to the scalar `else` below.
            int ei = 0;
            while (ei < n_main) {
              int b_ei[4];
              int b_n = 0;
              while (ei < n_main && b_n < 4) {
                auto& [ab, bl] = main_edges[ei];
                bool skip = ((ab == nz && bl == ns) && !do_reroot)
                    || (sector_mask && !(*sector_mask)[bl])
                    || (constrained && regraft_violates_constraint(bl, *cd))
                    || (!collapsed.empty() && collapsed[bl]);
                if (!skip) b_ei[b_n++] = ei;
                ++ei;
              }
              if (b_n == 0) break;

              double scores[4] = {best_candidate, best_candidate,
                                  best_candidate, best_candidate};
              if (b_n == 4) {
                if (has_na) {
                  indirect_na_iw_cached_flat_x4(
                      virtual_prelim.data(), clip_actives,
                      &vroot_cache[static_cast<size_t>(b_ei[0]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[1]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[2]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[3]) * tw],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[0]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[1]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[2]) * ds.n_blocks],
                      &below_actives_cache[
                          static_cast<size_t>(b_ei[3]) * ds.n_blocks],
                      ds, base_iw, iw_delta, best_candidate, scores);
                } else {
                  indirect_iw_cached_flat_x4(
                      virtual_prelim.data(),
                      &vroot_cache[static_cast<size_t>(b_ei[0]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[1]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[2]) * tw],
                      &vroot_cache[static_cast<size_t>(b_ei[3]) * tw],
                      ds, base_iw, iw_delta, best_candidate, scores);
                }
              } else {
                for (int k = 0; k < b_n; ++k) {
                  scores[k] = has_na
                      ? indirect_na_iw_length_cached(
                            virtual_prelim.data(), clip_actives,
                            &vroot_cache[static_cast<size_t>(b_ei[k]) * tw],
                            &below_actives_cache[
                                static_cast<size_t>(b_ei[k]) * ds.n_blocks],
                            ds, base_iw, iw_delta, best_candidate)
                      : indirect_iw_length_cached(
                            virtual_prelim.data(),
                            &vroot_cache[static_cast<size_t>(b_ei[k]) * tw],
                            ds, base_iw, iw_delta, best_candidate);
                }
              }

              n_evaluated += b_n;
              for (int k = 0; k < b_n; ++k) {
                double candidate = scores[k];
                if (candidate < best_candidate) {
                  best_candidate = candidate;
                  best_above = main_edges[b_ei[k]].first;
                  best_below = main_edges[b_ei[k]].second;
                  best_reroot_parent = sp;
                  best_reroot_child = sc;
                }
              }
            }
          } else {
            // === Scalar path (IW, or ratchet with upweight_mask) ===
            const int n_kept = static_cast<int>(kept_ei.size());
            for (int ki = 0; ki < n_kept; ++ki) {
              const int ei = kept_ei[ki];
              // Prefetch the NEXT kept vroot row.
              // At 180 tips vroot_cache is ~140 KB (L2); prefetch hides
              // the ~10-cycle L2 latency. No-op on small trees (L1-resident).
              if (ki + 2 < n_kept) {
                const size_t pf =
                    static_cast<size_t>(kept_ei[ki + 2]) * tree.total_words;
#if defined(__GNUC__) || defined(__clang__)
                __builtin_prefetch(&vroot_cache[pf], 0, 0);
#elif defined(_MSC_VER) && defined(TS_SIMD_SSE2)
                _mm_prefetch(reinterpret_cast<const char*>(&vroot_cache[pf]),
                             _MM_HINT_T0);
#endif
              }
              auto& [above, below] = main_edges[ei];
              // Skip predicate (nz/ns, sector, constraint, collapsed) already
              // applied once per clip via kept_ei — no per-candidate re-check.
              double candidate;
              if (has_na) {
                if (use_iw) {
                  candidate = indirect_na_iw_length_cached(
                      virtual_prelim.data(), clip_actives,
                      &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                      &below_actives_cache[
                          static_cast<size_t>(ei) * ds.n_blocks],
                      ds, base_iw, iw_delta, best_candidate);
                } else {
                  int extra = fitch_na_indirect_length_cached(
                      virtual_prelim.data(), clip_actives,
                      &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                      &below_actives_cache[
                          static_cast<size_t>(ei) * ds.n_blocks],
                      ds, cutoff);
                  candidate = divided_length + extra;
                }
              } else if (use_iw) {
                double iw_cutoff = best_candidate;
                candidate = indirect_iw_length_cached(
                    virtual_prelim.data(),
                    &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                    ds, base_iw, iw_delta, iw_cutoff);
              } else {
                int extra = fitch_indirect_length_cached(
                    virtual_prelim.data(),
                    &vroot_cache[static_cast<size_t>(ei) * tree.total_words],
                    ds, cutoff);
                candidate = divided_length + extra;
              }
              ++n_evaluated;
              if (candidate < best_candidate) {
                best_candidate = candidate;
                best_above = above;
                best_below = below;
                best_reroot_parent = sp;
                best_reroot_child = sc;
                cutoff = static_cast<int>(best_candidate - divided_length + 1);
              }
            }
          }
        }
        if (iw_timing) {
          t_rer_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::steady_clock::now() - _t_scan2).count();
          n_rer_t += n_evaluated - _e_rer;
        }
      }

      // L3b: restore the working edge_set/up buffers to the per-pass intact base
      // over exactly the nodes this clip patched, so the next clip starts from
      // base again.  The changed-node list IS the undo log (no separate save
      // buffer).  Both buffers are fully consumed by now (SPR scan read
      // edge_set_buf; the reroot vroot_cache copied its rows), and edge_set_up
      // must be restored too — the next patch reads ancestor up-messages from it
      // and assumes off-frontier nodes hold the base value.
      if (l3b_active) {
        for (int D : l3b_changed) {
          size_t db = static_cast<size_t>(D) * tw;
          std::memcpy(&edge_set_buf[db], &edge_set_base[db],
                      static_cast<size_t>(tw) * sizeof(uint64_t));
          std::memcpy(&edge_set_up[db], &up_base[db],
                      static_cast<size_t>(tw) * sizeof(uint64_t));
        }
        l3b_changed.clear();
      }

      // --- Phase 2: Restore original tree, verify best candidate ---
      // Restore states from pre-allocated undo (clip_undo_stack is empty)
      tree.restore_prealloc_undo();
      tree.spr_unclip();
      // Restore saved postorder (topology identical to pre-clip state)
      tree.postorder.assign(saved_postorder.begin(), saved_postorder.end());

      // DIAGNOSTIC (env TS_REVERT_CHECK): within a pass no move is accepted, so
      // the tree is invariant and the clip-undo restore above MUST equal the
      // pass-start snapshot (snap/state_snap, saved at 782-783).  Any mismatch =
      // clip-undo (spr_clip/unclip + saved_postorder) is not a perfect inverse,
      // leaving residue that accumulates across clips -- the latent bug the
      // abandonment edit exposed.  Reports which array diverges.
      if (revert_check) {
        bool topo = (tree.parent == snap.parent && tree.left == snap.left &&
                     tree.right == snap.right);
        size_t sb = state_snap.prelim.size() * sizeof(uint64_t);
        size_t cb = state_snap.local_cost.size() * sizeof(uint64_t);
        bool pre = sb == 0 || std::memcmp(tree.prelim.data(),
                     state_snap.prelim.data(), sb) == 0;
        bool fin = sb == 0 || std::memcmp(tree.final_.data(),
                     state_snap.final_.data(), sb) == 0;
        bool lc  = cb == 0 || std::memcmp(tree.local_cost.data(),
                     state_snap.local_cost.data(), cb) == 0;
        bool po  = (tree.postorder == state_snap.postorder);
        if (!(topo && pre && fin && lc && po))
          REprintf("CLIPUNDO-MISMATCH clip=%d topo=%d prelim=%d final=%d lc=%d post=%d\n",
                   clip_node, (int)topo, (int)pre, (int)fin, (int)lc, (int)po);
      }

      bool dominated = (best_candidate > best_score + eps) ||
                        (best_candidate > best_score - eps && !params.accept_equal);

      bool accepted = false;

      if (!dominated && best_above >= 0) {
        // Topology and state snapshot already saved at the top of the
        // while loop (optimization #7). No per-candidate save needed.

        bool ok = apply_tbr_move(tree, clip_node,
                                  best_reroot_parent, best_reroot_child,
                                  best_above, best_below);
        // Topology mutated; states no longer match best_score.
        score_fresh = false;

        // apply_tbr_move's own success flag (`ok`) is the functional check.
        // The full topology walk is debug-only paranoia (Tier 3a): compiled
        // out of release (NDEBUG) builds, where it cost ~2-3% on the per-accept
        // path. apply_tbr_move is trusted to produce a valid tree in release.
        bool topo_ok = ok;
#ifndef NDEBUG
        topo_ok = topo_ok && validate_topology(tree);
#endif
        if (!topo_ok) {
          restore_topology(tree, snap);
          state_snap.restore(tree);
          score_fresh = true;
          continue;
        }

        tree.build_postorder_prealloc(work_stack);

        // T-300: dirty-set incremental rescore for SPR moves.  The two
        // affected nodes after apply_tbr_move are nz (clip grandparent,
        // children changed: nx -> ns) and nx (regraft point, children
        // changed to {clip_node, below}).  fitch_dirty_downpass updates
        // every node on the union of paths nz->root and nx->root exactly
        // once in postorder; sums correctly with no shared-ancestor
        // ambiguity.  TBR moves with non-trivial rerooting and NA
        // datasets fall back to full_rescore.
        bool is_spr = (best_reroot_parent < 0 || best_reroot_parent == clip_node);
        double actual;
        if (is_spr && !has_na && incremental_ok) {
          int delta = fitch_dirty_downpass(tree, ds, nz, nx);
          fitch_dirty_uppass(tree, ds, nz, nx);
          if (use_iw) {
            std::fill(divided_steps.begin(), divided_steps.end(), 0);
            extract_char_steps(tree, ds, divided_steps);
            actual = compute_weighted_score(ds, divided_steps);
          } else {
            actual = best_score + static_cast<double>(delta);
          }
        } else if (is_spr && has_na && incremental_ok) {
          // T-300 NA variant: dirty-set Pass 1 + Pass 2 instead of full
          // rescore.  Pass 3 still runs over the full tree because it
          // populates internal down2 (read by extract_char_steps) and
          // counts NA-block steps directly.  Savings come from skipping
          // Pass 1 + Pass 2 on off-dirty nodes.
          fitch_na_dirty_downpass(tree, ds, nz, nx);
          fitch_na_dirty_uppass(tree, ds, nz, nx);
          int ew_total = fitch_na_pass3_score(tree, ds);
          if (use_iw) {
            std::fill(divided_steps.begin(), divided_steps.end(), 0);
            extract_char_steps(tree, ds, divided_steps);
            actual = compute_weighted_score(ds, divided_steps);
          } else {
            // EW score must include ds.ew_offset (topology-independent
            // step count added by fitch_score_ew on top of fitch_na_score)
            // — bug found 2026-05-19: omitting this produced systematic
            // diff=−3 against full_rescore (Vinther2008 offset = 3).
            actual = static_cast<double>(ew_total) + ds.ew_offset;
          }
        } else {
          // Non-trivial TBR rerooting, or a scoring mode whose incremental
          // delta is not exact (HSJ/XFORM, see incremental_ok): recompute the
          // authoritative score via score_tree().
          actual = full_rescore(tree, ds);
        }

        // DIAGNOSTIC (env TS_IW_SCANCHK): compare the scan's predicted
        // best_candidate against the authoritative post-apply score for EVERY
        // scorer.  Pure EW should be 0 (its indirect length is exact); a
        // mismatch under IW or NA flags the clip-internal omission (see
        // add_clip_internal_steps).  reroot=1 => TBR rerooting accept, 0 => SPR.
        // No-op unless the env var is set.
        if (iw_scanchk &&
            std::fabs(actual - best_candidate) > 1e-6) {
          REprintf("SCANCHK-MISMATCH mode=%s reroot=%d pred=%.5f actual=%.5f diff=%.5f\n",
                   has_na ? (use_iw ? "NA+IW" : "NA+EW") : (use_iw ? "IW" : "EW"),
                   is_spr ? 0 : 1, best_candidate, actual, actual - best_candidate);
        }

        // Post-hoc constraint validation: TBR rerooting can break
        // splits that were classified as UNCONSTRAINED during the
        // clip phase (the rerooting changes which constraint tips
        // end up on which side of the attachment edge).  Reject
        // any move that introduces a constraint violation.
        if (constrained) {
          map_constraint_nodes(tree, *cd);
          bool violation = false;
          for (int _s = 0; _s < cd->n_splits; ++_s) {
            if (cd->constraint_node[_s] < 0) {
              violation = true;
              break;
            }
          }
          if (violation) {
            restore_topology(tree, snap);
            state_snap.restore(tree);
            score_fresh = true;
            map_constraint_nodes(tree, *cd);
            compute_dfs_timestamps(tree, *cd);
            continue;
          }
        }

        // Compute topology hash for tabu checking
        uint64_t tree_hash = 0;
        if (tabu.active()) {
          tree_hash = hash_tree(tree);
        }

        if (actual < best_score - eps) {
          // Always accept strict improvements (but record in tabu)
          if (tabu.active()) tabu.insert(tree_hash);
          best_score = actual;
          score_fresh = true;
          ++n_accepted;
          hits = 1;
          accepted = true;
          keep_going = true;
          if (constrained) {
            compute_dfs_timestamps(tree, *cd);
          }
          if (collect_pool) {
            // Dedup on collapsed topology (matches the main-loop pool).  States
            // are valid here: the accepted move's incremental passes (or
            // full_rescore) left prelim/final_ coherent for the whole tree —
            // the same invariant the production compute_collapsed_flags call on
            // the post-accept path (below, gated on !collapsed.empty()) relies
            // on.  No-op vs plain add() when no edge is collapsible.
            if (enum_resolved) {
              collect_pool->add(tree, actual);
            } else {
              compute_collapsed_flags(tree, ds, collect_collapsed);
              collect_pool->add_collapsed(tree, actual, collect_collapsed);
            }
          }
        } else if (std::fabs(actual - best_score) <= eps
                   && params.accept_equal
                   && hits <= params.max_hits) {
          // Equal-score move: reject if tabu
          if (tabu.active() && tabu.contains(tree_hash)) {
            // Topology already visited — restore and skip.
            restore_topology(tree, snap);
            state_snap.restore(tree);
            score_fresh = true;
            // Re-sync constraint metadata to the restored topology: the
            // violation check above ran map_constraint_nodes() on the post-move
            // tree, so cd->constraint_node / DFS timestamps are stale after
            // restoration (same hazard handled in the !accepted path below).
            if (constrained) {
              map_constraint_nodes(tree, *cd);
              compute_dfs_timestamps(tree, *cd);
            }
            continue;
          }
          if (tabu.active()) tabu.insert(tree_hash);
          // Adopt `actual` so that best_score stays bit-exact for the new
          // topology (it is already within eps of the prior value).
          best_score = actual;
          score_fresh = true;
          ++hits;
          ++n_accepted;
          accepted = true;
          keep_going = true;
          if (constrained) {
            compute_dfs_timestamps(tree, *cd);
          }
          if (collect_pool) {
            // Dedup on collapsed topology (matches the main-loop pool).  States
            // are valid here: the accepted move's incremental passes (or
            // full_rescore) left prelim/final_ coherent for the whole tree —
            // the same invariant the production compute_collapsed_flags call on
            // the post-accept path (below, gated on !collapsed.empty()) relies
            // on.  No-op vs plain add() when no edge is collapsible.
            if (enum_resolved) {
              collect_pool->add(tree, actual);
            } else {
              compute_collapsed_flags(tree, ds, collect_collapsed);
              collect_pool->add_collapsed(tree, actual, collect_collapsed);
            }
          }
        }

        if (!accepted) {
          // Optimization #3: restore topology + states without full_rescore
          restore_topology(tree, snap);
          state_snap.restore(tree);
          score_fresh = true;
          // state_snap.restore() already restored postorder via memcpy
          // Re-sync constraint metadata to the restored topology.  When a
          // constrained move passes the violation check but fails the score
          // check, map_constraint_nodes() was already called for the
          // post-move tree; after restoration cd->constraint_node is stale
          // relative to the pre-move topology.  Without this re-mapping,
          // the next clip's classify_clip_constraints() can produce false-
          // positive or false-negative constraint violations.
          if (constrained) {
            map_constraint_nodes(tree, *cd);
            compute_dfs_timestamps(tree, *cd);
          }
        }
      }

      if (keep_going) {
        accepted_clip_size = clip_size;
        pass_candidates_evaluated += (n_evaluated - clip_evals_before);
        // Recompute collapsed regions after the accepted move (states are
        // valid from full_rescore in the accept path above).
        if (!collapsed.empty()) {
          compute_collapsed_flags(tree, ds, collapsed);
          refresh_collapsed_all_zero();  // lever #7 gate must not go stale
        }
        // Optimization #6: don't reshuffle after acceptance — the topology
        // changed near this clip, so re-trying the same ordering focuses
        // on the productive region.
        need_shuffle = false;

        // L3b base-incremental-update: bring the base + working buffers current
        // for the post-move tree WITHOUT the per-pass full recompute (the fixed
        // cost that drowned the clip-patch win in accept-dense passes).  prelim
        // is coherent here (the accept path ran dirty passes or full_rescore).
        // SPR moves only (local relocation of the clip subtree); TBR-reroot
        // accepts leave base_current false → full recompute at the next pass.
        if (l3b_active) {
          const bool spr_move =
              (best_reroot_parent < 0 || best_reroot_parent == clip_node);
          if (l3b_baseincr && spr_move) {
            update_base_after_spr_move(tree, ds, nz, nx, ns,
                best_above, best_below, up_base, edge_set_base,
                l3b_changed, l3b_worklist, l3b_base_tmp);
            if (l3b_oracle) {
              // The incrementally-updated base MUST equal a from-scratch
              // recompute of the post-move tree, full array.  Abort on mismatch.
              compute_insertion_edge_sets(tree, ds, l3b_oracle_es,
                                          l3b_oracle_up, l3b_oracle_pre);
              for (int D : l3b_oracle_pre) {
                if (D == tree.n_tip) continue;
                size_t db = static_cast<size_t>(D) * tw;
                if (std::memcmp(&edge_set_base[db], &l3b_oracle_es[db],
                        static_cast<size_t>(tw) * sizeof(uint64_t)) != 0) {
                  Rf_error("L3B-BASE-ORACLE mismatch clip=%d nz=%d nx=%d "
                           "above=%d below=%d node=%d (base != from-scratch)",
                           clip_node, nz, nx, best_above, best_below, D);
                }
              }
            }
            // Sync the working buffers to the updated base over the touched
            // nodes (buf == the OLD base after this clip's restore above).
            for (int D : l3b_changed) {
              size_t db = static_cast<size_t>(D) * tw;
              std::memcpy(&edge_set_buf[db], &edge_set_base[db],
                          static_cast<size_t>(tw) * sizeof(uint64_t));
              std::memcpy(&edge_set_up[db], &up_base[db],
                          static_cast<size_t>(tw) * sizeof(uint64_t));
            }
            l3b_base_current = true;
          } else {
            l3b_base_current = false;   // TBR-reroot accept → full recompute
          }
        }

        if (params.max_accepted_changes > 0
            && n_accepted >= params.max_accepted_changes) {
          keep_going = false;
        }
        break;
      }

      pass_candidates_evaluated += (n_evaluated - clip_evals_before);

      if (ts::check_interrupt()) break;
      ++clips_since_timeout_check;
      if (check_timeout && clips_since_timeout_check >= timeout_interval) {
        clips_since_timeout_check = 0;
        if (check_timeout()) { timed_out = true; break; }
      }
    }

    // Record per-pass diagnostics
    if (params.diagnostics) {
      diag_records.push_back({
        pass_index,
        /*productive=*/ accepted_clip_size > 0,
        accepted_clip_size,
        pass_clips_tried,
        pass_candidates_evaluated
      });
    }
    ++pass_index;

    if (params.max_accepted_changes > 0
        && n_accepted >= params.max_accepted_changes) {
      break;
    }
  }  // end inner convergence while

  // ----- outer reroot control -----
  if (!do_reroot || timed_out
      || (params.max_accepted_changes > 0
          && n_accepted >= params.max_accepted_changes)) break;

  // Direct in-pass root-edge enumeration for ALL scorers (EW fast additive;
  // IW/NA apply+rescore).  Clean => no improving move on ANY edge (the inner
  // loop certified the 2n-4 non-root edges) => true unrooted-TBR optimum.
  // TS_PHYS_REROOT forces the legacy physical-reroot sweep (validation ref).
  if (!phys_reroot) {
    // EW/IW: the indirect scan is exact, so the inner loop already certified
    // the 2n-4 non-root edges — only the root edge remains (try_root_edge_moves,
    // fast additive for EW / apply+rescore for IW).  NA: the indirect scan is
    // only approximate, so an EXACT full-neighbourhood sweep is required to
    // certify a true unrooted-TBR optimum (see exact_verify_sweep).
    bool improved = has_na
        ? exact_verify_sweep(tree, ds, best_score)
        : try_root_edge_moves(tree, ds, best_score, ew_directional);
    if (!improved) break;
    score_fresh = true;
    if (!collapsed.empty()) {
      if (collapse_aggr) compute_collapsed_flags_aggressive(tree, ds, collapsed);
      else               compute_collapsed_flags(tree, ds, collapsed);
      refresh_collapsed_all_zero();  // lever #7 gate must not go stale
    }
    continue;                              // re-descend from the improved tree
  }

  // Legacy physical-reroot sweep (TS_PHYS_REROOT=1): scorer-agnostic, known
  // complete; kept as the reference for validating the direct path above.
  if (!first_descent) {
    // Did the descent since the last re-root strictly improve?
    if (best_score < reroot_prev - eps) reroot_clean = 0;
    else ++reroot_clean;
    // A full cycle of distinct tip rootings with no improvement => optimal
    // under every rooting => a true unrooted-TBR local optimum.
    if (reroot_clean >= tree.n_tip) break;
  }
  first_descent = false;
  reroot_prev = best_score;
  reroot_at_tip(tree, reroot_tip);
  reroot_tip = (reroot_tip + 1) % tree.n_tip;
  best_score = full_rescore(tree, ds);   // root-invariant; refreshes states
  score_fresh = true;
  if (!collapsed.empty()) {
    compute_collapsed_flags(tree, ds, collapsed);
    refresh_collapsed_all_zero();  // lever #7 gate must not go stale (legacy reroot)
  }
  }  // end outer reroot for(;;)

  tree.prealloc_undo = nullptr;

  // States and best_score are kept in sync across every accepted move and
  // every state_snap.restore (`score_fresh` invariant). The trailing
  // full_rescore is therefore only needed if a code path left them stale —
  // it acts as a safety net.
  if (!score_fresh) {
    best_score = full_rescore(tree, ds);
  }

  bool converged = !(params.max_accepted_changes > 0
                     && n_accepted >= params.max_accepted_changes);

  // Accumulate candidate count into the dataset-level diagnostic counter
  // (one add per search call, not per candidate). See DataSet docs.
  ds.n_candidates_evaluated += n_evaluated;

  if (b2_ceiling && n_b2_eval > 0) {
    Rprintf("[B2_CEILING] aggressive-skippable SPR-regraft candidates: "
            "%lld / %lld evaluated = %.3f%% (work-weighted ceiling)\n",
            n_b2_aggr, n_b2_eval, 100.0 * n_b2_aggr / n_b2_eval);
  }
  if (iw_timing) {
    REprintf("IWT regime=%s clips=%.0f clip_us/clip=%.3f | "
             "SPR n=%.0f %.2fns | REROOT n=%.0f %.2fns | total %.1fms | "
             "sprmono flat=%lld cached=%lld iw=%lld/%.0f\n",
             use_iw ? "IW" : "EW", (double)n_clips_t,
             n_clips_t ? t_clip_ns / 1000.0 / (double)n_clips_t : 0.0,
             (double)n_spr_t, n_spr_t ? (double)t_spr_ns / (double)n_spr_t : 0.0,
             (double)n_rer_t, n_rer_t ? (double)t_rer_ns / (double)n_rer_t : 0.0,
             (t_clip_ns + t_spr_ns + t_rer_ns) / 1e6,
             spr_mono_flat_fired, spr_mono_cached_fired, spr_mono_iw_fired,
             (double)n_clips_t);
  }

  // L3b footprint cross-check: the mean patched-node fraction should track the
  // fp_ref measured in dev/profiling/l3b-footprint-482.md (a discovery that
  // over- or under-reaches would diverge).  Reported only when TS_L3B_STATS is
  // set; l3b_edges_sum is populated only under the oracle.
  if (l3b_active && l3b_patch_clips > 0 && std::getenv("TS_L3B_STATS")) {
    double meanChanged = (double)l3b_changed_sum / (double)l3b_patch_clips;
    double fp = l3b_edges_sum > 0
        ? (double)l3b_changed_sum / (double)l3b_edges_sum : -1.0;
    REprintf("[L3B_STATS] patch_clips=%lld mean_changed=%.1f "
             "mean_edges=%.1f fp_changed=%.3f\n",
             l3b_patch_clips, meanChanged,
             l3b_edges_sum > 0 ? (double)l3b_edges_sum / (double)l3b_patch_clips
                               : -1.0,
             fp);
  }

  return TBRResult{best_score, n_accepted, n_evaluated, n_zero_skipped,
                   converged, std::move(diag_records)};
}

} // namespace ts
