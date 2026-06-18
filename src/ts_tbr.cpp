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
static void reroot_at_tip(TreeState& tree, int t) {
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
  static thread_local std::vector<int> stack;
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
  static thread_local std::vector<int> stack;
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
  static thread_local std::vector<int> preorder;
  preorder.clear();
  {
    static thread_local std::vector<int> stack;
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
  std::vector<uint8_t> collapsed;
  if (!collect_pool) {
    compute_collapsed_flags(tree, ds, collapsed);
  }

  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;

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
  std::vector<uint64_t> edge_set_buf;

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
  int reroot_tip = 0;
  int reroot_clean = 0;            // consecutive reroots with no strict gain
  double reroot_prev = HUGE_VAL;
  bool first_descent = true;

  for (;;) {
   keep_going = true;
  while (keep_going && !timed_out) {
    keep_going = false;

    // Optimization #7: save state snapshot once per pass, not per candidate.
    // After a rejected move, state_snap.restore() returns the tree to exactly
    // the state saved here. The per-candidate save was redundant: consecutive
    // rejections all restore to the same state. Re-saving only happens when
    // the while loop restarts after an accepted move.
    save_topology(tree, snap);
    state_snap.save(tree);

    // Reset per-pass diagnostic counters
    pass_clips_tried = 0;
    pass_candidates_evaluated = 0;
    accepted_clip_size = 0;

    // Recompute subtree sizes (needed for smaller-subtree filtering
    // and for clip ordering strategies)
    compute_subtree_sizes(tree, subtree_sizes);

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

      // Optimization #2: skip clips of the larger subtree.
      // For each edge, only clip the side with fewer tips.
      int clip_size = subtree_sizes[clip_node];
      if (clip_size > tree.n_tip / 2) continue;

      // Skip collapsed edges: zero-length edge where clipping provably
      // cannot improve the score. Works for EW, IW, Profile, and NA.
      // Disabled during MPT enumeration (collapsed is empty).
      if (!collapsed.empty() && collapsed[clip_node]) {
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
        int delta = fitch_incremental_downpass(tree, ds, nz);
        fitch_incremental_uppass(tree, ds, nz);

        int nx_cost = 0;
        for (int b = 0; b < ds.n_blocks; ++b) {
          uint64_t lc = tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b];
          int nu = popcount64(lc);
          if (ds.blocks[b].upweight_mask) nu += popcount64(lc & ds.blocks[b].upweight_mask);
          nx_cost += ds.blocks[b].weight * nu;
        }
        divided_length = best_score + delta - nx_cost;
      }

      // For weighted scoring (IW or profile): precompute base score and deltas
      double base_iw = 0.0;
      if (use_iw) {
        std::fill(divided_steps.begin(), divided_steps.end(), 0);
        extract_char_steps(tree, ds, divided_steps);
        base_iw = compute_weighted_score(ds, divided_steps);
        precompute_weighted_delta(ds, divided_steps, iw_delta);
      }

      // Exact directional insertion edge sets for the pure-EW path; computed
      // once per clip from the current (clipped) main-tree downpass, then used
      // by both the SPR scan and the rerooting vroot cache below.
      if (ew_directional) {
        compute_insertion_edge_sets(tree, ds, edge_set_buf);
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
        if (!collapsed.empty() && collapsed[below]) {
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
            int cutoff = (best_candidate < HUGE_VAL)
                ? static_cast<int>(best_candidate - divided_length + 1)
                : INT_MAX;
            int extra = use_flat
                ? fitch_na_indirect_bounded_flat(clip_prelim, clip_actives,
                      tree, ds, above, below, cutoff)
                : fitch_na_indirect_length_bounded(clip_prelim, clip_actives,
                      tree, ds, above, below, cutoff);
            candidate = divided_length + extra;
          }
        } else if (use_iw) {
          candidate = indirect_iw_length_bounded(
              clip_prelim, tree, ds, above, below, base_iw, iw_delta,
              best_candidate);
        } else {
          int cutoff = (best_candidate < HUGE_VAL)
              ? static_cast<int>(best_candidate - divided_length + 1)
              : INT_MAX;
          // Exact directional cost: the edge set above `below` (= node_d) is
          // edge_set_buf[below], replacing the union-of-finals approximation.
          int extra = fitch_indirect_length_cached(
              clip_prelim, &edge_set_buf[static_cast<size_t>(below) * tw],
              ds, cutoff);
          candidate = divided_length + extra;
        }
        ++n_evaluated;
        if (candidate < best_candidate) {
          best_candidate = candidate;
          best_above = above;
          best_below = below;
          best_reroot_parent = -1;
          best_reroot_child = -1;
        }
      }

      // TBR candidates (rerooting) — with vroot cache (optimization #4)
      if (clip_node >= tree.n_tip) {
        compute_from_above(tree, ds, clip_node, from_above);
        collect_subtree_edges(tree, clip_node, sub_edges);

        // Precompute vroot for all main edges (optimization #4).  Pure EW uses
        // the exact directional edge set (vroot = edge_set_buf[below]); NA / IW
        // keep the union-of-finals form their cached scorers expect.
        int n_main = static_cast<int>(main_edges.size());
        if (ew_directional) {
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
            int ei = 0;
            while (ei < n_main) {
              int b_ei[4];
              int b_n = 0;
              while (ei < n_main && b_n < 4) {
                auto& [ab, bl] = main_edges[ei];
                bool skip = (ab == nz && bl == ns)
                    || (sector_mask && !(*sector_mask)[bl])
                    || (constrained && regraft_violates_constraint(bl, *cd))
                    || (!collapsed.empty() && collapsed[bl]);
                if (!skip) b_ei[b_n++] = ei;
                ++ei;
              }
              if (b_n == 0) break;

              int cutoff_b = (best_candidate < HUGE_VAL)
                  ? static_cast<int>(best_candidate - divided_length + 1)
                  : INT_MAX;
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
                }
              }
            }
          } else {
            // === Scalar path (IW, or ratchet with upweight_mask) ===
            for (int ei = 0; ei < n_main; ++ei) {
              // Prefetch vroot data for a future iteration.
              // At 180 tips vroot_cache is ~140 KB (L2); prefetch hides
              // the ~10-cycle L2 latency. No-op on small trees (L1-resident).
              if (ei + 2 < n_main) {
#if defined(__GNUC__) || defined(__clang__)
                __builtin_prefetch(
                    &vroot_cache[static_cast<size_t>(ei + 2) * tree.total_words],
                    0, 0);
#elif defined(_MSC_VER) && defined(TS_SIMD_SSE2)
                _mm_prefetch(reinterpret_cast<const char*>(
                    &vroot_cache[static_cast<size_t>(ei + 2) * tree.total_words]),
                    _MM_HINT_T0);
#endif
              }
              auto& [above, below] = main_edges[ei];
              if (above == nz && below == ns) continue;
              if (sector_mask && !(*sector_mask)[below]) continue;
              if (constrained && regraft_violates_constraint(below, *cd))
                continue;
              // Collapsed-region regraft merging (same as SPR loop above).
              if (!collapsed.empty() && collapsed[below]) {
                continue;
              }
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
                  int cutoff = (best_candidate < HUGE_VAL)
                      ? static_cast<int>(best_candidate - divided_length + 1)
                      : INT_MAX;
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
                int cutoff = (best_candidate < HUGE_VAL)
                    ? static_cast<int>(best_candidate - divided_length + 1)
                    : INT_MAX;
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
              }
            }
          }
        }
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
      if (std::getenv("TS_REVERT_CHECK")) {
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
          if (collect_pool) collect_pool->add(tree, actual);
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
          if (collect_pool) collect_pool->add(tree, actual);
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
        }
        // Optimization #6: don't reshuffle after acceptance — the topology
        // changed near this clip, so re-trying the same ordering focuses
        // on the productive region.
        need_shuffle = false;
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
  if (!collapsed.empty()) compute_collapsed_flags(tree, ds, collapsed);
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

  return TBRResult{best_score, n_accepted, n_evaluated, n_zero_skipped,
                   converged, std::move(diag_records)};
}

} // namespace ts
