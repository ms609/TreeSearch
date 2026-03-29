#include "ts_tbr.h"
#include "ts_fitch.h"
#include "ts_collapsed.h"
#include "ts_rng.h"
#include "ts_tabu.h"
#include "ts_splits.h"
#include <algorithm>
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

// --- Helpers (file-local) ---

static double full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

// Collect (parent, child) edge pairs reachable from root of main tree.
static void collect_main_edges(
    const TreeState& tree,
    std::vector<std::pair<int,int>>& edges)
{
  edges.clear();
  std::vector<int> stack;
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

  std::vector<int> stack;
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

  std::vector<int> preorder;
  {
    std::vector<int> stack;
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

static bool validate_topology(const TreeState& tree) {
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


// --- Main TBR search ---

TBRResult tbr_search(TreeState& tree, const DataSet& ds,
                     const TBRParams& params,
                     ConstraintData* cd,
                     const std::vector<bool>* sector_mask,
                     TreePool* collect_pool,
                     std::function<bool()> check_timeout) {
  double best_score = full_rescore(tree, ds);

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
  int n_targeted_calls = 0;
  int n_targeted_improved = 0;
  int hits = 1;
  const bool use_iw = std::isfinite(ds.concavity);
  // Floating-point tolerance for score equality
  const double eps = use_iw ? 1e-10 : 0.0;

  // Check if any block has inapplicable characters (for state snapshot)
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

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

  TopoSnapshot snap;
  bool keep_going = true;
  bool states_valid = true;  // track whether state arrays are current
  bool need_shuffle = true;  // optimization #6: defer reshuffle
  // Poll timeout every n_tip clips to avoid overhead on small trees
  // while ensuring responsiveness on large ones.
  const int timeout_interval = std::max(tree.n_tip, 50);
  int clips_since_timeout_check = 0;
  bool timed_out = false;

  while (keep_going && !timed_out) {
    keep_going = false;

    // Optimization #7: save state snapshot once per pass, not per candidate.
    // After a rejected move, state_snap.restore() returns the tree to exactly
    // the state saved here. The per-candidate save was redundant: consecutive
    // rejections all restore to the same state. Re-saving only happens when
    // the while loop restarts after an accepted move.
    save_topology(tree, snap);
    state_snap.save(tree);
    states_valid = true;

    // Optimization #6: only reshuffle when the previous pass found no
    // improvement. After an accepted move, retry with the same ordering
    // (the topology changed, so previously-failing clips may now succeed).
    if (need_shuffle) {
      std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);
    }
    need_shuffle = true;  // default: reshuffle next time (unless we accept)

    // Recompute subtree sizes for smaller-subtree filtering
    compute_subtree_sizes(tree, subtree_sizes);

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
          int extra = use_flat
              ? fitch_indirect_bounded_flat(clip_prelim, tree, ds,
                    above, below, cutoff)
              : fitch_indirect_length_bounded(clip_prelim, tree, ds,
                    above, below, cutoff);
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

        // Precompute vroot for all main edges (optimization #4)
        precompute_vroot_cache(tree, main_edges, vroot_cache);
        int n_main = static_cast<int>(main_edges.size());

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
        std::unordered_set<uint64_t> seen_vp_hashes;
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
          if (!seen_vp_hashes.insert(vp_hash).second) {
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

      bool dominated = (best_candidate > best_score + eps) ||
                        (best_candidate > best_score - eps && !params.accept_equal);

      bool accepted = false;

      if (!dominated && best_above >= 0) {
        // Topology and state snapshot already saved at the top of the
        // while loop (optimization #7). No per-candidate save needed.

        bool ok = apply_tbr_move(tree, clip_node,
                                  best_reroot_parent, best_reroot_child,
                                  best_above, best_below);

        if (!ok || !validate_topology(tree)) {
          restore_topology(tree, snap);
          state_snap.restore(tree);
          continue;
        }

        tree.build_postorder_prealloc(work_stack);
        double actual = full_rescore(tree, ds);

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
          ++n_accepted;
          hits = 1;
          accepted = true;
          keep_going = true;
          states_valid = true;
          if (constrained) {
            compute_dfs_timestamps(tree, *cd);
          }
          if (collect_pool) collect_pool->add(tree, actual);

          // Targeted sector refinement: after a strict improvement,
          // run TBR within the just-moved clip subtree to refine its
          // internal arrangement for the new context.
          if (params.targeted_sector
              && !collect_pool
              && !sector_mask
              && subtree_sizes[clip_node] >= params.targeted_min_size)
          {
            std::vector<bool> target_mask(tree.n_node, false);
            {
              std::vector<int> stk;
              stk.push_back(clip_node);
              while (!stk.empty()) {
                int nd = stk.back();
                stk.pop_back();
                target_mask[nd] = true;
                if (nd >= tree.n_tip) {
                  int ni = nd - tree.n_tip;
                  stk.push_back(tree.left[ni]);
                  stk.push_back(tree.right[ni]);
                }
              }
            }

            TBRParams inner_params;
            inner_params.max_hits = 1;
            inner_params.accept_equal = false;
            inner_params.targeted_sector = false;

            double pre_target = best_score;
            TBRResult inner = tbr_search(tree, ds, inner_params, cd,
                                         &target_mask);

            ++n_targeted_calls;
            if (inner.best_score < pre_target - eps) {
              best_score = inner.best_score;
              ++n_targeted_improved;
            }
          }
        } else if (std::fabs(actual - best_score) <= eps
                   && params.accept_equal
                   && hits <= params.max_hits) {
          // Equal-score move: reject if tabu
          if (tabu.active() && tabu.contains(tree_hash)) {
            // Topology already visited — restore and skip
            restore_topology(tree, snap);
            state_snap.restore(tree);
            continue;
          }
          if (tabu.active()) tabu.insert(tree_hash);
          ++hits;
          ++n_accepted;
          accepted = true;
          keep_going = true;
          states_valid = true;
          if (constrained) {
            compute_dfs_timestamps(tree, *cd);
          }
          if (collect_pool) collect_pool->add(tree, actual);
        }

        if (!accepted) {
          // Optimization #3: restore topology + states without full_rescore
          restore_topology(tree, snap);
          state_snap.restore(tree);
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

      if (ts::check_interrupt()) break;
      ++clips_since_timeout_check;
      if (check_timeout && clips_since_timeout_check >= timeout_interval) {
        clips_since_timeout_check = 0;
        if (check_timeout()) { timed_out = true; break; }
      }
    }

    if (params.max_accepted_changes > 0
        && n_accepted >= params.max_accepted_changes) {
      break;
    }
  }

  tree.prealloc_undo = nullptr;

  // Ensure state arrays match the final tree and return an authoritative score.
  // states_valid is always true at this point (it is initialized true and only
  // ever set back to true on acceptance), so both branches reduce to the same
  // call.  Simplified to a single assignment.
  best_score = full_rescore(tree, ds);

  bool converged = !(params.max_accepted_changes > 0
                     && n_accepted >= params.max_accepted_changes);

  return TBRResult{best_score, n_accepted, n_evaluated, n_zero_skipped,
                   converged, n_targeted_calls, n_targeted_improved};
}

} // namespace ts
