#include "ts_tbr.h"
#include "ts_fitch.h"
#include <algorithm>
#include <random>
#include <vector>
#include <climits>
#include <cmath>

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>

namespace ts {

// --- Helpers (file-local) ---

static int full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return fitch_score(tree, ds);
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


// --- Main TBR search ---

TBRResult tbr_search(TreeState& tree, const DataSet& ds,
                     const TBRParams& params) {
  double best_score = full_rescore(tree, ds);
  int n_accepted = 0;
  int n_evaluated = 0;
  int hits = 1;

  std::mt19937 rng(std::random_device{}());

  // Candidate clip nodes: all non-root nodes
  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;
    clip_candidates.push_back(node);
  }

  std::vector<std::pair<int,int>> main_edges;
  std::vector<std::pair<int,int>> sub_edges;

  // Temporary buffers
  std::vector<uint64_t> from_above(
      static_cast<size_t>(tree.n_node) * tree.total_words, 0);
  std::vector<uint64_t> virtual_prelim(tree.total_words);

  TopoSnapshot snap;
  bool keep_going = true;

  while (keep_going) {
    keep_going = false;
    std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

    for (int clip_node : clip_candidates) {
      if (tree.parent[clip_node] == tree.n_tip) continue;

      // --- Phase 1: Clip + indirect evaluation ---
      tree.spr_clip(clip_node);
      tree.build_postorder();

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;
      // Start at nz (grandparent), not ns — ns may be a tip.
      int nx = tree.clip_state.clip_parent;
      int delta = fitch_incremental_downpass(tree, ds, nz);
      fitch_incremental_uppass(tree, ds, nz);

      // nx was removed from the tree. Its local_cost (untouched by the
      // incremental pass, since nx is disconnected) must be subtracted.
      int nx_cost = 0;
      for (int b = 0; b < ds.n_blocks; ++b) {
        nx_cost += ds.blocks[b].weight * popcount64(
            tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b]);
      }
      double divided_length = best_score + delta - nx_cost;

      collect_main_edges(tree, main_edges);

      // Find best (reroot, regraft) combination
      int best_extra = INT_MAX;
      int best_above = -1, best_below = -1;
      int best_reroot_parent = -1, best_reroot_child = -1;

      // SPR candidates
      size_t clip_base = static_cast<size_t>(clip_node) * tree.total_words;
      const uint64_t* clip_prelim = &tree.prelim[clip_base];

      for (auto& [above, below] : main_edges) {
        if (above == nz && below == ns) continue;
        int extra = fitch_indirect_length(clip_prelim, tree, ds,
                                          above, below);
        ++n_evaluated;
        if (extra < best_extra) {
          best_extra = extra;
          best_above = above;
          best_below = below;
          best_reroot_parent = -1;
          best_reroot_child = -1;
        }
      }

      // TBR candidates (rerooting)
      if (clip_node >= tree.n_tip) {
        compute_from_above(tree, ds, clip_node, from_above);
        collect_subtree_edges(tree, clip_node, sub_edges);

        for (auto& [sp, sc] : sub_edges) {
          if (sp == clip_node) continue;

          fitch_join_states(
              &from_above[static_cast<size_t>(sc) * tree.total_words],
              &tree.prelim[static_cast<size_t>(sc) * tree.total_words],
              virtual_prelim.data(), ds);

          for (auto& [above, below] : main_edges) {
            if (above == nz && below == ns) continue;
            int extra = fitch_indirect_length(virtual_prelim.data(), tree,
                                              ds, above, below);
            ++n_evaluated;
            if (extra < best_extra) {
              best_extra = extra;
              best_above = above;
              best_below = below;
              best_reroot_parent = sp;
              best_reroot_child = sc;
            }
          }
        }
      }

      // --- Phase 2: Restore original tree, verify best candidate ---
      tree.spr_unclip();
      tree.build_postorder();

      double candidate_score = (best_extra == INT_MAX)
                              ? HUGE_VAL : divided_length + best_extra;
      bool dominated = (candidate_score > best_score) ||
                        (candidate_score == best_score && !params.accept_equal);

      bool accepted = false;

      if (!dominated && best_above >= 0) {
        save_topology(tree, snap);

        bool ok = apply_tbr_move(tree, clip_node,
                                  best_reroot_parent, best_reroot_child,
                                  best_above, best_below);

        if (!ok || !validate_topology(tree)) {
          restore_topology(tree, snap);
          tree.build_postorder();
          continue;
        }

        tree.build_postorder();
        double actual = full_rescore(tree, ds);

        if (actual < best_score) {
          best_score = actual;
          ++n_accepted;
          hits = 1;
          accepted = true;
          keep_going = true;
        } else if (actual == best_score && params.accept_equal
                   && hits <= params.max_hits) {
          ++hits;
          ++n_accepted;
          accepted = true;
          keep_going = true;
        }

        if (!accepted) {
          restore_topology(tree, snap);
          tree.build_postorder();
          full_rescore(tree, ds);
        }
      }

      if (keep_going) {
        if (params.max_accepted_changes > 0
            && n_accepted >= params.max_accepted_changes) {
          keep_going = false;
        }
        break;
      }

      R_CheckUserInterrupt();
    }

    if (params.max_accepted_changes > 0
        && n_accepted >= params.max_accepted_changes) {
      break;
    }
  }

  best_score = full_rescore(tree, ds);

  bool converged = !(params.max_accepted_changes > 0
                     && n_accepted >= params.max_accepted_changes);

  return TBRResult{best_score, n_accepted, n_evaluated, converged};
}

} // namespace ts
