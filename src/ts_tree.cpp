#include "ts_tree.h"
#include <algorithm>

namespace ts {

void TreeState::init_from_edge(
    const int* edge_parent, const int* edge_child,
    int n_edge, const DataSet& ds)
{
  n_tip = (n_edge / 2) + 1;
  n_internal = n_tip - 1;
  n_node = n_tip + n_internal;
  total_words = ds.total_words;
  n_blocks = ds.n_blocks;

  parent.assign(n_node, -1);
  left.assign(n_internal, -1);
  right.assign(n_internal, -1);

  // Edge matrix is 1-based; convert to 0-based.
  for (int i = 0; i < n_edge; ++i) {
    int p = edge_parent[i] - 1;
    int c = edge_child[i] - 1;
    parent[c] = p;
    int pi = p - n_tip;
    // Assign first child to left, second to right — matches tree_to_edge
    // output order, so R↔C++ round-trips preserve left/right assignment.
    if (left[pi] == -1) {
      left[pi] = c;
    } else {
      right[pi] = c;
    }
  }
  parent[n_tip] = n_tip;  // root is its own parent

  // Allocate state sets
  size_t state_size = static_cast<size_t>(n_node) * total_words;
  prelim.assign(state_size, 0ULL);
  final_.assign(state_size, 0ULL);
  down2.assign(state_size, 0ULL);
  subtree_actives.assign(state_size, 0ULL);

  // Local cost: one uint64_t per block per node (only internal meaningful)
  local_cost.assign(static_cast<size_t>(n_node) * n_blocks, 0ULL);

  load_tip_states(ds);
  build_postorder();
}

void TreeState::load_tip_states(const DataSet& ds) {
  for (int tip = 0; tip < n_tip; ++tip) {
    size_t base = static_cast<size_t>(tip) * total_words;
    for (int w = 0; w < total_words; ++w) {
      prelim[base + w] = ds.tip_states[base + w];
      // Tips: final = preliminary (observed states, unchanged by uppass)
      final_[base + w] = ds.tip_states[base + w];
    }
  }
  // Initialise tip subtree_actives: applicable states only (NA word = 0)
  for (int tip = 0; tip < n_tip; ++tip) {
    size_t base = static_cast<size_t>(tip) * total_words;
    for (int b = 0; b < ds.n_blocks; ++b) {
      int offset = ds.block_word_offset[b];
      if (ds.blocks[b].has_inapplicable) {
        subtree_actives[base + offset] = 0;  // state 0 (NA) always 0
        for (int s = 1; s < ds.blocks[b].n_states; ++s) {
          subtree_actives[base + offset + s] = ds.tip_states[base + offset + s];
        }
      }
    }
  }
}

void TreeState::build_postorder() {
  // Build postorder by recursive-style DFS from root.
  // Only includes nodes reachable from root — critical after SPR clip,
  // when the detached subtree must not appear in the traversal.
  postorder.clear();
  postorder.reserve(n_internal);

  // Two-pass DFS: first build preorder, then reverse to get postorder.
  // (Simpler and less error-prone than sentinel-based single-pass.)
  std::vector<int> preorder;
  preorder.reserve(n_internal);

  std::vector<int> stack;
  stack.push_back(n_tip);  // root

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    if (node < n_tip) continue;  // tip

    preorder.push_back(node);
    int ni = node - n_tip;
    // Push right first so left is processed first (standard preorder)
    stack.push_back(right[ni]);
    stack.push_back(left[ni]);
  }

  // Reverse preorder = postorder
  postorder.assign(preorder.rbegin(), preorder.rend());
}

// ---- NNI ----

TreeState::NNIUndo TreeState::nni_apply(int c, int which_child) {
  int p = parent[c];
  int ci = c - n_tip;
  int pi = p - n_tip;

  bool c_is_left = (left[pi] == c);
  int& sib_slot = c_is_left ? right[pi] : left[pi];
  int sib = sib_slot;

  int& c_child_slot = (which_child == 0) ? left[ci] : right[ci];
  int c_child = c_child_slot;

  sib_slot = c_child;
  c_child_slot = sib;
  parent[c_child] = p;
  parent[sib] = c;

  return NNIUndo{c, which_child};
}

void TreeState::nni_undo(const NNIUndo& undo) {
  nni_apply(undo.c, undo.which_child);
}

std::vector<int> TreeState::nni_edges() const {
  std::vector<int> edges;
  edges.reserve(n_internal);
  int root = n_tip;
  for (int c = n_tip; c < n_node; ++c) {
    if (c != root) {
      edges.push_back(c);
    }
  }
  return edges;
}

// ---- Helpers ----

int TreeState::sibling(int node) const {
  int p = parent[node];
  int pi = p - n_tip;
  return (left[pi] == node) ? right[pi] : left[pi];
}

void TreeState::save_node_state(int node) {
  NodeSnapshot snap;
  snap.node = node;

  size_t state_base = static_cast<size_t>(node) * total_words;
  snap.prelim_words.assign(prelim.begin() + state_base,
                           prelim.begin() + state_base + total_words);
  snap.final_words.assign(final_.begin() + state_base,
                          final_.begin() + state_base + total_words);

  size_t cost_base = static_cast<size_t>(node) * n_blocks;
  snap.local_costs.assign(local_cost.begin() + cost_base,
                          local_cost.begin() + cost_base + n_blocks);

  clip_undo_stack.push_back(std::move(snap));
}

void TreeState::restore_saved_states() {
  while (!clip_undo_stack.empty()) {
    const NodeSnapshot& snap = clip_undo_stack.back();
    int node = snap.node;

    size_t state_base = static_cast<size_t>(node) * total_words;
    std::copy(snap.prelim_words.begin(), snap.prelim_words.end(),
              prelim.begin() + state_base);
    std::copy(snap.final_words.begin(), snap.final_words.end(),
              final_.begin() + state_base);

    size_t cost_base = static_cast<size_t>(node) * n_blocks;
    std::copy(snap.local_costs.begin(), snap.local_costs.end(),
              local_cost.begin() + cost_base);

    clip_undo_stack.pop_back();
  }
}

void TreeState::reset_states(const DataSet& ds) {
  std::fill(prelim.begin(), prelim.end(), 0ULL);
  std::fill(final_.begin(), final_.end(), 0ULL);
  std::fill(down2.begin(), down2.end(), 0ULL);
  std::fill(subtree_actives.begin(), subtree_actives.end(), 0ULL);
  std::fill(local_cost.begin(), local_cost.end(), 0ULL);
  load_tip_states(ds);
}

// ---- SPR topology ----

void TreeState::spr_clip(int clip_node) {
  int nx = parent[clip_node];       // parent of clip node
  int nz = parent[nx];              // grandparent
  int ns = sibling(clip_node);      // sibling
  int nxi = nx - n_tip;

  // Save topology for undo
  clip_state.clip_node = clip_node;
  clip_state.clip_parent = nx;
  clip_state.clip_sibling = ns;
  clip_state.clip_grandpar = nz;
  clip_state.clip_is_left = (left[nxi] == clip_node);
  if (nz >= n_tip) {
    int nzi = nz - n_tip;
    clip_state.nx_is_left = (left[nzi] == nx);
  }

  // Detach: connect sibling directly to grandparent
  parent[ns] = nz;
  if (nz >= n_tip && nz != nx) {
    int nzi = nz - n_tip;
    if (left[nzi] == nx) {
      left[nzi] = ns;
    } else {
      right[nzi] = ns;
    }
  } else if (nz == nx) {
    // nx was root — ns becomes root's child... actually ns becomes root
    // This case: clip_node is child of root. After clip, the tree should
    // be re-rooted. For now, handle by making ns the new effective root
    // connection. Actually, if nx == nz, nx is root. We need to handle
    // this: ns takes over as root's child on the clip side.
    // In practice, for SPR we typically don't clip children of root,
    // but let's handle it:
    parent[ns] = nx;  // ns stays under root
    // nx's child slots: one was clip_node, the other was ns.
    // After clip, nx has ns as one child. We need to figure out the other.
    // Actually when clipping a child of root, the tree becomes: root has
    // one child (ns) and the clipped subtree. The "divided tree" is just
    // the main tree minus the clipped subtree, with root -> ns.
    // For simplicity, set nx's clip-side child to ns (already the case)
    // and leave the other slot as-is (it will be overwritten on regraft).
  }

  // Clear the undo stack (will be populated during incremental pass)
  clip_undo_stack.clear();
}

void TreeState::spr_unclip() {
  // Restore state sets first
  restore_saved_states();

  // Restore topology
  int nm = clip_state.clip_node;
  int nx = clip_state.clip_parent;
  int ns = clip_state.clip_sibling;
  int nz = clip_state.clip_grandpar;
  int nxi = nx - n_tip;

  // Reconnect nx between nz and ns
  if (nz >= n_tip && nz != nx) {
    int nzi = nz - n_tip;
    if (clip_state.nx_is_left) {
      left[nzi] = nx;
    } else {
      right[nzi] = nx;
    }
  }
  parent[nx] = nz;

  // Restore nx's children
  if (clip_state.clip_is_left) {
    left[nxi] = nm;
    right[nxi] = ns;
  } else {
    left[nxi] = ns;
    right[nxi] = nm;
  }
  parent[nm] = nx;
  parent[ns] = nx;
}

void TreeState::spr_regraft(int above, int below) {
  // Insert nx (clip parent) between above and below.
  int nx = clip_state.clip_parent;
  int nm = clip_state.clip_node;
  int nxi = nx - n_tip;

  // above -> below becomes above -> nx -> below, with nm as nx's other child.
  if (above >= n_tip) {
    int ai = above - n_tip;
    if (left[ai] == below) {
      left[ai] = nx;
    } else {
      right[ai] = nx;
    }
  }
  parent[nx] = above;

  // nx's children: nm (the clipped subtree) and below
  left[nxi] = nm;
  right[nxi] = below;
  parent[nm] = nx;
  parent[below] = nx;
}

void TreeState::spr_unregraft(int above, int below) {
  // Remove nx from between above and below.
  int nx = clip_state.clip_parent;

  // Reconnect above -> below directly
  if (above >= n_tip) {
    int ai = above - n_tip;
    if (left[ai] == nx) {
      left[ai] = below;
    } else {
      right[ai] = below;
    }
  }
  parent[below] = above;
}

} // namespace ts
