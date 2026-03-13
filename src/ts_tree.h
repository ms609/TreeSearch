#ifndef TS_TREE_H
#define TS_TREE_H

// TreeState: flat-array tree topology with per-node bit-packed state sets.
//
// Topology uses the same convention as build_postorder.h:
//   Tips:     0 .. n_tip-1
//   Internal: n_tip .. 2*n_tip-2
//   Root:     n_tip (parent[root] = root)
//
// left[] and right[] are indexed by (node - n_tip), i.e. only internal nodes.
//
// Per-node state arrays:
//   prelim:      preliminary (downpass) state sets
//   final_:      final (uppass) state sets — needed for indirect calculation
//   local_cost:  per-node per-block needs_union mask (which chars incurred a
//                step at this node). Used for length delta in incremental pass.
//
// All state arrays are node-major: all words for one node are contiguous.
//   prelim[node * total_words + block_word_offset[b] + s]
//   local_cost[node * n_blocks + b]

#include "ts_data.h"
#include <vector>
#include <cstdint>

namespace ts {

struct TreeState {
  int n_tip;
  int n_internal;      // = n_tip - 1
  int n_node;          // = 2 * n_tip - 1
  int total_words;     // sum of n_states across all blocks
  int n_blocks;        // number of character blocks

  // Topology (indexed by node id)
  std::vector<int> parent;      // [n_node] — parent[root] = root
  std::vector<int> left;        // [n_internal] — left child
  std::vector<int> right;       // [n_internal] — right child

  // --- Per-node state sets (flattened, node-major) ---

  // Preliminary (downpass) state sets.
  //   prelim[node * total_words + word_offset]
  std::vector<uint64_t> prelim;

  // Final (uppass) state sets. Computed from prelim + ancestor's final.
  //   final_[node * total_words + word_offset]
  std::vector<uint64_t> final_;

  // Per-node per-block local cost: the needs_union mask from the downpass.
  //   local_cost[node * n_blocks + b]
  // Only meaningful for internal nodes; tip entries are unused.
  std::vector<uint64_t> local_cost;

  // Postorder traversal sequence (internal nodes only, leaves to root)
  std::vector<int> postorder;

  // --- Undo buffer for incremental two-pass (Shortcut C) ---

  // Snapshot of one node's state during incremental pass.
  struct NodeSnapshot {
    int node;
    std::vector<uint64_t> prelim_words;   // total_words entries
    std::vector<uint64_t> final_words;    // total_words entries
    std::vector<uint64_t> local_costs;    // n_blocks entries
  };

  // Stack of snapshots saved during clip; restored on unclip.
  std::vector<NodeSnapshot> clip_undo_stack;

  // --- SPR clip state ---
  // Saved topology around the clip point for undo.
  struct SPRClipState {
    int clip_node;      // Nm: the node being clipped (subtree root)
    int clip_parent;    // Nx: parent of Nm before clip
    int clip_sibling;   // Ns: sibling of Nm before clip
    int clip_grandpar;  // Nz: grandparent (= parent of Nx)
    bool clip_is_left;  // was Nm the left child of Nx?
    bool nx_is_left;    // was Nx the left child of Nz?
    int divided_length; // length of divided tree after clip
  };
  SPRClipState clip_state;

  // ---- Initialization ----

  void init_from_edge(const int* edge_parent, const int* edge_child,
                      int n_edge, const DataSet& ds);
  void load_tip_states(const DataSet& ds);
  void build_postorder();

  // ---- NNI operations ----

  struct NNIUndo {
    int c;
    int which_child;
  };

  NNIUndo nni_apply(int c, int which_child);
  void nni_undo(const NNIUndo& undo);
  std::vector<int> nni_edges() const;

  // ---- SPR operations ----

  // Clip (prune) the subtree rooted at clip_node from the tree.
  // Connects clip_node's sibling directly to grandparent.
  // Does NOT update state sets — caller must run incremental pass.
  void spr_clip(int clip_node);

  // Undo a clip: restore the original topology.
  // Also restores state sets from clip_undo_stack.
  void spr_unclip();

  // Regraft the clipped subtree onto edge (above, below) where `above`
  // is the parent and `below` is the child in the target edge.
  // The clip's parent node (Nx) is re-inserted between above and below.
  void spr_regraft(int above, int below);

  // Undo a regraft: detach the re-inserted node and restore the edge.
  // This does NOT restore the original clip — use spr_unclip() for that.
  void spr_unregraft(int above, int below);

  // ---- Helpers ----

  // Save a node's state to the undo stack.
  void save_node_state(int node);

  // Restore all saved states from the undo stack (in reverse order).
  void restore_saved_states();

  // Get sibling of node under its parent.
  int sibling(int node) const;

  // Reset all state arrays (prelim, final_, local_cost) from current topology.
  // Reloads tip states and clears internal node states.
  void reset_states(const DataSet& ds);
};

} // namespace ts

#endif // TS_TREE_H
