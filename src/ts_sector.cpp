#include "ts_sector.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_wagner.h"
#include "ts_rng.h"
#include "ts_splits.h"
#include "ts_pool.h"

#include <algorithm>
#include <random>
#include <vector>
#include <cstring>
#include <R.h>

namespace ts {

// ---- From-above state computation for exact HTU ----
//
// Compute from_above[sector_root]: the Fitch state-set that the rest of
// the tree sends down to the sector boundary, EXCLUDING the sector's own
// contribution. Using this as the HTU state-set (instead of final_[parent])
// makes the sector score a better predictor of the full-tree impact.
//
// Algorithm: walk from root down to sector_root, computing fitch_join
// at each step. O(depth × total_words) — negligible.

static void compute_from_above_for_sector(
    const TreeState& tree, const DataSet& ds,
    int sector_root,
    std::vector<uint64_t>& from_above_out) {
  int tw = tree.total_words;
  from_above_out.resize(tw);

  // 1. Find path from sector_root to root (walk up via parent)
  std::vector<int> path;
  for (int cur = sector_root; ; cur = tree.parent[cur]) {
    path.push_back(cur);
    if (cur == tree.n_tip) break; // reached root
  }
  // path = [sector_root, ..., root]
  // Reverse to get [root, ..., sector_root]
  std::reverse(path.begin(), path.end());

  // 2. Seed: from_above[root] = all states (fully ambiguous)
  // All bits set for each state word within each block.
  std::vector<uint64_t> from_above_cur(tw);
  for (int b = 0; b < ds.n_blocks; ++b) {
    int off = ds.block_word_offset[b];
    uint64_t mask = ds.blocks[b].active_mask;
    for (int s = 0; s < ds.blocks[b].n_states; ++s) {
      from_above_cur[off + s] = mask;
    }
  }

  // 3. Walk down the path, computing from_above at each child step
  for (size_t i = 0; i + 1 < path.size(); ++i) {
    int node = path[i];
    int next = path[i + 1]; // child on the path

    // Find sibling of `next` under `node`
    int ni = node - tree.n_tip;
    int sib = (tree.left[ni] == next) ? tree.right[ni] : tree.left[ni];

    // from_above[next] = fitch_join(from_above[node], prelim[sib])
    // fitch_join: per-block, compute intersection; where empty, use union.
    const uint64_t* sib_prelim =
        &tree.prelim[static_cast<size_t>(sib) * tw];

    std::vector<uint64_t> new_from_above(tw);
    for (int b = 0; b < ds.n_blocks; ++b) {
      int off = ds.block_word_offset[b];
      int ns = ds.blocks[b].n_states;
      uint64_t any_isect = 0;
      for (int s = 0; s < ns; ++s) {
        any_isect |= (from_above_cur[off + s] & sib_prelim[off + s]);
      }
      uint64_t no_isect = ~any_isect & ds.blocks[b].active_mask;
      for (int s = 0; s < ns; ++s) {
        uint64_t isect = from_above_cur[off + s] & sib_prelim[off + s];
        uint64_t uni = from_above_cur[off + s] | sib_prelim[off + s];
        new_from_above[off + s] = (isect & any_isect) | (uni & no_isect);
      }
    }
    from_above_cur = std::move(new_from_above);
  }

  std::memcpy(from_above_out.data(), from_above_cur.data(),
              tw * sizeof(uint64_t));
}

// ---- Conflict-guided sector selection ----
//
// For each internal node (except root), compute how "conflicted" it is:
// the fraction of best-score pool trees that do NOT contain this node's split.
// Returns per-node values in [0, 1]: 0 = unanimous, 1 = absent from all pool
// trees. Tips and root get 0. When sft has <2 trees, returns all zeros.
//
// Also propagates max-descendant conflict upward: sector_conflict[node] =
// max conflict score among all nodes in node's subtree. This gives each
// eligible sector root a score reflecting the most uncertain region it contains.

static void compute_node_conflict(
    const TreeState& tree,
    const SplitFrequencyTable& sft,
    std::vector<double>& node_conflict,
    std::vector<double>& sector_conflict) {
  int nn = tree.n_node;
  node_conflict.assign(nn, 0.0);
  sector_conflict.assign(nn, 0.0);
  if (sft.n_trees < 2) return;

  int n_tip = tree.n_tip;
  int wps = (n_tip + 63) / 64;
  int trailing = n_tip % 64;
  uint64_t trail_mask = (trailing != 0) ? ((1ULL << trailing) - 1) : ~0ULL;

  // Build tip membership bitsets (same as compute_splits, but we also need
  // the node mapping so we do it in-place rather than calling compute_splits).
  size_t total = static_cast<size_t>(nn) * wps;
  std::vector<uint64_t> tip_bits(total, 0);

  for (int t = 0; t < n_tip; ++t) {
    tip_bits[static_cast<size_t>(t) * wps + t / 64] = 1ULL << (t % 64);
  }

  int root = n_tip;
  int root_right = tree.right[0];

  // Temporary buffer for canonicalized split
  std::vector<uint64_t> canon(wps);

  for (int node : tree.postorder) {
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* dst = &tip_bits[static_cast<size_t>(node) * wps];
    const uint64_t* lbits = &tip_bits[static_cast<size_t>(lc) * wps];
    const uint64_t* rbits = &tip_bits[static_cast<size_t>(rc) * wps];
    for (int w = 0; w < wps; ++w) {
      dst[w] = lbits[w] | rbits[w];
    }

    // Skip root and root_right (same exclusions as compute_splits)
    if (node == root || node == root_right) continue;

    // Check for trivial splits
    int count = 0;
    for (int w = 0; w < wps; ++w) {
      count += ts::popcount64(dst[w]);
    }
    if (count <= 1 || count >= n_tip - 1) continue;

    // Canonicalize: ensure bit 0 is clear
    bool flip = (dst[0] & 1ULL) != 0;
    for (int w = 0; w < wps; ++w) {
      canon[w] = flip ? ~dst[w] : dst[w];
    }
    canon[wps - 1] &= trail_mask;

    // Look up frequency in pool
    uint64_t sh = hash_single_split(canon.data(), wps);
    auto it = sft.freq.find(sh);
    int freq = (it != sft.freq.end()) ? it->second : 0;
    node_conflict[node] = 1.0 - static_cast<double>(freq) / sft.n_trees;
  }

  // Propagate max-descendant conflict upward (postorder)
  for (int t = 0; t < n_tip; ++t) sector_conflict[t] = 0.0;
  for (int node : tree.postorder) {
    int ni = node - n_tip;
    sector_conflict[node] = std::max({
        node_conflict[node],
        sector_conflict[tree.left[ni]],
        sector_conflict[tree.right[ni]]
    });
  }
}

// ---- Clade topology snapshot ----
// Saves only the internal nodes within a clade (sector) for fast undo.
// Much cheaper than copying the full tree's left/right/parent vectors.

struct CladeSnapshot {
  std::vector<int> internals;  // full-tree internal node indices
  std::vector<int> left;       // saved left[ni] for each
  std::vector<int> right;      // saved right[ni] for each
  // Parent links for all children (tips and internals) are restored
  // from left/right during restore.
};

static void save_clade(const TreeState& tree, int clade_root,
                        CladeSnapshot& snap) {
  snap.internals.clear();
  // DFS to collect internal nodes in the clade
  std::vector<int> stk;
  stk.push_back(clade_root);
  while (!stk.empty()) {
    int nd = stk.back();
    stk.pop_back();
    if (nd < tree.n_tip) continue;
    snap.internals.push_back(nd);
    int ni = nd - tree.n_tip;
    stk.push_back(tree.left[ni]);
    stk.push_back(tree.right[ni]);
  }
  // Save left/right for each internal node
  snap.left.resize(snap.internals.size());
  snap.right.resize(snap.internals.size());
  for (size_t i = 0; i < snap.internals.size(); ++i) {
    int ni = snap.internals[i] - tree.n_tip;
    snap.left[i] = tree.left[ni];
    snap.right[i] = tree.right[ni];
  }
}

static void restore_clade(TreeState& tree, const CladeSnapshot& snap) {
  for (size_t i = 0; i < snap.internals.size(); ++i) {
    int ni = snap.internals[i] - tree.n_tip;
    tree.left[ni] = snap.left[i];
    tree.right[ni] = snap.right[i];
    tree.parent[snap.left[i]] = snap.internals[i];
    tree.parent[snap.right[i]] = snap.internals[i];
  }
}

// ---- Clade helpers ----

int count_clade_tips(const TreeState& tree, int node) {
  if (node < tree.n_tip) return 1;
  std::vector<int> stack;
  stack.push_back(node);
  int count = 0;
  while (!stack.empty()) {
    int nd = stack.back();
    stack.pop_back();
    if (nd < tree.n_tip) {
      ++count;
    } else {
      int ni = nd - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }
  return count;
}

// Collect all nodes (tips + internals) in the clade rooted at `node`.
// Returns them in no particular order.
static void collect_clade_nodes(const TreeState& tree, int node,
                                std::vector<int>& tips,
                                std::vector<int>& internals) {
  std::vector<int> stack;
  stack.push_back(node);
  while (!stack.empty()) {
    int nd = stack.back();
    stack.pop_back();
    if (nd < tree.n_tip) {
      tips.push_back(nd);
    } else {
      internals.push_back(nd);
      int ni = nd - tree.n_tip;
      stack.push_back(tree.left[ni]);
      stack.push_back(tree.right[ni]);
    }
  }
}

// ---- Reduced dataset construction ----

ReducedDataset build_reduced_dataset(const TreeState& tree,
                                     const DataSet& ds,
                                     int sector_root) {
  ReducedDataset rd;
  rd.sector_root = sector_root;
  rd.n_htus = 1;

  // Collect nodes in the sector clade
  std::vector<int> clade_tips, clade_internals;
  collect_clade_nodes(tree, sector_root, clade_tips, clade_internals);
  rd.n_real_tips = static_cast<int>(clade_tips.size());

  int n_sector_tips = rd.n_real_tips + rd.n_htus; // +1 for HTU
  int n_sector_internal = n_sector_tips - 1;
  int n_sector_node = 2 * n_sector_tips - 1;

  // Build mappings: full → sector and sector → full
  rd.full_to_sector.assign(tree.n_node, -1);
  rd.sector_to_full.assign(n_sector_node, -1);

  // Real tips: sector indices 0..n_real_tips-1
  for (int i = 0; i < rd.n_real_tips; ++i) {
    int full_node = clade_tips[i];
    rd.full_to_sector[full_node] = i;
    rd.sector_to_full[i] = full_node;
  }

  // HTU pseudo-tip: sector index n_real_tips
  int htu_sector_idx = rd.n_real_tips;
  // The HTU represents "rest of tree" — mapped to parent of sector_root
  int htu_full_node = tree.parent[sector_root];
  rd.sector_to_full[htu_sector_idx] = htu_full_node;
  // Don't set full_to_sector for htu_full_node — it's not in the sector

  // Internal nodes: sector indices n_sector_tips..n_sector_tips+n_sector_internal-1
  // The sector root maps to the sector's root (n_sector_tips).
  // We need n_sector_internal = n_sector_tips - 1 internal nodes, but
  // clade_internals has rd.n_real_tips - 1 nodes (since the clade has
  // n_real_tips leaves and n_real_tips-1 internal nodes). We need one
  // more internal node to accommodate the HTU connection.
  //
  // Sector topology: sector_root becomes sector root (index n_sector_tips).
  // Its parent-side connection becomes the HTU. We insert a new internal
  // node connecting the HTU to sector_root's children.
  //
  // Actually: the sector has n_real_tips real leaves. The clade rooted at
  // sector_root has (n_real_tips - 1) internal nodes including sector_root.
  // In the reduced tree we have n_sector_tips = n_real_tips + 1 tips and
  // need n_sector_tips - 1 = n_real_tips internal nodes. That's one more
  // internal node than the clade has. This extra node is the sector root
  // (index n_sector_tips), whose children are the HTU and sector_root's
  // full-tree equivalent.
  //
  // Wait, let me reconsider the topology mapping:
  //
  // In the full tree:
  //   parent(sector_root) -- sector_root -- {clade}
  //
  // In the sector tree (rooted):
  //   sector_root_new
  //     ├─ HTU (representing parent side)
  //     └─ old_sector_root_equivalent
  //           ├─ left subtree
  //           └─ right subtree
  //
  // No — simpler: the sector root IS the new root, with the HTU as one
  // child and the clade's internal structure below the other child.
  //
  // But sector_root in the full tree has two children (left, right).
  // In the sector tree, the root must also have two children. So:
  //   - If we make the sector root's children = [HTU, subtree], we lose
  //     the binary structure of the clade.
  //
  // Correct approach: the sector tree root has the HTU on one side and
  // the full clade below sector_root on the other. Since sector_root
  // in the full tree has left and right children, we need:
  //
  //   new_root
  //     ├─ HTU
  //     └─ sector_root_copy
  //           ├─ left child (mapped)
  //           └─ right child (mapped)
  //
  // This requires n_real_tips internal nodes = clade_internals.size() + 1.
  // clade_internals has (n_real_tips - 1) entries, so we need exactly one
  // extra: the new_root.

  // Map clade internals (excluding sector_root which gets special treatment)
  // Actually, sector_root IS in clade_internals. Let's map all clade
  // internals to sector internal indices, then the new_root gets the last
  // slot.

  int sector_int_idx = n_sector_tips; // first available internal index
  int new_root = sector_int_idx;      // sector root index

  // The new_root's children: HTU and sector_root_equivalent
  // sector_root_equivalent = the clade's root, mapped to a sector internal
  ++sector_int_idx; // new_root takes n_sector_tips

  for (int i = 0; i < static_cast<int>(clade_internals.size()); ++i) {
    int full_node = clade_internals[i];
    rd.full_to_sector[full_node] = sector_int_idx;
    rd.sector_to_full[sector_int_idx] = full_node;
    ++sector_int_idx;
  }

  // new_root itself maps to... nothing in the full tree (it's synthetic)
  rd.sector_to_full[new_root] = -1; // synthetic root

  // ---- Build sector TreeState ----

  rd.subtree.n_tip = n_sector_tips;
  rd.subtree.n_internal = n_sector_internal;
  rd.subtree.n_node = n_sector_node;
  rd.subtree.total_words = ds.total_words;
  rd.subtree.n_blocks = ds.n_blocks;

  rd.subtree.parent.assign(n_sector_node, -1);
  rd.subtree.left.assign(n_sector_internal, -1);
  rd.subtree.right.assign(n_sector_internal, -1);

  // Set up root
  rd.subtree.parent[new_root] = new_root; // root is its own parent

  // new_root's children: HTU (left) and sector_root's mapped node (right)
  int sr_mapped = rd.full_to_sector[sector_root];
  int nr_i = new_root - n_sector_tips;
  rd.subtree.left[nr_i] = htu_sector_idx;
  rd.subtree.right[nr_i] = sr_mapped;
  rd.subtree.parent[htu_sector_idx] = new_root;
  rd.subtree.parent[sr_mapped] = new_root;

  // Map internal topology from full tree
  for (int full_nd : clade_internals) {
    int sec_nd = rd.full_to_sector[full_nd];
    int sec_ni = sec_nd - n_sector_tips;
    int full_ni = full_nd - tree.n_tip;

    int full_lc = tree.left[full_ni];
    int full_rc = tree.right[full_ni];

    int sec_lc = rd.full_to_sector[full_lc];
    int sec_rc = rd.full_to_sector[full_rc];

    rd.subtree.left[sec_ni] = sec_lc;
    rd.subtree.right[sec_ni] = sec_rc;
    rd.subtree.parent[sec_lc] = sec_nd;
    rd.subtree.parent[sec_rc] = sec_nd;
  }

  // ---- Build sector DataSet ----

  // Copy block structure from original
  rd.data.n_tips = n_sector_tips;
  rd.data.n_blocks = ds.n_blocks;
  rd.data.total_words = ds.total_words;
  rd.data.blocks = ds.blocks;
  rd.data.block_word_offset = ds.block_word_offset;
  rd.data.n_patterns = ds.n_patterns;
  rd.data.min_steps = ds.min_steps;
  rd.data.pattern_freq = ds.pattern_freq;
  rd.data.concavity = ds.concavity;
  rd.data.eff_k = ds.eff_k;
  rd.data.phi = ds.phi;

  // Copy scoring mode and simplification metadata
  rd.data.scoring_mode = ds.scoring_mode;
  rd.data.ew_offset = ds.ew_offset;
  rd.data.precomputed_steps = ds.precomputed_steps;
  rd.data.info_amounts = ds.info_amounts;
  rd.data.info_max_steps = ds.info_max_steps;
  rd.data.inapp_state = ds.inapp_state;

  // Build tip states for the sector
  size_t tip_state_size = static_cast<size_t>(n_sector_tips) * ds.total_words;
  rd.data.tip_states.assign(tip_state_size, 0ULL);

  // Real tips: copy from original dataset
  for (int i = 0; i < rd.n_real_tips; ++i) {
    int full_tip = rd.sector_to_full[i];
    size_t src_base = static_cast<size_t>(full_tip) * ds.total_words;
    size_t dst_base = static_cast<size_t>(i) * ds.total_words;
    for (int w = 0; w < ds.total_words; ++w) {
      rd.data.tip_states[dst_base + w] = ds.tip_states[src_base + w];
    }
  }

  // HTU pseudo-tip: use from_above[sector_root] — the Fitch state-set
  // that the rest of the tree sends down to the sector boundary, excluding
  // the sector's own contribution. This gives a better HTU approximation
  // than final_[parent], which circularly includes the sector's states.
  {
    std::vector<uint64_t> from_above_sr;
    compute_from_above_for_sector(tree, ds, sector_root, from_above_sr);

    size_t dst_base =
        static_cast<size_t>(htu_sector_idx) * ds.total_words;
    for (int w = 0; w < ds.total_words; ++w) {
      rd.data.tip_states[dst_base + w] = from_above_sr[w];
    }
  }

  // Allocate state arrays and load tip states
  size_t state_size = static_cast<size_t>(n_sector_node) * ds.total_words;
  rd.subtree.prelim.assign(state_size, 0ULL);
  rd.subtree.final_.assign(state_size, 0ULL);
  rd.subtree.down2.assign(state_size, 0ULL);
  rd.subtree.subtree_actives.assign(state_size, 0ULL);
  rd.subtree.local_cost.assign(
      static_cast<size_t>(n_sector_node) * ds.n_blocks, 0ULL);

  rd.subtree.load_tip_states(rd.data);
  rd.subtree.build_postorder();

  return rd;
}

// ---- Sector search ----

// Note: Sector trees use score_tree() which dispatches appropriately.
// For sectors with NA characters, the HTU subtree_actives approximation
// means the sector score is inexact, but the full-tree rescore after
// reinsertion catches any discrepancies.

// Search the reduced dataset and return the best score found.
// Modifies rd.subtree in place.
// Search the sector with TBR. The internal_ratchet_cycles parameter is
// reserved for future use (ratchet perturbation within sectors).
static double search_sector(ReducedDataset& rd, int /*internal_ratchet_cycles*/,
                            int max_hits) {
  int htu_idx = rd.n_real_tips;
  int root = rd.subtree.n_tip;
  int sr_mapped = rd.full_to_sector[rd.sector_root];

  // Save topology in case TBR disrupts root structure
  auto save_left = rd.subtree.left;
  auto save_right = rd.subtree.right;
  auto save_parent = rd.subtree.parent;

  double original_score = score_tree(rd.subtree, rd.data);

  TBRParams tp;
  tp.max_hits = max_hits;
  TBRResult tr = tbr_search(rd.subtree, rd.data, tp);

  // Verify root structure: HTU and sector_root_mapped must remain
  // direct children of the synthetic root. TBR can regraft onto root
  // edges, displacing nodes outside the clade — if so, discard result.
  int root_i = root - rd.subtree.n_tip;
  int root_lc = rd.subtree.left[root_i];
  int root_rc = rd.subtree.right[root_i];
  bool root_ok = (root_lc == htu_idx && root_rc == sr_mapped) ||
                 (root_lc == sr_mapped && root_rc == htu_idx);

  if (!root_ok) {
    rd.subtree.left = save_left;
    rd.subtree.right = save_right;
    rd.subtree.parent = save_parent;
    rd.subtree.build_postorder();
    return original_score;
  }

  return tr.best_score;
}

// ---- Reinsertion ----

// Reinsert the improved sector topology into the full tree.
// Only touches nodes within the sector clade.
static void reinsert_sector(TreeState& tree, const ReducedDataset& rd) {
  int n_sector_tips = rd.subtree.n_tip;
  int sector_root_mapped = rd.full_to_sector[rd.sector_root];

  // The sector tree's root (n_sector_tips) is synthetic — its right child
  // is the mapped sector_root. We only care about the subtree below
  // sector_root_mapped.
  //
  // Walk the sector tree below sector_root_mapped and write topology
  // back to the full tree.

  std::vector<int> stack;
  stack.push_back(sector_root_mapped);

  while (!stack.empty()) {
    int sec_nd = stack.back();
    stack.pop_back();

    if (sec_nd < n_sector_tips) continue; // sector tip — no children to map

    int full_nd = rd.sector_to_full[sec_nd];
    if (full_nd < 0) continue; // synthetic node (root)

    int sec_ni = sec_nd - n_sector_tips;
    int sec_lc = rd.subtree.left[sec_ni];
    int sec_rc = rd.subtree.right[sec_ni];

    // Map sector children to full tree nodes
    int full_lc = rd.sector_to_full[sec_lc];
    int full_rc = rd.sector_to_full[sec_rc];

    // Update full tree topology
    int full_ni = full_nd - tree.n_tip;
    tree.left[full_ni] = full_lc;
    tree.right[full_ni] = full_rc;
    tree.parent[full_lc] = full_nd;
    tree.parent[full_rc] = full_nd;

    stack.push_back(sec_lc);
    stack.push_back(sec_rc);
  }
}

// ---- XSS partitioning ----

// Partition the tree into approximately equal-sized non-overlapping sectors.
// Returns a vector of sector root node indices.
//
// O(n) algorithm: maintain unclaimed_below[] counts; when a sector is
// claimed, subtract from ancestors via a rootward walk. Each node is
// visited O(1) times in the main loop; rootward walks are O(height)
// each and there are at most n_partitions of them.
static std::vector<int> xss_partition(const TreeState& tree, int n_partitions) {
  std::vector<int> subtree_size(tree.n_node, 0);
  for (int i = 0; i < tree.n_tip; ++i) {
    subtree_size[i] = 1;
  }
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    subtree_size[node] = subtree_size[tree.left[ni]]
                       + subtree_size[tree.right[ni]];
  }

  int target = tree.n_tip / n_partitions;
  if (target < 4) target = 4;

  // unclaimed_below[node] = number of unclaimed tips in node's subtree.
  // Starts equal to subtree_size; reduced when descendant sectors are claimed.
  std::vector<int> unclaimed_below = subtree_size;

  int root = tree.n_tip;
  std::vector<int> sectors;

  for (int node : tree.postorder) {
    if (node == root) continue;

    if (unclaimed_below[node] >= target) {
      sectors.push_back(node);

      // Subtract claimed tips from all ancestors up to and including root
      int tips_claimed = unclaimed_below[node];
      for (int cur = tree.parent[node]; ; cur = tree.parent[cur]) {
        unclaimed_below[cur] -= tips_claimed;
        if (cur == root) break;
      }
      unclaimed_below[node] = 0;
    }
  }

  return sectors;
}

// ---- RSS ----

SectorResult rss_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd) {
  bool constrained = cd && cd->active && cd->has_posthoc;
  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  // Ensure full tree has current state sets
  double current_score = score_tree(tree, ds);

  SectorResult result;
  result.best_score = current_score;
  result.n_sectors_searched = 0;
  result.n_sectors_improved = 0;
  result.total_steps_saved = 0;

  int avg_size = (params.min_sector_size + params.max_sector_size) / 2;
  int n_picks = params.rss_picks_per_round;
  if (n_picks <= 0) {
    n_picks = std::max(1, 2 * tree.n_tip / std::max(1, avg_size));
  }

  // Precompute subtree sizes for sector selection
  std::vector<int> subtree_size(tree.n_node, 0);
  for (int i = 0; i < tree.n_tip; ++i) subtree_size[i] = 1;
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    subtree_size[node] = subtree_size[tree.left[ni]]
                       + subtree_size[tree.right[ni]];
  }

  // Collect eligible internal nodes (not root)
  std::vector<int> eligible;
  for (int node = tree.n_tip + 1; node < tree.n_node; ++node) {
    int sz = subtree_size[node];
    if (sz >= params.min_sector_size && sz <= params.max_sector_size) {
      eligible.push_back(node);
    }
  }

  if (eligible.empty()) {
    // No sectors of appropriate size; run global TBR and return
    TBRParams tp;
    tp.max_hits = params.internal_max_hits;
    TBRResult tr = tbr_search(tree, ds, tp);
    result.best_score = tr.best_score;
    result.total_steps_saved =
        static_cast<int>(current_score - tr.best_score);
    return result;
  }

  // Conflict-guided weighting: if pool frequency data is available,
  // bias random sector selection toward high-conflict regions.
  std::vector<double> pick_weights;
  bool use_weighted = false;
  if (params.split_freq && params.split_freq->n_trees >= 2) {
    std::vector<double> node_conf, sector_conf;
    compute_node_conflict(tree, *params.split_freq, node_conf, sector_conf);

    pick_weights.resize(eligible.size());
    double max_w = 0.0;
    for (size_t i = 0; i < eligible.size(); ++i) {
      // Blend: base weight 1.0 + conflict bonus (up to 3.0 extra).
      // Ensures even low-conflict sectors get some chance.
      pick_weights[i] = 1.0 + 3.0 * sector_conf[eligible[i]];
      if (pick_weights[i] > max_w) max_w = pick_weights[i];
    }
    // Only use weighted selection if there is meaningful variation
    use_weighted = (max_w > 1.5);
  }

  for (int pick = 0; pick < n_picks; ++pick) {
    int idx;
    if (use_weighted) {
      std::discrete_distribution<int> dist(pick_weights.begin(),
                                           pick_weights.end());
      idx = dist(rng);
    } else {
      idx = std::uniform_int_distribution<int>(
          0, static_cast<int>(eligible.size()) - 1)(rng);
    }
    int sector_root = eligible[idx];

    // State arrays are guaranteed valid: either from the initial
    // score_tree above, or from the previous iteration's acceptance
    // (which calls score_tree) or rejection (which also rescores).

    // Build reduced dataset
    ReducedDataset rd = build_reduced_dataset(tree, ds, sector_root);

    // Score the current sector topology
    double sector_current = score_tree(rd.subtree, rd.data);

    // Search the sector
    double sector_best = search_sector(rd, params.internal_ratchet_cycles,
                                       params.internal_max_hits);
    ++result.n_sectors_searched;

    bool improved = sector_best < sector_current;
    bool accept = improved ||
                  (params.accept_equal && sector_best == sector_current);

    if (accept && sector_best <= sector_current) {
      // Save only the sector clade's topology for potential undo
      CladeSnapshot snap;
      save_clade(tree, sector_root, snap);

      reinsert_sector(tree, rd);
      tree.build_postorder();
      double new_score = score_tree(tree, ds);

      // Post-hoc constraint check
      if (constrained && violates_constraint_posthoc(tree, *cd)) {
        restore_clade(tree, snap);
        tree.build_postorder();
        score_tree(tree, ds);
        continue;
      }

      if (new_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - new_score);
        result.best_score = new_score;
        ++result.n_sectors_improved;

        // Recompute subtree sizes and eligible list after topology change
        for (int i = 0; i < tree.n_tip; ++i) subtree_size[i] = 1;
        for (int node : tree.postorder) {
          int ni = node - tree.n_tip;
          subtree_size[node] = subtree_size[tree.left[ni]]
                             + subtree_size[tree.right[ni]];
        }
        eligible.clear();
        for (int node = tree.n_tip + 1; node < tree.n_node; ++node) {
          int sz = subtree_size[node];
          if (sz >= params.min_sector_size &&
              sz <= params.max_sector_size) {
            eligible.push_back(node);
          }
        }
        if (eligible.empty()) break;

        // Recompute conflict weights for new topology
        if (use_weighted) {
          std::vector<double> nc, sc;
          compute_node_conflict(tree, *params.split_freq, nc, sc);
          pick_weights.resize(eligible.size());
          for (size_t i = 0; i < eligible.size(); ++i) {
            pick_weights[i] = 1.0 + 3.0 * sc[eligible[i]];
          }
        }
      } else if (new_score == result.best_score && params.accept_equal) {
        // Equal score accepted — topology changed but score didn't
      } else {
        // HTU approximation caused full-tree score to worsen; revert
        restore_clade(tree, snap);
        tree.build_postorder();
        score_tree(tree, ds);
      }
    }

    if (ts::check_interrupt()) break;
  }

  // Global TBR after all sector picks
  {
    TBRParams tp;
    tp.max_hits = params.internal_max_hits;
    TBRResult tr = tbr_search(tree, ds, tp, cd);
    if (tr.best_score < result.best_score) {
      result.total_steps_saved +=
          static_cast<int>(result.best_score - tr.best_score);
      result.best_score = tr.best_score;
    }
  }

  return result;
}

// ---- XSS ----

SectorResult xss_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd) {
  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  double current_score = score_tree(tree, ds);

  SectorResult result;
  result.best_score = current_score;
  result.n_sectors_searched = 0;
  result.n_sectors_improved = 0;
  result.total_steps_saved = 0;

  bool constrained = cd && cd->active && cd->has_posthoc;

  for (int round = 0; round < params.xss_rounds; ++round) {
    double score_before_round = result.best_score;

    // Pick a random number of partitions around the target
    int n_parts = params.n_partitions;
    // Add some randomness: ±1
    if (params.n_partitions > 2) {
      int delta = std::uniform_int_distribution<int>(-1, 1)(rng);
      n_parts = std::max(2, params.n_partitions + delta);
    }

    // Partition the tree
    std::vector<int> sectors = xss_partition(tree, n_parts);

    // Search each sector
    for (int sector_root : sectors) {
      // Verify sector is still valid (topology may have changed)
      int sz = count_clade_tips(tree, sector_root);
      if (sz < 4) continue; // too small to be useful

      // State arrays are guaranteed valid: either from the initial
      // score_tree above, or from the previous sector's acceptance/
      // rejection (both paths call score_tree before continuing).

      ReducedDataset rd = build_reduced_dataset(tree, ds, sector_root);

      double sector_current = score_tree(rd.subtree, rd.data);
      double sector_best = search_sector(
          rd, params.internal_ratchet_cycles, params.internal_max_hits);
      ++result.n_sectors_searched;

      bool improved = sector_best < sector_current;
      bool accept = improved ||
                    (params.accept_equal && sector_best == sector_current);

      if (accept && sector_best <= sector_current) {
        // Save only the sector clade's topology for potential undo
        CladeSnapshot snap;
        save_clade(tree, sector_root, snap);

        reinsert_sector(tree, rd);
        tree.build_postorder();
        double new_score = score_tree(tree, ds);

        // Post-hoc constraint check
        if (constrained && violates_constraint_posthoc(tree, *cd)) {
          restore_clade(tree, snap);
          tree.build_postorder();
          score_tree(tree, ds);
          continue;
        }

        if (new_score < result.best_score) {
          result.total_steps_saved +=
              static_cast<int>(result.best_score - new_score);
          result.best_score = new_score;
          ++result.n_sectors_improved;
        } else if (new_score == result.best_score && params.accept_equal) {
          // Equal score accepted
        } else {
          // HTU approximation caused full-tree score to worsen; revert
          restore_clade(tree, snap);
          tree.build_postorder();
          score_tree(tree, ds);
        }
      }

      if (ts::check_interrupt()) break;
    }

    // Global TBR after each round of sectors
    {
      TBRParams tp;
      tp.max_hits = params.internal_max_hits;
      TBRResult tr = tbr_search(tree, ds, tp, cd);
      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
      }
    }

    // Adaptive: skip remaining rounds if this one found no improvement
    if (result.best_score >= score_before_round) break;
    if (ts::check_interrupt()) break;
  }

  return result;
}

// ---- CSS (Constrained Sectorial Search) ----
//
// Sector-restricted TBR on the full tree. No HTU approximation —
// scoring is exact against the full dataset.

// Build a sector mask: true for all nodes in the clade rooted at sector_root.
static void build_sector_mask(const TreeState& tree, int sector_root,
                               std::vector<bool>& mask) {
  mask.assign(tree.n_node, false);
  std::vector<int> stk;
  stk.push_back(sector_root);
  while (!stk.empty()) {
    int nd = stk.back();
    stk.pop_back();
    mask[nd] = true;
    if (nd >= tree.n_tip) {
      int ni = nd - tree.n_tip;
      stk.push_back(tree.left[ni]);
      stk.push_back(tree.right[ni]);
    }
  }
}

SectorResult css_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd) {
  double current_score = score_tree(tree, ds);

  SectorResult result;
  result.best_score = current_score;
  result.n_sectors_searched = 0;
  result.n_sectors_improved = 0;
  result.total_steps_saved = 0;

  int n_rounds = params.xss_rounds;
  if (n_rounds <= 0) n_rounds = 1;

  std::vector<bool> sector_mask;

  for (int round = 0; round < n_rounds; ++round) {
    double score_before_round = result.best_score;

    std::vector<int> sectors = xss_partition(tree, params.n_partitions);

    for (int sector_root : sectors) {
      int sz = count_clade_tips(tree, sector_root);
      if (sz < 4) continue;

      build_sector_mask(tree, sector_root, sector_mask);

      TBRParams tp;
      tp.max_hits = params.internal_max_hits;

      TBRResult tr = tbr_search(tree, ds, tp, cd, &sector_mask);
      ++result.n_sectors_searched;

      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
        ++result.n_sectors_improved;
      }

      if (ts::check_interrupt()) break;
    }

    // Global TBR after each round
    {
      TBRParams tp;
      tp.max_hits = params.internal_max_hits;
      TBRResult tr = tbr_search(tree, ds, tp, cd);
      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
      }
    }

    // Adaptive: skip remaining rounds if this one found no improvement
    if (result.best_score >= score_before_round) break;
    if (ts::check_interrupt()) break;
  }

  return result;
}

} // namespace ts
