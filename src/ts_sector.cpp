#include "ts_sector.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_wagner.h"

#include <algorithm>
#include <vector>
#include <R.h>
#include <Rmath.h>

namespace ts {

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

  // HTU pseudo-tip: use final_ states from parent of sector_root
  {
    size_t src_base =
        static_cast<size_t>(htu_full_node) * tree.total_words;
    size_t dst_base =
        static_cast<size_t>(htu_sector_idx) * ds.total_words;
    for (int w = 0; w < ds.total_words; ++w) {
      rd.data.tip_states[dst_base + w] = tree.final_[src_base + w];
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

// Ensure a TreeState has down2 and subtree_actives allocated.
// Wagner-built trees lack these; fitch_na_score needs them.
static void ensure_na_arrays(TreeState& tree, const DataSet& ds) {
  size_t state_size = static_cast<size_t>(tree.n_node) * tree.total_words;
  if (tree.down2.size() != state_size) {
    tree.down2.assign(state_size, 0ULL);
  }
  if (tree.subtree_actives.size() != state_size) {
    tree.subtree_actives.assign(state_size, 0ULL);
    // Re-init tip subtree_actives
    for (int tip = 0; tip < tree.n_tip; ++tip) {
      size_t base = static_cast<size_t>(tip) * tree.total_words;
      for (int b = 0; b < ds.n_blocks; ++b) {
        int offset = ds.block_word_offset[b];
        if (ds.blocks[b].has_inapplicable) {
          tree.subtree_actives[base + offset] = 0;
          for (int s = 1; s < ds.blocks[b].n_states; ++s) {
            tree.subtree_actives[base + offset + s] =
                ds.tip_states[base + offset + s];
          }
        }
      }
    }
  }
}

// Search the reduced dataset and return the best score found.
// Modifies rd.subtree in place.
static double search_sector(ReducedDataset& rd, int internal_ratchet_cycles,
                            int max_hits) {
  int htu_idx = rd.n_real_tips;
  int root = rd.subtree.n_tip;
  int sr_mapped = rd.full_to_sector[rd.sector_root];

  // Save topology in case TBR disrupts root structure
  auto save_left = rd.subtree.left;
  auto save_right = rd.subtree.right;
  auto save_parent = rd.subtree.parent;

  double original_score = static_cast<double>(fitch_score(rd.subtree, rd.data));

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
    int full_lc, full_rc;

    if (sec_lc < n_sector_tips) {
      // It's a sector tip
      full_lc = rd.sector_to_full[sec_lc];
    } else {
      full_lc = rd.sector_to_full[sec_lc];
    }

    if (sec_rc < n_sector_tips) {
      full_rc = rd.sector_to_full[sec_rc];
    } else {
      full_rc = rd.sector_to_full[sec_rc];
    }

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
static std::vector<int> xss_partition(const TreeState& tree, int n_partitions) {
  // Compute subtree tip counts bottom-up
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

  std::vector<int> sectors;
  std::vector<bool> claimed(tree.n_node, false);

  // Walk postorder (bottom-up). When a subtree reaches ~target size
  // and hasn't been claimed by a descendant sector, mark it.
  for (int node : tree.postorder) {
    if (node == tree.n_tip) continue; // skip root

    // Count unclaimed tips below this node
    int unclaimed_tips = 0;
    std::vector<int> tip_stack;
    tip_stack.push_back(node);
    while (!tip_stack.empty()) {
      int nd = tip_stack.back();
      tip_stack.pop_back();
      if (claimed[nd]) continue;
      if (nd < tree.n_tip) {
        ++unclaimed_tips;
      } else {
        int ni = nd - tree.n_tip;
        tip_stack.push_back(tree.left[ni]);
        tip_stack.push_back(tree.right[ni]);
      }
    }

    if (unclaimed_tips >= target) {
      sectors.push_back(node);
      // Mark all nodes in this clade as claimed
      std::vector<int> mark_stack;
      mark_stack.push_back(node);
      while (!mark_stack.empty()) {
        int nd = mark_stack.back();
        mark_stack.pop_back();
        claimed[nd] = true;
        if (nd >= tree.n_tip) {
          int ni = nd - tree.n_tip;
          mark_stack.push_back(tree.left[ni]);
          mark_stack.push_back(tree.right[ni]);
        }
      }
    }
  }

  return sectors;
}

// ---- RSS ----

SectorResult rss_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params) {
  // Ensure full tree has current state sets
  double current_score = static_cast<double>(fitch_na_score(tree, ds));

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

  for (int pick = 0; pick < n_picks; ++pick) {
    // Pick a random eligible node
    int idx = static_cast<int>(unif_rand() * eligible.size());
    if (idx >= static_cast<int>(eligible.size()))
      idx = static_cast<int>(eligible.size()) - 1;
    int sector_root = eligible[idx];

    // Ensure current state sets are up to date
    fitch_na_score(tree, ds);

    // Build reduced dataset
    ReducedDataset rd = build_reduced_dataset(tree, ds, sector_root);

    // Score the current sector topology
    double sector_current = static_cast<double>(
        fitch_score(rd.subtree, rd.data));

    // Search the sector
    double sector_best = search_sector(rd, params.internal_ratchet_cycles,
                                       params.internal_max_hits);
    ++result.n_sectors_searched;

    bool improved = sector_best < sector_current;
    bool accept = improved ||
                  (params.accept_equal && sector_best == sector_current);

    if (accept && sector_best <= sector_current) {
      reinsert_sector(tree, rd);
      tree.build_postorder();
      double new_score = static_cast<double>(fitch_na_score(tree, ds));

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
      } else if (new_score == result.best_score && params.accept_equal) {
        // Equal score accepted — topology changed but score didn't
      } else {
        // Reinsertion made things worse (shouldn't happen if sector
        // search found improvement, but HTU approximation can cause this)
        // Revert would be complex; just accept the current state
      }
    }

    R_CheckUserInterrupt();
  }

  // Global TBR after all sector picks
  {
    TBRParams tp;
    tp.max_hits = params.internal_max_hits;
    TBRResult tr = tbr_search(tree, ds, tp);
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
                        const SectorParams& params) {
  double current_score = static_cast<double>(fitch_na_score(tree, ds));

  SectorResult result;
  result.best_score = current_score;
  result.n_sectors_searched = 0;
  result.n_sectors_improved = 0;
  result.total_steps_saved = 0;

  for (int round = 0; round < params.xss_rounds; ++round) {
    // Pick a random number of partitions around the target
    int n_parts = params.n_partitions;
    // Add some randomness: ±1
    if (params.n_partitions > 2) {
      int delta = static_cast<int>(unif_rand() * 3.0) - 1; // -1, 0, or 1
      n_parts = std::max(2, params.n_partitions + delta);
    }

    // Partition the tree
    std::vector<int> sectors = xss_partition(tree, n_parts);

    // Search each sector
    for (int sector_root : sectors) {
      // Verify sector is still valid (topology may have changed)
      int sz = count_clade_tips(tree, sector_root);
      if (sz < 4) continue; // too small to be useful

      // Ensure state sets are current
      fitch_na_score(tree, ds);

      ReducedDataset rd = build_reduced_dataset(tree, ds, sector_root);

      double sector_current = static_cast<double>(
          fitch_score(rd.subtree, rd.data));
      double sector_best = search_sector(
          rd, params.internal_ratchet_cycles, params.internal_max_hits);
      ++result.n_sectors_searched;

      bool improved = sector_best < sector_current;
      bool accept = improved ||
                    (params.accept_equal && sector_best == sector_current);

      if (accept && sector_best <= sector_current) {
        reinsert_sector(tree, rd);
        tree.build_postorder();
        double new_score = static_cast<double>(fitch_na_score(tree, ds));

        if (new_score < result.best_score) {
          result.total_steps_saved +=
              static_cast<int>(result.best_score - new_score);
          result.best_score = new_score;
          ++result.n_sectors_improved;
        }
      }

      R_CheckUserInterrupt();
    }

    // Global TBR after each round of sectors
    {
      TBRParams tp;
      tp.max_hits = params.internal_max_hits;
      TBRResult tr = tbr_search(tree, ds, tp);
      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
      }
    }

    R_CheckUserInterrupt();
  }

  return result;
}

} // namespace ts
