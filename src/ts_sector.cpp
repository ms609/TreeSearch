#include "ts_sector.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_wagner.h"
#include "ts_rng.h"
#include "ts_splits.h"
#include "ts_pool.h"
#include "ts_drift.h"

#include <algorithm>
#include <random>
#include <vector>
#include <cstring>
#include <climits>
#include <cstdlib>
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

  // 3. Walk down the path, computing from_above at each child step.
  // `new_from_above` is allocated ONCE and swapped each step (O(1)) instead of
  // heap-allocated per step. Byte-identical: the inner loop overwrites every
  // state word each step, and any padding words (total_words > sum n_states)
  // start at 0 in both buffers and are never written, so they stay 0 — same as
  // the original fresh-zeroed allocation. (Sectorial micro-bank, T-S6c.)
  std::vector<uint64_t> new_from_above(tw);
  for (size_t i = 0; i + 1 < path.size(); ++i) {
    int node = path[i];
    int next = path[i + 1]; // child on the path

    // Find sibling of `next` under `node`
    int ni = node - tree.n_tip;
    int sib = (tree.left[ni] == next) ? tree.right[ni] : tree.left[ni];

    // from_above[next] = fitch_join(from_above[node], prelim[sib])
    // fitch_join: per-block, compute intersection; where empty, use union.
    // When tw == 0 (e.g. every character constant/autapomorphic under equal
    // weights, leaving zero Fitch blocks), tree.prelim is empty and the
    // block loop below never executes (ds.n_blocks == 0 too), so sib_prelim
    // is never dereferenced -- but taking its address still needs guarding
    // (address of element 0 of an empty vector is undefined behaviour;
    // aborts under _GLIBCXX_ASSERTIONS).
    const uint64_t* sib_prelim = (tw > 0)
        ? &tree.prelim[static_cast<size_t>(sib) * tw]
        : nullptr;

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
    std::swap(from_above_cur, new_from_above);
  }

  if (tw > 0) {
    std::memcpy(from_above_out.data(), from_above_cur.data(),
                tw * sizeof(uint64_t));
  }
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

#ifdef TS_AUDIT_PROBE
#include <cstdio>
// Audit #56: realized per-sector column-axis reduction counters.
long long g_sect_inf_chars = 0, g_sect_tot_chars = 0;
long long g_sect_fp_blocks = 0, g_sect_tot_blocks = 0;
unsigned long long g_sect_calls = 0;
#endif

// ---- Per-sector column-axis reduction (#56, EW only, opt-in TS_SECT_COLREDUCE) ----
// A character CONSTANT-within-{sector tips + HTU} (some single state shared by
// every tip) contributes 0 Fitch steps in every sector topology, so dropping it
// leaves every candidate's score UNCHANGED (ew_offset is untouched) -> the inner
// sector-search trajectory is byte-identical. Re-packing the informative
// survivors into fewer 64-char blocks shrinks the per-node block scan that
// dominates the inner-sector TBR -- especially the no-bail compute_insertion_
// edge_sets precompute (which scans every block at every node). Gated to the
// plain EW case (weight 1, no ratchet upweighting, no inapplicable) so the
// "0-step" invariance argument holds exactly.
static const bool kSectColReduce = []{
  const char* e = std::getenv("TS_SECT_COLREDUCE");
  return e != nullptr && e[0] == '1';
}();

static void reduce_sector_columns_ew(ReducedDataset& rd, int n_sector_tips) {
  DataSet& d = rd.data;
  if (d.scoring_mode != ScoringMode::EW || !d.all_weight_one) return;
  for (int b = 0; b < d.n_blocks; ++b)
    if (d.blocks[b].upweight_mask != 0 || d.blocks[b].has_inapplicable) return;

  const int otw = d.total_words;

  // 1. Collect informative survivors (n_states, old word offset, old bit).
  struct Surv { int n_states; int old_off; int old_bit; };
  std::vector<Surv> surv;
  int tot_active = 0;
  for (int b = 0; b < d.n_blocks; ++b) {
    const uint64_t active = d.blocks[b].active_mask;
    if (!active) continue;
    const int off = d.block_word_offset[b];
    const int nst = d.blocks[b].n_states;
    tot_active += popcount64(active);
    uint64_t constant = 0ULL;            // bit c set iff some state shared by ALL tips
    for (int s = 0; s < nst; ++s) {
      uint64_t all = ~0ULL;
      for (int t = 0; t < n_sector_tips; ++t)
        all &= d.tip_states[static_cast<size_t>(t) * otw + off + s];
      constant |= all;
    }
    uint64_t inf = active & ~constant;
    while (inf) {
      const int c = static_cast<int>(__builtin_ctzll(inf));
      inf &= inf - 1;
      surv.push_back({nst, off, c});
    }
  }
  if (surv.empty() ||
      static_cast<int>(surv.size()) >= tot_active) return;  // nothing to drop

  // 2. Group survivors by n_states (stable within group) -> one n_states/block.
  std::stable_sort(surv.begin(), surv.end(),
                   [](const Surv& a, const Surv& b) { return a.n_states < b.n_states; });

  // 3. New block layout: <=64 survivors per block, all of one n_states.
  std::vector<CharBlock> nb;
  std::vector<int> nbwo;
  int ntw = 0;
  for (size_t i = 0; i < surv.size();) {
    const int nst = surv[i].n_states;
    size_t j = i;
    while (j < surv.size() && surv[j].n_states == nst &&
           (j - i) < static_cast<size_t>(MAX_CHARS_PER_BLOCK)) ++j;
    const int nchar = static_cast<int>(j - i);
    CharBlock cb{};
    cb.n_chars = nchar;
    cb.n_states = nst;
    cb.weight = 1;
    cb.has_inapplicable = false;
    cb.active_mask = (nchar == 64) ? ~0ULL : ((1ULL << nchar) - 1ULL);
    cb.upweight_mask = 0;
    for (int k = 0; k < nchar; ++k) cb.pattern_index[k] = k;  // identity (EW)
    nbwo.push_back(ntw);
    ntw += nst;
    nb.push_back(cb);
    i = j;
  }
  const int nblk = static_cast<int>(nb.size());

  // 4. Repack tip_states: move each survivor's state-words to its new bit/offset.
  std::vector<uint64_t> nts(static_cast<size_t>(d.n_tips) * ntw, 0ULL);
  size_t si = 0;
  for (int blk = 0; blk < nblk; ++blk) {
    const int nst   = nb[blk].n_states;
    const int noff  = nbwo[blk];
    const int nchar = nb[blk].n_chars;
    for (int local = 0; local < nchar; ++local, ++si) {
      const int o_off = surv[si].old_off;
      const int o_bit = surv[si].old_bit;
      for (int t = 0; t < d.n_tips; ++t) {
        const size_t src = static_cast<size_t>(t) * otw + o_off;
        const size_t dst = static_cast<size_t>(t) * ntw + noff;
        for (int s = 0; s < nst; ++s) {
          const uint64_t bit = (d.tip_states[src + s] >> o_bit) & 1ULL;
          nts[dst + s] |= (bit << local);
        }
      }
    }
  }

  // 5. flat_blocks mirror the new blocks.
  std::vector<FlatBlock> nfb(nblk);
  for (int blk = 0; blk < nblk; ++blk) {
    nfb[blk].offset = nbwo[blk];
    nfb[blk].n_states = nb[blk].n_states;
    nfb[blk].active_mask = nb[blk].active_mask;
    nfb[blk].has_inapplicable = 0;
  }

  // 6. Commit. Per-pattern arrays (n_patterns/min_steps/precomputed_steps/
  //    pattern_freq) and ew_offset are left untouched: the EW indirect scorers
  //    read ONLY the block structure, and the dropped chars' 0 contribution keeps
  //    scores exact. (Extending past the EW gate would require remapping those
  //    per-pattern arrays over the survivors.)
#ifdef TS_AUDIT_PROBE
  static bool announced = false;
  if (!announced) {
    announced = true;
    Rprintf("COLREDUCE active: blocks %d->%d, words %d->%d\n",
            d.n_blocks, nblk, otw, ntw);
  }
#endif
  d.blocks = std::move(nb);
  d.block_word_offset = std::move(nbwo);
  d.flat_blocks = std::move(nfb);
  d.tip_states = std::move(nts);
  d.total_words = ntw;
  d.n_blocks = nblk;
  // CRITICAL: the TreeState carries its OWN stride fields (set from ds.* before
  // this call); every consumer of the reduced state arrays (load_tip_states,
  // fitch_downpass, compute_insertion_edge_sets, the inner tbr_search undo) reads
  // the stride from rd.subtree, NOT rd.data. Re-sync them or they index the
  // new-sized buffers with the old stride (heap OOB + wrong scores).
  rd.subtree.total_words = ntw;
  rd.subtree.n_blocks = nblk;
}

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
  rd.data.flat_blocks = ds.flat_blocks;
  rd.data.all_weight_one = ds.all_weight_one;
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

#ifdef TS_AUDIT_PROBE
  // Audit #56: a character is CONSTANT-within-sector iff some single state is
  // shared by ALL sector tips (incl. the HTU pseudo-tip) -> 0 Fitch steps in
  // every sector topology -> droppable (ranking-preserving). fp_blocks = blocks
  // the informative survivors re-pack into; drives the no-bail precompute saving
  // in compute_insertion_edge_sets (which scans all n_blocks per node).
  {
    const int tw_ = rd.data.total_words;
    long long inf_chars = 0, tot_chars = 0;
    for (int b = 0; b < rd.data.n_blocks; ++b) {
      uint64_t active = rd.data.blocks[b].active_mask;
      if (!active) continue;
      const int off = rd.data.block_word_offset[b];
      const int nst = rd.data.blocks[b].n_states;
      uint64_t constant = 0ULL;
      for (int s = 0; s < nst; ++s) {
        uint64_t all = ~0ULL;
        for (int t = 0; t < n_sector_tips; ++t)
          all &= rd.data.tip_states[static_cast<size_t>(t) * tw_ + off + s];
        constant |= all;
      }
      const uint64_t informative = active & ~constant;
      inf_chars += ts::popcount64(informative);
      tot_chars += ts::popcount64(active);
    }
    const long long fp_blocks = (inf_chars + 63) / 64;
    g_sect_inf_chars += inf_chars; g_sect_tot_chars += tot_chars;
    g_sect_fp_blocks += fp_blocks; g_sect_tot_blocks += rd.data.n_blocks;
    if ((++g_sect_calls % 2000ULL) == 0)
      std::fprintf(stderr,
        "SECT_REDUCE sectors=%llu inf/tot_chars=%.3f fp/tot_blocks=%.3f\n",
        (unsigned long long)g_sect_calls,
        g_sect_tot_chars ? (double)g_sect_inf_chars / (double)g_sect_tot_chars : 0.0,
        g_sect_tot_blocks ? (double)g_sect_fp_blocks / (double)g_sect_tot_blocks : 0.0);
  }
#endif

  // #56: optionally drop constant-within-sector columns + re-pack (EW only).
  if (kSectColReduce) reduce_sector_columns_ew(rd, n_sector_tips);

  // Allocate state arrays and load tip states (sized by rd.data, which the
  // reduction above may have shrunk; equals ds.* when the reduction is off).
  size_t state_size = static_cast<size_t>(n_sector_node) * rd.data.total_words;
  rd.subtree.prelim.assign(state_size, 0ULL);
  rd.subtree.final_.assign(state_size, 0ULL);
  rd.subtree.down2.assign(state_size, 0ULL);
  rd.subtree.subtree_actives.assign(state_size, 0ULL);
  rd.subtree.local_cost.assign(
      static_cast<size_t>(n_sector_node) * rd.data.n_blocks, 0ULL);

  rd.subtree.load_tip_states(rd.data);
  rd.subtree.build_postorder();

  return rd;
}

// Build a reduced dataset where deep sub-clades are COLLAPSED into composite
// terminals -- Goloboff 1999's reduced dataset ("internal nodes represented by
// their first-pass state sets"). The clade at sector_root is pruned at a frontier
// of ~target_tips nodes; each frontier node is a reduced "tip" (a real tip keeps
// its states; an internal sub-clade root takes its prelim/first-pass state set).
// The reduced internals are the skeleton ABOVE the frontier, so TBR/RAS rearranges
// the major sub-clades relative to one another -- the coarse-grained move TNT's
// sectsch makes that a fully-resolved contiguous sector cannot reach.
//
// reinsert_sector() works unchanged: a collapsed root is a sector TIP, so its
// full-tree sub-clade internals are never remapped (left intact) and the root is
// merely reparented to its new skeleton position. Scoring uses the composite
// state sets (approximate, like the HTU); the full-tree rescore after reinsertion
// is exact and gates acceptance.
static ReducedDataset build_reduced_dataset_collapsed(const TreeState& tree,
                                                      const DataSet& ds,
                                                      int sector_root,
                                                      int target_tips) {
  ReducedDataset rd;
  rd.sector_root = sector_root;
  rd.n_htus = 1;
  const int tw = ds.total_words;

  // Frontier: split the clade into ~target_tips sub-clades by repeatedly
  // expanding the frontier node with the largest subtree. Expanded nodes form
  // the skeleton (reduced internals, incl. sector_root); unexpanded frontier
  // nodes are collapsed composite terminals.
  std::vector<int> frontier;   // reduced tips: real tips + collapsed sub-clade roots
  std::vector<int> skeleton;   // reduced internals (expanded nodes)
  frontier.push_back(sector_root);
  while (static_cast<int>(frontier.size()) < target_tips) {
    int best_i = -1, best_sz = 1;
    for (int i = 0; i < static_cast<int>(frontier.size()); ++i) {
      int nd = frontier[i];
      if (nd < tree.n_tip) continue;             // real tip: not expandable
      int sz = count_clade_tips(tree, nd);
      if (sz > best_sz) { best_sz = sz; best_i = i; }
    }
    if (best_i < 0) break;                       // nothing left to expand
    int x = frontier[best_i];
    int xi = x - tree.n_tip;
    frontier[best_i] = tree.left[xi];
    frontier.push_back(tree.right[xi]);
    skeleton.push_back(x);
  }

  const int n_front = static_cast<int>(frontier.size());
  rd.n_real_tips = n_front;                       // reduced tips excluding the HTU
  const int n_sector_tips = n_front + rd.n_htus;  // + HTU
  const int n_sector_internal = n_sector_tips - 1;
  const int n_sector_node = 2 * n_sector_tips - 1;

  rd.full_to_sector.assign(tree.n_node, -1);
  rd.sector_to_full.assign(n_sector_node, -1);

  for (int i = 0; i < n_front; ++i) {             // frontier -> reduced tips
    rd.full_to_sector[frontier[i]] = i;
    rd.sector_to_full[i] = frontier[i];
  }
  const int htu_sector_idx = n_front;
  rd.sector_to_full[htu_sector_idx] = tree.parent[sector_root];

  const int new_root = n_sector_tips;             // synthetic root
  rd.sector_to_full[new_root] = -1;
  int next_internal = new_root + 1;
  for (int s : skeleton) {                        // skeleton -> reduced internals
    rd.full_to_sector[s] = next_internal;
    rd.sector_to_full[next_internal] = s;
    ++next_internal;
  }

  rd.subtree.n_tip = n_sector_tips;
  rd.subtree.n_internal = n_sector_internal;
  rd.subtree.n_node = n_sector_node;
  rd.subtree.total_words = ds.total_words;
  rd.subtree.n_blocks = ds.n_blocks;
  rd.subtree.parent.assign(n_sector_node, -1);
  rd.subtree.left.assign(n_sector_internal, -1);
  rd.subtree.right.assign(n_sector_internal, -1);

  const int sr_mapped = rd.full_to_sector[sector_root];
  const int nr_i = new_root - n_sector_tips;
  rd.subtree.left[nr_i] = htu_sector_idx;
  rd.subtree.right[nr_i] = sr_mapped;
  rd.subtree.parent[new_root] = new_root;
  rd.subtree.parent[htu_sector_idx] = new_root;
  rd.subtree.parent[sr_mapped] = new_root;

  for (int s : skeleton) {                        // map skeleton topology
    const int sec_s = rd.full_to_sector[s];
    const int s_i = s - tree.n_tip;
    const int sec_lc = rd.full_to_sector[tree.left[s_i]];
    const int sec_rc = rd.full_to_sector[tree.right[s_i]];
    rd.subtree.left[sec_s - n_sector_tips] = sec_lc;
    rd.subtree.right[sec_s - n_sector_tips] = sec_rc;
    rd.subtree.parent[sec_lc] = sec_s;
    rd.subtree.parent[sec_rc] = sec_s;
  }

  // Reduced DataSet (mirror build_reduced_dataset's block-structure copy).
  rd.data.n_tips = n_sector_tips;
  rd.data.n_blocks = ds.n_blocks;
  rd.data.total_words = ds.total_words;
  rd.data.blocks = ds.blocks;
  rd.data.block_word_offset = ds.block_word_offset;
  rd.data.flat_blocks = ds.flat_blocks;
  rd.data.all_weight_one = ds.all_weight_one;
  rd.data.n_patterns = ds.n_patterns;
  rd.data.min_steps = ds.min_steps;
  rd.data.pattern_freq = ds.pattern_freq;
  rd.data.concavity = ds.concavity;
  rd.data.eff_k = ds.eff_k;
  rd.data.phi = ds.phi;
  rd.data.scoring_mode = ds.scoring_mode;
  rd.data.ew_offset = ds.ew_offset;
  rd.data.precomputed_steps = ds.precomputed_steps;
  rd.data.info_amounts = ds.info_amounts;
  rd.data.info_max_steps = ds.info_max_steps;
  rd.data.inapp_state = ds.inapp_state;

  const size_t tip_state_size = static_cast<size_t>(n_sector_tips) * tw;
  rd.data.tip_states.assign(tip_state_size, 0ULL);
  for (int i = 0; i < n_front; ++i) {             // composite terminal states
    const int node = frontier[i];
    const size_t dst = static_cast<size_t>(i) * tw;
    // tw == 0 (zero Fitch blocks) leaves ds.tip_states/tree.prelim empty;
    // guard the address-of (UB on an empty vector) -- the copy loop below
    // is already a no-op (w < tw == 0), so src is never dereferenced.
    const uint64_t* src = (tw == 0) ? nullptr
        : (node < tree.n_tip)
        ? &ds.tip_states[static_cast<size_t>(node) * tw]   // real tip
        : &tree.prelim[static_cast<size_t>(node) * tw];    // collapsed sub-clade
    for (int w = 0; w < tw; ++w) rd.data.tip_states[dst + w] = src[w];
  }
  {                                               // HTU: rest-of-tree first pass
    std::vector<uint64_t> from_above_sr;
    compute_from_above_for_sector(tree, ds, sector_root, from_above_sr);
    const size_t dst = static_cast<size_t>(htu_sector_idx) * tw;
    for (int w = 0; w < tw; ++w) rd.data.tip_states[dst + w] = from_above_sr[w];
  }

  const size_t state_size = static_cast<size_t>(n_sector_node) * tw;
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

// Rebuild the sector's content topology from scratch by RAS Wagner: a random
// taxon ORDER with greedy (best-edge) PLACEMENT, keeping the HTU anchored at
// the synthetic root AND present throughout, so placements account for the
// rest-of-tree state it summarises. This is the operation
// TNT performs per sector (3 RAS+TBR restarts): it reaches sector topologies a
// TBR on the *existing* sector subtree cannot, because TBR only locally
// rearranges a tree the global TBR has already converged.
//
// The content is rooted at sr_mapped (= full_to_sector[sector_root]) — exactly
// the node reinsert_sector grafts from — and the synthetic root is set to
// (HTU, sr_mapped). So reinsert_sector and the root-structure check below work
// unchanged. Internal node ids for the content come from the free pool (every
// sector internal except new_root and sr_mapped); the pool size (n_real_tips-2)
// matches the number of insertions exactly.
static void build_ras_sector(ReducedDataset& rd, std::mt19937& rng) {
  TreeState& t = rd.subtree;
  const int n_real = rd.n_real_tips;
  if (n_real < 4) return;  // too small to rebuild; leave existing topology
  const int htu = n_real;                 // htu_sector_idx
  const int n_tip = t.n_tip;
  const int new_root = n_tip;             // synthetic root (= n_sector_tips)
  const int sr_mapped = rd.full_to_sector[rd.sector_root];
  const int tw = t.total_words;

  // Reset topology pointers (tip states stay loaded; only the tree changes).
  t.parent.assign(t.n_node, -1);
  t.left.assign(t.n_internal, -1);
  t.right.assign(t.n_internal, -1);

  // Free pool of internal ids: every internal except new_root and sr_mapped.
  std::vector<int> pool;
  pool.reserve(t.n_internal);
  for (int nd = n_tip; nd < t.n_node; ++nd) {
    if (nd != new_root && nd != sr_mapped) pool.push_back(nd);
  }

  // Random addition order over the real tips (Fisher-Yates).
  std::vector<int> order(n_real);
  for (int i = 0; i < n_real; ++i) order[i] = i;
  for (int i = n_real - 1; i > 0; --i) {
    int j = std::uniform_int_distribution<int>(0, i)(rng);
    std::swap(order[i], order[j]);
  }

  // Seed: anchor the HTU at the synthetic root and put the first two tips
  // under sr_mapped, so the HTU (carrying the rest-of-tree state) is present
  // before any scoring. Then score the seed so prelim/final_ are current for
  // the first greedy placement.
  const int nr_i = new_root - n_tip;
  t.left[nr_i] = htu;
  t.right[nr_i] = sr_mapped;
  t.parent[htu] = new_root;
  t.parent[sr_mapped] = new_root;
  t.parent[new_root] = new_root;
  const int sr_i = sr_mapped - n_tip;
  t.left[sr_i] = order[0];
  t.right[sr_i] = order[1];
  t.parent[order[0]] = sr_mapped;
  t.parent[order[1]] = sr_mapped;
  t.build_postorder();
  fitch_score(t, rd.data);

  // Add the remaining tips by GREEDY placement: each goes at the edge that
  // adds the fewest steps (Wagner), mirroring wagner_tree()'s inner loop.
  // Candidate edges are restricted to the subtree below sr_mapped (never a root
  // edge), so the HTU stays anchored at new_root and the content stays rooted
  // at sr_mapped -- keeping reinsert_sector and the root-structure check valid.
  // Placement uses the EW Fitch proxy even under IW/NA, exactly as wagner_tree;
  // search_sector()'s score_tree() is the authoritative scorer.
  int pool_idx = 0;
  std::vector<int> stack;
  // Exact DIRECTIONAL edge-set scoring (same fix as wagner_tree): the candidate
  // insertion cost is #chars where the tip downpass misses E(A,D) =
  // combine(prelim[D], up[D]) -- NOT the union-of-finals proxy, which undercounts
  // and degrades greedy placement.  Recomputed each step (the tree grows).  The
  // up-message at sub-sr_mapped nodes correctly carries the HTU (rest-of-tree)
  // state, so sector-internal placement accounts for the anchored context.
  std::vector<uint64_t> edge_set;
  // Caller-owned scratch reused across the insertion loop (size-ensured,
  // non-zeroing) so compute_insertion_edge_sets does not reallocate/zero its
  // up-message buffer and preorder list each step.
  std::vector<uint64_t> edge_set_up;
  std::vector<int> edge_set_pre;
  // When there are no Fitch words -- e.g. every character is constant or an
  // autapomorphy, leaving the equal-weights dataset with zero blocks (the
  // same all-uninformative case wagner_tree guards) -- rd.data.tip_states
  // and edge_set are empty and every insertion cost is identically zero.
  // Skip the edge-set precompute and the indirect length evaluation, which
  // would otherwise take the address of element 0 of an empty vector
  // (undefined behaviour; aborts under _GLIBCXX_ASSERTIONS/ASan).  The DFS
  // below still runs, so constraints are honoured and -- all costs being
  // equal -- the first legal edge is chosen (mirrors wagner_tree's fix).
  const bool have_words = rd.data.total_words > 0;
  for (int i = 2; i < n_real; ++i) {
    const int tip = order[i];
    const uint64_t* tip_prelim = have_words
        ? &rd.data.tip_states[static_cast<size_t>(tip) * tw]
        : nullptr;

    if (have_words) {
      compute_insertion_edge_sets(t, rd.data, edge_set, edge_set_up, edge_set_pre);
    }

    int best_above = -1, best_below = -1, best_extra = INT_MAX;
    stack.clear();
    stack.push_back(sr_mapped);
    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();
      if (node < n_tip) continue;  // tip -- no children to enumerate
      int ni = node - n_tip;
      int lc = t.left[ni];
      int rc = t.right[ni];
      if (lc >= 0) {
        int extra = have_words
            ? fitch_indirect_length_cached(
                  tip_prelim, &edge_set[static_cast<size_t>(lc) * tw],
                  rd.data, best_extra)
            : 0;
        if (extra < best_extra) { best_extra = extra; best_above = node; best_below = lc; }
        if (lc >= n_tip) stack.push_back(lc);
      }
      if (rc >= 0) {
        int extra = have_words
            ? fitch_indirect_length_cached(
                  tip_prelim, &edge_set[static_cast<size_t>(rc) * tw],
                  rd.data, best_extra)
            : 0;
        if (extra < best_extra) { best_extra = extra; best_above = node; best_below = rc; }
        if (rc >= n_tip) stack.push_back(rc);
      }
    }

    // sr_mapped always has >= 2 descendant edges once seeded; guard anyway.
    if (best_above < 0 || best_below < 0) {
      best_above = sr_mapped;
      best_below = t.left[sr_i];
    }

    const int new_internal = pool[pool_idx++];
    insert_tip_at_edge(t, tip, new_internal, best_above, best_below);
    wagner_incremental_rescore(t, rd.data, new_internal);
  }

  // Postorder for the subsequent score_tree()/TBR in search_sector().
  t.build_postorder();
}

// ---- TNT-style in-sector drift / combined analysis (godrift/gocomb) ----
//
// Small sectors are cheap to solve exhaustively with RAS+TBR; large sectors have
// more reach but need drift to escape their own local optima. This reuses the
// engine's drift_search() on the reduced sector dataset, PINNING the HTU pseudo-
// tip via a content-clade sector_mask (drift_solve_reduced) so every reduced result
// stays reinsertable -- an UNMASKED drift floats the HTU (re-roots the sector
// against the rest-of-tree summary) ~80% of the time and is then discarded on
// revert, which also loses the sector's TBR polish (net-negative). Gated behind
// sectorGoDrift/sectorGoComb size thresholds. Default OFF; env TS_SECT_DRIFT=0
// force-disables regardless of params (baseline arm).
//
// godrift: solve the sector with a single drift_search per RAS start.
// gocomb : run `combstarts` HTU-anchored RAS+drift starts and keep the best
//          (TNT's combined analysis is RAS+drift+FUSE; the fuse recombination
//          sub-step is DEFERRED — tree_fuse re-roots the reduced tree at tip 0
//          every round, which drops it out of the (HTU, sr_mapped) synthetic-root
//          layout that reinsert_sector requires, so a fused result is not
//          reinsertable without a root-agnostic reinsertion path. Validate
//          godrift first, then add the fuse step behind that work). Every
//          candidate kept here stays HTU-anchored, so reinsert_sector is unchanged.
//
// The env gates are read via FUNCTION-LOCAL statics inside search_sector (like
// _free_htu_probe), NOT namespace-scope statics: a namespace-scope dynamic
// initializer runs at DLL load, before R can Sys.setenv() after library(), so the
// kill switch / stats would silently ignore an env var set from R. Function-local
// statics initialise on the first search_sector call, after any Sys.setenv.
// Engagement instrumentation counters (drift solves / HTU-float reverts). A high
// revert rate means anchored drift is inert -- the signal that motivates pinning
// the HTU (sector_mask). The ++ is unsynchronised across parallel workers, which
// is acceptable for the opt-in serial trajectory-diff run.
static long long g_sect_drift_solves = 0;
static long long g_sect_drift_reverts = 0;

// Solve the reduced sector tree with drift in place; return its best score.
// `content_mask` PINS the HTU: it flags every node in the sector's content clade
// (real tips + content internals) true and the synthetic root, HTU pseudo-tip, and
// content-root node false, so drift's internal clips/regrafts stay strictly below
// sr_mapped and never re-root the reduced tree against the rest-of-tree summary.
// The reduced result is therefore always reinsertable (root_ok holds) -- unlike an
// unmasked drift, whose best landing tree floats the HTU ~80% of the time and is
// then discarded on revert (net-negative, since the fallback lost its TBR polish).
static double drift_solve_reduced(ReducedDataset& rd, const SectorParams& params,
                                  std::function<bool()>& check_timeout,
                                  const std::vector<bool>& content_mask) {
  DriftParams dp;
  dp.n_cycles  = params.sector_drift_cycles;
  dp.afd_limit = params.sector_drift_afd;
  dp.rfd_limit = params.sector_drift_rfd;
  dp.max_hits  = params.internal_max_hits;
  DriftResult dr = drift_search(rd.subtree, rd.data, dp, nullptr, check_timeout,
                                &content_mask);
  return dr.best_score;
}

// Search the reduced dataset and return the best score found.
// Modifies rd.subtree in place, leaving the best sector topology ready for
// reinsertion.
//
// `params.ras_starts` = total starts to try. Start 0 is TBR (or drift) on the
// EXISTING sector subtree (so ras_starts=1 reproduces the prior single-TBR
// behaviour exactly, and the existing topology is always a candidate floor).
// Starts 1.. are HTU-anchored random-addition restarts (RAS+TBR) — TNT's
// per-sector tactic. When the sector is large enough (sectorGoComb/sectorGoDrift),
// the per-start solver escalates to combined analysis or drift respectively.
static double search_sector(ReducedDataset& rd, const SectorParams& params,
                            std::mt19937& rng,
                            std::function<bool()> check_timeout = nullptr) {
  int ras_starts        = params.ras_starts;
  const int max_hits    = params.internal_max_hits;
  const int clip_order  = params.clip_order;
  const bool accept_equal = params.accept_equal;
  if (ras_starts < 1) ras_starts = 1;

  // Env gates, read once on first call (function-local -> after any Sys.setenv;
  // µs-scale ucrt getenv, so not per-pick). TS_SECT_DRIFT=0 force-disables godrift/
  // gocomb regardless of params (baseline arm); TS_SECT_DRIFT_STATS enables the
  // drift-solve / HTU-float-revert instrumentation.
  static const bool kSectDriftForceOff = []{
    const char* e = std::getenv("TS_SECT_DRIFT");
    return e != nullptr && e[0] == '0';
  }();
  static const bool kSectDriftStats =
      std::getenv("TS_SECT_DRIFT_STATS") != nullptr;

  // Solve-mode escalation by sector size (real tips, as count_clade_tips
  // measures max_sector_size). Both godrift and gocomb solve the sector with
  // drift; gocomb additionally raises the number of RAS+drift starts to
  // combstarts and takes precedence when both thresholds fire. (gocomb's fuse
  // recombination sub-step is deferred — see the header note above.)
  const int nrt = rd.n_real_tips;
  const bool comb_mode  = !kSectDriftForceOff &&
                          params.sector_go_comb > 0 && nrt >= params.sector_go_comb;
  const bool drift_mode = !kSectDriftForceOff && !comb_mode &&
                          params.sector_go_drift > 0 && nrt >= params.sector_go_drift;
  const bool use_drift = comb_mode || drift_mode;
  if (comb_mode && params.sector_comb_starts > ras_starts)
    ras_starts = params.sector_comb_starts;
  const int htu_idx = rd.n_real_tips;
  const int root = rd.subtree.n_tip;
  const int sr_mapped = rd.full_to_sector[rd.sector_root];
  const int root_i = root - rd.subtree.n_tip;

  // HTU-pinning mask for the drift path, rebuilt per start below (build_ras_sector
  // reassigns content-internal ids each restart, and sr_mapped's direct children
  // differ each restart).
  std::vector<bool> content_mask;

  double best_score = 0.0;
  bool have_best = false;
  std::vector<int> best_left, best_right, best_parent;

  // search_sector runs once per sector pick (1000s of times/search). std::getenv
  // is µs-scale on Windows/ucrt (linear env scan), so the per-pick D1-probe gate
  // is read ONCE into a static, not per pick/per start (T-S6c micro-bank).
  static const bool _free_htu_probe = std::getenv("TS_FREE_HTU_PROBE") != nullptr;

  // D1 confirm (env TS_FREE_HTU_PROBE): T0 sector's reduced length, baseline for
  // the floating-HTU free re-solve reported after the loop.
  double probe_orig = _free_htu_probe ? score_tree(rd.subtree, rd.data) : 0.0;

  for (int s = 0; s < ras_starts; ++s) {
    if (s > 0) {
      build_ras_sector(rd, rng);
      rd.subtree.build_postorder();
    }

    // Rebuild the HTU-pinning mask for this start's topology. Flag the content
    // clade true EXCEPT the synthetic root, the HTU pseudo-tip, the content root
    // sr_mapped, AND sr_mapped's two direct children: clipping a direct child of
    // sr_mapped splices sr_mapped out of the backbone (floats the HTU), so those
    // two nodes are barred from clips/regrafts. Everything deeper still rearranges
    // freely, including moves between the two halves. This cuts the HTU-float
    // (root_ok=false -> reverted) rate from ~80% (unmasked) to ~10-15%: the
    // residual is drift MID-SEARCH splicing a grandchild so that a fresh, still-
    // masked node becomes sr_mapped's child, which a per-start static mask can't
    // pre-exclude. Fully closing it needs a dynamic parent-in-mask clip gate in
    // tbr_search/drift_phase, whose mask is shared with css_search (clade root IN
    // its mask) -> deferred to avoid changing css semantics. Residual reverts are
    // safe (fall back to the anchored pre-drift topology; the root_ok check below
    // still guards reinsertion).
    if (use_drift) {
      content_mask.assign(rd.subtree.n_node, true);
      content_mask[htu_idx] = false;
      content_mask[root] = false;                     // synthetic root (new_root)
      content_mask[sr_mapped] = false;                // content root
      const int sr_i = sr_mapped - rd.subtree.n_tip;
      content_mask[rd.subtree.left[sr_i]] = false;    // sr_mapped's children:
      content_mask[rd.subtree.right[sr_i]] = false;   // barred to keep sr_mapped
    }

    // Snapshot the (valid, HTU-anchored) pre-TBR topology for revert.
    auto save_left = rd.subtree.left;
    auto save_right = rd.subtree.right;
    auto save_parent = rd.subtree.parent;

    double original_score = score_tree(rd.subtree, rd.data);

    // Solve the sector. Large sectors (>= sectorGoDrift/sectorGoComb) escalate to
    // drift to escape their own local optima (TNT `godrift`); otherwise RAS+TBR.
    double solved_score;
    if (use_drift) {
      solved_score = drift_solve_reduced(rd, params, check_timeout, content_mask);
    } else {
      TBRParams tp;
      tp.max_hits = max_hits;
      tp.clip_order = static_cast<ClipOrder>(clip_order);
      // Let the sector RE-SOLVE itself walk equal-length plateaus (TNT `equals`:
      // "accept equally good subtrees"), holding up to max_hits (sectorMaxHits)
      // equal trees.  Without this the only equal-move path was a coincidental
      // RAS-restart tie, which never fires on large sectors -> accept_equal inert
      // there.  best_score is unchanged (equal moves never worsen it); only the
      // returned topology differs, so reinsert can take the lateral step.
      tp.accept_equal = accept_equal;
      TBRResult tr = tbr_search(rd.subtree, rd.data, tp);
      solved_score = tr.best_score;
    }

    // Verify root structure: HTU and sector_root_mapped must remain direct
    // children of the synthetic root. TBR/drift can regraft onto root edges,
    // displacing nodes outside the clade — if so the reduced result is
    // unusable for reinsertion; revert to the (valid) pre-solve topology.
    int root_lc = rd.subtree.left[root_i];
    int root_rc = rd.subtree.right[root_i];
    bool root_ok = (root_lc == htu_idx && root_rc == sr_mapped) ||
                   (root_lc == sr_mapped && root_rc == htu_idx);

    // Engagement instrumentation (env TS_SECT_DRIFT_STATS): drift-solve count and
    // HTU-float revert count. A high revert rate means anchored drift is inert.
    if (use_drift && kSectDriftStats) {
      ++g_sect_drift_solves;
      if (!root_ok) ++g_sect_drift_reverts;
      if (g_sect_drift_solves == 1 || (g_sect_drift_solves % 20) == 0)
        REprintf("SECT_DRIFT solves=%lld reverts=%lld (%.1f%%)\n",
                 g_sect_drift_solves, g_sect_drift_reverts,
                 100.0 * g_sect_drift_reverts / g_sect_drift_solves);
    }

    // --- D1 warm test (env TS_FREE_HTU_PROBE): isolate the root_ok revert. ---
    // root_ok=false means the reduced solve's best move FLOATS the HTU (re-roots
    // the sector against the rest of the tree) -- the move discarded at the revert
    // below.  By the reduced = full - const invariance (const = rest-of-tree
    // standalone downpass length, independent of how the sector re-roots),
    // solved_score < original_score with root_ok=false PROVES a strictly
    // shorter FULL tree the anchored sectorial throws away.  GUARD: also reports
    // root_ok, so a null can be told apart from "TBR never floats the HTU"
    // (false-negative).  Run rasStarts=1 -> this is the warm T0 sector start.
    if (_free_htu_probe) {
      REprintf("REVERT sect=%d S=%d s=%d orig=%.0f tbr=%.0f root_ok=%d %s\n",
               rd.sector_root, rd.n_real_tips, s, original_score, solved_score,
               root_ok ? 1 : 0,
               (!root_ok && solved_score < original_score) ? "<<FLOAT-IMPROVES" : "");
    }

    double this_score;
    if (!root_ok) {
      rd.subtree.left = save_left;
      rd.subtree.right = save_right;
      rd.subtree.parent = save_parent;
      rd.subtree.build_postorder();
      this_score = original_score;
    } else {
      this_score = solved_score;
    }

    // Strictly-better always wins. With accept_equal, an equal-length RAS
    // rebuild (s>0) also REPLACES the kept topology -- this is the only way a
    // sector escapes onto a different equal-length arrangement (plateau walk),
    // which iterated sector picks then build a strict improvement from. At the
    // default ras_starts=1 there is no s>0, so this is a guaranteed no-op.
    if (ras_starts == 1) {
      // Single-start fast path (the default): rd.subtree already holds the
      // final (TBR-result or reverted) topology, so the best_* snapshot here
      // and the post-loop restore are a provable no-op round-trip — skip both.
      // reinsert_sector reads only left/right/parent (never postorder), so the
      // post-loop build_postorder is also unneeded. (T-S6c micro-bank.)
      best_score = this_score;
      have_best = true;
    } else {
      bool take = !have_best || this_score < best_score ||
                  (accept_equal && s > 0 && this_score == best_score);
      if (take) {
        best_score = this_score;
        best_left = rd.subtree.left;
        best_right = rd.subtree.right;
        best_parent = rd.subtree.parent;
        have_best = true;
      }
    }
  }

  // Restore the best topology found across starts, ready for reinsertion.
  // (Skipped at ras_starts==1: rd.subtree already holds it — see fast path.)
  if (ras_starts > 1) {
    rd.subtree.left = best_left;
    rd.subtree.right = best_right;
    rd.subtree.parent = best_parent;
    rd.subtree.build_postorder();
  }

  // D1 SCORING-ONLY CONFIRM (env TS_FREE_HTU_PROBE), NO reinsertion: does an
  // UNCONSTRAINED reduced search -- HTU = ordinary floating leaf among rd.data's
  // (S+1) tips -- find a LOWER reduced score than the anchored search (best_score)
  // or T0 (probe_orig)?  By the reduced = full - const invariance, free < anchored
  // PROVES a shorter FULL tree the anchored sectorial cannot reach (audit D1).  20
  // free RAS+TBR restarts so medium sectors reach their true optimum (free >= orig
  // on a LARGE sector may be cold-search weakness -- weigh the medium sectors).
  if (_free_htu_probe) {
    double free_min = HUGE_VAL;
    for (int fs = 0; fs < 20; ++fs) {
      TreeState ft;
      random_wagner_tree(ft, rd.data, nullptr);
      TBRParams ftp;
      ftp.max_hits = max_hits;
      ftp.clip_order = static_cast<ClipOrder>(clip_order);
      TBRResult ftr = tbr_search(ft, rd.data, ftp);
      if (ftr.best_score < free_min) free_min = ftr.best_score;
    }
    REprintf("FREEHTU sect=%d S=%d orig=%.0f anchored=%.0f free=%.0f %s\n",
             rd.sector_root, rd.n_real_tips, probe_orig, best_score, free_min,
             free_min < best_score ? "<<D1-CONFIRM" : "");
  }

  // Return: best score across starts
  return best_score;
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
  // Defensive: `n_partitions` reaches here straight from SearchControl, and a
  // value of 0 would make `tree.n_tip / n_partitions` (below) an integer
  // division by zero -- an uncatchable SIGFPE that kills the R session.
  // Treat any non-positive request as a single partition (the whole tree).
  if (n_partitions < 1) n_partitions = 1;

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
                        ConstraintData* cd,
                        std::function<bool()> check_timeout) {
  // Hoist the per-accept debug-trace gate (µs-scale ucrt getenv) to a static
  // (T-S6c micro-bank).
  static const bool _sect_debug = std::getenv("TS_SECT_DEBUG") != nullptr;
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

  // build_reduced_dataset() does not copy hierarchy_blocks, tip_labels,
  // n_orig_chars, hsj_alpha, or sankoff_* fields (T-303).  Sector-internal
  // scoring would silently degrade to Fitch-only.  Same class as T-275 guard.
  if (ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) {
    return result;
  }

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
    tp.clip_order = static_cast<ClipOrder>(params.clip_order);
    TBRResult tr = tbr_search(tree, ds, tp, nullptr, nullptr, nullptr,
                              check_timeout);
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

    // Build reduced dataset; collapse deep sub-clades into composite terminals
    // when collapse_target is set and the clade is larger (Goloboff 1999).
    int clade_sz_ = count_clade_tips(tree, sector_root);
    ReducedDataset rd =
        (params.collapse_target > 0 && clade_sz_ > params.collapse_target)
        ? build_reduced_dataset_collapsed(tree, ds, sector_root, params.collapse_target)
        : build_reduced_dataset(tree, ds, sector_root);
    const long long sector_cand0 = rd.data.n_candidates_evaluated;

    // Score the current sector topology
    double sector_current = score_tree(rd.subtree, rd.data);

    // Search the sector
    double sector_best = search_sector(rd, params, rng, check_timeout);
    ++result.n_sectors_searched;
    // Propagate reduced-dataset candidates to the parent (diagnostics).
    ds.n_candidates_evaluated += rd.data.n_candidates_evaluated - sector_cand0;

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
      if (_sect_debug)
        REprintf("  sect[%2d] red_cur=%.0f red_best=%.0f full_new=%.0f full_best=%.0f %s\n",
                 sector_root, sector_current, sector_best, new_score, result.best_score,
                 new_score < result.best_score ? "STRICT" :
                 (new_score == result.best_score ? "EQUAL-keep" : "WORSE-revert"));

      // Post-hoc constraint check
      if (constrained && violates_constraint_posthoc(tree, *cd)) {
        restore_clade(tree, snap);
        tree.build_postorder();
        score_tree(tree, ds);
        continue;
      }

      bool kept;
      if (new_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - new_score);
        result.best_score = new_score;
        ++result.n_sectors_improved;
        kept = true;
      } else if (new_score == result.best_score && params.accept_equal) {
        // Equal-length lateral move accepted (plateau walk): topology changed
        // but score did not.  This MUST refresh subtree_size / eligible just
        // like a strict improvement (shared block below).  Omitting it left a
        // STALE candidate set for every subsequent pick of the walk, so the
        // plateau walk drew sectors against the pre-move topology and made no
        // headway -- accept_equal was observably inert across the walk.
        kept = true;
      } else {
        // HTU approximation caused full-tree score to worsen; revert
        restore_clade(tree, snap);
        tree.build_postorder();
        score_tree(tree, ds);
        kept = false;
      }

      if (kept) {
        // Recompute subtree sizes and eligible list after the topology change
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
        // Re-sync constraint metadata to the updated topology.  The
        // global TBR cleanup at the end of rss_search passes `cd`
        // directly; stale constraint_node mapping after a sector
        // change would cause false-positive or false-negative
        // constraint violations for the first TBR clips.
        if (cd && cd->active) {
          map_constraint_nodes(tree, *cd);
          compute_dfs_timestamps(tree, *cd);
        }
      }
    }

    if (ts::check_interrupt() || (check_timeout && check_timeout())) break;
  }

  // Global TBR after all sector picks
  {
    TBRParams tp;
    tp.max_hits = params.internal_max_hits;
    tp.clip_order = static_cast<ClipOrder>(params.clip_order);
    TBRResult tr = tbr_search(tree, ds, tp, cd, nullptr, nullptr,
                              check_timeout);
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
                        ConstraintData* cd,
                        std::function<bool()> check_timeout) {
  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  double current_score = score_tree(tree, ds);

  SectorResult result;
  result.best_score = current_score;
  result.n_sectors_searched = 0;
  result.n_sectors_improved = 0;
  result.total_steps_saved = 0;

  // build_reduced_dataset() does not copy hierarchy_blocks, tip_labels,
  // n_orig_chars, hsj_alpha, or sankoff_* fields (T-303).  Sector-internal
  // scoring would silently degrade to Fitch-only.  Same class as T-275 guard.
  if (ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) {
    return result;
  }

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

      int clade_sz_ = count_clade_tips(tree, sector_root);
      ReducedDataset rd =
          (params.collapse_target > 0 && clade_sz_ > params.collapse_target)
          ? build_reduced_dataset_collapsed(tree, ds, sector_root, params.collapse_target)
          : build_reduced_dataset(tree, ds, sector_root);
      const long long sector_cand0 = rd.data.n_candidates_evaluated;

      double sector_current = score_tree(rd.subtree, rd.data);
      double sector_best = search_sector(rd, params, rng, check_timeout);
      ++result.n_sectors_searched;
      // Propagate reduced-dataset candidates to the parent (diagnostics).
      ds.n_candidates_evaluated += rd.data.n_candidates_evaluated - sector_cand0;

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
          // Re-sync constraint metadata to the updated topology.
          // Without this, the global TBR cleanup at the end of each
          // XSS round uses stale constraint_node mapping (same class
          // as T-278 / T-279 / T-280).
          if (cd && cd->active) {
            map_constraint_nodes(tree, *cd);
            compute_dfs_timestamps(tree, *cd);
          }
        } else if (new_score == result.best_score && params.accept_equal) {
          // Equal score accepted — topology changed; re-sync constraint.
          if (cd && cd->active) {
            map_constraint_nodes(tree, *cd);
            compute_dfs_timestamps(tree, *cd);
          }
        } else {
          // HTU approximation caused full-tree score to worsen; revert
          restore_clade(tree, snap);
          tree.build_postorder();
          score_tree(tree, ds);
        }
      }

      if (ts::check_interrupt() || (check_timeout && check_timeout())) break;
    }

    // Global TBR after each round of sectors
    {
      TBRParams tp;
      tp.max_hits = params.internal_max_hits;
      tp.clip_order = static_cast<ClipOrder>(params.clip_order);
      TBRResult tr = tbr_search(tree, ds, tp, cd, nullptr, nullptr,
                                check_timeout);
      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
      }
    }

    // Adaptive: skip remaining rounds if this one found no improvement
    if (result.best_score >= score_before_round) break;
    if (ts::check_interrupt() || (check_timeout && check_timeout())) break;
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
                        ConstraintData* cd,
                        std::function<bool()> check_timeout) {
  // No HSJ/XFORM guard needed here (cf. T-303 guards in rss_search/xss_search).
  // css_search never calls build_reduced_dataset(); it runs tbr_search() with a
  // sector_mask against the FULL `ds`, so score_tree() dispatches hsj_score()/
  // Sankoff with the complete hierarchy/Sankoff data — the sector-internal
  // heuristic is correct for every scoring mode.
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
      tp.clip_order = static_cast<ClipOrder>(params.clip_order);

      TBRResult tr = tbr_search(tree, ds, tp, cd, &sector_mask, nullptr,
                                check_timeout);
      ++result.n_sectors_searched;

      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
        ++result.n_sectors_improved;
      }

      if (ts::check_interrupt() || (check_timeout && check_timeout())) break;
    }

    // Global TBR after each round
    {
      TBRParams tp;
      tp.max_hits = params.internal_max_hits;
      tp.clip_order = static_cast<ClipOrder>(params.clip_order);
      TBRResult tr = tbr_search(tree, ds, tp, cd, nullptr, nullptr,
                                check_timeout);
      if (tr.best_score < result.best_score) {
        result.total_steps_saved +=
            static_cast<int>(result.best_score - tr.best_score);
        result.best_score = tr.best_score;
      }
    }

    // Adaptive: skip remaining rounds if this one found no improvement
    if (result.best_score >= score_before_round) break;
    if (ts::check_interrupt() || (check_timeout && check_timeout())) break;
  }

  return result;
}

} // namespace ts
