#include "ts_fuse.h"
#include "ts_fitch.h"
#include "ts_splits.h"
#include "ts_tbr.h"
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <vector>

namespace ts {

// ---------- Internal helpers ----------

// Re-root the tree so that tip 0 is a direct child of the root.
// This ensures every non-trivial split has was_flipped=false (tip 0 is
// never inside any non-root subtree), making split matching between
// differently-rooted trees consistent.
// Parsimony scores are rooting-invariant, so this is safe.
static void reroot_at_tip0(TreeState& tree) {
  int n_tip = tree.n_tip;
  int root = n_tip;

  if (tree.parent[0] == root) return;  // already a child of root

  // Collect path of internal nodes from tip 0's parent to root
  // (root-to-tip order after reversal).
  std::vector<int> path;  // tip-to-root order initially
  int cur = tree.parent[0];
  while (cur != root) {
    path.push_back(cur);
    cur = tree.parent[cur];
  }
  // path = [parent_of_0, ..., child_of_root_on_path], tip-to-root order
  // Reverse to root-to-tip order for processing
  std::reverse(path.begin(), path.end());
  // path = [child_of_root_on_path, ..., parent_of_0]

  int path_len = static_cast<int>(path.size());
  int root_ni = 0;  // root - n_tip

  // Sibling of path[0] under root (this gets "absorbed" by path[0])
  int root_other = (tree.left[root_ni] == path[0])
                       ? tree.right[root_ni]
                       : tree.left[root_ni];

  // Process each node on the path.  For each node, replace its child
  // that is "toward tip 0" with a new child from the root side:
  //   - path[0] absorbs root's other child
  //   - path[i>0] absorbs path[i-1] (the preceding node, now detached)
  for (int i = 0; i < path_len; ++i) {
    int node = path[i];
    int ni = node - n_tip;

    // The child of `node` that is toward tip 0:
    int toward_tip0 = (i + 1 < path_len) ? path[i + 1] : 0;

    // The replacement child coming from the root side:
    int replacement = (i == 0) ? root_other : path[i - 1];

    if (tree.left[ni] == toward_tip0) {
      tree.left[ni] = replacement;
    } else {
      tree.right[ni] = replacement;
    }
    tree.parent[replacement] = node;
  }

  // Finally, update root: children are tip 0 and the last path node
  int last_path = path[path_len - 1];
  tree.left[root_ni] = 0;
  tree.right[root_ni] = last_path;
  tree.parent[0] = root;
  tree.parent[last_path] = root;

  tree.build_postorder();
}

// Compute tip-membership bitsets for every node in the tree.
// Returns a flat vector of size n_node * words_per_split.
static std::vector<uint64_t> compute_tip_bits(const TreeState& tree) {
  int n_tip = tree.n_tip;
  int wps = (n_tip + 63) / 64;
  size_t total = static_cast<size_t>(tree.n_node) * wps;
  std::vector<uint64_t> bits(total, 0);

  // Tips: tip i has bit i set
  for (int t = 0; t < n_tip; ++t) {
    int word = t / 64;
    int bit  = t % 64;
    bits[static_cast<size_t>(t) * wps + word] = 1ULL << bit;
  }

  // Internal: union of children, postorder
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* dst = &bits[static_cast<size_t>(node) * wps];
    const uint64_t* lb = &bits[static_cast<size_t>(lc) * wps];
    const uint64_t* rb = &bits[static_cast<size_t>(rc) * wps];
    for (int w = 0; w < wps; ++w) {
      dst[w] = lb[w] | rb[w];
    }
  }
  return bits;
}

// Count the number of set bits (tips) in a split/bitset.
static int count_tips(const uint64_t* bits, int wps) {
  int c = 0;
  for (int w = 0; w < wps; ++w) {
    c += popcount64(bits[w]);
  }
  return c;
}

// Canonicalize a split: ensure bit 0 is unset (tip 0 in partition 0).
// Mask trailing bits in the last word.
static void canonicalize(uint64_t* s, int wps, int n_tips) {
  if (s[0] & 1ULL) {
    for (int w = 0; w < wps; ++w) s[w] = ~s[w];
  }
  int trailing = n_tips % 64;
  if (trailing != 0) {
    s[wps - 1] &= (1ULL << trailing) - 1;
  }
}

// Compare two split bitsets for equality.
static bool split_eq(const uint64_t* a, const uint64_t* b, int wps) {
  return std::memcmp(a, b, sizeof(uint64_t) * wps) == 0;
}

// A hash for a single split bitset (for use in hash maps).
static uint64_t hash_one_split(const uint64_t* s, int wps) {
  uint64_t h = 0;
  for (int w = 0; w < wps; ++w) {
    h ^= s[w] * (0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(w) * 0x517cc1b727220a95ULL);
  }
  // splitmix64 mix
  h ^= h >> 30;
  h *= 0xbf58476d1ce4e5b9ULL;
  h ^= h >> 27;
  h *= 0x94d049bb133111ebULL;
  h ^= h >> 31;
  return h;
}

// Info about one non-trivial split in a tree, including which node roots
// the corresponding clade.
struct SplitInfo {
  int node;             // internal node whose subtree defines this split
  int clade_size;       // number of tips in the clade
  std::vector<uint64_t> canonical;  // canonicalized split bitset (wps words)
  std::vector<uint64_t> raw;       // raw (uncanonicalized) tip bitset
  bool was_flipped;     // true if canonical != raw (tip 0 was in subtree)
  uint64_t hash;        // hash of the canonical split
};

// Build SplitInfo for all non-trivial, non-root splits in a tree.
static std::vector<SplitInfo> build_split_info(
    const TreeState& tree,
    const std::vector<uint64_t>& tip_bits)
{
  int n_tip = tree.n_tip;
  int wps = (n_tip + 63) / 64;
  int root = n_tip;
  int root_right = tree.right[0];

  std::vector<SplitInfo> infos;
  infos.reserve(n_tip - 3);

  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    if (node == root || node == root_right) continue;
    if (node < n_tip) continue;  // tips are not in postorder, but guard anyway

    const uint64_t* bits = &tip_bits[static_cast<size_t>(node) * wps];
    int sz = count_tips(bits, wps);
    if (sz <= 1 || sz >= n_tip - 1) continue;  // trivial split

    SplitInfo si;
    si.node = node;
    si.clade_size = sz;
    si.raw.assign(bits, bits + wps);
    si.was_flipped = (bits[0] & 1ULL) != 0;  // tip 0 in subtree → will flip
    si.canonical.assign(bits, bits + wps);
    canonicalize(si.canonical.data(), wps, n_tip);
    si.hash = hash_one_split(si.canonical.data(), wps);
    infos.push_back(std::move(si));
  }
  return infos;
}

// Collect all nodes in the subtree rooted at `root_node`.
// Separates tips and internal nodes.
static void collect_clade(const TreeState& tree, int root_node,
                          std::vector<int>& internals) {
  // DFS
  std::vector<int> stk;
  stk.push_back(root_node);
  while (!stk.empty()) {
    int node = stk.back();
    stk.pop_back();
    if (node < tree.n_tip) {
      // tip — no children
    } else {
      internals.push_back(node);
      int ni = node - tree.n_tip;
      stk.push_back(tree.left[ni]);
      stk.push_back(tree.right[ni]);
    }
  }
}

// Replace the subtree in `recipient` at `r_root` with the topology from
// `donor` at `d_root`. Both clades contain the same set of tips.
// After this call, recipient's topology is updated but state arrays
// (prelim, etc.) are NOT recomputed — caller must do that.
static void replace_subtree(TreeState& recipient, const TreeState& donor,
                            int r_root, int d_root) {
  int n_tip = recipient.n_tip;

  // Collect internal nodes in each clade
  std::vector<int> r_internals, d_internals;
  collect_clade(recipient, r_root, r_internals);
  collect_clade(donor, d_root, d_internals);

  // Build bijection: donor node → recipient node.
  // d_root must map to r_root; other internals get arbitrary pairing.
  // Remove roots from lists, pair explicitly, then pair the rest.
  //
  // NB: The pairing is topology-unaware — node indices get shuffled.
  // This is safe because we copy left/right through the mapping, so the
  // tree structure is preserved. However, any code that caches per-node
  // identity (e.g. for incremental scoring) must not assume node IDs are
  // stable across replace_subtree calls.
  std::vector<int> r_rest, d_rest;
  for (int n : r_internals) {
    if (n != r_root) r_rest.push_back(n);
  }
  for (int n : d_internals) {
    if (n != d_root) d_rest.push_back(n);
  }
  std::sort(r_rest.begin(), r_rest.end());
  std::sort(d_rest.begin(), d_rest.end());

  // Map: donor node → recipient node
  // Tips map to themselves.
  std::unordered_map<int, int> dtor;
  dtor[d_root] = r_root;
  for (size_t i = 0; i < d_rest.size(); ++i) {
    dtor[d_rest[i]] = r_rest[i];
  }

  // Copy topology from donor to recipient using the mapping
  for (size_t i = 0; i < d_internals.size(); ++i) {
    int d_node = d_internals[i];
    int r_node = dtor[d_node];
    int d_ni = d_node - n_tip;
    int r_ni = r_node - n_tip;

    int d_lc = donor.left[d_ni];
    int d_rc = donor.right[d_ni];

    // Map children: tips map to themselves, internals through dtor
    int r_lc = (d_lc < n_tip) ? d_lc : dtor[d_lc];
    int r_rc = (d_rc < n_tip) ? d_rc : dtor[d_rc];

    recipient.left[r_ni]  = r_lc;
    recipient.right[r_ni] = r_rc;
    recipient.parent[r_lc] = r_node;
    recipient.parent[r_rc] = r_node;
  }
  // r_root's parent stays unchanged (it's outside the clade)
}

// Create a topology-only copy of a TreeState with fresh state arrays.
// Tip states are loaded from the dataset; internal state arrays are zeroed.
// Avoids the cost of copying large prelim/final_/local_cost arrays that
// fitch_score will overwrite anyway.
static TreeState copy_topology(const TreeState& src, const DataSet& ds) {
  TreeState t;
  t.n_tip = src.n_tip;
  t.n_internal = src.n_internal;
  t.n_node = src.n_node;
  t.total_words = src.total_words;
  t.n_blocks = src.n_blocks;
  t.parent = src.parent;
  t.left = src.left;
  t.right = src.right;
  t.postorder = src.postorder;
  size_t state_size = static_cast<size_t>(t.n_node) * t.total_words;
  t.prelim.resize(state_size, 0ULL);
  t.final_.resize(state_size, 0ULL);
  t.down2.resize(state_size, 0ULL);
  t.subtree_actives.resize(state_size, 0ULL);
  t.local_cost.resize(static_cast<size_t>(t.n_node) * t.n_blocks, 0ULL);
  t.load_tip_states(ds);
  return t;
}

// Check if split `ancestor` is a strict superset of split `descendant`.
static bool is_ancestor_split(const uint64_t* ancestor,
                              const uint64_t* descendant, int wps) {
  // ancestor ⊃ descendant iff (descendant & ~ancestor) == 0
  // and they are not equal.
  bool is_superset = true;
  bool is_equal = true;
  for (int w = 0; w < wps; ++w) {
    if (descendant[w] & ~ancestor[w]) {
      is_superset = false;
      break;
    }
    if (ancestor[w] != descendant[w]) {
      is_equal = false;
    }
  }
  return is_superset && !is_equal;
}

// ---------- Main fuse algorithm ----------

FuseResult tree_fuse(TreeState& recipient, const DataSet& ds,
                     const TreePool& pool, const FuseParams& params) {
  FuseResult result;
  result.n_exchanges = 0;
  result.n_rounds = 0;

  // Re-root recipient at tip 0 for consistent split orientation.
  reroot_at_tip0(recipient);

  // Initial score
  double score = static_cast<double>(fitch_score(recipient, ds));

  const auto& entries = pool.all();
  int wps = (recipient.n_tip + 63) / 64;

  // Create tip-0-rooted copies of donor trees.
  // We need the re-rooted topology both for split computation (node → split
  // mapping) and for replace_subtree (clade topology extraction).
  std::vector<TreeState> donor_trees(entries.size());
  std::vector<std::vector<SplitInfo>> donor_splits(entries.size());
  for (int di = 0; di < static_cast<int>(entries.size()); ++di) {
    donor_trees[di] = copy_topology(entries[di].tree, ds);
    reroot_at_tip0(donor_trees[di]);
    std::vector<uint64_t> d_tip_bits = compute_tip_bits(donor_trees[di]);
    donor_splits[di] = build_split_info(donor_trees[di], d_tip_bits);
  }

  bool improved = true;
  while (improved && result.n_rounds < params.max_rounds) {
    improved = false;
    ++result.n_rounds;

    // Compute recipient's tip bits and split info (changes each round)
    std::vector<uint64_t> r_tip_bits = compute_tip_bits(recipient);
    std::vector<SplitInfo> r_splits = build_split_info(recipient, r_tip_bits);

    // Build a hash map from recipient split hash → index in r_splits
    std::unordered_multimap<uint64_t, int> r_split_map;
    for (int i = 0; i < static_cast<int>(r_splits.size()); ++i) {
      r_split_map.emplace(r_splits[i].hash, i);
    }

    for (int di = 0; di < static_cast<int>(entries.size()); ++di) {
      const TreeState& donor = donor_trees[di];
      const std::vector<SplitInfo>& d_splits = donor_splits[di];

      // Find shared splits: for each donor split, check if the recipient
      // has a matching one.
      struct SharedSplit {
        int r_node;      // clade root in recipient
        int d_node;      // clade root in donor
        int clade_size;  // number of tips
        int r_idx;       // index in r_splits (for ancestor checking)
      };

      std::vector<SharedSplit> shared;

      for (const auto& ds_info : d_splits) {
        auto range = r_split_map.equal_range(ds_info.hash);
        for (auto it = range.first; it != range.second; ++it) {
          int ri = it->second;
          if (split_eq(ds_info.canonical.data(),
                       r_splits[ri].canonical.data(), wps)) {
            // Both trees are re-rooted at tip 0, so all non-trivial
            // splits have was_flipped=false (tip 0 is always outside
            // the subtree). No flip-status filtering needed.

            SharedSplit ss;
            ss.r_node = r_splits[ri].node;
            ss.d_node = ds_info.node;
            ss.clade_size = ds_info.clade_size;
            ss.r_idx = ri;
            shared.push_back(ss);
            break;  // each donor split matches at most one recipient split
          }
        }
      }

      if (shared.empty()) continue;

      // Sort by clade size ascending (try smallest first)
      std::sort(shared.begin(), shared.end(),
                [](const SharedSplit& a, const SharedSplit& b) {
                  return a.clade_size < b.clade_size;
                });

      // Track which splits are "stale" (ancestor of an applied exchange)
      std::vector<bool> stale(shared.size(), false);

      for (size_t si = 0; si < shared.size(); ++si) {
        if (stale[si]) continue;

        const SharedSplit& ss = shared[si];

        // Trial: lightweight topology copy, replace subtree, rescore
        TreeState trial = copy_topology(recipient, ds);
        replace_subtree(trial, donor, ss.r_node, ss.d_node);
        trial.build_postorder();

        // Copy donor's tip prelim states for the tips in the clade
        // (tips haven't changed, so prelim for tips should still be correct
        //  from the recipient copy — they are the same tips)
        double new_score = static_cast<double>(fitch_score(trial, ds));

        bool accept = false;
        if (new_score < score) {
          accept = true;
        } else if (params.accept_equal && new_score == score) {
          // Check that the topology actually changed within the clade
          // (avoid infinite loops on identical clades).
          // Cheap check: compare left/right arrays — if any entry differs,
          // the topology changed.
          bool changed = false;
          for (int ni = 0; ni < recipient.n_internal && !changed; ++ni) {
            if (trial.left[ni] != recipient.left[ni] ||
                trial.right[ni] != recipient.right[ni]) {
              changed = true;
            }
          }
          if (changed) {
            accept = true;
          }
        }

        if (accept) {
          recipient = trial;
          score = new_score;
          ++result.n_exchanges;
          improved = true;

          // Mark ancestor splits as stale.
          // Must use raw (uncanonicalized) bitsets: canonical splits are
          // flipped when tip 0 is in the clade, which reverses subset
          // relationships and breaks ancestor detection.
          const uint64_t* exchanged_raw =
              r_splits[ss.r_idx].raw.data();
          for (size_t sj = si + 1; sj < shared.size(); ++sj) {
            if (stale[sj]) continue;
            if (is_ancestor_split(
                    r_splits[shared[sj].r_idx].raw.data(),
                    exchanged_raw, wps)) {
              stale[sj] = true;
            }
          }

          // After modifying the recipient, the r_splits and r_split_map
          // are stale for this donor. Break out and move to TBR cleanup.
          // (We'll recompute everything in the next round.)
          break;
        }
      }

      // If we found an improvement with this donor, break the donor loop
      // to run TBR and start a fresh round.
      if (improved) break;
    }

    if (improved) {
      // TBR search to clean up
      TBRParams tbr_params;
      tbr_params.accept_equal = false;
      tbr_params.max_hits = 1;
      TBRResult tbr_res = tbr_search(recipient, ds, tbr_params);
      score = tbr_res.best_score;
    }
  }

  result.best_score = score;
  return result;
}

} // namespace ts
