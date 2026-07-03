#include "ts_collapsed.h"
#include "ts_simd.h"
#include <cstring>

namespace ts {

void compute_collapsed_flags(
    const TreeState& tree,
    const DataSet& ds,
    std::vector<uint8_t>& collapsed) {

  collapsed.assign(tree.n_node, 0);

  // HSJ/XFORM store a clade's character support in ds.hierarchy_blocks (HSJ) or
  // ds.sankoff_* (XFORM); the hierarchy characters are zero-weighted out of the
  // standard ds.blocks[] (.NonHierarchyWeights) and dropped by simplify_patterns.
  // This kernel decides min-length-0 edges by iterating ONLY ds.blocks[], so a
  // clade supported solely by a hierarchy/Sankoff character looks unsupported
  // and would be wrongly contracted.  Until the collapse test folds in the HSJ
  // a(n)/p(n) DP or the Sankoff cost, disable collapse entirely for these modes
  // (all-zero flags == nothing collapses — the safe conservative outcome).
  // Falling back to the conservative flags is NOT sufficient: it is equally
  // blind to hierarchy/Sankoff support.  See red-team T-330.
  if (ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) return;

  // If all characters were simplified away (total_words == 0), every binary
  // resolution ties at the same score: no internal branch carries support,
  // so the correct collapsed answer is a star. Flag every internal,
  // non-root node's edge as collapsible (terminal/pendant edges never
  // collapse) — including root's children, since without them the tree
  // would stay partially resolved instead of collapsing to the star.
  if (tree.total_words == 0) {
    for (int c = tree.n_tip + 1; c < tree.n_node; ++c) collapsed[c] = 1;
    return;
  }

  // Detect whether any block has inapplicable characters.
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  const int nb = ds.n_blocks;
  const int tw = tree.total_words;
  const size_t word_bytes = static_cast<size_t>(tw) * sizeof(uint64_t);

  for (int c = 0; c < tree.n_node; ++c) {
    // Skip root (no parent edge).
    if (c == tree.n_tip) continue;
    int p = tree.parent[c];
    // Root's children can be clipped but removing one changes the root
    // structure; skip for safety.
    if (p == tree.n_tip) continue;

    int pi = p - tree.n_tip;  // internal index for left[]/right[]
    int s = (tree.left[pi] == c) ? tree.right[pi] : tree.left[pi];

    // --- Condition 1: zero standard-block cost at parent ---
    bool zero_std = true;
    for (int b = 0; b < nb && zero_std; ++b) {
      if (ds.blocks[b].has_inapplicable) continue;
      if (tree.local_cost[static_cast<size_t>(p) * nb + b])
        zero_std = false;
    }
    if (!zero_std) continue;

    // --- Condition 2: zero NA-block cost at parent ---
    if (has_na) {
      bool zero_na = true;
      int left = tree.left[pi];
      int right = tree.right[pi];
      size_t lb = static_cast<size_t>(left) * tw;
      size_t rb = static_cast<size_t>(right) * tw;
      size_t pb = static_cast<size_t>(p) * tw;

      for (int b = 0; b < nb && zero_na; ++b) {
        const CharBlock& blk = ds.blocks[b];
        if (!blk.has_inapplicable) continue;

        int k = blk.n_states;
        int off = ds.block_word_offset[b];

        // l_act = OR of applicable state words (states 1..k-1) for left
        uint64_t l_act = simd::or_reduce(
            &tree.subtree_actives[lb + off], k, 1);
        uint64_t r_act = simd::or_reduce(
            &tree.subtree_actives[rb + off], k, 1);

        if ((l_act & r_act) == 0) continue;  // auto-zero

        // Both subtrees have applicable tips — compute full condition.
        // ss_app = OR of applicable states at node p
        uint64_t ss_app = 0;
        for (int st = 1; st < k; ++st)
          ss_app |= tree.final_[pb + off + st];

        // any_isect = any D2 state intersection between children
        uint64_t any_isect = simd::any_hit_reduce(
            &tree.down2[lb + off], &tree.down2[rb + off], k);

        uint64_t needs_step =
            l_act & r_act & ~(ss_app & any_isect) & blk.active_mask;
        if (needs_step) zero_na = false;
      }
      if (!zero_na) continue;
    }

    // --- Condition 3: prelim[sibling] == prelim[parent] ---
    size_t sb = static_cast<size_t>(s) * tw;
    size_t pb = static_cast<size_t>(p) * tw;
    if (std::memcmp(&tree.prelim[sb], &tree.prelim[pb], word_bytes) != 0)
      continue;

    // --- Conditions 4–5 (NA only): down2 and subtree_actives preservation ---
    if (has_na) {
      if (std::memcmp(&tree.down2[sb], &tree.down2[pb], word_bytes) != 0)
        continue;
      if (std::memcmp(&tree.subtree_actives[sb],
                      &tree.subtree_actives[pb], word_bytes) != 0)
        continue;
    }

    collapsed[c] = 1;
  }
}

void compute_collapsed_flags_aggressive(
    const TreeState& tree,
    const DataSet& ds,
    std::vector<uint8_t>& collapsed) {

  // HSJ/XFORM: character support lives in ds.hierarchy_blocks / ds.sankoff_*,
  // which the marginal-MPR machinery below never reads (it iterates only
  // ds.blocks[]).  A clade supported solely by a hierarchy/Sankoff character
  // would be treated as unsupported and contracted.  Disable collapse for these
  // modes (all-zero flags).  Guarded independently of compute_collapsed_flags:
  // the has_na delegation at the bottom of this block only reaches it on NA
  // data, not the general HSJ/XFORM case.  See red-team T-330.
  if (ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) {
    collapsed.assign(tree.n_node, 0);
    return;
  }

  // Inapplicable characters: soft min-length-0 for the De-Laet 3-pass is not
  // derived; fall back to the conservative (exact) flags so NA datasets are
  // unaffected by the aggressive flag.
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }
  if (has_na) { compute_collapsed_flags(tree, ds, collapsed); return; }

  collapsed.assign(tree.n_node, 0);
  if (tree.total_words == 0) {
    for (int c = tree.n_tip + 1; c < tree.n_node; ++c) collapsed[c] = 1;
    return;
  }

  const int nb = ds.n_blocks;
  const int tw = tree.total_words;

  // We need the TRUE marginal MPR state sets (the set of states each node takes
  // in SOME most-parsimonious reconstruction): edge (p,c) has minimum possible
  // length 0 iff mpr[p] & mpr[c] != 0 for every character (criterion validated
  // bit-for-bit against a brute-force oracle, b2_minlength_oracle.R).  The
  // kernel's tree.final_ is a SIMPLER single-reconstruction final set (one
  // valid MPR, not the full marginal set), so we cannot use it here — we
  // recompute the marginal MPR sets from prelim + local_cost (both maintained
  // for scoring) via the Swofford & Maddison second pass:
  //   final[root] = prelim[root]
  //   for v (preorder, parent p, children L,R):
  //     if final[p] subset of prelim[v]:           final[v] = final[p]
  //     elif v formed by INTERSECTION (not union):  final[v] = prelim[v] | (final[p] & (prelim[L]|prelim[R]))
  //     else (v formed by union):                   final[v] = prelim[v] | final[p]
  // local_cost[v][b] holds the per-character `needs_union` mask, so a character
  // is an intersection node iff its bit is NOT set there.
  std::vector<uint64_t> mpr(static_cast<size_t>(tree.n_node) * tw);
  // tips & root: marginal MPR set == prelim (tips fixed; root's prelim is final)
  for (int v = 0; v < tree.n_node; ++v) {
    if (v >= tree.n_tip && v != tree.n_tip) continue;  // internal non-root: filled below
    size_t vb = static_cast<size_t>(v) * tw;
    std::memcpy(&mpr[vb], &tree.prelim[vb], static_cast<size_t>(tw) * sizeof(uint64_t));
  }
  // Preorder (root -> leaves) over internal non-root nodes; parent's mpr ready.
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int v = tree.postorder[i];
    if (v < tree.n_tip || v == tree.n_tip) continue;   // skip tips and root
    int p = tree.parent[v];
    int ni = v - tree.n_tip;
    int lc = tree.left[ni], rc = tree.right[ni];
    size_t vb = static_cast<size_t>(v) * tw, pb = static_cast<size_t>(p) * tw;
    size_t lb = static_cast<size_t>(lc) * tw, rb = static_cast<size_t>(rc) * tw;
    for (int b = 0; b < nb; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int off = ds.block_word_offset[b];
      int k = blk.n_states;
      uint64_t am = blk.active_mask;
      uint64_t needs_union =
          tree.local_cost[static_cast<size_t>(v) * nb + b] & am;
      // subset mask: chars where final[p] is a subset of prelim[v]
      uint64_t notsub = 0;
      for (int s = 0; s < k; ++s)
        notsub |= mpr[pb + off + s] & ~tree.prelim[vb + off + s];
      uint64_t subset = ~notsub & am;
      uint64_t intNonsub = (~needs_union) & (~subset) & am;  // intersection node, not subset
      uint64_t uniNonsub = needs_union & (~subset) & am;     // union node, not subset
      for (int s = 0; s < k; ++s) {
        uint64_t fp = mpr[pb + off + s];
        uint64_t pv = tree.prelim[vb + off + s];
        uint64_t kids = tree.prelim[lb + off + s] | tree.prelim[rb + off + s];
        uint64_t val = (fp & subset)
                     | ((pv | (fp & kids)) & intNonsub)
                     | ((pv | fp) & uniNonsub);
        mpr[vb + off + s] = val;
      }
    }
  }

  // Min-length-0 test per INTERNAL, non-root, non-root-child edge.  TNT
  // `collapse 3` collapses zero-length INTERNAL branches into polytomies;
  // terminal (pendant) edges are never collapsed, so we flag internal nodes
  // only (c >= n_tip).  This is both faithful to the operation and exactly
  // matches the brute-force oracle on internal edges (the terminal-edge marginal
  // test is confounded by the parsimony-uninformative characters the kernel
  // simplifies away — see b2_collapsed_kernel_validate.R).
  for (int c = tree.n_tip + 1; c < tree.n_node; ++c) {  // internal, skip root (= n_tip)
    int p = tree.parent[c];
    if (p == tree.n_tip) continue;          // root's children: skip (as conservative)
    size_t cb = static_cast<size_t>(c) * tw;
    size_t pb = static_cast<size_t>(p) * tw;
    bool collapsible = true;
    for (int b = 0; b < nb && collapsible; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;
      uint64_t acc = 0;
      for (int s = 0; s < k; ++s)
        acc |= mpr[pb + off + s] & mpr[cb + off + s];
      if ((acc & blk.active_mask) != blk.active_mask) collapsible = false;
    }
    if (collapsible) collapsed[c] = 1;
  }
}

void compute_collapsed_regions(
    const TreeState& tree,
    const DataSet& ds,
    CollapsedRegions& info) {

  // Step 1: compute per-node collapsed flags.
  compute_collapsed_flags(tree, ds, info.collapsed);

  const int n_node = tree.n_node;
  info.region_id.assign(n_node, -1);
  info.n_collapsed = 0;
  info.n_regions = 0;

  // Count collapsed edges.
  for (int c = 0; c < n_node; ++c) {
    if (info.collapsed[c]) ++info.n_collapsed;
  }
  if (info.n_collapsed == 0) return;

  // Step 2: assign region IDs.
  //
  // A collapsed edge connects child c to parent p (when collapsed[c] == 1).
  // Two nodes share a region if they are connected by a collapsed edge.
  //
  // Process internal nodes in REVERSE postorder (= preorder for trees):
  // parents visited before children. When a parent creates or joins a
  // region, its children inherit the same region_id.
  //
  // The root itself is never collapsed (no parent edge) and root's children
  // are excluded by compute_collapsed_flags(), so root always has
  // region_id == -1.
  const auto& po = tree.postorder;
  for (int idx = static_cast<int>(po.size()) - 1; idx >= 0; --idx) {
    int node = po[idx];
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    // Process left child
    if (info.collapsed[lc]) {
      if (info.region_id[node] >= 0) {
        // Parent already has a region — child joins it.
        info.region_id[lc] = info.region_id[node];
      } else {
        // Start a new region for both parent and child.
        int rid = info.n_regions++;
        info.region_id[node] = rid;
        info.region_id[lc] = rid;
      }
    }

    // Process right child
    if (info.collapsed[rc]) {
      if (info.region_id[node] >= 0) {
        info.region_id[rc] = info.region_id[node];
      } else {
        int rid = info.n_regions++;
        info.region_id[node] = rid;
        info.region_id[rc] = rid;
      }
    }
  }
}

} // namespace ts
