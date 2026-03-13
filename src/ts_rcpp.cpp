#include <Rcpp.h>
#include "ts_data.h"
#include "ts_tree.h"
#include "ts_fitch.h"
#include "ts_search.h"
#include "ts_tbr.h"
#include "ts_drift.h"
#include "ts_ratchet.h"
#include "ts_splits.h"
#include "ts_pool.h"
#include "ts_wagner.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_driven.h"

using namespace Rcpp;

namespace {

ts::DataSet make_dataset(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  int n_tips = tip_data.nrow();
  int n_patterns = tip_data.ncol();
  int n_tokens = contrast.nrow();
  int n_states = contrast.ncol();

  std::vector<std::string> level_strs(n_states);
  std::vector<const char*> level_ptrs(n_states);
  for (int i = 0; i < n_states; ++i) {
    level_strs[i] = as<std::string>(levels[i]);
    level_ptrs[i] = level_strs[i].c_str();
  }

  const int* min_steps_ptr = (min_steps.size() > 0) ? INTEGER(min_steps)
                                                     : nullptr;

  return ts::build_dataset(
      REAL(contrast), n_tokens, n_states,
      INTEGER(tip_data), n_tips, n_patterns,
      INTEGER(weight),
      level_ptrs.data(),
      min_steps_ptr,
      concavity);
}

// Convert TreeState topology back to R edge matrix (2-column, 1-based)
IntegerMatrix tree_to_edge(const ts::TreeState& tree) {
  int n_edge = 2 * (tree.n_tip - 1);
  IntegerMatrix edge(n_edge, 2);
  int row = 0;
  for (int node = tree.n_tip; node < tree.n_node; ++node) {
    int ni = node - tree.n_tip;
    // +1 to convert to R's 1-based indexing
    edge(row, 0) = node + 1;
    edge(row, 1) = tree.left[ni] + 1;
    ++row;
    edge(row, 0) = node + 1;
    edge(row, 1) = tree.right[ni] + 1;
    ++row;
  }
  return edge;
}

} // anonymous namespace

// [[Rcpp::export]]
double ts_fitch_score(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  return ts::score_tree(tree, ds);
}

// [[Rcpp::export]]
List ts_na_debug_char(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int target_pattern)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);
  int total_score = ts::fitch_na_score(tree, ds);

  // Find which block/bit contains target_pattern (0-based)
  int tgt = target_pattern - 1;  // R 1-based to C 0-based
  int tgt_block = -1, tgt_bit = -1;
  for (int b = 0; b < ds.n_blocks; ++b) {
    for (int c = 0; c < ds.blocks[b].n_chars; ++c) {
      if (ds.blocks[b].pattern_index[c] == tgt) {
        tgt_block = b;
        tgt_bit = c;
      }
    }
  }
  if (tgt_block < 0) return List::create(Named("error") = "pattern not found");

  const ts::CharBlock& blk = ds.blocks[tgt_block];
  int off = ds.block_word_offset[tgt_block];
  int k = blk.n_states;
  uint64_t mask = 1ULL << tgt_bit;

  // Per-node info
  int n_node = tree.n_node;
  IntegerVector node_id(n_node);
  CharacterVector prelim_str(n_node), final_str(n_node), down2_str(n_node);
  IntegerVector is_step(n_node);

  auto state_str = [&](const uint64_t* base) -> std::string {
    std::string s;
    for (int st = 0; st < k; ++st) {
      if (base[st] & mask) {
        if (!s.empty()) s += "/";
        if (blk.has_inapplicable && st == 0) s += "-";
        else s += std::to_string(st - (blk.has_inapplicable ? 1 : 0));
      }
    }
    return s.empty() ? "." : s;
  };

  for (int nd = 0; nd < n_node; ++nd) {
    node_id[nd] = nd + 1;
    size_t base = static_cast<size_t>(nd) * tree.total_words + off;
    prelim_str[nd] = state_str(&tree.prelim[base]);
    final_str[nd] = state_str(&tree.final_[base]);
    down2_str[nd] = state_str(&tree.down2[base]);
    is_step[nd] = 0;
  }

  // Determine which nodes had steps in Pass 3
  // Must match fitch_na_score: l_act & r_act & ~(ss_app & any_d2_isect)
  for (int node : tree.postorder) {
    if (!blk.has_inapplicable) continue;
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words + off;
    size_t lb = static_cast<size_t>(lc) * tree.total_words + off;
    size_t rb = static_cast<size_t>(rc) * tree.total_words + off;
    const uint64_t* F = &tree.final_[nb];
    const uint64_t* L2 = &tree.down2[lb];
    const uint64_t* R2 = &tree.down2[rb];

    uint64_t ss_app = 0;
    for (int s = 1; s < k; ++s) ss_app |= F[s];

    uint64_t any_d2_isect = 0;
    for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

    const uint64_t* la = &tree.subtree_actives[lb];
    const uint64_t* ra = &tree.subtree_actives[rb];
    uint64_t l_act = 0, r_act = 0;
    for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

    uint64_t needs_step = l_act & r_act
                        & ~(ss_app & any_d2_isect) & blk.active_mask;
    if (needs_step & mask) is_step[node] = 1;
  }

  // Parent info
  IntegerVector parent_id(n_node);
  IntegerVector left_child(tree.n_internal), right_child(tree.n_internal);
  for (int nd = 0; nd < n_node; ++nd) parent_id[nd] = tree.parent[nd] + 1;
  for (int ni = 0; ni < tree.n_internal; ++ni) {
    left_child[ni] = tree.left[ni] + 1;
    right_child[ni] = tree.right[ni] + 1;
  }

  return List::create(
    Named("total_score") = total_score,
    Named("block") = tgt_block + 1,
    Named("bit") = tgt_bit,
    Named("has_inapp") = blk.has_inapplicable,
    Named("n_states") = k,
    Named("node") = node_id,
    Named("parent") = parent_id,
    Named("left") = left_child,
    Named("right") = right_child,
    Named("prelim") = prelim_str,
    Named("final_state") = final_str,
    Named("down2") = down2_str,
    Named("is_step") = is_step,
    Named("n_tip") = tree.n_tip
  );
}

// [[Rcpp::export]]
List ts_na_char_steps(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  // Run the three-pass scoring
  int total_score = ts::fitch_na_score(tree, ds);

  int n_pat = tip_data.ncol();
  IntegerVector steps(n_pat, 0);
  int debug_std_total = 0, debug_na_total = 0;

  // Standard blocks: count from local_cost
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const ts::CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable) continue;
      uint64_t mask = tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      debug_std_total += blk.weight * ts::popcount64(mask);
      for (int c = 0; c < blk.n_chars; ++c) {
        if (mask & (1ULL << c)) {
          steps[blk.pattern_index[c]] += 1;
        }
      }
    }
  }

  // NA blocks: use same formula as Pass 3 (needs subtree_actives)
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      if (!ds.blocks[b].has_inapplicable) continue;
      const ts::CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      const uint64_t* D2 = &tree.down2[nb + off];

      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= D2[s];

      uint64_t any_d2_isect = 0;
      for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_d2_isect) & blk.active_mask;
      debug_na_total += blk.weight * ts::popcount64(needs_step);

      for (int c = 0; c < blk.n_chars; ++c) {
        if (needs_step & (1ULL << c)) {
          steps[blk.pattern_index[c]] += 1;
        }
      }
    }
  }

  // Block info
  int n_std = 0, n_na = 0;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) ++n_na; else ++n_std;
  }

  return List::create(
    Named("steps") = steps,
    Named("total_score") = total_score,
    Named("std_total") = debug_std_total,
    Named("na_total") = debug_na_total,
    Named("debug_sum") = debug_std_total + debug_na_total,
    Named("n_std_blocks") = n_std,
    Named("n_na_blocks") = n_na
  );
}

// [[Rcpp::export]]
List ts_debug_clip(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int clip_node_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  int whole_score = ts::fitch_na_score(tree, ds);
  int clip_node = clip_node_1based - 1;  // 0-based

  tree.spr_clip(clip_node);
  tree.build_postorder();

  int divided_score = ts::fitch_score(tree, ds);

  // Also check postorder size
  int postorder_size = static_cast<int>(tree.postorder.size());

  // Unclip
  tree.spr_unclip();
  tree.build_postorder();
  int restored_score = ts::fitch_score(tree, ds);

  return List::create(
    Named("whole_score") = whole_score,
    Named("divided_score") = divided_score,
    Named("restored_score") = restored_score,
    Named("postorder_size") = postorder_size,
    Named("n_internal") = tree.n_internal,
    Named("clip_parent_0based") = tree.clip_state.clip_parent,
    Named("clip_sibling_0based") = tree.clip_state.clip_sibling,
    Named("clip_grandpar_0based") = tree.clip_state.clip_grandpar
  );
}

// [[Rcpp::export]]
List ts_test_indirect(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int clip_node_1based,
    int above_1based,
    int below_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);
  ts::TreeState tree;
  tree.init_from_edge(&edge(0, 0), &edge(0, 1), edge.nrow(), ds);

  int whole_score = ts::fitch_na_score(tree, ds);
  int clip_node = clip_node_1based - 1;
  int above = above_1based - 1;
  int below = below_1based - 1;

  // Clip
  tree.spr_clip(clip_node);
  tree.build_postorder();
  tree.reset_states(ds);
  int main_score = ts::fitch_score(tree, ds);

  // Score clipped subtree
  int clip_score = 0;
  {
    std::vector<int> clip_stack;
    clip_stack.push_back(clip_node);
    std::vector<int> clip_preorder;
    while (!clip_stack.empty()) {
      int nd = clip_stack.back();
      clip_stack.pop_back();
      if (nd < tree.n_tip) continue;
      clip_preorder.push_back(nd);
      int ni = nd - tree.n_tip;
      clip_stack.push_back(tree.right[ni]);
      clip_stack.push_back(tree.left[ni]);
    }
    for (int j = static_cast<int>(clip_preorder.size()) - 1; j >= 0; --j) {
      int nd = clip_preorder[j];
      int ni = nd - tree.n_tip;
      int lc = tree.left[ni];
      int rc = tree.right[ni];
      for (int b = 0; b < ds.n_blocks; ++b) {
        const ts::CharBlock& blk = ds.blocks[b];
        int offset = ds.block_word_offset[b];
        clip_score += blk.weight * ts::fitch_downpass_node(
            &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset],
            &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset],
            &tree.prelim[static_cast<size_t>(nd) * tree.total_words + offset],
            blk.n_states, blk.active_mask);
      }
    }
  }

  int divided_length = main_score + clip_score;
  const uint64_t* clip_prelim =
      &tree.prelim[static_cast<size_t>(clip_node) * tree.total_words];
  int extra = ts::fitch_indirect_length(clip_prelim, tree, ds, above, below);
  int candidate = divided_length + extra;

  // Now regraft and get actual score
  tree.spr_regraft(above, below);
  tree.build_postorder();
  tree.reset_states(ds);
  int actual = ts::fitch_na_score(tree, ds);

  return List::create(
    Named("whole_score") = whole_score,
    Named("main_score") = main_score,
    Named("clip_score") = clip_score,
    Named("divided_length") = divided_length,
    Named("extra") = extra,
    Named("candidate") = candidate,
    Named("actual") = actual,
    Named("match") = (candidate == actual)
  );
}

// [[Rcpp::export]]
List ts_spr_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 20,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SearchResult result = ts::spr_search(tree, ds, maxHits);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score,
    Named("n_moves") = result.n_moves,
    Named("n_iterations") = result.n_iterations
  );
}

// [[Rcpp::export]]
List ts_tbr_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 1,
    bool acceptEqual = false,
    int maxChanges = 0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::TBRParams params;
  params.max_hits = maxHits;
  params.accept_equal = acceptEqual;
  params.max_accepted_changes = maxChanges;

  ts::TBRResult result = ts::tbr_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_accepted") = result.n_accepted,
    Named("n_evaluated") = result.n_evaluated,
    Named("converged") = result.converged
  );
}

// [[Rcpp::export]]
List ts_ratchet_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nCycles = 10,
    double perturbProb = 0.04,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::RatchetParams params;
  params.n_cycles = nCycles;
  params.perturb_prob = perturbProb;
  params.max_hits = maxHits;

  ts::RatchetResult result = ts::ratchet_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_cycles") = result.n_cycles_completed,
    Named("total_tbr_moves") = result.total_tbr_moves
  );
}

// [[Rcpp::export]]
List ts_drift_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nCycles = 10,
    int afdLimit = 3,
    double rfdLimit = 0.1,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::DriftParams params;
  params.n_cycles = nCycles;
  params.afd_limit = afdLimit;
  params.rfd_limit = rfdLimit;
  params.max_hits = maxHits;

  ts::DriftResult result = ts::drift_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_cycles_completed") = result.n_cycles_completed,
    Named("total_drift_moves") = result.total_drift_moves,
    Named("total_tbr_moves") = result.total_tbr_moves
  );
}

// [[Rcpp::export]]
List ts_wagner_tree(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector addition_order = IntegerVector(),
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  std::vector<int> order;
  if (addition_order.size() > 0) {
    // Convert from 1-based R indices to 0-based
    order.resize(addition_order.size());
    for (int i = 0; i < addition_order.size(); ++i) {
      order[i] = addition_order[i] - 1;
    }
  }

  ts::TreeState tree;
  ts::WagnerResult result = ts::wagner_tree(tree, ds, order);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score
  );
}

// [[Rcpp::export]]
List ts_random_wagner_tree(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  ts::WagnerResult result = ts::random_wagner_tree(tree, ds);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score
  );
}

// [[Rcpp::export]]
List ts_compute_splits(
    IntegerMatrix edge,
    int n_tip)
{
  // Build a minimal TreeState from the edge matrix (no character data needed).
  ts::TreeState tree;
  // We need a dummy DataSet just for init_from_edge dimensions.
  // Instead, manually initialize the tree topology.
  tree.n_tip = n_tip;
  tree.n_internal = n_tip - 1;
  tree.n_node = 2 * n_tip - 1;
  tree.total_words = 0;
  tree.n_blocks = 0;

  tree.parent.assign(tree.n_node, 0);
  tree.left.assign(tree.n_internal, 0);
  tree.right.assign(tree.n_internal, 0);

  // Parse edge matrix (R is 1-based)
  int n_edge = edge.nrow();
  for (int i = 0; i < n_edge; ++i) {
    int par = edge(i, 0) - 1;
    int child = edge(i, 1) - 1;
    tree.parent[child] = par;
    int ni = par - n_tip;
    if (tree.left[ni] == 0 && par != child) {
      // Check if left is unset (0 could be tip 0, but root's self-reference
      // won't use this path). We need a cleaner approach.
      // Actually: first child encountered goes left, second goes right.
      // Use a flag: if left == 0 and we haven't set it yet...
      // Simpler: track how many children we've seen per node.
      tree.left[ni] = child;
    } else {
      tree.right[ni] = child;
    }
  }
  // Fix: root's parent is itself
  tree.parent[n_tip] = n_tip;

  // Fix left/right assignment: re-parse properly
  // Reset and use a counter
  std::vector<int> child_count(tree.n_internal, 0);
  tree.left.assign(tree.n_internal, -1);
  tree.right.assign(tree.n_internal, -1);
  for (int i = 0; i < n_edge; ++i) {
    int par = edge(i, 0) - 1;
    int child = edge(i, 1) - 1;
    tree.parent[child] = par;
    int ni = par - n_tip;
    if (child_count[ni] == 0) {
      tree.left[ni] = child;
    } else {
      tree.right[ni] = child;
    }
    ++child_count[ni];
  }
  tree.parent[n_tip] = n_tip;

  tree.build_postorder();

  ts::SplitSet ss = ts::compute_splits(tree);

  // Return as a list of raw vectors (one per split)
  List result(ss.n_splits);
  for (int i = 0; i < ss.n_splits; ++i) {
    const uint64_t* s = ss.split(i);
    // Convert to an integer vector of tip indices in the split
    IntegerVector tips;
    for (int t = 0; t < n_tip; ++t) {
      int w = t / 64;
      int b = t % 64;
      if (s[w] & (1ULL << b)) {
        tips.push_back(t + 1);  // 1-based for R
      }
    }
    result[i] = tips;
  }

  return result;
}

// [[Rcpp::export]]
bool ts_trees_equal(
    IntegerMatrix edge1,
    IntegerMatrix edge2,
    int n_tip)
{
  auto build_tree = [&](IntegerMatrix edge) -> ts::TreeState {
    ts::TreeState tree;
    tree.n_tip = n_tip;
    tree.n_internal = n_tip - 1;
    tree.n_node = 2 * n_tip - 1;
    tree.total_words = 0;
    tree.n_blocks = 0;
    tree.parent.assign(tree.n_node, 0);
    tree.left.assign(tree.n_internal, -1);
    tree.right.assign(tree.n_internal, -1);

    int n_edge = edge.nrow();
    std::vector<int> child_count(tree.n_internal, 0);
    for (int i = 0; i < n_edge; ++i) {
      int par = edge(i, 0) - 1;
      int child = edge(i, 1) - 1;
      tree.parent[child] = par;
      int ni = par - n_tip;
      if (child_count[ni] == 0) {
        tree.left[ni] = child;
      } else {
        tree.right[ni] = child;
      }
      ++child_count[ni];
    }
    tree.parent[n_tip] = n_tip;
    tree.build_postorder();
    return tree;
  };

  ts::TreeState t1 = build_tree(edge1);
  ts::TreeState t2 = build_tree(edge2);

  ts::SplitSet s1 = ts::compute_splits(t1);
  ts::SplitSet s2 = ts::compute_splits(t2);

  return ts::splits_equal(s1, s2);
}

// [[Rcpp::export]]
List ts_pool_test(
    List edges,        // list of edge matrices
    NumericVector scores,
    int n_tip,
    int max_size = 100,
    double suboptimal = 0.0)
{
  auto build_tree = [&](IntegerMatrix edge) -> ts::TreeState {
    ts::TreeState tree;
    tree.n_tip = n_tip;
    tree.n_internal = n_tip - 1;
    tree.n_node = 2 * n_tip - 1;
    tree.total_words = 0;
    tree.n_blocks = 0;
    tree.parent.assign(tree.n_node, 0);
    tree.left.assign(tree.n_internal, -1);
    tree.right.assign(tree.n_internal, -1);

    int n_edge = edge.nrow();
    std::vector<int> child_count(tree.n_internal, 0);
    for (int i = 0; i < n_edge; ++i) {
      int par = edge(i, 0) - 1;
      int child = edge(i, 1) - 1;
      tree.parent[child] = par;
      int ni = par - n_tip;
      if (child_count[ni] == 0) {
        tree.left[ni] = child;
      } else {
        tree.right[ni] = child;
      }
      ++child_count[ni];
    }
    tree.parent[n_tip] = n_tip;
    tree.build_postorder();
    return tree;
  };

  ts::TreePool pool(max_size, suboptimal);

  int n_trees = edges.size();
  LogicalVector added(n_trees);
  for (int i = 0; i < n_trees; ++i) {
    IntegerMatrix e = as<IntegerMatrix>(edges[i]);
    ts::TreeState t = build_tree(e);
    added[i] = pool.add(t, scores[i]);
  }

  return List::create(
    Named("added") = added,
    Named("pool_size") = pool.size(),
    Named("best_score") = pool.best_score(),
    Named("hits_to_best") = pool.hits_to_best()
  );
}

// [[Rcpp::export]]
List ts_nni_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 20,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SearchResult result = ts::nni_search(tree, ds, maxHits);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.score,
    Named("n_moves") = result.n_moves,
    Named("n_iterations") = result.n_iterations
  );
}

// [[Rcpp::export]]
List ts_tree_fuse(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    List pool_edges,
    NumericVector pool_scores,
    bool accept_equal = false,
    int max_rounds = 10,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  // Build recipient tree
  ts::TreeState recipient;
  recipient.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  // Build pool from R-side inputs
  ts::TreePool pool(static_cast<int>(pool_edges.size()), 1e18);
  for (int i = 0; i < pool_edges.size(); ++i) {
    IntegerMatrix pe = as<IntegerMatrix>(pool_edges[i]);
    ts::TreeState t;
    t.init_from_edge(&pe(0, 0), &pe(0, 1), pe.nrow(), ds);
    double sc = pool_scores[i];
    pool.add(t, sc);
  }

  ts::FuseParams params;
  params.accept_equal = accept_equal;
  params.max_rounds = max_rounds;

  ts::FuseResult result = ts::tree_fuse(recipient, ds, pool, params);

  return List::create(
    Named("edge") = tree_to_edge(recipient),
    Named("score") = result.best_score,
    Named("n_exchanges") = result.n_exchanges,
    Named("n_rounds") = result.n_rounds
  );
}

// [[Rcpp::export]]
List ts_sector_diag(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int sector_root_1based)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  int full_score = ts::fitch_na_score(tree, ds);

  int sector_root = sector_root_1based - 1;
  int clade_size = ts::count_clade_tips(tree, sector_root);

  ts::ReducedDataset rd = ts::build_reduced_dataset(tree, ds, sector_root);

  int sector_score = ts::fitch_score(rd.subtree, rd.data);

  return List::create(
    Named("full_score") = full_score,
    Named("sector_root") = sector_root,
    Named("clade_size") = clade_size,
    Named("n_sector_tips") = rd.subtree.n_tip,
    Named("n_sector_nodes") = rd.subtree.n_node,
    Named("sector_score") = sector_score,
    Named("sector_total_words") = rd.data.total_words,
    Named("sector_n_blocks") = rd.data.n_blocks
  );
}

// [[Rcpp::export]]
List ts_rss_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int minSectorSize = 6,
    int maxSectorSize = 50,
    bool acceptEqual = false,
    int rssPicks = 0,
    int ratchetCycles = 6,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SectorParams params;
  params.min_sector_size = minSectorSize;
  params.max_sector_size = maxSectorSize;
  params.accept_equal = acceptEqual;
  params.rss_picks_per_round = rssPicks;
  params.internal_ratchet_cycles = ratchetCycles;
  params.internal_max_hits = maxHits;

  ts::SectorResult result = ts::rss_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_sectors_searched") = result.n_sectors_searched,
    Named("n_sectors_improved") = result.n_sectors_improved,
    Named("total_steps_saved") = result.total_steps_saved
  );
}

// [[Rcpp::export]]
List ts_xss_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int nPartitions = 4,
    int xssRounds = 3,
    bool acceptEqual = false,
    int ratchetCycles = 6,
    int maxHits = 1,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  ts::SectorParams params;
  params.n_partitions = nPartitions;
  params.xss_rounds = xssRounds;
  params.accept_equal = acceptEqual;
  params.internal_ratchet_cycles = ratchetCycles;
  params.internal_max_hits = maxHits;

  ts::SectorResult result = ts::xss_search(tree, ds, params);

  return List::create(
    Named("edge") = tree_to_edge(tree),
    Named("score") = result.best_score,
    Named("n_sectors_searched") = result.n_sectors_searched,
    Named("n_sectors_improved") = result.n_sectors_improved,
    Named("total_steps_saved") = result.total_steps_saved
  );
}

// [[Rcpp::export]]
List ts_driven_search(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxReplicates = 100,
    int targetHits = 10,
    int tbrMaxHits = 1,
    int ratchetCycles = 10,
    double ratchetPerturbProb = 0.04,
    int driftCycles = 6,
    int driftAfdLimit = 3,
    double driftRfdLimit = 0.1,
    int xssRounds = 3,
    int xssPartitions = 4,
    int sectorMinSize = 6,
    int sectorMaxSize = 50,
    int fuseInterval = 3,
    bool fuseAcceptEqual = false,
    int poolMaxSize = 100,
    double poolSuboptimal = 0.0,
    IntegerVector min_steps = IntegerVector(),
    double concavity = R_PosInf)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels,
                                min_steps, concavity);

  ts::DrivenParams params;
  params.max_replicates = maxReplicates;
  params.target_hits = targetHits;
  params.tbr_max_hits = tbrMaxHits;
  params.ratchet_cycles = ratchetCycles;
  params.ratchet_perturb_prob = ratchetPerturbProb;
  params.drift_cycles = driftCycles;
  params.drift_afd_limit = driftAfdLimit;
  params.drift_rfd_limit = driftRfdLimit;
  params.xss_rounds = xssRounds;
  params.xss_partitions = xssPartitions;
  params.sector_min_size = sectorMinSize;
  params.sector_max_size = sectorMaxSize;
  params.fuse_interval = fuseInterval;
  params.fuse_accept_equal = fuseAcceptEqual;
  params.pool_max_size = poolMaxSize;
  params.pool_suboptimal = poolSuboptimal;

  ts::TreeState best_tree;
  ts::DrivenResult result = ts::driven_search(best_tree, ds, params);

  if (result.pool_size == 0) {
    // No replicates completed (e.g. max_replicates=0) — no valid tree
    return List::create(
      Named("edge") = IntegerMatrix(0, 2),
      Named("score") = result.best_score,
      Named("replicates") = result.replicates_completed,
      Named("hits_to_best") = result.hits_to_best,
      Named("pool_size") = result.pool_size
    );
  }

  return List::create(
    Named("edge") = tree_to_edge(best_tree),
    Named("score") = result.best_score,
    Named("replicates") = result.replicates_completed,
    Named("hits_to_best") = result.hits_to_best,
    Named("pool_size") = result.pool_size
  );
}
