#include <Rcpp.h>
#include "ts_data.h"
#include "ts_tree.h"
#include "ts_fitch.h"
#include "ts_search.h"
#include "ts_tbr.h"

using namespace Rcpp;

namespace {

ts::DataSet make_dataset(
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
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

  return ts::build_dataset(
      REAL(contrast), n_tokens, n_states,
      INTEGER(tip_data), n_tips, n_patterns,
      INTEGER(weight),
      level_ptrs.data());
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
int ts_fitch_score(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

  ts::TreeState tree;
  tree.init_from_edge(
      &edge(0, 0), &edge(0, 1),
      edge.nrow(), ds);

  return ts::fitch_score(tree, ds);
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

  int whole_score = ts::fitch_score(tree, ds);
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

  int whole_score = ts::fitch_score(tree, ds);
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
        clip_score += ts::fitch_downpass_node(
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
  int actual = ts::fitch_score(tree, ds);

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
    int maxHits = 20)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

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
    int maxChanges = 0)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

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
List ts_nni_search(
    IntegerMatrix edge,
    NumericMatrix contrast,
    IntegerMatrix tip_data,
    IntegerVector weight,
    CharacterVector levels,
    int maxHits = 20)
{
  ts::DataSet ds = make_dataset(contrast, tip_data, weight, levels);

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
