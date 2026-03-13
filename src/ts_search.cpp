#include "ts_search.h"
#include "ts_fitch.h"
#include <algorithm>
#include <random>
#include <vector>

#include <Rcpp.h>

// R API for interrupt checking
#include <R.h>
#include <Rinternals.h>

namespace ts {

// ---- NNI search (Phase 1) ----

SearchResult nni_search(TreeState& tree, const DataSet& ds, int maxHits) {
  double best_score = fitch_score(tree, ds);
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  std::vector<int> edges = tree.nni_edges();
  int n_edges = static_cast<int>(edges.size());

  std::mt19937 rng(std::random_device{}());

  bool keep_going = true;
  while (keep_going) {
    keep_going = false;
    std::shuffle(edges.begin(), edges.end(), rng);

    for (int ei = 0; ei < n_edges; ++ei) {
      int c = edges[ei];

      for (int which = 0; which < 2; ++which) {
        auto undo = tree.nni_apply(c, which);
        tree.build_postorder();
        double new_score = fitch_score(tree, ds);
        ++n_iterations;

        if (new_score < best_score) {
          best_score = new_score;
          ++n_moves;
          hits = 1;
          keep_going = true;
          goto nni_next_pass;
        } else if (new_score == best_score) {
          ++hits;
          if (hits <= maxHits) {
            ++n_moves;
            keep_going = true;
            goto nni_next_pass;
          }
        }

        tree.nni_undo(undo);
        tree.build_postorder();
      }
    }

    nni_next_pass:
    R_CheckUserInterrupt();
  }

  return SearchResult{best_score, n_moves, n_iterations};
}

// ---- SPR search (Phase 2) ----

// Collect edges in the main (divided) tree reachable from root.
static void collect_destination_edges(
    const TreeState& tree,
    std::vector<std::pair<int,int>>& destinations)
{
  destinations.clear();

  std::vector<int> stack;
  stack.push_back(tree.n_tip);  // root

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    if (node < tree.n_tip) continue;

    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    destinations.push_back({node, lc});
    destinations.push_back({node, rc});

    stack.push_back(lc);
    stack.push_back(rc);
  }
}

// Score the tree from scratch: reset states, run full two-pass.
static int full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return fitch_score(tree, ds);
}

SearchResult spr_search(TreeState& tree, const DataSet& ds, int maxHits) {
  double best_score = full_rescore(tree, ds);
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  std::mt19937 rng(std::random_device{}());

  // Enumerate candidate clip nodes: any non-root node.
  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;  // root
    clip_candidates.push_back(node);
  }

  std::vector<std::pair<int,int>> destinations;
  bool keep_going = true;

  while (keep_going) {
    keep_going = false;
    std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

    for (int clip_node : clip_candidates) {
      // Skip children of root (see Lessons learned #2)
      if (tree.parent[clip_node] == tree.n_tip) continue;

      // --- Clip phase ---
      tree.spr_clip(clip_node);
      tree.build_postorder();

      // Score the divided tree (main + clipped parts)
      tree.reset_states(ds);
      int main_score = fitch_score(tree, ds);

      // Score the clipped subtree separately via postorder DFS
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
            const CharBlock& blk = ds.blocks[b];
            int offset = ds.block_word_offset[b];
            clip_score += fitch_downpass_node(
                &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset],
                &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset],
                &tree.prelim[static_cast<size_t>(nd) * tree.total_words + offset],
                blk.n_states, blk.active_mask);
          }
        }
      }

      double divided_length = main_score + clip_score;

      // Pointer to clipped subtree's basal prelim states
      const uint64_t* clip_prelim =
          &tree.prelim[static_cast<size_t>(clip_node) * tree.total_words];

      // --- Rearrangement phase: exact indirect calculation ---
      collect_destination_edges(tree, destinations);

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;

      bool accepted = false;
      for (auto& [above, below] : destinations) {
        if (above == nz && below == ns) continue;  // skip original position

        // Exact indirect calculation (Goloboff 1996): union-based virtual root
        int extra = fitch_indirect_length(clip_prelim, tree, ds, above, below);
        double candidate_score = divided_length + extra;
        ++n_iterations;

        if (candidate_score < best_score) {
          // Strict improvement: regraft and update states
          tree.spr_regraft(above, below);
          tree.build_postorder();
          best_score = full_rescore(tree, ds);
          ++n_moves;
          hits = 1;
          accepted = true;
          keep_going = true;
          break;
        } else if (candidate_score == best_score && hits <= maxHits) {
          // Plateau move
          tree.spr_regraft(above, below);
          tree.build_postorder();
          best_score = full_rescore(tree, ds);
          ++hits;
          ++n_moves;
          accepted = true;
          keep_going = true;
          break;
        }
      }

      if (!accepted) {
        tree.spr_unclip();
      }

      tree.build_postorder();

      if (keep_going) break;
    }

    R_CheckUserInterrupt();
  }

  // Final rescore for safety
  best_score = full_rescore(tree, ds);

  return SearchResult{best_score, n_moves, n_iterations};
}

} // namespace ts
