#include "ts_search.h"
#include "ts_fitch.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <Rcpp.h>

// R API for interrupt checking
#include <R.h>
#include <Rinternals.h>

namespace ts {

// ---- NNI search (Phase 1) ----

SearchResult nni_search(TreeState& tree, const DataSet& ds, int maxHits) {
  double best_score = score_tree(tree, ds);
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  std::vector<int> edges = tree.nni_edges();
  int n_edges = static_cast<int>(edges.size());

  // Seed from R's RNG for reproducibility with set.seed()
  GetRNGstate();
  std::mt19937 rng(static_cast<unsigned>(unif_rand() * 4294967295.0));
  PutRNGstate();

  bool keep_going = true;
  while (keep_going) {
    keep_going = false;
    std::shuffle(edges.begin(), edges.end(), rng);

    for (int ei = 0; ei < n_edges; ++ei) {
      int c = edges[ei];

      for (int which = 0; which < 2; ++which) {
        auto undo = tree.nni_apply(c, which);
        tree.build_postorder();
        double new_score = score_tree(tree, ds);
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
static double full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

SearchResult spr_search(TreeState& tree, const DataSet& ds, int maxHits) {
  double best_score = full_rescore(tree, ds);
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  // Seed from R's RNG for reproducibility with set.seed()
  GetRNGstate();
  std::mt19937 rng(static_cast<unsigned>(unif_rand() * 4294967295.0));
  PutRNGstate();

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

      tree.reset_states(ds);
      double main_score = score_tree(tree, ds);

      // Score clipped subtree separately
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
            clip_score += blk.weight * fitch_downpass_node(
                &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset],
                &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset],
                &tree.prelim[static_cast<size_t>(nd) * tree.total_words + offset],
                blk.n_states, blk.active_mask);
          }
        }
      }

      double divided_length = main_score + clip_score;

      const uint64_t* clip_prelim =
          &tree.prelim[static_cast<size_t>(clip_node) * tree.total_words];

      // IW: precompute base IW score and marginal deltas
      const bool use_iw = std::isfinite(ds.concavity);
      double base_iw = 0.0;
      std::vector<int> div_steps;
      std::vector<double> iw_del;
      if (use_iw) {
        div_steps.assign(ds.n_patterns, 0);
        iw_del.resize(ds.n_patterns, 0.0);
        // Main tree steps from local_cost
        for (int nd : tree.postorder) {
          for (int b = 0; b < ds.n_blocks; ++b) {
            uint64_t mask = tree.local_cost[
                static_cast<size_t>(nd) * tree.n_blocks + b];
            while (mask) {
              int c = ts::ctz64(mask);
              div_steps[ds.blocks[b].pattern_index[c]] += 1;
              mask &= mask - 1;
            }
          }
        }
        // Clip subtree steps (from prelim computed by fitch_downpass_node)
        {
          std::vector<int> cstack;
          cstack.push_back(clip_node);
          while (!cstack.empty()) {
            int nd = cstack.back();
            cstack.pop_back();
            if (nd < tree.n_tip) continue;
            int ni = nd - tree.n_tip;
            int lc = tree.left[ni];
            int rc = tree.right[ni];
            for (int b = 0; b < ds.n_blocks; ++b) {
              const CharBlock& blk = ds.blocks[b];
              int offset = ds.block_word_offset[b];
              const uint64_t* ls =
                  &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
              const uint64_t* rs =
                  &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
              uint64_t any_isect = 0;
              for (int s = 0; s < blk.n_states; ++s)
                any_isect |= (ls[s] & rs[s]);
              uint64_t nu = ~any_isect & blk.active_mask;
              while (nu) {
                int c = ts::ctz64(nu);
                div_steps[blk.pattern_index[c]] += 1;
                nu &= nu - 1;
              }
            }
            cstack.push_back(lc);
            cstack.push_back(rc);
          }
        }
        base_iw = compute_iw(ds, div_steps);
        precompute_iw_delta(ds, div_steps, iw_del);
      }

      // --- Rearrangement phase: screen with indirect calc, verify accepts ---
      collect_destination_edges(tree, destinations);

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;

      bool accepted = false;
      for (auto& [above, below] : destinations) {
        if (above == nz && below == ns) continue;

        double candidate_score;
        if (use_iw) {
          candidate_score = indirect_iw_length(clip_prelim, tree, ds,
                                               above, below, base_iw, iw_del);
        } else {
          int extra = fitch_indirect_length(clip_prelim, tree, ds, above, below);
          candidate_score = divided_length + extra;
        }
        ++n_iterations;

        bool dominated = (candidate_score > best_score) ||
                         (candidate_score == best_score && hits > maxHits);
        if (dominated) continue;

        // Candidate passes screen — regraft and verify with full rescore
        tree.spr_regraft(above, below);
        tree.build_postorder();
        double actual = full_rescore(tree, ds);

        if (actual < best_score) {
          best_score = actual;
          ++n_moves;
          hits = 1;
          accepted = true;
          keep_going = true;
          break;
        } else if (actual == best_score && hits <= maxHits) {
          ++hits;
          ++n_moves;
          accepted = true;
          keep_going = true;
          break;
        }

        // Not accepted — undo regraft
        tree.spr_unregraft(above, below);
      }

      if (!accepted) {
        tree.spr_unclip();
      }

      tree.build_postorder();

      if (keep_going) break;
    }

    R_CheckUserInterrupt();
  }

  best_score = full_rescore(tree, ds);

  return SearchResult{best_score, n_moves, n_iterations};
}

} // namespace ts
