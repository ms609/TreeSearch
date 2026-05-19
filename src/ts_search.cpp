#include "ts_search.h"
#include "ts_fitch.h"
#include "ts_collapsed.h"
#include "ts_rng.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <random>
#include <vector>

#include <Rcpp.h>

#include <R.h>
#include <Rinternals.h>

namespace ts {

// ---- Helpers (file-local) ----

static double full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}

// Compute the number of tips in the subtree below each node.
static void compute_subtree_sizes(const TreeState& tree,
                                  std::vector<int>& sizes) {
  sizes.assign(tree.n_node, 0);
  for (int i = 0; i < tree.n_tip; ++i) sizes[i] = 1;
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    sizes[node] = sizes[tree.left[ni]] + sizes[tree.right[ni]];
  }
}

// ---- NNI search ----

SearchResult nni_search(TreeState& tree, const DataSet& ds, int maxHits,
                        std::function<bool()> check_timeout) {
  double best_score = score_tree(tree, ds);
  // No informative characters: all trees have the same score.
  if (ds.total_words == 0) return {best_score, 0, 0};
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  std::vector<int> edges = tree.nni_edges();
  int n_edges = static_cast<int>(edges.size());

  // Detect inapplicable characters
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  // Poll timeout every n_tip edges (capped at 50)
  const int timeout_interval = std::max(1, std::min(tree.n_tip, 50));

  bool keep_going = true;
  while (keep_going) {
    keep_going = false;
    std::shuffle(edges.begin(), edges.end(), rng);

    for (int ei = 0; ei < n_edges; ++ei) {
      int c = edges[ei];

      for (int which = 0; which < 2; ++which) {
        auto undo = tree.nni_apply(c, which);
        ++n_iterations;

        double new_score;
        if (has_na) {
          // NA datasets: fall back to full rescore (three-pass needed)
          tree.build_postorder();
          new_score = score_tree(tree, ds);
        } else {
          // Incremental downpass: O(depth × C) instead of O(n × C)
          tree.clip_undo_stack.clear();
          int delta = fitch_incremental_downpass(tree, ds, c);
          if (std::isfinite(ds.concavity)) {
            // Weighted (IW or profile): integer EW delta cannot be added
            // to a float weighted score.  After the chain walk,
            // local_cost is correct for the whole tree (NNI only changes
            // children at edge c; off-chain nodes retain valid local_cost
            // from the score_tree at function entry), so extract per-
            // pattern step counts and dispatch by ds.scoring_mode.
            std::vector<int> char_steps(ds.n_patterns, 0);
            extract_char_steps(tree, ds, char_steps);
            new_score = compute_weighted_score(ds, char_steps);
          } else {
            new_score = best_score + delta;
          }
        }

        if (new_score < best_score) {
          best_score = new_score;
          // Update postorder + final_ arrays for subsequent iterations
          tree.build_postorder();
          tree.clip_undo_stack.clear();
          fitch_uppass(tree, ds);
          ++n_moves;
          hits = 1;
          keep_going = true;
          goto nni_next_pass;
        } else if (new_score == best_score) {
          ++hits;
          if (hits <= maxHits) {
            tree.build_postorder();
            tree.clip_undo_stack.clear();
            fitch_uppass(tree, ds);
            ++n_moves;
            keep_going = true;
            goto nni_next_pass;
          }
        }

        tree.nni_undo(undo);
        if (!has_na) {
          // Restore prelim/local_cost by re-scoring the original topology
          tree.clip_undo_stack.clear();
          fitch_incremental_downpass(tree, ds, c);
          tree.clip_undo_stack.clear();
        } else {
          tree.build_postorder();
        }
      }

      // Periodic timeout check
      if (check_timeout && (ei % timeout_interval == 0) && check_timeout()) {
        keep_going = false;
        goto nni_done;
      }
    }

    nni_next_pass:
    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  nni_done:
  // Authoritative final score
  tree.build_postorder();
  best_score = full_rescore(tree, ds);

  return SearchResult{best_score, n_moves, n_iterations};
}

// ---- SPR search ----

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

SearchResult spr_search(TreeState& tree, const DataSet& ds, int maxHits,
                        std::function<bool()> check_timeout) {
  double best_score = full_rescore(tree, ds);
  // No informative characters: all trees have the same score.
  if (ds.total_words == 0) return {best_score, 0, 0};
  int n_moves = 0;
  int n_iterations = 0;
  int hits = 1;

  const bool use_iw = std::isfinite(ds.concavity);
  const double eps = use_iw ? 1e-10 : 0.0;

  // Detect inapplicable characters
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;  // root
    clip_candidates.push_back(node);
  }

  // Collapsed flags: edges that provably cannot yield an improvement
  // (clip skipping + regraft merging).
  std::vector<uint8_t> collapsed;
  compute_collapsed_flags(tree, ds, collapsed);

  std::vector<std::pair<int,int>> destinations;

  // Pre-allocate IW buffers
  std::vector<int> div_steps;
  std::vector<double> iw_del;
  if (use_iw) {
    div_steps.resize(ds.n_patterns, 0);
    iw_del.resize(ds.n_patterns, 0.0);
  }

  // Subtree sizes for smaller-subtree filtering
  std::vector<int> subtree_sizes;

  bool keep_going = true;
  bool need_shuffle = true;
  const int timeout_interval = std::max(tree.n_tip, 50);
  int clips_since_timeout_check = 0;
  bool timed_out = false;

  while (keep_going && !timed_out) {
    keep_going = false;

    // Deferred reshuffling: only reshuffle when previous pass found nothing
    if (need_shuffle) {
      std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);
    }
    need_shuffle = true;

    // Recompute subtree sizes for smaller-subtree filtering
    compute_subtree_sizes(tree, subtree_sizes);

    for (int clip_node : clip_candidates) {
      if (tree.parent[clip_node] == tree.n_tip) continue;

      // Skip collapsed edges (zero-length, provably unimprovable).
      if (!collapsed.empty() && collapsed[clip_node])
        continue;

      // Smaller-subtree filtering: skip clips of the larger half
      if (subtree_sizes[clip_node] > tree.n_tip / 2) continue;

      // Save clip subtree's actives before clipping (for NA indirect)
      const uint64_t* clip_actives = nullptr;
      std::vector<uint64_t> clip_actives_buf;
      if (has_na) {
        size_t clip_sa_base =
            static_cast<size_t>(clip_node) * tree.total_words;
        clip_actives_buf.assign(
            tree.subtree_actives.begin() + clip_sa_base,
            tree.subtree_actives.begin() + clip_sa_base + tree.total_words);
        clip_actives = clip_actives_buf.data();
      }

      // --- Clip phase: incremental scoring (matches TBR pattern) ---
      tree.spr_clip(clip_node);
      tree.build_postorder();

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;
      int nx = tree.clip_state.clip_parent;

      double divided_length;
      if (has_na) {
        fitch_na_incremental_downpass(tree, ds, nz);
        fitch_na_incremental_uppass(tree, ds, nz);
        divided_length = static_cast<double>(fitch_na_pass3_score(tree, ds));
      } else {
        int delta = fitch_incremental_downpass(tree, ds, nz);
        fitch_incremental_uppass(tree, ds, nz);

        int nx_cost = 0;
        for (int b = 0; b < ds.n_blocks; ++b) {
          uint64_t lc = tree.local_cost[static_cast<size_t>(nx) * tree.n_blocks + b];
          int nu = popcount64(lc);
          if (ds.blocks[b].upweight_mask) nu += popcount64(lc & ds.blocks[b].upweight_mask);
          nx_cost += ds.blocks[b].weight * nu;
        }
        divided_length = best_score + delta - nx_cost;
      }

      const uint64_t* clip_prelim =
          &tree.prelim[static_cast<size_t>(clip_node) * tree.total_words];

      // IW: precompute base score and marginal deltas
      double base_iw = 0.0;
      if (use_iw) {
        std::fill(div_steps.begin(), div_steps.end(), 0);
        extract_char_steps(tree, ds, div_steps);
        base_iw = compute_weighted_score(ds, div_steps);
        precompute_weighted_delta(ds, div_steps, iw_del);
      }

      // --- Rearrangement phase: screen with bounded indirect calc ---
      collect_destination_edges(tree, destinations);
      // Partial shuffle: seed bound with diverse sample
      {
        int ne = static_cast<int>(destinations.size());
        int k = std::min(20, ne);
        for (int i = 0; i < k; ++i) {
          std::uniform_int_distribution<int> dist(i, ne - 1);
          std::swap(destinations[i], destinations[dist(rng)]);
        }
      }

      double best_candidate = HUGE_VAL;
      int best_above = -1, best_below = -1;

      for (auto& [above, below] : destinations) {
        if (above == nz && below == ns) continue;

        // Collapsed-region regraft merging: skip interior collapsed edges.
        if (!collapsed.empty() && collapsed[below])
          continue;

        double candidate_score;
        if (has_na) {
          if (use_iw) {
            candidate_score = indirect_na_iw_length_bounded(
                clip_prelim, clip_actives, tree, ds, above, below,
                base_iw, iw_del, best_candidate);
          } else {
            int cutoff = (best_candidate < HUGE_VAL)
                ? static_cast<int>(best_candidate - divided_length + 1)
                : INT_MAX;
            int extra = fitch_na_indirect_length_bounded(
                clip_prelim, clip_actives, tree, ds, above, below, cutoff);
            candidate_score = divided_length + extra;
          }
        } else if (use_iw) {
          candidate_score = indirect_iw_length_bounded(
              clip_prelim, tree, ds, above, below, base_iw, iw_del,
              best_candidate);
        } else {
          int cutoff = (best_candidate < HUGE_VAL)
              ? static_cast<int>(best_candidate - divided_length + 1)
              : INT_MAX;
          int extra = fitch_indirect_length_bounded(
              clip_prelim, tree, ds, above, below, cutoff);
          candidate_score = divided_length + extra;
        }
        ++n_iterations;

        if (candidate_score < best_candidate) {
          best_candidate = candidate_score;
          best_above = above;
          best_below = below;
        }
      }

      // --- Verify best candidate with full rescore ---
      bool dominated = (best_candidate > best_score + eps) ||
                       (best_candidate > best_score - eps
                        && hits > maxHits);

      bool accepted = false;

      if (!dominated && best_above >= 0) {
        tree.spr_regraft(best_above, best_below);
        tree.build_postorder();
        double actual = full_rescore(tree, ds);

        if (actual < best_score - eps) {
          best_score = actual;
          ++n_moves;
          hits = 1;
          accepted = true;
          keep_going = true;
        } else if (std::fabs(actual - best_score) <= eps
                   && hits <= maxHits) {
          ++hits;
          ++n_moves;
          accepted = true;
          keep_going = true;
        }

        if (!accepted) {
          tree.spr_unregraft(best_above, best_below);
        }
      }

      bool regraft_was_rescored = !accepted && !dominated && best_above >= 0;

      if (!accepted) {
        tree.spr_unclip();
      }

      tree.build_postorder();

      // full_rescore during rejected regraft overwrites all state arrays.
      // spr_unclip only restores states along the clip-to-root path.
      // Recompute so indirect scoring is correct for subsequent clips.
      if (regraft_was_rescored) {
        full_rescore(tree, ds);
      }

      if (keep_going) {
        // Recompute collapsed flags after the accepted move.
        compute_collapsed_flags(tree, ds, collapsed);
        // Deferred reshuffling: don't reshuffle after acceptance
        need_shuffle = false;
        break;
      }

      if (ts::check_interrupt()) { keep_going = false; break; }
      ++clips_since_timeout_check;
      if (check_timeout && clips_since_timeout_check >= timeout_interval) {
        clips_since_timeout_check = 0;
        if (check_timeout()) { timed_out = true; break; }
      }
    }

    if (ts::check_interrupt()) break;
  }

  best_score = full_rescore(tree, ds);

  return SearchResult{best_score, n_moves, n_iterations};
}

} // namespace ts
