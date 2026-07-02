#include "ts_nni_perturb.h"
#include "ts_constraint.h"
#include "ts_tbr.h"
#include "ts_fitch.h"
#include "ts_rng.h"

#include <algorithm>
#include <random>
#include <unordered_set>
#include <R.h>

namespace ts {

namespace {

void copy_topology(TreeState& dst, const TreeState& src) {
  dst.parent = src.parent;
  dst.left   = src.left;
  dst.right  = src.right;
}

} // anonymous namespace


int random_nni_perturb(TreeState& tree, double fraction) {
  std::vector<int> edges = tree.nni_edges();
  std::mt19937 rng = ts::make_rng();
  std::shuffle(edges.begin(), edges.end(), rng);

  std::bernoulli_distribution coin(fraction);
  std::uniform_int_distribution<int> which_dist(0, 1);

  // Track nodes involved in applied swaps. An edge `e` conflicts with
  // a prior swap if `e` or `parent[e]` is in the touched set.
  // This ensures no two adjacent edges are both swapped.
  std::unordered_set<int> touched;
  int n_applied = 0;

  for (int c : edges) {
    if (touched.count(c) || touched.count(tree.parent[c])) continue;
    if (!coin(rng)) continue;

    tree.nni_apply(c, which_dist(rng));
    touched.insert(c);
    touched.insert(tree.parent[c]);
    ++n_applied;
  }

  if (n_applied > 0) {
    tree.build_postorder();
  }

  return n_applied;
}


NNIPerturbResult nni_perturb_search(
    TreeState& tree, const DataSet& ds,
    const NNIPerturbParams& params,
    ConstraintData* cd,
    std::function<bool()> check_timeout)
{
  // Initial TBR to establish baseline
  TBRParams search_params;
  search_params.accept_equal = false;
  search_params.max_accepted_changes = 0;
  search_params.max_hits = params.max_hits;
  search_params.tabu_size = params.tabu_size;

  TBRResult initial = tbr_search(tree, ds, search_params, cd,
                                  nullptr, nullptr, check_timeout);

  double best_score = initial.best_score;
  int total_moves = initial.n_accepted;
  int n_escapes = 0;

  TreeState best_tree = tree;

  int cycles_completed = 0;

  for (int cycle = 0; cycle < params.n_cycles; ++cycle) {
    // 1. Perturb topology: random compatible NNI swaps
    int n_swaps = random_nni_perturb(tree, params.perturb_fraction);

    if (n_swaps == 0) {
      ++cycles_completed;
      continue;
    }

    // Repair constraint violations from blind NNI perturbation, then
    // re-sync constraint metadata for the (now repaired) topology.
    // update_constraint must be called even when impose_constraint is
    // skipped — cd->constraint_node and DFS timestamps are stale after
    // the topology change from random_nni_perturb.
    if (cd) {
      if (cd->active) impose_constraint(tree, *cd);
      update_constraint(tree, *cd);
    }

    // Rescore after perturbation (+ repair)
    tree.reset_states(ds);
    score_tree(tree, ds);

    // 2. TBR to new local optimum on original landscape
    TBRResult tbr_result = tbr_search(tree, ds, search_params, cd,
                                       nullptr, nullptr, check_timeout);
    total_moves += tbr_result.n_accepted;

    // impose_constraint() (above) is a heuristic that can fail to
    // converge on nested/complex constraints (it may cycle without
    // fully repairing all splits). TBR cannot self-repair a violated
    // constraint — regraft_violates_constraint() rejects all moves once
    // a split is unmapped. So a still-violating tree can carry an
    // improved score through TBR unchallenged. Verify before capturing,
    // mirroring the fused-tree check in ts_driven.cpp (search
    // "fused_ok").
    bool accept = tbr_result.best_score < best_score;
    if (accept && cd && cd->active) {
      map_constraint_nodes(tree, *cd);
      for (int _s = 0; _s < cd->n_splits; ++_s) {
        if (cd->constraint_node[_s] < 0) { accept = false; break; }
      }
    }

    if (accept) {
      best_score = tbr_result.best_score;
      best_tree = tree;
      ++n_escapes;
    } else {
      // Revert to best known tree
      copy_topology(tree, best_tree);
      tree.build_postorder();
      tree.reset_states(ds);
      // Re-sync constraint metadata after topology revert.
      // Same bug class as T-278 (TBR), T-279 (drift), F-015 (ratchet).
      if (cd) update_constraint(tree, *cd);
    }

    ++cycles_completed;

    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  // Ensure tree holds the best result
  if (cycles_completed > 0) {
    copy_topology(tree, best_tree);
    tree.build_postorder();
    tree.reset_states(ds);
  }

  return NNIPerturbResult{
    best_score,
    cycles_completed,
    total_moves,
    n_escapes
  };
}

} // namespace ts
