#include "ts_ratchet.h"
#include "ts_tbr.h"
#include "ts_fitch.h"

#include <random>
#include <algorithm>
#include <R.h>           // R_CheckUserInterrupt
#include <Rmath.h>       // unif_rand

namespace ts {

namespace {

// Save active_mask for every block so we can restore after perturbation.
std::vector<uint64_t> save_active_masks(const DataSet& ds) {
  std::vector<uint64_t> saved(ds.n_blocks);
  for (int b = 0; b < ds.n_blocks; ++b) {
    saved[b] = ds.blocks[b].active_mask;
  }
  return saved;
}

void restore_active_masks(DataSet& ds,
                          const std::vector<uint64_t>& saved) {
  for (int b = 0; b < ds.n_blocks; ++b) {
    ds.blocks[b].active_mask = saved[b];
  }
}

// Zero random character bits in active_mask across all blocks.
void perturb_characters(DataSet& ds, double prob, std::mt19937& rng) {
  std::bernoulli_distribution coin(prob);
  for (auto& blk : ds.blocks) {
    uint64_t mask = blk.active_mask;
    for (int i = 0; i < blk.n_chars; ++i) {
      if (coin(rng)) {
        mask &= ~(uint64_t(1) << i);
      }
    }
    blk.active_mask = mask;
  }
}

// Copy topology from src to dst (same-sized trees).
void copy_topology(TreeState& dst, const TreeState& src) {
  dst.parent = src.parent;
  dst.left   = src.left;
  dst.right  = src.right;
}

} // anonymous namespace

RatchetResult ratchet_search(TreeState& tree, DataSet& ds,
                             const RatchetParams& params) {
  // Initial TBR to get a baseline
  TBRParams search_params;
  search_params.accept_equal = false;
  search_params.max_accepted_changes = 0;
  search_params.max_hits = params.max_hits;

  TBRResult initial = tbr_search(tree, ds, search_params);

  double best_score = initial.best_score;
  int total_moves = initial.n_accepted;

  // Save best topology
  TreeState best_tree = tree;

  // Perturbation TBR params
  int perturb_max_changes = std::max(20,
                              std::min(200, tree.n_tip / 8));
  TBRParams perturb_params;
  perturb_params.accept_equal = true;
  perturb_params.max_accepted_changes = perturb_max_changes;
  perturb_params.max_hits = 1;

  // Seed from R's RNG for reproducibility with set.seed()
  std::mt19937 rng(static_cast<unsigned>(unif_rand() * 4294967295.0));

  int cycles_completed = 0;
  for (int cycle = 0; cycle < params.n_cycles; ++cycle) {
    // 1. Perturbation phase
    std::vector<uint64_t> saved_masks = save_active_masks(ds);
    perturb_characters(ds, params.perturb_prob, rng);

    TBRResult perturb_result = tbr_search(tree, ds, perturb_params);
    total_moves += perturb_result.n_accepted;

    // 2. Search phase: restore original weights, standard TBR
    restore_active_masks(ds, saved_masks);
    TBRResult search_result = tbr_search(tree, ds, search_params);
    total_moves += search_result.n_accepted;

    if (search_result.best_score < best_score) {
      best_score = search_result.best_score;
      best_tree = tree;
    } else {
      // Reset to best known tree
      copy_topology(tree, best_tree);
      tree.build_postorder();
      tree.reset_states(ds);
    }

    ++cycles_completed;
    R_CheckUserInterrupt();
  }

  // Ensure tree holds the best result
  if (cycles_completed > 0) {
    copy_topology(tree, best_tree);
    tree.build_postorder();
    tree.reset_states(ds);
  }

  return RatchetResult{best_score, cycles_completed, total_moves};
}

} // namespace ts
