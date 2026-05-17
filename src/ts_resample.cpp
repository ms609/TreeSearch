#include "ts_resample.h"
#include "ts_fitch.h"
#include "ts_tree.h"
#include "ts_rng.h"

#include <R.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace ts {

// ---- Resampling (jackknife / bootstrap) ----

ResampleResult resample_search(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* original_weights,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const ResampleParams& params,
    const double* info_amounts_r,
    int info_max_steps,
    ConstraintData* cd,
    bool xpiwe,
    double xpiwe_r,
    double xpiwe_max_f,
    const int* obs_count_r)
{
  // Expand original weights into a flat character index
  // (each pattern p appears original_weights[p] times).
  //
  // Guard: n_total_chars is cast to int and used as an array-index bound.
  // If sum(weights) > INT_MAX the cast overflows to a negative value and
  // the subsequent array access is undefined behaviour (SIGSEGV).  Check
  // via size_t arithmetic before allocating and error out cleanly.
  {
    size_t total_chars = 0;
    for (int p = 0; p < n_patterns; ++p) {
      if (original_weights[p] < 0) {
        Rf_error("TreeSearch: character weight[%d] = %d is negative",
                 p, original_weights[p]);
      }
      total_chars += static_cast<size_t>(original_weights[p]);
    }
    if (total_chars > static_cast<size_t>(INT_MAX)) {
      Rf_error("TreeSearch: sum of character weights (%zu) exceeds INT_MAX.\n"
               "  Reduce options(\"TreeSearch.fractional.scale\") or set\n"
               "  weights to smaller values before calling Resample().",
               total_chars);
    }
  }

  std::vector<int> char_index;
  for (int p = 0; p < n_patterns; ++p) {
    for (int w = 0; w < original_weights[p]; ++w) {
      char_index.push_back(p);
    }
  }

  int n_total_chars = static_cast<int>(char_index.size());
  std::vector<int> new_weights(n_patterns, 0);

  ts::rng_state_begin();
  if (params.bootstrap) {
    // Sample WITH replacement
    for (int i = 0; i < n_total_chars; ++i) {
      int j = static_cast<int>(ts::thread_safe_unif() * n_total_chars);
      if (j >= n_total_chars) j = n_total_chars - 1;
      new_weights[char_index[j]]++;
    }
  } else {
    // Jackknife: sample WITHOUT replacement (Fisher-Yates partial shuffle)
    int n_keep = static_cast<int>(
        std::ceil(params.jack_proportion * n_total_chars));
    if (n_keep < 1) n_keep = 1;
    if (n_keep >= n_total_chars) n_keep = n_total_chars - 1;

    // Partial Fisher-Yates: shuffle first n_keep elements
    std::vector<int> indices(char_index);
    for (int i = 0; i < n_keep; ++i) {
      int j = i + static_cast<int>(ts::thread_safe_unif() * (n_total_chars - i));
      if (j >= n_total_chars) j = n_total_chars - 1;
      std::swap(indices[i], indices[j]);
    }
    for (int i = 0; i < n_keep; ++i) {
      new_weights[indices[i]]++;
    }
  }
  ts::rng_state_end();

  // Build a DataSet with the resampled weights
  DataSet ds = build_dataset(
      contrast_r, n_tokens, n_states,
      tip_data_r, n_tips, n_patterns,
      new_weights.data(),
      levels_r,
      min_steps_r,
      concavity,
      info_amounts_r,
      info_max_steps,
      xpiwe,
      xpiwe_r,
      xpiwe_max_f,
      obs_count_r);

  // Run driven search
  TreePool pool(params.search.pool_max_size, params.search.pool_suboptimal);
  DrivenResult dr = driven_search(pool, ds, params.search, cd);

  // Extract result
  ResampleResult result;
  result.n_tip = n_tips;
  result.score = dr.best_score;

  if (pool.size() > 0) {
    const TreeState& best = pool.best().tree;
    int n_edge = 2 * (best.n_tip - 1);
    result.edge_parent.resize(n_edge);
    result.edge_child.resize(n_edge);
    int row = 0;
    for (int node = best.n_tip; node < best.n_node; ++node) {
      int ni = node - best.n_tip;
      result.edge_parent[row] = node + 1;
      result.edge_child[row] = best.left[ni] + 1;
      ++row;
      result.edge_parent[row] = node + 1;
      result.edge_child[row] = best.right[ni] + 1;
      ++row;
    }
  }

  return result;
}

// ---- Successive Approximations ----

SAResult successive_approximations(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* original_weights,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const SAParams& params,
    const double* info_amounts_r,
    int info_max_steps,
    ConstraintData* cd,
    bool xpiwe,
    double xpiwe_r,
    double xpiwe_max_f,
    const int* obs_count_r)
{
  SAResult result;
  result.n_tip = n_tips;
  result.score = -1.0;
  result.sa_iterations = 0;
  result.converged = false;

  int n_internal = n_tips - 1;

  // SA weights start at 1.0 for all patterns
  std::vector<double> sa_weights(n_patterns, 1.0);

  // Track previous iteration's per-pattern steps for convergence
  std::vector<int> prev_steps(n_patterns, -1);

  for (int iter = 0; iter < params.max_sa_iter; ++iter) {
    if (ts::check_interrupt()) break;

    // Build effective weights: original_weight * sa_weight (rounded to int)
    // SA weighting modifies the relative importance, but we need integer
    // weights for build_dataset. Multiply original_weight * sa_weight and round.
    // To preserve precision, we scale so the minimum nonzero weight is 1.
    std::vector<double> effective_d(n_patterns);
    for (int p = 0; p < n_patterns; ++p) {
      effective_d[p] = original_weights[p] * sa_weights[p];
    }

    // Find minimum nonzero effective weight for scaling
    double min_nz = 1e18;
    for (int p = 0; p < n_patterns; ++p) {
      if (effective_d[p] > 1e-12 && effective_d[p] < min_nz) {
        min_nz = effective_d[p];
      }
    }
    if (min_nz < 1e-12) min_nz = 1.0;

    std::vector<int> effective_weights(n_patterns);
    for (int p = 0; p < n_patterns; ++p) {
      effective_weights[p] = std::max(0,
          static_cast<int>(std::round(effective_d[p] / min_nz)));
    }

    // Build dataset with effective weights
    DataSet ds = build_dataset(
        contrast_r, n_tokens, n_states,
        tip_data_r, n_tips, n_patterns,
        effective_weights.data(),
        levels_r,
        min_steps_r,
        concavity,
        info_amounts_r,
        info_max_steps,
        xpiwe,
        xpiwe_r,
        xpiwe_max_f,
        obs_count_r);

    // Run driven search
    TreePool pool(params.search.pool_max_size, params.search.pool_suboptimal);
    driven_search(pool, ds, params.search, cd);

    ++result.sa_iterations;

    if (pool.size() == 0) break;

    // Get the best tree
    const TreeState& best_tree_state = pool.best().tree;

    // Extract per-pattern step counts using EW scoring
    // Rebuild dataset with original weights (weight=1 per pattern) to
    // get unweighted per-character steps for the SA reweighting formula
    DataSet ds_ew = build_dataset(
        contrast_r, n_tokens, n_states,
        tip_data_r, n_tips, n_patterns,
        original_weights,
        levels_r,
        min_steps_r,
        HUGE_VAL);

    // Re-initialize the tree against the EW dataset (different block layout)
    // Extract edge matrix first, then rebuild
    int n_edge = 2 * (best_tree_state.n_tip - 1);
    std::vector<int> ep(n_edge), ec(n_edge);
    {
      int row = 0;
      for (int node = best_tree_state.n_tip;
           node < best_tree_state.n_node; ++node) {
        int ni = node - best_tree_state.n_tip;
        ep[row] = node + 1;  // 1-based for init_from_edge
        ec[row] = best_tree_state.left[ni] + 1;
        ++row;
        ep[row] = node + 1;
        ec[row] = best_tree_state.right[ni] + 1;
        ++row;
      }
    }
    TreeState tree_copy;
    tree_copy.init_from_edge(ep.data(), ec.data(), n_edge, ds_ew);

    double ew_score;
    if (ds_ew.blocks.size() > 0 &&
        std::any_of(ds_ew.blocks.begin(), ds_ew.blocks.end(),
                    [](const CharBlock& b) { return b.has_inapplicable; })) {
      ew_score = static_cast<double>(fitch_na_score(tree_copy, ds_ew));
    } else {
      ew_score = static_cast<double>(fitch_score(tree_copy, ds_ew));
    }

    // Add back topology-independent steps for correct EW total
    result.score = ew_score + ds_ew.ew_offset;

    // Extract per-pattern step counts (reduced by simplification)
    std::vector<int> char_steps(n_patterns, 0);
    extract_char_steps(tree_copy, ds_ew, char_steps);

    // Add back precomputed_steps for SA reweighting — the ratio
    // p_i = steps / (n_internal - 1) should reflect total steps
    // including topology-independent autapomorphies/singletons.
    for (int p = 0; p < n_patterns; ++p) {
      if (!ds_ew.precomputed_steps.empty()) {
        char_steps[p] += ds_ew.precomputed_steps[p];
      }
    }

    // Check convergence: same steps as previous iteration
    if (char_steps == prev_steps) {
      result.converged = true;
      result.edge_parent = ep;
      result.edge_child = ec;
      break;
    }

    prev_steps = char_steps;

    // Reweight: w_i = (p_i)^(-k) - 1
    // where p_i = steps_i / (n_internal - 1)
    // Characters with 0 steps get maximum weight.
    for (int p = 0; p < n_patterns; ++p) {
      if (char_steps[p] <= 0) {
        // Perfect character: maximum weight
        sa_weights[p] = std::pow(1.0 / (n_internal - 1), -params.k) - 1.0;
      } else {
        double p_i = static_cast<double>(char_steps[p])
                     / static_cast<double>(n_internal - 1);
        if (p_i >= 1.0) {
          sa_weights[p] = 0.0;
        } else {
          sa_weights[p] = std::pow(p_i, -params.k) - 1.0;
        }
      }
      if (sa_weights[p] < 0.0) sa_weights[p] = 0.0;
    }

    // Save edge matrix (reuse ep/ec from above; may be overwritten next iter)
    result.edge_parent = ep;
    result.edge_child = ec;
  }

  return result;
}

} // namespace ts
