#ifndef TS_RESAMPLE_H
#define TS_RESAMPLE_H

// Resampling and successive approximations for the C++ search engine.
//
// - Jackknife: subsample characters (without replacement)
// - Bootstrap: resample characters (with replacement)
// - Successive approximations: iterative reweighting (Farris 1969)

#include "ts_data.h"
#include "ts_driven.h"
#include "ts_constraint.h"
#include <vector>

namespace ts {

// ---- Resampling (jackknife / bootstrap) ----

struct ResampleParams {
  // Resampling
  bool bootstrap = false;       // false = jackknife, true = bootstrap
  double jack_proportion = 2.0 / 3.0; // jackknife: proportion to keep

  // Driven search params for each replicate
  DrivenParams search;
};

struct ResampleResult {
  // Edge matrix (flattened, 2 columns) for the best tree found
  std::vector<int> edge_parent;
  std::vector<int> edge_child;
  int n_tip;
  double score;
};

// Run one resampled search iteration.
// Modifies weights (jackknife or bootstrap), rebuilds DataSet, runs search.
// `original_weights` is the pattern-frequency vector (length n_patterns).
// Returns the best tree from the resampled search.
ResampleResult resample_search(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* original_weights,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const ResampleParams& params,
    const double* info_amounts_r = nullptr,
    int info_max_steps = 0,
    ConstraintData* cd = nullptr);

// ---- Successive Approximations ----

struct SAParams {
  double k = 3.0;               // SA weighting constant (>= 1)
  int max_sa_iter = 20;         // maximum SA iterations
  DrivenParams search;          // driven search params per iteration
};

struct SAResult {
  // Best tree from the final iteration
  std::vector<int> edge_parent;
  std::vector<int> edge_child;
  int n_tip;
  double score;                 // EW parsimony score of best tree
  int sa_iterations;            // number of SA iterations completed
  bool converged;               // true if weights stabilized
};

// Run successive approximations search.
// Iteratively reweights characters based on fit, searching under reweighted
// parsimony until the optimal tree stabilizes.
SAResult successive_approximations(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* original_weights,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const SAParams& params,
    const double* info_amounts_r = nullptr,
    int info_max_steps = 0,
    ConstraintData* cd = nullptr);

} // namespace ts

#endif // TS_RESAMPLE_H
