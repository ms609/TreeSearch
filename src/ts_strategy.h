#ifndef TS_STRATEGY_H
#define TS_STRATEGY_H

// Adaptive starting-tree strategy selection via Thompson sampling.
//
// Each replicate draws its starting-tree strategy from a probability
// distribution maintained by a Beta-Bernoulli multi-armed bandit.
// The reward signal is whether the replicate hit the pool's best score.
//
// Strategy arms:
//   WAGNER_RANDOM       - Random addition-order Wagner (baseline)
//   WAGNER_GOLOBOFF     - Goloboff (2014) non-ambiguous-char bias
//   WAGNER_ENTROPY      - State-specificity bias
//   RANDOM_TREE         - Purely random topology (no character data)
//   POOL_RATCHET        - Ratchet perturbation from a pool tree
//   POOL_NNI_PERTURB    - NNI perturbation from a pool tree
//
// Pool-based arms are only available once the pool has ≥2 best-score trees.
// RANDOM_TREE starts with a pessimistic prior (Beta(1,2)) reflecting the
// expectation that it's usually worse, but letting data override.
//
// On new best score, all counts are decayed by 0.5× to discount stale evidence.
//
// Reference: Thompson (1933), "On the likelihood that one unknown probability
//            exceeds another in view of the evidence of two samples."

#include <random>
#include <array>
#include <algorithm>

namespace ts {

enum class StartStrategy : int {
  WAGNER_RANDOM     = 0,
  WAGNER_GOLOBOFF   = 1,
  WAGNER_ENTROPY    = 2,
  RANDOM_TREE       = 3,
  POOL_RATCHET      = 4,
  POOL_NNI_PERTURB  = 5,
  N_STRATEGIES      = 6
};

inline const char* strategy_name(StartStrategy s) {
  switch (s) {
    case StartStrategy::WAGNER_RANDOM:    return "wag_rand";
    case StartStrategy::WAGNER_GOLOBOFF:  return "wag_golob";
    case StartStrategy::WAGNER_ENTROPY:   return "wag_entropy";
    case StartStrategy::RANDOM_TREE:      return "rand_tree";
    case StartStrategy::POOL_RATCHET:     return "pool_ratch";
    case StartStrategy::POOL_NNI_PERTURB: return "pool_nni";
    default:                              return "unknown";
  }
}

constexpr int N_STRAT = static_cast<int>(StartStrategy::N_STRATEGIES);

// Returns true if this strategy requires a populated pool.
inline bool strategy_needs_pool(StartStrategy s) {
  return s == StartStrategy::POOL_RATCHET ||
         s == StartStrategy::POOL_NNI_PERTURB;
}

// Returns true if this is a Wagner-based strategy.
inline bool strategy_is_wagner(StartStrategy s) {
  return s == StartStrategy::WAGNER_RANDOM ||
         s == StartStrategy::WAGNER_GOLOBOFF ||
         s == StartStrategy::WAGNER_ENTROPY;
}

class StrategyTracker {
public:
  // Initialise with default priors.
  // RANDOM_TREE gets Beta(1,2) = pessimistic prior.
  // All others get Beta(1,1) = uniform.
  StrategyTracker() {
    for (int i = 0; i < N_STRAT; ++i) {
      alpha_[i] = 1.0;
      beta_[i] = 1.0;
      attempts_[i] = 0;
      successes_[i] = 0;
    }
    // Pessimistic prior for random tree
    beta_[static_cast<int>(StartStrategy::RANDOM_TREE)] = 2.0;
  }

  // Select a strategy via Thompson sampling.
  // `pool_available` gates pool-based arms.
  // `rng` is the caller's RNG (not R's — safe for parallel use).
  StartStrategy select(std::mt19937& rng, bool pool_available) const {
    double best_sample = -1.0;
    StartStrategy best = StartStrategy::WAGNER_RANDOM;

    for (int i = 0; i < N_STRAT; ++i) {
      auto s = static_cast<StartStrategy>(i);
      if (strategy_needs_pool(s) && !pool_available) continue;

      // Sample from Beta(alpha, beta) via two Gamma draws
      std::gamma_distribution<double> ga(alpha_[i], 1.0);
      std::gamma_distribution<double> gb(beta_[i], 1.0);
      double x = ga(rng);
      double y = gb(rng);
      double theta = (x + y > 0.0) ? x / (x + y) : 0.5;

      if (theta > best_sample) {
        best_sample = theta;
        best = s;
      }
    }
    return best;
  }

  // Update after a replicate completes.
  // `hit_best` = true if the replicate matched or improved the pool best score.
  void update(StartStrategy strategy, bool hit_best) {
    int i = static_cast<int>(strategy);
    attempts_[i]++;
    if (hit_best) {
      successes_[i]++;
      alpha_[i] += 1.0;
    } else {
      beta_[i] += 1.0;
    }
  }

  // Decay all counts when the best score improves.
  // Old evidence is stale because the landscape effectively changed.
  void decay(double factor = 0.5) {
    for (int i = 0; i < N_STRAT; ++i) {
      // Decay toward prior (don't let alpha/beta go below 1.0)
      alpha_[i] = std::max(1.0, 1.0 + (alpha_[i] - 1.0) * factor);
      beta_[i]  = std::max(1.0, 1.0 + (beta_[i]  - 1.0) * factor);
    }
  }

  // Accessors for diagnostics / reporting
  int attempts(StartStrategy s) const { return attempts_[static_cast<int>(s)]; }
  int successes(StartStrategy s) const { return successes_[static_cast<int>(s)]; }
  double alpha(StartStrategy s) const { return alpha_[static_cast<int>(s)]; }
  double beta_param(StartStrategy s) const { return beta_[static_cast<int>(s)]; }

  // Pre-compute a round-robin strategy sequence for parallel dispatch.
  // Returns a vector of length `n_replicates` cycling through all arms.
  // Pool-based arms are only included from `pool_available_after` onwards.
  static std::vector<StartStrategy> round_robin(
      int n_replicates, int pool_available_after = 3)
  {
    std::vector<StartStrategy> seq(n_replicates);
    // Fresh-only arms for early replicates
    const StartStrategy fresh[] = {
      StartStrategy::WAGNER_RANDOM,
      StartStrategy::WAGNER_GOLOBOFF,
      StartStrategy::WAGNER_ENTROPY,
      StartStrategy::RANDOM_TREE,
    };
    constexpr int n_fresh = 4;

    // All arms for later replicates
    const StartStrategy all[] = {
      StartStrategy::WAGNER_RANDOM,
      StartStrategy::WAGNER_GOLOBOFF,
      StartStrategy::WAGNER_ENTROPY,
      StartStrategy::RANDOM_TREE,
      StartStrategy::POOL_RATCHET,
      StartStrategy::POOL_NNI_PERTURB,
    };
    constexpr int n_all = N_STRAT;

    int fresh_idx = 0, all_idx = 0;
    for (int r = 0; r < n_replicates; ++r) {
      if (r < pool_available_after) {
        seq[r] = fresh[fresh_idx % n_fresh];
        ++fresh_idx;
      } else {
        seq[r] = all[all_idx % n_all];
        ++all_idx;
      }
    }
    return seq;
  }

private:
  std::array<double, N_STRAT> alpha_;  // Beta distribution alpha param
  std::array<double, N_STRAT> beta_;   // Beta distribution beta param
  std::array<int, N_STRAT> attempts_;
  std::array<int, N_STRAT> successes_;
};

} // namespace ts

#endif // TS_STRATEGY_H
