#include "ts_driven.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_drift.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_wagner.h"
#include "ts_splits.h"
#include "ts_rng.h"

#include <R.h>
#include <chrono>
#include <functional>

namespace ts {

namespace {

// Build a ProgressInfo snapshot. pool may be null (during single-replicate).
ProgressInfo make_progress(int rep, const DrivenParams& params,
                           const TreePool* pool,
                           const char* phase, double elapsed,
                           double phase_score) {
  ProgressInfo pi;
  pi.replicate = rep;
  pi.max_replicates = params.max_replicates;
  pi.best_score = pool && pool->size() > 0 ? pool->best_score() : 1e18;
  pi.hits_to_best = pool ? pool->hits_to_best() : 0;
  pi.target_hits = params.target_hits;
  pi.pool_size = pool ? pool->size() : 0;
  pi.phase = phase;
  pi.elapsed_seconds = elapsed;
  pi.phase_score = phase_score;
  return pi;
}

} // anonymous namespace

// --- Single-replicate pipeline ---

ReplicateResult run_single_replicate(
    DataSet& ds,
    const DrivenParams& params,
    ConstraintData* cd,
    std::function<bool()> check_timeout,
    int verbosity)
{
  ReplicateResult result;
  result.interrupted = false;

  bool tree_large_enough_for_sectors =
      ds.n_tips >= 2 * params.sector_min_size;

  // 1. Random addition sequence Wagner tree (try multiple, keep best)
  random_wagner_tree(result.tree, ds, cd);
  double best_wag = score_tree(result.tree, ds);
  for (int ws = 1; ws < params.wagner_starts; ++ws) {
    TreeState trial;
    random_wagner_tree(trial, ds, cd);
    double trial_score = score_tree(trial, ds);
    if (trial_score < best_wag) {
      result.tree = std::move(trial);
      best_wag = trial_score;
    }
  }

  if (verbosity >= 2) {
    Rprintf("  Wagner tree score: %.1f%s\n", best_wag,
            params.wagner_starts > 1 ? " (best of multiple starts)" : "");
  }

  // 2. TBR to local optimum
  {
    TBRParams tp;
    tp.tabu_size = params.tabu_size;
    tbr_search(result.tree, ds, tp, cd);
  }
  if (verbosity >= 2) {
    Rprintf("  TBR score: %.1f\n", score_tree(result.tree, ds));
  }

  if (ts::check_interrupt() || check_timeout()) {
    result.interrupted = true;
    result.score = score_tree(result.tree, ds);
    return result;
  }

  // 3. Sectorial search (XSS + RSS) if tree is large enough
  if (tree_large_enough_for_sectors) {
    SectorParams sp;
    sp.min_sector_size = params.sector_min_size;
    sp.max_sector_size = params.sector_max_size;
    sp.internal_ratchet_cycles = 0;
    sp.internal_max_hits = 1;

    // XSS: systematic partitioning
    sp.n_partitions = params.xss_partitions;
    sp.xss_rounds = params.xss_rounds;
    xss_search(result.tree, ds, sp, cd);

    if (verbosity >= 2) {
      Rprintf("  XSS score: %.1f\n", score_tree(result.tree, ds));
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // RSS: random sector picks
    if (params.rss_rounds > 0) {
      for (int rr = 0; rr < params.rss_rounds; ++rr) {
        rss_search(result.tree, ds, sp, cd);

        if (ts::check_interrupt() || check_timeout()) {
          result.interrupted = true;
          result.score = score_tree(result.tree, ds);
          return result;
        }
      }
      if (verbosity >= 2) {
        Rprintf("  RSS score: %.1f\n", score_tree(result.tree, ds));
      }
    }

    // CSS: sector-restricted TBR on full tree (exact scoring)
    if (params.css_rounds > 0) {
      SectorParams css_sp;
      css_sp.n_partitions = params.css_partitions;
      css_sp.xss_rounds = params.css_rounds;
      css_sp.internal_max_hits = 1;
      css_search(result.tree, ds, css_sp, cd);

      if (verbosity >= 2) {
        Rprintf("  CSS score: %.1f\n", score_tree(result.tree, ds));
      }

      if (ts::check_interrupt() || check_timeout()) {
        result.interrupted = true;
        result.score = score_tree(result.tree, ds);
        return result;
      }
    }
  }

  if (ts::check_interrupt() || check_timeout()) {
    result.interrupted = true;
    result.score = score_tree(result.tree, ds);
    return result;
  }

  // 4. Ratchet perturbation to escape local optima
  {
    RatchetParams rp;
    rp.n_cycles = params.ratchet_cycles;
    rp.perturb_prob = params.ratchet_perturb_prob;
    rp.max_hits = params.tbr_max_hits;
    rp.perturb_mode = static_cast<PerturbMode>(params.ratchet_perturb_mode);
    rp.perturb_max_moves = params.ratchet_perturb_max_moves;
    rp.adaptive = params.ratchet_adaptive;
    rp.tabu_size = params.tabu_size;
    ratchet_search(result.tree, ds, rp, cd);
  }
  if (verbosity >= 2) {
    Rprintf("  Ratchet score: %.1f\n", score_tree(result.tree, ds));
  }

  if (ts::check_interrupt() || check_timeout()) {
    result.interrupted = true;
    result.score = score_tree(result.tree, ds);
    return result;
  }

  // 5. Drifting (suboptimal + equal-score exploration)
  if (params.drift_cycles > 0) {
    DriftParams dp;
    dp.n_cycles = params.drift_cycles;
    dp.afd_limit = params.drift_afd_limit;
    dp.rfd_limit = params.drift_rfd_limit;
    dp.max_hits = params.tbr_max_hits;
    dp.tabu_size = params.tabu_size;
    drift_search(result.tree, ds, dp, cd);

    if (verbosity >= 2) {
      Rprintf("  Drift score: %.1f\n", score_tree(result.tree, ds));
    }
  }

  if (ts::check_interrupt() || check_timeout()) {
    result.interrupted = true;
    result.score = score_tree(result.tree, ds);
    return result;
  }

  // 6. Final TBR polish
  {
    TBRParams tp;
    tp.tabu_size = params.tabu_size;
    tbr_search(result.tree, ds, tp, cd);
  }
  if (verbosity >= 2) {
    Rprintf("  Final TBR score: %.1f\n", score_tree(result.tree, ds));
  }

  result.score = score_tree(result.tree, ds);
  return result;
}

// --- Full driven search (serial) ---

DrivenResult driven_search(TreePool& pool, DataSet& ds,
                           const DrivenParams& params,
                           ConstraintData* cd) {
  DrivenResult result;
  result.best_score = 1e18;
  result.replicates_completed = 0;
  result.hits_to_best = 0;
  result.pool_size = 0;
  result.timed_out = false;

  if (params.max_replicates <= 0) {
    result.best_score = -1.0;
    return result;
  }

  bool use_timeout = params.max_seconds > 0.0;
  auto start_time = std::chrono::steady_clock::now();

  auto elapsed = [&]() -> double {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_time).count();
  };

  auto check_timeout = [&]() -> bool {
    if (!use_timeout) return false;
    return elapsed() >= params.max_seconds;
  };

  bool has_callback = static_cast<bool>(params.progress_callback);

  // Helper: report progress via callback or Rprintf fallback
  auto report = [&](const char* phase, int min_verbosity,
                    double phase_score, int rep_1based) {
    if (params.verbosity < min_verbosity) return;
    if (has_callback) {
      ProgressInfo pi = make_progress(rep_1based, params, &pool,
                                       phase, elapsed(), phase_score);
      params.progress_callback(pi);
    }
    // Rprintf fallback is handled inline at callsites (only for
    // messages that don't map cleanly to a single Rprintf call)
  };

  for (int rep = 0; rep < params.max_replicates; ++rep) {
    int rep1 = rep + 1;

    if (params.verbosity >= 1) {
      if (has_callback) {
        // Callback gets a "replicate_start" event; not a separate phase
        // — the per-replicate "replicate" event is emitted after pool add.
      } else {
        Rprintf("Replicate %d/%d", rep1, params.max_replicates);
        if (pool.size() > 0) {
          Rprintf(" (best: %.1f, pool: %d, hits: %d)",
                  pool.best_score(), pool.size(), pool.hits_to_best());
        }
        Rprintf("\n");
      }
    }

    // Run the single-replicate pipeline
    ReplicateResult rep_result = run_single_replicate(
        ds, params, cd, check_timeout, params.verbosity);

    if (rep_result.interrupted) {
      if (rep_result.score < 1e18) {
        pool.add(rep_result.tree, rep_result.score);
      }
      result.timed_out = true;
      goto finish;
    }

    // Add to pool
    pool.add(rep_result.tree, rep_result.score);

    ++result.replicates_completed;

    // Report end of replicate
    report("replicate", 1, rep_result.score, rep1);

    // Periodic tree fusing
    if ((rep + 1) % params.fuse_interval == 0 && pool.size() >= 2) {
      int hits_before = pool.hits_to_best();
      double best_before = pool.best_score();

      TreeState fused = pool.best().tree;
      FuseParams fp;
      fp.accept_equal = params.fuse_accept_equal;
      fp.max_rounds = 10;
      tree_fuse(fused, ds, pool, fp);

      double fused_score = score_tree(fused, ds);

      bool fuse_ok = true;
      if (cd && cd->active) {
        fuse_ok = !violates_constraint_posthoc(fused, *cd);
      }
      if (fuse_ok) pool.add(fused, fused_score);

      if (fused_score < best_before) {
        pool.set_hits_to_best(0);
        if (params.verbosity >= 1) {
          if (has_callback) {
            report("fuse", 1, fused_score, rep1);
          } else {
            Rprintf("  Fuse improved: %.1f -> %.1f\n",
                    best_before, fused_score);
          }
        }
      } else {
        pool.set_hits_to_best(hits_before);
      }
    }

    // Convergence check
    if (pool.hits_to_best() >= params.target_hits) {
      if (params.verbosity >= 1) {
        if (!has_callback) {
          Rprintf("Converged: %d hits to best score %.1f\n",
                  pool.hits_to_best(), pool.best_score());
        }
      }
      break;
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.timed_out = true;
      goto finish;
    }
  }

finish:
  result.hits_to_best = pool.hits_to_best();
  result.pool_size = pool.size();

  if (pool.size() > 0) {
    result.best_score = pool.best_score();
  } else {
    result.best_score = -1.0;
  }

  // Final "done" callback
  if (has_callback && params.verbosity >= 1) {
    ProgressInfo pi = make_progress(result.replicates_completed, params,
                                     &pool, "done", elapsed(),
                                     result.best_score);
    pi.replicate = result.replicates_completed;
    params.progress_callback(pi);
  } else if (result.timed_out && params.verbosity >= 1) {
    Rprintf("Timeout reached (%.1f s)\n", params.max_seconds);
  }

  return result;
}

} // namespace ts
