#include "ts_driven.h"
#include "ts_fitch.h"
#include "ts_search.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_drift.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_wagner.h"
#include "ts_splits.h"
#include "ts_rng.h"

#include <R.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
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
    int verbosity,
    TreeState* starting_tree)
{
  ReplicateResult result;
  result.interrupted = false;

  bool tree_large_enough_for_sectors =
      ds.n_tips >= 2 * params.sector_min_size;

  using PhClock = std::chrono::steady_clock;
  auto ph_start = PhClock::now();
  auto ph_lap = [&]() {
    auto now = PhClock::now();
    double ms = std::chrono::duration<double, std::milli>(now - ph_start).count();
    ph_start = now;
    return ms;
  };

  // 1. Starting tree: use provided tree or build Wagner
  double best_wag;
  if (starting_tree) {
    result.tree = *starting_tree;
    best_wag = score_tree(result.tree, ds);
  } else {
    random_wagner_tree(result.tree, ds, cd);
    best_wag = score_tree(result.tree, ds);
    for (int ws = 1; ws < params.wagner_starts; ++ws) {
      TreeState trial;
      random_wagner_tree(trial, ds, cd);
      double trial_score = score_tree(trial, ds);
      if (trial_score < best_wag) {
        result.tree = std::move(trial);
        best_wag = trial_score;
      }
    }
  }

  result.timings.wagner_ms = ph_lap();
  if (verbosity >= 2) {
    if (starting_tree) {
      Rprintf("  Starting tree score: %.5g [%.0f ms]\n", best_wag,
              result.timings.wagner_ms);
    } else {
      Rprintf("  Wagner tree score: %.5g [%.0f ms]%s\n", best_wag,
              result.timings.wagner_ms,
              params.wagner_starts > 1 ? " (best of multiple starts)" : "");
    }
  }

  // 2. Hill-climbing to local optimum (SPR→TBR escalation)
  if (params.spr_first) {
    spr_search(result.tree, ds, 1);
  }
  {
    TBRParams tp;
    tp.tabu_size = params.tabu_size;
    tbr_search(result.tree, ds, tp, cd);
  }
  result.timings.tbr_ms = ph_lap();
  if (verbosity >= 2) {
    Rprintf("  TBR score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
            result.timings.tbr_ms);
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

    result.timings.xss_ms = ph_lap();
    if (verbosity >= 2) {
      Rprintf("  XSS score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
              result.timings.xss_ms);
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
      result.timings.rss_ms = ph_lap();
      if (verbosity >= 2) {
        Rprintf("  RSS score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
                result.timings.rss_ms);
      }
    }

    // CSS: sector-restricted TBR on full tree (exact scoring)
    if (params.css_rounds > 0) {
      SectorParams css_sp;
      css_sp.n_partitions = params.css_partitions;
      css_sp.xss_rounds = params.css_rounds;
      css_sp.internal_max_hits = 1;
      css_search(result.tree, ds, css_sp, cd);

      result.timings.css_ms = ph_lap();
      if (verbosity >= 2) {
        Rprintf("  CSS score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
                result.timings.css_ms);
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
    ratchet_search(result.tree, ds, rp, cd, check_timeout);
  }
  result.timings.ratchet_ms = ph_lap();
  if (verbosity >= 2) {
    Rprintf("  Ratchet score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
            result.timings.ratchet_ms);
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
    drift_search(result.tree, ds, dp, cd, check_timeout);

    result.timings.drift_ms = ph_lap();
    if (verbosity >= 2) {
      Rprintf("  Drift score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
              result.timings.drift_ms);
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
  result.timings.final_tbr_ms = ph_lap();
  if (verbosity >= 2) {
    Rprintf("  Final TBR score: %.5g [%.0f ms]\n", score_tree(result.tree, ds),
            result.timings.final_tbr_ms);
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

  // Cancel file: read path from environment variable (set by Shiny app).
  std::string cancel_path;
  {
    const char* cancel_env = std::getenv("TREESEARCH_CANCEL_FILE");
    if (cancel_env && cancel_env[0] != '\0') cancel_path = cancel_env;
  }
  auto check_cancel = [&]() -> bool {
    if (cancel_path.empty()) return false;
    FILE* cf = std::fopen(cancel_path.c_str(), "r");
    if (cf) {
      std::fclose(cf);
      return true;
    }
    return false;
  };

  auto elapsed = [&]() -> double {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_time).count();
  };

  // Combines timeout + cancel-file check so both are evaluated inside
  // run_single_replicate() (which only receives check_timeout).
  auto check_timeout = [&]() -> bool {
    if (use_timeout && elapsed() >= params.max_seconds) return true;
    return check_cancel();
  };

  bool has_callback = static_cast<bool>(params.progress_callback);

  // Helper: report progress via callback or Rprintf fallback.
  // Callbacks are ALWAYS invoked when present (regardless of verbosity)
  // so that Shiny progress polling works at verbosity=0.
  auto report = [&](const char* phase, int min_verbosity,
                    double phase_score, int rep_1based) {
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
          Rprintf(" (best: %.5g, pool: %d, hits: %d)",
                  pool.best_score(), pool.size(), pool.hits_to_best());
        }
        Rprintf("\n");
      }
    }

    // Use starting tree for replicate 0 if provided
    TreeState* start_ptr = nullptr;
    TreeState start_tree;
    if (rep == 0 && params.start_n_edge > 0 &&
        static_cast<int>(params.start_edge.size()) >= 2 * params.start_n_edge) {
      const int* edge_parent = params.start_edge.data();
      const int* edge_child = params.start_edge.data() + params.start_n_edge;
      start_tree.init_from_edge(edge_parent, edge_child,
                                params.start_n_edge, ds);
      start_ptr = &start_tree;
    }

    // Run the single-replicate pipeline
    ReplicateResult rep_result = run_single_replicate(
        ds, params, cd, check_timeout, params.verbosity, start_ptr);

    result.timings += rep_result.timings;

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
    if (params.fuse_interval > 0 &&
        (rep + 1) % params.fuse_interval == 0 && pool.size() >= 2) {
      auto fuse_start = std::chrono::steady_clock::now();

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
        report("fuse", 1, fused_score, rep1);
        if (params.verbosity >= 1 && !has_callback) {
          Rprintf("  Fuse improved: %.5g -> %.5g\n",
                  best_before, fused_score);
        }
      } else {
        pool.set_hits_to_best(hits_before);
      }

      auto fuse_end = std::chrono::steady_clock::now();
      result.timings.fuse_ms +=
          std::chrono::duration<double, std::milli>(fuse_end - fuse_start).count();
    }

    // Convergence check
    if (pool.hits_to_best() >= params.target_hits) {
      if (params.verbosity >= 1) {
        if (!has_callback) {
          Rprintf("Converged: %d hits to best score %.5g\n",
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

  // Capture hits_to_best BEFORE MPT enumeration.  Only main-loop
  // replicates are independent; MPT enumeration discovers variant
  // topologies from existing pool trees and should not inflate the count.
  result.hits_to_best = pool.hits_to_best();

  // 7. MPT enumeration: TBR plateau walk from each pool tree to discover
  //    additional equal-score topologies.  Each replicate contributes only
  //    one tree, so different TBR-connected islands are only discovered if
  //    different replicates landed on them.  We enumerate from each seed
  //    tree to explore its island, stopping when the pool is full.
  if (pool.size() > 0 && pool.size() < pool.max_size && !result.timed_out) {
    TBRParams tp;
    tp.accept_equal = true;
    tp.tabu_size = params.tabu_size > 0 ? params.tabu_size : 100;

    // Snapshot current pool entries as seeds (new trees discovered during
    // enumeration of seed i become additional seeds for later iterations).
    int seed_idx = 0;
    while (seed_idx < pool.size() && pool.size() < pool.max_size) {
      TreeState enum_tree = pool.all()[seed_idx].tree;
      // Budget remaining capacity across remaining seeds
      tp.max_hits = std::max(10, (pool.max_size - pool.size()) * 2);
      tbr_search(enum_tree, ds, tp, cd, nullptr, &pool);
      ++seed_idx;
    }
    if (params.verbosity >= 2) {
      Rprintf("MPT enumeration: %d trees in pool\n", pool.size());
    }
  }

  // result.hits_to_best already set before MPT enumeration
  result.pool_size = pool.size();

  if (pool.size() > 0) {
    result.best_score = pool.best_score();
  } else {
    result.best_score = -1.0;
  }

  // Final "done" callback (always fired when callback exists)
  if (has_callback) {
    ProgressInfo pi = make_progress(result.replicates_completed, params,
                                     &pool, "done", elapsed(),
                                     result.best_score);
    pi.replicate = result.replicates_completed;
    params.progress_callback(pi);
  } else if (result.timed_out && params.verbosity >= 1) {
    Rprintf("Timeout reached (%.5g s)\n", params.max_seconds);
  }

  return result;
}

} // namespace ts
