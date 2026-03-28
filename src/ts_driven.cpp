#include "ts_driven.h"
#include "ts_fitch.h"
#include "ts_search.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_nni_perturb.h"
#include "ts_drift.h"
#include "ts_temper.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_pool.h"
#include "ts_constraint.h"
#include "ts_wagner.h"
#include "ts_splits.h"
#include "ts_prune_reinsert.h"
#include "ts_rng.h"

#include <R.h>
#include <Rmath.h>
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
    TreeState* starting_tree,
    const SplitFrequencyTable* split_freq,
    StartStrategy strategy,
    const TreePool* pool)
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

  // NNI warmup per Wagner start is skipped when constraints are active
  // because nni_search() does not enforce topological constraints.
  bool nni_wagner = params.nni_first && (!cd || !cd->active);

  // 1. Starting tree: dispatch on StartStrategy.
  //
  // All arms build a fresh tree from scratch, ensuring each replicate is
  // an independent sample from the landscape (basin coverage).
  //
  // When nni_first is true, NNI-optimize each start before selecting the
  // best. This finds better starting basins for TBR (O(n) per pass).
  double best_wag;
  if (starting_tree) {
    // User-supplied starting tree
    result.tree = *starting_tree;
    best_wag = score_tree(result.tree, ds);
  } else {

    // Build first start according to selected strategy
    switch (strategy) {
      case StartStrategy::WAGNER_GOLOBOFF: {
        ts::BiasedWagnerParams bwp;
        bwp.bias        = ts::WagnerBias::GOLOBOFF;
        bwp.temperature = params.wagner_bias_temp;
        biased_wagner_tree(result.tree, ds, bwp, cd);
        break;
      }
      case StartStrategy::WAGNER_ENTROPY: {
        ts::BiasedWagnerParams bwp;
        bwp.bias        = ts::WagnerBias::ENTROPY;
        bwp.temperature = params.wagner_bias_temp;
        biased_wagner_tree(result.tree, ds, bwp, cd);
        break;
      }
      case StartStrategy::RANDOM_TREE:
        if (cd && cd->active) {
          random_constrained_tree(result.tree, ds, *cd);
        } else {
          random_topology_tree(result.tree, ds);
        }
        break;
      default:  // WAGNER_RANDOM (and pool-based fallback)
        random_wagner_tree(result.tree, ds, cd);
        break;
    }
    best_wag = score_tree(result.tree, ds);
    if (nni_wagner) {
      auto nr = nni_search(result.tree, ds, 0, check_timeout);
      best_wag = nr.score;
    }
    // Additional Wagner starts: always random-order for basin diversity
    for (int ws = 1; ws < params.wagner_starts; ++ws) {
      TreeState trial;
      random_wagner_tree(trial, ds, cd);
      double trial_score = score_tree(trial, ds);
      if (nni_wagner) {
        auto nr = nni_search(trial, ds, 0, check_timeout);
        trial_score = nr.score;
      }
      if (trial_score < best_wag) {
        result.tree = std::move(trial);
        best_wag = trial_score;
      }
    }
  }

  result.timings.wagner_ms = ph_lap();
  if (verbosity >= 2) {
    if (starting_tree) {
      // User-supplied or warm-start tree
      Rprintf("  Starting tree score: %.5g [%.0f ms]\n", best_wag,
              result.timings.wagner_ms);
    } else {
      Rprintf("  %s%s tree score: %.5g [%.0f ms]%s\n",
              strategy_name(strategy),
              params.nni_first ? "+NNI" : "",
              best_wag, result.timings.wagner_ms,
              params.wagner_starts > 1 ? " (best of multiple starts)" : "");
    }
  }

  // 2. Hill-climbing to local optimum.
  // When NNI is active and unconstrained, each Wagner start was already
  // NNI-optimized, and SPR is skipped (NNI→TBR outperforms NNI→SPR→TBR).
  // When constrained, NNI was skipped above; fall back to SPR warmup.
  if (!nni_wagner && params.spr_first) {
    spr_search(result.tree, ds, 1, check_timeout);
  }
  {
    TBRParams tp;
    tp.tabu_size = params.tabu_size;
    tbr_search(result.tree, ds, tp, cd, nullptr, nullptr, check_timeout);
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

  // Steps 3–6: outer cycle loop.
  // Each outer cycle runs [XSS + RSS + CSS → Ratchet → NNI-perturb →
  // Drift → TBR polish] once.  With outer_cycles = 1 and
  // max_outer_resets = 0 (the defaults), this is exactly one pass through
  // the pipeline.  With outer_cycles > 1 we interleave fresh XSS passes
  // after each ratchet/drift escape, matching TNT's xmult pattern
  // (Goloboff 1999 §2.3).
  //
  // When max_outer_resets > 0 (or -1 = unlimited), a cycle that improves
  // the score resets the counter to exploit the new basin.  Resets are
  // capped to avoid runaway expansion on datasets with many small
  // incremental improvements (late resets have diminishing returns:
  // empirically <1 step/s vs >40 steps/s for the first cycle).
  //
  // Perturbation cycles are divided evenly among outer cycles so that the
  // total compute budget is approximately unchanged.
  const int n_outer = std::max(1, params.outer_cycles);
  const int max_resets = params.max_outer_resets;  // 0=none, -1=unlimited
  int resets_used = 0;
  // Ceiling division: each outer cycle gets at least 1 ratchet cycle.
  const int ratchet_per = std::max(1,
      (params.ratchet_cycles + n_outer - 1) / n_outer);
  const int drift_per = (params.drift_cycles == 0) ? 0 : std::max(1,
      (params.drift_cycles + n_outer - 1) / n_outer);
  const int nni_perturb_per = (params.nni_perturb_cycles == 0) ? 0 : std::max(1,
      (params.nni_perturb_cycles + n_outer - 1) / n_outer);
  const int prune_reinsert_per = (params.prune_reinsert_cycles == 0) ? 0 :
      std::max(1,
      (params.prune_reinsert_cycles + n_outer - 1) / n_outer);

  int outer = 0;
  while (outer < n_outer) {
    const double score_before_cycle = score_tree(result.tree, ds);
    // Outer-cycle label for verbose output (only shown when n_outer > 1)
    auto outer_label = [&](const char* phase) -> std::string {
      if (n_outer <= 1) return phase;
      return std::string(phase) + " [cycle " + std::to_string(outer + 1) + "]";
    };

    // 3. Sectorial search (XSS + RSS + CSS) if tree is large enough
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

      result.timings.xss_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g [%.0f ms total]\n",
                outer_label("XSS").c_str(),
                score_tree(result.tree, ds), result.timings.xss_ms);
      }

      if (ts::check_interrupt() || check_timeout()) {
        result.interrupted = true;
        result.score = score_tree(result.tree, ds);
        return result;
      }

      // RSS: random sector picks (conflict-guided when pool data available)
      if (params.rss_rounds > 0) {
        sp.split_freq = split_freq;
        for (int rr = 0; rr < params.rss_rounds; ++rr) {
          rss_search(result.tree, ds, sp, cd);

          if (ts::check_interrupt() || check_timeout()) {
            result.interrupted = true;
            result.score = score_tree(result.tree, ds);
            return result;
          }
        }
        result.timings.rss_ms += ph_lap();
        if (verbosity >= 2) {
          Rprintf("  %s score: %.5g [%.0f ms total]\n",
                  outer_label("RSS").c_str(),
                  score_tree(result.tree, ds), result.timings.rss_ms);
        }
      }

      // CSS: sector-restricted TBR on full tree (exact scoring)
      if (params.css_rounds > 0) {
        SectorParams css_sp;
        css_sp.n_partitions = params.css_partitions;
        css_sp.xss_rounds = params.css_rounds;
        css_sp.internal_max_hits = 1;
        css_search(result.tree, ds, css_sp, cd);

        result.timings.css_ms += ph_lap();
        if (verbosity >= 2) {
          Rprintf("  %s score: %.5g [%.0f ms total]\n",
                  outer_label("CSS").c_str(),
                  score_tree(result.tree, ds), result.timings.css_ms);
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
      rp.n_cycles = ratchet_per;
      rp.perturb_prob = params.ratchet_perturb_prob;
      rp.max_hits = params.tbr_max_hits;
      rp.perturb_mode = static_cast<PerturbMode>(params.ratchet_perturb_mode);
      rp.perturb_max_moves = params.ratchet_perturb_max_moves;
      rp.adaptive = params.ratchet_adaptive;
      rp.tabu_size = params.tabu_size;
      ratchet_search(result.tree, ds, rp, cd, check_timeout);
    }
    result.timings.ratchet_ms += ph_lap();
    if (verbosity >= 2) {
      Rprintf("  %s score: %.5g [%.0f ms total]\n",
              outer_label("Ratchet").c_str(),
              score_tree(result.tree, ds), result.timings.ratchet_ms);
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 4a. Post-ratchet sectorial search (T-257)
    // After ratchet perturbation, the tree is in a new basin. A second
    // sectorial pass can exploit local improvements in this basin before
    // TBR polish, approximating TNT's interleaved sectorial pattern.
    if (params.post_ratchet_sectorial && tree_large_enough_for_sectors) {
      SectorParams sp;
      sp.min_sector_size = params.sector_min_size;
      sp.max_sector_size = params.sector_max_size;
      sp.internal_ratchet_cycles = 0;
      sp.internal_max_hits = 1;

      if (params.xss_rounds > 0) {
        sp.n_partitions = params.xss_partitions;
        sp.xss_rounds = params.xss_rounds;
        xss_search(result.tree, ds, sp, cd);
        result.timings.xss_ms += ph_lap();
      }

      if (ts::check_interrupt() || check_timeout()) {
        result.interrupted = true;
        result.score = score_tree(result.tree, ds);
        return result;
      }

      if (params.rss_rounds > 0) {
        sp.split_freq = split_freq;
        for (int rr = 0; rr < params.rss_rounds; ++rr) {
          rss_search(result.tree, ds, sp, cd);
          if (ts::check_interrupt() || check_timeout()) {
            result.interrupted = true;
            result.score = score_tree(result.tree, ds);
            return result;
          }
        }
        result.timings.rss_ms += ph_lap();
      }

      if (params.css_rounds > 0) {
        SectorParams css_sp;
        css_sp.n_partitions = params.css_partitions;
        css_sp.xss_rounds = params.css_rounds;
        css_sp.internal_max_hits = 1;
        css_search(result.tree, ds, css_sp, cd);
        result.timings.css_ms += ph_lap();

        if (ts::check_interrupt() || check_timeout()) {
          result.interrupted = true;
          result.score = score_tree(result.tree, ds);
          return result;
        }
      }

      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g\n",
                outer_label("PostRatch-XSS").c_str(),
                score_tree(result.tree, ds));
      }
    }

    // 4b. NNI perturbation (topology-space escape)
    // `cd` is passed through to nni_perturb_search(), which calls
    // impose_constraint() after perturbation; safe under constraints.
    if (nni_perturb_per > 0) {
      NNIPerturbParams np;
      np.n_cycles = nni_perturb_per;
      np.perturb_fraction = params.nni_perturb_fraction;
      np.max_hits = params.tbr_max_hits;
      np.tabu_size = params.tabu_size;
      nni_perturb_search(result.tree, ds, np, cd, check_timeout);

      result.timings.nni_perturb_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g [%.0f ms total]\n",
                outer_label("NNI-perturb").c_str(),
                score_tree(result.tree, ds), result.timings.nni_perturb_ms);
      }
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 5. Drifting (suboptimal + equal-score exploration)
    if (drift_per > 0) {
      DriftParams dp;
      dp.n_cycles = drift_per;
      dp.afd_limit = params.drift_afd_limit;
      dp.rfd_limit = params.drift_rfd_limit;
      dp.max_hits = params.tbr_max_hits;
      dp.tabu_size = params.tabu_size;
      drift_search(result.tree, ds, dp, cd, check_timeout);

      result.timings.drift_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g [%.0f ms total]\n",
                outer_label("Drift").c_str(),
                score_tree(result.tree, ds), result.timings.drift_ms);
      }
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 5b. SA perturbation (multi-cycle PCSA with best-tree restart).
    // Each cycle: perturb best tree via scheduled SA cooling -> TBR
    // reconverge -> keep if improved (T-207).
    if (params.anneal_cycles > 0) {
      AnnealParams ap;
      ap.t_start = params.anneal_t_start;
      ap.t_end = params.anneal_t_end;
      ap.n_phases = params.anneal_phases;
      ap.moves_per_phase = params.anneal_moves_per_phase;

      double best_sa_score = score_tree(result.tree, ds);
      TreeState best_sa_tree = result.tree;

      for (int cyc = 0; cyc < params.anneal_cycles; ++cyc) {
        if (cyc > 0) result.tree = best_sa_tree;

        anneal_search(result.tree, ds, ap, cd, check_timeout);

        // TBR reconverge after SA perturbation
        {
          TBRParams tp;
          tp.tabu_size = params.tabu_size;
          tbr_search(result.tree, ds, tp, cd, nullptr, nullptr, check_timeout);
        }

        double cyc_score = score_tree(result.tree, ds);
        if (cyc_score < best_sa_score - 1e-10) {
          best_sa_score = cyc_score;
          best_sa_tree = result.tree;
        }

        if (ts::check_interrupt() || check_timeout()) break;
      }

      result.tree = best_sa_tree;
      result.timings.anneal_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g [%.0f ms total]\n",
                outer_label("SA").c_str(),
                score_tree(result.tree, ds), result.timings.anneal_ms);
      }
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 5c. Taxon pruning-reinsertion (T-266).
    // Drop a fraction of leaves, TBR-optimize the backbone, then greedily
    // re-add the dropped taxa and TBR-polish.  Complementary to the ratchet
    // (which perturbs weights) and NNI-perturbation (which perturbs topology).
    if (prune_reinsert_per > 0) {
      PruneReinsertParams prp;
      prp.n_cycles = prune_reinsert_per;
      prp.drop_fraction = params.prune_reinsert_drop;
      prp.selection = static_cast<PruneSelection>(
          params.prune_reinsert_selection);
      prp.tbr_max_hits = params.tbr_max_hits;
      prp.tabu_size = params.tabu_size;
      prp.tbr_max_moves       = params.prune_reinsert_tbr_moves;
      prp.tbr_full_max_moves  = params.prune_reinsert_full_moves;
      prp.nni_full            = (params.prune_reinsert_nni != 0);

      prune_reinsert_search(result.tree, ds, prp, cd, split_freq,
                            check_timeout);

      result.timings.prune_reinsert_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g [%.0f ms total]\n",
                outer_label("PruneRI").c_str(),
                score_tree(result.tree, ds),
                result.timings.prune_reinsert_ms);
      }
    }

    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 6. TBR polish after each outer cycle.
    // Restores local optimality after drift (which accepts suboptimal moves),
    // and seeds the next cycle's XSS from a clean local optimum.
    {
      TBRParams tp;
      tp.tabu_size = params.tabu_size;
      tbr_search(result.tree, ds, tp, cd, nullptr, nullptr, check_timeout);
    }
    result.timings.final_tbr_ms += ph_lap();
    if (verbosity >= 2) {
      Rprintf("  %s score: %.5g [%.0f ms total]\n",
              outer_label("TBR").c_str(),
              score_tree(result.tree, ds), result.timings.final_tbr_ms);
    }

    // Check cancel/timeout after TBR so a stop during TBR is detected
    // here rather than forcing the caller to run MPT enumeration.
    if (ts::check_interrupt() || check_timeout()) {
      result.interrupted = true;
      result.score = score_tree(result.tree, ds);
      return result;
    }

    // 6b. Intra-replicate fusing (T-258).
    // After TBR polish, fuse the current tree against pool donors.
    // The pool is read-only; the fused tree replaces the current replicate
    // tree.  This approximates TNT's within-replicate fusing pattern.
    // tree_fuse() runs TBR internally after each improvement round, so
    // no extra TBR pass is needed here.
    if (params.intra_fuse && pool && pool->size() >= 1) {
      FuseParams fp;
      fp.accept_equal = params.fuse_accept_equal;
      fp.max_rounds = 3;  // brief: just grab low-hanging improvements
      tree_fuse(result.tree, ds, *pool, fp);

      // Rebuild state arrays: tree_fuse may have modified the topology
      // and the internal TBR uses a separate scoring path.
      result.tree.build_postorder();
      result.tree.reset_states(ds);

      result.timings.fuse_ms += ph_lap();
      if (verbosity >= 2) {
        Rprintf("  %s score: %.5g\n",
                outer_label("Intra-fuse").c_str(),
                score_tree(result.tree, ds));
      }
    }

    // If this cycle improved the score and resets are allowed, reset the
    // counter to exploit the new basin.
    const double score_after_cycle = score_tree(result.tree, ds);
    const bool improved = score_after_cycle < score_before_cycle - 1e-8;
    const bool can_reset = (max_resets < 0) ||
                           (max_resets > 0 && resets_used < max_resets);
    if (improved && can_reset) {
      outer = 0;
      ++resets_used;
      if (verbosity >= 2) {
        Rprintf("  Outer cycle improved score (%.5g -> %.5g); resetting"
                " (%d/%s)\n",
                score_before_cycle, score_after_cycle, resets_used,
                max_resets < 0 ? "inf" : std::to_string(max_resets).c_str());
      }
    } else {
      ++outer;
      if (improved && verbosity >= 2) {
        Rprintf("  Outer cycle improved score (%.5g -> %.5g);"
                " reset cap reached (%d)\n",
                score_before_cycle, score_after_cycle, max_resets);
      }
    }
  } // end outer loop

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
  result.n_topologies_at_best = 0;
  result.last_improved_rep = 0;
  result.timed_out = false;
  result.consensus_stable = false;
  result.perturb_stop = false;

  // Perturbation-count stopping rule (T-187).
  int unsuccessful_reps = 0;
  const int perturb_stop_limit = (params.perturb_stop_factor > 0)
      ? ds.n_tips * params.perturb_stop_factor : 0;

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

  // Two-phase timeout (T-202):
  // Main loop exits at budget × (1 - enum_fraction), reserving the
  // remainder for MPT enumeration.
  const double enum_frac = std::max(0.0, std::min(params.enum_time_fraction, 0.5));
  const double main_deadline = params.max_seconds * (1.0 - enum_frac);
  const double full_deadline = params.max_seconds;

  auto check_timeout = [&]() -> bool {
    if (use_timeout && elapsed() >= main_deadline) return true;
    return check_cancel();
  };

  auto check_enum_timeout = [&]() -> bool {
    if (use_timeout && elapsed() >= full_deadline) return true;
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

  // Adaptive search level state.
  // We make a mutable copy of the params to adjust per-replicate cycles.
  // Base values are the originally configured ratchet/drift cycles.
  const int base_ratchet_cycles = params.ratchet_cycles;
  const int base_drift_cycles = params.drift_cycles;
  DrivenParams adaptive_params = params;  // mutable copy for adaptive level

  // Adaptive starting-tree strategy (T-190).
  StrategyTracker strategy_tracker;
  // Seed a dedicated RNG for Thompson sampling from R's RNG.
  GetRNGstate();
  std::mt19937 bandit_rng(static_cast<unsigned>(unif_rand() * 4294967295.0));
  PutRNGstate();

  // Cross-replicate consensus constraint tightening.
  // When enabled and no user constraint is supplied, the strict consensus
  // of pool trees is enforced as topological constraints for subsequent
  // replicates. Cleared whenever a new best score is found.
  bool use_auto_constraint = params.consensus_constrain && (!cd || !cd->active);
  ConstraintData auto_cd;  // built from pool consensus; reused across reps
  double auto_cd_best_score = 1e18;  // score when auto_cd was last built

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

    // Adaptive level: adjust ratchet/drift cycles based on hit rate.
    // Uses a sliding window of recent replicates to compute hit rate.
    if (params.adaptive_level && rep > 0 && pool.size() > 0) {
      // Simple hit rate: hits_to_best / replicates_completed
      double hit_rate = (result.replicates_completed > 0)
          ? static_cast<double>(pool.hits_to_best()) /
            result.replicates_completed
          : 0.0;

      // Scale factor: high hit rate → reduce effort, low → increase
      double scale;
      if (hit_rate > 0.7) {
        scale = 0.5;     // easy landscape: halve effort
      } else if (hit_rate > 0.4) {
        scale = 0.75;    // moderate: reduce slightly
      } else if (hit_rate < 0.15) {
        scale = 1.5;     // hard landscape: increase effort
      } else {
        scale = 1.0;     // use base values
      }

      adaptive_params.ratchet_cycles = std::max(
          1, static_cast<int>(base_ratchet_cycles * scale));
      adaptive_params.drift_cycles = std::max(
          0, static_cast<int>(base_drift_cycles * scale));

      if (params.verbosity >= 2) {
        Rprintf("  Adaptive level: hit_rate=%.2f, scale=%.2f "
                "(ratchet=%d, drift=%d)\n",
                hit_rate, scale,
                adaptive_params.ratchet_cycles,
                adaptive_params.drift_cycles);
      }
    }

    // Adaptive ratchet taper (T-182): reduce perturbation probability as
    // the pool stabilizes.  High hit rate = stable pool = gentler perturbation
    // for finer local exploration.  Resets to base when a new best is found.
    if (params.ratchet_taper && rep > 0 && pool.size() > 0) {
      double stability = (result.replicates_completed > 0)
          ? static_cast<double>(pool.hits_to_best()) /
            result.replicates_completed
          : 0.0;
      double taper_factor = std::max(
          params.ratchet_taper_floor,
          1.0 - params.ratchet_taper_strength * stability);
      adaptive_params.ratchet_perturb_prob =
          params.ratchet_perturb_prob * taper_factor;

      if (params.verbosity >= 2 && !has_callback) {
        Rprintf("  Ratchet taper: stability=%.2f, prob=%.3f (base=%.3f)\n",
                stability, adaptive_params.ratchet_perturb_prob,
                params.ratchet_perturb_prob);
      }
    }

    // Use adaptive_params when any per-replicate adaptation is active
    const DrivenParams& rep_params =
        (params.adaptive_level || params.ratchet_taper)
        ? adaptive_params : params;

    // Conflict-guided sector selection: compute pool split frequencies
    // for RSS weighting (only useful when pool has ≥2 best-score trees).
    const SplitFrequencyTable* sft_ptr = nullptr;
    SplitFrequencyTable sft;
    if (pool.size() >= 2) {
      sft = pool.compute_split_frequencies();
      if (sft.n_trees >= 2) sft_ptr = &sft;
    }

    // Consensus constraint tightening: build/update auto-constraints
    ConstraintData* rep_cd = cd;  // default: user-supplied constraint
    if (use_auto_constraint &&
        result.replicates_completed >= params.consensus_constrain_min_reps &&
        pool.size() >= 3) {
      // Rebuild if pool changed score (meaning old constraints may be wrong)
      if (pool.best_score() < auto_cd_best_score || !auto_cd.active) {
        auto_cd.active = false;  // clear old constraints
        auto_cd_best_score = pool.best_score();

        int n_unan = 0, wps = 0;
        auto bits = pool.extract_consensus_splits(n_unan, wps);
        if (n_unan > 0) {
          auto_cd = build_constraint_from_bitsets(
              bits.data(), n_unan, wps, ds.n_tips);
          if (params.verbosity >= 2 && !has_callback) {
            Rprintf("  Auto-constraint: %d consensus splits locked\n",
                    n_unan);
          }
        }
      }
      if (auto_cd.active) {
        rep_cd = &auto_cd;
      }
    }

    // Select starting-tree strategy for this replicate.
    StartStrategy rep_strategy = StartStrategy::WAGNER_RANDOM;
    if (start_ptr) {
      // User-supplied starting tree for rep 0 — strategy is moot
    } else if (params.adaptive_start) {
      rep_strategy = strategy_tracker.select(bandit_rng);
    } else if (params.wagner_bias != 0) {
      // Legacy fixed-bias mode
      rep_strategy = static_cast<StartStrategy>(params.wagner_bias);
    }

    if (params.verbosity >= 2 && params.adaptive_start && !has_callback) {
      Rprintf("  Strategy: %s\n", strategy_name(rep_strategy));
    }

    // Run the single-replicate pipeline
    ReplicateResult rep_result = run_single_replicate(
        ds, rep_params, rep_cd, check_timeout, params.verbosity, start_ptr,
        sft_ptr, rep_strategy, &pool);

    result.timings += rep_result.timings;

    // Compute collapsed flags for collapsed-topology pool dedup.
    // Trees that differ only in zero-length resolutions are treated
    // as duplicates, improving pool diversity (Goloboff & Farris 2001).
    std::vector<uint8_t> rep_collapsed;
    compute_collapsed_flags(rep_result.tree, ds, rep_collapsed);

    if (rep_result.interrupted) {
      if (rep_result.score < 1e18) {
        pool.add_collapsed(rep_result.tree, rep_result.score, rep_collapsed);
      }
      result.timed_out = true;
      goto finish;
    }

    // Add to pool with collapsed-topology dedup
    double prev_best = pool.best_score();
    pool.add_collapsed(rep_result.tree, rep_result.score, rep_collapsed);
    bool score_improved = pool.best_score() < prev_best;
    if (score_improved) {
      result.last_improved_rep = rep1;
      unsuccessful_reps = 0;
    } else {
      ++unsuccessful_reps;
    }

    // Update strategy bandit (T-190)
    if (params.adaptive_start) {
      bool hit_best = (rep_result.score <= pool.best_score());
      strategy_tracker.update(rep_strategy, hit_best);
      if (score_improved) {
        strategy_tracker.decay(0.5);
      }
    }

    ++result.replicates_completed;
    result.replicate_scores.push_back(rep_result.score);

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

      // Check and repair constraint violations on fused tree.
      // impose_constraint() is heuristic — verify the repair succeeded
      // and discard the tree if it still violates.
      bool fused_ok = true;
      if (cd && cd->active) {
        map_constraint_nodes(fused, *cd);
        bool viol = false;
        for (int _s = 0; _s < cd->n_splits; ++_s) {
          if (cd->constraint_node[_s] < 0) { viol = true; break; }
        }
        if (viol) {
          impose_constraint(fused, *cd);
          fused.build_postorder();
          fused.reset_states(ds);
          fused_score = score_tree(fused, ds);
          // Verify repair succeeded
          map_constraint_nodes(fused, *cd);
          for (int _s = 0; _s < cd->n_splits; ++_s) {
            if (cd->constraint_node[_s] < 0) { fused_ok = false; break; }
          }
        }
      }
      if (fused_ok) {
        std::vector<uint8_t> fused_collapsed;
        compute_collapsed_flags(fused, ds, fused_collapsed);
        pool.add_collapsed(fused, fused_score, fused_collapsed);
      }

      if (fused_ok && fused_score < best_before) {
        pool.set_hits_to_best(0);
        result.last_improved_rep = rep1;
        unsuccessful_reps = 0;  // fuse found a better score; reset perturb-stop counter
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

    // Consensus stability check
    if (params.consensus_stable_reps > 0 && pool.size() >= 2) {
      int unchanged = pool.update_consensus_stability();
      if (unchanged >= params.consensus_stable_reps) {
        if (params.verbosity >= 1) {
          if (!has_callback) {
            Rprintf("Consensus stable for %d replicates (score %.5g, "
                    "pool %d trees)\n",
                    unchanged, pool.best_score(), pool.size());
          }
        }
        result.consensus_stable = true;
        break;
      }
    }

    // Convergence check (hit count)
    if (pool.hits_to_best() >= params.target_hits) {
      if (params.verbosity >= 1) {
        if (!has_callback) {
          Rprintf("Converged: %d hits to best score %.5g\n",
                  pool.hits_to_best(), pool.best_score());
        }
      }
      break;
    }

    // Perturbation-count stopping rule (T-187)
    if (perturb_stop_limit > 0 && unsuccessful_reps >= perturb_stop_limit) {
      if (params.verbosity >= 1) {
        if (!has_callback) {
          Rprintf("Stopped: %d consecutive unsuccessful replicates "
                  "(limit %d = %d tips x %d)\n",
                  unsuccessful_reps, perturb_stop_limit,
                  ds.n_tips, params.perturb_stop_factor);
        }
      }
      result.perturb_stop = true;
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
  if (pool.size() > 0 && pool.size() < pool.max_size) {
    TBRParams tp;
    tp.accept_equal = true;
    tp.tabu_size = params.tabu_size > 0 ? params.tabu_size : 100;

    // Snapshot current pool entries as seeds (new trees discovered during
    // enumeration of seed i become additional seeds for later iterations).
    int seed_idx = 0;
    while (seed_idx < pool.size() && pool.size() < pool.max_size) {
      if (check_enum_timeout()) break;
      TreeState enum_tree = pool.all()[seed_idx].tree;
      // Budget remaining capacity across remaining seeds
      tp.max_hits = std::max(10, (pool.max_size - pool.size()) * 2);
      tbr_search(enum_tree, ds, tp, cd, nullptr, &pool, check_enum_timeout);
      ++seed_idx;
    }
    if (params.verbosity >= 2) {
      Rprintf("MPT enumeration: %d trees in pool (%.1f s)\n",
              pool.size(), elapsed());
    }
  }

  // result.hits_to_best already set before MPT enumeration
  result.pool_size = pool.size();
  result.n_topologies_at_best = pool.count_at_best();

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

  // Populate per-strategy diagnostics
  if (params.adaptive_start) {
    for (int i = 0; i < N_STRAT; ++i) {
      auto s = static_cast<StartStrategy>(i);
      result.strategy_attempts[i] = strategy_tracker.attempts(s);
      result.strategy_successes[i] = strategy_tracker.successes(s);
    }
  }

  return result;
}

} // namespace ts
