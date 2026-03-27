#include "ts_parallel.h"
#include "ts_collapsed.h"
#include "ts_constraint.h"
#include "ts_rng.h"
#include "ts_fitch.h"
#include "ts_fuse.h"
#include "ts_tbr.h"

#include <R.h>
#include <Rmath.h>

#include <thread>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <string>

namespace ts {

// --- ThreadSafePool ---

void ThreadSafePool::fuse_round(DataSet& ds, const DrivenParams& params,
                                ConstraintData* cd) {
  std::lock_guard<std::mutex> lock(mu_);
  if (pool_.size() < 2) return;

  int hits_before = pool_.hits_to_best();
  double best_before = pool_.best_score();

  TreeState fused = pool_.best().tree;
  FuseParams fp;
  fp.accept_equal = params.fuse_accept_equal;
  fp.max_rounds = 10;
  tree_fuse(fused, ds, pool_, fp);

  double fused_score = score_tree(fused, ds);

  bool fused_ok = true;
  if (cd && cd->active &&
      violates_constraint_posthoc(fused, *cd)) {
    impose_constraint(fused, *cd);
    fused.build_postorder();
    fused.reset_states(ds);
    fused_score = score_tree(fused, ds);
    // Verify repair succeeded — impose_constraint is heuristic
    map_constraint_nodes(fused, *cd);
    for (int s = 0; s < cd->n_splits; ++s) {
      if (cd->constraint_node[s] < 0) { fused_ok = false; break; }
    }
  }
  if (fused_ok) {
    std::vector<uint8_t> fused_collapsed;
    compute_collapsed_flags(fused, ds, fused_collapsed);
    pool_.add_collapsed(fused, fused_score, fused_collapsed);
  }

  if (fused_ok && fused_score < best_before) {
    pool_.set_hits_to_best(0);
  } else {
    pool_.set_hits_to_best(hits_before);
  }
}

void ThreadSafePool::extract_into(TreePool& out) {
  // No lock needed — called after all threads joined
  const auto& entries = pool_.all();
  for (const auto& e : entries) {
    out.add(e.tree, e.score);
  }
  // Propagate the actual independent-hit count.  The add() calls above
  // only count one hit per distinct topology; the internal pool tracks
  // the true number of independent replicate hits (including duplicates
  // that matched best_score but were deduped).
  out.set_hits_to_best(pool_.hits_to_best());
}

// --- Worker thread function ---

namespace {

struct WorkerContext {
  const DataSet* ds_prototype;
  const DrivenParams* params;
  const ConstraintData* cd_prototype;

  ThreadSafePool* shared_pool;
  std::atomic<bool>* stop_flag;
  std::atomic<int>* replicates_claimed;
  std::atomic<int>* replicates_done;

  // Pre-generated seeds (one per replicate)
  const std::vector<unsigned>* seeds;

  // Pre-computed strategy sequence for round-robin (T-190)
  const std::vector<StartStrategy>* strategies;

  // Per-thread timing accumulator (index = thread_id)
  PhaseTimings* thread_timings;
  int thread_id;

  // Per-thread score accumulator (index = thread_id)
  std::vector<double>* thread_scores;
};

void worker_thread(WorkerContext ctx) {
  // Make thread-local copies of mutable data
  DataSet ds_local = *ctx.ds_prototype;

  ConstraintData cd_local;
  ConstraintData* cd_ptr = nullptr;
  if (ctx.cd_prototype && ctx.cd_prototype->active) {
    cd_local = *ctx.cd_prototype;
    cd_ptr = &cd_local;
  }

  // Set up thread-local RNG
  std::mt19937 local_rng(0);  // will be re-seeded per replicate
  ts::thread_rng = &local_rng;
  ts::thread_stop_flag = ctx.stop_flag;

  auto check_timeout_noop = []() -> bool { return false; };
  // Timeout is handled by the main thread setting stop_flag

  while (true) {
    int rep = ctx.replicates_claimed->fetch_add(1, std::memory_order_relaxed);
    if (rep >= ctx.params->max_replicates) break;
    if (ctx.stop_flag->load(std::memory_order_relaxed)) break;

    // Seed RNG for this replicate
    local_rng.seed((*ctx.seeds)[rep]);

    // Use starting tree for replicate 0 if provided
    TreeState* start_ptr = nullptr;
    TreeState start_tree;
    if (rep == 0 && ctx.params->start_n_edge > 0 &&
        static_cast<int>(ctx.params->start_edge.size()) >=
            2 * ctx.params->start_n_edge) {
      const int* edge_parent = ctx.params->start_edge.data();
      const int* edge_child =
          ctx.params->start_edge.data() + ctx.params->start_n_edge;
      start_tree.init_from_edge(edge_parent, edge_child,
                                ctx.params->start_n_edge, ds_local);
      start_ptr = &start_tree;
    }

    // Strategy for this replicate (round-robin when adaptive, else default)
    StartStrategy rep_strat = StartStrategy::WAGNER_RANDOM;
    if (ctx.strategies && rep < static_cast<int>(ctx.strategies->size())) {
      rep_strat = (*ctx.strategies)[rep];
    }

    // Run the replicate pipeline (verbosity=0 for parallel)
    // pool=nullptr: intra-fuse disabled in parallel mode (between-replicate
    // fusing via ThreadSafePool::fuse_round() is already active)
    ReplicateResult rep_result = run_single_replicate(
        ds_local, *ctx.params, cd_ptr, check_timeout_noop, 0, start_ptr,
        nullptr, rep_strat, nullptr);

    if (ctx.stop_flag->load(std::memory_order_relaxed)) break;

    // Accumulate phase timings for this thread
    ctx.thread_timings[ctx.thread_id] += rep_result.timings;

    // Add to shared pool with collapsed-topology dedup
    std::vector<uint8_t> rep_collapsed;
    compute_collapsed_flags(rep_result.tree, ds_local, rep_collapsed);
    ctx.shared_pool->add_collapsed(rep_result.tree, rep_result.score,
                                   rep_collapsed);

    // Record per-replicate score for Chao1 coverage estimation
    ctx.thread_scores[ctx.thread_id].push_back(rep_result.score);

    ctx.replicates_done->fetch_add(1, std::memory_order_relaxed);

    // Check convergence
    if (ctx.shared_pool->hits_to_best() >= ctx.params->target_hits) {
      ctx.stop_flag->store(true, std::memory_order_relaxed);
      break;
    }

    // Periodic fuse
    int done = ctx.replicates_done->load(std::memory_order_relaxed);
    if (ctx.params->fuse_interval > 0 &&
        done > 0 && done % ctx.params->fuse_interval == 0
        && ctx.shared_pool->size() >= 2) {
      ctx.shared_pool->fuse_round(ds_local, *ctx.params, cd_ptr);
    }
  }

  // Clean up thread-local pointers
  ts::thread_rng = nullptr;
  ts::thread_stop_flag = nullptr;
}

} // anonymous namespace

// --- Parallel driven search ---

DrivenResult parallel_driven_search(
    TreePool& pool_out,
    const DataSet& ds_prototype,
    const DrivenParams& params,
    const ConstraintData* cd,
    int n_threads)
{
  DrivenResult result;
  result.best_score = 1e18;
  result.replicates_completed = 0;
  result.hits_to_best = 0;
  result.pool_size = 0;
  result.n_topologies_at_best = 0;
  result.last_improved_rep = 0;  // not tracked in parallel (replicates out of order)
  result.timed_out = false;
  result.consensus_stable = false;
  result.perturb_stop = false;

  if (params.max_replicates <= 0) {
    result.best_score = -1.0;
    return result;
  }

  // Auto-detect thread count
  if (n_threads <= 0) {
    n_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads <= 1) n_threads = 2;  // at least 2 if auto
    n_threads = std::min(n_threads, params.max_replicates);
  }
  n_threads = std::max(1, std::min(n_threads, params.max_replicates));

  // Pre-generate RNG seeds from R (must be done on main thread)
  std::vector<unsigned> seeds(params.max_replicates);
  GetRNGstate();
  for (int i = 0; i < params.max_replicates; ++i) {
    seeds[i] = static_cast<unsigned>(unif_rand() * 4294967295.0);
  }
  PutRNGstate();

  // Set up shared state
  ThreadSafePool shared_pool(params.pool_max_size, params.pool_suboptimal);
  std::atomic<bool> stop_flag(false);
  std::atomic<int> replicates_claimed(0);
  std::atomic<int> replicates_done(0);

  // Perturbation-count stopping rule (T-187, parallel path).
  const int perturb_stop_limit = (params.perturb_stop_factor > 0)
      ? ds_prototype.n_tips * params.perturb_stop_factor : 0;
  double last_known_best = 1e18;
  int reps_at_last_improvement = 0;

  // Prepare worker context
  WorkerContext ctx;
  ctx.ds_prototype = &ds_prototype;
  ctx.params = &params;
  ctx.cd_prototype = cd;
  ctx.shared_pool = &shared_pool;
  ctx.stop_flag = &stop_flag;
  ctx.replicates_claimed = &replicates_claimed;
  ctx.replicates_done = &replicates_done;
  ctx.seeds = &seeds;

  // Pre-compute round-robin strategy sequence for adaptive start (T-190)
  std::vector<StartStrategy> strategies;
  if (params.adaptive_start) {
    strategies = StrategyTracker::round_robin(params.max_replicates);
  }
  ctx.strategies = params.adaptive_start ? &strategies : nullptr;

  // Two-phase timeout (T-202): main loop exits early, reserving time for
  // MPT enumeration.
  bool use_timeout = params.max_seconds > 0.0;
  auto start_time = std::chrono::steady_clock::now();
  const double enum_frac = std::max(0.0,
                                     std::min(params.enum_time_fraction, 0.5));
  const double main_deadline = params.max_seconds * (1.0 - enum_frac);
  const double full_deadline = params.max_seconds;

  // Cancel file: read path from environment variable (set by Shiny app).
  // If the file exists, the search should stop.
  std::string cancel_path;
  {
    const char* cancel_env = std::getenv("TREESEARCH_CANCEL_FILE");
    if (cancel_env && cancel_env[0] != '\0') cancel_path = cancel_env;
  }

  // Per-thread timing and score accumulators
  std::vector<PhaseTimings> thread_timings(n_threads);
  std::vector<std::vector<double>> thread_scores(n_threads);

  // Spawn worker threads
  std::vector<std::thread> workers;
  workers.reserve(n_threads);
  for (int t = 0; t < n_threads; ++t) {
    ctx.thread_timings = thread_timings.data();
    ctx.thread_scores = thread_scores.data();
    ctx.thread_id = t;
    workers.emplace_back(worker_thread, ctx);
  }

  // Main thread: poll for interrupt and timeout
  int last_stab_done = 0;  // replicates_done at last consensus check
  while (true) {
    // Sleep briefly to avoid spinning
    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    // Check if all workers are done
    if (replicates_done.load(std::memory_order_relaxed) >= params.max_replicates
        || stop_flag.load(std::memory_order_relaxed)) {
      break;
    }

    // Check if all replicates have been claimed (workers finishing up)
    if (replicates_claimed.load(std::memory_order_relaxed) >= params.max_replicates) {
      // Workers are still running their last replicate — keep waiting
      // but less aggressively
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
      if (stop_flag.load(std::memory_order_relaxed)) break;
      // Don't spin forever — check if replicates_done caught up
      if (replicates_done.load(std::memory_order_relaxed)
          >= replicates_claimed.load(std::memory_order_relaxed) - n_threads) {
        // Allow some in-flight replicates
      }
    }

    // Check user interrupt (R API — main thread only)
    try {
      R_CheckUserInterrupt();
    } catch (...) {
      stop_flag.store(true, std::memory_order_relaxed);
      break;
    }

    // Check timeout (main loop deadline, reserving time for MPT enum)
    if (use_timeout) {
      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - start_time).count();
      if (elapsed >= main_deadline) {
        stop_flag.store(true, std::memory_order_relaxed);
        result.timed_out = true;
        break;
      }
    }

    // Check cancel file (set by Shiny app or external process)
    if (!cancel_path.empty()) {
      FILE* cf = std::fopen(cancel_path.c_str(), "r");
      if (cf) {
        std::fclose(cf);
        stop_flag.store(true, std::memory_order_relaxed);
        result.timed_out = true;
        break;
      }
    }

    // Consensus stability check (parallel path).
    // Only check when new replicates have completed; otherwise the
    // unchanged counter increments on idle polls (every 200 ms) and
    // can trigger premature termination with slow replicates.
    if (params.consensus_stable_reps > 0) {
      int done_now = replicates_done.load(std::memory_order_relaxed);
      auto st = shared_pool.status();
      if (st.pool_size >= 2 && done_now > last_stab_done) {
        last_stab_done = done_now;
        int unchanged = shared_pool.update_consensus_stability();
        if (unchanged >= params.consensus_stable_reps) {
          stop_flag.store(true, std::memory_order_relaxed);
          result.consensus_stable = true;
          if (params.verbosity >= 1) {
            Rprintf("Consensus stable for %d replicates (score %.5g, "
                    "pool %d trees)\n",
                    unchanged, st.best_score, st.pool_size);
          }
          break;
        }
      }
    }


    // Perturbation-count stopping rule (T-187, parallel path)
    if (perturb_stop_limit > 0) {
      int done = replicates_done.load(std::memory_order_relaxed);
      double cur_best = shared_pool.best_score();
      if (cur_best < last_known_best) {
        last_known_best = cur_best;
        reps_at_last_improvement = done;
      }
      if (done - reps_at_last_improvement >= perturb_stop_limit) {
        stop_flag.store(true, std::memory_order_relaxed);
        result.perturb_stop = true;
        if (params.verbosity >= 1) {
          Rprintf("Stopped: %d consecutive unsuccessful replicates "
                  "(limit %d = %d tips x %d)\n",
                  done - reps_at_last_improvement, perturb_stop_limit,
                  ds_prototype.n_tips, params.perturb_stop_factor);
        }
        break;
      }
    }
    // Progress reporting
    if (params.verbosity >= 1) {
      auto st = shared_pool.status();
      int done = replicates_done.load(std::memory_order_relaxed);
      Rprintf("[%d threads] Replicates: %d/%d | Best: %.5g | Pool: %d | Hits: %d\n",
              n_threads, done, params.max_replicates,
              st.best_score, st.pool_size, st.hits_to_best);
    }
  }

  // Signal stop to all workers
  stop_flag.store(true, std::memory_order_relaxed);

  // Join all worker threads
  for (auto& w : workers) {
    if (w.joinable()) w.join();
  }

  // Sum per-thread timings; merge per-thread replicate scores
  for (int t = 0; t < n_threads; ++t) {
    result.timings += thread_timings[t];
    for (double s : thread_scores[t]) {
      result.replicate_scores.push_back(s);
    }
  }

  // Extract results
  shared_pool.extract_into(pool_out);

  // Capture hits_to_best BEFORE MPT enumeration — only main-loop
  // replicates are independent; MPT enumeration should not inflate count.
  result.replicates_completed = replicates_done.load();
  result.hits_to_best = pool_out.hits_to_best();

  // MPT enumeration: TBR plateau walk from each pool tree (serial, main
  // thread).  T-202: always runs (even after timeout), subject to the
  // reserved enum time budget.
  auto elapsed_now = [&]() -> double {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_time).count();
  };
  auto check_enum_timeout = [&]() -> bool {
    if (use_timeout && elapsed_now() >= full_deadline) return true;
    if (!cancel_path.empty()) {
      FILE* cf = std::fopen(cancel_path.c_str(), "r");
      if (cf) { std::fclose(cf); return true; }
    }
    return false;
  };

  if (pool_out.size() > 0 && pool_out.size() < pool_out.max_size) {
    TBRParams tp;
    tp.accept_equal = true;
    tp.tabu_size = 100;
    std::unique_ptr<ConstraintData> cd_local;
    ConstraintData* cd_ptr = nullptr;
    if (cd && cd->active) {
      cd_local = std::make_unique<ConstraintData>(*cd);
      cd_ptr = cd_local.get();
    }
    int seed_idx = 0;
    while (seed_idx < pool_out.size() && pool_out.size() < pool_out.max_size) {
      if (check_enum_timeout()) break;
      TreeState enum_tree = pool_out.all()[seed_idx].tree;
      tp.max_hits = std::max(10, (pool_out.max_size - pool_out.size()) * 2);
      tbr_search(enum_tree, ds_prototype, tp, cd_ptr, nullptr, &pool_out,
                 check_enum_timeout);
      ++seed_idx;
    }
  }

  // result.replicates_completed and result.hits_to_best already set
  // before MPT enumeration (above).
  result.pool_size = pool_out.size();
  result.n_topologies_at_best = pool_out.count_at_best();
  if (pool_out.size() > 0) {
    result.best_score = pool_out.best_score();
  } else {
    result.best_score = -1.0;
  }

  if (params.verbosity >= 1) {
    if (result.timed_out) {
      Rprintf("Timeout reached (%.1f s)\n", params.max_seconds);
    } else if (result.hits_to_best >= params.target_hits) {
      Rprintf("Converged: %d hits to best score %.5g (%d replicates)\n",
              result.hits_to_best, result.best_score,
              result.replicates_completed);
    }
  }

  return result;
}

// --- Parallel resampling ---

std::vector<ResampleResult> parallel_resample(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* original_weights,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const ResampleParams& params,
    int n_replicates,
    int n_threads,
    const double* info_amounts_r,
    int info_max_steps,
    ConstraintData* cd,
    bool xpiwe,
    double xpiwe_r,
    double xpiwe_max_f,
    const int* obs_count_r)
{
  if (n_threads <= 0) {
    n_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads <= 1) n_threads = 2;
    n_threads = std::min(n_threads, n_replicates);
  }
  n_threads = std::max(1, std::min(n_threads, n_replicates));

  // Pre-generate seeds (main thread)
  std::vector<unsigned> seeds(n_replicates);
  GetRNGstate();
  for (int i = 0; i < n_replicates; ++i) {
    seeds[i] = static_cast<unsigned>(unif_rand() * 4294967295.0);
  }
  PutRNGstate();

  std::vector<ResampleResult> results(n_replicates);
  std::atomic<int> next_rep(0);
  std::atomic<bool> stop_flag(false);

  auto worker = [&](int /*thread_id*/) {
    // Set up thread-local RNG
    std::mt19937 local_rng(0);
    ts::thread_rng = &local_rng;
    ts::thread_stop_flag = &stop_flag;

    while (true) {
      int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
      if (rep >= n_replicates) break;
      if (stop_flag.load(std::memory_order_relaxed)) break;

      local_rng.seed(seeds[rep]);

      results[rep] = resample_search(
          contrast_r, n_tokens, n_states,
          tip_data_r, n_tips, n_patterns,
          original_weights, levels_r, min_steps_r,
          concavity, params,
          info_amounts_r, info_max_steps, cd,
          xpiwe, xpiwe_r, xpiwe_max_f, obs_count_r);
    }

    ts::thread_rng = nullptr;
    ts::thread_stop_flag = nullptr;
  };

  // Spawn workers
  std::vector<std::thread> workers;
  workers.reserve(n_threads);
  for (int t = 0; t < n_threads; ++t) {
    workers.emplace_back(worker, t);
  }

  // Main thread polls for interrupt
  while (true) {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    if (next_rep.load(std::memory_order_relaxed) >= n_replicates) break;
    if (stop_flag.load(std::memory_order_relaxed)) break;

    try {
      R_CheckUserInterrupt();
    } catch (...) {
      stop_flag.store(true, std::memory_order_relaxed);
      break;
    }
  }

  stop_flag.store(true, std::memory_order_relaxed);
  for (auto& w : workers) {
    if (w.joinable()) w.join();
  }

  return results;
}

} // namespace ts
