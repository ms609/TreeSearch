#include "ts_parallel.h"
#include "ts_rng.h"
#include "ts_fitch.h"
#include "ts_fuse.h"

#include <R.h>
#include <Rmath.h>

#include <thread>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <functional>

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

  bool fuse_ok = true;
  if (cd && cd->active) {
    fuse_ok = !violates_constraint_posthoc(fused, *cd);
  }
  if (fuse_ok) pool_.add(fused, fused_score);

  if (fused_score < best_before) {
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

    // Run the replicate pipeline (verbosity=0 for parallel)
    ReplicateResult rep_result = run_single_replicate(
        ds_local, *ctx.params, cd_ptr, check_timeout_noop, 0);

    if (ctx.stop_flag->load(std::memory_order_relaxed)) break;

    // Add to shared pool
    ctx.shared_pool->add(rep_result.tree, rep_result.score);
    ctx.replicates_done->fetch_add(1, std::memory_order_relaxed);

    // Check convergence
    if (ctx.shared_pool->hits_to_best() >= ctx.params->target_hits) {
      ctx.stop_flag->store(true, std::memory_order_relaxed);
      break;
    }

    // Periodic fuse
    int done = ctx.replicates_done->load(std::memory_order_relaxed);
    if (done > 0 && done % ctx.params->fuse_interval == 0
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
  result.timed_out = false;

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

  // Timeout setup
  bool use_timeout = params.max_seconds > 0.0;
  auto start_time = std::chrono::steady_clock::now();

  // Spawn worker threads
  std::vector<std::thread> workers;
  workers.reserve(n_threads);
  for (int t = 0; t < n_threads; ++t) {
    workers.emplace_back(worker_thread, ctx);
  }

  // Main thread: poll for interrupt and timeout
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

    // Check timeout
    if (use_timeout) {
      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - start_time).count();
      if (elapsed >= params.max_seconds) {
        stop_flag.store(true, std::memory_order_relaxed);
        result.timed_out = true;
        break;
      }
    }

    // Progress reporting
    if (params.verbosity >= 1) {
      auto st = shared_pool.status();
      int done = replicates_done.load(std::memory_order_relaxed);
      Rprintf("[%d threads] Replicates: %d/%d | Best: %.1f | Pool: %d | Hits: %d\n",
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

  // Extract results
  shared_pool.extract_into(pool_out);

  result.replicates_completed = replicates_done.load();
  result.hits_to_best = pool_out.hits_to_best();
  result.pool_size = pool_out.size();
  if (pool_out.size() > 0) {
    result.best_score = pool_out.best_score();
  } else {
    result.best_score = -1.0;
  }

  if (params.verbosity >= 1) {
    if (result.timed_out) {
      Rprintf("Timeout reached (%.1f s)\n", params.max_seconds);
    } else if (result.hits_to_best >= params.target_hits) {
      Rprintf("Converged: %d hits to best score %.1f (%d replicates)\n",
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
    ConstraintData* cd)
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
          info_amounts_r, info_max_steps, cd);
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
