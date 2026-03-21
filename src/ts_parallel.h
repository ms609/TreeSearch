#ifndef TS_PARALLEL_H
#define TS_PARALLEL_H

// Thread-parallel driven search and resampling.
//
// Inter-replicate parallelism: each search replicate runs on its own
// thread with a private DataSet + ConstraintData copy. A shared
// ThreadSafePool collects results.
//
// Design:
// - Main thread: pre-generates RNG seeds, spawns workers, polls for
//   R interrupt and timeout, joins workers, returns results.
// - Worker threads: claim replicates via atomic counter, run pipeline,
//   add to shared pool.
// - No R API calls from worker threads (RNG, interrupt, output all
//   handled via ts_rng.h thread-local indirection).

#include "ts_driven.h"
#include "ts_pool.h"
#include "ts_data.h"
#include "ts_constraint.h"
#include "ts_resample.h"

#include <mutex>
#include <atomic>
#include <vector>

namespace ts {

// --- Thread-safe pool ---

class ThreadSafePool {
public:
  ThreadSafePool(int max_size, double suboptimal)
    : pool_(max_size, suboptimal) {}

  bool add(const TreeState& tree, double score) {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.add(tree, score);
  }

  PoolEntry best() const {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.best();
  }

  double best_score() const {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.best_score();
  }

  int hits_to_best() const {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.hits_to_best();
  }

  void set_hits_to_best(int n) {
    std::lock_guard<std::mutex> lock(mu_);
    pool_.set_hits_to_best(n);
  }

  int size() const {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.size();
  }

  // Get status snapshot (score + hits) atomically
  struct PoolStatus {
    double best_score;
    int hits_to_best;
    int pool_size;
  };

  PoolStatus status() const {
    std::lock_guard<std::mutex> lock(mu_);
    PoolStatus s;
    s.best_score = pool_.best_score();
    s.hits_to_best = pool_.hits_to_best();
    s.pool_size = pool_.size();
    return s;
  }

  // Perform fuse round under lock. Caller provides thread-local ds/cd.
  void fuse_round(DataSet& ds, const DrivenParams& params,
                  ConstraintData* cd);

  // Move contents into output pool (called after all threads joined)
  void extract_into(TreePool& out);

  // Update consensus stability and return consecutive-unchanged count.
  int update_consensus_stability() {
    std::lock_guard<std::mutex> lock(mu_);
    return pool_.update_consensus_stability();
  }

  // Direct access to underlying pool (only safe when no threads running)
  const TreePool& pool() const { return pool_; }

private:
  TreePool pool_;
  mutable std::mutex mu_;
};

// --- Parallel driven search ---

DrivenResult parallel_driven_search(
    TreePool& pool_out,
    const DataSet& ds_prototype,
    const DrivenParams& params,
    const ConstraintData* cd,
    int n_threads);

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
    const double* info_amounts_r = nullptr,
    int info_max_steps = 0,
    ConstraintData* cd = nullptr);

} // namespace ts

#endif // TS_PARALLEL_H
