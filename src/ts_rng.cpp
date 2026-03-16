#include "ts_rng.h"

#include <R.h>
#include <Rmath.h>

namespace ts {

thread_local std::mt19937* thread_rng = nullptr;
thread_local const std::atomic<bool>* thread_stop_flag = nullptr;

std::mt19937 make_rng() {
  if (thread_rng) {
    // Parallel mode: seed from thread-local RNG
    return std::mt19937((*thread_rng)());
  }
  // Serial mode: seed from R's RNG
  GetRNGstate();
  unsigned seed = static_cast<unsigned>(unif_rand() * 4294967295.0);
  PutRNGstate();
  return std::mt19937(seed);
}

bool check_interrupt() {
  if (thread_stop_flag) {
    return thread_stop_flag->load(std::memory_order_relaxed);
  }
  // Serial mode: R_CheckUserInterrupt() longjmps on interrupt
  R_CheckUserInterrupt();
  return false;
}

double thread_safe_unif() {
  if (thread_rng) {
    // Map mt19937 output to [0, 1)
    return std::uniform_real_distribution<double>(0.0, 1.0)(*thread_rng);
  }
  return unif_rand();
}

void rng_state_begin() {
  if (!thread_rng) {
    GetRNGstate();
  }
}

void rng_state_end() {
  if (!thread_rng) {
    PutRNGstate();
  }
}

} // namespace ts
