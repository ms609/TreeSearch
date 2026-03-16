#ifndef TS_RNG_H
#define TS_RNG_H

// Thread-safe RNG and interrupt checking for parallel search.
//
// In serial mode (thread_rng == nullptr), these helpers fall back to
// R's GetRNGstate()/unif_rand()/PutRNGstate() and R_CheckUserInterrupt().
//
// In parallel mode, worker threads set thread_rng to their own mt19937
// and thread_stop_flag to a shared atomic. No R API calls are made from
// worker threads.

#include <random>
#include <atomic>

namespace ts {

// Thread-local RNG pointer. When non-null, make_rng() draws seeds from
// this instead of R's unif_rand(). Set by parallel worker threads.
extern thread_local std::mt19937* thread_rng;

// Thread-local stop flag pointer. When non-null, check_interrupt()
// tests this instead of calling R_CheckUserInterrupt(). Set by parallel
// worker threads.
extern thread_local const std::atomic<bool>* thread_stop_flag;

// Return a seeded std::mt19937.
// - Serial mode: seeds from R's unif_rand() (with Get/PutRNGstate)
// - Parallel mode: seeds from *thread_rng
std::mt19937 make_rng();

// Check for user interrupt or stop signal.
// - Serial mode: calls R_CheckUserInterrupt() (may longjmp)
// - Parallel mode: checks *thread_stop_flag; returns true if set
// Returns true only in parallel mode when stop is requested.
// In serial mode, returns false (R_CheckUserInterrupt longjmps on
// interrupt rather than returning a value).
bool check_interrupt();

// Direct random draws for algorithms that use unif_rand() directly
// (e.g., Fisher-Yates in Wagner tree and resampling).
// - Serial mode: calls unif_rand() (caller must manage RNG state)
// - Parallel mode: draws from *thread_rng (normalised to [0, 1))
double thread_safe_unif();

// Begin RNG state access (serial mode: calls GetRNGstate())
// In parallel mode this is a no-op.
void rng_state_begin();

// End RNG state access (serial mode: calls PutRNGstate())
// In parallel mode this is a no-op.
void rng_state_end();

} // namespace ts

#endif // TS_RNG_H
