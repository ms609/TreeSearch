#ifndef TS_HEARTBEAT_H
#define TS_HEARTBEAT_H

// Intra-phase progress heartbeat.
//
// Phase-boundary reporting alone leaves a user watching total silence for as long
// as a single phase runs.  On a 182-tip, 420-character matrix with inapplicable
// tokens throughout, one measured TBR phase ran 582 s and one 3-cycle ratchet
// 549 s -- 96% of a 1173 s replicate spent in two phases that each print only on
// completion.  This lets those phases say "still working, here's the score".
//
// Thread safety.  Rprintf is an R API call and must never be made from a worker
// thread.  Every entry point here is a no-op unless the caller is the R main
// thread, tested via ts::thread_stop_flag (null only in serial mode -- see
// ts_rng.h).  The parallel path reports from its coordinating thread instead
// (ts_parallel.cpp), which is already the only thread that touches R there.
//
// State is deliberately plain globals rather than thread_local: only the main
// thread ever writes them, workers return before reading, and TreeSearch has a
// history of MinGW emutls corruption from thread_local objects (PR #253).

#include "ts_rng.h"

namespace ts {

namespace hb {
// True between heartbeat_begin() and heartbeat_end().  Written only by the main
// thread; worker threads never read it (they return at the thread test first).
extern bool active;

// Out-of-line slow path: decides whether enough time has passed, and prints.
void tick(const char* label, double score, int stride);
}  // namespace hb

// Arm the heartbeat for a serial search.  Disabled at verbosity 0, and by
// setting TS_HEARTBEAT_SECONDS=0.  No-op when called from a worker thread: the
// parallel path reports from its coordinator instead.
//
// The default cadence differs by destination, taken from TS_HEARTBEAT_SECONDS if
// set: 30 s to a terminal, where the line overwrites itself and costs nothing to
// leave on; 120 s to a batch log, where every heartbeat is a permanent line.
void heartbeat_begin(int verbosity);

// The cadence heartbeat_begin() would choose for this destination, so the
// parallel coordinator can poll on the same schedule without duplicating the
// TS_HEARTBEAT_SECONDS parsing.  Returns 0 when the heartbeat is disabled.
double heartbeat_interval(bool isTty);

// Reset the in-phase timer and clear any open line, so the next phase's
// heartbeat reports time spent in *that* phase rather than since the search
// began.  Call at each phase boundary.
void heartbeat_phase(const char* label);

// Suspend/resume reporting around a search whose score is NOT on the user's
// objective -- notably the ratchet's perturbed TBR, which runs on a reweighted
// matrix and can report a score far below the true optimum.  Emitting those
// would make the score appear to jump around meaninglessly.  Nests via a depth
// counter, so an inner suspend cannot resume an outer one.
void heartbeat_suspend();
void heartbeat_resume();

// RAII form of the above.
struct HeartbeatSuspend {
  HeartbeatSuspend() { heartbeat_suspend(); }
  ~HeartbeatSuspend() { heartbeat_resume(); }
  HeartbeatSuspend(const HeartbeatSuspend&) = delete;
  HeartbeatSuspend& operator=(const HeartbeatSuspend&) = delete;
};

// Disarm.  Safe to call unconditionally, including when never armed.
void heartbeat_end();

// Scope guard, so that an early return or an interrupt-unwound stack cannot
// leave the heartbeat armed with a half-written line on the console.
struct HeartbeatScope {
  explicit HeartbeatScope(int verbosity) { heartbeat_begin(verbosity); }
  ~HeartbeatScope() { heartbeat_end(); }
  HeartbeatScope(const HeartbeatScope&) = delete;
  HeartbeatScope& operator=(const HeartbeatScope&) = delete;
};

// Call from inside a long-running phase.  `label` names the phase ("TBR",
// "Ratchet"); `score` is the phase's current best.  `stride` amortizes the clock
// read: only every `stride`-th call consults the clock, so pass a large stride
// from tight loops (per-clip) and 1 from coarse ones (per-ratchet-cycle).
inline void heartbeat(const char* label, double score, int stride = 64) {
  // Worker threads must not reach the R API.  This test first, so `hb::active`
  // is never read off-thread.
  if (thread_stop_flag != nullptr) return;
  if (!hb::active) return;
  hb::tick(label, score, stride);
}

}  // namespace ts

#endif  // TS_HEARTBEAT_H
