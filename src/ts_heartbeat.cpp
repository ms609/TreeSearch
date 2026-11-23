#include "ts_heartbeat.h"

// <R.h> rather than just <R_ext/Print.h>: R_FlushConsole is not declared by the
// latter.  This is the same include ts_parallel.cpp uses for the same call.
#include <R.h>
#include <chrono>
#include <cstdlib>

// Detect a terminal without calling fileno(stdout): referencing `stdout` from
// compiled code draws a CRAN NOTE.  Same approach as ts_parallel.cpp.
#ifdef _WIN32
  #include <io.h>
  #define TS_HB_ISATTY()  (_isatty(1) != 0)
#else
  #include <unistd.h>
  #define TS_HB_ISATTY()  (isatty(1) != 0)
#endif

namespace ts {

namespace hb {

bool active = false;

namespace {

using Clock = std::chrono::steady_clock;

double interval_ = 0;
long long counter_ = 0;
Clock::time_point last_;
Clock::time_point phase_start_;
// True when an unterminated \r line is on the console, so it can be cleared
// before the phase prints its own summary.
bool lineOpen_ = false;
bool isTty_ = false;
// >0 while inside a search whose score is not on the user's objective.
int suspendDepth_ = 0;

void CloseLine() {
  if (!lineOpen_) return;
  // Blank the line, then return to column 0 so the next write starts clean.
  Rprintf("\r%-72s\r", "");
  lineOpen_ = false;
}

}  // namespace

void tick(const char* label, double score, int stride) {
  if (suspendDepth_ > 0) return;
  if (stride < 1) stride = 1;
  if (++counter_ % stride != 0) return;

  const Clock::time_point now = Clock::now();
  const double sinceLast =
      std::chrono::duration_cast<std::chrono::duration<double>>(now - last_)
          .count();
  if (sinceLast < interval_) return;
  last_ = now;

  const double inPhase =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          now - phase_start_).count();

  if (isTty_) {
    // Overwrite in place: a 10-minute phase would otherwise scroll the console.
    Rprintf("\r    %s: best %.5g, %.0fs in phase      ", label, score, inPhase);
    R_FlushConsole();
    lineOpen_ = true;
  } else {
    // Batch log (Rscript, SLURM): \r would collapse the file into one unreadable
    // line, so emit discrete lines instead.
    Rprintf("    %s: best %.5g, %.0fs in phase\n", label, score, inPhase);
  }
}

void PhaseBegin() {
  CloseLine();
  phase_start_ = Clock::now();
  last_ = phase_start_;
  counter_ = 0;
}

void Finish() { CloseLine(); }

}  // namespace hb

double heartbeat_interval(bool isTty) {
  // A terminal line overwrites itself, so it can be cheap and frequent; a batch
  // log keeps every line, so slow it down to stay readable over a multi-hour run.
  double interval = isTty ? 30.0 : 120.0;
  if (const char* env = std::getenv("TS_HEARTBEAT_SECONDS")) {
    if (*env != '\0') {
      char* end = nullptr;
      const double parsed = std::strtod(env, &end);
      // Reject junk rather than silently treating it as 0 (which would disable
      // the heartbeat and look like the feature is broken).
      if (end != env && *end == '\0' && parsed >= 0) interval = parsed;
    }
  }
  return interval > 0 ? interval : 0;
}

void heartbeat_begin(int verbosity) {
  if (thread_stop_flag != nullptr) return;  // worker thread: never arm
  hb::active = false;
  if (verbosity < 1) return;

  hb::isTty_ = TS_HB_ISATTY();
  const double interval = heartbeat_interval(hb::isTty_);
  if (interval <= 0) return;  // explicitly disabled

  hb::interval_ = interval;
  hb::lineOpen_ = false;
  hb::PhaseBegin();
  hb::active = true;
}

void heartbeat_phase(const char* /*label*/) {
  if (thread_stop_flag != nullptr || !hb::active) return;
  hb::PhaseBegin();
}

void heartbeat_suspend() {
  if (thread_stop_flag != nullptr) return;
  // Clear any in-place line now: the suspended search may run for a while, and a
  // stale score should not sit on the console implying it is current.
  hb::CloseLine();
  ++hb::suspendDepth_;
}

void heartbeat_resume() {
  if (thread_stop_flag != nullptr) return;
  if (hb::suspendDepth_ > 0) --hb::suspendDepth_;
}

void heartbeat_end() {
  if (thread_stop_flag != nullptr) return;
  hb::Finish();
  hb::active = false;
}

}  // namespace ts
