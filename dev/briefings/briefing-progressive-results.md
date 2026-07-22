# Briefing: Progressive Search Result Display

**Task:** T-129  
**Author:** Agent A  
**Date:** 2026-03-19

---

## Summary

**Recommendation: implement progress file polling using the existing C++ callback infrastructure.**

This approach reuses the cancel-file pattern already in place, requires minimal new
code on both the C++ and Shiny sides, and gives users real-time per-replicate feedback
without architectural changes or streaming intermediates.

Do **not** stream partial tree results to the UI mid-search. The benefits are marginal
and the implementation cost is high (see below).

---

## Context

Searches are invoked via `ExtendedTask` wrapping `future::future({MaximizeParsimony(...)})`.
The future runs in a separate R process — no reactive communication until it resolves.
The only currently-visible progress signal is a static "Searching…" notification and
a frozen output panel.

For short searches (<10 s), this is not a problem. For long searches (large datasets,
many replicates, long timeouts), users have no feedback that the search is alive or
making progress.

---

## What Already Exists

### C++ side: `progress_callback` (ts_driven.h / ts_rcpp.cpp)

`DrivenParams::progress_callback` is an `std::function<void(const ProgressInfo&)>`
already called after every phase and every replicate. The `ProgressInfo` struct carries:

```
replicate       — current replicate (1-based)
max_replicates  — configured max
best_score      — pool best so far
hits_to_best    — independent discoveries of best
target_hits     — convergence target
pool_size       — trees currently in pool
phase           — "replicate", "done", "tbr", "ratchet", etc.
elapsed_seconds — wall time since search start
phase_score     — score after this phase
```

The Rcpp bridge (`ts_rcpp.cpp` lines 1360–1375) already accepts an optional R function
and wraps it into this callback. **The infrastructure is complete — it just isn't used
by MaximizeParsimony() yet.**

### Shiny side: file-polling pattern (mod_search.R profile prep)

`profilePrepTask` already uses:
1. `tempfile()` progress path passed to the background task
2. `invalidateLater(500)` observer polling the file every 500ms
3. Notification update on each poll

This is exactly the right pattern for search progress too.

---

## Recommended Approach: Progress File Polling

### How it works

1. Before invoking `searchTask`, create a `progressPath` temp file.
2. In the background future, after `MaximizeParsimony()` sets `TREESEARCH_CANCEL_FILE`,
   also set a `TREESEARCH_PROGRESS_FILE` environment variable.
3. In `MaximizeParsimony()` (R level), if `TREESEARCH_PROGRESS_FILE` is set, pass an R
   function as `progressCallback` to `ts_driven_search()`. This callback writes
   a single line to the file after each replicate: `{rep} {max_rep} {best_score} {hits}`.
4. The main Shiny process polls `progressPath` every 500ms. On each poll, update the
   notification text.

### What the user sees

Currently:  
> `Searching (50 runs, k=6, 2 threads)…`

With progress polling:  
> `Searching… Rep 15/50 | Best: 42 | 3 hits`

Or with elapsed time:  
> `Searching… Rep 15/50 | Best: 42 | 3 hits | 8.2s elapsed`

When targetHits is reached before maxReplicates, this naturally shows convergence:  
> `Searching… Rep 23/50 | Best: 42 | 5/5 hits ✔ (wrapping up…)`

### Implementation path

**R package changes (MaximizeParsimony.R):**
- Check `Sys.getenv("TREESEARCH_PROGRESS_FILE")` before calling ts_driven_search
- If set, construct a `progressCallback` function that `writeLines()` to the file
  on `phase == "replicate"` events only (skip phase-level noise)
- ~20 lines R

**Shiny changes (mod_search.R):**
- Add `progressPath <- tempfile()` and pass it to the future alongside `cancelPath`
- Set `Sys.setenv(TREESEARCH_PROGRESS_FILE = progressPath)` in the future alongside
  the cancel env var
- Add `invalidateLater(500)` observer (mirroring the profile prep observer) that
  reads and parses the progress file
- On read, update the `r$searchNotification` message text
- ~30 lines R, no new UI elements

**C++ changes:** None required. The existing `progress_callback` / `progressCallback`
infrastructure handles everything.

**Total estimated effort:** ~2–3 hours.

---

## What NOT to Build: Partial Tree Streaming

A common wish is to "show best trees so far" during a search. This sounds appealing
but has significant problems:

### The pool is not intermediate-result-safe

The internal `TreePool` accumulates trees across replicates. At any mid-search point,
the pool contains a **subset of replicates' local optima** — not the final MPT set.
Trees in the pool at rep 15/50 may be suboptimal relative to the final result; the
tree topology at the "current best score" may not even survive MPT enumeration.

Displaying these trees as search results would be misleading: users might interpret
them as MPTs, save them, or make decisions based on incomplete evidence.

### R-level chunking doesn't help

Splitting the search into multiple short `MaximizeParsimony()` calls (each returning
a partial result) is tempting but:
- The pool and search state don't persist across calls (each call starts fresh)
- The quality/time tradeoff from very short searches is poor (no ratchet convergence)
- This is exactly what the "Continue search" button already provides at the user level

If users want intermediate trees, "Continue search" with small `maxReplicates` already
achieves this.

### The right display is convergence status

What users actually need to know mid-search is not *which trees* are in the pool, but:
- Is the search still running? (alive check)
- Is the best score improving? (convergence progress)
- How many hits to the best score so far? (convergence confidence)

All three are available from `ProgressInfo` at the replicate level with no new C++
work.

---

## Secondary Improvement: Elapsed Timer (Trivially Easy)

Even without the C++ callback, a simple elapsed-time counter can be added with zero
package changes:

```r
# In mod_search.R, near the searchInProgress observer:
observe({
  req(r$searchInProgress)
  invalidateLater(1000)  # fire every second
  elapsed <- as.integer(difftime(Sys.time(), r$searchStartTime, units = "secs"))
  # Update notification text: "Searching… (42s elapsed)"
})
```

This costs ~10 lines and prevents "is the app frozen?" uncertainty. However, it
provides no information about progress — a 5-minute search with a frozen score
gives the user no convergence signal. The file-polling approach is clearly superior.

---

## Decision Matrix

| Approach | Effort | Value | Verdict |
|----------|--------|-------|---------|
| Elapsed timer only | ~10 lines | Low — no convergence info | Not worth it alone |
| **Progress file polling** | **~50 lines** | **High — reps + score + hits** | **✅ Recommended** |
| Partial tree streaming | ~200+ lines + arch changes | Low — misleads user | ✗ Do not build |
| R-level chunking | ~150+ lines + pool state | Medium — duplicates "Continue" | ✗ Redundant |

---

## Concrete Task Proposal

File as **T-141** (P3):

> **Shiny: Per-replicate search progress display**  
> Use existing `progress_callback` / `TREESEARCH_PROGRESS_FILE` env var pattern
> (mirrors cancel file + profile prep). MaximizeParsimony() writes rep/score/hits
> to file on each replicate. mod_search.R polls every 500ms during search.
> Result: notification updates from static "Searching…" to live "Rep 15/50 | Best: 42 | 3 hits".
> No C++ changes needed.  
> Estimate: ~2–3 hours.
