# Parallel fan-out for the converse-constraint Bremer engine.
#
# Each clade's converse-constraint search is fully independent of the others, so
# the dominant cost of `Bremer(method = "constraint")` -- N sequential
# `MaximizeParsimony()` searches -- parallelises cleanly ACROSS clades.  It must
# parallelise across clades and NOT within a single search: every converse
# search runs serially (`nThreads = 1`) because the negative-constraint pool
# guard (`TreePool::set_forbidden`) lives on the serial search path only.
#
# `.BremerConverseScores()` runs the N per-clade tasks -- optionally over a
# user-supplied `parallel` cluster -- and collects their results.  The per-clade
# task (`processFn`, i.e. `.BremerConstraint()`'s `processConverse`) runs one
# converse search, maps the engine's "no tree found" sentinel to NA and verifies
# the returned tree lacks the clade, so this helper stays agnostic to that logic
# and can be validated in isolation against the enumeration oracle.
#
# Determinism: each clade is seeded independently (seeds drawn once, in order,
# from the caller's RNG stream).  Seeding per CLADE rather than per WORKER makes
# the result independent of how tasks are scheduled onto workers, so a parallel
# run reproduces a serial one under the same `set.seed()` -- serial
# `MaximizeParsimony()` is itself reproducible under `set.seed()`.
#
# @param n Number of clades (tasks).
# @param processFn `function(i)` returning the (scalar) Bremer contribution of
#   clade `i` -- a length, or `NA_real_`.  The closure carries the reference,
#   dataset and scoring/search arguments.
# @param cl A `parallel` cluster whose workers have TreeSearch loaded, or `NULL`
#   for a seeded serial fan-out.
# @param seeds Optional integer per-clade RNG seeds; drawn reproducibly from the
#   caller's stream when `NULL`.
# @return A numeric vector of length `n`, one entry per clade.
.BremerConverseScores <- function(n, processFn, cl = NULL, seeds = NULL) {
  if (n == 0L) {
    return(numeric(0))
  }
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, n)
  }

  # One task per clade.  Seed inside the task so the outcome depends only on the
  # clade index, never on which worker (or in what order) runs it.
  worker <- function(i) {
    set.seed(seeds[[i]])
    processFn(i)
  }

  scores <- if (is.null(cl)) {
    lapply(seq_len(n), worker)
  } else {
    if (!inherits(cl, "cluster")) {
      stop("`cl` must be a `parallel` cluster or NULL.")
    }
    parallel::parLapply(cl, seq_len(n), worker)
  }

  vapply(scores, function(s) {
    if (length(s) == 0L) NA_real_ else as.numeric(s)[[1L]]
  }, double(1))
}
