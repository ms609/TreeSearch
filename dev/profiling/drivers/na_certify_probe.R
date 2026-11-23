#!/usr/bin/env Rscript
# Where does NA certification actually fire under the SHIPPED recipe?
#
# The 97.7% figure in dev/profiling/na-exact-verify-dominates.md was measured
# through ts_tbr_search / ts_ratchet_search, whose TBRParams/RatchetParams leave
# tabu_size at 0.  do_reroot -- the gate on exact_verify_sweep -- requires
# tabu_size == 0, and the shipped presets set tabuSize = 100 (default) / 200
# (thorough).  So before gating certification we must know which call sites
# reach it in production at all: a null wall result from a flag that never
# fired is indistinguishable from a flag that fired and bought nothing.
#
# Reports MaximizeParsimony()'s `naDiag` (n_evs = certifications executed,
# n_evs_skipped = certifications the gate removed) with the sector / fuse /
# ratchet phases toggled off one at a time, which attributes the calls to a
# phase without needing per-site instrumentation.
#
# Usage: Rscript dev/profiling/drivers/na_certify_probe.R [lib] [dataset] [cells]
#   cells: comma-separated subset of the labels below.  The `tabu0` cell is
#   O(n^3)-expensive on large matrices -- drop it there.

args  <- commandArgs(trailingOnly = TRUE)
lib   <- if (length(args) >= 1L) args[[1]] else ".agent-p0"
key   <- if (length(args) >= 2L) args[[2]] else "Vinther2008"
pick  <- if (length(args) >= 3L) strsplit(args[[3]], ",", fixed = TRUE)[[1]] else NULL

library("TreeSearch", lib.loc = lib)
Sys.setenv(TS_NA_TIMING = "1")

e <- new.env()
load("data/inapplicable.phyData.rda", envir = e)
allDat <- get(ls(e)[1], envir = e)
dat <- allDat[[key]]
stopifnot(!is.null(dat))

runCell <- function(label, noCertify, ...) {
  if (noCertify) Sys.setenv(TS_NA_NOCERTIFY = "1") else Sys.unsetenv("TS_NA_NOCERTIFY")
  set.seed(1)
  t0 <- Sys.time()
  # Preset knobs go through DOTS, never `control =`: control= alongside
  # strategy= discards the preset (control-clobbers-strategy-preset).
  res <- MaximizeParsimony(dat, strategy = "default", verbosity = 0,
                           nThreads = 1L, maxReplicates = 2L, ...)
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  d <- attr(res, "naDiag")
  data.frame(label = label, noCertify = noCertify,
             score = attr(res, "score"),
             wall = round(wall, 2),
             nEvs = if (is.null(d)) NA_real_ else d$n_evs,
             nSkipped = if (is.null(d)) NA_real_ else d$n_evs_skipped,
             tEvsMs = if (is.null(d)) NA_real_ else round(d$t_evs_ms),
             stringsAsFactors = FALSE)
}

cells <- list(
  list("full-default",    list()),
  list("noSector",        list(xssRounds = 0L, rssRounds = 0L, cssRounds = 0L)),
  list("noFuse",          list(fuseInterval = 0L)),
  list("noSector+noFuse", list(xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
                               fuseInterval = 0L)),
  list("noRatchet",       list(ratchetCycles = 0L)),
  list("tabu0",           list(tabuSize = 0L))
)

if (!is.null(pick)) {
  cells <- cells[vapply(cells, function(cl) cl[[1]] %in% pick, logical(1))]
  stopifnot(length(cells) > 0L)
}

out <- do.call(rbind, lapply(cells, function(cl) {
  do.call(rbind, lapply(c(FALSE, TRUE), function(nc) {
    do.call(runCell, c(list(label = cl[[1]], noCertify = nc), cl[[2]]))
  }))
}))

cat(sprintf("\n== %s (%d taxa) ==\n", key, length(dat)))
print(out, row.names = FALSE)
cat("\nnEvs > 0 => certification fires under that config.\n")
cat("nSkipped > 0 => the certify_unrooted gate reached a live call site.\n")
