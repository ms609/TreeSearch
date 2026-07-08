# Convert off-corpus .nex matrices to verified phyDat .rds for the addseq
# confirmation array (bench_addseq_cell.R with TS_EXTRADIR). Run ONCE,
# single-tenant. For each entry it reads the .nex, saves a phyDat .rds into
# TS_EXTRADIR, and runs a single MaximizeParsimony replicate to confirm the
# matrix loads and scores under the native engine without error.
#
# Note: whether "-" resolves to inapplicable or missing does NOT affect the
# A/B validity -- all three wagnerBias arms see the identical phyDat, so the
# mode-vs-mode comparison is internally consistent regardless. The loaded
# levels are printed for transparency (na-validation-alignment-gotcha).
#
# Env:
#   TS_LIB        -- addseq probe lib (lib.loc for TreeSearch)
#   TS_EXTRADIR   -- output directory for the .rds objects
#   TS_EXTRA_SPEC -- "name=path;name=path" matrices to convert
# Local test: TS_EXTRA_SPEC="lobo=dev/benchmarks/data-na/lobo.nex" Rscript ...

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})

extradir <- Sys.getenv("TS_EXTRADIR", "dev/benchmarks/extra_addseq")
spec <- Sys.getenv("TS_EXTRA_SPEC", "")
if (!nzchar(spec)) stop("TS_EXTRA_SPEC not set (expected \"name=path;name=path\")")
dir.create(extradir, showWarnings = FALSE, recursive = TRUE)

entries <- strsplit(trimws(spec), ";")[[1]]
ok <- character(0)
for (e in entries) {
  kv <- strsplit(trimws(e), "=")[[1]]
  name <- trimws(kv[[1]]); path <- trimws(kv[[2]])
  cat(sprintf("\n=== %s <- %s ===\n", name, path))
  if (!file.exists(path)) { cat("  MISSING .nex, skipping\n"); next }
  res <- tryCatch({
    d <- ReadAsPhyDat(path)
    cat(sprintf("  tips=%d chars=%d levels={%s}\n",
                length(d), sum(attr(d, "weight")),
                paste(attr(d, "levels"), collapse = ",")))
    saveRDS(d, file.path(extradir, paste0(name, ".rds")))
    set.seed(1)
    r <- suppressWarnings(MaximizeParsimony(d,
      control = SearchControl(wagnerBias = 5L, wagnerBiasTemp = 0.3,
                              adaptiveStart = FALSE, wagnerStarts = 1L),
      maxReplicates = 1L, targetHits = 999L, maxSeconds = 0,
      nThreads = 1L, verbosity = 0L))
    cat(sprintf("  smoke score (1 rep, bias=5): %g -> saved %s.rds\n",
                attr(r, "score"), name))
    TRUE
  }, error = function(err) { cat("  ERROR:", conditionMessage(err), "\n"); FALSE })
  if (isTRUE(res)) ok <- c(ok, name)
}
cat(sprintf("\nConverted %d/%d: %s\n", length(ok), length(entries),
            paste(ok, collapse = ", ")))
