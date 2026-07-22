# profvis + Rprof triage for the standard-Fitch search call (R vs C++).
# Confirms there is no per-replicate R loop to port: the search is a single
# .Call, so R overhead is one-time data prep only.
LIBDIR <- "dev/profiling/.vtune-lib-20260616051420"
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))
suppressMessages(library(profvis))

raw <- inapplicable.phyData[["Zhu2013"]]
m <- PhyDatToMatrix(raw, ambigNA = FALSE); m[m == "-"] <- "?"
dataset <- MatrixToPhyDat(m)

# Warm up (compile/JIT, load DLL) outside measurement
invisible(suppressWarnings(MaximizeParsimony(dataset, maxReplicates = 1L,
          nThreads = 1L, strategy = "auto", verbosity = 0L)))

p <- profvis::profvis({
  set.seed(1)
  invisible(suppressWarnings(MaximizeParsimony(dataset, maxReplicates = 6L,
            nThreads = 1L, strategy = "auto", verbosity = 0L)))
})
htmlwidgets::saveWidget(p, "dev/profiling/drivers/fitch-tnt-profvis.html",
                        selfcontained = FALSE)

# Numeric R-vs-native split via Rprof
tf <- tempfile()
Rprof(tf, interval = 0.005, line.profiling = FALSE)
set.seed(1)
invisible(suppressWarnings(MaximizeParsimony(dataset, maxReplicates = 6L,
          nThreads = 1L, strategy = "auto", verbosity = 0L)))
Rprof(NULL)
s <- summaryRprof(tf)
cat("\n=== Top self-time (by.self) ===\n")
print(utils::head(s$by.self, 12))
cat("\n total.time:", s$sampling.time, "s\n")
