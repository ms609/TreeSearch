# Tier-2 behaviour-neutrality check: same seed + nThreads=1 => the dedup table
# must produce byte-identical search trajectories => identical scores on EVERY
# dataset and BOTH scoring paths (NA three-pass on raw data; standard Fitch
# after '-'->'?'). Run on after- and before-libs; scores must match exactly.
LIBDIR <- normalizePath(Sys.getenv("TREESEARCH_VTUNE_LIB"), winslash = "/")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))

sets <- c("Vinther2008", "Longrich2010", "Sansom2010", "Aria2015",
          "Dikow2009", "Zhu2013")

for (nm in sets) {
  raw <- inapplicable.phyData[[nm]]
  # NA three-pass path (raw inapplicables)
  set.seed(42)
  r_na <- suppressWarnings(MaximizeParsimony(raw, maxReplicates = 3L,
            nThreads = 1L, strategy = "default", verbosity = 0L))
  # Standard-Fitch path ('-' -> '?')
  m <- PhyDatToMatrix(raw, ambigNA = FALSE); m[m == "-"] <- "?"
  std <- MatrixToPhyDat(m)
  set.seed(42)
  r_std <- suppressWarnings(MaximizeParsimony(std, maxReplicates = 3L,
            nThreads = 1L, strategy = "default", verbosity = 0L))
  cat(sprintf("%-14s  NA=%-8s  STD=%-8s\n",
              nm, attr(r_na, "score"), attr(r_std, "score")))
}
