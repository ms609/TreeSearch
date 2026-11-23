# Mission-wide A/B for the getenv hoist + T-S6c levers: full MaximizeParsimony
# (ratchet + sectorial + TBR), NOT isolated sectorial. Confirms the per-clip
# getenv finding is cross-cutting (every tbr_search clip, mission-wide), and that
# score is unchanged (byte-identical levers). Run with TREESEARCH_VTUNE_LIB.
LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB", unset = ".agent-sect")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))
ds_name <- Sys.getenv("TS_DATASET", unset = "Zhu2013")
reps    <- as.integer(Sys.getenv("TS_REPS", unset = "3"))
seed    <- as.integer(Sys.getenv("TS_SEED", unset = "1"))

raw <- inapplicable.phyData[[ds_name]]
m <- PhyDatToMatrix(raw, ambigNA = FALSE); m[m == "-"] <- "?"
dataset <- MatrixToPhyDat(m)

set.seed(seed)
t0 <- proc.time()
res <- suppressWarnings(MaximizeParsimony(dataset, maxReplicates = reps,
                                          nThreads = 1L, strategy = "thorough",
                                          verbosity = 0L))
el <- (proc.time() - t0)["elapsed"]
cat(sprintf("MISSION %s reps=%d seed=%d : elapsed=%.2f s score=%s\n",
            ds_name, reps, seed, el, attr(res, "score")))
