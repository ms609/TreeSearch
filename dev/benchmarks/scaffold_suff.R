#!/usr/bin/env Rscript
# Mission A angle #1, SCAFFOLD-SUFFICIENCY (advisor: payoff-before-feasibility).
# Constrain MaximizeParsimony on nested subsets of the 1943 deep-backbone splits
# (deepest-k of the ~35 min-side>=100 splits), RAS-fill the rest, matched budget,
# measure hit-rate to <=1944. Degradation curve = the discriminator:
#   need ~all 35 to route in  -> blind-find ~= full search -> clique route DEAD (cheap).
#   5-10 suffice              -> enumeration pipeline worth building; target size known.
# NOT circular: uses true splits to measure ROUTABILITY (a ceiling), not to claim a
# blind method finds them. Self-check: FULL-tree constraint MUST return 1943 (honored).
# Prior: BSS K_ceiling=323 leans negative for small k, but that froze clades+reoptimised;
# split-constrain + RAS-build is a different operation.
suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(ape)}))
task <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
NEX <- "/nobackup/pjjg18/lgsweep/matrices/project5432.nex"
FLOOR <- "/nobackup/pjjg18/floors/project5432_tnt_floor_one.tre"
OUT <- "/nobackup/pjjg18/reeval/scaffold_out"; dir.create(OUT, showWarnings = FALSE)

pd <- ReadAsPhyDat(NEX); tips <- names(pd); nt <- length(tips)
TL <- function(t) TreeLength(Preorder(t), pd, concavity = Inf, inapplicable = "missing")
tnt <- ReadTntTree(FLOOR, tipLabels = tips); if (inherits(tnt, "multiPhylo")) tnt <- tnt[[1]]
tnt <- Preorder(tnt); s1943 <- TL(tnt)
stopifnot(abs(s1943 - 1943) <= 3)                       # floor sanity (self-check the anchor)

spl <- as.Splits(tnt, tipLabels = tips); tis <- TipsInSplits(spl); ms <- pmin(tis, nt - tis)
allNodes <- as.integer(names(tis)); deepMask <- ms >= 100
deepNodes <- allNodes[deepMask]; deepByDepth <- deepNodes[order(ms[deepMask], decreasing = TRUE)]
nDeep <- length(deepByDepth)

LEV <- c(0, 5, 10, 20, nDeep); SEEDS <- 12
grid <- expand.grid(k = LEV, seed = seq_len(SEEDS)); ng <- nrow(grid)
if (task <= ng) { k <- grid$k[task]; seed <- grid$seed[task]; full <- FALSE
} else { full <- TRUE; k <- NA; seed <- task - ng }               # tasks ng+1.. = FULL self-check
set.seed(2000 + task)

if (full) { ctree <- tnt; ck <- "FULL"
} else if (k == 0) { ctree <- NULL; ck <- "RAS0"
} else { keep <- deepByDepth[seq_len(min(k, nDeep))]; ctree <- CollapseNode(tnt, setdiff(allNodes, keep)); ck <- paste0("k", k) }

STRAT <- "thorough"; R <- 8L; CAP <- 1800
t0 <- Sys.time()
res <- tryCatch(MaximizeParsimony(pd, constraint = ctree, concavity = Inf, inapplicable = "missing",
    strategy = STRAT, maxReplicates = R, targetHits = NULL, maxSeconds = CAP,
    nThreads = 1L, verbosity = 0L, collapse = FALSE), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.null(res)) { cat("NULL result task", task, "\n"); quit(save = "no") }
trees <- if (inherits(res, "phylo")) list(res) else res
sc <- vapply(trees, TL, numeric(1)); best <- min(sc); bt <- Preorder(trees[[which.min(sc)]])

## constraint-held check: fraction of constrained splits present in the best tree
held <- NA_real_
if (!is.null(ctree)) {
  sf <- tryCatch(SplitFrequency(ctree, structure(list(bt), class = "multiPhylo")), error = function(e) NULL)
  if (!is.null(sf)) held <- mean(!is.na(sf) & sf >= 1)
}
rp <- attr(res, "replicate_scores"); if (is.null(rp)) rp <- attr(res, "scores"); if (is.null(rp)) rp <- NA
cat(sprintf("task=%d %s seed=%d best=%.0f hit1944=%d hit1943=%d held=%.2f el=%.0fs reps=%s\n",
    task, ck, seed, best, best <= 1944, best <= 1943, held, el, paste(round(rp), collapse = ";")))
write.csv(data.frame(task = task, level = ck, k = ifelse(full, nDeep + 1, k), seed = seed,
    best = best, hit1944 = (best <= 1944), hit1943 = (best <= 1943), held = held,
    elapsed = el, ntip_deep = nDeep, repscores = paste(round(rp), collapse = ";")),
    file.path(OUT, sprintf("cell_%03d.csv", task)), row.names = FALSE)
