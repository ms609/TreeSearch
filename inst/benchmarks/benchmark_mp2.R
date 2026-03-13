# Benchmark: MaximizeParsimony2 (C++ driven search) vs MaximizeParsimony (R loop)
#
# Compares wall-clock time and best score found on a selection of datasets
# from inapplicable.phyData, using equal-weight Fitch parsimony throughout.

library(TreeSearch)
library(TreeTools)

data("inapplicable.phyData")

#' Convert inapplicable tokens to fully ambiguous for pure Fitch EW scoring
#' @param ds A phyDat object
#' @return The modified phyDat with "-" treated as "?"
strip_inapp <- function(ds) {
  cont <- attr(ds, "contrast")
  lvls <- attr(ds, "levels")
  dash_col <- which(lvls == "-")
  if (length(dash_col) == 0L) return(ds)
  # Tokens that code for "-": make them fully ambiguous over applicable states
  has_dash <- cont[, dash_col] == 1
  app_cols <- setdiff(seq_len(ncol(cont)), dash_col)
  cont[has_dash, app_cols] <- 1
  # Drop the "-" state column
  cont <- cont[, -dash_col, drop = FALSE]
  attr(ds, "contrast") <- cont
  attr(ds, "levels") <- lvls[-dash_col]
  ds
}

bench_datasets <- c(
  "Vinther2008",    # 23 tips, 50 chars
  "Asher2005",      # 23 tips, 125 chars
  "Wortley2006",    # 37 tips, 105 chars
  "Wills2012",      # 55 tips, 87 chars
  "Agnarsson2004",  # 62 tips, 225 chars
  "Dikow2009"       # 88 tips, 204 chars
)

results <- data.frame(
  dataset = character(), tips = integer(), patterns = integer(),
  mp2_score = numeric(), mp1_score = numeric(), score_diff = numeric(),
  mp2_time = numeric(), mp1_time = numeric(), speedup = numeric(),
  stringsAsFactors = FALSE
)

for (nm in bench_datasets) {
  ds <- strip_inapp(inapplicable.phyData[[nm]])
  n_tip <- length(ds)
  n_pat <- attr(ds, "nr")
  cat("\n---", nm, "(", n_tip, "tips,", n_pat, "pat) ---\n")

  # --- MaximizeParsimony2 (C++ driven search) ---
  set.seed(6218)
  t2 <- system.time({
    r2 <- MaximizeParsimony2(ds, verbosity = 0L)
  })
  s2 <- TreeLength(r2[[1]], ds)

  # --- MaximizeParsimony (R loop) ---
  set.seed(6218)
  t1 <- system.time({
    r1 <- MaximizeParsimony(ds, ratchIter = 7L, tbrIter = 2L,
                            maxHits = n_tip * 1.8, maxTime = 5,
                            verbosity = 0L)
  })
  s1 <- TreeLength(r1[[1]], ds)

  cat("  MP2:", s2, sprintf("(%.2fs, %d reps)", t2["elapsed"],
      attr(r2, "replicates")),
      "  MP1:", s1, sprintf("(%.2fs)", t1["elapsed"]),
      "  diff:", s2 - s1, "\n")

  results <- rbind(results, data.frame(
    dataset = nm, tips = n_tip, patterns = n_pat,
    mp2_score = s2, mp1_score = s1, score_diff = s2 - s1,
    mp2_time = t2["elapsed"], mp1_time = t1["elapsed"],
    speedup = t1["elapsed"] / t2["elapsed"],
    stringsAsFactors = FALSE
  ))
}

cat("\n\n=== SUMMARY ===\n")
print(results, row.names = FALSE)
