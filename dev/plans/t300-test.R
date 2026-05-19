library(TreeSearch, lib.loc = "C:/Users/pjjg18/AppData/Local/Temp/ts_t300_lib")
library(TreeTools)

set.seed(7531)

# Check whether congreveLamsdellMatrices has inapplicable chars
data("congreveLamsdellMatrices", package = "TreeSearch")
ds <- congreveLamsdellMatrices[[1]]
has_na <- "-" %in% unique(unlist(as.character(ds)))
cat("congreveLamsdellMatrices[[1]] has '-' tokens:", has_na, "\n")
cat("levels:", paste(attr(ds, "levels"), collapse = ","), "\n")
cat("contrast rows:", nrow(attr(ds, "contrast")), "\n")

# Build a pure-EW dataset (no NA) for sure: take Vinther2008 and replace - with ?
data("inapplicable.phyData", package = "TreeSearch")
ds_v <- inapplicable.phyData[["Vinther2008"]]
cat("\nVinther2008 levels:", paste(attr(ds_v, "levels"), collapse = ","), "\n")

# Build a synthetic pure-EW dataset using a random small matrix
mat <- matrix(sample(c("0","1"), 30 * 25, replace = TRUE),
              nrow = 30, ncol = 25,
              dimnames = list(paste0("t", 1:30), paste0("c", 1:25)))
ds_ew <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))
cat("\nPure EW synthetic dataset: n_taxa=", length(ds_ew),
    " n_chars=", attr(ds_ew, "nr"), "\n", sep = "")

cat("\n=== Pure EW TBR test ===\n")
set.seed(8421)
sink_path <- tempfile(fileext = ".txt")
sink(sink_path, split = TRUE)
result <- MaximizeParsimony(
  ds_ew,
  maxReplicates = 5L,
  targetHits = 20L,
  verbosity = 1L,
  nThreads = 1L
)
sink()
log_text <- readLines(sink_path)
n_debug <- sum(grepl("DEBUG_RESCORE", log_text))
cat("DEBUG_RESCORE mismatch lines emitted (EW):", n_debug, "(expect 0)\n")
if (n_debug > 0) cat(grep("DEBUG_RESCORE", log_text, value = TRUE), sep = "\n")
cat("Final score:", attr(result, "score"), "\n")

cat("\n=== Pure EW IW test (concavity=10) ===\n")
set.seed(8421)
sink_path2 <- tempfile(fileext = ".txt")
sink(sink_path2, split = TRUE)
result2 <- MaximizeParsimony(
  ds_ew,
  maxReplicates = 5L,
  targetHits = 20L,
  verbosity = 1L,
  nThreads = 1L,
  concavity = 10
)
sink()
log_text2 <- readLines(sink_path2)
n_debug2 <- sum(grepl("DEBUG_RESCORE", log_text2))
cat("DEBUG_RESCORE mismatch lines emitted (IW):", n_debug2, "(expect 0)\n")
if (n_debug2 > 0) cat(grep("DEBUG_RESCORE", log_text2, value = TRUE), sep = "\n")
cat("Final IW score:", attr(result2, "score"), "\n")
