# corpus_state_heterogeneity.R -- quantify the packing lever's applicability.
#
# TS pads EVERY block to the dataset's GLOBAL applicable-state count
# (ts_data.cpp:214-250: blk.n_states = total_app_states for all blocks). TNT
# bins characters by state-count (2-bit/4-bit fields), so a binary char costs
# 2 bits even in a 9-state matrix. The packing lever (re-bin + per-block state
# alphabet) is worth something ONLY if real datasets have within-matrix state
# HETEROGENEITY: a low-state majority dragged up to a high global n_states.
#
# For each local dataset: global applicable n_states, per-character distinct
# applicable-state count, and a first-order "packing headroom":
#   current words/node   = n_blocks * global_n_states
#   binned  words/node   = sum over 64-char blocks (chars sorted by state-count)
#                          of (max distinct-state-count in block)
#   headroom = 1 - binned/current   (fraction of state-words removable)
# This upper-bounds the per-word cost a state-count-aware packing could recover.

suppressMessages({
  library(TreeSearch, lib.loc = ".agent-disc")
  library(TreeTools)
})

dsdir <- system.file("datasets", package = "TreeSearch", lib.loc = ".agent-disc")
files <- list.files(dsdir, pattern = "\\.nex$", full.names = TRUE)

analyse <- function(f) {
  nm <- sub("\\.nex$", "", basename(f))
  phy <- tryCatch(ReadAsPhyDat(f), error = function(e) NULL)
  if (is.null(phy) || !inherits(phy, "phyDat")) return(NULL)
  at <- attributes(phy)
  contrast <- at$contrast                     # n_tokens x n_states
  levels   <- at$levels
  gstates  <- length(levels)                  # global applicable states (incl "-"?)
  # per-character ALPHABET = states the char can actually take, EXCLUDING the
  # all-states "?" token (which maps to the global mask under the current global
  # remap but semantically = the char's own alphabet; representing it as the local
  # alphabet is Fitch-EXACT and needs fewer state-planes). Counting "?" as all
  # global states (the prior bug) inflated every missing-data char to n_states.
  mat <- do.call(rbind, phy)                  # n_tips x n_patterns (token indices)
  w   <- at$weight
  tok_ns   <- rowSums(contrast)               # states each token spans
  fullmask <- which(tok_ns == gstates)        # the all-states "?" token(s)
  perpat <- integer(ncol(mat))
  for (p in seq_len(ncol(mat))) {
    toks <- unique(mat[, p]); toks <- toks[!(toks %in% fullmask)]
    if (!length(toks)) { perpat[p] <- 1L; next }
    perpat[p] <- sum(colSums(contrast[toks, , drop = FALSE]) > 0)  # char's own alphabet
  }
  # expand by weight to characters
  sc <- rep(perpat, times = w)
  nChar <- length(sc)
  # current words/node (global width, all blocks)
  nBlk_cur <- ceiling(nChar / 64)
  cur <- nBlk_cur * gstates
  # binned: sort chars by state-count, 64/block, block width = max in block
  scs <- sort(sc)
  binned <- 0L
  i <- 1L
  while (i <= nChar) {
    j <- min(i + 63L, nChar)
    binned <- binned + max(scs[i:j])
    i <- j + 1L
  }
  data.frame(dataset = nm, nTip = length(phy), nChar = nChar,
             global_ns = gstates,
             med_char_ns = stats::median(sc), mean_char_ns = round(mean(sc), 2),
             pct_binary = round(100 * mean(sc <= 2), 1),
             pct_le4 = round(100 * mean(sc <= 4), 1),
             cur_words = cur, binned_words = binned,
             headroom = round(1 - binned / cur, 3))
}

rows <- lapply(files, function(f) tryCatch(analyse(f), error = function(e) NULL))
res <- do.call(rbind, Filter(Negate(is.null), rows))
res <- res[order(-res$headroom), ]
options(width = 200)
cat(sprintf("Analysed %d datasets\n\n", nrow(res)))
print(res, row.names = FALSE)
cat("\n=== summary ===\n")
cat(sprintf("datasets with global_ns >= 5: %d\n", sum(res$global_ns >= 5)))
cat(sprintf("median packing headroom (all): %.3f\n", median(res$headroom)))
cat(sprintf("median headroom | global_ns>=5: %.3f\n",
            median(res$headroom[res$global_ns >= 5])))
cat(sprintf("datasets headroom > 0.30: %d (%.0f%%)\n",
            sum(res$headroom > 0.30), 100 * mean(res$headroom > 0.30)))
cat(sprintf("datasets headroom > 0.50: %d (%.0f%%)\n",
            sum(res$headroom > 0.50), 100 * mean(res$headroom > 0.50)))
saveRDS(res, "dev/profiling/reeval/corpus_state_heterogeneity.rds")
