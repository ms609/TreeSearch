# Tier-1 kill test (sectorial pattern-reduction lever):
# For each mission dataset, build a good tree, decompose into clades of size 6-50
# (the actual sector regime), and for each clade measure:
#   - n informative-within-sector patterns (>=2 distinct observed states among
#     the sector's tips; constant/single-state contribute a fixed cost offset)
#   - the COMPACT BLOCK COUNT after re-packing survivors (64 chars/block)
#     vs the original block count.
# Decision: if a typical size-25 Zanol sector still spans the SAME #blocks after
# re-pack, the lever is dead (no fewer SIMD blocks scanned -> no throughput win).

suppressMessages({
  library(TreeSearch)
  library(TreeTools)
})

CHARS_PER_BLOCK <- 64L

analyse <- function(nex, label, nstart = 1L) {
  dat <- ReadAsPhyDat(nex)
  ntip <- length(dat)
  # decode phyDat to a tip x char matrix of integer state codes; ambiguous/?
  # become NA (treated as "does not pin a state").
  cont <- attr(dat, "contrast")
  alev <- attr(dat, "allLevels")
  idx  <- attr(dat, "index")          # maps site -> pattern
  w    <- attr(dat, "weight")         # pattern frequency
  npat <- length(w)
  # For each token (row of contrast), how many states it is compatible with:
  ntok_states <- rowSums(cont > 0)
  # A token "pins" a single observed state only if it maps to exactly 1 state.
  # token -> the single state index (or NA if ambiguous/missing)
  single_state <- ifelse(ntok_states == 1L,
                         max.col(cont, ties.method = "first"), NA_integer_)
  # tip x pattern matrix of token codes
  M <- matrix(unlist(dat, use.names = FALSE), nrow = ntip, byrow = TRUE)
  # map token code -> pinned state (NA if ambiguous)
  Mstate <- matrix(single_state[M], nrow = ntip)

  # A good tree: parsimony-ratchet-free quick search start
  tr <- TreeTools::NJTree(dat)
  tr <- RootTree(tr, 1)
  # node clade sizes
  ce <- tr$edge
  nnode <- tr$Nnode
  # descendant tip count per node via postorder
  desc <- phangorn::Descendants(tr, seq_len(ntip + nnode), type = "tips")
  sizes <- lengths(desc)

  internal <- (ntip + 1):(ntip + nnode)
  sect_nodes <- internal[sizes[internal] >= 6 & sizes[internal] <= 50]

  orig_blocks <- ceiling(npat / CHARS_PER_BLOCK)

  rows <- lapply(sect_nodes, function(nd) {
    tips <- desc[[nd]]
    sub <- Mstate[tips, , drop = FALSE]   # |tips| x npat pinned states
    # informative within sector: >= 2 distinct non-NA states present
    n_distinct <- apply(sub, 2L, function(col) {
      u <- unique(col[!is.na(col)])
      length(u)
    })
    informative <- n_distinct >= 2L
    n_inf <- sum(informative)
    compact_blocks <- ceiling(n_inf / CHARS_PER_BLOCK)
    data.frame(node = nd, sector_tips = length(tips),
               n_inf = n_inf, npat = npat,
               orig_blocks = orig_blocks,
               compact_blocks = compact_blocks,
               block_reduction = orig_blocks - compact_blocks)
  })
  df <- do.call(rbind, rows)
  df$dataset <- label
  df
}

ds_dir <- "C:/Users/pjjg18/GitHub/TreeSearch/inst/datasets"
targets <- list(
  Zanol = file.path(ds_dir, "Zanol2014.nex"),
  Zhu   = file.path(ds_dir, "Zhu2013.nex"),
  Giles = file.path(ds_dir, "Giles2015.nex"),
  Wortley = file.path(ds_dir, "Wortley2006.nex")
)

all <- do.call(rbind, lapply(names(targets), function(nm) {
  cat("=== ", nm, " ===\n")
  tryCatch(analyse(targets[[nm]], nm), error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n"); NULL
  })
}))

cat("\n\n==== SUMMARY: block-count reduction by dataset ====\n")
for (nm in unique(all$dataset)) {
  d <- all[all$dataset == nm, ]
  cat(sprintf("\n%s: npat=%d, orig_blocks=%d, n_sectors=%d\n",
              nm, d$npat[1], d$orig_blocks[1], nrow(d)))
  cat(sprintf("  sector_tips: min=%d med=%.0f max=%d\n",
              min(d$sector_tips), median(d$sector_tips), max(d$sector_tips)))
  cat(sprintf("  n_inf (informative within sector): min=%d med=%.0f max=%d\n",
              min(d$n_inf), median(d$n_inf), max(d$n_inf)))
  cat(sprintf("  compact_blocks: min=%d med=%.1f max=%d\n",
              min(d$compact_blocks), median(d$compact_blocks), max(d$compact_blocks)))
  cat(sprintf("  block_reduction (orig-compact): med=%.1f  | sectors with reduction>=1: %d/%d (%.0f%%)\n",
              median(d$block_reduction),
              sum(d$block_reduction >= 1), nrow(d),
              100 * mean(d$block_reduction >= 1)))
  # focus on the ~size-25 sectors
  mid <- d[d$sector_tips >= 18 & d$sector_tips <= 32, ]
  if (nrow(mid) > 0)
    cat(sprintf("  [size 18-32 sectors, n=%d]: med n_inf=%.0f, med compact_blocks=%.1f, med reduction=%.1f\n",
                nrow(mid), median(mid$n_inf), median(mid$compact_blocks), median(mid$block_reduction)))
}

saveRDS(all, "C:/Users/pjjg18/GitHub/TreeSearch/dev/profiling/sector_pattern_killtest.rds")
