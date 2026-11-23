# Count distinct COLLAPSED topologies among the saved Zhu MPTs, using the
# engine's own zero-length-branch collapse criterion (ts_collapsed_flags_debug,
# the same one add_collapsed uses).  This is the TNT-"collapse zero-length
# branches"-comparable count.
suppressMessages(devtools::load_all(".", quiet = TRUE))
trees <- readRDS("dev/profiling/zhu_probe_trees.rds")
nex <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013/zhu2013_orig.nex"
dat <- TreeTools::ReadAsPhyDat(nex)
dat <- TreeSearch:::.GapsAsMissing(dat)   # search used inapplicable = "missing"

at <- attributes(dat)
contrast <- at$contrast
tip_data <- matrix(unlist(dat, use.names = FALSE), nrow = length(dat), byrow = TRUE)
weight <- as.integer(at$weight)
levels <- at$levels
labs <- names(dat)

# Canonical collapsed Newick from engine parent[]/collapsed[] arrays.
collapsed_key <- function(tree) {
  tr <- TreeTools::Preorder(TreeTools::RenumberTips(tree, labs))
  if (tr[["edge"]][1L, 2L] > length(labs)) tr <- TreeTools::RootTree(tr, 1L)
  fl <- ts_collapsed_flags_debug(tr[["edge"]], contrast, tip_data, weight,
                                 levels, FALSE)
  parent <- fl$parent          # 0-based, root is self-parent
  collapsed <- fl$collapsed
  n_tip <- fl$n_tip
  n_node <- fl$n_node
  # children lists
  kids <- vector("list", n_node)
  root <- which(seq_len(n_node) - 1L == parent)[1] - 1L  # self-parent (0-based)
  for (i in seq_len(n_node) - 1L) {
    p <- parent[i + 1L]
    if (i != root) kids[[p + 1L]] <- c(kids[[p + 1L]], i)
  }
  # recurse, splicing collapsed internal nodes into their parent
  emit <- function(node) {
    if (node < n_tip) return(labs[node + 1L])     # tip
    ch <- kids[[node + 1L]]
    parts <- character(0)
    for (c in ch) {
      if (c >= n_tip && collapsed[c + 1L]) {
        # collapsed internal edge: splice this node's children up
        parts <- c(parts, vapply(kids[[c + 1L]], emit, character(1)))
      } else {
        parts <- c(parts, emit(c))
      }
    }
    paste0("(", paste(sort(parts), collapse = ","), ")")
  }
  emit(root)
}

resolved_key <- function(tree)
  ape::write.tree(TreeTools::SortTree(ape::unroot(tree)))

res <- vapply(trees, resolved_key, character(1))
col <- vapply(trees, collapsed_key, character(1))
cat("trees:", length(trees), "\n")
cat("distinct RESOLVED topologies :", length(unique(res)), "\n")
cat("distinct COLLAPSED topologies:", length(unique(col)), "\n")
# how many internal edges collapse, per tree (distribution)
