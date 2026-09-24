suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
bare <- "dev/benchmarks/tnt_bare"
nm <- Sys.getenv("DS", "Zanol2014")
phy <- readRDS(file.path(bare, paste0(nm, ".phy.rds")))
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))      # the 1271 set (hold 1000)
first <- sub("[*]$", "", L[2])                                # first tree, drop trailing '*'
writeLines(c("tread 'single T0 = tree1 of best set'", paste0(first, ";"), "proc-;"),
           file.path(bare, paste0(nm, ".t0single.tre")))
t <- ReadTntTree(file.path(bare, paste0(nm, ".t0single.tre")))
if (inherits(t, "multiPhylo")) t <- t[[1]]
t <- RootTree(t, t$tip.label[1])
cat(sprintf("%s single-tree T0 score (TreeLength) = %.0f tips=%d\n",
            nm, TreeLength(t, phy), length(t$tip.label)))
