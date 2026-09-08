# tbr_gate.R -- pre-flight gate #1: TS->TNT tree round-trip length identity.
source("dev/benchmarks/tbr_shared_start_lib.R")

d  <- prepareDataset("Zanol2014")
t0 <- ape::read.tree(file.path(T0_DIR, "Zanol2014.tre"))
cat("TreeLength(T0) =", TreeLength(t0, d$phy), "\n\n")

script <- c("mxram 1024;", "taxname=;", "proc data.tnt;",
            paste0("tread ", ToTntTree(t0), ";"),
            "length;", "quit;")
out <- RunTnt(d$phy, script, tag = "gate")
cat("--- RAW TNT OUTPUT ---\n")
cat(out, sep = "\n")
