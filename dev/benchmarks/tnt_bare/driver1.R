source("dev/benchmarks/tnt_bare/harness.R")

# Sanity: read T0 fresh and report TNT's own score for it.
chk <- run_tnt(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
                 "proc tee.tre;", "score;", "quit;"))
cat("---- start-tree score lines (raw) ----\n")
cat(paste0("  ", trimws(grep("score|Tree|length|1271|1275", chk, ignore.case = TRUE,
           value = TRUE))), sep = "\n")

cat("\n\n======== EXPERIMENT BATCH 1 ========\n")
# Bare-bones: strip global TBR, equal-acceptance, fuse, drift
bare <- run_config("noglobal noequals nofuse godrift 9999", rounds = 12, label = "BARE")
print_config(bare)

# TNT default (no settings changed) for reference
def  <- run_config("", rounds = 12, label = "DEFAULT")
print_config(def)
