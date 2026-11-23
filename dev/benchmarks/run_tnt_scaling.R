setwd("C:/Users/pjjg18/GitHub/TreeSearch")
source("dev/benchmarks/bench_tnt_settings.R")

# Quick smoke-test: project691 x sect+fuse x seed=1
cat("--- Smoke test: project691 / sect+fuse / seed=1 ---\n")
info <- export_nexus_dataset("project691")
cat(sprintf("Exported: %dt %dc\n", info$ntip, info$nchar))

sc1 <- write_phase1_script("project691.tnt", seed = 1L, timeout_s = 60L)
r1  <- run_tnt(sc1, hard_timeout_s = 90L)
cat(sprintf("Phase1 seed=1: score=%g  wall=%.1fs\n", r1$score, r1$wall_s))

if (!is.na(r1$score)) {
  B <- r1$score
  cfg <- CONFIGS[["sect+fuse"]]
  sc2 <- write_survey_script("project691.tnt", cfg, B, seed = 1L, timeout_s = 60L)
  r2  <- run_tnt(sc2, hard_timeout_s = 90L)
  ttt <- parse_ttt(r2$raw, B)
  reached <- isTRUE(!is.na(r2$score) && r2$score <= B + 1e-6)
  ttb <- if (!is.na(ttt$ttb_s) && ttt$ttb_s > 0) ttt$ttb_s else r2$wall_s
  cat(sprintf("Phase2 sect+fuse: score=%g  reached=%s  TTT=%.1fs\n",
              r2$score, reached, ttb))
}

cat("\n--- Launching full scaling survey ---\n")
results <- tnt_scaling_full()
message("Scaling survey complete. ", nrow(results), " rows.")
