# NNI survey: measure batch-NNI feasibility
#
# For each dataset, builds Wagner trees and surveys all NNI candidates to
# count how many moves improve the score. This measures the theoretical
# payoff of batch/simultaneous NNI at different search stages.
#
# Usage: Rscript dev/benchmarks/bench_nni_survey.R <lib_path>

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else stop("Usage: Rscript bench_nni_survey.R <lib_path>")
.libPaths(c(lib_path, .libPaths()))

pkg_name <- basename(lib_path)
agent_letter <- sub(".*-", "", pkg_name)
renamed <- paste0("TreeSearch.", agent_letter)
library(renamed, character.only = TRUE)
if (is.null(.Internal(getRegisteredNamespace("TreeSearch"))))
  .Internal(registerNamespace("TreeSearch", asNamespace(renamed)))

library(TreeTools)

prepare_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(dataset)
  )
}

build_wagner <- function(ds, seed) {
  set.seed(seed)
  TreeSearch:::ts_wagner_tree(ds$contrast, ds$tip_data, ds$weight, ds$levels)
}

run_survey <- function(edge_mat, ds) {
  TreeSearch:::ts_nni_survey(
    edge_mat, ds$contrast, ds$tip_data, ds$weight, ds$levels
  )
}

run_nni <- function(edge_mat, ds, maxHits = 20L) {
  TreeSearch:::ts_nni_search(
    edge_mat, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = maxHits
  )
}

analyze_survey <- function(survey) {
  deltas <- survey$delta
  n_candidates <- length(deltas)
  n_improving <- sum(deltas < 0)
  n_equal <- sum(deltas == 0)

  edge_ids <- survey$edge
  best_per_edge <- tapply(deltas, edge_ids, min)
  n_edges_improving <- sum(best_per_edge < 0)

  total_improvement <- -sum(deltas[deltas < 0])
  best_improvement <- if (n_improving > 0) -min(deltas) else 0L

  data.frame(
    base_score = survey$base_score,
    n_edges = survey$n_edges,
    n_candidates = n_candidates,
    n_improving = n_improving,
    n_equal = n_equal,
    n_edges_improving = n_edges_improving,
    total_improvement = total_improvement,
    best_single_improvement = best_improvement,
    pct_edges_improving = round(100 * n_edges_improving / survey$n_edges, 1)
  )
}

# All standard Fitch datasets (no inapplicable-dominant ones)
DATASETS <- c(
  "Vinther2008",    # 23 tips
  "Griswold1999",   # 43 tips
  "Eklund2004",     # 54 tips
  "Agnarsson2004",  # 62 tips
  "Zhu2013",        # 75 tips
  "Giles2015",      # 78 tips
  "Dikow2009"       # 88 tips
)

SEEDS <- c(1742L, 5281L, 8093L, 3647L, 9210L)

cat("=== NNI Survey: Batch-NNI Feasibility ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

all_wagner <- list()
all_converged <- list()

for (nm in DATASETS) {
  ds_raw <- TreeSearch::inapplicable.phyData[[nm]]
  if (is.null(ds_raw)) { cat("SKIP:", nm, "\n"); next }
  ds <- prepare_ts_data(ds_raw)
  n_tips <- ds$n_taxa

  cat(sprintf("\n--- %s (%d tips, %d edges) ---\n", nm, n_tips, n_tips - 2L))

  for (seed in SEEDS) {
    # Stage 1: Wagner tree
    wagner <- build_wagner(ds, seed)
    survey_w <- run_survey(wagner$edge, ds)
    info_w <- analyze_survey(survey_w)
    info_w$dataset <- nm
    info_w$n_tips <- n_tips
    info_w$seed <- seed
    info_w$stage <- "wagner"

    cat(sprintf("  seed=%d  Wagner: score=%d, %d/%d edges improving (total delta=%d, best=%d)\n",
                seed, as.integer(info_w$base_score),
                info_w$n_edges_improving, info_w$n_edges,
                info_w$total_improvement, info_w$best_single_improvement))

    all_wagner <- c(all_wagner, list(info_w))

    # Stage 2: After NNI convergence (maxHits=20, full plateau search)
    nni_result <- run_nni(wagner$edge, ds, maxHits = 20L)
    survey_c <- run_survey(nni_result$edge, ds)
    info_c <- analyze_survey(survey_c)
    info_c$dataset <- nm
    info_c$n_tips <- n_tips
    info_c$seed <- seed
    info_c$stage <- "nni_converged"
    info_c$nni_moves <- nni_result$n_moves
    info_c$nni_iterations <- nni_result$n_iterations

    cat(sprintf("           NNI converged: score=%d (%d moves, %d iter), %d improving edges\n",
                as.integer(info_c$base_score),
                nni_result$n_moves, nni_result$n_iterations,
                info_c$n_edges_improving))

    all_converged <- c(all_converged, list(info_c))
  }
}

wagner_df <- do.call(rbind, all_wagner)
converged_df <- do.call(rbind, all_converged)

cat("\n\n========================================\n")
cat("=== SUMMARY: Wagner Tree Surveys ===\n")
cat("========================================\n\n")

for (nm in unique(wagner_df$dataset)) {
  sub <- wagner_df[wagner_df$dataset == nm, ]
  csub <- converged_df[converged_df$dataset == nm, ]
  cat(sprintf("%s (%d tips, %d NNI edges):\n", nm, sub$n_tips[1], sub$n_edges[1]))
  cat(sprintf("  Wagner scores:    %d-%d (median %d)\n",
              min(as.integer(sub$base_score)),
              max(as.integer(sub$base_score)),
              as.integer(median(sub$base_score))))
  cat(sprintf("  Improving edges:  %d-%d (median %.0f, %.0f%% of edges)\n",
              min(sub$n_edges_improving), max(sub$n_edges_improving),
              median(sub$n_edges_improving),
              median(sub$pct_edges_improving)))
  cat(sprintf("  Total delta:      %d-%d steps (median %d)\n",
              min(sub$total_improvement), max(sub$total_improvement),
              as.integer(median(sub$total_improvement))))
  cat(sprintf("  Best single move: %d-%d steps\n",
              min(sub$best_single_improvement),
              max(sub$best_single_improvement)))
  cat(sprintf("  NNI-converged:    score %d-%d (%d-%d moves)\n\n",
              min(as.integer(csub$base_score)),
              max(as.integer(csub$base_score)),
              min(csub$nni_moves), max(csub$nni_moves)))
}

cat("\n=== Key Finding: Batch Size (improving edges on Wagner trees) ===\n")
cat(sprintf("%-15s %5s %10s %10s %10s %10s\n",
            "Dataset", "Tips", "Med.Batch", "Max.Batch", "%Edges", "Med.Delta"))
for (nm in unique(wagner_df$dataset)) {
  sub <- wagner_df[wagner_df$dataset == nm, ]
  cat(sprintf("%-15s %5d %10.0f %10d %9.0f%% %10d\n",
              nm, sub$n_tips[1],
              median(sub$n_edges_improving),
              max(sub$n_edges_improving),
              median(sub$pct_edges_improving),
              as.integer(median(sub$total_improvement))))
}
