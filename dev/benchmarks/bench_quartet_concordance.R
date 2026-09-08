# Benchmark and profile quartet_concordance
# T-298: Profile flat array vs NumericMatrix, int vs double
#
# Usage: Rscript dev/benchmarks/bench_quartet_concordance.R
# Or interactively for profvis output.
#
# Install from source first (tarball build per AGENTS.md).

library(TreeSearch)

set.seed(42)

make_inputs <- function(n_taxa, n_splits, n_chars, n_states = 4) {
  splits <- matrix(sample(c(TRUE, FALSE), n_taxa * n_splits, replace = TRUE),
                   nrow = n_taxa, ncol = n_splits)
  # IntegerMatrix with NAs (~5% missing)
  chars <- matrix(sample(c(0:(n_states - 1), NA_integer_),
                         n_taxa * n_chars, replace = TRUE, prob = c(rep(0.95/n_states, n_states), 0.05)),
                  nrow = n_taxa, ncol = n_chars)
  list(splits = splits, chars = chars)
}

sizes <- list(
  small  = list(n_taxa = 25,  n_splits = 23,  n_chars = 50),
  medium = list(n_taxa = 100, n_splits = 98,  n_chars = 200),
  large  = list(n_taxa = 300, n_splits = 298, n_chars = 400)
)

cat("=== Timing quartet_concordance ===\n")
for (nm in names(sizes)) {
  sz <- sizes[[nm]]
  inp <- make_inputs(sz$n_taxa, sz$n_splits, sz$n_chars)
  t <- system.time(
    for (i in seq_len(10)) TreeSearch:::quartet_concordance(inp$splits, inp$chars)
  )
  cat(sprintf("%-8s (%3d taxa, %3d splits, %3d chars): %.3fs / call\n",
              nm, sz$n_taxa, sz$n_splits, sz$n_chars, t[["elapsed"]] / 10))
}

# --- profvis profile of the medium case ---
if (requireNamespace("profvis", quietly = TRUE)) {
  inp <- make_inputs(100, 98, 200)
  p <- profvis::profvis({
    for (i in seq_len(50)) TreeSearch:::quartet_concordance(inp$splits, inp$chars)
  })
  print(p)
} else {
  message("Install profvis for flame graph: install.packages('profvis')")
}
