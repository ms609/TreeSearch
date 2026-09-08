# Validation gate: C++ LS scorer vs phangorn::nnls.tree / designTree.
suppressMessages({
  library(ape)
  pkgload::load_all(".", quiet = TRUE)
})
stopifnot(requireNamespace("phangorn", quietly = TRUE))

fit_cpp <- function(tree, D, method) {  # method 0=OLS, 1=NNLS
  D <- D[tree$tip.label, tree$tip.label]
  TreeSearch:::ts_ls_fit(tree$edge, D, NULL, as.integer(method))
}

# cophenetic of the rooted tree carrying fitted lengths (root edges are
# returned as (full, 0) so the rooted cophenetic already equals the unrooted).
cophen_cpp <- function(tree, fit) {
  t2 <- tree
  t2$edge.length <- fit$edge_length
  cophenetic(t2)
}

report <- function(label, ok, extra = "") {
  cat(sprintf("%-46s %s %s\n", label, if (ok) "PASS" else "**FAIL**", extra))
  invisible(ok)
}

set.seed(42)
all_ok <- TRUE

for (n in c(5, 6, 8, 12)) {
  # Raw rtree keeps standard ape numbering (root = n+1) which both
  # build_topology_tree and phangorn accept; TreeTools::Preorder sets an
  # `order` attribute that trips phangorn's reorder().
  tree <- rtree(n, br = function(k) runif(k, 0.1, 2))
  labs <- tree$tip.label

  ## ---- additive matrix from this very tree ----
  D <- cophenetic(tree)[labs, labs]

  for (meth in c(0, 1)) {
    mname <- if (meth == 0) "OLS" else "NNLS"
    fit <- fit_cpp(tree, D, meth)
    coph <- cophen_cpp(tree, fit)[labs, labs]
    add_ok <- max(abs(coph - D)) < 1e-6 && fit$rss < 1e-8
    all_ok <- report(sprintf("n=%2d additive  %-4s recovers D, RSS~0", n, mname),
                     add_ok, sprintf("(rss=%.2e, maxerr=%.2e)", fit$rss, max(abs(coph - D)))) && all_ok
  }

  ## ---- phangorn NNLS oracle on the SAME additive matrix ----
  ph <- phangorn::nnls.tree(D, tree, method = "unrooted")
  rss_ph <- attr(ph, "RSS"); if (is.null(rss_ph)) rss_ph <- 0
  coph_ph <- cophenetic(ph)[labs, labs]
  fit_nnls <- fit_cpp(tree, D, 1)
  coph_me <- cophen_cpp(tree, fit_nnls)[labs, labs]
  match_ok <- max(abs(coph_me - coph_ph)) < 1e-6
  all_ok <- report(sprintf("n=%2d additive  NNLS fitted-dist == phangorn", n),
                   match_ok, sprintf("(maxdiff=%.2e)", max(abs(coph_me - coph_ph)))) && all_ok

  ## ---- non-additive matrix ----
  Dn <- as.matrix(as.dist(matrix(0, n, n)))
  rv <- runif(n * (n - 1) / 2, 0.5, 3)
  Dn[lower.tri(Dn)] <- rv; Dn <- Dn + t(Dn)
  dimnames(Dn) <- list(labs, labs)

  # phangorn NNLS oracle
  ph_n <- phangorn::nnls.tree(Dn, tree, method = "unrooted")
  rss_ph_n <- attr(ph_n, "RSS")
  if (is.null(rss_ph_n)) {                          # quadprog branch: recompute
    rss_ph_n <- sum((Dn[labs, labs][lower.tri(Dn)] -
                     cophenetic(ph_n)[labs, labs][lower.tri(Dn)])^2)
  }
  fit_n <- fit_cpp(tree, Dn, 1)
  coph_n <- cophen_cpp(tree, fit_n)[labs, labs]
  rss_me_n <- sum((Dn[lower.tri(Dn)] - coph_n[lower.tri(coph_n)])^2)
  nnls_ok <- abs(rss_me_n - rss_ph_n) < 1e-5 * max(1, rss_ph_n) &&
             max(abs(coph_n - cophenetic(ph_n)[labs, labs])) < 1e-5
  all_ok <- report(sprintf("n=%2d non-add   NNLS RSS == phangorn", n),
                   nnls_ok, sprintf("(me=%.5f ph=%.5f)", rss_me_n, rss_ph_n)) && all_ok

  # OLS vs direct normal-equation solve on the unrooted design
  ut <- unroot(tree)
  X <- as.matrix(phangorn::designTree(ut))
  dm <- Dn[ut$tip.label, ut$tip.label]
  y <- dm[lower.tri(dm)]
  beta <- solve(crossprod(X), crossprod(X, y))
  rss_direct <- sum((y - X %*% beta)^2)
  fit_o <- fit_cpp(tree, Dn, 0)
  ols_ok <- abs(fit_o$rss - rss_direct) < 1e-5 * max(1, rss_direct)
  all_ok <- report(sprintf("n=%2d non-add   OLS RSS == direct solve", n),
                   ols_ok, sprintf("(me=%.5f direct=%.5f)", fit_o$rss, rss_direct)) && all_ok

  # OLS RSS must be <= NNLS RSS (unconstrained)
  mono_ok <- fit_o$rss <= rss_me_n + 1e-8
  all_ok <- report(sprintf("n=%2d non-add   OLS RSS <= NNLS RSS", n), mono_ok) && all_ok
  cat("\n")
}

cat(if (all_ok) "ALL VALIDATION CHECKS PASSED\n" else "SOME CHECKS FAILED\n")
