# Chain to the user-level profile.
#
# R reads EITHER ./.Rprofile OR ~/.Rprofile -- never both (see ?Startup).  The
# mere existence of this file therefore silently suppresses the user profile,
# and on CI that profile is where `r-lib/actions/setup-r` writes
# `options(repos = c(RSPM = <Posit binary repo>, CRAN = ...), Ncpus,
# HTTPUserAgent)`.  Without this chaining, every CI job resolves 100% of its
# dependencies from source off the default CRAN mirror -- ~10 min of compiling
# per leg -- because Posit Package Manager never enters `getOption("repos")`
# and the User-Agent it needs to serve binaries is never set.
local({
  user <- path.expand("~/.Rprofile")
  if (file.exists(user) &&
      normalizePath(user, mustWork = FALSE) !=
        normalizePath("./.Rprofile", mustWork = FALSE)) {
    try(sys.source(user, envir = globalenv()))
  }
})

# Patch Rcpp::compileAttributes() to write LF line endings, preventing CRLF
# churn in RcppExports files on Windows.  Rcpp always deletes and rewrites
# R/RcppExports.R from scratch (so it is always CRLF on Windows without this
# patch), while src/RcppExports.cpp uses writeIfChanged (so it is only rewritten
# when content changes).  Patching here intercepts every caller — devtools,
# pkgbuild, RStudio auto-build — not just our manual wrapper script.
local({
  .normalise_lf <- function(pkgdir) {
    for (f in file.path(pkgdir, c("src/RcppExports.cpp", "R/RcppExports.R"))) {
      if (!file.exists(f)) next
      lines <- readLines(f, warn = FALSE)
      con <- file(f, open = "wb")
      writeLines(lines, con = con, sep = "\n")
      close(con)
    }
  }

  .patch_rcpp <- function(...) {
    env <- tryCatch(asNamespace("Rcpp"), error = function(e) NULL)
    if (is.null(env) || !exists("compileAttributes", envir = env, inherits = FALSE))
      return()
    orig <- get("compileAttributes", envir = env, inherits = FALSE)
    if (isTRUE(attr(orig, ".lf_patched"))) return()
    patched <- function(pkgdir = ".", verbose = getOption("verbose")) {
      result <- orig(pkgdir = pkgdir, verbose = verbose)
      .normalise_lf(pkgdir)
      invisible(result)
    }
    attr(patched, ".lf_patched") <- TRUE
    try(suppressWarnings(unlockBinding("compileAttributes", env)), silent = TRUE)
    assign("compileAttributes", patched, envir = env)
    try(suppressWarnings(lockBinding("compileAttributes", env)), silent = TRUE)
  }

  setHook(packageEvent("Rcpp", "onLoad"), .patch_rcpp)
  if (isNamespaceLoaded("Rcpp")) .patch_rcpp()
})
