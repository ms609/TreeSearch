if (requireNamespace("spelling", quietly = TRUE)) {
  # spell_check_test() warns "Failed to find package source directory" when run
  # from covr's temp install path; that warning is harmless, so catch and
  # suppress it rather than letting it become an error under error-on="warning".
  withCallingHandlers(
    spelling::spell_check_test(vignettes = TRUE, error = TRUE,
                               skip_on_cran = TRUE),
    warning = function(w) {
      if (grepl("find package source", conditionMessage(w), fixed = FALSE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
