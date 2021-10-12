# Suppress "NOTE: Nothing imported from Rdpack":
#' @importFrom Rdpack reprompt
.onUnload <- function (libpath) {
  library.dynam.unload("TreeSearch", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you cleared GitHub issues for this release milestone?",
    "Have you checked the Vignettes for sanity?"
  )
}


# Additional checks:
#
# codemetar::write_codemeta()
# 
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
#
# devtools::check_win_devel(); rhub::check_for_cran()
#
#
# tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# tools::checkRdaFiles('R') - set optimal compression in `data-raw`
