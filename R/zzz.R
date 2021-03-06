.onUnload <- function (libpath) {
  library.dynam.unload("TreeSearch", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you cleared GitHub issues for this release milestone?",
    "Have you checked the Vignettes for sanity?",
    "Have you checked pkgdown::build_reference_index()?",
    "Have you updated the version number in NEWS & DESCRIPTION?"
  )
}



# Additional steps:
#
# Propogate changes in README.md to R/TreeTools-package.R


# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
#
# devtools::check_win_devel(); rhub::check_for_cran()
# rhub::check_with_valgrind() # runs the build and check on Linux, in valgrind to find memory leaks and pointer errors.
# rhub::check_with_sanitizers() # runs all package package tests, examples and vignettes with Address Sanitizer and Undefined Behavior Sanitizer.
#
# revdepcheck::revdep_check()
#
# codemetar::write_codemeta()
#
# tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# tools::checkRdaFiles('R') - set optimal compression in `data-raw`
