# Run with:
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/all.R
# Package must be installed first.
testthat::test_local()
devtools::run_examples()
tools::buildVignettes(dir = ".")
