# Run with:
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/tests.R
# Package must be installed first.
testthat::test_local()
