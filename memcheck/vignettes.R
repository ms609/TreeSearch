# Run with:
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/vignettes.R
# Package must be installed first.
tools::buildVignettes(dir = ".")
