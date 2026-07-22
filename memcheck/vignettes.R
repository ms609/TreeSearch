# Run with:
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/vignettes.R
# Package must be installed first.
pkgdown::build_articles(preview = FALSE)
