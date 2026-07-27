Rcpp::compileAttributes()
for (f in c("src/RcppExports.cpp", "R/RcppExports.R")) {
  if (file.exists(f)) {
    lines <- readLines(f, warn = FALSE)
    con <- file(f, open = "wb")
    writeLines(lines, con = con, sep = "\n")
    close(con)
  }
}
source("check_init.R")
