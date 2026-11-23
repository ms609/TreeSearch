# Compare arg counts between TreeSearch-init.c and RcppExports.cpp

# Parse TreeSearch-init.c
init_lines <- readLines("src/TreeSearch-init.c")
init_pattern <- '[{]"(_TreeSearch_\\w+)".*,\\s*(\\d+)[}]'
init_matches <- regmatches(init_lines, regexec(init_pattern, init_lines))
init_matches <- init_matches[lengths(init_matches) > 0]
init_df <- data.frame(
  name = vapply(init_matches, `[`, "", 2),
  init_args = as.integer(vapply(init_matches, `[`, "", 3)),
  stringsAsFactors = FALSE
)

# Parse RcppExports.cpp
export_lines <- readLines("src/RcppExports.cpp")
export_pattern <- "RcppExport SEXP (_TreeSearch_\\w+)[(]([^)]*)[)]"
export_matches <- regmatches(export_lines, regexec(export_pattern, export_lines))
export_matches <- export_matches[lengths(export_matches) > 0]
export_df <- data.frame(
  name = vapply(export_matches, `[`, "", 2),
  export_args = vapply(export_matches, function(m) {
    params <- trimws(m[3])
    if (nchar(params) == 0) return(0L)
    length(strsplit(params, ",")[[1]])
  }, integer(1)),
  stringsAsFactors = FALSE
)

cat("init.c entries:", nrow(init_df), "\n")
cat("RcppExports.cpp entries:", nrow(export_df), "\n\n")

# Merge and compare
merged <- merge(init_df, export_df, by = "name", all = TRUE)

# Mismatches in shared entries
mis <- merged[!is.na(merged$init_args) & !is.na(merged$export_args) &
              merged$init_args != merged$export_args, ]
if (nrow(mis) > 0) {
  cat("ARG COUNT MISMATCHES:\n")
  print(mis, row.names = FALSE)
} else {
  cat("All shared entries: arg counts match.\n")
}
cat("\n")

# In init.c but not RcppExports.cpp
manual <- merged[is.na(merged$export_args), ]
if (nrow(manual) > 0) {
  cat("Manual entries (init.c only, not in RcppExports.cpp):", nrow(manual), "\n")
  print(manual[, c("name", "init_args")], row.names = FALSE)
}
cat("\n")

# In RcppExports.cpp but missing from init.c
missing_reg <- merged[is.na(merged$init_args), ]
if (nrow(missing_reg) > 0) {
  cat("MISSING from init.c (in RcppExports.cpp but not registered):\n")
  print(missing_reg[, c("name", "export_args")], row.names = FALSE)
}
