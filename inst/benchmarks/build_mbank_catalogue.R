#!/usr/bin/env Rscript
# Build a catalogue of MorphoBank matrices from the neotrans corpus.
#
# Scans neotrans/inst/matrices/*.nex, attempts to parse each as phyDat,
# and records metadata (ntax, nchar, patterns, missing%, inapplicable%).
#
# Output: inst/benchmarks/mbank_catalogue.csv
#
# Run from the TreeSearch source root:
#   Rscript inst/benchmarks/build_mbank_catalogue.R
#
# Or from inst/benchmarks/:
#   Rscript build_mbank_catalogue.R

library(TreeTools)

# --- Path resolution ---
find_neotrans_dir <- function() {
  candidates <- c(
    file.path(getwd(), "..", "neotrans", "inst", "matrices"),
    file.path(getwd(), "..", "..", "neotrans", "inst", "matrices"),
    file.path(dirname(getwd()), "neotrans", "inst", "matrices")
  )
  for (d in candidates) {
    d <- normalizePath(d, mustWork = FALSE)
    if (dir.exists(d)) return(d)
  }
  stop("Cannot find neotrans/inst/matrices/. ",
       "Run from TreeSearch source root or inst/benchmarks/.")
}

find_output_dir <- function() {
  candidates <- c(
    file.path(getwd(), "inst", "benchmarks"),
    getwd()
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "bench_datasets.R"))) return(d)
  }
  # Fall back to inst/benchmarks if it exists
  d <- file.path(getwd(), "inst", "benchmarks")
  if (dir.exists(d)) return(d)
  stop("Cannot find inst/benchmarks/ directory.")
}

neotrans_dir <- find_neotrans_dir()
output_dir <- find_output_dir()

cat("Neotrans matrices dir:", neotrans_dir, "\n")
cat("Output dir:", output_dir, "\n")

# --- Find all .nex files ---
nex_files <- list.files(neotrans_dir, pattern = "\\.nex$",
                        full.names = TRUE, recursive = FALSE)
cat("Found", length(nex_files), ".nex files\n")

# --- Parse each file and collect metadata ---
characterize_phyDat <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  lvls <- at$levels
  n_taxa <- length(dataset)
  n_patterns <- length(at$weight)
  n_chars <- sum(at$weight)
  n_states <- ncol(contrast)

  inapp_idx <- which(lvls == "-")
  n_app_states <- n_states - length(inapp_idx)

  td <- matrix(unlist(dataset, use.names = FALSE),
               nrow = n_taxa, byrow = TRUE)
  total_cells <- n_taxa * n_patterns

  n_inapp <- 0L
  n_missing <- 0L
  has_inapp <- length(inapp_idx) > 0
  for (i in seq_len(nrow(contrast))) {
    is_inapp <- has_inapp && contrast[i, inapp_idx] > 0.5
    cols_check <- setdiff(seq_len(n_states), inapp_idx)
    is_all <- length(cols_check) > 0 && all(contrast[i, cols_check] > 0.5)
    count <- sum(td == i)
    if (is_inapp && !is_all) n_inapp <- n_inapp + count
    if (is_all) n_missing <- n_missing + count
  }

  list(
    ntax = n_taxa,
    nchar = n_chars,
    n_patterns = n_patterns,
    n_states = n_app_states,
    pct_missing = round(100 * n_missing / total_cells, 1),
    pct_inapp = round(100 * n_inapp / total_cells, 1)
  )
}

results <- vector("list", length(nex_files))

for (i in seq_along(nex_files)) {
  f <- nex_files[i]
  bname <- basename(f)

  # Extract project ID and matrix index
  if (grepl("^project", bname, ignore.case = TRUE)) {
    proj_num <- as.integer(sub("^project(\\d+).*", "\\1", bname,
                                ignore.case = TRUE))
    # Multi-matrix index: "project1037 (2).nex" -> 2
    if (grepl("\\(\\d+\\)", bname)) {
      mat_idx <- as.integer(sub(".*\\((\\d+)\\).*", "\\1", bname))
    } else {
      mat_idx <- NA_integer_
    }
    source_type <- "morphobank"
  } else if (grepl("^syab", bname, ignore.case = TRUE)) {
    proj_num <- NA_integer_
    mat_idx <- NA_integer_
    source_type <- "syab"
  } else {
    proj_num <- NA_integer_
    mat_idx <- NA_integer_
    source_type <- "other"
  }

  # Unique key for this matrix
  key <- sub("\\.nex$", "", bname, ignore.case = TRUE)
  key <- gsub(" ", "_", key)

  # Assign split
  if (!is.na(proj_num) && proj_num %% 5 == 0) {
    split <- "validation"
  } else {
    split <- "training"
  }

  # Try to parse
  row <- list(
    key = key,
    filename = bname,
    project_id = proj_num,
    matrix_idx = mat_idx,
    source_type = source_type,
    split = split,
    ntax = NA_integer_,
    nchar = NA_integer_,
    n_patterns = NA_integer_,
    n_states = NA_integer_,
    pct_missing = NA_real_,
    pct_inapp = NA_real_,
    parse_ok = FALSE,
    error_message = ""
  )

  tryCatch({
    pd <- ReadAsPhyDat(f)
    chars <- characterize_phyDat(pd)
    row$ntax <- chars$ntax
    row$nchar <- chars$nchar
    row$n_patterns <- chars$n_patterns
    row$n_states <- chars$n_states
    row$pct_missing <- chars$pct_missing
    row$pct_inapp <- chars$pct_inapp
    row$parse_ok <- TRUE
  }, error = function(e) {
    row$error_message <<- conditionMessage(e)
  }, warning = function(w) {
    # Warnings during parsing are common (e.g. "Duplicate taxon names")
    # Try to continue
    tryCatch({
      pd <- suppressWarnings(ReadAsPhyDat(f))
      chars <- characterize_phyDat(pd)
      row$ntax <<- chars$ntax
      row$nchar <<- chars$nchar
      row$n_patterns <<- chars$n_patterns
      row$n_states <<- chars$n_states
      row$pct_missing <<- chars$pct_missing
      row$pct_inapp <<- chars$pct_inapp
      row$parse_ok <<- TRUE
      row$error_message <<- paste("WARNING:", conditionMessage(w))
    }, error = function(e2) {
      row$error_message <<- paste("WARNING:", conditionMessage(w),
                                   "; ERROR:", conditionMessage(e2))
    })
  })

  results[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)

  if (i %% 50 == 0 || i == length(nex_files)) {
    cat(sprintf("  [%d/%d] %s\n", i, length(nex_files), bname))
  }
}

catalogue <- do.call(rbind, results)

# --- Dedup: flag near-duplicate multi-file matrices ---
# Multi-file projects (e.g. "project1037 (1).nex", "project1037 (2).nex") often
# contain the same character data with minor taxon-sampling differences. We flag
# redundant copies so the benchmark loader can exclude them by default.
#
# Method: for each project with multiple usable files, load all matrices,
# compute pairwise character identity on shared taxa, and greedily keep the
# largest (most taxa) representative from each cluster of >=95% identical pairs.

usable_mask <- catalogue$parse_ok & !is.na(catalogue$ntax) & catalogue$ntax >= 20
catalogue$dedup_drop <- FALSE

usable_multi <- catalogue[usable_mask & !is.na(catalogue$matrix_idx), ]
if (nrow(usable_multi) > 0) {
  usable_multi$project <- sub("_\\(\\d+\\)$", "", usable_multi$key)
  proj_counts <- table(usable_multi$project)
  multi_projects <- names(proj_counts[proj_counts >= 2])

  cat(sprintf("\nDedup: checking %d multi-file projects (%d matrices)...\n",
              length(multi_projects),
              sum(usable_multi$project %in% multi_projects)))

  drop_keys <- character(0)

  for (proj in multi_projects) {
    rows <- usable_multi[usable_multi$project == proj, ]
    keys <- rows$key
    mats <- list()
    for (j in seq_len(nrow(rows))) {
      fpath <- file.path(neotrans_dir, rows$filename[j])
      tryCatch({
        mats[[rows$key[j]]] <- suppressWarnings(ReadAsPhyDat(fpath))
      }, error = function(e) NULL)
    }
    if (length(mats) < 2) next

    # Build pairwise character-identity matrix
    mk <- names(mats)
    identity_mat <- matrix(NA_real_, length(mk), length(mk),
                           dimnames = list(mk, mk))
    for (a in seq_len(length(mk) - 1)) {
      for (b in (a + 1):length(mk)) {
        taxa_a <- names(mats[[mk[a]]])
        taxa_b <- names(mats[[mk[b]]])
        common <- intersect(taxa_a, taxa_b)
        # Require >=80% taxon overlap with the smaller matrix
        if (length(common) < 0.8 * min(length(taxa_a), length(taxa_b))) next
        mat_a <- as.matrix(mats[[mk[a]]])[common, , drop = FALSE]
        mat_b <- as.matrix(mats[[mk[b]]])[common, , drop = FALSE]
        if (ncol(mat_a) != ncol(mat_b)) next
        identity_mat[mk[a], mk[b]] <- mean(mat_a == mat_b, na.rm = TRUE)
        identity_mat[mk[b], mk[a]] <- identity_mat[mk[a], mk[b]]
      }
    }

    # Greedy dedup: sort by ntax desc, keep first, drop near-dups
    sorted_keys <- rows$key[order(-rows$ntax, -rows$nchar)]
    kept <- character(0)
    for (k in sorted_keys) {
      is_dup <- FALSE
      for (kk in kept) {
        ci <- identity_mat[k, kk]
        if (!is.na(ci) && ci >= 0.95) { is_dup <- TRUE; break }
      }
      if (is_dup) drop_keys <- c(drop_keys, k)
      else kept <- c(kept, k)
    }
  }

  catalogue$dedup_drop[catalogue$key %in% drop_keys] <- TRUE
  cat(sprintf("Dedup: flagged %d near-duplicate matrices for exclusion.\n",
              length(drop_keys)))
}

# --- Summary ---
cat("\n=== Catalogue Summary ===\n")
cat("Total files scanned:", nrow(catalogue), "\n")
cat("Parse OK:", sum(catalogue$parse_ok), "\n")
cat("Parse failed:", sum(!catalogue$parse_ok), "\n")
cat("\nAfter ntax >= 20 filter:\n")
usable <- catalogue$parse_ok & !is.na(catalogue$ntax) & catalogue$ntax >= 20
cat("  Usable (before dedup):", sum(usable), "\n")
cat("  Dedup dropped:", sum(usable & catalogue$dedup_drop), "\n")
usable_dedup <- usable & !catalogue$dedup_drop
cat("  Usable (after dedup):", sum(usable_dedup), "\n")
cat("  Training:", sum(usable_dedup & catalogue$split == "training"), "\n")
cat("  Validation:", sum(usable_dedup & catalogue$split == "validation"), "\n")

cat("\nSize tiers (after dedup):\n")
usable_cat <- catalogue[usable_dedup, ]
usable_cat$tier <- cut(usable_cat$ntax,
                        breaks = c(0, 30, 60, 120, Inf),
                        labels = c("Small(20-30)", "Medium(31-60)",
                                   "Large(61-120)", "XLarge(121+)"))
print(table(usable_cat$split, usable_cat$tier))

cat("\nParse failures:\n")
if (any(!catalogue$parse_ok)) {
  fails <- catalogue[!catalogue$parse_ok, c("key", "error_message")]
  for (j in seq_len(nrow(fails))) {
    cat(sprintf("  %s: %s\n", fails$key[j],
                substr(fails$error_message[j], 1, 80)))
  }
}

# --- Save ---
out_path <- file.path(output_dir, "mbank_catalogue.csv")
write.csv(catalogue, out_path, row.names = FALSE)
cat("\nCatalogue written to:", out_path, "\n")
