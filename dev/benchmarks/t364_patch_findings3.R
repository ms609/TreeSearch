#!/usr/bin/env Rscript
# Third surgical patch: attach the two scope caveats to the T-384 row's magnitude
# figures -- (a) every budget-96 number is small+medium only, (b) the budget-24
# `no0` comparison is budget-bound.  Both change what the figures may be cited for,
# not what they conclude.  Same match-exactly-once discipline.
f <- Sys.getenv("FINDINGS", "C:/Users/pjjg18/GitHub/TreeSearch/dev/red-team/findings.md")
raw <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

old <- "`no0` 1.053, pooled **1.211** (19/8, 16 >10%, p = 0.052). "
new <- paste0(
  "`no0` 1.053, pooled **1.211** (19/8, 16 >10%, p = 0.052). ",
  "**Two scope caveats that bound what these may be cited for.** (i) Every ",
  "`maxReplicates = 96` figure here -- including the 2.531x and the merged arm's ",
  "0.967 below -- is **small+medium tiers ONLY** (14 matrices, 27 pairs, the `in0` ",
  "rows just 13); large/xlarge were not run at that budget, so none of them is a ",
  "corpus-wide number. (ii) The budget-24 `no0` comparison is **budget-bound and ",
  "near-vacuous**: all three arms ran 24/24/24 replicates there (exhaustion ",
  "64/65/67%), so wall is roughly equal by construction; `in0` at budget 24 was ",
  "genuinely convergence-limited (16/24/17.5), and the budget-96 run un-binds `no0` ",
  "(31.5/29 replicates, not exhausted, merged/pre-fix 0.969), which is the ",
  "informative figure for that stratum. ")

n <- gregexpr(old, raw, fixed = TRUE)[[1]]
if (identical(n[1], -1L)) stop("PATTERN NOT FOUND -- did patch_findings2.R run?")
if (length(n) != 1L) stop("matched ", length(n), " times, expected 1")
raw <- sub(old, new, raw, fixed = TRUE)

con <- file(f, open = "wb")            # binary: preserve LF
writeLines(raw, con, sep = "\n", useBytes = TRUE)
close(con)
cat("patch 3 applied\n")
