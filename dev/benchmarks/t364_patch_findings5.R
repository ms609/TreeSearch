#!/usr/bin/env Rscript
# Fifth surgical patch: patch_findings2.R's prevalence sentence also carried the
# 7-matrix LOCAL pilot rates.  Replace with the 21-matrix corpus figures, which
# additionally expose a 0.0 floor on the common geometry -- the useful number for
# anyone estimating exposure.
f <- Sys.getenv("FINDINGS", "C:/Users/pjjg18/GitHub/TreeSearch/dev/red-team/findings.md")
raw <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
old <- "gives **6-11.5% complement-rooted** (median 8.5%, reproducing this row's 9.5%); tip 0 in the small side gives **71-81%**."
new <- "gives **median 6.0%, range [0.0-11.5]** complement-rooted -- this row's 9.5% sits inside that range, and the **0.0 floor means some matrices never trigger it at all**; tip 0 in the small side gives **median 81.2%, range [71.0-94.0]** (21 matrices x 200 identical addition orders per arm)."
n <- gregexpr(old, raw, fixed = TRUE)[[1]]
if (identical(n[1], -1L)) stop("PATTERN NOT FOUND")
if (length(n) != 1L) stop("matched ", length(n), " times, expected 1")
raw <- sub(old, new, raw, fixed = TRUE)
con <- file(f, open = "wb"); writeLines(raw, con, sep = "\n", useBytes = TRUE); close(con)
cat("patch 5 applied\n")
