files <- Sys.glob("/nobackup/pjjg18/reeval/scaffold_out/cell_*.csv")
d <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
d$level <- factor(d$level, levels = c("RAS0","k5","k10","k20","k25","k26","FULL"))
cat(sprintf("cells=%d\n\n", nrow(d)))
cat("SCAFFOLD-SUFFICIENCY: hit-rate to <=1944 (and <=1943) by #deep splits fixed\n")
cat(sprintf("%-6s %4s %8s %8s %8s %8s %7s\n","level","n","hit1944","hit1943","medBest","minBest","held"))
for (lv in levels(d$level)) { x <- d[d$level==lv,]; if(!nrow(x)) next
  cat(sprintf("%-6s %4d %8.2f %8.2f %8.0f %8.0f %7.2f\n", lv, nrow(x),
    mean(x$hit1944), mean(x$hit1943), median(x$best), min(x$best), mean(x$held, na.rm=TRUE))) }
cat("\nper-cell best by level:\n")
for (lv in levels(d$level)) { x <- d[d$level==lv,]; if(!nrow(x)) next
  cat(sprintf("  %-5s: %s\n", lv, paste(sort(x$best), collapse=" "))) }
