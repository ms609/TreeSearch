#!/usr/bin/env Rscript
# Resolve VTune MinGW-DWARF `func@0xADDR` hotspots to demangled C++ names.
#
# VTune 2026's CSV reporter leaves MinGW DWARF symbols unparsed (`func@0xADDR`)
# even when the DLL carries `.debug_info`. This joins the hotspot addresses to
# `nm -C` output by address (nearest function start <= addr) so the flat
# self-time profile is readable. Recurs every profiling round on this toolchain.
#
# Inputs (dump these first):
#   vtune -report hotspots -r <res> -group-by function -format=csv \
#         -csv-delimiter="`t" > hs.tsv          # TAB delimiter (names have commas)
#   nm -C --defined-only <dll> > nm.txt
# Usage: Rscript resolve_syms.R <hs.tsv> <nm.txt> [topN]

args <- commandArgs(trailingOnly = TRUE)
hs_path <- args[1]; nm_path <- args[2]
topN <- if (length(args) >= 3) as.integer(args[3]) else 30L

hs <- read.delim(hs_path, check.names = FALSE, stringsAsFactors = FALSE)
names(hs)[names(hs) == "CPU Time"] <- "cpu"
hs$addr <- suppressWarnings(as.numeric(hs[["Start Address"]]))

nm <- readLines(nm_path, warn = FALSE)
mm <- regmatches(nm, regexec("^([0-9a-fA-F]{8,})[ \t]+(\\S)[ \t]+(.*)$", nm))
ok <- lengths(mm) == 4
syma <- as.numeric(paste0("0x", vapply(mm[ok], `[`, "", 2)))
symn <- vapply(mm[ok], `[`, "", 4)
o <- order(syma); syma <- syma[o]; symn <- symn[o]

resolve <- function(a) {
  if (is.na(a)) return(NA_character_)
  i <- findInterval(a, syma)
  if (i < 1) return(NA_character_)
  # Return: name
  symn[i]
}
hs$resolved <- vapply(hs$addr, resolve, "")
isTS <- hs$Module == "TreeSearch.dll"
hs$name <- ifelse(isTS & !is.na(hs$resolved), hs$resolved,
                  ifelse(grepl("^func@", hs$Function),
                         paste0(hs$Module, "!", hs$Function), hs$Function))

total <- sum(hs$cpu, na.rm = TRUE)
ts_cpu <- sum(hs$cpu[isTS], na.rm = TRUE)
agg <- aggregate(cpu ~ name + Module, hs, sum)
agg$pct <- 100 * agg$cpu / total
agg <- agg[order(-agg$cpu), ]

cat(sprintf("Total CPU time (all modules): %.3f s\n", total))
cat(sprintf("TreeSearch.dll self CPU:      %.3f s (%.1f%% of total)\n\n",
            ts_cpu, 100 * ts_cpu / total))
cat(sprintf("%7s %7s  %s\n", "self_s", "pct", "function [module]"))
for (i in seq_len(min(topN, nrow(agg)))) {
  cat(sprintf("%7.3f %6.1f%%  %s  [%s]\n",
              agg$cpu[i], agg$pct[i], agg$name[i], agg$Module[i]))
}
