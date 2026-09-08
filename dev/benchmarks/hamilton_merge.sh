#!/bin/bash
# Merge per-cell partial CSVs from the array into one panel CSV.
# Submit with --dependency=afterany:<arrayjob> (afterany so partial failures
# still merge what succeeded). Validate row count == expected grid size.
# DISPATCH-UNTESTED template.
#SBATCH --job-name=ts-merge
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH --time=0:10:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/merge_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/merge_%j.err

module load r/4.5.1
P=/nobackup/$USER/TreeSearch/panel_partials
O=/nobackup/$USER/TreeSearch/panel_results
mkdir -p "$O"
Rscript -e "f<-list.files('$P',pattern='cell_.*csv\$',full.names=TRUE); \
  d<-do.call(rbind,lapply(f,read.csv)); \
  write.csv(d,file.path('$O','panel.csv'),row.names=FALSE); \
  cat(sprintf('%d rows from %d cells -> %s/panel.csv\n', nrow(d), length(f), '$O'))"
