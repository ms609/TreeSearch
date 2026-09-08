# SLURM array cell for the basin-diversity harness.
# Runs ONE grid cell (one INDEP seed, or one PROD variant/seed) and writes a
# partial .rds.  Mirrors dev/benchmarks/bench_cell.R.
#
# Local test:   BD_MODE=cell Rscript dev/benchmarks/basin_cell.R 0
# Array task:   SLURM_ARRAY_TASK_ID supplies the 0-based cell index.
#
# Grid size = BD_INDEP_K + 2 * BD_PROD_SEEDS  (set --array=0-(N-1) to match).
Sys.setenv(BD_MODE = "cell")
source(file.path("dev", "benchmarks", "basin_diversity.R"), local = FALSE)
.bdMain()  # source() suppresses the auto-run guard, so invoke explicitly
