---
name: hamilton-hpc
description: Run benchmarks, builds, and batch R jobs on the Hamilton HPC cluster (Durham University). Covers SSH connection via R's ssh package, SLURM job submission, module loading, TreeSearch installation, and result retrieval. Use when the user mentions Hamilton, HPC, remote benchmarks, SLURM jobs, or running anything on the cluster.
---

# Hamilton HPC — Remote Job Execution

## Connection

Hamilton is accessed via the R `ssh` package (the local machine has no
CLI SSH config for Hamilton). The session object is typically stored in
the R global environment as `session`.

Connection credentials may be set in `.Renviron`:

| Variable | Purpose |
|----------|---------|
| `sshLogin` | `user@host` (e.g. `pjjg18@hamilton8.dur.ac.uk`) |
| `sshKey` | Path to private key (relative to `.Renviron` location) |
| `sshPass` | Passphrase hint / lookup key |

Try environment variables first, falling back to hardcoded defaults:

```r
library(ssh)
login <- Sys.getenv("sshLogin", "pjjg18@hamilton8.dur.ac.uk")
key   <- Sys.getenv("sshKey", "")
session <- if (nzchar(key)) {
  ssh_connect(login, keyfile = key)
} else {
  ssh_connect(login)
}
ssh_info(session)  # verify connection
```

If the session has gone stale (`ssh_info()` errors), reconnect using
the same pattern.

### Running commands

`ssh_exec_internal()` runs in a minimal shell **without** module
environments. Commands that need `Rscript`, `gcc`, etc. will fail with
status 127. For anything that requires loaded modules, use SLURM job
submission instead.

```r
# Simple commands work fine
out <- ssh_exec_internal(session, "ls /nobackup/pjjg18/ts-bench/")
cat(rawToChar(out$stdout))

# This will FAIL (status 127) — Rscript not in PATH:
# ssh_exec_internal(session, "Rscript -e '1+1'")
```

### File transfer

```r
# Upload
scp_upload(session, "local_script.R", "/nobackup/pjjg18/ts-bench/local_script.R")

# Download
scp_download(session, "/nobackup/pjjg18/ts-bench/results/output.csv", "local_output.csv")
```

## Working directory layout

All TreeSearch benchmark work lives under `/nobackup/pjjg18/ts-bench/`:

| Path | Purpose |
|------|---------|
| `TreeSearch/` | Git clone of the repo (used for builds) |
| `lib-baseline/` | R library with deps + baseline TreeSearch build |
| `lib-optimized/` | R library with alternative TreeSearch build (A/B testing) |
| `results/` | SLURM logs and CSV output files |
| `data/` | Additional datasets |
| `*.R` | Benchmark and build R scripts |
| `*.sh` | SLURM job scripts |

## Environment modules

Hamilton uses `module load` for toolchains. Every SLURM script must
include:

```bash
module load r/4.5.1
module load gcc/14.2
```

Pin single-threaded execution to avoid BLAS/OpenMP contention with
TreeSearch's own threading:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

## SLURM job submission

### Writing a job script

```bash
#!/bin/bash
#SBATCH --job-name=ts-bench
#SBATCH --output=/nobackup/pjjg18/ts-bench/results/my_job.log
#SBATCH --error=/nobackup/pjjg18/ts-bench/results/my_job.err
#SBATCH -n 1
#SBATCH --time=0:15:00
#SBATCH --mem=4000M
#SBATCH -p shared

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript /nobackup/pjjg18/ts-bench/my_script.R
```

Typical resource requests:
- **Benchmarks (30s budget, 5 seeds):** `--time=0:05:00 --mem=4000M -n 1`
- **Benchmarks (120s budget, 5 seeds):** `--time=0:15:00 --mem=4000M -n 1`
- **Package builds:** `--time=0:30:00 --mem=4000M -n 4`
- **Dependency installs:** `--time=0:30:00 --mem=4000M -n 4`

### Submitting and monitoring from R

```r
# Upload script then submit
scp_upload(session, "my_job.sh", "/nobackup/pjjg18/ts-bench/my_job.sh")
out <- ssh_exec_internal(session, "cd /nobackup/pjjg18/ts-bench && sbatch my_job.sh")
cat(rawToChar(out$stdout))  # "Submitted batch job 12345678"

# Parameterized submission (pass variables via --export)
cmd <- sprintf(
  "cd /nobackup/pjjg18/ts-bench && sbatch --export=VAR1=%d,VAR2=%d -o results/run_%d.log -e results/run_%d.err my_job.sh",
  val1, val2, val1, val1
)
out <- ssh_exec_internal(session, cmd)

# Monitor queue
out <- ssh_exec_internal(session, "squeue -u pjjg18 --format='%i %j %T %M'")
cat(rawToChar(out$stdout))

# Read results after completion
out <- ssh_exec_internal(session, "cat /nobackup/pjjg18/ts-bench/results/my_job.log")
cat(rawToChar(out$stdout))
```

### Polling for completion

There is no built-in callback. Poll `squeue` or check for output files:

```r
# Poll until a specific job finishes
poll_job <- function(session, job_id, interval = 30) {
  repeat {
    out <- ssh_exec_internal(session, sprintf("squeue -j %s -h 2>/dev/null", job_id))
    if (nchar(rawToChar(out$stdout)) == 0) return(TRUE)
    Sys.sleep(interval)
  }
}
```

## Building TreeSearch on Hamilton

The repo is cloned at `/nobackup/pjjg18/ts-bench/TreeSearch/`.

### Single-branch rebuild (fast, cached objects)

When rebuilding the same branch after a small change, in-place install
reuses cached `.o` files — only changed translation units recompile:

```bash
BENCH=/nobackup/pjjg18/ts-bench
REPO=$BENCH/TreeSearch
LIB=$BENCH/lib-baseline

cd $REPO
git pull origin <branch>
R CMD INSTALL --library=$LIB "$REPO"
```

### A/B profiling rebuild (tarball, clean objects)

When switching branches to build two variants (e.g. `lib-baseline` vs
`lib-optimized`), you **must** use the tarball method. In-place install
compiles `.o` files directly in `src/`, so stale objects from branch A
contaminate the build of branch B. The tarball copies sources to a temp
directory, keeping the repo clean.

```bash
BENCH=/nobackup/pjjg18/ts-bench
REPO=$BENCH/TreeSearch
LIB=$BENCH/lib-baseline   # or lib-optimized

cd $REPO
git checkout <branch>
git pull origin <branch>

rm -f src/*.o src/*.so
TMPDIR=$(mktemp -d)
(cd "$TMPDIR" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$REPO")
R CMD INSTALL --library=$LIB "$TMPDIR"/TreeSearch_*.tar.gz
rm -rf "$TMPDIR"
```

The `rm -f src/*.o src/*.so` before each build is belt-and-suspenders:
the tarball method already isolates compilation, but this ensures no
stale artifacts are accidentally picked up if someone later does an
in-place build.

## Installing R package dependencies

Submit `install_deps.R` as a SLURM job (needs network access from
compute node):

```r
lib <- '/nobackup/pjjg18/ts-bench/lib-baseline'
.libPaths(c(lib, .libPaths()))
install.packages(
  c('Rcpp', 'ape', 'cli', 'fastmatch', 'lifecycle', 'Rdpack',
    'phangorn', 'TreeTools', 'PlotTools', 'TreeDist', 'Quartet',
    'Rogue', 'protoclust', 'abind'),
  lib = lib,
  repos = 'https://cloud.r-project.org',
  Ncpus = 4L
)
```

## Writing benchmark R scripts

Benchmark scripts load TreeSearch from a specified library:

```r
lib_dir <- "/nobackup/pjjg18/ts-bench/lib-baseline"
.libPaths(c(lib_dir, .libPaths()))
library(TreeSearch)
library(TreeTools)
```

Standard benchmark datasets live in the repo clone:
- `TreeSearch/dev/benchmarks/mbank_X30754.nex` — 180 taxa, 418 patterns
  (large-tree benchmark)
- `TreeSearch/inst/` — bundled datasets (Vinther2008, etc.)

Write CSV results to `/nobackup/pjjg18/ts-bench/results/` for later
download.

## Hardware

- **CPU:** AMD EPYC 7702 64-Core (2 sockets, 128 cores, no SMT)
- **Clock:** ~2.0 GHz base
- **Cache:** 32K L1d/L1i, 512K L2, 16M L3 per CCX
- **ISA:** x86_64 with AVX2 support
- **Partition:** `shared` (default for single-node jobs)
