# Fast iteration for profiling/benchmark work

How to keep the edit→measure loop in seconds–minutes instead of 10+ minutes.
Validated 2026-06-16 (Windows / Rtools45 / R-devel, i7-10700 8C/16T). Numbers below
are measured on this box.

## TL;DR loop

| Step | Command | Stop condition | Wall-clock |
|---|---|---|---|
| Build (C++ edit) | `Rscript dev/build-fast.R [lib]` | — | **~3s** (1 TU) / 0.5s no-op |
| Build (header edit) | `rm -f src/*.o && Rscript dev/build-fast.R [lib]` | — | **~19s** (ccache+`-j8`) |
| Smoke (every edit) | `Rscript dev/benchmarks/bench_smoke.R` | `maxReplicates` | **~0.2s**, exit 1 on regress |
| Iterate gate (pre-commit) | `Rscript dev/benchmarks/bench_iterate.R` | `maxReplicates` | ~1–2 min |
| Batch panel | `Rscript dev/benchmarks/bench_parallel.R` | `maxReplicates` | ~5–7× serial |
| Hamilton panel | `hamilton_*` job-array chain (below) | `maxReplicates` | ~10–15 min incl. queue |

Two iron rules (the reason this works):
1. **Replicate-bounded, never `maxSeconds`, for any candidate/score signal.**
   `candidates_evaluated` is only deterministic when the run stops on replicate count
   (or target hits). Under `maxSeconds`, replicates-completed (hence candidates) depend on
   machine load — non-reproducible. Verified: serial vs parallel give *identical* candidate
   counts when replicate-bounded.
2. **Measure with the pool drained.** Any candidate/timing number must run alone
   (8 physical cores, memory-bandwidth-bound Fitch). The process pool / SLURM array is for
   *batch* panels only, never for the single authoritative measurement.

## One-time machine setup (global, reversible)

Already applied on this box. To reproduce elsewhere:

1. **ccache** (object cache). `ccache --max-size=20G` (3GB thrashes — it was 99.96% full).
2. **`~/.R/Makevars.win`** — wires ccache + parallel make into *every* build while
   preserving `-O2`. On Windows this file SHADOWS `~/.R/Makevars` (no merge), so the `-O2`
   flags are reproduced in it; do not drop them or timings silently become invalid:
   ```make
   MAKEFLAGS = -j8
   CCACHE = ccache
   CC = $(CCACHE) gcc
   CXX = $(CCACHE) g++
   CXX11 = $(CCACHE) g++
   CXX14 = $(CCACHE) g++
   CXX17 = $(CCACHE) g++
   CXXFLAGS = -g -O2 -Wall -mfpmath=sse -msse2 -mstackrealign
   PKG_CXXFLAGS =
   ```
   Revert by deleting the file.

## Build: `dev/build-fast.R`

`compileAttributes()` (guards the stale-`RcppExports` trap) + `pkgbuild::compile_dll(debug=FALSE)`
(incremental, **-O2** — recompiles only changed TUs) + hot-swap the DLL into the install lib
(`.agent-p0` by default) so benchmarks pick it up without a full `R CMD INSTALL`.

- **C++-only edits:** `Rscript dev/build-fast.R`. ~3s for one `.cpp`.
- **Header edits:** R's Makefile has no header-dependency tracking (a `.h` change recompiles
  nothing), so `rm -f src/*.o` first, then build-fast — ccache + `-j8` makes the forced full
  rebuild ~19s, not ~90s.
- **R / roxygen / `[[Rcpp::export]]` signature changes:** need a full `R CMD INSTALL`
  (build-fast only rebuilds C++). For a release/CI build, full install + tests.
- **-O0 fast builds** (`compile_dll(debug=TRUE)`, ~2–4× faster compile) are for
  *correctness/logic* loops only — tag them to a throwaway lib and NEVER benchmark them.
  Candidate *counts* are opt-invariant (safe at -O0); *timing/throughput/profiles* are not.

Do not have an R session holding the target DLL open during a hot-swap (Windows file lock).

## Measurement tiers

- **`bench_smoke.R`** — breakage tripwire. One R process, 3 tiny datasets, `maxReplicates=4`,
  seed 1. Compares score (must be unchanged) + candidates (±5%) vs `smoke_baseline.csv`;
  exit 1 on regress. NOT a ship gate — tiny datasets don't exercise sectorial search.
  Re-baseline with `SMOKE_WRITE_BASELINE=1`.
- **`bench_iterate.R`** — the real lever gate. Gap panel, fixed `maxReplicates`, 2–3 seeds,
  median candidates + score + gap-to-TNT. A win = lower median candidates at equal/better
  score. ~0.7% seed spread on candidates (vs the ±2–4 step score lottery) → 2–3 seeds suffice.
- **`bench_parallel.R`** — PSOCK pool batch runner (`conc = cores − TS_HEADROOM`, OMP=1 per
  worker). For batch panels/sweeps only. Raise `TS_HEADROOM` while anything else runs.
- **`bench_tnt_headtohead.R`** — full TNT comparison (validate tier); `bench_phase_yield.R` —
  per-phase wall-clock shares.

## Hamilton job-array (parallel validation panels)

Convert the serial in-R panel loop to one SLURM task per `(dataset × seed)` cell. Chain:
build once → array → merge.

```bash
# from a Hamilton login node (scripts scp'd to /nobackup/$USER/TreeSearch/scripts, CRLF-stripped)
bid=$(sbatch --parsable hamilton_build_once.sh)                       # compile ONCE into shared $LIB
aid=$(sbatch --parsable --dependency=afterok:$bid hamilton_panel_array.sh)  # one cell per array task
sbatch --dependency=afterany:$aid hamilton_merge.sh                   # rbind partials -> panel.csv
```

- `bench_cell.R` runs one cell from `$SLURM_ARRAY_TASK_ID` (replicate-bounded), writing one
  partial CSV; locally testable: `Rscript dev/benchmarks/bench_cell.R 0`.
- Array tasks NEVER build (the `afterok` build job populates a read-only `$LIB`).
- Collapses a ~4.5h serial panel to ~10–15 min incl. queue (not the optimistic 3–5 min:
  count SLURM queue + ~2–5s NFS R-startup per cell + the build/merge jobs).
- **Blocker for the authoritative wall-clock head-to-head:** there is NO 64-bit Linux TNT on
  Hamilton (no `module load tnt`, none under `/nobackup`). Candidate-count comparison is
  bitness-independent (valid on local 32-bit TNT); only the wall-clock ratio needs a 64-bit
  TNT staged under `/nobackup` and fed to `bench_tnt_headtohead.R` via `TNT_EXE`. The
  Hamilton `*_hamilton.sh` SLURM templates here are **dispatch-untested** (cell logic is
  validated locally); smoke them in `test.q` (5-min partition) before a full run.

## Planned refinement: a true candidate-budget stop

The cleanest iterate signal is "score at a fixed candidate budget" / "candidates to reach a
target score", which removes even the ~0.7% per-seed replicate-spend variance. Needs a small
C++ change: a `max_candidates` guard (and optional `target_score` early-exit) checked at the
replicate boundary in `src/ts_driven.cpp` (alongside the timeout at ~`:627`, the timed-out
set at ~`:897`/`:1043`, and the target-hits check at ~`:1006`), plumbed through
`runtimeConfig` → `ts_rcpp.cpp` → `SearchControl`. Until then, fixed `maxReplicates` is the
deterministic signal.
