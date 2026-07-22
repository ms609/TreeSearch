# Profiling focus areas

Ranked by `(estimated wall-time share) × (remaining fixability)`. Areas at the
memory-bandwidth ceiling or already optimised drop to the bottom but remain
visible so the rotation knows to skip them.

Signals used to build this list:
- Phase distribution from `.positai/expertise/profiling.md` (Zhu2013, 75 t,
  thorough preset, post-T-261/T-262/T-263, 2026-03-27).
- Hot Rcpp entries grepped from `src/ts_rcpp.cpp`.
- Active and parked profiling tasks in `to-do.md` (T-274, T-298, T-300,
  S-PROF round 7).
- File mtimes vs S-PROF last-run date (2026-05-12).

Statuses: `NEW` (never profiled), `PROFILED` (profiled, no fix yet),
`OPTIMISED` (fix shipped — re-profile if `src/` changes), `AT-LIMIT` (no
further wins — skip unless code changes), `SKIPPED` (out of rotation).

| #  | Area                                | Files                                                                  | Why hot                                                                                       | Last known cost                                  | Last profiled | Status     |
|----|-------------------------------------|------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|--------------------------------------------------|---------------|------------|
| 1  | NNI-perturb in driven pipeline      | `src/ts_nni_perturb.cpp`, `src/ts_driven.cpp` (perturb call sites)     | Disabled in thorough preset via T-274 (`nniPerturbCycles=0L` in `R/MaximizeParsimony.R`) — code on path only when caller sets `nni_perturb_per > 0` | T-274 filed; disabled at R level; re-evaluate only if default changes | 2026-05-18    | AT-LIMIT   |
| 2  | Ratchet inner loop                  | `src/ts_ratchet.cpp`, `src/ts_tbr.cpp` (called from ratchet)           | 62 % of inner-loop search time (verbosity=2, Zhu2013 thorough, 2026-05-18); TBR dominates (perturbation overhead < 2 %) | 2.80 s/rep median (Zhu2013 thorough ×1 rep, nThreads=1); T-300 (`full_rescore`) is pending fix | 2026-05-18    | PROFILED   |
| 3  | RSS / sector search                 | `src/ts_sector.cpp`, `src/ts_prune_reinsert.cpp`                       | THROUGHPUT at-limit by inheritance (R6 2026-06-20): ~96 % is inner+global tbr_search (at-limit kernel); sector scaffolding ≤2 %. Banked T-S6c byte-identical ~2.8 %; T-S6d per-clip getenv ~22 % (TBR-wide). Efficiency axis (work-to-target) untouched. | inner tbr_search-dominated; see findings R6 | 2026-06-20    | AT-LIMIT   |
| 4  | TBR full-rescore at acceptance      | `src/ts_tbr.cpp:1138` (`full_rescore` after every accepted move)       | T-300 RESOLVED — dirty-set incremental rescore landed for SPR accept (EW path `fitch_dirty_*`, NA path `fitch_na_dirty_*`); GHA-green; 15.2 % wall-time speedup on Zhu2013 (3.88 s → 3.29 s) confirmed via dev/profiling/t300_na_bench.R 2026-05-19 | resolved | 2026-05-19    | DONE       |
| 5  | quartet_concordance.cpp allocation  | `src/quartet_concordance.cpp`                                          | T-298 active PR #242 — matrix allocation hoist already benchmarked; re-profile after merge    | hoist-fix in flight                              | 2026-05-12    | PROFILED   |
| 6  | CSS / XSS sector pipeline           | `src/ts_sector.cpp`, `src/ts_simplify.cpp` (`ts_simplify_diag` entry)  | Same verdict as #3 (R6): XSS uses search_sector (=RSS scaffolding ≤2 % + inner tbr_search); CSS uses sector-masked tbr_search directly — both inner-tbr-dominated ⇒ THROUGHPUT at-limit by inheritance. T-S6c levers + T-S6d getenv apply to all three modes. | inner tbr_search-dominated | 2026-06-20    | AT-LIMIT   |
| 7  | Hierarchical resampling parallelism | `src/ts_resample.cpp`, `R/Resample.R`                                  | HSJ/XFORM hierarchical resampling 2-thread speedup 1.1× (vs Brazeau 2.5×) — serial R loop     | known limitation (2026-03-19 Agent A)            | —             | NEW        |
| 8  | Simulated-annealing phase           | `src/ts_temper.cpp`                                                    | 7.4 % at 180 t, 14 % hit rate, 0.8 steps/s — efficiency far below ratchet (4.5) or XSS (13.8) | 1241 ms/rep at 180 t (T-179)                     | —             | NEW        |
| 9  | MaddisonSlatkin solver              | `src/MaddisonSlatkin.cpp`                                              | Hash-map infrastructure was 53 % of DLL time; T-151/T-152 raised but check if shipped         | ~1.4 s (6 %) gain estimated cold-cache           | 2026-03-19    | PROFILED   |
| 10 | Wagner tree construction            | `src/ts_wagner.cpp`                                                    | < 0.1 % of search time on all datasets ≤ 88 t                                                 | 300–1400 µs / tree                               | 2026-03-18    | AT-LIMIT   |
| 11 | Per-candidate indirect scoring      | `src/ts_driven.cpp`, `src/ts_fitch.cpp`, `src/ts_fitch_na_incr.h`      | At memory-bandwidth ceiling (~23 ns / 75 tips). T-075 confirmed no further wins.              | 23 ns / candidate (75 t)                         | 2026-03-18    | AT-LIMIT   |
| 12 | R-loop search engine (`MaxParsi`)   | `R/MaximizeParsimony.R`, `R/TreeSearch.R`                              | < 0.5 % of wall time (Rprof, 2026-03-18). R is a passenger, not a bottleneck on hot path.     | < 0.5 %                                          | 2026-03-18    | AT-LIMIT   |
| 13 | **Standard-Fitch path (TNT-parity)** | `src/ts_tbr.cpp`, `src/ts_fitch.cpp`, `src/ts_simd.h`, `src/ts_tree.cpp` | The `-`→`?` path (has_na=FALSE) the TNT benchmark uses; ~20× faster/rep than NA. tbr_search self 25 %, SIMD 21 %, uppass 13 %, per-clip bookkeeping 18 % | Zhu2013 627 in 0.56 s/rep; total DLL 2.70 s/8rep | 2026-06-16 (r3) | PROFILED |
| 13a| → per-clip bookkeeping              | `ts_tbr.cpp` (build_postorder, collect_*, compute_from_above), `ts_tree.cpp` | postorder rebuilt every clip+accept (5.2 %) + incremental down/uppass (6.4 %) — TNT minimises exactly this | — | 2026-06-16 | NEW (top code lever) |
| 13b| → SIMD reduce / uppass arithmetic   | `ts_simd.h`, `ts_fitch.cpp:54`                                        | any_hit_reduce 21 % (compiler-optimal); uppass scalar loop 1.22× only | — | 2026-06-16 | AT-LIMIT |

## Notes on the ranking

- Area #1 (NNI-perturb) is now AT-LIMIT — disabled in the thorough preset
  via T-274 (`nniPerturbCycles=0L`); the code path is only active when a
  caller explicitly sets that parameter, and profiling dead code is wasteful.
- Areas #2–#3 are the live wins: Ratchet has the largest absolute share
  (now ~60–70 % after NNI-perturb disabled) with no per-line profile yet;
  RSS grew 16× without a profile pass to confirm cost source.
- Area #4 (T-300 lazy rescore) is PARKED in `to-do.md` but stays in the
  rotation because the path is well-understood and the predicted gain is
  large; rerun after T-300 lands to verify.
- Areas #5 and #9 are PROFILED — skip unless their files change.
- Areas #10–#12 are AT-LIMIT — recorded so the rotation skips them.
- Per the skill's `init` rules, IO setup, config readers, and CLI parsers
  are SKIPPED (not listed) and only profiled on explicit `/profile <name>`.

## Profvis-required cases

Per the skill's tool guide, areas with non-trivial R surface area on the hot
path must go through `profvis` before VTune so we don't miss a `[Port]`
finding:

- Area #7 (hierarchical resampling) — bulk of the loop is in R.
- Area #12 (R-loop search) — already known to be < 0.5 %, but rerun profvis
  if `R/MaximizeParsimony.R` or `R/TreeSearch.R` change to confirm.

All other areas are pure C++ on the hot path and can go straight to VTune
after the dry-run.
