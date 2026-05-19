# Profiling baselines

Snapshot of the workloads used by `/profile` and the timings they currently
produce on this machine. Refreshed each round for the area profiled.
`/profile regress` reruns these drivers and flags any > 10 % slowdown.

Machine context belongs alongside each driver entry (CPU, cores used,
power profile) so future regressions can be compared apples-to-apples.

## Driver baselines

| Driver                          | Dataset / N    | Bare wall (s) | Top hotspot (mod=TreeSearch.dll) | % share | Recorded   | Machine note          |
|---------------------------------|----------------|---------------|----------------------------------|---------|------------|-----------------------|
| dev/profiling/drivers/ratchet.R | Zhu2013 / 1 rep thorough nThreads=1 | 2.80 (median of 3) | ts_driven_search (>95 %; no VTune; from profvis) | >95 % | 2026-05-18 | Windows 10 i-series, R-devel, .vtune-lib debug build |

## End-to-end reference timings (from `.positai/expertise/profiling.md`)

Kept here for cross-check against driver-level numbers — *not* a substitute
for the project benchmark suite.

| Dataset       | Tips | Chars | Median wall (s) | Score   | Preset          | Recorded   |
|---------------|------|-------|-----------------|---------|-----------------|------------|
| Vinther2008   | 23   | 57    | 0.42            | 79      | sprint          | 2026-03-19 |
| Agnarsson2004 | 62   | 242   | 1.79            | 778     | default         | 2026-03-19 |
| Zhu2013       | 75   | 253   | 3.17            | 648–666 | thorough        | 2026-03-19 |
| Dikow2009     | 88   | 220   | 4.90            | 1612–14 | thorough        | 2026-03-19 |
| mbank_X30754  | 180  | 425   | 17.3 / rep      | 1202    | large, 30 s bud | 2026-03-26 |
