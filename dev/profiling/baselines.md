# Profiling baselines

Snapshot of the workloads used by `/profile` and the timings they currently
produce on this machine. Refreshed each round for the area profiled.
`/profile regress` reruns these drivers and flags any > 10 % slowdown.

Machine context belongs alongside each driver entry (CPU, cores used,
power profile) so future regressions can be compared apples-to-apples.

## Driver baselines

| Driver                              | Dataset / N    | Bare wall (s) | Top hotspot (mod=TreeSearch.dll) | % share | Recorded   | Machine note          |
|-------------------------------------|----------------|---------------|----------------------------------|---------|------------|-----------------------|
| dev/profiling/drivers/ratchet.R     | Zhu2013 / 1 rep thorough nThreads=1 | 2.80 (median of 3) | ts_driven_search (>95 %; no VTune; from profvis) | >95 % | 2026-05-18 | Windows 10 i-series, R-devel, .vtune-lib debug build |
| dev/profiling/drivers/tbr-rescore.R | Zhu2013 / 12 ratchet reps nCycles=12 nThreads=1 | 3.9 | ts::fitch_na_score (full_rescore path via callstack) | 18.2 % | 2026-05-19 | Windows 10 EARTHSCI-PJJG18, 2.904 GHz 16-core, R-devel, .vtune-lib-20260519061049 (HEAD c504ea87) |
| dev/profiling/t300_na_bench.R       | Zhu2013 / 12 ratchet reps nCycles=12 nThreads=1 | 3.29 (median of 5, score 647) | (post T-300 NA dirty-set; fitch_na_pass3_score expected dominant; not VTune-attributed yet) | n/a | 2026-05-19 | Same machine, HEAD 5b210fdd; 15.2 % wall-time speedup vs c504ea87 baseline (3.88 s median of 3) |
| dev/profiling/drivers/fitch-tnt.R   | Zhu2013 `-`→`?` / 8 reps auto→thorough nThreads=1 | 5.57 (0.56 s/rep, score 627) | ts::tbr_search (orchestration self) | 25.1 % | 2026-06-16 | Same machine, HEAD 841eead3, .vtune-lib-20260616052323 (-O2 -g). STANDARD-Fitch path (has_na=FALSE) — TNT-parity objective; TNT 1.6 = 624 |

### Round 3 top hotspots (TreeSearch.dll, Zhu2013 75t **standard-Fitch** `-`→`?`, total 2.70 s)

Different path from rounds 1-2 (NA): no `fitch_na_*`; flat/x4 kernels.
Names via `nm` (VTune CSV shows `func@0x…`; image base 0x2cc1a0000 stable).

| Rank | Function                              | Self time | % of DLL | Notes |
|------|---------------------------------------|-----------|----------|-------|
| 1    | ts::tbr_search (orchestration)        | 0.678 s   | 25.1 %   | candidate-loop control + collapsed/sector vector<bool> bit-tests + inlined scoring |
| 2    | ts::simd::any_hit_reduce_avx2         | 0.392 s   | 14.5 %   | 2-op Fitch reduce; AT-LIMIT (compiler-optimal, disasm-confirmed) |
| 3    | ts::uppass_node                       | 0.357 s   | 13.2 %   | incremental uppass; scalar update loop; AT-LIMIT (1.22× micro-bench, nil for 2-state) |
| 4    | ts::simd::any_hit_reduce3_avx2        | 0.171 s   | 6.3 %    | 3-op reduce (SPR bounded) |
| 5    | ts::TreeState::build_postorder_prealloc | 0.141 s | 5.2 %    | O(n) per clip + per accept — top per-clip-bookkeeping target |
| 6    | ts::fitch_incremental_downpass        | 0.110 s   | 4.1 %    | per clip |
| 7    | ts::fitch_indirect_bounded_flat       | 0.109 s   | 4.0 %    | SPR candidate scoring (flat) |
| 8    | ts::hash_tree                         | 0.078 s   | 2.9 %    | pool/tabu dedup |
| 8    | ts::fitch_indirect_length_cached      | 0.078 s   | 2.9 %    | scalar cached (MIXED-ratchet perturbed sub-search) |
| 8    | ts::validate_topology                 | 0.077 s   | 2.9 %    | per-accept DFS sanity check (allocates 2 vectors/call) |

### Round 2 top-5 hotspots (TreeSearch.dll, Zhu2013 75t ratchet)

| Rank | Function                         | Self time | % of DLL | Notes |
|------|----------------------------------|-----------|----------|-------|
| 1    | ts::fitch_na_score               | 0.585 s   | 18.2 %   | Full Fitch pass (called via full_rescore → tbr_search, confirmed by callstack) |
| 2    | ts::simd::any_hit_reduce_avx2    | 0.309 s   | 9.6 %    | SIMD candidate hit reduction — inner evaluation loop |
| 3    | ts::tbr_search (residual)        | 0.297 s   | 9.3 %    | Control-flow overhead not attributed to child callees |
| 4    | ts::fitch_na_pass3_score         | 0.281 s   | 8.8 %    | Incremental scoring uppass (candidate evaluation) |
| 5    | ts::fitch_na_incremental_uppass  | 0.110 s   | 3.4 %    | Incremental uppass after candidate topology |

**full_rescore (ts_tbr.cpp:1138) total = fitch_na_score + load_tip_states = 0.617 s = 19.2 % of DLL CPU time**
Note: prior S-PROF round 7 estimate was 28 %; measured 19.2 % (see Round 2 log for context).
Context: virtually all of fitch_na_score flows through line 1138 (acceptance path), not line 563 (entry call), because ratchet-driven TBR accepts ~100–200 moves per sub-optimal restart vs 1 entry call.

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
