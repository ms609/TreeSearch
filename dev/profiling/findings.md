# Profiling findings

One row per verified optimisation opportunity, in `to-do.md` paste-ready
format. A finding only lands here if an isolated `std::chrono` micro-bench
reproduces the predicted delta.

Tags:
- `[Port]` — R loop on the hot path that should move to C++.
- `[Optimise]` — C++ change with verified expected speedup.
- `[AT-LIMIT]` — function is at a hardware ceiling; record so the rotation
  skips it in future rounds.

| ID-suggest | P? | Status | Depends | Headline | Detail (% time, mechanism, verified Δ, micro-bench path) |
|------------|----|--------|---------|----------|---------------------------------------------------------|
| T-300 | P1 | OPEN | — | [Optimise] `full_rescore` after accepted TBR move (ts_tbr.cpp:1138): replace with incremental rescore | 19.2 % of TreeSearch.dll CPU time (0.617 s / 3.211 s); confirmed via VTune hotspots 2026-05-19. `ts::fitch_na_score` 18.2 % + `load_tip_states` 1.0 %. Per-accepted-move full rescore is O(n_nodes × n_chars); incremental uppass/downpass path already exists (fitch_na_pass3_score + fitch_na_incremental_uppass). Not micro-bench verified (T-300 implementation in flight by parallel agent). |
