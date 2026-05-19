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
| (none yet — scaffold round)                                            |
