# Agent E — Progress Log

## Current Task
- **Task:** T-260 — Per-evaluation overhead profiling (VTune)
- **Status:** IN PROGRESS. Building with symbols, writing driver, collecting hotspots.

### T-260 — Per-evaluation overhead profiling
- Goal: profile TBR per-evaluation overhead to find why TreeSearch evaluates 1.5–3.6× fewer rearrangements/s than TNT
- VTune 2025.10 at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest`
- Steps: build with symbols → driver script → VTune hotspot collection → analysis

### Previous: T-255 — Reduce drift from presets — PARKED
- GHA 23591874696 in progress. Fixed test-ts-anneal.R:106 (7dc2ed96).

### Previous: S-RED focus 7 — recently landed code — DONE
- Reviewed T-258, PR #230, T-255. No bugs. 624 targeted tests pass.
