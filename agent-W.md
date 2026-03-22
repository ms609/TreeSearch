# Agent W — XPIWE Feature Branch

**Worktree**: `C:/Users/pjjg18/GitHub/TS-Xpiwe`
**Branch**: `feature/xpiwe`
**Library**: `.agent-W/` (inside `../TS-Xpiwe/`)
**Status**: IDLE — ready to begin T-156

## Resolved design decisions

| Question | Decision |
|---|---|
| Formula | `e / (k + m + e)`, i.e. `eff_k = k + global_min[p]`. Verify against Goloboff (2014) §2 before first commit, but proceed — TNT + secondary sources agree. |
| Resample / SA in Phase 1? | **Yes — include.** R-side change is trivial; ships a consistent API. |
| Silent-ignore EW + `extended_iw`? | **Yes — silent ignore.** No warning; parameter is vacuously correct. |

## Task log

_No tasks completed yet._
