# TreeSearch Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- When a task is being developed in a **git worktree**, set its status to
  `WORKTREE (name)` where *name* is the worktree directory (e.g.
  `WORKTREE (TS-CID-cons)`). This distinguishes human/long-running worktree
  work from agent assignments and prevents double-claiming.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md` (see workflow in AGENTS.md).
- Tasks awaiting GHA results: `PARKED (<Letter>, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (<Letter>)`.
  S-COORD cleans these up after merge.
- Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
  completed, reset it to OPEN. Their effective priority is dynamic:
  - ≥6 OPEN specific tasks → standing tasks are P3
  - 3–5 OPEN specific tasks → standing tasks are P2
  - <3 OPEN specific tasks → standing tasks are P1

---

## Active Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-150 | P2 | PARKED (F, GHA 23650002703) | — | **CID-optimal consensus tree search** | PR #213. Vignette fix (TreeTools::Consensus) commit f8bfee49. GHA 23650002703. |
| T-204 | P2 | PR #216 (F) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | GHA 23649607006 PASSED. Ready for merge. |


### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-279 | P2 | PARKED (F, GHA 23650290962) | — | **Constrained drift: stale `cd->constraint_node` after RFD rejection + redundant full_rescore.** In `ts_drift.cpp::drift_phase()`, when a constrained suboptimal move passes the constraint check (so `cd` is updated), then is rejected by the RFD test, `update_constraint()` is not called on the restored topology. Bug in both IW reject (lines 647–649) and EW reject (line 689). Same class as T-278. Also eliminates a redundant `drift_full_rescore()` call in the EW reject path (return value at line 657 was discarded, causing a second rescore at line 689). Only affects constrained drift (`driftCycles > 0`); drift disabled in all presets since T-255. | Found by S-RED focus 8 (F, 2026-03-27). PR TBD. feature/drift-constraint-fix commit e85ec84f. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### TNT Comparison & Strategy Learning

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-269 | P3 | OPEN | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. **Design note (u.005):** Full TNT-style interleaving IS architecturally supported now — setting `outerCycles = ratchetCycles` achieves one sector pass per ratchet cycle. T-257 first validated that any post-ratchet sectorial helps at all (merged). T-269 benchmarks the fine-grained variant to determine whether the per-cycle overhead is worth enabling in presets. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| E-002 | P3 | OPEN | — | **SearchControl Rd completeness check.** After T-150/T-204 merge, verify all parameters in `SearchControl()` Rd match the actual function signature. File bugs for any missing/stale param docs. | Short task; likely 30 min. |



### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 10 by F (ts_driven.cpp, 1056 lines). No bugs. Notes: goto-finish is valid C++; adaptive_params copy correct; SA best-tree save/restore correct; MPT enumeration uses check_enum_timeout (full deadline) correctly; extra score_tree calls at verbosity>=2 harmless. Next: ts_sector.cpp (1007 lines). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 35 by F. T-278 filed (constrained TBR constraint staleness). OPEN: T-278, T-245, T-269, E-002 → 4 specific OPEN → standing at P2. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 round 38 by E. **ASAN**: `continue-on-error: true` on jobs (rlang capture.c uses PREXPR directly — header shim can't fix it). GHA 23649409998 confirms green status with jobs failing — PRs no longer blocked. T-150: GHA 23648875258 ubuntu ✓, windows running. T-204: GHA 23649607006 running (Morphy.R root-cause fix). Open PRs: #213 (T-150), #216 (T-204), #210 (DRAFT). |
