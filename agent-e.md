# Agent E — Progress Log

## Current Task
- **Status:** IDLE (T-289 parked SLURM 16608629; S-RED area 4 done 2026-03-27)

### S-RED Area 4 — Parallelism & RNG — DONE (2026-03-27)

Reviewed ts_rng.h/.cpp (110 lines) and ts_parallel.h/.cpp (732 lines).
ts_driven.cpp covered in E-003 (see below).

**No bugs found.** Thread safety correct throughout.

Observations (non-bugs):
- fuse_round holds pool mutex across entire tree_fuse() call (O(n) TBR
  exchanges). Workers block for full fuse duration. Performance only.
- Multiple workers may trigger fuse_round at the same `replicates_done`
  checkpoint due to relaxed read races. Redundant fuse (harmless).
- Lines 323-325 in main polling loop: empty if-block, dead code.
- Verbosity Rprintf acquires pool mutex via status(). If fuse_round holds
  the lock, interrupt/timeout polling is delayed by fuse duration.
- ts_rng.h serial/parallel dispatch verified correct in all paths.

### S-RED Focus 4 — ts_driven.cpp review — DONE (2026-03-27)

Reviewed ts_driven.cpp (1054 lines) and ts_driven.h (322 lines) in full.
Focus areas per AGENTS.md: cross-replicate constraint tightening, outer
cycle loop, and features added since T-189.

**Bugs fixed (committed to cpp-search):**

1. **`unsuccessful_reps` not reset on fuse improvement** (`ts_driven.cpp`
   line ~923). When inter-replicate fusing found a better score, the
   perturb-stop counter was not cleared. Meanwhile `last_improved_rep`
   *was* updated by fuse. This inconsistency could cause `perturb_stop`
   to fire prematurely when fusing is still productive. Low severity
   (factor defaults 0; limit = n_tips × factor is high when enabled),
   but logically wrong. Fixed by adding `unsuccessful_reps = 0;` in the
   fuse-improvement branch.

2. **`DrivenResult::perturb_stop` flag missing** (`ts_driven.h`).
   T-276 ("print convergence summary") explicitly lists perturb_stop as
   a convergence indicator. Added the field and set it at the stopping
   site in `driven_search()`.

3. **Stale NNI-perturb comment** (step 4b, `ts_driven.cpp`). Opening
   sentence said "Skip when constraints are active" but the code passes
   `cd` through `nni_perturb_search()` and has been safe under constraints
   for several tasks. Replaced with accurate one-liner.

**Other observations (no fix needed):**

- `consensus_constrain = true` with 0 unanimous splits calls
  `extract_consensus_splits()` every replicate (performance, not
  correctness). consensus_constrain defaults to false; low priority.
- `timed_out = true` is set for both timeout and user interrupt — no
  distinction in DrivenResult. Acceptable for now; T-276 can note this
  in summary text ("search interrupted/timed out").
- `score_tree()` called at top of each outer cycle for improvement
  comparison — minor overhead, by design.
- MPT enumeration uses user constraint `cd` only (not auto_cd) — by
  design: enumeration should be unconstrained.
- Outer cycle reset logic correct; `score_before_cycle` / `score_after_cycle`
  correctly bound the improvement check.
- Adaptive level, ratchet taper, consensus constraint tightening, and
  adaptive start bandit logic all look correct.

### Previous work
### ASAN vector OOB fix — DONE (2026-03-26)
- Root cause: total_words == 0 when all characters are parsimony-uninformative
- Fix: early returns in TBR, SPR, NNI, drift, ratchet, collapsed-flags
- Commit: 6505803f on cpp-search

### S-COORD Round 27 — DONE
- Fixed R 4.1 `%||%` compat bug in `test-ts-anneal.R` (58fc2552)

### T-265 — RESOLVED (scoring method confound)
### Previous: S-RED Focus 8, T-261+T-262, T-255, T-260
