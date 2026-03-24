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


| T-150 | P2 | WORKTREE (TS-CID-cons) | — | **CID-optimal consensus tree search** | Use SPR/TBR/NNI + Ratchet to find consensus trees minimising CID to input trees. Needs: collapse/resolve moves for non-binary trees, CID batch scorer, `CIDBootstrap` (resample input trees). See `briefing-cid-consensus.md`. PoC validated in TreeDist. |




### Parallel Tempering (Objective 17)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-190 | P2 | ASSIGNED (C) | — | **Stochastic TBR + Boltzmann acceptance.** `stochastic_tbr_phase()` in `ts_temper.h/.cpp`: random clip+regraft, indirect scoring, temperature-scaled acceptance. Building block for hot chains. | Plan: `.positai/plans/2026-03-24-parallel-tempering.md` |
| T-191 | P2 | OPEN | T-190 | **Multi-chain parallel tempering framework.** N ChainStates with temperature ladder, round structure (K moves then swap), Metropolis swap criterion, best-score promotion. | |
| T-192 | P2 | OPEN | T-191 | **Pipeline integration.** Wire into `run_single_replicate()` as new phase. `SearchParams` fields, `SearchControl()` R API, timing diagnostics. | |
| T-193 | P2 | OPEN | T-192 | **Benchmark evaluation.** Compare PT vs current pipeline at equal wall-clock. 180t dataset + additional large matrices. Time-adjusted expected best. | |

### Bugs from S-RED Focus 10

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-196 | P2 | ASSIGNED (G) | — | **[Bug] `extract_divided_steps` wrong for NA+IW.** Four static copies (ts_tbr.cpp, ts_search.cpp, ts_drift.cpp, ts_temper.cpp) read `local_cost` for all blocks, but NA blocks need the three-pass correction from `extract_char_steps`. IW candidate screening uses wrong base score + deltas. Impact: conservative (final `score_tree()` always correct), but suboptimal move selection for `inapplicable="bgs"` + finite `concavity`. Fix: add NA block loop (matching ts_fitch.cpp:420-465) to each copy, or extract a shared function. | Found by S-RED focus 10. |
| T-197 | P3 | ASSIGNED (B) | — | **[Bug] `concavity = 0` produces NaN.** `precompute_iw_delta` computes `0/(0+0)=NaN` when `k=0` and `e=0`. R-side `TreeLength` also has `0/0` at line 132. No input validation rejects `k ≤ 0`. Fix: validate `concavity > 0` in R layer (`MaximizeParsimony`, `TreeLength`, `AdditionTree`), or handle `k=0` as limit case in C++. | Found by S-RED focus 10. |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-177 | P1 | ASSIGNED (Human+AI) | — | **Bug fix: mid-TBR/SPR timeout.** Pass `check_timeout` callback into `tbr_search()` and `spr_search()` so they can bail out mid-pass. At 180+ tips, a single TBR convergence takes ~13s; without mid-pass timeout, `maxSeconds` can overrun significantly. | **Implemented**, building and testing. Added `std::function<bool()>` param to `tbr_search`/`spr_search`; threaded through driven, ratchet, drift calls. 1762 Tier 2 tests pass. |
| T-179 | P2 | OPEN | T-177, T-178 | **Large-tree strategy preset.** Add a preset for ≥120 tips (or ≥100 tips with ≥200 patterns): NNI→SPR→TBR escalation, scaled ratchet/drift cycles, tuned sector sizes. Benchmark on 180-taxon dataset. | May need `nniFirst` param in `SearchControl()`. |


| T-182 | P3 | OPEN | — | **Adaptive ratchet perturbation probability.** Start aggressive (~40%) and taper by hit rate as pool stabilizes. Extend `adaptive_level` infrastructure to also scale `perturb_prob`. Risk: premature convergence if tapering too fast. | Also consider IQ-TREE-style tree-size scaling: 50% at ≤50 tips, 5% at ≥400 tips. Two orthogonal axes: hit-rate tapering (within run) and size scaling (across datasets). |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** Use pool consensus as backbone constraint during Wagner construction for later replicates (after ≥N diverse trees). Partially randomised addition order. Concern: run independence — mitigate by only activating late. | Constraint infrastructure already exists (`consensus_constrain`). |

| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop search after `nTip × K` unsuccessful perturbations (no score improvement). Dataset-adaptive convergence criterion complementary to time-based `maxSeconds` and consensus-stability stopping. IQ-TREE uses `nSeq × 100`. | From T-185 IQ-TREE review. Needs counter in `driven_search()` loop tracking perturbations since last improvement. |



### Alternative Inapplicable-Handling Algorithms — Phase 3: Integration & polish

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|



### Shiny Search UX Improvements

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|







### Infrastructure

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-195 | P3 | OPEN | — | **GHA benchmark workflow.** Create `agent-benchmark.yml` for performance regression testing on GHA. Adapt `bench_regression.R` to accept CLI args (`--datasets`, `--budget`). Upload results as artifacts. | See `.positai/plans/2026-03-24-0625-remote-compute-integration-for-multi-agent-workflows.md` §1b. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-19 by B (focus 7: Shiny module wiring). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-19 by A (round 3). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-23 by G (round 9). |
