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
| T-197 | P3 | ASSIGNED (D) | — | **[Bug] `concavity = 0` produces NaN.** `precompute_iw_delta` computes `0/(0+0)=NaN` when `k=0` and `e=0`. R-side `TreeLength` also has `0/0` at line 132. No input validation rejects `k ≤ 0`. Fix: validate `concavity > 0` in R layer (`MaximizeParsimony`, `TreeLength`, `AdditionTree`), or handle `k=0` as limit case in C++. | Found by S-RED focus 10. |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|



| T-190 | P2 | ASSIGNED (A) | — | **Adaptive starting-tree strategy mixing (bandit).** Replace fixed `wagnerBias` with per-replicate strategy selection using Thompson sampling. See cold-start brief below. | Subsumes the "mix strategies" aspect of T-183. Benchmark on 180-taxon mbank_X30754. |
| T-182 | P3 | OPEN | — | **Adaptive ratchet perturbation probability.** Start aggressive (~40%) and taper by hit rate as pool stabilizes. Extend `adaptive_level` infrastructure to also scale `perturb_prob`. Risk: premature convergence if tapering too fast. | Also consider IQ-TREE-style tree-size scaling: 50% at ≤50 tips, 5% at ≥400 tips. Two orthogonal axes: hit-rate tapering (within run) and size scaling (across datasets). |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** Use pool consensus as backbone constraint during Wagner construction for later replicates (after ≥N diverse trees). Partially randomised addition order. Concern: run independence — mitigate by only activating late. | Constraint infrastructure already exists (`consensus_constrain`). |

| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop search after `nTip × K` unsuccessful perturbations (no score improvement). Dataset-adaptive convergence criterion complementary to time-based `maxSeconds` and consensus-stability stopping. IQ-TREE uses `nSeq × 100`. | From T-185 IQ-TREE review. Needs counter in `driven_search()` loop tracking perturbations since last improvement. |











### T-190 Cold-Start Brief: Adaptive Starting-Tree Strategy Mixing

**Problem.** The current driven search uses a fixed starting-tree strategy
(`wagnerBias`) for all replicates. This is suboptimal because:
- Different datasets favour different strategies; picking the wrong one wastes
  replicates.
- Even when one strategy dominates, diversity of starting basins improves
  exploration. Goloboff (2014 §3.3) found random starting trees sometimes
  necessary — they access basins that Wagner trees systematically miss.

**Idea.** Treat starting-tree strategy as a multi-armed bandit. Each replicate
draws its strategy from a probability distribution. After each replicate,
update the probabilities based on whether the replicate hit the best score.
Use Thompson sampling (Beta-Bernoulli) for the explore/exploit tradeoff.

**Strategy arms.**

| Arm | Type | Description | Available |
|-----|------|-------------|-----------|
| WAGNER_RANDOM | Fresh | Random addition-order Wagner tree | Always |
| WAGNER_GOLOBOFF | Fresh | Goloboff (2014) non-ambiguous-char-biased Wagner | Always |
| WAGNER_ENTROPY | Fresh | State-specificity-biased Wagner | Always |
| RANDOM_TREE | Fresh | Purely random topology (no character data used) | Always |
| POOL_RATCHET | Pool-based | Pick random pool tree → heavy ratchet → TBR | Pool ≥1 tree |
| POOL_NNI_PERTURB | Pool-based | Pick random pool tree → NNI perturbation → TBR | Pool ≥1 tree |

Fresh-start arms: build tree → NNI warmup → TBR → [XSS → ratchet → drift → TBR].
Pool-based arms: perturb pool tree → TBR → [XSS → ratchet → drift → TBR].
Pool-based arms skip Wagner + initial descent, investing time in escape instead.
Give RANDOM_TREE a lower prior (e.g. Beta(1,2)) so it starts with ~33%
probability rather than 50%, reflecting the expectation that it's usually
worse but allowing the data to override.

**Thompson sampling.**
```
For each arm i: maintain successes[i], failures[i] (init: 1,1 or 1,2)
To select: sample theta[i] ~ Beta(successes[i], failures[i]); pick argmax.
After replicate: if hit best score → successes[arm]++; else failures[arm]++.
On new best score: decay all counts by ×0.5 (landscape changed, old evidence
is stale).
Pool-based arms excluded when pool is empty (zero-armed at rep 0).
```

**Key code locations.**
- `ts_driven.h`: `DrivenParams` struct (add `adaptive_start` bool + enum)
- `ts_driven.cpp`: `run_single_replicate()` (accept strategy; build
  starting tree accordingly), `driven_search()` loop (strategy selection)
- `ts_wagner.h/cpp`: `WagnerBias` enum, `random_wagner_tree()`,
  `biased_wagner_tree()`
- `build_postorder.h`: `random_tree()` — random topology (needs
  conversion to TreeState format)
- `ts_nni_perturb.h/cpp`: `random_nni_perturb()` — for POOL_NNI_PERTURB
- `ts_ratchet.h`: `ratchet_search()` — for POOL_RATCHET
- `ts_pool.h`: `TreePool` — for selecting a random pool tree
- `ts_parallel.cpp`: worker_thread loop (round-robin assignment for
  parallel path)
- `R/SearchControl.R`: add `adaptiveStart` param
- `R/MaximizeParsimony.R`: pass through to C++; default on for thorough/large
- `inst/benchmarks/bench_framework.R`: benchmarking infrastructure

**Implementation plan.**
1. Add `StartStrategy` enum to `ts_driven.h` (6 arms).
2. Add `StrategyTracker` class: Thompson sampling, per-arm Beta params,
   `select()`, `update()`, `decay()`, `set_pool_available()`.
3. Modify `run_single_replicate()` to accept a `StartStrategy` param.
   - Fresh arms: switch on strategy to set `wagner_bias` or build random tree.
   - For RANDOM_TREE: implement `random_topology_tree(TreeState&, DataSet&)`
     that builds a random topology and scores it (adapt `random_tree()` from
     `build_postorder.h`).
   - Pool-based arms: caller passes a pre-perturbed tree as `starting_tree`.
4. Modify `driven_search()` serial loop:
   - Before each rep: `tracker.select()` → strategy.
   - If pool-based and pool non-empty: extract random pool tree, apply
     perturbation (ratchet or NNI), pass as starting_tree.
   - After rep: `tracker.update(strategy, hit_best)`.
   - On new best score: `tracker.decay(0.5)`.
5. Modify `parallel_driven_search()`: assign strategies round-robin across
   replicates (pre-computed array). No adaptation (too complex for marginal
   benefit in parallel).
6. R API: `SearchControl(adaptiveStart = TRUE)` → maps to
   `DrivenParams.adaptive_start`. Default TRUE for thorough/large presets.
7. Benchmark: 180-taxon mbank_X30754, 5 seeds × {fixed-random, fixed-Goloboff,
   adaptive} × 60s budget. Compare best score, time-to-best, per-arm hit rates.
8. Tests: unit test StrategyTracker (selection, update, decay, pool-gating);
   integration test that adaptive search runs without error.

**Risks.**
- Attribution noise: perturbation may dominate over starting tree, making
  bandit signal weak. Mitigation: this is fine — the main value is diversity
  from the start, not learning per se.
- Pool-based strategies early in search: pool has only 1 mediocre tree.
  Mitigation: pool-based arms auto-excluded until pool has ≥2 best-score trees.
- Parallel adaptation: not attempted; round-robin is sufficient.

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-24 by B (focus 10: Profile & IW scoring). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-19 by A (round 3). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-23 by G (round 9). |
