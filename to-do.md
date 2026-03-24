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
| T-150 | P2 | WORKTREE (TS-CID-cons) | — | **CID-optimal consensus tree search** | PR #213 open to cpp-search. |
| T-204 | P2 | PR #216 (B) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | On `feature/native-search`. GHA run 23495097795. |

### Parallel Tempering (Objective 17)

**Note:** These tasks were originally misnumbered T-190–T-193, colliding
with existing completed IDs. Renumbered to T-198–T-201 in S-COORD round 10.

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-198 | P2 | DONE (C), needs PR | — | **Stochastic TBR + Boltzmann acceptance.** `stochastic_tbr_phase()` in `ts_temper.h/.cpp`. | Was T-190 (PT). On `feature/parallel-temper`. |
| T-199 | P2 | DONE (C), needs PR | T-198 | **Multi-chain parallel tempering framework.** N ChainStates with temperature ladder, Metropolis swaps. | Was T-191 (PT). On `feature/parallel-temper`. |
| T-200 | P2 | DONE (C), needs PR | T-199 | **Pipeline integration.** Wired into `run_single_replicate()` as new phase. | Was T-192 (PT). On `feature/parallel-temper`. |
| T-201 | P2 | DONE (C), needs PR | T-200 | **Benchmark evaluation.** PT vs current pipeline at equal wall-clock. | Was T-193 (PT). On `feature/parallel-temper`. |

### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-196 | P2 | DONE (G), on PT branch | — | **[Bug] `extract_divided_steps` wrong for NA+IW.** Four static copies read `local_cost` for NA blocks instead of three-pass correction. Conservative (final `score_tree()` always correct), but suboptimal move selection. | Found by S-RED focus 10. Fix committed on `feature/parallel-temper` (`6dc28a2`); will arrive with PT PR. |
| T-202 | P2 | ASSIGNED (B) | — | **[Bug] MPT enumeration skipped on timeout.** `!result.timed_out` guard in `driven_search()` and `parallel_driven_search()` prevents MPT plateau walk after timeout. At ≥40 tips, timeout is the normal exit, so pool often has only 1–4 trees. Fix: remove the `!result.timed_out` guard; also remove `check_timeout` from MPT enum's TBR call. | Investigated by F; root cause confirmed. See `agent-f.md`. |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-177 | P1 | ASSIGNED (Human+AI) | — | **Bug fix: mid-TBR/SPR timeout.** | Implemented, building and testing. 1762 Tier 2 tests pass. |
| T-179 | P2 | DONE (G), needs PR | T-177, T-178 | **Large-tree strategy preset.** For ≥120 tips. | On `feature/parallel-temper`. Commit `fab1e52c`. |
| T-203 | P2 | ASSIGNED (G) | — | **Simulated annealing for large trees.** Single-chain linear cooling schedule using `stochastic_tbr_phase()`. Replace drift in `large` preset. | On `feature/anneal` (TS-anneal worktree). |
| T-190 | P2 | ASSIGNED (A) | — | **Adaptive starting-tree strategy mixing (bandit).** Thompson sampling over 6 strategy arms. | See cold-start brief below. |
| T-182 | P3 | OPEN | — | **Adaptive ratchet perturbation probability.** Taper by hit rate as pool stabilizes. | |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** | Constraint infrastructure exists (`consensus_constrain`). |
| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop after `nTip × K` unsuccessful perturbations. | From T-185 IQ-TREE review. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-24 by B (focus 10: Profile & IW scoring). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-20 by A (round 3 + MaddisonSlatkin opts). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-24 by D (round 10). |

---

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
- `dev/benchmarks/bench_framework.R`: benchmarking infrastructure

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
