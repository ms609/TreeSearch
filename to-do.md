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

### Parallel Tempering / SA (Objective 17)

T-198–T-201 (PT core) are on PR #215. T-199 evaluation (agent-c) found
Boltzmann PT is broken for parsimony but PCSA (post-convergence SA with
best-tree restart) is highly effective under EW at 125+ tips. See
`TS-PTeval/dev/pt_t199_findings.md`.

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-198–201 | P2 | PR #215 (C) | — | **PT core + pipeline integration.** Boltzmann PT disabled by default. | On `feature/parallel-temper`. |
| T-207 | P2 | PR #222 (C) | — | **SA perturbation phase in `run_single_replicate()`.** Multi-cycle PCSA integrated. | On `feature/pt-eval` (TS-PTeval). Includes T-210 fix. GHA 23509475416 PASS. |

### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-196 | P2 | PR #215 (M) | — | **[Bug] `extract_divided_steps` wrong for NA+IW.** Four static copies read `local_cost` for NA blocks instead of three-pass correction. Conservative (final `score_tree()` always correct), but suboptimal move selection. | Found by S-RED focus 10. Fix committed on `feature/parallel-temper` (`6dc28a2`); arrives with PT PR #215. |
| T-208 | P2 | PR #219 (G) | — | **[Bug] `random_topology_tree` ignores constraints.** When `adaptiveStart=true` AND constraints active, bandit selects RANDOM_TREE → constraint violation. | PR #219 open (WAGNER_RANDOM fallback). |

| T-210 | P2 | PR #222 (C) | — | **[Bug] SA doesn't save best-found topology.** Fix: `anneal_search` tracks/restores best tree at phase boundaries. | On `feature/pt-eval` (TS-PTeval). In T-207 PR #222. |
| T-211 | P2 | ASSIGNED (C) | — | **[Bug] Stale `final_` in temper candidate scoring.** Same stale-score pattern as SPR: cached score not refreshed after topology changes, biasing candidate selection. | Race: A released, C investigating. |

### Testing & Constraint Handling

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-212 | P2 | OPEN | — | **Test `random_constrained_tree` under RANDOM_TREE strategy.** Add testthat test that forces `StartStrategy::RANDOM_TREE` with active constraints via `ts_driven_search`. Verify constraint satisfaction and score validity. Exercises the parallel round-robin path indirectly. | `Rf_error` posthoc assertion at ts_wagner.cpp:976 is worker-thread-unsafe (S-RED focus 4). |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-179 | P2 | PR #215 (M) | — | **Large-tree strategy preset.** For ≥120 tips. | On `feature/parallel-temper`. Commit `fab1e52c`. Arrives with PT PR #215. |

| T-182 | P3 | PR #221 (G) | — | **Adaptive ratchet perturbation probability.** Taper by hit rate as pool stabilizes. | On `feature/adaptive-ratchet`. GHA 23508899686 PASS. |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** | Constraint infrastructure exists (`consensus_constrain`). |
| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop after `nTip × K` unsuccessful perturbations. | From T-185 IQ-TREE review. |


### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-24 by A (focus 4: parallelism & RNG, Rf_error-on-worker noted). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-24 by E (supplement: outer cycle reset analysis, T-206 filed). Round 4 by G (re-baseline). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-24 by A (round 13: 7 open PRs, all agents idle, project in holding pattern). |


