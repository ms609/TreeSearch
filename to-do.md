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
| T-207 | P2 | ASSIGNED (C) | — | **SA perturbation phase in `run_single_replicate()`.** Add multi-cycle PCSA (SA+TBR with best-tree restart) after drift, before final TBR. EW ≥100 tips. GHA-test. | On `feature/pt-eval` (TS-PTeval). |

### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-196 | P2 | PR #215 (M) | — | **[Bug] `extract_divided_steps` wrong for NA+IW.** Four static copies read `local_cost` for NA blocks instead of three-pass correction. Conservative (final `score_tree()` always correct), but suboptimal move selection. | Found by S-RED focus 10. Fix committed on `feature/parallel-temper` (`6dc28a2`); arrives with PT PR #215. |
| T-208 | P2 | OPEN | — | **[Bug] `random_topology_tree` ignores constraints.** When `adaptiveStart=true` (thorough preset) AND constraints active, bandit can select RANDOM_TREE → starting tree violates constraint → TBR blocks all constraint-relevant moves (`cn=-1`). Could return constraint-violating trees. | Found by S-RED focus 11. Fix: fall back to WAGNER_RANDOM when `strategy == RANDOM_TREE && cd && cd->active` in `run_single_replicate()` (ts_driven.cpp:111). |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-179 | P2 | PR #215 (M) | — | **Large-tree strategy preset.** For ≥120 tips. | On `feature/parallel-temper`. Commit `fab1e52c`. Arrives with PT PR #215. |
| T-206 | P3 | PARKED (E, GHA 23503681296) | — | **Outer cycle reset cap / minimum-Δ gate.** `outerCycles=1` repeats until no improvement; late cycles yield <1 step/s. Add reset cap or minimum-Δ threshold. Also fix misleading comment at `ts_driven.cpp:180`. | On `feature/outer-cap-t206`. Retry after test fix. |
| T-182 | P3 | ASSIGNED (G) | — | **Adaptive ratchet perturbation probability.** Taper by hit rate as pool stabilizes. | On `feature/adaptive-ratchet`. |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** | Constraint infrastructure exists (`consensus_constrain`). |
| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop after `nTip × K` unsuccessful perturbations. | From T-185 IQ-TREE review. |


### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-24 by F (focus 11: T-190/T-202/XPIWE merge review). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-24 by E (supplement: outer cycle reset analysis, T-206 filed). Round 4 by G (re-baseline). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-24 by A (round 11: cleanup T-190/T-177/T-202, pipeline refresh). |


