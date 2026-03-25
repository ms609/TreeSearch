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
| T-210 | P2 | PR #222 (C) | — | **[Bug] SA doesn't save best-found topology.** Fix: `anneal_search` tracks/restores best tree at phase boundaries. | On `feature/pt-eval` (TS-PTeval). In T-207 PR #222. |
### Testing & Constraint Handling

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-212 | P2 | PARKED (B, GHA 23543892219) | — | **Test `random_constrained_tree` under RANDOM_TREE strategy.** Tests on cpp-search. GHA 23528636505 failed (59 inapplicable + constraint failures). Needs re-dispatch now that T-214 is done. | Was blocked by T-214 (now complete). |

| T-235 | P3 | ASSIGNED (B) | — | **[Bug] SPR search stale state arrays after rejected regraft.** In `spr_search()`, when a candidate passes indirect screening but fails `full_rescore`, `spr_unclip()` only partially restores states (clip-to-root), leaving other nodes with regrafted-topology states. Degrades screening for subsequent clips (NA/IW datasets only). Final score always correct. Fix: add `full_rescore()` after `spr_unclip()` on the rejection path. | Found by S-RED focus 2. SPR rarely used (all presets disable `sprFirst`). |

### Large-Tree Scaling & Search Optimization (Objective 15)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-179 | P2 | PR #215 (M) | — | **Large-tree strategy preset.** For ≥120 tips. | On `feature/parallel-temper`. Commit `fab1e52c`. Arrives with PT PR #215. |

| T-182 | P3 | PR #221 (G) | — | **Adaptive ratchet perturbation probability.** Taper by hit rate as pool stabilizes. | On `feature/adaptive-ratchet`. GHA 23508899686 PASS. |
| T-183 | P3 | OPEN | — | **Pool-seeded Wagner / consensus backbone.** | Constraint infrastructure exists (`consensus_constrain`). |
| T-187 | P3 | OPEN | — | **Perturbation-count stopping rule.** Stop after `nTip × K` unsuccessful perturbations. | From T-185 IQ-TREE review. |

### Shiny App

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|

| T-232 | P2 | PARKED (D, GHA 23543699366) | — | **[Shiny] "Tips to show" input bounces back on decrement.** Clicking "down" arrow resets to previous value (e.g. 54 for Sun dataset). | From a.013. Fix: `isolate(input$keepNTips)` in `UpdateKeepNTipsRange`. |

| T-226 | P2 | OPEN | — | **[Shiny] "Trees in sequence" connect mode — review/remove.** May not make sense under new C++ search engine (no meaningful replicate ordering). | From a007. Design question. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-25 by C (focus 2: Search topology invariants). Found T-235 (SPR stale state). |
| S-PROF | dyn | BLOCKED: Do not run this task until 2026-03-26 | — | **Standing: Performance profiling** | Last run: 2026-03-24 by E (supplement: outer cycle reset analysis, T-206 filed). Round 4 by G (re-baseline). |
| S-COORD | dyn | ASSIGNED (D) | — | **Standing: Coordination review** | Round 19. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-25 by B. Open PRs: #210 (cpp-search→main, MERGEABLE), #213/#215/#216/#221/#222 (feature branches, all failing GHA or CONFLICTING). #224 closed. |

