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



### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-242 | P1 | PARKED (C, GHA 23545987517) | — | **[Bug?] Agnarsson2004 IW search quality regression.** 230 runs, only 5 hit best score (2% hit rate). User reports "1 trees in memory: 1 sampled, each with score 50.1872 (k = 5.62)". May indicate search regression or IW landscape difficulty. | From a.20. Investigate whether this is a genuine regression or expected IW behaviour. |



### Shiny App

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-232 | P2 | PARKED (D, GHA 23543699366†) | — | **[Shiny] "Tips to show" input bounces back on decrement.** Clicking "down" arrow resets to previous value (e.g. 54 for Sun dataset). | From a.013. Fix committed. GHA failed pre-T-214 fix; awaiting re-validation on current HEAD (run 23547582438 queued). |
| T-240 | P2 | PARKED (D, GHA 23544604214†) | — | **[Shiny] Pool suboptimal filter not applied when changed mid-search.** After search, changing "Keep if suboptimal by" from ≤6 to ≤2 or ≤0 doesn't filter existing pool trees. | From a.17. Fix committed. GHA failed pre-T-214; awaiting re-validation. |
| T-239 | P3 | PARKED (D, GHA 23545538742†) | — | **[Shiny] Cluster consensus: highlight edges unique to a cluster.** Heatmap colouring: emphasize "unique to cluster" vs "in 5/6 clusters". Agnarsson (6 clusters) is testbed. | From a.07. Feature committed. GHA failed pre-T-214; awaiting re-validation. |
| T-241 | P3 | PARKED (D, GHA 23545261957†) | — | **[Shiny] Show cluster assignment next to tree selector.** Add "(cluster X)" in cluster colour after "Tree to plot" label. | From a.19. Feature committed. GHA failed pre-T-214; awaiting re-validation. |

### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-25 by F (focus 6: R↔C++ interface). Clean — no bugs. |
| S-PROF | dyn | BLOCKED: Do not run this task until 2026-03-26 | — | **Standing: Performance profiling** | Last run: 2026-03-24 by E (supplement: outer cycle reset analysis, T-206 filed). Round 4 by G (re-baseline). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-25 by F (round 20). Cleaned stale entries, updated PR refs. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-25 by F (round 20 triage). Open PRs: #216 (native-search), #213 (CID-consensus). #226/#227 merged, #215/#222 CLOSED. #210 (draft cpp-search→main). #178/#106 stale — consider closing. |

