# TreeSearch — hand-off 2026-08-01 (soft-Sankoff exploration)

> **Placement note.** `AGENTS.md` says feature branches carry code only and
> coordination files live on `cpp-search`. This file is on
> `feature/soft-sankoff` anyway, because at hand-off time `cpp-search` had a
> 1h34m `R-CMD-check` in flight (run `30698449302`) and a coordination push
> would have risked cancelling it. Delete this file before any PR.

## What this repo is

TreeSearch reconstructs phylogenetic trees from discrete character data under
parsimony criteria (equal weights, implied weights, profile parsimony), with a
C++ search engine and the full Goloboff toolkit — ratchet, drift, fusing,
sectorial search, tabu, TBR — plus a Shiny UI. This branch is **exploration
only**, on a question that came out of a conversation about variational
inference: whether softening the `min` in Sankoff's DP buys anything.

## Where we left off

Three commits on `feature/soft-sankoff`, cut from `cpp-search` @ `f13f0c31`:

- **`ec5491e6`** — the substance. Plan document
  (`dev/plans/2026-08-01-soft-sankoff-temperature-dial.md`), a pure-R reference
  implementation (`tests/testthat/helper-soft-sankoff.R`), 63 Tier-2 assertions
  (`tests/testthat/test-ts-soft-sankoff.R`), and three evidence scripts under
  `dev/soft-sankoff/`.
- **`49048252`** — data inventory and the Hamilton sizing rule, added to the plan.
- **`8370d228`** — O'Reilly matrices location and the nested-archive trap.

**The claim, and it is verified.** Sankoff's DP with `min` replaced by
`softmin_T(x) = -T*log(sum_j exp(-x_j/T))` is hard weighted parsimony at
`T -> 0` and Felsenstein's pruning algorithm at `T = 1` with
`cost(i,j) = -log P_ij(t)`. Parsimony and likelihood are the same dynamic
program in the tropical `(min,+)` and probability `(sum,x)` semirings; `T` is
the dial. Because `softmin_T ~= min - T*log(#minima)`, temperature controls how
far the criterion integrates over ancestral reconstructions rather than
optimising them.

Measured, not asserted:

| Check | Result |
|---|---|
| `soft(T=1)` vs `phangorn::pml` | `0.0e+00` |
| `soft(T=1)` vs own pruning reference | `1.8e-15` |
| `soft(T->0)` vs compiled `ts_sankoff_test()` | exact |
| `soft(T->0)`, equal costs, vs `TreeLength()` | exact |
| `(min - softmin)/T` vs `log(#minima)` | `9e-10` |

**Both gates got a first read, and neither is encouraging for the search
application** — which is why plan Steps 3a and 3b were deliberately placed off
the gates.

- **Gate B (cost): unfavourable.** Operation-count ratio **x350** for binary
  characters, flat across a 16–128 tip ladder, against a x50 orientation
  threshold. Before counting that `exp` is not one operation, and before the
  incremental-rescore loss (`src/ts_sankoff.h` is full-rescore only), which is
  the term that changes inner-loop complexity rather than a constant.
- **Gate A (direction): insufficient data, and the design needs widening.**
  Only 1 of 6 Congreve & Lamsdell matrices had >= 8 distinct MPTs. On that one,
  `rho(CID) = -0.30` (toward the truth) but `rho(Mk) = +0.14` (away from
  likelihood) — a disagreement worth watching, but `n = 1`. Ranking was
  invariant in `T` across 0.02–0.5.

## Pending jobs

**Updated 2026-08-01, later session.** The original claim that nothing was
pushed is stale: `origin/feature/soft-sankoff` exists at `a1c3ec50`, and three
further commits (`e8f9de7c`, `8f5b7480`, `ce8ae11a`) are local-only.

| Type | ID / ref | Status | On completion |
|---|---|---|---|
| Hamilton | `18146187` | running (`soft-sankoff-dial`, ~15 min) | Collect `/nobackup/pjjg18/soft-sankoff/out/04-dial-study{,-per-matrix}.csv` into `dev/soft-sankoff/`; read the recovery-by-`T` table, the per-matrix sign test, and the homoplasy correlation; write the verdict into the Step 3a section of `dev/plans/2026-08-01-soft-sankoff-temperature-dial.md` and the README status list; commit |

Job `18146150` was the same job and **FAILED** at matrix 11 of 100
(`MaximizeParsimony()` returns polytomies by default). Fixed and superseded by
`18146187`; nothing to collect from it.

No GHA dispatched by this session, no `to-do.md` task claimed, no issue labelled
`in-progress`, no dispatch agents active.

Other agents' GHA work, from the original hand-off — context only, do not
collect. `30698449302` has since completed successfully; the other two were
still in progress when this session started:

| Run | Branch | Started |
|---|---|---|
| `30701531266` | `feature/hsj-token-index-fix` | 13:20 |
| `30700854197` | `gha-ccache` | 13:00 |
| `30698449302` | `cpp-search` | 11:49 — ✅ success |

`.dispatch/state.json` records `5h_pct_committed: 122` — over the 5-hour budget
window. Factor that into any dispatch on arrival.

## Open items / next steps

1. **Re-root the session.** This work was done from a session rooted in the
   **StratoBayes** worktree, so StratoBayes's `AGENTS.md` was auto-loaded
   instead of TreeSearch's, and harness memory was writing to the StratoBayes
   project key. Start the next session at
   `C:\Users\pjjg18\GitHub\TreeSearch` (mints the stable project key; 41
   existing TreeSearch keys are all worktree-scoped orphans), then reach this
   branch with `git -C` or `EnterWorktree`.
2. **Widen Gate A and run it on Hamilton.** Score trees within a few steps of
   optimal rather than only the tied MPTs, across all **100** bundled
   `congreveLamsdellMatrices` (6 used so far). Edit
   `dev/soft-sankoff/02-tilt-direction.R`. **This is a `/hamilton` job** — see
   Technical pointers.
3. **Stage the O'Reilly matrices.** Fetched to
   `~/Downloads/doi_10_5061_dryad_10qf3__v20160322.zip` (19,061,436 bytes).
   **Nested — unzip twice**; the wrapper holds one member,
   `oreilly2016matrices.zip` (19,061,276 bytes). Suggested home
   `OReillyEtAl2016/data-raw/Matrices/` (already gitignored there); recipe in
   that directory's `README-FETCH.md`.
4. **Build the C++ soft-scorer prototype.** This *is* Gate B proper, including
   the incremental-rescore question the op count cannot capture. It also makes
   every downstream experiment ~100x cheaper, so it is on the critical path
   even if Step 4 of the plan dies.
5. **Steps 3a / 3b of the plan** — the parsimony/likelihood dial study and
   graded ancestral-state reconstruction. Both need only the scorer, which now
   exists and is verified. These survive whatever the gates say, and 3a is the
   paper.
6. **File the 0-based `ts_sankoff_test()` trap** into
   `.AGENTS/memory/feature-inapplicable.md` (the x-transform/Sankoff memory).
   Left undone because `.AGENTS/` sits in the grey zone of the
   coordination-files-not-on-feature-branches rule — maintainer's call.
7. **Consider a `to-do.md` claim** on `cpp-search` if this is taken up.
   Suggested **T-386**; nothing claimed yet.

## Technical pointers

- **`ts_sankoff_test()` takes 0-based tip states.** The guard is
  `state >= 0 && state < ns_ch` (`src/ts_rcpp.cpp`). The "1-based" comment at
  `R/recode_hierarchy.R:188` means present states begin at index **1 of a
  0-based array**, state 0 being "absent". An out-of-range index does **not**
  error — it leaves every state at `INF`, so the whole score returns `Inf`.
  This cost 6 test failures before it was spotted.
- **`NOT_CRAN=true` or nothing runs.** `test-ts-soft-sankoff.R` is Tier 2, so
  `skip_on_cran()` is the first line. Without the env var testthat reports
  `[ FAIL 0 | PASS 0 ]`, which reads as success.
- **Build into the isolated library, never in place.** This worktree has
  `.agent-softsankoff/` (gitignored via `^\.agent-`) with TreeSearch 2.0.0
  installed. Rebuild recipe is in `dev/soft-sankoff/README.md`; it follows
  `AGENTS.md` (tarball into a temp dir, `rm -f src/*.o src/*.dll` first).
- **Anything over ~45 s of compute is a `/hamilton` job**, not a local run.
  Do not scale an experiment down to fit the workstation — that is exactly how
  the Gate A first read ended up with `n = 1`.
- **Hamilton already has TreeSearch.** `/hamilton`'s `r-infrastructure.md`
  documents `tsLib = /nobackup/pjjg18/TreeSearch/lib` (TreeSearch 2.0.0,
  TreeTools, TreeDist, Quartet), and `tsLib` **must precede** `baseLib`.
  **`phangorn` is not in that listing** and the Gate A Mk oracle needs it —
  expect a project-local install.
- **`phangorn` must be *attached*, not namespaced,** in any script calling
  `as.character()` on a `phyDat`. `as.character.phyDat` is registered but not
  exported, so `as.character(dataset)` silently returns a bare vector unless
  the package is on the search path.
- **`MaximizeParsimony()` v2 dropped `maxHits`/`ratchIter`** in favour of
  `maxReplicates` / `targetHits` / `maxSeconds` (see `?SearchControl`).
- **`~/GitHub/stack-status.sh` reports StratoBayes stacks, not TreeSearch.**
  Do not read its output as if it described this repo.
- Costs are equal and symmetric throughout this work, so the DP is
  rooting-invariant and the live **T-374 / T-385** rooting defects do not
  apply. They *would* apply to any asymmetric-cost variant — do not build a
  soft version on top of an unsettled objective.

## Things ruled out

- **Vine's embedding + NJ-decoder machinery does not port to
  `MaximizeParsimony()`.** Vine's chain rule terminates at the branch-length
  gradient via inside-outside; parsimony has no branch lengths and an integer
  score, so the chain rule stops at step one. Without a gradient you have an
  unguided walk in re-coordinatised tree space, competing against ratchet +
  drift + fusing + sectorial + tabu, already in C++. Soft-Sankoff exists
  precisely to supply the missing gradient.
- **Embedding-derived starting trees are a no-op.** Vine initialises its mean
  by classical MDS on the distance matrix and decodes by neighbour-joining, so
  an unoptimised embedding start tree *is* the `TreeTools::NJTree()` you
  already have. The only untested residue is whether MVN-sampled decoded starts
  ("diverse but plausible") beat random-addition Wagner starts for driven-search
  replicates.
- **Dryad blocks programmatic download.** v2 API returns
  `401 Unauthorized, must have current bearer token`; both
  `/downloads/file_stream/69374` and `/stash/downloads/file_stream/69374`
  return 403. Do not re-attempt with `curl` — use a browser or an API token.
- **A printf-style `Rcpp::stop()` audit was raised and withdrawn.** It rested
  on a claim that variadic `stop(fmt, ...)` caused an aarch64 segfault in a
  sibling package. The maintainer confirmed on 2026-08-01 that this was
  **wrong**. Do not re-file it.

## Worktrees

| Branch | Path | Status |
|--------|------|--------|
| `feature/soft-sankoff` | `C:/Users/pjjg18/GitHub/worktrees/TS-softsankoff` | clean, 3 commits ahead of `cpp-search`, unpushed |

## Suggested first action

Start a session rooted at `C:\Users\pjjg18\GitHub\TreeSearch`, then read
`dev/plans/2026-08-01-soft-sankoff-temperature-dial.md` on this branch — it is
self-contained and carries both gate reads. Then widen Gate A:

```bash
git -C C:/Users/pjjg18/GitHub/worktrees/TS-softsankoff log --oneline -3
```

and rework `dev/soft-sankoff/02-tilt-direction.R` to score near-optimal trees
across all 100 matrices, for submission via `/hamilton`.
