# TreeSearch Multi-Agent Development Notes

## Current phase: bug-fixing / pre-release (as of 2026-03-20)

The project is in a **bug-fixing and stabilisation phase** with the goal of
shipping the package. Agents should:

- Monitor `to-do.md` as usual for task selection.
- **Prioritise bug fixes, test failures, documentation issues, and R CMD check
  problems** over new functionality.
- **Do not implement new features on `cpp-search` or `main`.**  
  Feature work is allowed only on dedicated `feature/<name>` branches, and
  only when the task is explicitly labelled as a feature and has been approved
  for active development.
- When in doubt, prefer a conservative fix (minimal diff, no API changes) over
  an ambitious refactor.

This phase ends when a clean `R CMD check` (0 errors, 0 warnings) is confirmed
and the maintainer signals readiness to tag a release.

---

## Validation workflow — GHA first (mandatory)

**Use GitHub Actions for all validation:** R CMD check, full test suites,
and benchmarks. Local builds are for **targeted iteration only** (editing
code → building → running one or two specific test files to check your
change). Never run a full test suite or R CMD check locally.

### GHA dispatch (primary validation path)

```bash
# Push your branch and dispatch checks
git push -u origin feature/<name>
cd ..
bash gha-dispatch.sh agent-check.yml feature/<name>

# Poll for results
bash gha-poll.sh <run_id>
```

Park the task while waiting and pick up another (see root `AGENTS.md` for
parking protocol). Do **not** block waiting for GHA results.

### Local builds (targeted iteration only)

Multiple agents share the same `src/` directory. In-place `R CMD INSTALL .`
compiles `.o` files and links the DLL directly in `src/`, causing races.

**Always build via tarball** so compilation happens in an isolated temp
directory:

```bash
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-X "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
```

Key points:
- `rm -f src/*.o src/*.dll` **must** precede every build — stale artifacts slow traversal and corrupt DLLs.
- Build into an agent-specific `$TMPBUILD` outside the source tree — avoids tarball collision when multiple agents build concurrently.
- `--no-resave-data` skips unnecessary `.rda` re-saving (not needed for dev installs).

Run **targeted** tests only:
```bash
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='test-ts-foo')"
```

**Never** use `R CMD INSTALL --library=.agent-X .` (in-place build).

**Never** install to the default library. On Windows, a loaded DLL locks
the file and blocks other agents.

**Never** use `devtools::load_all()` or `pkgbuild::compile_dll()` — these
target a shared temp location and will conflict.

**Never** run full test suites or R CMD check locally — use GHA.

## Build failure recovery

### Debug `.o` contamination

`roxygen2::roxygenise()` (default mode) calls `pkgbuild::compile_dll(debug=TRUE)`,
which leaves debug `.o` files in `src/`. Subsequent `R CMD INSTALL` reuses them,
producing a DLL that crashes at runtime (exit code 127/139).

**Fix:** `rm -f src/*.o src/*.dll` then rebuild.

**Prevention:** Never use bare `roxygen2::roxygenise()`. To regenerate docs:
```bash
Rscript -e ".libPaths(c('.agent-X', .libPaths())); roxygen2::roxygenise(load_code = roxygen2::load_installed)"
```

### DLL lock

If `R CMD INSTALL` fails with "Access is denied", another R process has the
DLL loaded. Kill it or wait, then retry.

### `TreeSearch-init.c` arg count mismatch

After `Rcpp::compileAttributes()`, **always** run `Rscript check_init.R` to
verify arg counts match between `RcppExports.cpp` and `TreeSearch-init.c`.

### Quick recovery

```bash
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-X "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
Rscript check_init.R
```

## CPU limits — max 2 cores per agent

Use `nThreads = 2L` at most in tests/benchmarks. Never `nThreads = 0L`
(auto-detect). Use `-j2` at most for make.

## Shared files — coordination rules

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` are modified by every agent.
**Append only** — add new entries at the end. Do not reformat or reorder.

### `concavity` sentinel

Rcpp can't translate `R_PosInf`. All Rcpp-exported functions use
`concavity = -1.0` as the C++ default (sentinel for "equal weights / Inf").
Conversion `if (concavity < 0) concavity = HUGE_VAL` happens at three
gateway points in `ts_rcpp.cpp`: `make_dataset()`, `ts_resample_search()`,
`ts_successive_approx()`.

### `src/Makevars.win`

**Never leave a `src/Makevars.win` in place.** Debug/PGO/UBSan flags cause
crashes or miscompilation. Delete after any profiling session.

### `src/TreeSearch-win.def`

**Keep this file.** It explicitly exports `R_init_TreeSearch` for Windows
DLL builds. Without it, the default `nm | sed` pipeline generates a
`tmp.def` that truncates long C++ mangled symbols, causing linker failures
or corrupt DLLs (especially under `pkgbuild::compile_dll(debug=TRUE)`).

## Branch structure

```
main              ← stable, taggable; receives only reviewed bug fixes
  └─ cpp-search   ← integration branch; all feature work merges here
       ├─ feature/cid-consensus
       ├─ feature/hsj-polish
       └─ feature/<name>   (one per major feature)
```

### Rules

- **`main`**: bug fixes and release tags only. No experiments.
- **`cpp-search`**: integration target. **Agents must not merge directly to
  `cpp-search`.** All code changes go through PRs reviewed by the human.
  Coordination-only commits (agent logs, to-do.md updates) may be pushed
  directly.
- **`feature/*`**: branch from `cpp-search`; contain **code changes only**.
  Each feature branch is owned by a single agent at a time.

### Worktree discipline

Each worktree directory is **locked to its branch**. The mapping may change
over time (e.g. when a feature merges and the worktree is reassigned), but
**only the human updates this mapping**.

**Hard rules:**

1. **Never `git checkout <branch>` on any worktree.** This silently replaces
   another agent's (or the human's) working tree and can cause data loss.
2. **Never `git merge` into a worktree's branch** unless that is the explicit
   task you were assigned (e.g. "pull cpp-search into feature/X").
3. To read files from another branch, use `git show <branch>:<path>`.

**Current worktree mapping** (run `git worktree list` to verify):

| Directory | Branch | Purpose |
|-----------|--------|---------|
| `TreeSearch-a` | `cpp-search` | Main source dir; integration branch |
| `TS-anneal` | `feature/anneal` | Simulated annealing for large trees (T-203) |
| `TS-CID-cons` | `feature/cid-consensus` | CID consensus feature |
| `TS-MadSlat` | `feature/madslatkin-profiling` | Mad-Slatkin profiling |
| `TS-ParsSim` | `feature/parssim-ambiguous` | Parsimony simulation |
| `TS-PTeval` | `feature/pt-eval` | Parallel tempering evaluation |
| `TS-NativeSearch` | `feature/native-search` | Native scorer decoupling (T-204) |
| `TS-TNT-bench` | `feature/tnt-bench` | TNT benchmarking |
| `TS-Xpiwe` | `feature/xpiwe` | Extended implied weighting |

There is **no permanent worktree for `cpp-search`**. To commit to
`cpp-search`, use one of these approaches (in order of preference):

**Option A — temporary worktree (safest):**
```bash
git worktree add ../TS-tmp-<Letter> cpp-search
# work in ../TS-tmp-<Letter>, commit, push
git worktree remove ../TS-tmp-<Letter>
```

**Option B — single-file coordination commit (from a feature worktree):**
```bash
git show cpp-search:agent-X.md > agent-X.md   # read current version
# edit the file
git stash                                      # stash any code changes
git checkout cpp-search -- agent-X.md          # stage from cpp-search
cp agent-X.md.edited agent-X.md               # apply your edit
git add agent-X.md && git commit -m "chore: agent X progress note"
git push origin cpp-search
git checkout HEAD -- agent-X.md                # restore feature version
git stash pop                                  # restore code work
```

**Never use Option B for multi-file edits.** Use a temporary worktree.

### Coordination files live on `cpp-search` only

`to-do.md`, `issues.md`, `agent-X.md`, `completed-tasks.md`, `coordination.md`,
and `AGENTS.md` are **never committed on feature branches**. When an agent
working on a feature branch needs to log progress or claim a task, they commit
those changes directly to `cpp-search` (coordination-only commit), keeping the
feature branch clean.

To read coordination files while on a feature branch without switching:
```bash
git show cpp-search:to-do.md
git show cpp-search:agent-X.md
```

### Shared files at merge time

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` use the existing append-only
convention — merge conflicts resolve cleanly by keeping both appended blocks.
`DESCRIPTION` (Collate field) and `NAMESPACE` require a manual merge pass;
this is expected and should be done carefully at feature-merge time.

### Feature branch lifecycle

1. `git checkout cpp-search && git checkout -b feature/<name>`
   Optionally create a worktree: `git worktree add ../TS-<name> feature/<name>`
2. Claim task on `cpp-search`'s `to-do.md` (coordination commit).
3. Do all code work on `feature/<name>`. Use local targeted tests only
   during iteration; use GHA for full validation.
4. When ready: push and dispatch GHA checks:
   ```bash
   git push -u origin feature/<name>
   bash gha-dispatch.sh agent-check.yml feature/<name>
   ```
5. On GHA success, open a PR:
   ```bash
   gh pr create --base cpp-search --head feature/<name> \
     --title "T-nnn: <description>" --body "Agent <Letter>. ..."
   ```
6. Set `to-do.md` status to `PR #N (<Letter>)`. Move on.
7. Human reviews and merges the PR.
8. After merge, clean up:
   ```bash
   git worktree remove ../TS-<name>  # if worktree was used
   git branch -d feature/<name>
   git push origin --delete feature/<name>
   ```

---

## Multi-agent workflow protocol

### Worktree tasks

Tasks with status `WORKTREE (name)` are actively developed in a dedicated git
worktree (e.g. `C:/Users/pjjg18/GitHub/TS-CID-cons`). **Do not claim or
modify these tasks.** They are reserved for the human developer working in
that worktree. To mark a task as in-flight on a worktree, set its status to
`WORKTREE (name)` where *name* matches the worktree directory basename.

### Assignment

On `/assign X`:

1. Read `agent-X.md`. If a task is already in-progress, resume it.
2. **Triage `aXXX.md` / `a.XXX` bug reports** (see "Shiny bug report intake" below):
   a. List all `a[0-9]*.md` **and** `a.[0-9]*` files in the project root
      (excluding any `*.claimed-*.md` / `*.claimed-*` files).
   b. For each file, check its size first. **Skip files shorter than
      20 characters** (likely mid-edit — the human may still be typing).
      Do not rename or touch these files; leave them for a later pass.
   c. For files ≥20 characters, attempt `mv aXXX.md aXXX.claimed-X.md`.
      If the rename fails (file gone or access denied), another agent
      claimed it — skip.
   d. Create a `to-do.md` task under `### Shiny App` for each valid
      report. Assign the next available `T-nnn` ID, a priority based on
      severity (default P2), and tag as `[Shiny]`. Use the file's content
      as the description/notes.
   e. Delete the `aXXX.claimed-X.md` file once the `to-do.md` entry is
      written.
   f. Repeat for all claimed files before moving on. **Do not start
      working a task until all pending reports are triaged.**
3. Check `issues.md` **before** `to-do.md`:
   a. If `issues.md` contains any unclaimed issues (blocks whose first line
      does **not** start with `CLAIMED`), **claim the bottom-most unclaimed
      issue** by prepending `CLAIMED (X):` to its first line.
   b. Triage the claimed issue: determine what needs doing, then add one or
      more discrete tasks to `to-do.md` (assign appropriate IDs and
      priorities — issues may be P0). Begin work on the first task.
   c. Once the `to-do.md` tasks are created, delete the entire issue block
      (including its `---` separator) from `issues.md`.
   d. **While `issues.md` still has unclaimed issues, triaging them takes
      priority over picking up existing `to-do.md` tasks** (an issue may
      contain a P0).
4. If `issues.md` is empty or all issues are already claimed, claim the next
   OPEN task from `to-do.md` as before.

Set `CONVERSATIONSUMMARY` to `Agent X: <task description>`.

> **Concurrency guard (issues.md):** Only the bottom-most *unclaimed* issue
> may be claimed. Because agents always target the bottom and mark it
> `CLAIMED (X)` immediately, two agents will never parse the same issue. If
> an agent sees the bottom issue is already `CLAIMED`, it moves up to the
> next unclaimed one.
>
> **Concurrency guard (aXXX.md / a.XXX):** Atomic rename (`mv aXXX.md
> aXXX.claimed-X.md`) ensures exactly one agent wins each file. NTFS
> rename is atomic; losers see "file not found" and skip.

### During work

- Update `agent-X.md` after every significant step (crash-recovery record).
- All work uses `.agent-X/` as library directory.
- **All builds, tests, and benchmarks in bash subprocesses** — never in the
  RStudio R session.
- **Use GHA for validation** (full test suites, R CMD check, benchmarks).
  Local builds are for targeted iteration only (build + run 1–2 test files).
  See "Validation workflow" section above.

### On task completion

1. **Delete** the task row from `to-do.md`. If the task was the last open
   row in a section/group, delete the section header too.
2. **Append** a summary row to `completed-tasks.md` under the current date
   heading (create a new `## YYYY-MM-DD` heading if needed):
   `| T-nnn | Short description | X | Brief notes |`
3. Set `agent-X.md` to IDLE.
4. Append a brief entry to this file documenting what changed.
5. Update `coordination.md` if strategic objectives are affected.
6. Take next task.

### Shiny bug report intake

The human files Shiny app bugs as individual files in the project root.
Naming convention:

`a.010`, `a.011`, … (dot-separated, no extension)

Each file contains a free-text bug description. The human's workflow is:
create file → write bug → save → never touch the file again.

**Agent responsibility:** Triage all pending bug report files into
`to-do.md` at the start of every `/assign` (step 2 in Assignment above).

**Claim protocol:**
```bash
# List unclaimed reports (both conventions)
ls a.[0-9]* 2>/dev/null | grep -v 'claimed'

# Skip short files (< 20 chars) — don't rename, don't touch
wc -c < a.010  # check size first

# Claim atomically (rename)
mv a.010 a.010.claimed-X

# Read, triage into to-do.md, then delete
cat a.001.claimed-X
# ... create to-do.md entry ...
rm a.001.claimed-X
```

**Skip guard:** Files shorter than 20 characters are likely mid-edit.
Do **not** rename them — just leave them in place for a later pass.
(Renaming and renaming back triggers RStudio "file moved" dialogs.)

**to-do.md placement:** All Shiny bugs go under `### Shiny App`. Create
the section if it doesn't exist. Default priority P2 unless the content
clearly indicates higher severity (crash = P1, cosmetic = P3).

**Shiny bug fixes** are committed directly to `cpp-search` (they are
bug fixes in `inst/Parsimony/`, not feature work, so no feature branch
is needed). Use the temporary-worktree approach if changes span multiple
files.

### Standing tasks

| ID | Type | Expertise file |
|----|------|---------------|
| S-RED | Red-team review | `.positai/expertise/red-team.md` |
| S-PROF | Performance profiling | `.positai/expertise/profiling.md` |
| S-COORD | Coordination review | `.positai/expertise/coordination.md` |

Priority: P3 when ≥6 OPEN tasks, P2 when 3–5, P1 when <3.

### Key files

| File | Purpose |
|------|---------|
| `a.XXX` | Individual Shiny bug reports (agents triage → `to-do.md`, then delete) |
| `issues.md` | Human-entered issues (agents triage → `to-do.md`) |
| `to-do.md` | Task queue (active/open tasks only) |
| `completed-tasks.md` | Archive of completed tasks |
| `coordination.md` | Strategic plan |
| `agent-X.md` | Agent progress log |
| `AGENTS.md` | Conventions + architecture reference |
| `.positai/expertise/*.md` | Standing task methodology |

## Test file conventions

All `tests/testthat/test-ts-*.R` files must use `TreeSearch:::` to call
internal C++ bridge functions. Define short local wrappers for readability.

Shared helpers are in `tests/testthat/helper-ts.R` (`make_ts_data()`,
`ts_score()`, `validate_result()`, `skip_extended()`).

**Never use `%in%` on Splits objects in test files** — S3 dispatch fails
in the cloned namespace created by `test_check()`. Use `as.logical()`
matrix comparison instead.

### Test tiering

Every new `test-ts-*.R` file must be assigned to one of three tiers.
See `tests/testing-strategy.md` for the full rationale.

| Tier | Guard | When it runs | Use for |
|------|-------|-------------|---------|
| 1 — CRAN | none | always (CRAN + CI + local) | Fast (< ~2 s) API and data-structure unit tests |
| 2 — CI | `skip_on_cran()` at **file level** (first executable line) | CI + local | C++ engine correctness, scoring, search algorithms |
| 3 — Extended | `skip_extended()` at **file level** | `TREESEARCH_EXTENDED_TESTS=true` only | Stress tests, benchmarks, timing measurements |

**Default for new `test-ts-*` files: Tier 2.** Add `skip_on_cran()` as the
very first executable line (before any helpers or `test_that()` calls):

```r
# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
```

Use Tier 3 only for tests that take > ~10 s or are sensitive to machine load.

## R source file ordering

`DESCRIPTION` has an explicit `Collate:` field. When adding a new `.R` file,
**update the Collate field** — otherwise R sources alphabetically, which can
break if one file's top-level code depends on a later file.

## Documentation checks (mandatory)

After any change to a function signature or roxygen block, run:

```r
devtools::check_man()
```

After writing or updating documentation prose, also run:

```r
spelling::spell_check_package()
```

Both should be clean before committing. These are fast and catch issues
(`check_man` catches Rd parse errors, cross-ref failures, `\usage` mismatches;
`spell_check_package` catches typos in `@description`/`@details`/`@param` text).

References are added using Rdpack's \insertCite{}, with
\insertAllCited{} in the references section.

## Algorithm vignette (mandatory updates)

`vignettes/search-algorithm.Rmd` documents the search algorithm for
publication. **Any change that modifies search behaviour** — new heuristics,
parameter tuning, scoring methods, stopping criteria, pool management, or
rearrangement operators — **must be accompanied by an update to this
vignette.**

- Published techniques: add a short summary and `@Key` citation.
- Novel contributions: describe the algorithm in enough detail for a reader
  to understand the design and rationale. Include empirical results where
  available (e.g. benchmark deltas).
- New references: add `@article{Key, ...}` to `inst/REFERENCES.bib`.

The vignette uses pandoc-style `@Key` citations (same as the other
vignettes), not Rdpack `\insertCite{}`.

## Architecture reference

### R-level API

| Function | Engine | Purpose |
|----------|--------|---------|
| `MaximizeParsimony()` | C++ driven search | Primary search (EW, IW, profile, constraints) |
| `Morphy()` | R-loop + MorphyLib | Legacy search (custom stopping, per-iteration callbacks) |
| `MaximizeParsimony2()` | — | Deprecated alias for `MaximizeParsimony()` |
| `Resample()` | C++ | Jackknife/bootstrap resampling |
| `SuccessiveApproximations()` | C++ | Successive approximations weighting |
| `TreeLength()` | C++ `ts_fitch_score` | Score one or more trees |
| `FastCharacterLength()` | C++ `ts_char_steps` | Per-character step counts |
| `AdditionTree()` | C++ `ts_wagner_tree` | Wagner tree construction |
| `RandomTreeScore()` | C++ (phyDat) or MorphyLib (morphyPtr) | Score a random tree |
| `TaxonInfluence()` | C++ via `MaximizeParsimony()` | Per-taxon search |
| `SearchControl()` | — | Expert parameter constructor for `MaximizeParsimony()` |
| `ParsSim()` | Pure R | Simulate datasets under parsimony (EW/IW/profile) |

`MaximizeParsimony()` has a backward-compatibility shim: passing old
Morphy-style parameters (`ratchIter`, `tbrIter`, etc.) triggers a deprecation
warning and delegates to `Morphy()`. Scheduled for removal in 2028.

### Driven search pipeline per replicate

1. Random Wagner tree → optional SPR → TBR to local optimum
2. XSS sectorial search (if tree large enough)
3. RSS random sectorial search
4. CSS constrained sectorial search
5. Ratchet perturbation to escape local optima
6. Drift search (accept suboptimal moves)
7. Final TBR polish
8. Add to pool
9. Fuse against pool (every `fuse_interval` replicates)

Post-search: TBR plateau enumeration from all pool seeds to find MPTs.

### Strategy presets (auto-selected by `NTip` and signal density)

| Preset | Condition | Key settings |
|--------|-----------|-------------|
| sprint | ≤30 tips | 3 ratchet (4%), 0 drift, XSS only, NNI-first, consensus-stop 3 |
| default | 31–64 tips; or ≥65 tips with <100 char patterns | 12 ratchet (25%, 5 moves), 2 drift (AFD 5, RFD 0.15), XSS+RSS, consensus-stop 3, Wagner×3, NNI-first, adaptive level |
| thorough | 65–119 tips with ≥100 char patterns | 20 ratchet (25%, 5 moves, adaptive), 5 NNI-perturb, 12 drift (AFD 5, RFD 0.15), XSS+RSS+CSS, consensus-stop 3, Wagner×3, NNI-first, outerCycles=2 |
| large | ≥120 tips with ≥100 char patterns | 12 ratchet (25%, 5 moves, adaptive), 0 NNI-perturb, 0 drift, 1 SA cycle (T=20→0, 5 phases), XSS(3)+RSS(2)+CSS(1), consensus-stop 2, Wagner×1 biased (Goloboff 2014), NNI-first, outerCycles=1, tbrMaxHits=1, sectorMaxSize=100 |

`sprint`/`default`/`thorough` set `consensusStableReps = 3`; `large` sets
`consensusStableReps = 2` (faster convergence detection since replicates at
120+ tips take 30–90s each).

**Large preset design rationale (T-179, 2026-03-24):** At 180 tips, each TBR
convergence takes ~5–7s, making phases like NNI-perturbation (~5.5s/cycle) and
drift (~4s/cycle) extremely expensive. Systematic benchmarking on mbank_X30754
(180t, 418p) showed that reducing cycle counts (12 ratchet, 4 drift, no NNI-perturb)
with outerCycles=1 and a single biased Wagner start outperforms the thorough
preset by 4–7 steps (median) at 30–60s budgets and ties at 120s, while
consistently completing more replicates.

**Post-T-206 Hamilton HPC baselines (2026-03-26, EPYC 7702, 5 seeds):**
30s median=1202 (range 1189–1214), 60s median=1190 (1190–1202), 120s
median=1185 (1171–1189). Per-replicate median 17.3s (cf. ~60s pre-T-206).
The 65–74 step improvement over pre-T-206 Intel baselines is primarily
from the outer cycle reset cap (maxOuterResets=0), not hardware.
Phase distribution: TBR 43.6%, Ratchet 32.2%, SA 7.4% (14% hit rate,
0.8 steps/s — least productive phase). T-248 benchmarked annealCycles
0/1/3: AC=1 (400ms/rep, 40% hit rate) is most cost-effective; AC=3
(1370ms/rep, 21% hit rate) showed no significant score gain (p>0.5,
n=5 seeds). Large preset reduced to annealCycles=1.

All presets set `nniFirst = TRUE` (NNI warmup before TBR) and
`sprFirst = FALSE` (SPR is counterproductive when NNI is active —
empirically NNI→TBR outperforms NNI→SPR→TBR). With `nniFirst`, each
Wagner start is NNI-optimized before selection (best of 3 NNI-local optima
rather than 3 raw Wagner scores). `default` also enables `adaptiveLevel =
TRUE` (scale ratchet/drift by hit rate); `thorough` omits it because high
base cycle counts already cover hard landscapes.

**Ratchet perturbation tuning (2026-03-22)**: Systematic profiling across
all 14 benchmark datasets showed the previous 4% perturbation probability
was far too gentle. With 253 characters (Zhu2013), 4% zeroes only ~10
characters — insufficient to reshape the landscape. Increasing to 25%
with fewer perturbed TBR moves (5 instead of auto=20) improves median
scores by 3–7 steps on hard datasets while completing fewer but more
productive replicates. 9/14 datasets improved, 4 unchanged, 1 marginal at
10s budget (resolves at 20s). The key insight: the perturbed-phase TBR
should be short (the landscape is warped, so extensive search on it is
wasteful), but the perturbation itself should be aggressive enough to
meaningfully displace the tree from its current basin of attraction.

Signal-density gate: datasets with few character patterns (<100) have flat
parsimony landscapes where intensive search adds no benefit.

### Adaptive sectorial search

XSS and CSS use **adaptive early-exit**: after each round of sector searches
+ global TBR polish, if the overall best score did not improve, remaining
rounds are skipped. This avoids wasting ~7% of replicate time on datasets
where sectorial search is unproductive (e.g. Dikow2009). On productive
datasets (e.g. Zhu2013), the early exit never fires.

### Conflict-guided RSS

RSS uses **conflict-guided sector selection**: before each replicate's RSS
phase, `driven_search()` computes a `SplitFrequencyTable` from the pool's
best-score trees. Within `rss_search()`, each internal node's "conflict
score" is `1 − (fraction of pool trees containing that split)`.
Max-descendant conflict is propagated upward, and eligible sector roots
are sampled via `std::discrete_distribution` with weight `1 + 3 × conflict`.
Falls back to uniform selection when the pool has <2 best-score trees or
when conflict variation is negligible.

### Consensus-stability stopping

After each replicate, if `consensus_stable_reps > 0` (default 3 in all
presets), the pool's strict consensus hash is compared to the previous
replicate's. If unchanged for `consensus_stable_reps` consecutive
replicates, the search terminates early. `compute_consensus_hash()` uses
XOR of per-split FNV-1a hashes for O(pool × splits) cost.

### Adaptive search level

When `adaptive_level = true`, ratchet and drift cycle counts are scaled
each replicate based on the cumulative hit rate:
- hit_rate > 0.7 → 0.5× (easy landscape)
- hit_rate > 0.4 → 0.75×
- hit_rate < 0.15 → 1.5× (hard landscape)
- else → 1.0×

### TBR zero-length clip skipping + regraft merging (collapsed flags)

`compute_collapsed_flags()` (`ts_collapsed.h/.cpp`) identifies edges where
clipping provably cannot improve score. Checks 5 conditions: (1) zero
standard-block cost at parent, (2) zero NA-block cost at parent, (3) prelim
preservation (`prelim[sibling] == prelim[parent]`), (4) down2 preservation
(NA), (5) subtree_actives preservation (NA). Works for EW, IW, Profile,
and NA-aware scoring. Integrated into TBR, SPR, and drift search.
Disabled during MPT enumeration (equal-score topologies may exist).
Recomputed after every accepted move.

**Regraft merging** (Goloboff 1996): within a collapsed region (connected
set of nodes linked by zero-length edges), all regraft positions yield the
same full score. Only boundary edges (entering the region) are evaluated;
interior collapsed edges are skipped via `if (collapsed[below]) continue`.
TBR, SPR, and drift all use this. The `CollapsedRegions` struct exists in
the header but callers use `compute_collapsed_flags()` directly (the
`region_id` field is unused — only the boolean flag array matters).

**Collapsed-topology pool dedup**: `compute_collapsed_splits()` in
`ts_splits.cpp` produces the split set excluding collapsed edges. Two
binary trees differing only in zero-length resolutions produce the same
collapsed split set → treated as duplicates by `TreePool::add_collapsed()`.
Both serial (`driven_search`) and parallel (`ThreadSafePool`) paths use
collapsed dedup.

**Benchmark results** (2026-03-22, 4 standard datasets, 3 seeds each):
Skip rate = 0% on all datasets (Vinther2008 23t, Agnarsson2004 62t,
Zhu2013 75t, Dikow2009 88t). Near-optimal trees in these morphological
datasets have negligible zero-length edges. Overhead from flag computation
is negligible. Score equivalence confirmed (enabled vs disabled produce
identical best scores). Benefit expected on sparse/synthetic data.

### C++ module map

| Module | Header/Source | Purpose |
|--------|--------------|---------|
| Fitch scoring | `ts_fitch.h/.cpp` | Downpass, uppass, incremental, indirect |
| NA scoring | `ts_fitch_na.h` | Three-pass inapplicable algorithm (Brazeau et al. 2019) |
| NA incremental | `ts_fitch_na_incr.h` | Incremental NA-aware scoring for TBR/drift |
| SIMD | `ts_simd.h` | SSE2/NEON portability layer for bit-parallel ops |
| Data | `ts_data.h/.cpp` | `DataSet`, `CharBlock`, `build_dataset`, simplification |
| Tree | `ts_tree.h/.cpp` | `TreeState`, topology manipulation, `PreallocUndo` |
| Constraint | `ts_constraint.h/.cpp` | Topological constraint enforcement |
| TBR | `ts_tbr.h/.cpp` | TBR search (with sector_mask for CSS) |
| SPR/NNI | `ts_search.h/.cpp` | SPR and NNI search (standalone, not in driven pipeline) |
| Ratchet | `ts_ratchet.h/.cpp` | Perturbation (zero/upweight/mixed, adaptive) |
| Drift | `ts_drift.h/.cpp` | Accept suboptimal moves within AFD/RFD limits |
| Wagner | `ts_wagner.h/.cpp` | Greedy addition tree (incremental scoring, NA-aware) |
| Sectorial | `ts_sector.h/.cpp` | RSS (conflict-guided), XSS, CSS; from-above HTU |
| Fuse | `ts_fuse.h/.cpp` | Tree fusing (in-place exchange) |
| Pool | `ts_pool.h/.cpp` | Dedup, eviction, consensus hash, split frequency table |
| Splits | `ts_splits.h/.cpp` | Bipartition computation, comparison, `hash_single_split()` |
| Driven | `ts_driven.h/.cpp` | Multi-replicate orchestrator |
| Resample | `ts_resample.h/.cpp` | Jackknife, bootstrap, successive approximations |
| Parallel | `ts_parallel.h/.cpp` | `std::thread` inter-replicate parallelism |
| RNG | `ts_rng.h/.cpp` | Thread-safe RNG (`thread_local` dispatch) |
| Simplify | `ts_simplify.h/.cpp` | Character compression and uninformativeness checks |
| Collapsed | `ts_collapsed.h/.cpp` | Zero-length edge detection for clip skipping |
| NNI perturb | `ts_nni_perturb.h/.cpp` | Stochastic NNI-perturbation (IQ-TREE-style topology escape) |
| Rcpp bridge | `ts_rcpp.cpp` | All Rcpp-exported functions |

### Scoring modes

`ScoringMode` enum in `ts_data.h`: `EW`, `IW`, `PROFILE`.
- **EW**: standard Fitch parsimony
- **IW**: implied weights via `e/(k+e)` where `e = steps - min_steps`
- **PROFILE**: lookup in `info_amounts` table (structurally identical to IW pipeline)

Profile mode sets `ds.concavity = 1.0` (finite sentinel) so existing
`isfinite()` checks activate the weighted pipeline without code duplication.

### Parallelism design

- `std::thread` (not OpenMP) to avoid R memory allocator conflicts
- Per-thread: `DataSet` copy, `ConstraintData` copy, `std::mt19937` RNG
- Shared: `ThreadSafePool` (mutex-guarded), atomic stop flag
- Main thread: pre-generates seeds from R's RNG, polls
  `R_CheckUserInterrupt()` and timeout every 200ms
- Worker threads make no R API calls — `ts_rng.h` provides `thread_local`
  dispatch (null → R API for serial; set → thread-local for parallel)

### Scoring notes

- `.h` file changes (`ts_fitch_na.h`, `ts_fitch_na_incr.h`) may require
  `touch src/ts_fitch.cpp` before rebuild if the build system doesn't track
  header dependencies.
- Incremental scoring is a **screening heuristic** for candidate selection;
  `full_rescore()` / `score_tree()` is always authoritative.
- See `.positai/expertise/fitch-scoring.md` for detailed invariants:
  uppass correctness proof, NA staleness analysis, `upweight_mask` audit.

### Constraint enforcement

- `build_constraint()` reads R split matrix with **column-major** indexing:
  `split_matrix[s + n_splits * t]`.
- Wagner uses LCA-based constraint mapping (`wagner_map_constraint_nodes`)
  since splits aren't fully present during incremental construction.
- Wagner has a posthoc retry loop (up to 100 random addition orders) as a
  safety net for edge cases.

## Exported Rcpp functions

All registered in `ts_rcpp.cpp` and `TreeSearch-init.c`. Run
`Rscript check_init.R` to verify consistency.

| Function | Module | Purpose |
|----------|--------|---------|
| `ts_fitch_score` | ts_fitch | Score a tree |
| `ts_char_steps` | ts_rcpp | Per-pattern step counts (with simplification offsets) |
| `ts_na_debug_char` | ts_fitch_na | Per-node debug for a single pattern |
| `ts_na_char_steps` | ts_fitch_na | Per-pattern step counts (raw, no offsets) |
| `ts_debug_clip` | ts_fitch | Debug SPR clip/regraft |
| `ts_test_indirect` | ts_fitch | Debug indirect length |
| `ts_nni_search` | ts_search | NNI hill-climbing |
| `ts_spr_search` | ts_search | SPR hill-climbing |
| `ts_tbr_search` | ts_tbr | TBR with plateau exploration |
| `ts_ratchet_search` | ts_ratchet | Ratchet perturbation |
| `ts_drift_search` | ts_drift | Drift search |
| `ts_wagner_tree` | ts_wagner | Wagner tree (specified addition order) |
| `ts_random_wagner_tree` | ts_wagner | Wagner tree (random order) |
| `ts_compute_splits` | ts_splits | Bipartition splits from edge matrix |
| `ts_trees_equal` | ts_splits | Compare two trees |
| `ts_pool_test` | ts_pool | Pool deduplication test |
| `ts_tree_fuse` | ts_fuse | Fuse two trees |
| `ts_sector_diag` | ts_sector | Sectorial search diagnostics |
| `ts_rss_search` | ts_sector | Random Sectorial Search |
| `ts_xss_search` | ts_sector | Exclusive Sectorial Search |
| `ts_driven_search` | ts_driven | Full driven search |
| `ts_resample_search` | ts_resample | One jackknife/bootstrap replicate |
| `ts_successive_approx` | ts_resample | Successive approximations |
| `ts_parallel_resample` | ts_parallel | Batch resample with parallelism |
| `ts_bench_tbr_phases` | ts_rcpp | TBR phase timing diagnostic |

## MorphyLib deprecation status

Migration plan in `inst/deprecation/morphy-migration.md`.

**Already migrated to C++:** `MaximizeParsimony`, `AdditionTree`, `Resample`,
`SuccessiveApproximations`, `TreeLength`, `CharacterLength`,
`FastCharacterLength`, `RandomTreeScore`, `TaxonInfluence`.

**Still using MorphyLib:** Legacy search functions (`Ratchet`, `Jackknife`,
`MorphyBootstrap`, `CustomSearch`), R-level tree rearrangement functions.
These are candidates for deprecation rather than migration.

## Shiny app (`inst/Parsimony/`)

Fully modularized from monolithic `app.R` into Shiny modules:
- `global.R` — library calls, constants, helpers, colours, citations, module UI instantiation
- `ui.R` — `fluidPage(...)` definition using module UI elements
- `server.R` — `AppState()` + module wiring + `ShowConfigs` observer + `onStop()`
- `server/app_state.R` — `AppState()` typed `reactiveValues()` constructor
- `server/logging.R` — session logging infrastructure
- `server/mod_*.R` — 7 Shiny modules (`NS()`/`moduleServer()`)

**All server logic now lives in modules.** The old `events.R` has been
dissolved; its `ShowConfigs` function and `plotFormat` observer are inlined
in `server.R` (they operate on top-level DOM elements).

**Modules:**
- `mod_references.R` — references panel (no state)
- `mod_downloads.R` — all 8 download handlers
- `mod_data.R` — data loading + tree management (9 returned reactives).
  Uses `cb_ref` forward-reference env for circular deps with consensus module.
- `mod_clustering.R` — clustering analysis + tree distances (5 returned reactives)
- `mod_search.R` — search engine, scoring, weighting.
  Owns ExtendedTask, search config modal, result accumulation.
- `mod_treespace.R` — tree space visualization + plot settings (14 returned reactives)
- `mod_consensus.R` — consensus plotting, character mapping, stability/rogue analysis,
  concordance, cluster consensus, main plot dispatch, plot logging (1327 lines).
  Returns `MainPlot`, `RCode`, `UpdateKeepNTipsRange`,
  `UpdateDroppedTaxaDisplay`, `UpdateOutgroupInput`.

**Important:** Server source files are in `server/` NOT `R/`. Shiny 1.5+
auto-sources all `.R` files in an app's `R/` directory at startup (before
any session exists), which crashes on references to `output`/`input`/`session`.

Test suite: `NOT_CRAN=true` required for shinytest2 integration tests.
Run from `inst/Parsimony/`:
```bash
NOT_CRAN=true Rscript -e "testthat::test_dir('tests/testthat')"
```
`setup.R` loads `library(shinytest2)` for `AppDriver` availability.

**Important:** Integration tests trigger `pkgbuild::compile_dll(debug=TRUE)`
via `load_all()`. `src/TreeSearch-win.def` prevents linker failures from
corrupted auto-generated `tmp.def` on Windows.

Module tests: `test-mod-references.R` (4), `test-mod-data.R` (9),
`test-mod-clustering.R` (12), `test-mod-treespace.R` (5),
`test-mod-downloads.R` (11), `test-mod-search.R` (28),
`test-mod-consensus.R` (9).
Integration tests: `test-app-smoke.R` (3), `test-Distribution.R` (13),
`test-SearchLog.R` (4), `test-ViewChars.R` (12). Total: 110 assertions.

## Version and CRAN status

- **Version**: 2.0.0 (major bump for new `MaximizeParsimony()` API)
- **R CMD check**: 0 ERRORs, 0 WARNINGs, 1 NOTE (R 4.5.2 internal bug)
- **Test suite**: ~9200 R-level + 1859 ts-* + 128 ParsSim + 37 MaddisonSlatkin + 49 recode-hierarchy pass

## Key design decisions (reference)

1. **PreallocUndo** (`ts_tree.h`): Pre-allocated flat buffers for TBR/drift
   undo stack. Uses `grow()` to dynamically expand when capacity exceeded
   (NA uppass saves both internal nodes and tips). Initial capacity `3 * n_node`.

2. **TBR symmetry breaking** (`ts_tbr.cpp`): FNV-1a hash deduplication of
   `virtual_prelim` vectors to skip redundant rerooting evaluations.

3. **Bounded indirect scoring**: All search modules use `_bounded` variants
   that bail out when accumulated score exceeds best candidate.

4. **Profile parsimony**: Reuses IW indirect pipeline unchanged; only delta
   precomputation differs. `ds.concavity = 1.0` sentinel activates weighted
   path. Max 2 informative states per character; inapplicable → ambiguous.

5. **MPT enumeration**: Post-search TBR plateau walk from all pool seeds.
   `tbr_search()` accepts optional `TreePool* collect_pool` parameter.

6. **All-ambiguous phyDat guard**: `TreeLength()` and `MaximizeParsimony()`
   check for `levels = NULL` / 0-column contrast matrix before calling C++.

7. **From-above HTU for sectorial search** (`ts_sector.cpp`):
   `compute_from_above_for_sector()` computes `from_above[sector_root]` —
   the Fitch state-set the rest of the tree sends *down* to the sector
   boundary, excluding the sector's own contribution. Used instead of
   `final_[parent]` in `build_reduced_dataset()`. O(depth × total_words).

8. **Split frequency table** (`ts_pool.h/.cpp`): `SplitFrequencyTable` maps
   per-split FNV-1a hash → occurrence count across best-score pool trees.
   Used by conflict-guided RSS to weight sector selection. The same FNV-1a
   hash (`hash_single_split()` in `ts_splits.h`) is used by consensus
   hashing and split frequency counting — must stay consistent.

9. **Consensus-stability hash** (`ts_pool.cpp`): XOR of FNV-1a hashes of
   splits present in ALL best-score trees. Updated after each replicate.
   Hash collision false-matches are conservative (over-count stability).

10. **Diversity-aware pool eviction** (`ts_pool.cpp`): When the pool is full
    and a new tree ties the worst score, the entry most similar to the new
    tree (most shared splits, counted via per-split FNV-1a hash set
    membership) is evicted. This maintains topological diversity in the pool,
    improving fusing effectiveness. Falls back to arbitrary worst entry when
    the new tree is strictly better.

11. **Cross-replicate consensus constraint tightening** (`ts_driven.cpp`):
    When `consensus_constrain = true` and no user constraint is supplied,
    after ≥5 replicates, unanimous pool splits are extracted and enforced
    as topological constraints via `build_constraint_from_bitsets()`. The
    TBR/SPR search then avoids breaking established consensus clades.
    Constraints are cleared and rebuilt whenever the best score changes.
    Sector/fuse operations do not enforce auto-constraints (no posthoc
    DataSet is built).

## Alternative inapplicable-handling algorithms (in progress)

Plan: `.positai/plans/2026-03-19-0643-alternative-inapplicable-handling-algorithms.md`

Adding HSJ (Hopkins & St. John 2021) and step-matrix/x-transformation
(Goloboff et al. 2021) scoring as alternatives to the existing Brazeau
et al. (2019) three-pass algorithm. Both require an explicit character
hierarchy specification.

### New files

| File | Purpose | Status |
|------|---------|--------|
| `R/CharacterHierarchy.R` | `CharacterHierarchy` S3 class, `validate_hierarchy()`, `hierarchy_from_names()`, `hierarchy_chars()`, `hierarchy_controlling()`, `non_hierarchy_weights()` | Complete, 34 tests passing |
| `tests/testthat/test-CharacterHierarchy.R` | Unit tests for hierarchy specification + weight partitioning | Complete |
| `src/ts_hsj.h` | `HierarchyBlock` struct (with `absent_state`), `hsj_score()` declaration, `partition_weights()` | Complete |
| `src/ts_hsj.cpp` | `partition_weights()`, `fitch_label_char()` (with uppass), `score_hierarchy_block()`, `hsj_score()` | Complete (full-rescore only; not wired to search pipeline) |
| `src/ts_sankoff.h` | `SankoffChar`, `SankoffData` structs, `sankoff_score()`, `sankoff_score_char()`, `sankoff_uppass()` | Complete |
| `src/ts_sankoff.cpp` | Sankoff downpass, uppass, root forcing | Complete |
| `R/recode_hierarchy.R` | `recode_hierarchy()`: x-transformation recoding (Goloboff et al. 2021) | Complete, 49 tests |
| `tests/testthat/test-recode-hierarchy.R` | Unit tests for recode_hierarchy() | Complete |
| `inst/REFERENCES.bib` | Added `Goloboff2021b` entry | Complete |

### Modified files

| File | Change |
|------|--------|
| `DESCRIPTION` | Added `CharacterHierarchy.R` to Collate field |
| `R/MaximizeParsimony.R` | Added `hierarchy`, `inapplicable`, `hsj_alpha` params with validation; non-brazeau methods currently `stop()` with "not yet implemented" |
| `src/ts_data.h` | Added `inapp_state` field to `DataSet` (for HSJ) |
| `src/ts_data.cpp` | Populate `inapp_state` in `build_dataset()` |

### Design decisions

- `hierarchy` is a **separate argument** to `MaximizeParsimony()` (not a phyDat attribute)
- `inapplicable` and `hsj_alpha` are **top-level args** alongside `concavity`
- Default `hsj_alpha = 1.0`
- IW + hierarchy and Profile + hierarchy: **deferred**
- Constraint interaction: **ignored** for now
- Resampling: **hierarchical** — resample top-level chars; when a controlling primary is sampled, also resample within its block; recurse for nested hierarchies

### Resampling with hierarchy (T-124)

`Resample()` now accepts `hierarchy`, `inapplicable`, and `hsj_alpha`
parameters. When `inapplicable != "brazeau"`, resampling is hierarchy-aware:

- **Resampling units**: each non-hierarchy character = 1 unit; each
  top-level hierarchy block (primary + all dependents) = 1 atomic unit.
- **Jackknife**: retain `proportion` of units without replacement.
- **Bootstrap**: sample `n_units` units with replacement (blocks can be
  duplicated).
- Per replicate: `.HierarchicalResampleWeights()` computes pattern weights
  for non-hierarchy chars and per-block sample counts. `.ResampleHierarchy()`
  calls `ts_driven_search` per replicate with filtered HSJ blocks or xform
  chars.
- **No C++ changes**: reuses existing `ts_driven_search` HSJ/xform infrastructure.
- **Parallelism**: serial R loop over replicates (C++ inter-search parallelism
  via `nThreads` still available within each replicate). Adding C++-level
  inter-replicate parallelism is a future optimization.

### Remaining work (Phase 1c–f)

1. ~~Pass `absent_token`, `n_tokens` from R to C++~~ **Done** (T-115): `absent_state` in HierarchyBlock, `inapp_state` in DataSet.
2. ~~Partition original characters into hierarchy vs non-hierarchy sets~~ **Done** (T-115): `partition_weights()` (C++) and `non_hierarchy_weights()` (R).
3. ~~Implement `hsj_score()` core algorithm in `ts_hsj.cpp`~~ **Done** (T-116): `fitch_label_char()` + `score_hierarchy_block()` + `hsj_score()`.
4. ~~Add Rcpp bridge function for HSJ scoring in `ts_rcpp.cpp`~~ **Done** (T-116): `ts_hsj_score()` registered in init.c.
5. ~~R-side marshalling~~ **Done** (T-116): `build_tip_labels()`, `hierarchy_to_blocks()`.
6. ~~Remove placeholder `stop()`~~ **Done** (T-117): HSJ wired into `score_tree()` dispatch, `ts_driven_search()` bridge, `MaximizeParsimony()`. End-to-end test against paper examples (T-118)

### Key algorithm notes (HSJ)

- Paper's Algorithm 1 initializes `a(l) = p(l) = 0` for all leaves. This is
  incorrect for enforcing observed leaf states. Correct initialization:
  leaf with primary absent → `a(l) = 0, p(l) = INF`; primary present →
  `a(l) = INF, p(l) = 0`. Verified against hand-computed example.
- `score_hierarchy_block()` operates per hierarchy block. Non-hierarchy
  characters use standard Fitch. Total = Fitch(non-hierarchy) + Σ HSJ(blocks).
- Secondary character labels at internal nodes from Fitch first-pass
  (inapplicable treated as a separate state).
- HSJ is full-rescore only (no incremental variant). Performance mitigation:
  candidate screening via Fitch, full HSJ only for promising candidates.

### Phase 2 (step-matrix) — Complete (end-to-end functional)

Sankoff engine (`ts_sankoff.h/.cpp`) implements downpass, uppass, root forcing.
R-level `recode_hierarchy()` combines primary + secondaries into composite
step-matrix character with asymmetric costs (gain:loss = n+1:1). Multistate
secondaries supported (state count = ∏k_i + 1). Nested hierarchies deferred.
Integration complete: `ScoringMode::XFORM` in `score_tree()` dispatches
Fitch(non-hierarchy) + Sankoff(recoded). `MaximizeParsimony()` accepts
`inapplicable = "xform"`. End-to-end search verified.

### Stochastic NNI-perturbation (T-186)

`ts_nni_perturb.h/.cpp` implements a topology-space escape mechanism inspired
by IQ-TREE's `doRandomNNIs()` (Nguyen et al. 2015). Complementary to the
weight-perturbation ratchet: the ratchet reshapes the objective function, while
NNI-perturbation directly displaces the tree topology.

**Algorithm:** Collect all internal NNI edges. For each edge (with probability
`perturb_fraction`, default 0.5), apply a random NNI swap — but skip edges
adjacent to already-swapped edges (two NNIs conflict if their edges share an
endpoint). Track touched nodes in a hash set. After all compatible swaps,
rebuild postorder and full rescore, then TBR to a new local optimum. Repeat
for `n_cycles`.

**Pipeline placement:** Between ratchet (phase 4) and drift (phase 5) in
`run_single_replicate()`. Disabled by default (`nniPerturbCycles = 0`).
Enabled in the `thorough` preset (5 cycles, 0.5 fraction).

**R API:** `SearchControl(nniPerturbCycles, nniPerturbFraction)`.
Timings reported as `nni_perturb_ms`.

### Biased Wagner addition (T-188, 2026-03-23)

`biased_wagner_tree()` (`ts_wagner.h/.cpp`) samples the taxon-addition
order from a softmax distribution weighted by informativeness score rather
than purely at random. Two criteria available:

- **GOLOBOFF** (bias=1): `score[t]` = number of non-ambiguous characters
  for taxon t. Ref: Goloboff 2014 (*Extended implied weighting*) §3.3.
- **ENTROPY** (bias=2): `score[t]` = Σ_c (n_states_c − |state set for t|).

**R API:** `SearchControl(wagnerBias = 0L, wagnerBiasTemp = 0.3)`.
Applied only to the first of `wagnerStarts` starts; remaining starts
use random order to preserve basin diversity.

**Benchmark results** (2026-03-23, 14 standard + crico-174):
- Wagner→TBR gap reduction: ~80% at 174t (random: 1356 steps, Goloboff: 244)
- Score improvement after TBR convergence: ~22 steps at 174t; 1–2 steps at ≤88t
- Anomalous slight regression at 75–100t (Giles2015, Zanol2014); T=0.3
  stochastic is safer than pure greedy (T=0)

### Outer search cycle loop (T-189, 2026-03-23)

`outer_cycles` in `SearchParams` / `outerCycles` in `SearchControl()`.
Wraps steps 3–6 of `run_single_replicate()` in a configurable outer
loop: [XSS+RSS+CSS → Ratchet → NNI-perturb → Drift → TBR] × N.
Ratchet/NNI-perturb/drift cycles are divided evenly among N outer cycles
(ceiling division, minimum 1 per cycle).

`outerCycles = 1` (default) is bit-for-bit identical to the previous
linear pipeline. `thorough` preset defaults to `outerCycles = 2`.

**Pattern:** Matches TNT's `xmult` interleaving (Goloboff 1999 §2.3):
after each ratchet/drift escape, a fresh XSS pass exploits the new
topology before the next perturbation round.

**Citation tracking:** see `papers.md` in project root for full
reference list used in optimization work.

### NNI in the driven pipeline

`nni_search()` in `ts_search.cpp` is implemented but **never called** in the
driven pipeline. At ≤88 tips, NNI is strictly redundant — TBR subsumes it and
completes in <1s per pass, so there's nothing to save.

**At 180 tips, NNI becomes essential.** TBR evaluates O(n²) candidates per
pass (358 clips × 356 regrafts × rerooting ≈ millions of evaluations),
and a single convergence from Wagner takes many minutes. NNI evaluates O(n)
candidates per pass (178 edges × 2 swaps = 356), roughly 1000× cheaper.
Most improvements during initial descent are NNI-reachable.

Proposed escalation: NNI → SPR → TBR, gated on `n_tip > ~100`.
See T-178 in `to-do.md`.

**Empirical comparison at 180 tips** (mbank_X30754, 3 seeds, EW):

| Strategy | Median score | Median time |
|----------|:-----------:|:-----------:|
| TBR alone | 1427 | 13.6s |
| SPR→TBR | 1360 | 13.1s |
| **NNI→TBR** | **1326** | **6.8s** |
| NNI→SPR→TBR | 1369 | 8.8s |

NNI→TBR wins on both score AND time (~2× faster, ~100 steps better than
TBR alone). SPR intermediate step adds time without benefit at this scale.
The NNI descent path leads TBR to better basins of attraction.

**Recommendation:** `nniFirst = TRUE` (always on — NNI costs ~1.5s at
180 tips, negligible at ≤88 tips). Replace `sprFirst` with `nniFirst`
for large trees (n_tip > ~80), or just always run NNI since the overhead
is negligible. SPR warmup is counterproductive at 180 tips.

**Metric note:** When comparing strategies, the right metric is
**time-adjusted expected best** — the expected minimum score from
k = budget / time_per_rep independent replicates, since multi-start search
keeps the best tree. The median measures typical quality, but a strategy
with high variance and occasional excellent finds can dominate if it gets
enough draws. Bootstrap estimation: sample k scores with replacement, take
the min, repeat 5000×, take the mean.

**Time-adjusted expected best (5 seeds, EW):**

| Budget | 88t: TBR | 88t: NNI→SPR→TBR | 180t: TBR | 180t: NNI→SPR→TBR |
|--------|:--------:|:-----------------:|:---------:|:-----------------:|
| 20s | 1617 | 1619 (+2) | 1388 | 1278 (−110) |
| 60s | 1617 | 1619 (+2) | 1348 | 1253 (−95) |
| 120s | 1617 | 1619 (+2) | 1337 | 1247 (−90) |

At ≤88 tips: NNI has a consistent but negligible 2-step penalty (within
noise of MPT enumeration). At 180 tips: NNI saves 90–110 steps. The
crossover is between 88 and 180 tips. No reactive per-run switching
needed — a simple always-on NNI warmup policy is optimal.

**Escalation strategy:** NNI→TBR (skip SPR) is simplest and nearly
optimal at all sizes. NNI→SPR→TBR adds ~7 steps at 180 tips for
negligible extra time, but adds complexity. SPR alone (without NNI)
is counterproductive at 180 tips (15s vs 7s, worse scores).

### Large-tree scaling issues (discovered 2026-03-23)

The 180-taxon `mbank_X30754` dataset (425 chars, 374 informative patterns,
40% missing, 20% inapplicable) exposed:

1. ~~**`maxTime` triggers Morphy delegation.**~~ **Fixed (T-184)**:
   `maxTime` is now intercepted before the Morphy shim check and
   mapped to `maxSeconds` with a deprecation warning. The C++ driven
   search is used correctly. (Note: initial 180-taxon benchmarking
   mistakenly used `maxTime`, which in an older version delegated to
   `Morphy()` — ~10× slower. All results since T-184 use `maxSeconds`.)
2. **C++ TBR convergence at 180 tips takes ~13s** (Wagner ~2560 → local
   optimum ~1420). NNI warmup (~1.5s) followed by TBR reduces this to
   ~7s while finding better scores. T-178 filed.
3. **Strategy presets assume replicate time O(seconds).** At 180 tips,
   a single replicate (TBR + XSS + ratchet + drift) takes ~60-100s.
   Cycle counts and fuse intervals need recalibration for large trees.

### Benchmarking methodology notes

**Early vs late search:** The character of the search changes over time.
Early replicates are dominated by initial descent quality (Wagner → local
optimum); late replicates test ratchet/drift escape effectiveness. At ≤88
tips, 20s gives 10–40 replicates spanning both regimes. At 180 tips, 20s
doesn't complete one replicate. A warm-start benchmark (T-180) would
isolate the escape-effectiveness question.

**Generalization to large trees:** All 14 existing benchmark datasets are
≤88 tips. Algorithmic choices (e.g. TBR vs NNI warmup, ratchet cycle counts)
that are optimal at 88 tips may be suboptimal at 180+. The 180-taxon dataset
should be added to the benchmark suite as a separate tier (T-181).

**`maxTime` confound (2026-03-23):** Initial 180-taxon testing used
`maxTime` (legacy Morphy parameter), which silently delegated to the
R-loop `Morphy()` engine. The C++ driven search (via `maxSeconds`) is
~10× faster at 180 tips. All subsequent profiling used the C++ path.

**180-taxon baseline (C++ driven search, EW, single replicate):**
- Wagner (best of 3): ~2560 steps, 16ms
- NNI convergence: ~1600 steps, 1.5s
- TBR convergence: ~1330 steps, 7s (from NNI-optimal start)
- XSS: additional ~60 steps improvement, 5s
- Total single replicate: ~25s (before ratchet/drift)

## Search optimization roadmap

Plan: `.positai/plans/2026-03-21-search-optimizations.md`

Ranked by priority:
1. ~~Consensus-guided sector targeting~~ — **Done**: RSS weighted by
   pool split conflict scores
2. ~~Diverse pool maintenance~~ — **Done**: evict most-similar entry on ties
3. ~~Cross-replicate constraint tightening~~ — **Done**: opt-in via
   `consensusConstrain = TRUE`
4. ~~Collapsed-tree clip skipping~~ — **Done**: zero-length edges skipped
   in TBR, SPR, and drift. Benchmark shows 0% skip rate on standard
   morphological datasets (Vinther2008, Agnarsson2004, Zhu2013, Dikow2009)
   because near-optimal trees have few zero-length edges. Negligible
   overhead. Benefit expected primarily on sparse/synthetic data.
   Full polytomy search remains post-2.0.0.
5. ~~Collapsed-region regraft merging + pool dedup~~ — **Done**: within
   collapsed regions (connected zero-length edges), only the boundary
   regraft position is evaluated (Goloboff 1996). Collapsed-topology
   pool dedup treats trees differing only in zero-length resolutions as
   duplicates. Parallel path also uses collapsed dedup. Diversity-aware
   pool eviction selects most-similar entry on ties.
6. ~~Strategy preset tuning~~ — **Done**: `default` preset now uses
   `wagnerStarts=3`, `sprFirst=TRUE`, `adaptiveLevel=TRUE`; `thorough`
   preset uses `sprFirst=TRUE`.
7. ~~Ratchet perturbation tuning~~ — **Done**: perturbation probability
   increased from 4% to 25%, perturbed TBR moves reduced from auto=20
   to 5, ratchet cycles increased from 5 to 10 (default) and kept at
   20 (thorough). Drift cycles increased from 2 to 4 with wider
   acceptance (AFD 5, RFD 0.15). Validated on all 14 datasets.
8. ~~Biased Wagner addition~~ — **Done** (T-188): `wagnerBias`
   (0=RANDOM, 1=GOLOBOFF, 2=ENTROPY) + `wagnerBiasTemp` in
   `SearchControl()`. First Wagner start uses biased addition order;
   remaining starts use random for basin diversity. Benchmarked on
   14 standard + crico-174 datasets: 80% Wagner→TBR gap reduction at
   174t; ~1–2 steps improvement at ≤88t (negligible). Goloboff 2014
   §3.3.
9. ~~Outer search cycle loop~~ — **Done** (T-189): `outerCycles` in
   `SearchControl()`. Wraps [XSS+RSS+CSS → Ratchet → NNI-perturb →
   Drift → TBR] in configurable outer loop; cycles divided evenly.
   `thorough` preset defaults to `outerCycles=2`. Matches TNT xmult
   pattern. Goloboff 1999 §2.3. Needs benchmarking to validate
   improvement on standard datasets.
10. ~~Drift MPT diversity experiment~~ — **Done** (T-254): Drift
    (`driftCycles=2`) provides zero score benefit, zero MPT enumeration
    benefit, zero topological diversity benefit, and costs 10–22% of
    replicates. On Wortley2006, no-drift consistently finds 4 MPTs vs
    1–3 with drift. On Geisler2001/Zhu2013, mean pairwise RF is
    identical (7.3 vs 7.4; 11.6 vs 10.2). Drift delays consensus
    stability without improving the answer. **Recommendation:** set
    `driftCycles=0` in default and thorough presets (T-255).

## Benchmarks and profiling

### MorphoBank external benchmark corpus

The neotrans repo (`../neotrans/inst/matrices/`) contains ~800 MorphoBank
NEXUS matrices. These complement the 14 bundled datasets and 1 large-tree
dataset for broader strategy validation.

**Catalogue:** `dev/benchmarks/mbank_catalogue.csv` (659 usable matrices
after ntax≥20 filter and dedup). Regenerate with
`Rscript dev/benchmarks/build_mbank_catalogue.R`.

**Train/validation split:** Matrices whose MorphoBank project number is
divisible by 5 are **validation** (124 matrices, ~19%). All others are
**training** (535 matrices). The 7 `syab*` files (non-MorphoBank, from a
Systematic Biology paper) are always training.

**Dedup:** Multi-file projects with ≥95% character identity on shared taxa
(≥80% taxon overlap) are flagged `dedup_drop = TRUE`. Greedy selection keeps
the largest matrix per redundancy cluster. 24 near-duplicates excluded.

**IMPORTANT:** Validation results must **never** be used to guide strategy
tuning. They confirm generalization only. This is a one-way door.

**Fixed 25-matrix training sample:** `MBANK_FIXED_SAMPLE` in
`bench_datasets.R` — 7 small, 7 medium, 7 large, 4 xlarge. Selected via
max-min distance on standardized features. Do not modify. Used by
`benchmark_mbank_sample()`.

**Key functions** (in `dev/benchmarks/bench_datasets.R`):
- `load_mbank_catalogue()` — loads metadata CSV (excludes dedup by default)
- `load_mbank_sample(cat, n, seed, split)` — stratified random sample
- `load_mbank_datasets(cat, keys)` — load specific matrices by key

**Benchmark runners** (in `dev/benchmarks/bench_framework.R`):
- `benchmark_mbank_sample()` — fixed 25-matrix training sample (routine)
- `benchmark_mbank_sweep(split)` — full training or validation sweep
- `benchmark_mbank_validation()` — validation sweep with prominent warning

**TNT comparison suite** lives in `../TS-TNT-bench/`. Key files:
- `dev/benchmarks/bench_tnt_compare.R` — runner (smoke/medium/full)
- `dev/benchmarks/tnt_comparison.qmd` — Quarto report (HTML output)
- `dev/benchmarks/.tnt-bench/` — staging dir for TNT I/O
- Requires TNT 1.6 at `C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe`

Benchmark scripts in `dev/benchmarks/`. Key files:
- `bench_regression.R` — CI regression test (score quality + timing bounds)
- `bench_framework.R` — Dataset × strategy × replicate grid
- `strategies.md` — Strategy space documentation

Profiling baselines in `.positai/expertise/profiling.md`. Phase distribution
(default strategy, EW, 2026-03-22, post-ratchet tuning):

| Phase | Zhu2013 (75t) | Dikow2009 (88t) | Agnarsson2004 (62t) |
|-------|:---:|:---:|:---:|
| Ratchet | 27% | 39% | 35% |
| Drift | 32% | 28% | 24% |
| TBR | 31% | 21% | 19% |
| XSS | 8% | 7% | 8% |
| RSS | 2% | 3% | 9% |
| Final TBR | 2% | 3% | 4% |

Per-candidate indirect scoring is at memory-throughput limit (~23 ns at
75 tips).

## VTune driver scripts — dry-run first

**Always test a VTune driver script with plain `Rscript` before launching
VTune.** Software-sampling overhead can be 5–20×; if the bare script takes
30 s, VTune may need 10 min. Target < 5 s bare run for a lite driver.

MaddisonSlatkin is exponential in tip count — even n=20 with k=3 can take
seconds per call. Use small n (≤ 15 for k=3, ≤ 12 for k=4, ≤ 9 for k=5)
and few iterations for VTune drivers.

### Ratchet tuning validation (2026-03-22)

Full 14-dataset comparison, optimized vs original defaults (10s budget,
3 seeds each). Median scores shown (lower is better).

> **Metric note:** Median per-replicate score is adequate for comparing
> parameter changes on a fixed pipeline (same time-per-rep). For comparing
> strategies with different time costs (e.g. NNI→TBR vs TBR), use
> **time-adjusted expected best** instead — see "Metric note" under
> "NNI in the driven pipeline".

| Dataset | Tips | Original | Optimized | Delta |
|---------|:---:|:---:|:---:|:---:|
| Longrich2010 | 20 | 131 | 131 | 0 |
| Vinther2008 | 23 | 79 | 79 | 0 |
| Sansom2010 | 23 | 189 | 189 | 0 |
| DeAssis2011 | 33 | 64 | 64 | 0 |
| Aria2015 | 35 | 143 | 143 | 0 |
| Wortley2006 | 37 | 494 | 491 | +3 |
| Griswold1999 | 43 | 408 | 407 | +1 |
| Schulze2007 | 52 | 165 | 164 | +1 |
| Eklund2004 | 54 | 442 | 441 | +1 |
| Agnarsson2004 | 62 | 778 | 778 | 0 |
| Zanol2014 | 74 | 1338 | 1331 | +7 |
| Zhu2013 | 75 | 649 | 650 | −1 |
| Giles2015 | 78 | 720 | 716 | +4 |
| Dikow2009 | 88 | 1614 | 1614 | 0 |

Zhu2013 marginal regression at 10s resolves at 20s (median 649→644).
At 20s with 5 seeds: Zhu2013 645/643, Giles2015 712/710, Dikow2009
1611/1611 (all improvements).
