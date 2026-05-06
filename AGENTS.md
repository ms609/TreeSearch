# TreeSearch Multi-Agent Development Notes

Always check the contents of `.AGENTS` for memories and policies relevant to
the task you have been assigned.

Update memory files with anything relevant you learn. But keep them lean.

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

### Coordination files live on `cpp-search` only

`to-do.md`, `u.nnn`, `agent-X.md`, `completed-tasks.md`, `coordination.md`,
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
     --title "<Letter>-nnn: <description>" --body "Agent <Letter>. ..."
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

> **Task IDs:** New tasks use `<Letter>-nnn` format (e.g. `A-042`), where
> `<Letter>` is your agent letter and `nnn` is your personal counter
> (tracked in `agent-<letter>.md`). Existing `T-nnn` IDs in `to-do.md`,
> `completed-tasks.md`, PRs, and git log are valid and need not be renamed.
> Before adding or removing rows in `to-do.md`, acquire the lock:
> `bash ../../todo-lock.sh . acquire` / `bash ../../todo-lock.sh . release`.

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
      report. Assign a `<Letter>-nnn` ID (your letter + incremented
      counter; see "Task IDs" in the master AGENTS.md), a priority based
      on severity (default P2), and tag as `[Shiny]`. Use the file's
      content as the description/notes.
   e. Delete the `aXXX.claimed-X.md` file once the `to-do.md` entry is
      written.
   f. Repeat for all claimed files before moving on. **Do not start
      working a task until all pending reports are triaged.**
3. **Triage user issues** (`u.*` files) before `to-do.md`. See the
   parent `../AGENTS.md` "User issue files" section for the full protocol
   (scan → claim via rename → triage → delete). While untriaged issues
   remain, triaging takes priority over `to-do.md` tasks (an issue may
   contain a P0).
4. **Check `remote-jobs.md`** for retrievable results. If a job is listed
   as complete (or past its expected duration), retrieve and process the
   results before claiming a new task.
5. If no untriaged issues or pending remote results, claim the next OPEN
   task from `to-do.md`.

Set `CONVERSATIONSUMMARY` to `Agent X: <task description>`.

> **Concurrency guard (u.nnn / a.XXX):** Atomic rename (`mv u.001
> u.001.claimed-X`) ensures exactly one agent wins each file. NTFS
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
   `| <Letter>-nnn | Short description | X | Brief notes |`
3. Set `agent-X.md` to IDLE.
4. Append a brief entry to this file documenting what changed.
5. Update `coordination.md` if strategic objectives are affected.
6. Take next task.

### Shiny bug report intake

The human files Shiny app bugs as individual files in the project root.
Naming convention:

`u.010`, `u.011`, … (dot-separated, no extension)

Each file contains a free-text bug description. The human's workflow is:
create file → write bug → save → never touch the file again.

**Agent responsibility:** Triage all pending bug report files into
`to-do.md` at the start of every `/assign` (step 2 in Assignment above).

**Claim protocol:**
```bash
# List unclaimed reports (both conventions)
ls u.[0-9]* 2>/dev/null | grep -v 'claimed'

# Skip short files (< 20 chars) — don't rename, don't touch
wc -c < u.010  # check size first

# Claim atomically (rename)
mv u.010 u.010.claimed-X

# Read, triage into to-do.md, then delete
cat u.001.claimed-X
# ... create to-do.md entry ...
rm u.001.claimed-X
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
| `u.nnn` | User issue files (agents triage → `to-do.md`, then delete) |
| `to-do.md` | Task queue (active/open tasks only) |
| `remote-jobs.md` | Pending async jobs (Hamilton SLURM, long GHA) — check at `/assign` |
| `completed-tasks.md` | Archive of completed tasks |
| `coordination.md` | Strategic plan |
| `agent-X.md` | Agent progress log |
| `AGENTS.md` | Conventions + workflow reference |
| `.positai/expertise/*.md` | Standing task methodology |

---

## User-level skills

The user has skills installed at `~/.positai/skills/`.
**Load these with the `skill` tool before starting relevant work:**

| Skill name | When to load |
|------------|-------------|
| `r-package-profiling` | Profiling, benchmarking, VTune, A/B comparison, hotspot analysis |
| `hamilton-hpc` | Hamilton HPC, SLURM jobs, SSH, remote benchmarking |

Example: `skill(skill: "hamilton-hpc")` before any Hamilton dispatch work.

---

## Mandatory checks

Run these before committing whenever the trigger applies:

| Trigger | Command |
|---------|---------|
| Function signature or roxygen block changed | `Rscript -e "devtools::check_man()"` |
| Documentation prose changed | `Rscript -e "spelling::spell_check_package()"` |
| `Rcpp::compileAttributes()` run | `Rscript check_init.R` (verifies `ts_rcpp.cpp` / `TreeSearch-init.c` arg counts) |
| Search behaviour changed (heuristics, scoring, stopping, pool) | Update `vignettes/search-algorithm.Rmd` |

Full details: `.AGENTS/memory/r-package-conventions.md`.

---

## MorphyLib deprecation status

Migration plan in `inst/deprecation/morphy-migration.md`.

**Already migrated to C++:** `MaximizeParsimony`, `AdditionTree`, `Resample`,
`SuccessiveApproximations`, `TreeLength`, `CharacterLength`,
`FastCharacterLength`, `RandomTreeScore`, `TaxonInfluence`.

**Still using MorphyLib:** Legacy search functions (`Ratchet`, `Jackknife`,
`MorphyBootstrap`, `CustomSearch`), R-level tree rearrangement functions.
These are candidates for deprecation rather than migration.

## Version and CRAN status

- **Version**: 2.0.0 (major bump for new `MaximizeParsimony()` API)
- **R CMD check**: 0 ERRORs, 0 WARNINGs, 1 NOTE (R 4.5.2 internal bug)
- **Test suite**: ~9200 R-level + 1859 ts-* + 128 ParsSim + 37 MaddisonSlatkin + 49 recode-hierarchy pass

---

## Technical reference

Load the relevant `.AGENTS/memory/` file before starting work in that area:

| Memory file | Load when... |
|-------------|--------------|
| `architecture.md` | Editing `src/ts_*.cpp`/`.h`, adding Rcpp exports, reviewing R-level API or key design decisions |
| `benchmarking.md` | Running benchmarks, doing VTune profiling, interpreting phase-distribution or Brazeau/Fitch results |
| `feature-inapplicable.md` | Working on HSJ, x-transform/Sankoff, `inapplicable=` parameter, or `CharacterHierarchy` |
| `r-package-conventions.md` | Adding `.R` files to `Collate:`, writing roxygen docs, updating vignettes |
| `search-algorithms.md` | Researching NNI warmup, biased Wagner, outer cycles, large-tree behaviour, or the search optimization history |
| `search_strategy.md` | Understanding the driven pipeline, strategy presets, adaptive search, collapsed-flag optimization |
| `shiny_app.md` | Working on `inst/Parsimony/`, Shiny modules, or app tests |
| `testing.md` | Adding or modifying `tests/testthat/test-ts-*.R`, choosing test tiers, writing helpers |
