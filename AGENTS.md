# TreeSearch Multi-Agent Development Notes

Always check the contents of `.AGENTS` for memories and policies relevant to
the task you have been assigned.

Update memory files with anything relevant you learn. But keep them lean.

## Current phase: bug-fixing / pre-release (as of 2026-03-29)

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

Before dispatching, run `spelling::spell_check_package()` (or a targetted `spell_check_files()`).
GHA will fail on spelling errors.
If any "errors" can be avoided (e.g. by spelling out acronyms or wrapping in 
\acronym{}; by hyphenating compound words), reword. Add false positives 
to `inst/WORDLIST`.

Once confirmed, dispatch GHA with:

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
  R CMD INSTALL --library=.agent-<id> "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
```

Key points:
- `rm -f src/*.o src/*.dll` **must** precede every build — stale artifacts slow traversal and corrupt DLLs.
- Build into an agent-specific `$TMPBUILD` outside the source tree — avoids tarball collision when multiple agents build concurrently.
- `--no-resave-data` skips unnecessary `.rda` re-saving (not needed for dev installs).

Run **targeted** tests only:
```bash
Rscript -e "library(TreeSearch, lib.loc='.agent-<id>'); testthat::test_dir('tests/testthat', filter='test-ts-foo')"
```

**Never** use `R CMD INSTALL --library=.agent-<id> .` (in-place build).

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
Rscript -e ".libPaths(c('.agent-<id>', .libPaths())); roxygen2::roxygenise(load_code = roxygen2::load_installed)"
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
  R CMD INSTALL --library=.agent-<id> "$TMPBUILD"/TreeSearch_*.tar.gz && \
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

`to-do.md`, `u.nnn`, `completed-tasks.md`, `coordination.md`,
and `AGENTS.md` are **never committed on feature branches**. When a dispatched
agent working on a feature branch needs to claim a task or update coordination
files, they commit those changes directly to `cpp-search` (coordination-only
commit), keeping the feature branch clean.

To read coordination files while on a feature branch without switching:
```bash
git show cpp-search:to-do.md
git show cpp-search:coordination.md
```

### Shared files at merge time

`src/ts_rcpp.cpp` and `src/TreeSearch-init.c` use the existing append-only
convention — merge conflicts resolve cleanly by keeping both appended blocks.
`DESCRIPTION` (Collate field) and `NAMESPACE` require a manual merge pass;
this is expected and should be done carefully at feature-merge time.

### Feature branch lifecycle

1. `git checkout cpp-search && git checkout -b feature/<name>`
   Optionally create a worktree: `git worktree add ../worktrees/TS-<name> feature/<name>`
   **Never** switch the main `./TreeSearch` checkout away from `cpp-search` (or a
   feature branch actively being worked). Worktrees must always live under `../worktrees/`.
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
     --title "T-nnn: <description>" --body "Dispatched agent <id>. ..."
   ```
6. Set `to-do.md` status to `PR #N (<id>)`. Move on.
7. Human reviews and merges the PR.
8. After merge, clean up:
   ```bash
   git worktree remove ../worktrees/TS-<name>  # if worktree was used
   git branch -d feature/<name>
   git push origin --delete feature/<name>
   ```

---

## Multi-agent workflow protocol

> **Task IDs:** New tasks use `T-nnn` format. Existing `T-nnn`, `<Letter>-nnn`
> IDs in `to-do.md`, `completed-tasks.md`, PRs, and git log are valid and need
> not be renamed.
> Before adding or removing rows in `to-do.md`, acquire the lock:
> `bash ../../todo-lock.sh . acquire` / `bash ../../todo-lock.sh . release`.

### Dispatcher model

Agents are launched by the dispatcher (`dispatch.sh`) and receive an ephemeral
ID of the form `d1`, `d2`, etc. The dispatcher:

1. Reads `.dispatch/state.json` to determine which tasks are already in-flight.
2. Selects a task (via the Haiku ranker or an explicit task ID).
3. Mints a new agent ID and updates `to-do.md` to `ASSIGNED (d1)`.
4. Spawns a `claude -p` subprocess whose brief is loaded from
   `dev/dispatch/agent-brief.md`.
5. Logs output to `.dispatch/logs/<id>-<task>.log`.

**Agents do not edit `.dispatch/state.json` directly.** State is written only
by `dispatch.sh checkin`.

#### Starting a session

The user (or a parent dispatched session) calls:

```bash
bash dispatch.sh allocate <budget>      # e.g. 5%/5h, 2%/wk, 15m
bash dispatch.sh task <T-ID> [budget]   # explicit task; budget optional
```

#### Session start protocol (every dispatched agent, before claiming work)

1. **Resume check:** read your brief — if a resume action is recorded, execute
   it now.
2. **Triage user reports** (`a.*` and `u.*` files) — see "User report intake"
   below for the full claim protocol:
   a. List all `a.[0-9]*` **and** `u.[0-9]*` files in the project root
      (excluding any `*.claimed-*` files).
   b. For each file, check its size first. **Skip files shorter than
      20 characters** (likely mid-edit — the human may still be typing).
      Do not rename or touch these files; leave them for a later pass.
   c. For files ≥20 characters, claim atomically:
      `mv a.010 a.010.claimed-<id>` (or `mv u.010 u.010.claimed-<id>`).
      If the rename fails, another agent claimed it — skip.
   d. Create a `to-do.md` entry. **`a.*` files** → `### Shiny App`, tag
      `[Shiny]`. **`u.*` files** → section matching content (search bug,
      docs issue, etc.). Default priority P2; crash = P1, cosmetic = P3.
   e. Delete the `.claimed-<id>` file once the `to-do.md` entry is written.
   f. Repeat for all files before moving on. **Do not start working a
      task until all pending reports are triaged.** (An issue may be P0.)
3. **Check `remote-jobs.md`** for retrievable results. If a job is listed
   as complete (or past its expected duration), retrieve and process the
   results before claiming a new task.
4. If no untriaged issues or pending remote results, proceed with the
   assigned task.

> **Concurrency guard:** Atomic rename (`mv a.010 a.010.claimed-<id>` or
> `mv u.001 u.001.claimed-<id>`) ensures exactly one agent wins each file.
> NTFS rename is atomic; losers see "file not found" and skip.

### Worktree tasks

Tasks with status `WORKTREE (name)` are actively developed in a dedicated git
worktree under `C:/Users/pjjg18/GitHub/worktrees/` (e.g.
`../worktrees/TS-CID-cons`). **Do not claim or modify these tasks.** They are
reserved for the human developer working in that worktree. To mark a task as
in-flight on a worktree, set its status to `WORKTREE (name)` where *name*
matches the worktree directory basename.

> **Worktree rule:** Worktrees must **always** be created under `../worktrees/`
> (i.e. `C:/Users/pjjg18/GitHub/worktrees/<name>`). **Never** create a worktree
> directly inside `../` alongside the main checkout, and **never** switch the
> main `C:/Users/pjjg18/GitHub/TreeSearch` directory to a different branch using
> `git checkout` — it must remain on `cpp-search` (or the current feature branch
> being actively developed). Use a worktree instead.

### During work

- All work uses `.agent-<id>/` as library directory (e.g. `.agent-d1/`).
- **All builds, tests, and benchmarks in bash subprocesses** — never in the
  RStudio R session.
- **Use GHA for validation** (full test suites, R CMD check, benchmarks).
  Local builds are for targeted iteration only (build + run 1–2 test files).
  See "Validation workflow" section above.

### On task completion

1. **Delete** the task row from `to-do.md`. If the task was the last open
   row in a section/group, delete the section header too.
2. **`completed-tasks.md` is decision-only — not an archive.** For a routine
   fix, the commit/PR *is* the record; do **not** add a row. Add a row **only**
   when the task closes without a routine fix — a **not-a-bug determination, a
   superseded/ruled-out design, or a negative experimental result** whose
   reasoning a future agent would otherwise re-investigate. When you do, append
   one row to the matching section with the terminal decision + a pointer to the
   write-up (e.g. `dev/benchmarks/*.md`). Keep it to a line or two; the detail
   lives in the linked file, not the row.
3. Update `coordination.md` if strategic objectives are affected.
4. Run `bash dispatch.sh checkin <id> --done`.

### Parking (waiting for GHA / Hamilton / human review)

When the dispatched agent must stop and wait for an external event:

```bash
bash dispatch.sh checkin <id> \
  --kind=<gha|hamilton|human> \
  --ref=<run-id-or-ref> \
  --eta=<iso-datetime> \
  --resume="<one-sentence next action>"
```

Then exit cleanly. The dispatcher's `reap` subcommand surfaces parked agents
once their ETA has passed. `to-do.md` status flips to `PARKED (<id>, <kind> <ref>)`.

### User report intake (`a.*` / `u.*`)

The human files reports as individual files in the project root.
Each file contains a free-text description. The human's workflow is:
create file → write → save → never touch again.

**Naming convention:**
- `a.###` — app (Shiny) bug. Always routes to `### Shiny App` in `to-do.md`.
- `u.###` — general user issue (search quality, docs, API, etc.). Route by content.

**Agent responsibility:** Triage all pending files into `to-do.md` at
the start of every dispatched session (step 2 in Session start protocol above).

**Claim protocol:**
```bash
# List unclaimed reports
ls a.[0-9]* u.[0-9]* 2>/dev/null | grep -v 'claimed'

# Skip short files (< 20 chars) — don't rename, don't touch
wc -c < a.010  # check size first

# Claim atomically (rename)
mv a.010 a.010.claimed-d1

# Read, triage into to-do.md, then delete
cat a.010.claimed-d1
# ... create to-do.md entry ...
rm a.010.claimed-d1
```

**Skip guard:** Files shorter than 20 characters are likely mid-edit.
Do **not** rename them — leave in place for a later pass.
(Renaming and renaming back triggers RStudio "file moved" dialogs.)

**Shiny (`a.*`) fixes** are committed directly to `cpp-search` (bug
fixes in `inst/Parsimony/`, no feature branch needed). Use a temporary
worktree if changes span multiple files.

### Standing tasks

| ID | Type | Expertise file |
|----|------|---------------|
| S-RED | Red-team review | `dev/expertise/red-team.md` |
| S-PROF | Performance profiling | `dev/expertise/profiling.md` |
| S-COORD | Coordination review | `dev/expertise/coordination.md` |

Priority: P3 when ≥6 OPEN tasks, P2 when 3–5, P1 when <3.

### Key files

| File | Purpose |
|------|---------|
| `a.###` | App (Shiny) bug reports → triage to `### Shiny App`, then delete |
| `u.###` | General user issue reports → triage to matching section, then delete |
| `to-do.md` | Task queue (active/open tasks only) |
| `remote-jobs.md` | Pending async jobs (Hamilton SLURM, long GHA) — check at session start |
| `completed-tasks.md` | Decision-only log: not-a-bug / superseded / negative-result closures. `grep` before reopening a closed task; don't archive routine fixes here |
| `coordination.md` | Strategic plan |
| `AGENTS.md` | Conventions + workflow reference |
| `.dispatch/state.json` | Live dispatcher state (active agents, check-ins, budget tally) |
| `dev/expertise/*.md` | Standing-task methodology references |
| `dev/dispatch/ranker.txt` | Haiku ranker prompt template used by `dispatch.sh` |
| `dev/dispatch/agent-brief.md` | Spawned-agent system prompt template used by `dispatch.sh` |

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
