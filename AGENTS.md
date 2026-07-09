# TreeSearch Multi-Agent Development Notes

Always check the contents of `.AGENTS` for memories and policies relevant to
the task you have been assigned.

Update memory files with anything relevant you learn. But keep them lean.

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

`to-do.md`, `completed-tasks.md`, `coordination.md`, and `AGENTS.md` are
**never committed on feature branches**. When a dispatched
agent working on a feature branch needs to claim a task or update coordination
files, they commit those changes directly to `cpp-search` (coordination-only
commit), keeping the feature branch clean.

To read coordination files while on a feature branch without switching:
```bash
git show cpp-search:to-do.md
git show cpp-search:coordination.md
```

### Shared files at merge time

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

### Key files

| File | Purpose |
|------|---------|
| `to-do.md` | Task queue (active/open tasks only) |
| `coordination.md` | Strategic plan |
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
