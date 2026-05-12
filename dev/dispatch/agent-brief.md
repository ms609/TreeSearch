# Dispatcher Agent Brief

You are agent **{{AGENT_ID}}**, assigned to task **{{TASK_ID}}**: `{{TASK_ROW}}`

## Budget & Model

- **Budget**: {{BUDGET_MINUTES}} minutes (stay within this slice; if work won't fit, do a sub-step and check in)
- **Model assigned**: {{MODEL}}
- **Effort level**: {{EFFORT}}
- **Resume action** (if parked): {{RESUME_HINT}}

## Workflow

1. **Startup intake** (before claiming work):
   - Triage any new user reports (`a.*` and `u.*` files in project root)
   - Check `remote-jobs.md` for pending async results
   - See AGENTS.md for full protocols

2. **Read conventions**: See `AGENTS.md` for:
   - Build/test/branch rules (GHA-first validation, tarball builds, `.agent-{{AGENT_ID}}/` isolation)
   - Shared-file coordination (append-only for `ts_rcpp.cpp`, `TreeSearch-init.c`)
   - Feature branch lifecycle and mandatory pre-commit checks
   - Multi-agent workflow (worktree reserved tasks, user-report claim protocol)

3. **Worktree rule**: If you need a worktree, create it under `../worktrees/TS-<name>`.
   **Never** switch the main `C:/Users/pjjg18/GitHub/TreeSearch` checkout to a
   different branch — it must stay on `cpp-search` (or the current feature branch).

4. **Build isolation**: Use `.agent-{{AGENT_ID}}/` as the install library
   ```bash
   SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
     rm -f src/*.o src/*.dll && \
     (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
     R CMD INSTALL --library=.agent-{{AGENT_ID}} "$TMPBUILD"/TreeSearch_*.tar.gz && \
     rm -rf "$TMPBUILD"
   ```

5. **Validation via GHA** (never run full test suites or R CMD check locally):
   - Push your branch: `git push -u origin feature/<name>`
   - Dispatch checks: `bash gha-dispatch.sh agent-check.yml feature/<name>`
   - Poll results: `bash gha-poll.sh <run_id>` (from another agent slice; don't block)

6. **Exit protocol**:

   **When blocking on external wait** (GHA, Hamilton, human review):
   ```bash
   bash dispatch.sh checkin {{AGENT_ID}} \
     --kind=<gha|hamilton|human|other> \
     --ref=<id> \
     --eta=<iso-8601-datetime> \
     --resume="<one-sentence next action>"
   ```
   Exit cleanly. The dispatcher will park this task and resume when the ETA passes.

   **When complete**:
   - Update `to-do.md` (delete task row; create new sections if needed)
   - Append summary row to `completed-tasks.md` under today's date
   - Call:
     ```bash
     bash dispatch.sh checkin {{AGENT_ID}} --done
     ```
   - The dispatcher will mark the agent slot as free.

## Budget discipline

If the work won't fit in {{BUDGET_MINUTES}} minutes:
1. Do a **meaningful sub-step** (fix one bug, implement one small feature, resolve one blocker)
2. Check in with a resume action: `bash dispatch.sh checkin {{AGENT_ID}} --kind=other --eta=<next> --resume="<next step>"`
3. Exit cleanly rather than blowing the budget

## Tools

- `.AGENTS/memory/` — technical references (architecture, testing, benchmarking, conventions)
- `todo-lock.sh` — lock protocol for coordinating `to-do.md` changes
- `gha-dispatch.sh` / `gha-poll.sh` — GitHub Actions integration
- Claude Code skills — use `skill(skill: "hamilton-hpc")` for Hamilton SLURM, `skill(skill: "r-package-profiling")` for profiling
