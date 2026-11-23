# `/next-issue` config — TreeSearch

Read by `~/.claude/skills/next-issue/SKILL.md`. Only what varies per repo lives
here; the doctrine is in the skill.

| Key | Value |
|-----|-------|
| `base_branch` | **`cpp-search`**, not `main` — the integration branch, on `origin` = `agent-issues/TreeSearch`. **Never `upstream` (`ms609/TreeSearch`)**: push URL is `no-push-use-gha`, and `ms609-agent` is read-only there by design. Work reaches `ms609` only when the human syncs the fork. |
| `grouping` | `labels`, prefix `area:` (1–15). Fall back to `hot_files` for an unlabelled issue. |
| `exclude_label` | `deferred` is **not** the maintainer-call filter — it means "not now", not "needs a human decision". Exclude by judgment per the skill until `needs-decision` exists here. |
| `identity` | `agent` — `ms609-agent`, token env var `CLAUDE_GH_TOKEN`. Table in `~/.claude/CLAUDE.md`. |

## `hot_files` — never split across parallel chips

Derive the current list from the tranche's `area:N` label plus
`dev/red-team/focus-areas.md`, which maps each area to its files. The standing
coordination-critical ones:

```
src/Makevars.win          # shared; coordination rules in AGENTS.md
src/TreeSearch-win.def    # shared; append-only
src/TreeSearch-init.c     # arg-count mismatches are a known build failure mode
```

## `pr_command`

```bash
GH_TOKEN=$CLAUDE_GH_TOKEN gh pr create --base cpp-search --head <branch> --reviewer ms609 --body-file <file>
```

`--base cpp-search` is not optional — the fork's default branch governs whether
`Fixes #N` closes on merge. Reads need no token prefix. Never `export GH_TOKEN`.

## `build` / `test`

GHA is the primary validation path. Local builds are for targeted iteration:

```bash
ID="a$(printf %s "$BRANCH" | sha1sum | cut -c1-5)"
bash build-agent.sh TreeSearch-a "$ID"
bash test-agent.sh  TreeSearch-a "$ID" [filter]
```

Derive the id from the branch, never the session id — a resumed chip gets a new
session id and would silently build a second tree.

Targeted filter always — never a full suite.

## `branch_rule`

One agent owns a `feature/*` branch at a time; **never `git checkout` a branch
you do not own** — use `git worktree add`. **Never `git stash`**: the stash is
repo-wide across worktrees and can apply another agent's work into yours; use a
patch file. Coordination files live on `cpp-search` only.

## `extra`

- **Max 2 cores per agent.**
- Debug `.o` contamination and DLL locks are the two recurring build failures —
  recovery steps in `AGENTS.md` → *Build failure recovery*.
- Prefer `TreeTools` over `ape` (see `AGENTS.md`).
- `dev/red-team/log.md`, `findings-archive.md` and `migration-map*.tsv` are
  **frozen**. Never add a row to any of them, and never to
  `escalation-backlog.md` — findings are issues.
