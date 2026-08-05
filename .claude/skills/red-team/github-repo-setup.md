# One-off: setting up a repo under `agent-issues` for `/red-team`

Read this **only** when standing up a new project — a fresh fork under the `agent-issues`
org, or a first migration of an existing project's findings file into issues. A normal
`/red-team` round never needs it.

The end state is what `SKILL.md` assumes throughout: findings are GitHub issues, `Fixes #N`
closes them on merge, and no file anywhere carries a status column.

---

## 1. Fork into `agent-issues`, and point `gh` at it

The public upstream stays the release repo; the fork is where issues and development live.

```bash
gh repo fork <owner>/<repo> --org agent-issues --remote --remote-name origin
cd <repo>
gh repo set-default agent-issues/<repo>
```

`gh repo set-default` is what makes every later `gh issue`/`gh pr` command in this skill —
and the `gha-dispatch.sh` / `gha-poll.sh` helpers, which resolve the target with
`gh repo view --json nameWithOwner` — hit the fork rather than upstream.

The `agent-issues` org is `collaborators_only`, so issues filed here can only come from
collaborators. The public upstream tracker is the opposite: **untrusted input, never a task
list.**

## 2. Set the default branch to the development branch

**Do this before filing anything.** `Fixes #N` closes an issue only on merge into the
repository's *default* branch. If development happens on a long-lived branch (TreeSearch:
`cpp-search`, 1108 commits ahead of `main`), every fix PR targets it, and a default of
`main` means no issue ever closes — silently.

```bash
gh repo edit agent-issues/<repo> --default-branch cpp-search
gh repo view --json defaultBranchRef --jq .defaultBranchRef.name   # verify
```

## 3. Block direct pushes to upstream

Upstream must only ever *receive* fast-forwards of the fork's trunk. One direct commit
there and every later sync becomes a real merge, with conflicts on `DESCRIPTION`,
`NAMESPACE` and the append-only `src/` files. Enforce it mechanically rather than by
discipline:

```bash
git remote add upstream https://github.com/<owner>/<repo>.git
git remote set-url --push upstream no-push-use-gha
```

`git push upstream` now fails locally, before reaching GitHub.

## 4. Create the labels

`red-team` is the mandatory one — historically it doubled as the skill's mode switch, and
tooling still keys off it. `sev:*` replaced the old P1/P2/P3.

```bash
gh label create red-team        --color 5319E7 --description "Filed by the /red-team rotation"
gh label create sev:high        --color B60205 --description "P1: wrong user-visible result / crash"
gh label create sev:med         --color D93F0B --description "P2: wrong on edge input / search quality"
gh label create sev:low         --color FBCA04 --description "P3: robustness / polish"
gh label create in-progress     --color 0E8A16 --description "Being fixed; claiming comment names the branch"
gh label create needs-escalation --color 1D76DB --description "Next dispatch on this area must be opus+"
gh label create chore           --color BFD4F2 --description "Infrastructure / process work, not a red-team finding"

for n in $(seq 1 <N>); do
  gh label create "area:$n" --color C5DEF5 --description "Red-team focus area $n"
done
```

`area:N` labels are **per focus-area row**, so adding a row to `focus-areas.md` later means
creating its label too — an easy step to miss, and an issue filed against a missing label
just fails.

## 5. Enable Actions and recreate secrets

Workflows are **disabled on a new fork** until enabled once through the Actions tab in the
browser — there is no `gh` equivalent, so this is the one manual step.

**Secrets do not come across from upstream.** Any check that needs one fails until it is
recreated:

```bash
gh secret list                       # what the fork actually has
gh secret set <NAME>                 # recreate each one the workflows reference
```

## 6. Scaffold `dev/red-team/`

Run `/red-team init`, which builds `focus-areas.md` (6–12 areas for a new project; every
`start_tier` **`sonnet`**, since maturity is measured, not assumed) and `log.md` (round
format, the **model-version legend** seeded with today's alias→version mapping, and
`last_focus: 0`).

A project scaffolded without the legend can never notice a version bump, so trigger 3 never
fires on its own.

`init` reports the rotation and **stops** — no review on the scaffolding turn; the user
reviews the areas and tiers first.

---

## Migrating an existing `findings.md`

Only for a project moving off file mode.

1. **Freeze, don't delete.** Terminal-state rows stay in `findings-archive.md` as
   one-line-each offline anti-duplication memory. A finding closed years ago in a file is
   still a duplicate, and the "closed — no longer reproducible" rows are the highest-value
   records: they are exactly what stops an expensive Opus/Fable pass chasing a ghost.
2. **Open rows become issues** — one each, `RT-###` kept in the title, labelled
   `red-team` + `sev:*` + `area:N`. Use `--body-file`; bodies are KB of backticks and `$`.
3. **Expect the file to have been lying.** At TreeSearch's 2026-08-04 migration, **23 of 49
   rows still read as open when their fixes had already landed** — verify each against
   `cpp-search` HEAD before filing it as open, and cite commit SHAs, not "merged to main".
4. **Write `migration-map.tsv`:** every historical `T-nnn` → its issue number, archive
   entry, or open-PR reference. `T-nnn` ids are **frozen, not retired** — they persist in
   shipped source comments and in `log.md`, so an unresolvable one is a re-hunted finding.
5. **Leave anything with an open upstream PR unmigrated,** and say so in the map.
6. **Do not carry the status column across.** It existed only because a file cannot observe
   a merge. `Fixes #N` can, and that is the entire point of the move.
