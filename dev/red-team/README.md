# `dev/red-team/` — the red-team rotation's state

Everything the `/red-team` skill reads and writes. The skill itself is user-level and shared
across projects (`~/.claude/skills/red-team/`), rewritten 2026-08-05 around GitHub issues —
its dual-mode `findings.md` machinery is retired to `~/.claude/skills-retired/red-team/`, the
only surviving record of file mode. **This directory is the project-local memory that makes
it work.** If a file
here goes stale, the next round pays for it in wasted finder tokens — that is not
hypothetical, it is what the 2026-07-27 rounds measured.

## The files

| File | What it is | Who writes it |
|------|------------|---------------|
| [`focus-areas.md`](focus-areas.md) | The rotation table: 15 numbered areas, the files each owns, its `start_tier`, and its key questions. Built once, edited rarely. **Adding a row** also needs, and nothing currently automates: a matching `area:N` GitHub label (`gh label create area:N --description "Red-team focus area N"`), and recomputing `N` in `log.md`'s rotation-formula comment (see RT12-01). | A round, when it finds the scope row wrong |
| [`log.md`](log.md) | Append-only, **newest first**. One entry per round (`area` / `reviewed_by` / `date` / `tier` / `yield` / `notes`), the **model-version legend** at the top, and `last_focus:` at the very bottom. | Every round |
| **GitHub issues** in [`agent-issues/TreeSearch`](https://github.com/agent-issues/TreeSearch/issues?q=label%3Ared-team) | **OPEN verified findings live here since 2026-08-04**, labelled `red-team` + `sev:high\|med\|low` + `area:N`. Status is GitHub state, so it cannot drift from merge state. Filed *after* verification; trivial issues are fixed inline and noted in `log.md` instead. | A round files; a merged `Fixes #N` closes |
| [`findings-archive.md`](findings-archive.md) | **FROZEN 2026-08-04.** Terminal-state findings from the file era, one compressed line each. **Offline anti-duplication memory, not a trophy case** — the one thing the tracker doesn't provide. | Nobody; it is closed to new rows |
| [`migration-map.tsv`](migration-map.tsv) | Every historical `T-nnn` → its issue number, archive entry, or open-PR reference. `T-nnn` ids are **frozen, not retired**: they persist in shipped source comments and in `log.md`. | Written once, at migration |
| [`migration-map-todo.tsv`](migration-map-todo.tsv) | Same idea, for the pre-tracker `T-nnn` ids that became `task` issues (#27+) rather than `red-team` findings — kept separate because it maps a different label family. | Written once, at migration |
| [`escalation-backlog.md`](escalation-backlog.md) | Seams that are re-eligible *now* but are not next in rotation — chiefly ones reopened by a model-version bump. Split into Open / Resolved-history. | `revisit`, and rounds that leave a residual |
| `proofs/` | Written derivations backing a specific finding (e.g. `union-construct-lower-bound.md`). | Whoever needs one |
| `heavy-tests/` | Standalone harnesses too slow for the test suite (e.g. `impose_validity/`, the T-327/T-333 constraint-repair driver). | Whoever needs one |
| `reviews/` | Per-review working notes and repro scripts (e.g. `cpp-search-sect-colreduce/`, which backs open T-335). | Whoever needs one |
| `union-of-finals-bound-proof.md`, `union-construct-gate.R` | **Do not move these into `proofs/`.** They live at top level because things outside this directory cite them *by path* — most importantly a comment in shipped source, `src/ts_fitch.cpp:421`, plus `dev/plans/2026-07-14-mission-B-kernel-8x-goloboff-screen.md`, `dev/profiling/mission-b-goloboff-gates.md` and `proofs/union-construct-lower-bound.md` itself. Relocating them for neatness breaks a source-code reference. | math-prover lane; leave in place |

Working artifacts under `proofs/`, `heavy-tests/` and `reviews/` are **live** as long as the
finding they back is open. Never sweep them without checking the **open issue list** first.

**Globbing gotcha (found 2026-08-04 doing the area-12 file-coverage diff):** `R/` holds one
lowercase-extension file, `R/pp_info_extra_step.r` — a case-sensitive `ls R/*.R` / `Glob`
pattern silently skips it. Any future scope-coverage diff should glob `R/*.[Rr]`, not `R/*.R`.

**Current state (2026-08-04, post-migration):** **24 open issues** and **61 archived rows** in
`findings-archive.md`. Two findings (former T-395, T-396) were not migrated because they have
an open upstream PR.

Do not maintain the open count, or its severity breakdown, by hand — both are now a query (a
hand-kept breakdown drifted within the same round it was written: an area-12 audit on
2026-08-04 found this line reading 6/5/13 against a true 6/4/14):

```bash
gh issue list --repo agent-issues/TreeSearch --label red-team --state open --json number --jq length
gh issue list --repo agent-issues/TreeSearch --label red-team --state open --json labels \
  --jq '[.[].labels[].name | select(startswith("sev:"))] | group_by(.) | map({(.[0]): length}) | add'
```

The 2026-07-27 history is still worth knowing before reading an empty high-severity column as a
finished job: a `tidy` pass took the table to 11 open with neither a P1 nor a P2 — the first
time that had happened — and one opus round on a seam whose files had been touched only
*incidentally* by four other-area commits put both straight back.

## The finding lifecycle

Since 2026-08-04 the lifecycle is **GitHub's**, not a status column:

```
open  →  open + `in-progress` (claimed; the claiming comment names the branch)
      →  closed by `Fixes #N` on merge into the fork's default branch
      →  closed as "not planned" (+ `wontfix`) when adjudicated not worth fixing
```

**Do not reintroduce a status column.** It existed only because a file cannot observe a merge.
`Fixes #N` can, and that is what removed the drift described below — so re-adding it re-adds
the bug.

A **closed issue is the archive**: it stays searchable forever, which is the whole point of
archive-never-delete.

### Three rules this directory has actually been burned by

**1. `Fixes #N` only fires on merge into the *default* branch.** The fork's default branch is
`cpp-search` for exactly this reason. Merging a fix into any other branch closes nothing, and
the issue silently stays open.

**2. "Landed" means present in `cpp-search` HEAD — not merged to `main`.** Development happens
on `cpp-search`, which is **1108 commits ahead of `main`** (and `main` carries one commit
`cpp-search` lacks). Generic lifecycle text saying "merged to main" would make nothing ever
closable here. Cite **commit SHAs**.

**3. Archive, never delete — and never trust a hand-maintained status.** `findings.md` carried a
"when a finding lands, remove the row" instruction until 2026-07-27; two rows (T-331, T-333)
were deleted under it and had to be reconstructed from git log. Worse, at the 2026-08-04
migration **23 of 49 rows still read as open when their fixes had already landed**, only four
carrying the "awaiting `tidy` archive" marker. Both failures share one cause: keeping status in
a file that no merge can update. That is why findings moved to issues.

## Running it

```bash
claude -p "/red-team"
```

- `/red-team` — next area in rotation, at its earned tier · `/red-team <N> [tier]` — force one
- `/red-team status` — what's been done, what's next
- `/red-team revisit` — a rung's model version moved; re-mine seams it makes re-eligible
- `/red-team tidy` — this housekeeping pass. No finder, no verifier, no new findings, and it
  does **not** advance `last_focus:`
- `/red-team deep` — Workflow fan-out. Costs far more; explicit opt-in only

Durable lessons that outlive any one round live in [`../expertise/red-team.md`](../expertise/red-team.md).
