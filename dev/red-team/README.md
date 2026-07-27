# `dev/red-team/` — the red-team rotation's state

Everything the `/red-team` skill reads and writes. The skill itself is user-level and shared
across projects; **this directory is the project-local memory that makes it work.** If a file
here goes stale, the next round pays for it in wasted finder tokens — that is not
hypothetical, it is what the 2026-07-27 rounds measured.

## The files

| File | What it is | Who writes it |
|------|------------|---------------|
| [`focus-areas.md`](focus-areas.md) | The rotation table: 13 numbered areas, the files each owns, its `start_tier`, and its key questions. Built once, edited rarely. | A round, when it finds the scope row wrong |
| [`log.md`](log.md) | Append-only, **newest first**. One entry per round (`area` / `reviewed_by` / `date` / `tier` / `yield` / `notes`), the **model-version legend** at the top, and `last_focus:` at the very bottom. | Every round |
| [`findings.md`](findings.md) | **OPEN** verified findings only. Filed *after* verification; trivial issues are fixed inline and noted in `log.md` instead. | A round adds; `tidy` never files |
| [`findings-archive.md`](findings-archive.md) | Terminal-state findings, one compressed line each. **Anti-duplication memory, not a trophy case.** | `tidy` only |
| [`escalation-backlog.md`](escalation-backlog.md) | Seams that are re-eligible *now* but are not next in rotation — chiefly ones reopened by a model-version bump. Split into Open / Resolved-history. | `revisit`, and rounds that leave a residual |
| `proofs/` | Written derivations backing a specific finding (e.g. `union-construct-lower-bound.md`). | Whoever needs one |
| `heavy-tests/` | Standalone harnesses too slow for the test suite (e.g. `impose_validity/`, the T-327/T-333 constraint-repair driver). | Whoever needs one |
| `reviews/` | Per-review working notes and repro scripts (e.g. `cpp-search-sect-colreduce/`, which backs open T-335). | Whoever needs one |
| `union-of-finals-bound-proof.md`, `union-construct-gate.R` | **Do not move these into `proofs/`.** They live at top level because things outside this directory cite them *by path* — most importantly a comment in shipped source, `src/ts_fitch.cpp:421`, plus `dev/plans/2026-07-14-mission-B-kernel-8x-goloboff-screen.md`, `dev/profiling/mission-b-goloboff-gates.md` and `proofs/union-construct-lower-bound.md` itself. Relocating them for neatness breaks a source-code reference. | math-prover lane; leave in place |

Working artifacts under `proofs/`, `heavy-tests/` and `reviews/` are **live** as long as the
finding they back is open. Never sweep them without checking `findings.md` first.

**Current state:** **11 open findings** — 0 P1, 0 P2, 11 P3 — and **36 archived**, as of the
`tidy` pass of 2026-07-27. That is the first time the open table has held neither a P1 nor a
P2. (Count is maintained by `tidy`; if it disagrees with `findings.md`, trust the file.)

## The finding lifecycle

```
open  →  open (fix PR #N) | in-review (PR #N)  →  fixed (<sha>) | closed (Round R) | wontfix  →  archived
```

- **`open`** — verified, unfixed.
- **`open (fix PR #N)` / `in-review (PR #N)`** — a fix is proposed, not yet landed.
- **`fixed (<sha>)`** — the fix is in the tree. **`closed (Round R)`** — no longer
  reproducible (superseded by an unrelated change, or inert on inspection); kept as an
  anti-dup signal. **`wontfix`** — adjudicated not worth fixing.
- **`archived`** — relocated to `findings-archive.md`, one line, provenance preserved.

### Two rules that this directory has actually been burned by

**1. "Landed" means present in `cpp-search` HEAD — not merged to `main`.** Development happens
on `cpp-search`, which is ~180 commits ahead of `main`. The skill's generic lifecycle text
says `fixed (PR #N)` = "merged to main"; taken literally here, nothing would ever be
archivable. Cite **commit SHAs**.

**2. Archive, never delete.** A resolved finding *moved* to `findings-archive.md` still stops a
future finder re-hunting it. A resolved finding *deleted* does not — and the deletion leaves no
trace. `findings.md` carried a "when a finding lands, remove the row" instruction until
2026-07-27; two rows (T-331, T-333) were deleted under it and had to be reconstructed from
git log, and five more went stale because flipping a status was nobody's job. Only the `tidy`
pass moves rows; a rotation round only ever adds them.

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
