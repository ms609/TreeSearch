---
description: Run a focused red-team review of TreeSearch — one focus area per invocation, at a model tier earned by recorded yield, severity-matched verification, findings filed as GitHub issues in agent-issues/TreeSearch.
when_to_use: When the user says /red-team, asks for a red-team review, code audit, bug hunt, or wants the next round of adversarial review. Each invocation reviews ONE focus area and rotates to the next next time.
---

# /red-team skill

You are the **orchestrator** of the red-team rotation. The bug-hunting and verification
doctrine lives in two agents — `red-team-finder` and `red-team-verifier` — so this file is
pure orchestration: which area, at what tier, verified how, recorded where. You never
re-read the finder/verifier doctrine; you hand each a short dynamic brief and parse what
comes back.

Findings are **GitHub issues** in `agent-issues/TreeSearch` (`gh` already defaults to it),
labelled `red-team` + `sev:*` + `area:N`. Rotation state is `dev/red-team/`. There is no
`findings.md`: it was deleted at the 2026-08-04 migration because a file cannot observe a
merge, and 23 of its 49 rows read as open when their fixes had already landed. The
file-mode doctrine survives only in the retired user-level copy at
`~/.claude/skills-retired/red-team/SKILL.md`; do not reintroduce it here.

Setting up a *new* repo under `agent-issues` so this skill can run against it is a one-off
that most sessions never touch: [`github-repo-setup.md`](github-repo-setup.md).

The goal is **issues fixed per token spent**, not issues found in the abstract. Depth over
breadth: one focused review that finds a real bug beats a broad sweep that confirms "all
green." Spend cheap capability where it yields; reserve premium capability for where a
missed bug is catastrophic.

---

## Model tiers — the core idea

Bug-finding is a **capability cliff**, not a smooth gradient. A cheaper model doesn't just
find *fewer* bugs — it is *blind* to whole classes of subtle ones. Two consequences:

- **Maturity is measured, not assumed.** A cheap first pass *measures* an area: if Sonnet
  keeps finding glaring bugs, it's immature — keep mining cheaply. When Sonnet runs dry,
  that does **not** mean the area is clean; it means Sonnet is tapped out.
- **An empty cheap sweep is a signal to escalate, never that the area is clean.** "Sonnet
  found nothing" on a heavily-revised numerical core is meaningless — the remaining bugs
  are below Sonnet's cliff by construction.

The finder ladder is **`sonnet` → `opus` → `fable`** (Agent-tool `model` values). Haiku is
*not* a finder — too far below the cliff; it is the cheap **verifier**.

| Tier | `model` | When it's the right finder |
|------|---------|----------------------------|
| 1 | `sonnet` | Default first pass for every area. A bug fixed at Sonnet prices is the best issues-per-token outcome. |
| 2 | `opus` | After Sonnet's seam runs dry, or for an area `focus-areas.md` has *proven* mature. |
| 3 | `fable` | Top escalation — most capable, ~2× Opus cost. Reserve for the highest-magnitude core (Fitch/NA scoring, topology invariants, constraint machinery) where a missed bug is catastrophic and Opus's seam has run dry. |

**Three escalation triggers** (any fires):

1. **Dry seam** — the last visit at the current tier found nothing new. Escalate one step —
   **but check *Model versions* first: a version bump is a cheaper step than a rung bump,
   and takes precedence.**
2. **High-severity signal** — a finder flags something high-severity it couldn't pin down.
   Escalate *immediately* to confirm and sweep adjacent code; don't wait for the seam to
   run dry.
3. **Rung version bump** — the model behind the tier alias has been superseded since the
   dry verdict was recorded. Re-visit at the **same rung on the new version** before
   escalating a rung or honouring a dormant record.

**While a seam keeps yielding, re-visit at the same tier with a fresh agent.** A second
agent with no memory of the first takes a different angle and catches overlooked things
even at the same price — fresh-angle recall is the cheapest recall there is.

**`start_tier`** lives in `focus-areas.md`, one per area, and encodes *measured* maturity
from the round history — not a guess. A newly-added area starts at `sonnet` and is
annotated **UNMEASURED / no inherited maturity**. So does a **row that grows**: a scope row
that gains files does not inherit the dry verdicts earned before it grew, because those
verdicts were about the old file list.

---

## Model versions — a tier is a rung, not a model

`sonnet` / `opus` / `fable` / `haiku` are Agent-tool **aliases**. Each names a rung; the
*version* it resolves to changes under you. You cannot request a version — the enum has no
version field — so the skill **records what ran**, it never pins.

This matters because the tier doctrine rests on "this model ran dry here." That verdict is
evidence about **one model version**, not about the rung forever. Left unversioned, the
tracker silently converts "Opus 4.8 found nothing" into "opus is tapped out" — retiring an
area, or paying for a `fable` escalation, on stale evidence.

- **Never pin; always take the newest.** Pass `model: "opus"` (not a dated id), and keep
  the `red-team-finder` / `red-team-verifier` profile frontmatter on bare aliases. A dated
  id anywhere in that chain freezes a rung at a superseded version.
- **Reconcile the legend at round start** (step 3 of *Normal run*) against the
  model-version legend at the top of `dev/red-team/log.md`. Moved? Add a legend row and
  fire trigger 3 **before** dispatching — otherwise a bump is only ever noticed by a human.
- **Stamp every round with the version that ran:**
  `## Round R — Area N (<name>) — <tier> (<Model Version>) — <date>`. If you cannot confirm
  what the alias resolved to, record `opus (version unconfirmed)` rather than asserting.
- **Backward-looking verdicts are version-scoped.** Write `dormant at <rung>-<version>`.
- **Forward-looking routing stays unversioned.** "Escalate to opus" names the rung.
- **Version bump before rung bump:** `opus-4.8 dry → opus-5 → fable`.
- **A version bump reopens dormancy** — automatically. Reopen it (label the relevant open
  issue `needs-escalation`, or add a row to `escalation-backlog.md` if no issue exists yet)
  rather than waiting for rotation to reach it.
- **Two independent dry verdicts is the real dormancy bar.**
- **A revisit must change its angle, not just its model** — brief the finder to attack the
  *prior round's derivation*, pasted in as claims to break.

The alias→version table is **project-local**, in `dev/red-team/log.md`. Never hardcode one
here.

---

## Severity-matched verification

The finder reports everything to maximise recall; the verifier restores precision. Routing
is **by severity**:

- **Low-severity OR low-confidence** → one **`haiku`** `red-team-verifier` batch (cheap
  refute-or-confirm). This filter is what licenses report-everything.
- **High-severity (any confidence)** → a `red-team-verifier` at **peer-or-higher tier**
  (`model` ≥ finder's tier). A subtle scoring or invariant claim is as hard to *verify* as
  to *find*; a cheap false-refute would drop the crown-jewel bug. Never cheap-verify a
  high-severity claim.

Verification is **sequential** (finder first, then verifier passes) — the verifier is a
*filter*, not a parallel finder.

---

## Arguments

`/red-team` — run the next area in rotation, at its earned tier.
`/red-team <N>` — force area #N.
`/red-team <N> <tier>` — force area #N at an explicit tier, overriding the earned tier for one round.
`/red-team init` — (re)scaffold `dev/red-team/` and rebuild the focus-area list.
`/red-team status` — read `log.md`; summarise rounds done, each area's tier and yield, next up.
`/red-team revisit` — **model-version sweep, no rotation advance.** Update the legend, reopen every `dry`/`dormant`/`retired` verdict scoped to a superseded version, re-visit at the same rung with a fresh-angle brief. Does not touch `last_focus:`.
`/red-team tidy` — **housekeeping, no finder.** Scope-coverage diff, link integrity, artifact sweep, map refresh. See *Tidy pass*.
`/red-team deep` — **explicit opt-in.** Run the current area through a Workflow that fans out dimension-finders and atomises verification. Costs far more.

---

## Where findings live

| Label | Meaning |
|-------|---------|
| `red-team` | Every finding this skill files. Mandatory. |
| `sev:high` / `sev:med` / `sev:low` | Severity. |
| `area:N` | **Which area owns the code**, not which round found it. An issue may carry several — 10 of 24 open issues did at migration — and `--label needs-escalation,area:N` ANDs correctly, so a multi-labelled issue is visible to each of its areas. |
| `in-progress` | Someone is fixing it; the claiming comment names the branch. |
| `needs-escalation` | **A tier flag and nothing else:** this area's next dispatch must be `opus`+. See *Escalation tracking* for what must NOT be encoded here. |

File with `--body-file`, never `--body "$string"` — bodies run to several KB of backticks,
quotes and `$`, and shell quoting will mangle them:

```bash
gh issue create --title "RT-437: <one-line claim>" \
  --label red-team,sev:high,area:3 --body-file <tmpfile>
```

Keep the `RT-###` prefix. IDs are cited from `log.md`, PR bodies and prior findings; a
finding that cannot be found by its ID gets re-hunted. Allocate the next ID from the **max
across** open issues *and* the frozen `findings-archive.md`.

**Anti-duplication spans both stores.** Before filing, search
`gh issue list --search "RT- <keyword>" --state all` *and* grep `findings-archive.md`. A
finding closed years ago in a file is still a duplicate.

### Files (all under the project root)

```
dev/red-team/
  README.md             # Dir map + the finding lifecycle
  focus-areas.md        # Rotation table — the files each area owns, its start_tier, key questions
  log.md                # Per-round notes, newest first; model-version legend at top; `last_focus:` at the bottom
  escalation-backlog.md # Residuals that are re-eligible now but not next in rotation
  findings-archive.md   # FROZEN 2026-08-04 — file-era terminal findings, one line each. Offline anti-dup memory
  migration-map*.tsv    # Historical T-nnn -> issue number / archive entry
  proofs/ heavy-tests/ reviews/   # Working artifacts backing specific findings — live while their finding is open
```

Two lifecycle rules this directory has been burned by:

- **`Fixes #N` fires only on merge into the default branch**, which is why the fork's
  default is `cpp-search`. Merging a fix anywhere else closes nothing, silently.
- **"Landed" means present in `cpp-search` HEAD, not merged to `main`.** Cite commit SHAs.

---

## Escalation tracking — the tier channel, and what must not go through it

`needs-escalation` encodes exactly one ask: *dispatch this area at `opus`+*. Keep it that
way.

- **A specific issue drives a tier escalation** (the common case): label that issue
  `needs-escalation` alongside its `area:N`. When the issue closes, it drops out of the
  query automatically — no second copy of "is this still open" to drift. This replaced
  hand-maintained rows that had gone stale for ~20 rounds.
- **A residual whose ask is *not* more capability stays prose** in
  `escalation-backlog.md` — cross-area routing ("someone owning area 9 should look"),
  sequencing ("harness first, then #18/#19"), or a soft signal not yet filed. **Do not
  promote these to `needs-escalation`.** It changes no decision when the area already
  starts at `opus`, and worse, a label hit is a short-circuit: it would suppress reading
  the backlog row that holds the actual ask. Labelling a sequencing question converts it
  into a tier answer and then hides the question. (This happened on 2026-08-04 to #18/#19
  and was reverted the same round.)
- **A *work-shape* verdict — "the next visit should not be a finder" — is not an escalation
  at all.** It gates dispatch; see step 5 of *Normal run*.
- `escalation-backlog.md` keeps the two things a label can't hold: unfiled soft signals,
  and the resolved-history narrative. It is not a duplicate ledger of issue state.

---

## Normal run

1. Read `focus-areas.md` and the **bottom** of `log.md` for `last_focus:`.
2. Next area: `(last_focus mod N) + 1`, where `N` is the **current row count** in
   `focus-areas.md` — recompute it, never trust a number written into prose. A stale `N`
   made areas 11–13 mathematically unreachable for a month.
3. **Determine the tier** (override if the user passed one):
   - **Reconcile the model-version legend** against what each alias resolves to *now*.
     Moved? Add a legend row, reopen the verdicts the bump makes re-eligible, let that
     drive routing below.
   - **Check for an open tier escalation:** `gh issue list --label needs-escalation,area:N
     --state open`. **Then read `escalation-backlog.md` for area N regardless of the
     result** — a label hit is not licence to skip it; the non-tier residuals that gate
     *how* to spend the round live only there.
   - Read the area's `start_tier` and its most recent stamped `tier: <rung> (<version>)` +
     `yield:`.
   - Never visited → `start_tier`.
   - Last visit **yielded** → same tier, fresh agent.
   - Last visit **empty** and the rung's **version has moved on** → same tier, newer model,
     fresh-angle brief. *This precedes rung escalation.*
   - Last visit **empty** at the current version → escalate one rung.
   - Already `fable` and empty → note `dormant at fable-<version>` and rotate on.
   - A prior round raised a **high-severity signal** → escalate immediately.
4. **Assemble the brief.** Read recent `log.md` entries for the area, and run
   `gh issue list --label red-team,area:N --state open --json number,title` — the area's
   open issues go into the brief. Without them a fresh finder can spend its whole budget
   re-investigating a mechanism an open issue already describes; the anti-duplication
   search at *filing* time prevents a duplicate filing, not duplicate investigation. (A
   verifier independently rediscovered T-400/#16 this way.)
5. **Work-shape gate — check before dispatching, and honour it.** If the most recent
   `log.md` entry for this area records an explicit verdict that the next visit should
   **not** be a finder — "NEXT VISIT: NOT another finder — a bounded exhaustive harness"
   (area 13, 2026-07-03), "the next visit should NOT be a finder", a wall-matched A/B owed
   (area 10, 2026-08-03) — then **stop and report that to the user instead of
   auto-dispatching.** Say what is owed, and offer either to do that work or to skip the
   area and rotate on. Do not launch a finder against a standing verdict not to; a tier
   decision cannot answer a work-shape question.
6. **Launch one finder** via the Agent tool: `subagent_type: "red-team-finder"`, `model`
   from step 3, `description: "Red-team area N (<tier>): <name>"`, `prompt`: the finder
   brief below.
7. **Verify** before filing — the finder returns *candidates*:
   - Low-sev / low-confidence → one `red-team-verifier` at `model: "haiku"`.
   - High-sev → a fresh `red-team-verifier` at `model:` ≥ the finder's tier.
   - Drop REFUTED; keep REAL with the verdict.
8. Append the round to `log.md`: area, `tier: <rung> (<version>)`, `yield:` (count of
   *confirmed* findings), verifier verdicts, escalation decision, and any work-shape verdict
   for the next visit — **state that unmissably**, since step 5 is what reads it back.
   File still-open confirmed findings as issues.
9. Report to the user (see *Reporting back*).

**Do not file an issue for a finding whose fix already landed this same round.** A
`gh issue create` immediately followed by `gh issue close` is notification noise and a
wasted number, with no advantage over the round's `log.md` entry — which is already the
grep-able anti-dup memory. Allocate the `RT-###` id (so it stays globally unique for later
cross-reference), describe it in `log.md`, and stop. `gh issue create` is for findings that
need someone else's attention.

---

## Briefs

Standing doctrine lives in the agent profiles; you pass only the per-round specifics.

**Finder brief** (to `red-team-finder`):

```
Focus area #N (<name>) on the project at <absolute path>. You are running at the <tier> tier.

## Files in scope
<list from focus-areas.md row>

## Key questions
<from focus-areas.md row>

## Prior rounds on this area (DO NOT REPEAT)
<most recent 1–3 log entries for area N — tier used, what was checked/found/ruled out>

## Open issues already filed against this area (DO NOT RE-INVESTIGATE)
<gh issue list --label red-team,area:N --state open — number + title each>
Anything here is known. Finding it again is not a finding. Adjacent mechanisms these
issues do NOT cover are fair game — say so explicitly if you go there.
```

On a **version-bump revisit** (trigger 3), append:

```
## This is a version-bump revisit
A prior visit at this same tier ran DRY and recorded these conclusions, reached by an
earlier, now-superseded model at this tier. Treat each as a CLAIM TO BREAK:

<paste the prior round's "clean by derivation" conclusions verbatim>

Re-deriving one independently and reaching the same answer is a useful result — say so.
But do not re-read the same files the same way: prioritise arguments resting on an
unstated assumption, an informal proof sketch, or "structurally unreachable" reasoning.
```

**Verifier brief** (to `red-team-verifier`):

```
Verify these red-team findings for the project at <absolute path>.

## Findings to verify
<paste the candidate rows: id | severity | confidence | title | file:line + detail>
```

**Return contracts** — Finder → a *Round summary* (incl. `Seam status: still yielding | ran
dry`), candidate rows `suggested-id | severity | confidence | [Bug/Perf] title | file:line
— detail`, a *High-severity signals* block, *Notes for next reviewer*. Verifier → rows
`id | REAL | REFUTED | one-line verdict + reproduction`.

---

## Deep mode (`/red-team deep` — explicit opt-in)

Default `/red-team` is sequential and user-paced. Deep mode is the exception: a **Workflow**
that fans out one `red-team-finder` per concern (numerical conditioning, cache coherence,
API/contract, edge cases, concurrency), each reading the full scope in its own context;
atomises verification (one `red-team-verifier` per finding, `pipeline`, severity-matched);
and synthesises confirmed findings into issues. This is a *thoroughness* play, not a
cheapness one — N× the reading cost.

Requires the Workflow tool's explicit opt-in. Script lives at
`dev/red-team/deep.workflow.js` — author it on first `deep` use, then reuse it.

Do NOT fan out parallel **finders** in default mode: the user drives the cadence to react
between rounds.

---

## Tidy pass (`tidy` — housekeeping, no finder)

Reconciles what is filed against reality and restructures for findability. It **never files
or fixes a bug** — that is the rotation's job — and it does **not** advance `last_focus:` or
touch tier/yield state.

1. **Scope-coverage diff — the highest-value step.** Diff every file in the tree against
   the union of `focus-areas.md`'s scope rows. Files owned by no area are never reviewed at
   any tier, at any point in the rotation, and this has cost real findings twice: three
   findings filed in two unowned `R/` files a month after the gap was flagged, and a
   `sev:high` `TreeLength()` out-of-bounds write in an unowned family. **Glob `R/*.[Rr]`,
   not `R/*.R`** — `R/pp_info_extra_step.r` is lowercase and is silently skipped by a
   case-sensitive pattern. Add unowned files to the right row, annotated **UNMEASURED / no
   inherited maturity**; a subsystem too large or too distinct to fold in is a new area,
   which means `gh label create area:<N>` and a recomputed `N`.
2. **Reconcile issues against merge state.** List open `red-team` issues whose linked PR
   merged without closing them (a PR body missing `Fixes #N`) and close those by hand.
   Clear `in-progress` from any issue whose claiming branch no longer exists. Sweep
   `escalation-backlog.md` for rows that now name a specific issue *and whose ask is a tier
   escalation* — promote those to a `needs-escalation` label and drop the row. Leave
   sequencing and cross-area-routing rows as prose (see *Escalation tracking*).
3. **Link integrity.** Confirm every `proofs/`, `heavy-tests/`, `reviews/` and `[[memory]]`
   cross-ref still resolves. Artifacts backing an **open** issue are live — never sweep
   them. Some paths are cited from shipped source (`src/ts_fitch.cpp:421` cites
   `union-of-finals-bound-proof.md`); relocating them for neatness breaks the reference.
4. **Artifact hygiene.** `.gitignore` covers run outputs (`heavy-tests/*.log`, experiment
   `*.rds`, `*-results/`) without catching committed fixtures (`RT-area*-*.rds`). Remove
   stray dumps.
5. **Keep the map current.** Refresh `README.md`. Never hand-maintain a count — a
   hand-kept severity breakdown drifted within the same round it was written. Use the
   query:
   ```bash
   gh issue list --label red-team --state open --json number --jq length
   ```
   Over ~150 open findings is a genuine backlog signal — surface it, don't hide it.

Run every ~10 rounds, when the open count exceeds ~200, or on demand.

---

## Reporting back to the user

```
Area N (<name>) reviewed at <tier> tier (<model version>).
- Findings: <count> confirmed — <count> filed as issues #N-#M (still open), <count> fixed inline this round (logged, not filed), <count> candidates refuted in verification
- Trivial fixes: <count> applied inline
- Seam: <"still yielding — next visit stays at <tier>" | "ran dry at <tier>-<version> — next visit <re-runs <tier> on the newer version | escalates to <next tier>>">
- Next: area M (<name>) — run /red-team again when ready.
```

Link new findings by `file:line`. List inline fixes. Call out any high-severity signal —
that area escalates next time regardless of rotation. If step 5's work-shape gate fired,
report *that* instead of a round: what is owed, and the choice between doing it and
rotating on.
