# Red-team escalation backlog — TreeSearch

Seams that are **re-eligible now** but are not next in the rotation. The `/red-team` skill
reopens entries here when a rung's model version moves on (`/red-team revisit`); a normal
rotation round may also pull an item from here when it reaches that area.

This file does **not** advance `last_focus:` and does not replace `log.md` — each item is a
pointer back to the round that recorded the residual. Full context stays in `log.md`; only
the *ask* and the *rung* live here.

**Rung rule.** Version bump before rung bump: `opus-4.8 dry → opus-5 → fable`. A dry verdict
is evidence about the version that ran, not about the rung. See the model-version legend at
the top of `log.md`.

---

## Open — five version-bump reopenings, plus one orphaned cross-area residual

The **five table rows** were reopened 2026-07-27 by the `opus` 4.8 → 5 version bump. Each was
recorded on 2026-07-24 (or earlier) as **"escalate to FABLE"**; the opus bump makes the same-rung
step available and cheaper, so each is re-queued at **`opus` (Opus 5) with a fresh-angle brief**.
Fable remains the step after Opus 5 *also* runs dry. The residuals themselves are unchanged —
only the rung aimed at them.

**Item 5(a), written up below the table, is a different animal** — not a version-bump reopening
but a residual that one area raised and routed to another, then vanished for six weeks because no
mechanism reads a cross-area note. Added by the area-9 round of 2026-07-27.

Per the skill's *Model versions* rule, a revisit **must change its angle, not just its
model**: brief the finder to attack the prior round's *derivation* (paste its "clean by
derivation" conclusions in as claims to break), not to re-read the same files the same way.

| # | Area | Rung | The residual (verbatim ask) | Recorded | Notes |
|---|------|------|-----------------------------|----------|-------|
| 1 | **1 — Fitch scoring correctness** | `opus` (Opus 5) | **Absolute** correctness of the shared `ts_fitch_combine` / full scorer. The L3b oracle proves incremental == from-scratch (*relative* equivalence) and is structurally blind to a bug both sides share; needs an **independent reference-scorer cross-check**. | log.md, area 1, 2026-07-24 (dry at opus-4.8) | Highest-value item in this file. Scope is ONLY this residual — area 1's *parked leads* (monomorphized EW/IW SPR scans, packing-in-its-own-right, IW x4-reroot, NA three-pass vs a reference NA scorer) are ordinary fresh-agent work, not this escalation. |
| 2 | **2 — Search topology invariants** | `opus` (Opus 5) | **Two** residuals, each with its own right check — carry BOTH. **(#1, oracle-blind, higher value):** absolute correctness of `compute_insertion_edge_sets` after the `00d73d6a` zero-fill removal; the oracle cannot see it (both sides now call the same zero-fill-skipped from-scratch) ⇒ confirmable only by score/reach-equivalence vs an **independent baseline** on homoplasy-rich data. **(#2, oracle-visible):** value-correctness of `update_base_after_spr_move` across many consecutive accepts on large real data (n_tip≥150, default-ON) ⇒ the Hamilton L3b oracle CAN confirm this. | log.md, area 2, 2026-07-24 (dry at opus-4.8 on the reject-invariant angle) | A static-only retirement of #1 is exactly what the escalation exists to double-check — do **not** re-aim only at #2. Both checks are Hamilton-class ([[feedback-no-local-heavy-compute]]). Blast radius of #1 is bounded to REACH by the full-rescore firewall, never a wrong returned score. |
| 3 | **3 — Ratchet & perturbation** | `opus` (Opus 5) | Exact numeric correctness of `compute_from_above_for_sector` / `build_ras_sector` heuristic proxies. | log.md, area 3, 2026-07-24 (dry at opus-4.8) | **THIN residual** — the whole new surface is opt-in/heuristic behind the full-dataset-rescore firewall, so it is reach-only, never a returned-score risk. Also a standing DOWNTIER/retire candidate, but *hold the retire decision until this Opus 5 pass*: three of area 3's four dry rounds are pre-tier and version-unrecorded, so only one counts toward the dormancy bar. |
| 4 | **5 — Data pipeline & simplification** | `opus` (Opus 5) | The ONE untested corner: does `TS_PACK_LOCAL` (default-ON) interact with `hierarchy_blocks` / Sankoff (HSJ/XFORM) tip data, which the Rcpp bridge sets **after** `build_dataset`? The packed `tip_states` layout vs HSJ tip data was never directly verified. | log.md, area 5, 2026-07-24 (dry at opus-4.8, evidenced empirically) | **THIN** — area 3's round showed rss/xss SKIP HSJ/XFORM, so this is likely moot; it is the only corner the 2026-07-24 trace did not reach. |
| 5 | **10 — Profile/IW/HSJ/XFORM kernels** | `opus` (Opus 5) | **(a) THE INAPPLICABLE IW MIN-STEPS CONVENTION — carry this first, it is the one item here with a hard numeric repro** (see the dedicated row below the table). **(b)** Un-spent area-10 surface: profile delta capping (`7cff7870`), `e/(k+e)` delta, `precomputed_steps` offset, `info_amounts` capping, `concavity=1.0` profile sentinel, DAT-002 `obs==0` XPIWE division reachability, and the OPEN clipped-subtree IW-screening follow-up. | (a) log.md, **area 9, 2026-06-16** (raised there, routed here, never picked up — see below); (b) log.md, area 10, 2026-06-16 (signal-resolution round; finder not spent) | **Weakest dry-evidence in this file.** Both cited dry rounds (2026-05-19, 2026-05-26) are pre-tier and version-unrecorded ⇒ area 10 has **zero** version-scoped dry verdicts. This is a first-versioned measurement, not a re-mine — and area 10 has never had a proper finder sweep. |

### Item 5(a) in full — the orphaned inapplicable-IW signal

Recorded here **2026-07-27** by the area-9 round, which found it had been sitting in a log entry
and nowhere else for six weeks. Unlike the five rows above, this is **not** a version-bump
reopening — see the note in *Not in this backlog* about why it belongs in this file anyway.

The 2026-06-16 **area-9** round raised a high-severity **cross-component** signal, explicitly
routed it to "a scoring area (≥opus)", and area 10 has not been visited since. It was never
filed in `findings.md`, never archived, never entered in `to-do.md`, and was not in this file
until now — so every mechanism that would have resurfaced it was bypassed at once. `tidy`
cannot catch this class: it reconciles *filed* findings against reality, and a residual that
only ever existed as log prose is invisible to it.

**The signal, verbatim in substance.** On Vinther2008 (inapplicable-bearing), at `concavity = 10`:

| Quantity | Value |
|---|---|
| Wagner/`AdditionTree` kernel NA+IW score | **3.003497** |
| the documented IW reference formula `Σ (cs−ml)/((cs−ml)+k)·w`, and `test-iw-scoring.R` | **3.003497** (exact match) |
| `TreeLength(tree, pd, concavity = 10)` | **2.974744** |
| EW step totals, both sides | 96 == 96 (agree) |

So the divergence is **purely the IW per-character minimum on inapplicable characters** — about
7 characters. Example given: char 23, `cs = 3`, `ml = 1` ⇒ reference `2/12 = 0.166667`,
`TreeLength` `0.16573`. On Lobo (also NA-bearing) the two agree, so it is
inapplicable-*pattern* dependent, not NA-dependent.

**The ask.** Settle which convention is correct for the minimum steps of an inapplicable
character under implied weights, then make the two paths agree. The 2026-06-16 conclusion was
that the *kernel* is correct per the package's stated contract and needs no change — the root
cause is R-layer `TreeLength` / `MinimumLength` / `CharacterLength`. Treat that as the prior,
not as settled: it was a round's judgement call, not a verified derivation.

**Why it matters beyond a rounding difference.** If `TreeLength` is the user-facing ground
truth, then `MaximizeParsimony` under NA+IW optimises and reports a subtly *different objective*
than `TreeLength` for affected datasets — two numbers for one tree, with no indication which the
user should believe.

**Test coverage is the reason it stayed invisible:** `test-iw-scoring.R` does assert
`TreeLength == reference`, but **only on Lobo** — the one dataset of the two where they agree.
The divergent case is uncovered, so no test will ever fail because of this.

**One arithmetic pointer, computed here from the logged numbers — a lead, not a finding.**
Inverting `x/(x+10) = 0.16573` gives `x ≈ 1.9865`, not 2: the effective `cs−ml` behind
`TreeLength`'s value is **fractional**, which points at `MinimumLength` returning a non-integer
for inapplicable characters rather than at the IW formula itself. Also note char 23 is **not
representative** — its per-character gap is `0.000937`, while the whole-tree gap is `0.028753`,
~31× larger, so other divergent characters dominate and should be identified before theorising.

**Scope gap this exposes, worth fixing in a future area-12 round:** the named root-cause files
(`TreeLength` / `MinimumLength` / `CharacterLength`, R-layer IW scoring) are **not in area 10's
scope row**, which lists only `src/ts_fitch.cpp` IW paths, `src/ts_data.cpp`, `src/ts_hsj.*` and
`src/ts_sankoff.*`. They are not in any other area's row either. This is a concrete instance of
the `R/*.R` coverage gap area 12 flagged on 2026-07-03 but never diffed.

### Item 6 — `expand_and_reinsert` ignores its `ConstraintData*` (UNVERIFIED lead)

Recorded 2026-07-27 by the area-9 round. **Cross-area, so it is here rather than only in
`log.md`** — the rule this file adopted in the same round, applied to its first case: area 9
found it, areas 3 and 13 own it, and a note left only in area 9's log entry would be read by
nobody.

`src/ts_prune_reinsert.cpp:353` — `expand_and_reinsert(…, ts::ConstraintData* cd)` takes a
constraint pointer and **never references it**, which is why it shows up as a pre-existing
`-Wunused-parameter` under `g++ -Wall -Wextra`. If that function genuinely re-inserts tips
without consulting the constraint, it is a [`T-324`](findings.md)-shaped gap on a different
path — reinsertion producing a violating tree that only the downstream posthoc check might
catch.

**Status: UNVERIFIED, and it must not be treated as a finding until someone reads the
function.** Three things a reader should settle: whether reinsertion is genuinely
constraint-blind or the constraint is enforced by a caller; whether any capture path can retain
a violating tree; and whether the parameter is simply vestigial (in which case deleting it is
the fix, and the warning goes away).

Two reasons it is worth someone's time rather than a shrug. It surfaced from a **compiler
warning, not from reading** — nobody has read this function's constraint handling, so its
silence is not evidence. And the same function is already the subject of
[`T-366`](findings.md) (mixed-regime `prelim`), so a reader is going in there anyway; settling
both in one pass costs barely more than settling one.

### Not in this backlog (deliberately)

- **Area 4 (Parallelism & RNG), 6 (R↔C++), 7 (Shiny), 8 (Tests), 9 (Wagner), 11 (Collapse),
  12 (Meta), 13 (Constraints)** — all recorded **still yielding** at their last visit, so they
  escalate nothing; they re-visit at the same rung with a fresh agent when rotation reaches
  them. A version bump does not reopen a yielding seam (there is no dormancy to reopen).
- **Area 13's next visit** is a *bounded exhaustive harness*, not a finder at any rung
  (log.md, area 13, 2026-07-03) — a work-shape decision the version bump does not change.
- **The `sonnet` 4.6 → 5 and `fable` bumps of the same day reopened NOTHING** — checked by
  `tidy` 2026-07-27, recorded here so `/red-team revisit` need not redo it. A version bump only
  reopens a *dry* verdict, and every sonnet-tier area (6, 7, 8, 12) was recorded **still
  yielding** at its last visit; area 6 has in any case already escalated past sonnet (its
  2026-07-24 round ran opus + fable dual-tier and yielded T-339/T-340, both since fixed). So
  the opus bump is the only one with re-eligible seams behind it, and they are the five above.
- **Area 8's residual from the 2026-07-27 round is deliberately NOT here.** It is a
  *within-area, same-rung* next-round note (extend mutation testing to the rest of the
  per-finding regression corpus; settle whether `test-ts-hsj.R:352`'s smoke test would catch a
  *wrong* `absent_state` rather than only a crash), so `log.md` is its home.
- **Two classes of residual belong in this file, not one.** The original class is a seam
  re-eligible at a **different rung** than rotation would give it (rows 1–5's version-bump
  reopenings). Item **5(a)** adds the second: a **cross-area-routed** residual, where the area
  that *found* it is not the area that *owns* it. That routing is precisely what makes it get
  lost — a within-area note is read by the next round on that area, but a note saying "some
  other area should look at this" is read by nobody, because `log.md` entries are consulted
  per-area. Six weeks of invisibility is the evidence. **The test for whether a residual belongs
  here is therefore "would the next round that should act on it actually read it?", not "is a
  rung change involved?"**

---

## Resolved history (one line each)

*(Nothing resolved yet — this file was created 2026-07-27 with the opus 4.8 → 5 bump.)*
