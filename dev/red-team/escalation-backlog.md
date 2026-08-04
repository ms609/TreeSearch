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

## Open — five version-bump reopenings, plus one cross-area residual

The **five table rows** were reopened 2026-07-27 by the `opus` 4.8 → 5 version bump. Each was
recorded on 2026-07-24 (or earlier) as **"escalate to FABLE"**; the opus bump makes the same-rung
step available and cheaper, so each is re-queued at **`opus` (Opus 5) with a fresh-angle brief**.
Fable remains the step after Opus 5 *also* runs dry. The residuals themselves are unchanged —
only the rung aimed at them.

**Item 5(a) was retired 2026-07-28** — the bug claim it carried had already been refuted, in
`log.md`, by the receiving area on the same day it was raised. It is written up in *Resolved
history*; two narrow sub-items survive and are noted under item 5(a) below. It was the first
entry of the second class this file admits (a **cross-area-routed** residual rather than a
version-bump reopening), and retiring it does not retire that class — item 6 is one.

Per the skill's *Model versions* rule, a revisit **must change its angle, not just its
model**: brief the finder to attack the prior round's *derivation* (paste its "clean by
derivation" conclusions in as claims to break), not to re-read the same files the same way.

| # | Area | Rung | The residual (verbatim ask) | Recorded | Notes |
|---|------|------|-----------------------------|----------|-------|
| 1 | **1 — Fitch scoring correctness** | `opus` (Opus 5) | **Absolute** correctness of the shared `ts_fitch_combine` / full scorer. The L3b oracle proves incremental == from-scratch (*relative* equivalence) and is structurally blind to a bug both sides share; needs an **independent reference-scorer cross-check**. | log.md, area 1, 2026-07-24 (dry at opus-4.8) | Highest-value item in this file. Scope is ONLY this residual — area 1's *parked leads* (monomorphized EW/IW SPR scans, packing-in-its-own-right, IW x4-reroot, NA three-pass vs a reference NA scorer) are ordinary fresh-agent work, not this escalation. |
| 2 | **2 — Search topology invariants** | `opus` (Opus 5) | **Two** residuals, each with its own right check — carry BOTH. **(#1, oracle-blind, higher value):** absolute correctness of `compute_insertion_edge_sets` after the `00d73d6a` zero-fill removal; the oracle cannot see it (both sides now call the same zero-fill-skipped from-scratch) ⇒ confirmable only by score/reach-equivalence vs an **independent baseline** on homoplasy-rich data. **(#2, oracle-visible):** value-correctness of `update_base_after_spr_move` across many consecutive accepts on large real data (n_tip≥150, default-ON) ⇒ the Hamilton L3b oracle CAN confirm this. | log.md, area 2, 2026-07-24 (dry at opus-4.8 on the reject-invariant angle) | A static-only retirement of #1 is exactly what the escalation exists to double-check — do **not** re-aim only at #2. Both checks are Hamilton-class ([[feedback-no-local-heavy-compute]]). Blast radius of #1 is bounded to REACH by the full-rescore firewall, never a wrong returned score. |
| 3 | **3 — Ratchet & perturbation** | `opus` (Opus 5) | Exact numeric correctness of `compute_from_above_for_sector` / `build_ras_sector` heuristic proxies. | log.md, area 3, 2026-07-24 (dry at opus-4.8) | **THIN residual** — the whole new surface is opt-in/heuristic behind the full-dataset-rescore firewall, so it is reach-only, never a returned-score risk. Also a standing DOWNTIER/retire candidate, but *hold the retire decision until this Opus 5 pass*: three of area 3's four dry rounds are pre-tier and version-unrecorded, so only one counts toward the dormancy bar. |
| 4 | **5 — Data pipeline & simplification** | `opus` (Opus 5) | ~~Does `TS_PACK_LOCAL` (default-ON) interact with `hierarchy_blocks` / Sankoff (HSJ/XFORM) tip data?~~ **SCORING-SIDE HALF CLOSED 2026-07-28** by the area-10 round: `TS_PACK_LOCAL=0` vs default is bit-identical for PROFILE, and HSJ/XFORM are identical across the full 2×2 with `TS_CHAR_ORDER`, with a reasoned argument so it need not be re-run — profile depends only on `char_steps[]`, and HSJ reads `ds.tip_labels` while XFORM reads `ds.sankoff_tip_costs`, **never** `ds.tip_states`, so plane relabelling cannot reach either. **What remains for area 5 is only the tip-data-construction half**: whether the Rcpp bridge, which sets HSJ/Sankoff tip data *after* `build_dataset`, is itself consistent with the packed layout. | log.md, area 5, 2026-07-24 (dry at opus-4.8, evidenced empirically); scoring half closed log.md, area 10, 2026-07-28 | **THINNER THAN BEFORE** — area 3's round showed rss/xss SKIP HSJ/XFORM, and the scoring side is now positively verified. |
| 5 | **10 — Profile/IW/HSJ/XFORM kernels** | — | **SPENT 2026-07-28 — REMOVE ON NEXT `tidy`.** The Opus 5 round ran and yielded **4 P1s + 2 P2s + 4 P3s** (T-373…T-382). Item (a) was retired pre-dispatch as already-refuted (see Resolved history). Item (b) is fully worked: every static residual (profile delta capping, `e/(k+e)`, `precomputed_steps` offset, `info_amounts` capping, `concavity=1.0` sentinel, DAT-002 `obs==0` reachability) was **independently re-derived and confirmed clean**, and the char-ordering × per-pattern-indexing question came back clean on a proven non-identity permutation. **One item did NOT get reached and stays open: the clipped-subtree IW-screening follow-up.** | log.md, area 10, 2026-07-28 | The dry-evidence caveat proved to be the right call and then some: area 10 had zero version-scoped dry verdicts, and the reason was structural — `ts_hsj.*`/`ts_sankoff.*` were added to the scope row on 2026-07-03 but **no finder had ever read them at any tier**, so both prior dry rounds were dry about a *different half* of the area. Every P1 came from that unread half. Area 10 escalates next round anyway on T-374's open design question. |

### Item 5(a) — RETIRED 2026-07-28. The bug claim was already refuted; two sub-items survive.

**Do not carry the bug claim into a finder brief.** It was resolved on 2026-06-16 by the
area-10 signal-resolution round, and the area-9 round of 2026-07-27 resurrected it because it
read the area-9 log entry that *raised* the signal but not the area-10 entry that *killed* it.
Full write-up in *Resolved history* at the bottom of this file; the original text is preserved
there so the anti-duplication value is not lost.

The two things that genuinely survive, both carried into the 2026-07-28 area-10 round:

1. **A test-coverage gap — real, but much narrower than the 2026-07-27 write-up claimed.** That
   write-up asserted "test coverage is the reason it stayed invisible: `test-iw-scoring.R` does
   assert `TreeLength == reference`, but only on Lobo." **That is wrong on both halves, checked
   2026-07-28.** `test-iw-scoring.R:20-21` passes `extended_iw = FALSE`, so it pins the *plain-IW*
   path and never exercises XPIWE at all — Lobo is beside the point. And a dedicated
   `test-ts-xpiwe.R` exists that the write-up did not know about: it pins the XPIWE formula
   against a hand-computed reference (`:62-91`, tolerance 1e-10, computing `f`, `eff_k`, `phi`
   from first principles), asserts XPIWE == IW on complete data (`:38-48` — the identity property
   derived above), and asserts they differ with missing entries (`:50-60`).

   **What is genuinely uncovered:** every XPIWE reference test uses `make_missing_data()`
   (`:18-32`), whose only non-observed token is `?`. **No test pins XPIWE against a reference on
   inapplicable-bearing (`-`) input** — which is exactly the input class the original signal was
   observed on. The substantive open question is narrow and answerable: does `.ObsCount` treat
   `-` as observed, and is that the intended convention for `f = 1 + r·(nTaxa − obs)/obs`? Adding
   a `-`-bearing case to `test-ts-xpiwe.R` is a small, well-defined job.
2. **A scope gap (area 12's, not area 10's).** `TreeLength` / `MinimumLength` /
   `CharacterLength` and the R-layer IW scoring are in **no area's scope row**. This remains a
   concrete instance of the `R/*.R` coverage gap area 12 flagged on 2026-07-03 and never diffed.

**The process lesson is the durable part, and it is the opposite of what the 2026-07-27 entry
concluded.** That round diagnosed "no mechanism reads a cross-area note" and built this file to
fix it. That diagnosis was half right: the routing mechanism *was* missing, but in this instance
the receiving area **had already acted** — the resolution existed, in `log.md`, under the very
area the signal was routed to. The actual failure was reading only the entry that raised a
signal and not searching for a later entry that resolved it. **So the rule this file needs is
not only "route cross-area residuals here" but "before re-queuing any residual, grep `log.md`
for later rounds on the *receiving* area."** A backlog that re-queues refuted work is worse than
no backlog, because it spends premium-tier budget with the authority of a tracked item.

The signal's full original text — the Vinther2008 numbers, the `0.16573` char-23 example, and
the "fractional `cs−ml`" lead — is preserved verbatim in *Resolved history*, together with the
arithmetic that retires it. Keeping it readable is the point: a future finder that re-derives
those same numbers should be able to find, in one grep, why they are not a bug.

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

### Item 7 — area 13 gets two filed constraint findings from an area-11 round, one of them P1

Recorded 2026-08-04 by the area-11 round. **Cross-area class** (the second one this file admits):
area 11 found them, area 13 owns them, and area 13's *recorded next-visit plan predates them*.

Receiving-area check done as this file requires: area 13's most recent round is **2026-07-03**,
and nothing later in `log.md` touches either finding. So this is genuinely open, not a re-queue
of resolved work.

**The two findings.** [`T-402`](findings.md) (**P1**) — a `constraint` is silently ignored when
the caller supplies a violating start via `tree =`; the search freezes on it, reports a
better-than-constrained score, and *evicts* every compliant tree other replicates find.
[`T-403`](findings.md) (P2) — the "enforced splits are protected from collapse" promise is an
exact-match test with no access to `consZero`, so under the **default** `collapse = TRUE` the
returned trees can violate the constraint outright (20/20 seeds).

**The ask is a sequencing decision, not a review.** Area 13's next visit was recorded as *"a
bounded exhaustive harness, not a finder"* (the `topology_spr` / `build_postorder`-guard
equivalence). That plan is orthogonal to these two and still stands on its merits — but it was
set when area 13 had no filed P1. Whoever takes area 13 next should decide explicitly which
comes first and record the reason, rather than defaulting to the older note.

**Two things to read before patching anything in this class**, both already in the rows:

1. **A verify-and-revert gate of the T-390/T-391 shape does not fix T-402.** `nni_perturb`
   snapshots the violating start *before* repair and then rejects the repaired legal tree for
   scoring worse, so the illegal score is an unbeatable baseline. Gating the pool capture alone
   is worse than useless: the pool empties at `maxReplicates = 1` and
   `R/MaximizeParsimony.R:1682-1684` returns the user's violating start anyway. The fix has to
   act at the `startEdge` boundary.
2. **T-324's row was amended on 2026-08-04** because its repair claim was over-optimistic in
   exactly this regime. T-402 and T-324 share T-324's downstream half verbatim (ungated pool
   capture, no downstream filter), so they should be fixed together — with T-402's deterministic
   8-taxon repro as the standing regression test for the shared half. **T-402 does not settle
   T-324's own reachability question**, and neither row should be read as if it does.

**One part is a maintainer adjudication, not a fixer's call** (same shape as T-396): whether the
`startEdge` boundary should *repair* a violating start, or *reject* it with an error/warning.
Both satisfy the contract; they differ in whether `tree =` stays usable as a warm start under a
constraint, which is a user-facing design choice.

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
  reopenings). The second is a **cross-area-routed** residual, where the area that *found* it is
  not the area that *owns* it — a within-area note is read by the next round on that area, but a
  note saying "some other area should look at this" is read by nobody, because `log.md` entries
  are consulted per-area. Item **6** is the live instance. **The test for whether a residual
  belongs here is "would the next round that should act on it actually read it?", not "is a rung
  change involved?"**
- **Corollary added 2026-07-28, and it is the harder half: entering a residual here requires
  checking that it is still open.** Item 5(a) was the founding example of the cross-area class
  and it turned out to be *already resolved* — refuted the same day it was raised, by the very
  area it was routed to, in an entry sitting in `log.md` under that area's own heading. Nobody
  read it, because the round that re-queued the residual read the entry that *raised* the signal
  and stopped there. **So: before adding or re-queuing any residual, grep `log.md` for later
  rounds on the RECEIVING area, not just the raising one.** A backlog that re-queues refuted work
  is worse than no backlog — it spends premium-tier budget carrying the authority of a tracked
  item, and 5(a) was flagged "carry this first, it is the one item here with a hard numeric
  repro."

---

## Resolved history (one line each)

- **Item 5(a) — inapplicable IW min-steps convention. RETIRED 2026-07-28** as a bug claim; it was
  never live. Resolved 2026-06-16 by the area-10 signal-resolution round (plain-IW vs XPIWE);
  re-queued in error 2026-07-27; arithmetic independently re-confirmed 2026-07-28. Two sub-items
  survive: a narrow XPIWE test gap (no `-`-bearing reference case — the write-up's broader
  coverage claim does **not** survive, `test-ts-xpiwe.R` already pins the formula) and area 12's
  R-layer scope gap. Full write-up below.

### Item 5(a) in full — why it is retired, and the original text

**Why it was re-queued.** The 2026-06-16 **area-9** round raised this as a high-severity
cross-component signal and routed it to "a scoring area (≥opus)". The 2026-07-27 area-9 round
found it in a log entry and nowhere else, concluded that area 10 "has not been visited since",
and created this file to stop such residuals getting lost. **But area 10 *had* been visited —
the same day, 2026-06-16, in a signal-resolution round that refuted this exact signal** with
three exact equalities. The re-queue read the entry that raised the signal and not the one that
killed it.

**The refutation (2026-06-16, area 10).** It is **plain-IW vs XPIWE**, not a min-steps
convention. `TreeLength()` defaults to `extended_iw = TRUE` (XPIWE — Goloboff 2014 Extension-3
missing-data correction, `R/tree_length.R:144-156`: `f = 1 + r·(nTaxa−obs)/obs`, `eff_k = k/f`,
`phi = (1+eff_k)/(1+k)`, `fit = h/(h+eff_k)`, `Σ fit·w·phi`), whereas the area-9 cross-check
called the kernel via `ts_fitch_score(..., min_steps, concavity)` — plain IW, no XPIWE args. Two
different objectives by construction. The round further confirmed (i) `TreeLength(extended_iw =
FALSE)` == kernel plain IW exactly, and (ii) kernel XPIWE == an independent `TreeLength()` XPIWE
rescore of a real `MaximizeParsimony(concavity = 10)` result, exactly — so production optimises
and reports the *same* objective in both modes. There is no optimise-vs-report mismatch.

**The arithmetic, re-confirmed independently 2026-07-28.** The 2026-07-27 re-queue treated the
"non-rational" `0.16573` as its strongest evidence, reasoning that inverting `x/(x+10) = 0.16573`
gives `x ≈ 1.9865` rather than 2, and that a fractional `cs−ml` points at `MinimumLength`. That
inversion assumes plain IW. Under XPIWE with `h = 2`, `k = 10`, `eff_k ≈ 9.3`:

```
fit·phi = (2 / (2 + 9.3)) · ((1 + 9.3) / (1 + 10)) = 0.176991 · 0.936364 = 0.165728
```

against the logged `0.16573` — agreement to five significant figures, with `h` an exact
integer 2. The non-rational value is the `phi`/`eff_k` scaling, exactly as the 2026-06-16 round
said. Note also that at `eff_k = 10` (no missing data ⇒ `f = 1`) the expression collapses to
`(2/12)·(11/11) = 0.166667` = the plain-IW reference **exactly** — i.e. the XPIWE correction is
the *identity* when nothing is missing. That is the mechanism behind the one fact the re-queue
found most suspicious: Lobo agrees and Vinther2008 does not.

**Original text, preserved verbatim for anti-duplication.** On Vinther2008
(inapplicable-bearing), at `concavity = 10`:

| Quantity | Value |
|---|---|
| Wagner/`AdditionTree` kernel NA+IW score | **3.003497** |
| the documented IW reference formula `Σ (cs−ml)/((cs−ml)+k)·w`, and `test-iw-scoring.R` | **3.003497** (exact match) |
| `TreeLength(tree, pd, concavity = 10)` | **2.974744** |
| EW step totals, both sides | 96 == 96 (agree) |

> So the divergence is **purely the IW per-character minimum on inapplicable characters** — about
> 7 characters. Example given: char 23, `cs = 3`, `ml = 1` ⇒ reference `2/12 = 0.166667`,
> `TreeLength` `0.16573`. On Lobo (also NA-bearing) the two agree, so it is
> inapplicable-*pattern* dependent, not NA-dependent.
>
> **One arithmetic pointer, computed here from the logged numbers — a lead, not a finding.**
> Inverting `x/(x+10) = 0.16573` gives `x ≈ 1.9865`, not 2: the effective `cs−ml` behind
> `TreeLength`'s value is **fractional**, which points at `MinimumLength` returning a non-integer
> for inapplicable characters rather than at the IW formula itself.

The EW row is the tell that was available at the time and went unused: **96 == 96**. A
min-steps-convention bug on inapplicable characters would move the EW step totals too. Only a
weighting correction leaves EW identical and IW divergent.

**A second claim in the original that does not survive checking.** The write-up asserted that
thin test coverage was "the reason it stayed invisible." `test-ts-xpiwe.R` — which the write-up
never mentions — pins the XPIWE formula against a first-principles hand computation and encodes
the very identity property that explains the Lobo/Vinther2008 asymmetry. The coverage was there;
it simply was not looked for. The narrow, genuine gap that remains (no `-`-bearing XPIWE
reference case) is recorded in the Open section above.
