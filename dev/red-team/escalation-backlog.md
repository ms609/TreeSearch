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

## Open — reopened 2026-07-27 by the `opus` 4.8 → 5 version bump

All five items below were recorded on 2026-07-24 (or earlier) as **"escalate to FABLE"**.
The opus bump makes the same-rung step available and cheaper, so each is re-queued at
**`opus` (Opus 5) with a fresh-angle brief**. Fable remains the step after Opus 5 *also* runs
dry. The residuals themselves are unchanged — only the rung aimed at them.

Per the skill's *Model versions* rule, a revisit **must change its angle, not just its
model**: brief the finder to attack the prior round's *derivation* (paste its "clean by
derivation" conclusions in as claims to break), not to re-read the same files the same way.

| # | Area | Rung | The residual (verbatim ask) | Recorded | Notes |
|---|------|------|-----------------------------|----------|-------|
| 1 | **1 — Fitch scoring correctness** | `opus` (Opus 5) | **Absolute** correctness of the shared `ts_fitch_combine` / full scorer. The L3b oracle proves incremental == from-scratch (*relative* equivalence) and is structurally blind to a bug both sides share; needs an **independent reference-scorer cross-check**. | log.md, area 1, 2026-07-24 (dry at opus-4.8) | Highest-value item in this file. Scope is ONLY this residual — area 1's *parked leads* (monomorphized EW/IW SPR scans, packing-in-its-own-right, IW x4-reroot, NA three-pass vs a reference NA scorer) are ordinary fresh-agent work, not this escalation. |
| 2 | **2 — Search topology invariants** | `opus` (Opus 5) | **Two** residuals, each with its own right check — carry BOTH. **(#1, oracle-blind, higher value):** absolute correctness of `compute_insertion_edge_sets` after the `00d73d6a` zero-fill removal; the oracle cannot see it (both sides now call the same zero-fill-skipped from-scratch) ⇒ confirmable only by score/reach-equivalence vs an **independent baseline** on homoplasy-rich data. **(#2, oracle-visible):** value-correctness of `update_base_after_spr_move` across many consecutive accepts on large real data (n_tip≥150, default-ON) ⇒ the Hamilton L3b oracle CAN confirm this. | log.md, area 2, 2026-07-24 (dry at opus-4.8 on the reject-invariant angle) | A static-only retirement of #1 is exactly what the escalation exists to double-check — do **not** re-aim only at #2. Both checks are Hamilton-class ([[feedback-no-local-heavy-compute]]). Blast radius of #1 is bounded to REACH by the full-rescore firewall, never a wrong returned score. |
| 3 | **3 — Ratchet & perturbation** | `opus` (Opus 5) | Exact numeric correctness of `compute_from_above_for_sector` / `build_ras_sector` heuristic proxies. | log.md, area 3, 2026-07-24 (dry at opus-4.8) | **THIN residual** — the whole new surface is opt-in/heuristic behind the full-dataset-rescore firewall, so it is reach-only, never a returned-score risk. Also a standing DOWNTIER/retire candidate, but *hold the retire decision until this Opus 5 pass*: three of area 3's four dry rounds are pre-tier and version-unrecorded, so only one counts toward the dormancy bar. |
| 4 | **5 — Data pipeline & simplification** | `opus` (Opus 5) | The ONE untested corner: does `TS_PACK_LOCAL` (default-ON) interact with `hierarchy_blocks` / Sankoff (HSJ/XFORM) tip data, which the Rcpp bridge sets **after** `build_dataset`? The packed `tip_states` layout vs HSJ tip data was never directly verified. | log.md, area 5, 2026-07-24 (dry at opus-4.8, evidenced empirically) | **THIN** — area 3's round showed rss/xss SKIP HSJ/XFORM, so this is likely moot; it is the only corner the 2026-07-24 trace did not reach. |
| 5 | **10 — Profile/IW/HSJ/XFORM kernels** | `opus` (Opus 5) | Un-spent area-10 surface: profile delta capping (`7cff7870`), `e/(k+e)` delta, `precomputed_steps` offset, `info_amounts` capping, `concavity=1.0` profile sentinel, DAT-002 `obs==0` XPIWE division reachability, and the OPEN clipped-subtree IW-screening follow-up. | log.md, area 10, 2026-06-16 (signal-resolution round; finder not spent) | **Weakest dry-evidence in this file.** Both cited dry rounds (2026-05-19, 2026-05-26) are pre-tier and version-unrecorded ⇒ area 10 has **zero** version-scoped dry verdicts. This is a first-versioned measurement, not a re-mine — and area 10 has never had a proper finder sweep. |

### Not in this backlog (deliberately)

- **Area 4 (Parallelism & RNG), 6 (R↔C++), 7 (Shiny), 8 (Tests), 9 (Wagner), 11 (Collapse),
  12 (Meta), 13 (Constraints)** — all recorded **still yielding** at their last visit, so they
  escalate nothing; they re-visit at the same rung with a fresh agent when rotation reaches
  them. A version bump does not reopen a yielding seam (there is no dormancy to reopen).
- **Area 13's next visit** is a *bounded exhaustive harness*, not a finder at any rung
  (log.md, area 13, 2026-07-03) — a work-shape decision the version bump does not change.

---

## Resolved history (one line each)

*(Nothing resolved yet — this file was created 2026-07-27 with the opus 4.8 → 5 bump.)*
