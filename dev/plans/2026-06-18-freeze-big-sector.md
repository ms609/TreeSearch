# Freeze-big sector reduction breaks the ratchet-off null (TNT selectem mechanism)

Date 2026-06-18. Branch `claude/selectem-diversity` (worktree `C:/Users/pjjg18/GitHub/TS-selectem`,
off `cpp-search`). Env-gated, default path byte-identical. NOT yet on cpp-search.

## Question
Can TreeSearch's RSS sectorial, with **global ratchet/drift OFF**, escape the canonical
frozen T0 the way TNT does? (TNT reaches its scores with ratchet OFF — so matching this is a
**wall-clock** lever, not only quality. See [[tnt-sectorial-recipe]].)

## Target (define_target.R; ratchet-off TNT `mult`+`sectsch=rss`, canonical hold-1000 T0)
| dataset | n | T0 | TNT ratchet-off sectorial | escape |
|---|---|---|---|---|
| Zanol2014  | 74 | 1271 | 1261 | −10 |
| Wortley2006| 37 |  485 |  480 | −5 |
| Zhu2013    | 75 |  631 |  624 | −7 |
| Giles2015  | 78 |  672 |  670 | −2 |
TreeSearch baseline from the SAME T0 (ts_arms.R, base/coll30): **~0** (Zhu −2 only). Stuck.

## Mechanism found: freeze-big sector reduction
The existing collapse (`build_reduced_dataset_collapsed`) is **break-big**: expands the largest
sub-clade until `target_tips` units → surviving composites are SMALL leftovers → no large
movable units → null (coll30). TNT `selectem` (Goloboff 1999 App.1) is **freeze-big**: keep
tips individual, FREEZE whole sub-clades (≥0.8·cap) into single composite terminals (random
order) until ≤cap units. Relocating one composite = transplanting a multi-taxon clade as a
**unit** — a large-radius move single-step hill-climbing cannot reach by moving fragments one
at a time. Implemented `build_reduced_dataset_freeze` (env `TS_FREEZE_COLLAPSE`; cap/thresh =
`TS_FREEZE_CAP`/`TS_FREEZE_THRESH` knobs; `TS_FREEZE_RANDOM` = random vs deterministic order;
shares `assemble_reduced` with break-big, refactor verified byte-identical).

Critical tuning: a USEFUL reduction (many tips + ONE big composite) needs a HIGH freeze
threshold (only freeze near-cap sub-clades, so a freeze *overshoots* cap). Low threshold →
degenerate (one near-whole composite, e.g. clade=74→units=2 maxcomp=73). cap=33 thr=28 worked
for n≈74; must scale with n.

## Result (ts_arms.R, canonical T0, ratchet/drift OFF, rss-only, 30 rounds, 20 picks, seeds 1-3)
| dataset | target | base/coll30 | freezeHT**det** (H2) | freezeHT**rand** (H1) |
|---|---|---|---|---|
| Zanol | −10 | 0 / 0 | **−4** | **−5** |
| Zhu   | −7  | 0 / −2 | **−2** | **−3** |
(cap33 thr28; Wortley/Giles need n-scaled cap — see freezeScaled run.)

## H1 vs H2 — ANSWERED
Deterministic high-threshold freeze ALREADY breaks the null (Zanol −4, Zhu −2). Randomisation
adds ~+1 step (−5, −3). So per the pre-registered ablation: **H2 (large movable units) is the
PRIMARY lever; per-pass diversity (H1) is a secondary ~+1 increment.** The user's "per-pass
diversity" intuition is real but not the main thing — the missing ingredient was the structural
move type (frozen-clade-as-unit), available even deterministically.

## n-scaled run (negative field = pct of n; cap .45n thr .38n band [.42n,.99n])
| dataset | target | freezeScalDet | freezeScaled(rand) |
|---|---|---|---|
| Zanol | −10 | −4 | −5 |
| Zhu   | −7  | 0  | −1 |
| Wortley | −5 | 0 | 0 |
| Giles | −2  | 0 | 0 |
n-scaling did NOT unlock Wortley/Giles and made Zhu noisier (signal is 1-5 steps, seed-sensitive).

## BUDGET CONFOUND — freeze framing COLLAPSES (advisor-caught)
The "freeze breaks the null" claim was confounded: the null (base/coll30) ran at default
~5 picks x 15 rounds (~75 sector searches); the freeze arms ran 20 picks x 30 rounds (~600).
8x more search. Decisive budget-matched run (ALL at 20 picks x 30 rounds, seeds 1-3):
| dataset | target | base20 (small [6,50]) | coll30_20 (break-big) | freezeHT (freeze-big) |
|---|---|---|---|---|
| Zanol | −10 | 0 | **−4** | −5 |
| Zhu   | −7  | 0 | **−4** | −3 |
| Wortley | −5 | 0 | 0 | 0 |
| Giles | −2  | 0 | 0 | 0 |
**coll30 (plain break-big) reaches −4 at matched budget == freeze (±1, freeze WORSE on Zhu).**
So freeze-big / "large movable units" (H2) and "per-pass diversity" (H1) are NOT the lever —
they add nothing over the pre-existing break-big collapse. `build_reduced_dataset_freeze` adds
no value; KEEP IT OUT of production.

## What ACTUALLY breaks the ratchet-off null (corrected)
- `base20` (small-clade [6,50] selection) finds 0 even at 600 searches → small-clade sectorial
  is the dead end, regardless of budget.
- `coll30_20` (LARGE-clade [31,99] selection + collapse + ras3 + 600 searches) finds −4 on
  Zanol AND Zhu → **large-clade selection + sufficient budget** is the partial lever (~half the
  gap). The earlier coll30 "null" (−2 Zhu) was itself a LOW-BUDGET artifact.
- This also implies the memory's "RAS-multistart on large sectors = null" was likely a budget
  artifact (it was at default low budget). [[tnt-sectorial-recipe]] selection-quality verdicts
  need re-reading through the budget lens.
- Still 0 on Wortley/Giles even at high budget — genuinely unresponsive from this T0.
- (isolation run b42fxd23o: does large-selection ALONE escape, or need collapse/ras3? — fill in)

## Clean surviving finding
The H1/H2 ablation (freezeHTdet −4 vs rand −5, both 30x20) is internally clean BUT moot now
that coll30 (deterministic break-big) also = −4: freeze machinery is within noise of break-big.
The honest result: **TreeSearch's existing large-clade collapse sectorial DOES partially escape
the frozen T0 (~half the Zanol/Zhu gap) once given TNT-like budget (≈20 picks/round, ≈30
rounds) — but only on large-clade datasets, and the remaining gap + Wortley/Giles are open.**
