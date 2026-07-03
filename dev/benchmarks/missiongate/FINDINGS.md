# Sectorial-escape mission gate — VERDICT: not a wall-clock lever (retires to budget axis)

Date 2026-07-03. Branch `claude/sect-missiongate` (off `cpp-search` c74ee6e6).
Hamilton array 17755546 (12 cells = 4 datasets x 3 TNT-T0 seeds), 5 sectorial-RNG
replicates per config => n=15 per (dataset, config). Clean uncontended cores.

## Question (the mission gate)
The distinct, non-budget residual of the TNT sectorial-escape thread (after D1/HTU-
float was refuted and the "geometry bulk" was shown to be mostly the budget confound
goofy-cannon owns) is: does the TNT-`selectem`-style config (large clade band +
sub-clade collapse + rasStarts=3 + rssPicks) escape a frozen shared TNT T0 at a
budget CHEAPER than the ratchet it would displace? Quality is closed WITH ratchet;
sectorial-escape only becomes a wall-clock win if it can escape cheaply enough to
displace ratchet (66-83% of production wall). Config-first (existing knobs via the
new `rssPicks` plumbing), ratchet/drift/xss/fuse OFF, from the identical TNT `mult` T0.

## Result: NO — 0/4 datasets escape cheaper than ratchet
`ratchet_ref` = default strategy, maxReplicates=1 (single ratchet start), same T0.
Escape = score - start (more negative = better). Wall = median seconds.

| dataset | target | ratchet_ref (esc @ wall) | best sectorial at <=ratchet wall | best sectorial ANY wall |
|---|---|---|---|---|
| Wortley2006 | -5 | **-8 @ 0.09s** | -7 | -7 (never reaches -8) |
| Zanol2014   | -13 | **-8 @ 0.54s** | -3 (@<=0.54s) | -7 @ 3.20s (6x wall) |
| Giles2015   | -4 | -3 @ 0.35s | -3 (tie) | -4 @ 0.64s (2x wall) |
| Zhu2013     | -8 | -1 @ 0.32s (noisy) | -1 (tie) | -6 @ 0.64s (2x wall) |

- **Wortley + Zanol: ratchet flatly Pareto-dominates.** A single ratchet start beats
  every sectorial config on BOTH score and wall. On Wortley ratchet reaches -8 (below
  the -5 target) in 0.09s; no sectorial config reaches -8 at all.
- **Giles + Zhu: the selectem `ras=3` machinery is NOT inert** — it yields reliable
  median escape (Giles -4 = target; Zhu -6) where plain baseline sectorial gets ~0 and
  a single ratchet start is a coin-flip (Zhu median -1, best -10). BUT it needs 2-7x
  ratchet's wall to do so. Not a wall win.

## What the data says about the sub-mechanisms
- **rasStarts=3 (the wagner-fixed RAS sector rebuild) is the only selectem ingredient
  that moves the needle** — on Giles/Zhu the ras=3 rows dominate the ras=1/baseline rows
  in median escape. Large-clade band + collapse + high rssPicks WITHOUT ras=3 ~ baseline.
- **Picks (rssPicks) and rounds buy escape only by buying budget** (monotone more
  candidates -> more escape, at proportional wall). This is the budget axis, not a
  distinct geometry lever — consistent with the freeze-big budget-confound finding.
- **Nothing reaches the TNT sectsch target reliably at competitive wall** (Zanol target
  -13; best sectorial median -7 at 6x wall; even best-of-15 only -13/-14 on lucky draws).

## Verdict
The selectem sectorial escape is REAL but Pareto-dominated by ratchet on wall-clock: it
cannot escape cheaper than the ratchet it would replace on any of the 4 datasets, so it
is not a wall-clock lever and cannot displace ratchet. The residual observation — ras=3
sectorial is more RELIABLE (lower variance) than a single ratchet start on Zhu/Giles — is
a variance/reliability point on the BUDGET axis (multi-start vs single-start), which is
goofy-cannon's territory, not a distinct algorithmic win. The true walk-up-selection src
routine is NOT justified: the config-first proxy already captures large-band+collapse and
still loses on wall by 2-6x; even a perfect walk-up would have to erase that deficit,
which nothing here suggests. THREAD RETIRES to the budget axis.

Plumbing kept: `rssPicks` SearchControl knob (exposes existing SectorParams field,
default 0 = byte-identical) — reusable for any future sequential-picks experiment.
