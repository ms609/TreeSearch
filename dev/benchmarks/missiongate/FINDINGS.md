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

## First, the diagnostic result: the MECHANISM gap is CLOSED
This thread was seeded by `tnt-outperformance-is-diagnostic`: from an identical T0,
ours historically escaped **0** while TNT sectsch escaped +3..+11. That per-iteration
gap is now GONE. Best-of-15 sectorial escape vs the TNT sectsch target:

| dataset | TNT target | ours esc_best (was 0) | reaches target? |
|---|---|---|---|
| Zanol2014   | -13 | **-13** | yes |
| Wortley2006 | -5  | **-15** | yes (overshoots) |
| Zhu2013     | -8  | **-11** | yes (beats) |
| Giles2015   | -4  | **-5**  | yes (beats) |

So ours now escapes too — the wagner-fixed RAS sector rebuild (`93071cae`) + rasStarts>=3
supply the barrier-crossing move the audit chased under the wrong mechanism (D1). The
diagnostic clue is RESOLVED: not a missing capability. NB reachable (best-of-15), NOT
reliable (Zanol median only -7) — closure is of the mechanism, not of variance.

## Then, the wall-clock question: is it a LEVER? NO — 0/4 escape cheaper than ratchet
This is NOT ROI-gating the TNT win (that win is resolved above); it asks whether the
now-working sectorial escape can DISPLACE ratchet on wall-clock. Lead on Mcand
(machine-portable, per `policy-in-replicates-not-seconds`); competitive-end walls are
sub-second where R/setup overhead swamps ratios. `ratchet_ref` = default strategy,
maxReplicates=1 (single ratchet start), same T0. Escape = score - start.

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
- **Escape-per-candidate settles it (machine-portable).** Zanol (hard target): ratchet
  esc_med -8 @ 6.0 Mcand = ~-1.33 steps/Mcand; best sectorial -7 @ 27.9 Mcand = ~-0.25
  steps/Mcand. Ratchet is ~5x more candidate-efficient at escaping — the same ordering
  the (noisier) sub-second wall ratios show, on the axis that doesn't drift with hardware.

## Verdict: diagnostic clue RESOLVED; wall-lever RULED OUT; thread retires
Two separate conclusions, both clean:
1. **Mechanism (the diagnostic obligation): RESOLVED.** Ours now escapes to the TNT
   sectsch targets (best-of-15, all 4 datasets) where it historically escaped 0. The
   wagner-fix + rasStarts>=3 rebuild is the barrier-crosser. No missing capability remains.
2. **Wall-clock lever: RULED OUT.** The now-working escape is Pareto-dominated by ratchet
   (~5x less candidate-efficient on the hard Zanol target; loses on 3/4, needs 2x wall on
   the 4th), so it cannot displace ratchet — not a time-to-optimum win.

The walk-up-selection src routine is NOT justified: the config proxy already REACHES the
escape targets, so the limiter is wall/candidate-efficiency, not sampling geometry —
walk-up would only matter if far more candidate-efficient, and ratchet already crushes
sectorial on escape-per-candidate. Leave it un-written.

HANDOFF (one line, goofy-cannon's call): rasStarts>=3 sectorial is lower-variance than a
SINGLE ratchet start on ratchet-unreliable datasets (Zhu median -6 vs ratchet -1) — a
single-vs-multi-start BUDGET/recipe question, not a distinct geometry lever. Worth a
thorough-recipe slot iff the recipe work wants it; not pursued here.

Plumbing kept: `rssPicks` SearchControl knob (exposes existing SectorParams field,
default 0 = byte-identical) — reusable for any future sequential-picks experiment.
