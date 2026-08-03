# Does more effort recover a better score on tough matrices?

Panel 3 (array 18143234, 60/60 cells at pin `db5c6578`). Six hard
inapplicable matrices x 10 seeds x {`effort` 0/+1/+2} x {certify, gate}.
Scored on **score**, not binary attainment. Raw cells:
`dev/profiling/na-hardtail-hamilton.csv`.

Panels 1 and 2 never tested this. Panel 2's `targetHits` arm could not act where
it mattered (see the retraction in `na-certify-gate.md`), and a first draft of
this panel repeated the mistake in mirror image — raising only `maxReplicates`
was inert on Aria2015, where `targetHits` binds at 39 replicates. Only moving
both escalates every dataset, so the arms are the shipped `effort` dial.

## Yes — but only where there is headroom, and only for the first notch

Floor attainment, share of 10 seeds reaching that matrix's best known score:

| matrix | certify e0 | e1 | e2 | gate e0 | e1 | e2 |
|---|---|---|---|---|---|---|
| Zanol2014 | **0.5** | **0.9** | 0.9 | 0.0 | 0.3 | 0.3 |
| Aria2015 | 0.9 | 1.0 | 1.0 | 0.7 | 1.0 | 1.0 |
| Wortley2006 | 1.0 | 1.0 | 1.0 | 0.8 | 0.9 | 0.9 |
| Aguado2009 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| Geisler2001 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| Zhu2013 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |

Zanol2014 is the case the hypothesis describes: at `effort = 0` half the seeds
stop at 1312, at `+1` nine of ten reach 1311. Per seed, four runs (3, 6, 9, 10)
move 1312 → 1311 and none move the other way.

Direction never reverses anywhere: **5 of 60 certified cells improved from e0 to
e2 and 0 worsened; 10 of 60 gated cells improved and 1 worsened.** But the effect
lives almost entirely in one matrix, so treat this as descriptive. At matrix
level it is 2 of 6 better, 0 worse — a sign test there is n = 2 and says nothing.
That thinness was predicted and is the honest limit of the bundled NA corpus:
under the real shipped recipe, only Zanol2014 is still hard.

## The second notch is inert

`effort = +2` matched `+1` on **every matrix, both certification states**, while
doing measurably more work — Zanol2014 ran 798 replicates against 500, for
14 493 s against 8 747 s. Returns die immediately after the first notch.

Aria2015 also confirms the **dead notch** predicted from the smoke test: `+1` and
`+2` both ran exactly 24 replicates, because it is `targetHits`-bound and rungs
1–4 never raise the hit target. The notch is not merely low-yield there, it is
structurally incapable of acting.

## Certification is worth more than 10x the replicates

| Zanol2014 arm | attainment | wall (median) |
|---|---|---|
| `certify_e0` | 0.5 | 1 795 s |
| `certify_e1` | **0.9** | 8 747 s |
| `gate_e2` | 0.3 | 8 235 s |

`certify_e0` beats `gate_e2` on reach at **a fifth of the wall**, and `gate` never
exceeds 0.3 on this matrix at any effort. **This closes the gating question on the
hard tail: do not gate.** Panel 1's matched-wall win was real but does not survive
here — the wall gating frees cannot be spent to buy back what certification finds,
even at ten times the replicates.

## The actionable finding is about `.AutoRung`, not about `effort`

`effort = +1` on Zanol2014 means rung 4 (`large`) — which `auto` reserves for
>= 120 tips. Zanol2014 has **74**. So the automatic choice is one notch too low on
this matrix, and the reason is that `.AutoRung` reads size and character count
only, while difficulty has a second axis it cannot see (`two-axis-difficulty`:
score-reach vs island-completeness, predicted by homoplasy as much as size).

That is a `campaign-recipes` question — a difficulty predictor for the rung — and
it is worth more than any further tuning of the ladder above rung 4, which this
panel shows is flat.

## Scoping note that matters for reading panel 2

Panel 2 pinned `strategy = "default"` on every matrix, so on the 65–119-tip
matrices it measured a preset weaker than the shipped `auto` (which selects
`thorough` there). Panel 3's `effort = 0` is the genuine default, and under it
five of these six matrices saturate — which is why panel 2 showed headroom on
Aguado2009 and Geisler2001 that has since closed.
