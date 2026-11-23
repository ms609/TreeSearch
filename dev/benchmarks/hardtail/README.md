# project5432 hard-tail — best-known trees & the shared-floor result

Two score-verified trees for the Mission-A epitome instance **project5432**
(482 taxa, 189 characters, ~44% missing; `.../lgsweep/matrices/project5432.nex`,
a.k.a. `m01`). Committed here so the best-known artifacts can never go missing
again (the lesson of the retracted "1939", which was a sectorial-search *trace
phantom*, never a real tree — see below).

| file | score | reached by | verified |
|---|---|---|---|
| `project5432_best_1942_tnt.tre` | **1942** | TNT regeneration (`xmult`, fine blocks) — one lucky draw | canonical re-score = 1942 ✓ (2026-07-22) |
| `project5432_ts_reach_1943.tre` | **1943** | TreeSearch, aggressive perturbation swarm (below), seed s15 | canonical re-score = 1943 ✓ (2026-07-22) |

## Scoring regime (get this wrong and every number is noise)

The canonical objective for project5432 is **unordered Fitch, gaps-as-missing**:

```r
TreeLength(Preorder(t), pd, concavity = Inf, inapplicable = "missing")
```

This equals TNT `ccode -` (unordered) **exactly** — the floor tree scores 1943
under both. TNT `ccode +` (ordered/additive) scores 2475, a different regime;
do not compare across it. Both trees here carry **taxon-name** tip labels (not
the 0-indexed integer labels TNT floor trees ship with — a recurring
alignment gotcha; relabel to matrix order before scoring integer-labelled
trees). Both are binary (score with `collapse = FALSE` upstream, else
`TreeLength` rejects the collapsed MPT).

## The headline: 1943 is a *shared* TS/TNT floor

The best-known ladder for project5432 is **1942 (stored) / 1943 (both engines'
working floor) / 1939 (does not exist — retracted phantom)**:

- **1942** is a *rare* draw. A 28 h heavy TNT hammer (16 seeds, ratchet
  50–120 + drift, `hits 999`) floored at **1943 across every seed**; the 1942
  tree came from a single lucky regeneration run, not a repeatable recipe.
- **1943** is now reached by *both* engines' working search. TNT reaches it
  reliably; TreeSearch reaches it too (this directory's `_ts_reach_1943` tree).
- **1939** was never a real tree — it was a sectorial-search *bestScore-column
  trace phantom* (a sector estimate against a frozen backbone), unreproduced by
  every method. `best` over all 10 000 saved trees of the run that "showed"
  1939 = 1942.

**⇒ The residual 1943 → 1942 gap is a rarity / luck-and-speed gap, not a
TreeSearch pathology.** Both engines sit in the same 1943-quality basin regime;
1942 is a rare escape for *either*. This corrects the older framing that read
the gap as a TS-specific mechanism failure ("TS 0/20 cold vs TNT 4/12 → 1943").

## How the TS 1943 was reached (recorded observationally)

Hamilton array **17943743** (`ts_warm.sh`), 16 seeds, ~48 h to the 47:59 wall.
A block loop: block 1 is a **cold** `MaximizeParsimony(strategy="thorough")`;
each later block carries the running best forward (`tree = best`) and re-runs
the same heavy config, saving on every improvement (so the wall never costs
progress). ~57 blocks/seed. Config (deltas from `thorough`):

```
ratchetCycles=40  ratchetPerturbMaxMoves=0 (auto/deep)  ratchetAdaptive=TRUE
driftCycles=25  postRatchetSectorial=TRUE  stallEscalateFactor=1.5
adaptiveStart=TRUE  intraFuse=TRUE  poolSuboptimal=3
xss/rss/css 5/3/2  sectorMax=80  rasStarts=3  outerCycles=2  maxOuterResets=3
```

Result: **floor 1943, 6 of 16 seeds** (s5, s7, s8, s11, s15, s16); remaining
seeds at 1944. What is and isn't claimed:

- **Block-1 cold still floors ~1945–1947.** It is the *sustained iteration*
  under this stack that descends to 1943 — the cold single-start floor is
  unchanged, consistent with the campaign's cold-start diagnosis.
- **The differentiator is the perturbation *stack*, not restart volume.** An
  earlier 1870-replicate breadth run scored 0 hits below the floor; raw restart
  count is not what moves it. Because `tree=` warms only replicate 1 (7/8 reps
  per block are cold RAS), this run cannot separate warm-seed accumulation from
  repeated perturbed restarts — so no single mechanism is named beyond "the
  aggressive perturbation stack, iterated."
- **1943 is the new best TS-*generated* score** for project5432, one step below
  the prior TS-side best of 1944 (`floors/project5432_ts_1944.tre`).

## Footnote — the "should we re-run?" decision

This swarm ran on the Hamilton lib **TreeSearch 2.0.0** (built 2026-07-08),
which lags `cpp-search`: it lacks the sector-drift levers
(`sectorGoDrift`/`sectorDriftCycles`, stripped from the config), and predates
merged engine work (exact SPR #273, L3b lever #6, packing, cutoff). So this
records *"TS 2.0.0 + this stack → 1943"*, not *"the TS engine cannot reach
1942."* Settling that would mean rebuilding current `cpp-search` on Hamilton
and re-running the same config with the sector-drift levers restored (~2-day
16-seed array). Held for a go decision — not launched.

Provenance & sibling instances (where the merged engine already reaches the
TNT floor cold): campaign hub in memory `campaign-mission-a-5432`, and
`dev/plans/2026-07-14-mission-A-cold-start-basin-capture.md`.
