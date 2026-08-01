# Gating NA certification: `TBRParams::certify_unrooted`

Follows `na-exact-verify-dominates.md`, which measured `exact_verify_sweep` at
**97.7% of `tbr_search` wall** on inapplicable data and proposed a
`certify_unrooted` flag as lever (a). This note is the mechanism, two corrections
to the premise, and the gate.

Status **2026-07-29: mechanism landed opt-in; panel not yet returned.** The
default is unchanged and byte-identical to `0a9e7c24` (verified: 3 matrices x
2 seeds x tabuSize {100, 0}, identical scores *and* identical topology sets).

## Mechanism

`TBRParams::certify_unrooted` (default `true`) says whether *this* caller needs a
certified optimum. Cleared by every internal sub-search whose tree is never
reported: `ts_ratchet.cpp` (both cycle searches and the baseline), `ts_drift.cpp`,
`ts_sector.cpp` (6 sites), `ts_fuse.cpp`, `ts_nni_perturb.cpp`,
`ts_prune_reinsert.cpp`, and the anneal reconverge in `ts_driven.cpp`. Left at
the default for `driven_search`'s step-2 and step-6 whole-tree polishes and for
the diagnostic Rcpp entry points.

Honoured only when **`TS_NA_NOCERTIFY`** is set. That is deliberate: skipping
certification returns worse trees, so the default flips only if the panel below
earns it, and until then the merge is provably a no-op.

On the skip path `best_score` is re-synced with an explicit `full_rescore`.
`exact_verify_sweep` did that on both its exits, and a score that drifts from its
returned topology (the `ts_tbr.cpp:660` bug class) would read as a quality
regression in a floor-attainment panel when it is really a reporting bug.
Pinned by `tests/testthat/test-ts-na-certify.R`.

## Correction 1: under the shipped presets the certifier is a SECTOR cost

`do_reroot` (`ts_tbr.cpp:1843`) requires `params.tabu_size == 0`. The shipped
presets set `tabuSize = 100` (`default`) and `200` (`thorough`); only `sprint`
sets 0. The 97.7% was measured through `ts_tbr_search` / `ts_ratchet_search`,
whose `TBRParams` / `RatchetParams` leave `tabu_size` at 0 — **not** the shipped
recipe.

Attributed by phase-ablation (`dev/profiling/drivers/na_certify_probe.R`,
`MaximizeParsimony(strategy = "default")`, 2 replicates, `naDiag$n_evs`):

| Vinther2008 (23t) | wall | sweeps | evs ms |
|---|---|---|---|
| full default | 0.34 s | 4 | 62 |
| `xssRounds/rssRounds/cssRounds = 0` | 0.17 s | **0** | 0 |
| `fuseInterval = 0` | 0.29 s | 4 | 60 |
| `ratchetCycles = 0` | 0.25 s | 4 | 52 |
| `tabuSize = 0` | 1.01 s | **40** | 673 |

Removing the sector phase takes the sweep count to zero; removing fuse or the
ratchet does not touch it. So under `default`:

* the sector sub-searches (and, when a fuse improves, the fuse cleanup) are the
  only call sites that reach the certifier — they leave `tabu_size` at 0;
* the ratchet, drift and both driven polishes never reach it;
* therefore **the tree a replicate reports has never been certified**, while
  trees nobody reports are certified repeatedly. The spend is in the wrong place,
  which is a finding independent of whether gating pays.

`tabuSize` is consequently an arm dimension of the panel, not a detail: at 0 the
sweep fires at every whole-tree convergence (40 calls vs 4).

## Correction 2: certification is not a formality, it is the only exact step

Under Brazeau's three-pass the indirect scan is approximate, so the sweep does
not merely *prove* optimality — it *finds* improvers the scan cannot see.
`naDiag$n_evs_improved` counts those, and it is routinely non-zero.

Zanol2014 (74t, 64% NA), `default` preset, 2 replicates, seed 1:

| arm | score | wall | sweeps |
|---|---|---|---|
| certify (today) | **1321** | 43.3 s | 20 (16.3 s) |
| `TS_NA_NOCERTIFY=1` | **1326** | 2.74 s | 0 (14 skipped) |
| no sector at all | 1325 | 0.83 s | 0 |

Five steps for 40 s. So the lever is a genuine quality-for-speed trade and
**must** be judged on floor attainment, never on wall alone — and never on any
metric whose numerator is improvements found (`reach-not-escape-rate`).

## Third arm: certify the reported result once

Given correction 1, the coherent change is not only "certify less" but "certify
in the right place". **`TS_NA_FINAL_CERTIFY`** (also default off) runs one
certification per replicate, on the tree that replicate contributes to the pool
(`run_single_replicate`, after the outer-cycle loop; forces `tabu_size = 0`
because `do_reroot` gates on it; skipped under constraints, where `do_reroot` is
off anyway).

This is also the answer to the delegated question *"should a cleared flag still
certify when the tree is a new global best?"* — **no**. `tbr_search` has no view
of the caller's incumbent, so that would mean threading a threshold through
`TBRParams` for a case one end-of-replicate pass covers better. The related
question *"should the ratchet's final cycle certify?"* is **no**, and is
mis-posed: `ratchet_search` returns `best_tree` across cycles and discards a
non-improving cycle's tree outright, so certifying the last cycle certifies a
tree that may be thrown away. `driven_search`'s last convergence keeps the
default (it is the reported one, and `n_outer` is typically 1).

## The gate

`dev/benchmarks/bench_na_certify_cell.R` (cell = dataset x seed x tabuSize; all
arms on one node), `hamilton_na_certify_{build,array}.sh`,
`na_certify_analyze.R`.

* **Panel**: all **30** bundled `inapplicable.phyData` matrices, **native** (never
  recoded `-`→`?`), 20–88 tips. Five seeds. `strategy = "default"`, EW
  (`concavity = Inf`, the `MaximizeParsimony` default; `TS_CONCAVITY=10` gives the
  IW regime the profile used). `targetHits`/`stopPatience` neutralised so the
  budget alone ends a run. 30 x 5 x 2 = 300 array tasks.
  Four matrices — what the profile used — cannot support a paired test: a sign
  test at n = 4 bottoms out at p = 0.125.
* **Arms**: `A_certify` (today) / `B_gate` / `C_final` at a fixed 8-replicate
  budget, then `B_gate_mw` / `C_final_mw` at `maxSeconds = ` arm A's *own*
  measured wall, re-run in the same cell. The matched-wall arms exist because a
  wall win is spendable on replicates, and a fixed-replicate comparison cannot
  see that.
* **Metrics**: floor attainment (share of seeds reaching the best score any arm
  found for that matrix) aggregated to **one number per matrix** before any
  paired test — `(matrix x seed)` pairing is pseudo-replication and has
  manufactured p = 0.0007 from an effect that was p = 0.98 at matrix level. Wall
  as per-matrix **median** + sign count + count >10% slower, never the mean.
  Plus `n_evs` / `n_evs_skipped` per arm: **a null result with
  `n_evs_skipped == 0` means the flag never fired**, which given correction 1 is
  the likelier reading, and is a different finding from "it fired and bought
  nothing".
* **Verdict rule**: an arm ships only if floor attainment does not regress.
* Training data only. The MorphoBank validation split is sealed and informs
  nothing here; the bundled inapplicable corpus is not part of that split.

### Local pre-read (2 matrices x 2 seeds, `tabuSize = 100`, 4 replicates)

Suggestive only — n = 2 matrices, and the sample is deliberately mid-size.

| cell | A_certify | B_gate | C_final | B_gate_mw | C_final_mw |
|---|---|---|---|---|---|
| Wills2012 s1 | 273 / 14.4 s | 273 / 0.8 s | 273 / 5.4 s | 273 (111 reps) | 273 (12 reps) |
| Wills2012 s2 | 273 / 23.1 s | 273 / 0.7 s | 273 / 3.5 s | 273 (181 reps) | 273 (22 reps) |
| Griswold1999 s1 | **407** / 8.6 s | 408 / 0.9 s | 409 / 2.0 s | **407** (64 reps) | **407** (17 reps) |
| Griswold1999 s2 | 407 / 10.5 s | 407 / 0.6 s | 407 / 1.9 s | 407 (100 reps) | 407 (22 reps) |

The shape to test at scale: at a fixed replicate budget the gate costs reach
(Griswold1999 s1: 407 → 408/409); at matched wall it buys 16–45x the replicates
and recovers it. That is the `lever-b-oversearch` pattern — certification is
over-search — but two matrices cannot establish it, and the large NA matrices
(Zanol2014, Zhu2013, Dikow2009, Giles2015) are exactly where the fixed-replicate
gap was largest and are absent from this pre-read.

Note `C_final` is not a strict improvement on `B_gate` at fixed replicates
(Griswold1999 s1: 409 vs 408): the final certification moves the tree, which
changes the pool and every later replicate's trajectory. It is a different
trajectory, not a superset.

### Reading the results — four things not to get wrong

1. **The pre-read's shape is least likely to hold where the lever matters.** These
   are 43–55 tips; the fixed-replicate gap was 5 steps on Zanol2014 (74t) against
   1–2 here. Certification cost is O(n³) while the replicates a wall win buys
   scale inversely, so on 74–88 tips arm B gets *fewer* extra replicates and has
   *more* reach to recover. Read the large-matrix rows first. If `B_gate_mw` fails
   to recover on Zanol2014 / Zhu2013 / Dikow2009 / Giles2015, that is the verdict
   whatever the 20–40-tip matrices say, and the honest conclusion is a
   size-conditioned default, not a global flip.
2. **`tabuSize = 0` is a different search, not "default plus certification".** It
   also disables the tabu list, changing plateau exploration for *every* arm.
   That is fine for the A-vs-B contrast within a tabu level, which is why the
   analyser segments by it — but the `tabu = 0` rows bound the lever's size, they
   are not evidence about `default`-preset behaviour.
3. **`nEvsImproved` is not "improvers the gate cost you".** For gated arms it
   counts only the final-certify improvers (arm C) or driven-polish improvers
   (arm B at `tabu = 0`). What certification *would* have found in a run that
   never ran it is unobservable by construction.
4. **Check the matched-wall budget audit before the verdict.** `maxSeconds` is
   polled at replicate boundaries, so a `*_mw` arm can overshoot; an arm handed
   more wall than arm A is not comparable to it. The analyser prints per-matrix
   overshoot and warns above 1.25x. The `run_single_replicate` budget-spent guard
   that prevents a certification starting past the deadline landed *after* the
   panel array was built, so the panel's `C_final_mw` rows are the unguarded
   behaviour — the audit is how that gets caught rather than assumed away.

The panel build is the pre-rebase tree (`b9bc14d6`); the mechanism is unchanged by
the rebase onto `12a5866d`, which pulled in unrelated T-373/T-378 commits.

## PANEL RESULT (2026-07-30, array 18080528, 300/300 COMPLETED)

30 native matrices x 5 seeds x tabuSize {100, 0}; unit of replication = MATRIX
(n = 30); wall = per-matrix median; two-sided exact sign test.

**Budget audit passed, in the conservative direction.** Every matched-wall arm
came in at **0.90-0.97x** arm A's wall — an *under*-spend, because `maxSeconds`
is polled at replicate boundaries so a replicate that would overrun never starts.
The `*_mw` arms were therefore given ~9% **less** wall than the arm they beat, so
their wins are floors, not artefacts. No overshoot anywhere.

**The gate fired everywhere it was supposed to.** At `tabuSize = 100`, `B_gate`
executed **0** sweeps against 5 526 skipped — under the shipped preset the gate
removes *all* certification. `A_certify` ran 6 957 sweeps of which **1 087 (16%)
found a real improver** (`nEvsImproved`), which is why removing them costs reach.

### tabuSize = 100 — the shipped `default` preset

| arm | attain Δ | better / worse | p | wall (median ratio) |
|---|---|---|---|---|
| `B_gate` | −0.100 | 3 / 12 | **0.035** | 0.052 (**19x faster**) |
| `C_final` | −0.033 | 5 / 8 | 0.581 | 0.294 (3.4x) |
| `B_gate_mw` | **+0.167** | **10 / 1** | **0.012** | 0.904 |
| `C_final_mw` | +0.087 | 8 / 1 | 0.039 | 0.923 |

### tabuSize = 0 — `sprint`, and the configuration the 97.7% was profiled in

| arm | attain Δ | better / worse | p | wall (median ratio) |
|---|---|---|---|---|
| `B_gate` | −0.140 | 0 / 12 | **0.00049** | 0.160 |
| `C_final` | −0.107 | 0 / 10 | **0.0020** | 0.165 |
| `B_gate_mw` | +0.013 | 4 / 2 | 0.687 | 0.906 |
| `C_final_mw` | +0.033 | **5 / 0** | 0.063 | 0.907 |

### Verdict: certification is OVER-SEARCH — do not flip the default

**At equal replicates the gate regresses reach at both tabu levels** — but the
two levels differ in kind, and the distinction matters. At `tabu = 0` the
regression is **uniform**: 0 matrices better, 12 worse. At `tabu = 100` it is a
**net, not a direction**: 3 matrices *improved* with certification removed
(Aguado2009, Aria2015, Conrad2008, all +0.2) against 12 that lost. So under the
shipped preset the gate already trades in both directions at fixed replicates,
and the sign test is reporting a balance rather than a one-way mechanism. 12 > 3
at p = 0.035, so the rule (`ships only if floor attainment does not regress`)
still says `B_gate` does not ship as a default — but do not narrate it as
"removing improvers costs reach everywhere", because at `tabu = 100` that is not
what the data show.

`C_final` is not significant either way, but its point estimate is also negative,
so it is "not shown to regress" rather than "shown not to regress" — not enough
for a regression-averse default (`auto-vs-thorough-objective`).

**At equal wall the same gate wins decisively.** Under the production preset
`B_gate_mw` improves floor attainment on 10 of 30 matrices and loses on 1, while
spending 10% less wall than the arm it beats. The gains land exactly where the
engine is weakest:

| matrix | A_certify | B_gate_mw |
|---|---|---|
| Aguado2009 | 0.0 | **1.0** |
| Geisler2001 | 0.2 | **1.0** |
| Zhu2013 | 0.2 | **1.0** |
| Dikow2009 | 0.4 | **1.0** |
| Liljeblad2008 | 0.4 | **1.0** |
| Aria2015 | 0.6 | **1.0** |
| Wortley2006 | 0.2 | 0.6 |

So the finding is not "certification is wasteful" — it demonstrably finds
improvers 16% of the time. It is that **certification is over-search**: the wall
it consumes buys more reach when spent on replicates instead. That is
`lever-b-oversearch` reappearing on the NA path.

**Which regime production is in decides the default.** `MaximizeParsimony`
defaults to `maxSeconds = 0` with `maxReplicates = 96` and `targetHits` stopping
— a *replicate*-bounded budget, i.e. the regime where the gate loses. So the
opt-in default stays.

**And "set `maxSeconds` and gate" is NOT the arm that won.** `B_gate_mw` ran with
`maxReplicates = 1000`; at arm A's wall it completed a median of **101–226
replicates on 23 of the 30 matrices** (Zhu2013 226, Wills2012 222, Zanol2014 173,
Dikow2009 147). Today's default cap of 96 would bind on all 23, leaving the freed
wall unspent, so a user who merely sets `maxSeconds` gets a replicate-capped run
rather than the arm measured here. The recommendation is therefore
**`TS_NA_NOCERTIFY=1` *together with* a raised `maxReplicates`**, for wall-bounded
runs. Making that the default is a recipe change belonging to `campaign-recipes`,
gated on its own panel — not a flag flip here.

`C_final_mw` is the exception: it never exceeded **33** replicates on any matrix,
so the 96 cap never binds for it. It buys less (+0.087 vs +0.167) but it works
with the shipped budget unchanged.

**The one honest counter-example is Zanol2014**, the hardest matrix in the corpus:
`B_gate_mw` is −0.2 at both tabu levels, i.e. certification pays there even at
matched wall. `C_final_mw` recovers it (+0.2 at `tabu = 0`, level at 100) and is
the only arm in the panel with **no regression observed anywhere** (5 better /
0 worse at `tabu = 0` — p = 0.063 on 5 discordant matrices, so suggestive rather
than established). If a size- or difficulty-conditioned default is ever built,
`C_final` is the arm for the hard tail.

*Scale caveat:* the per-block floor is "best any arm found **in that block**", and
Zanol2014's differs between them — 1311 at `tabu = 0`, 1312 at `tabu = 100`. Its
0.2 in the two blocks is therefore not the same achievement, and the blocks'
attainment numbers for that matrix are not directly comparable to each other.

## PANEL 2 RESULT — the shipped stopping rules (2026-08-01, array 18126933)

150/150 cells (30 matrices x 5 seeds x 6 arms), all at pin `644b5b10`. Every
stopping rule left exactly as shipped: `maxSeconds = 0`, `maxReplicates = 96`,
`targetHits = max(10, ntax/5)`. Raw cells in
`dev/profiling/na-certify-stop-hamilton.csv`.

| arm | attain Δ | better / worse | p | wall (median ratio) |
|---|---|---|---|---|
| `B_gate` | −0.047 | 0 / 4 | 0.125 | 0.060 (**17x faster**) |
| `C_final` | −0.013 | 1 / 2 | **1.000** | 0.322 (3.1x) |
| `B_gate_reps` | −0.033 | 1 / 3 | 0.625 | 0.074 |
| **`A_hits3`** | **0.000** | **0 / 0** | — | **2.582 (2.6x SLOWER)** |
| `B_gate_hits3` | −0.040 | 0 / 4 | 0.125 | 0.128 |

### `targetHits` escalation is pure cost — 0 of 30 matrices improved

Tripling the hit target changed floor attainment on **not one matrix**, while
costing **2.58x the wall** (26 of 30 slower, 22 of them by >10%, p = 6e-05).

It is not that the runs ignored it: only 33 of 150 cells were identical in score
*and* replicate count, so runs genuinely went longer — they just never found
anything better. The mechanism is the disjoint-population argument, now measured
at corpus scale rather than on one cell:

* on **hard** matrices the replicate cap binds before the hit target is reached,
  so raising it changes nothing (`hitCapBound` rises 22% → 47% between
  `A_default` and `A_hits3` — the extra demand just pushes more runs into the
  cap);
* on **easy** matrices it does add replicates, but those matrices already attain
  1.0, so there is nothing left to find.

The extra work therefore lands exactly where it cannot help. **`targetHits` is
not an effort knob for reach on inapplicable data.**

### Gating still costs reach here, so the default still does not flip

`B_gate` is 0 better / 4 worse — not significant at n = 30 (p = 0.125) but
uniformly negative, and the losses are concentrated on the hard tail:
Zanol2014 −0.6, Wortley2006 −0.4, Aria2015 −0.2, Zhu2013 −0.2. Panel 1's verdict
survives contact with the real stopping rules.

### `C_final` is the arm worth pursuing

Statistically indistinguishable from baseline (−0.013, 1 better / 2 worse,
p = 1.000) at **a third of the wall**, and it *gains* on Aguado2009 (+0.4). Its
certifications are also the productive ones: 11 000 sweeps of which **4 349
(40%) found a real improver**, against `A_default`'s 9 555 of 44 889 (21%) —
certifying the tree a replicate actually reports is four times cheaper and twice
as likely to pay.

Wall on the hard matrices, median seconds:

| matrix | `A_default` | `C_final` | `B_gate` |
|---|---|---|---|
| Zanol2014 | 1947 | 628 | 86 |
| Dikow2009 | 1437 | 365 | 69 |
| Giles2015 | 1294 | 351 | 46 |
| Zhu2013 | 1256 | 330 | 37 |
| Aguado2009 | 1008 | 245 | 39 |

**But it regresses on Zanol2014 (0.6 → 0.0)**, the hardest matrix in the corpus
and the same one panel 1 flagged. By the stated rule that is a regression, so
`C_final` is a candidate for a future default and not a change to make now. What
would settle it is the hard tail specifically — Zanol2014, Zhu2013, Wortley2006,
Aguado2009 at more seeds — rather than another corpus-wide sweep, since 24 of the
30 matrices are saturated at 1.0 in every arm and can only dilute the signal.

## What to switch on, in practical terms

**Presets today: nothing changes.** Certification stays on in `default` and
`thorough`. Neither panel-1 budget is the shipped one, and a preset default may
not be flipped on a regime that was not measured.

**`thorough`: keep certification on, permanently.** Not "pending panel 2" — its
objective is reaching the global optimum, not wall (`auto-vs-thorough-objective`),
and Zanol2014, the hard tail, is precisely where certification still pays at
matched wall.

**`default` / `auto`: the candidate is a PAIR, not a flag.** Gating alone is the
arm that *loses*. The shipped `maxReplicates = 96` is what stops the freed wall
being spendable, so any flip must be `TS_NA_NOCERTIFY` **plus** a raised replicate
cap. Panel 2 (`bench_na_certify_stop_cell.R`, array 18096945) decides it under the
real stopping rules.

For a user on inapplicable data **today**:

```r
# Wall-limited, want the best tree within a fixed time: the +0.167 arm.
Sys.setenv(TS_NA_NOCERTIFY = "1")
MaximizeParsimony(dat, maxSeconds = 600, maxReplicates = 500)
```

Want the best tree regardless of time: change nothing. **Do not set
`TS_NA_NOCERTIFY=1` on its own** — without the raised cap that is the arm that
regressed.

## Does anything need adjusting when `targetHits` rises?

**No compensating adjustment, and raising it helps rather than hurts** — but it is
not the knob that pays for gating.

`targetHits` defaults to `max(10, ntax/5)` and stops the run once the best score
has been hit that many times. Raising it means more replicates, which moves the
run toward the many-replicate regime where the gate won. It is already the shipped
idiom for "this matrix is hard": `.IwRatchetDepth()` reads `targetHits /
defaultHits` as a user escalation signal and deepens the IW ratchet by it.

Two cautions:

* **`targetHits` counts hits against the CURRENT best, not the true optimum.** So
  it cannot rescue a uniformly weaker search that plateaus one step high: raising
  it just buys more confirmations of the same wrong score. What protects against
  that is more *independent* replicates. **`maxReplicates` is the knob that pays
  for gating; `targetHits` is not.**
* **Raising `targetHits` without raising `maxReplicates` makes the run bump the
  96 cap instead of reaching its hit target** — and the more expensive the
  replicate, the sooner that happens. This is the concrete reason the two knobs
  are coupled on NA data: certification is what makes the replicate budget
  unusable.

Panel 2's `A_hits3` / `B_gate_hits3` arms measure exactly this at 3x the default
hit target.

### Harvesting the panel

Build `18080527` (COMPLETED, gate smoke `n_evs_skipped = 8`); array `18080528`
(300 tasks), submitted 2026-07-29.

```bash
ssh hamilton8.dur.ac.uk 'ls /nobackup/$USER/TreeSearch/na_certify_partials | wc -l'
```

When complete (300 files), pull and analyse:

```bash
scp -r hamilton8.dur.ac.uk:/nobackup/pjjg18/TreeSearch/na_certify_partials dev/benchmarks/
Rscript dev/benchmarks/na_certify_analyze.R dev/benchmarks/na_certify_partials
```

Read in this order: (1) the matched-wall budget audit — a warning there
invalidates the `*_mw` rows; (2) `nEvsSkipped` per arm — zero means the flag never
fired and the arm says nothing; (3) the large-matrix rows of the floor-attainment
table; (4) only then the paired summary.

## What this does not do

Pruning the sweep remains the larger, separate project.
`dev/red-team/union-of-finals-bound-proof.md` shows the deployed union-of-finals
is a sound *upper* bound (the in-code comment is backwards), so it cannot prune;
its step S3 shows the union of *true MPR* finals would be a valid lower bound,
and the engine does not compute those. Gating is the tractable half.
