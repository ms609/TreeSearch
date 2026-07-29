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
