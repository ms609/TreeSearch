# A/B verdict — bandit credit for user-warm-started replicates

**Date:** 2026-07-27 · **Harness:** `dev/benchmarks/ab_bandit_warmstart.R` ·
**Raw:** `dev/benchmarks/ab_bandit_warmstart.csv` (local only —
`dev/benchmarks/*.csv` is gitignored, so the tables below are the durable record)

The harness was run twice by accident; the CSV holds **replication 2 only**
(the second run overwrote it).  Scores were identical between replications, so
nothing is lost there; replication 1's wall figures survive only in the quoted
table below.

## The change

`adaptiveStart = TRUE` updated the Thompson-sampling bandit for *every*
replicate, including those handed a starting tree via `tree = `.  Such a
replicate never reaches the strategy switch in `run_single_replicate()`
(`ts_driven.cpp`, `if (starting_tree)` short-circuits it), so the default
`WAGNER_RANDOM` was credited for a Wagner start that was never built.
Reseeded (`POOL_RESEED`) reps were already excluded on exactly this reasoning.

Latent for replicate 0 alone until `tree = <multiPhylo>` began warm-starting one
replicate per supplied tree — a k-tree pool then feeds k phantom credits.

`decay()` was hoisted out of the same guard: it measures how stale the
accumulated evidence is, not which arm ran, so it still fires when a
warm-started rep improves the best score.  Leaving it inside would have made
this A/B measure two changes at once.

## Design

Only `adaptiveStart = TRUE` **and** `tree =` together can differ, and only when
spare **cold** replicates remain after the pool is consumed — those are the reps
whose arm is drawn from the (previously polluted) posterior.  With
`nRep == nPool` the fix is provably invisible.  Hence `nPool = 5`, `nRep = 20`.

3 datasets × 5 seeds, serial (`nThreads = 1`).  One starting pool per
(dataset, seed), generated once with the bandit **off** and shared verbatim by
both arms, so no phase-1 noise enters the comparison.  Both builds were compiled
from the same working tree, differing only in the guard (the `before` library
was built by temporarily reverting it in place — *not* from HEAD, which lacks
the multiPhylo feature entirely and would make the comparison meaningless).

**Pre-registered decision rule** (fixed before results were seen): adopt unless
time-to-optimum regresses.  A tie adopts, because the phantom credit is wrong on
its own terms.

## Results

Scores are deterministic and were identical across two independent
replications; wall-clock differs only by timing noise.

| Measure | Outcome |
| --- | --- |
| Score (paired, n = 15) | 13 tied, 1 after-better, 1 after-worse; mean −0.067 |
| Wall, replication 1 | after faster 9/15, mean **−1.48 s** |
| Wall, replication 2 | after faster 11/15, mean **−1.95 s** |

**Mechanism (the high-signal readout).**  Share of bandit attempts going to
`wag_rand`, summed over all 15 runs:

| Build | `wag_rand` | `wag_golob` | `wag_entropy` | `rand_tree` | total |
| --- | --- | --- | --- | --- | --- |
| before | **181 (60.3 %)** | 46 | 35 | 38 | 300 |
| after | **62 (25.9 %)** | 47 | 70 | 60 | 239 |

Uniform would be 25 %.  The phantom credits do not merely add k counts — they
compound: an early unearned success raises `wag_rand`'s posterior, which wins
more Thompson draws, which earns more real credit.  The clearest case is
`Agnarsson2004` seed 2, where `before` recorded `20,0,0,0` — the bandit locked
onto `wag_rand` for all 15 cold replicates and never sampled another arm.  After
the fix the same cell explores (`7,0,9,0`).

`before` totals 300 = 15 runs × 20 reps (every rep votes); `after` totals
239 ≈ 15 × 15 cold reps (the shortfall is reps that ended early).

## Verdict

**ADOPT.**  No time-to-optimum regression detected: wall trends slightly in
`after`'s favour in both replications, and score is a 13/15 tie with a 1–1
split.  The wall measurements were taken on a machine shared with another
active session (see scope note), so they support "no regression" but should
**not** be read as a speed-up.  The mechanism result
is unambiguous and is the real evidence: arm selection is no longer dominated by
credit for work that never happened.

## Power and scope — read before citing this

- **Score is underpowered here.**  `Agnarsson2004` and `Griswold1999` returned a
  single score in every cell (778, 407), so they carry *no* discriminating
  power on reach; only `Wortley2006` varied (482–484).  The score column is
  evidence of *no regression*, not of equivalence.
- **The wall column was measured under CPU contention.**  Another session was
  editing and rebuilding this checkout across the same window.  Score and the
  attempts distribution are deterministic and immune; wall is not.  Contention
  hits both arms (they are interleaved per seed) and the direction held across
  two replications, so the "no regression" reading stands — but a clean
  machine is needed before quoting any speed-up.
- Serial only; the parallel path uses a stateless round-robin and has no bandit
  (see below), but was not benchmarked.
- 3 datasets × 5 seeds is a pilot.  A powered reach comparison is Hamilton
  scale; this harness takes `nSeed` as its first argument and needs only the two
  libraries to run there.
- Cold-start invariance (no `tree =`) is verified separately and exactly by
  `dev/benchmarks/ab_bandit_coldstart_invariance.R`: **6/6 cells identical**
  across 2 datasets × 3 seeds, matching on score *and* on the full per-arm
  attempt vector — not merely on the final score.  This is the empirical
  backing for the NEWS claim that `adaptiveStart` alone is unaffected.
- Both libraries were built from the working tree as of ~07:25 on 2026-07-27.
  A concurrent session subsequently edited `ts_rcpp.cpp` and
  `R/MaximizeParsimony.R`; the A/B remains a valid paired comparison (the two
  arms differ *only* in the guard), but its absolute numbers predate those
  edits.  The test suite was re-run against the later source and passes.
- One run of the invariance harness died with SIGSEGV in the **parent** R
  process — which never loads TreeSearch — while a concurrent `R CMD INSTALL`
  was in flight.  It did not reproduce, and the re-run was clean 6/6.  Recorded
  here rather than omitted; if it recurs outside a concurrent install, it wants
  a proper look (cf. the open `ts_ratchet_search` multistate crash).

## `ts_parallel.cpp`: no mirror change needed

Its only contact with `StrategyTracker` is the **static, stateless**
`round_robin()` factory (`seq[r] = r % N_STRAT`, `ts_strategy.h`).  There is no
tracker instance, no `update()`, no `decay()` on that path — no bandit state
exists there to corrupt.  A warm-started rep does consume a round-robin slot,
which shifts (but does not bias) the arm cycle over the cold reps; that is
cosmetic and was left alone.
