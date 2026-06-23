# Briefing for #40: add PAIRWISE FUSE as a candidate strategy

**To:** the composition agent (#40), currently testing recipes against the full
training corpus.
**From:** the B1 / architecture-audit thread (2026-06-22).
**Ask:** add "pairwise fuse + suboptimal pool" as a new strategy *variant* in your
panel and score it against the corpus, alongside the existing presets. It is a
default-OFF, byte-identical-when-off prototype, so adding it cannot regress
anything you already measure.

---

## 1. What it is (one paragraph)

TreeSearch's tree-fusing (`tree_fuse`) recombines clades between pool trees, but
in the production driven search it does **zero productive recombination**: the
recipient is hard-wired to `pool.best()` (the already-optimised global best) and
acceptance requires a *strict* improvement, so a clade swap from a suboptimal
donor can never strictly improve it (root cause:
`dev/plans/2026-06-22-tree-fuse-zero-exchanges.md`). The **pairwise-fuse**
prototype fixes this: it also fuses into a few *suboptimal* pool recipients
(selected by score, worst-first, up to 4), the Goloboff-1999 configuration where
a better donor clade *can* strictly improve a weaker tree. Validated functional:
on 88-tip Dikow it recombines real suboptimal trees down to the optimum
(`Fuse improved: 1607 -> 1606`), and the default (OFF) path is byte-identical to
current behaviour (all `MaximizeParsimony` tests pass).

It is gated behind env flag **`TS_FUSE_PAIRWISE=1`** and **requires
`poolSuboptimal > 0`** (else the pool holds only tied-best trees, there are no
suboptimal recipients, and it silently degrades to the best-only no-op).

Code (uncommitted in worktree `condescending-hawking-50a5fc`):
`src/ts_driven.cpp` ~983 and `src/ts_parallel.cpp` ~51.

---

## 2. Honest evidence so far — prior is LOW, but the corpus is the right test

I measured it on the **5-dataset mission roster** (Wortley/Eklund/Zanol/Zhu/
Giles/Dikow, all ≤88 tips) and it is **neither a quality nor a robust efficiency
lever there**:

- **Quality:** plain multistart already reaches the TNT optimum on every roster
  dataset, so there is nothing above the optimum for recombination to close.
  Pairwise fuse adds 0 to the best score reached (wall-clock-matched A/B, n=16).
- **Efficiency (time-to-optimum, n=48, wall-clock-matched):** Dikow median
  4.61s (ON) vs 6.06s (OFF), faster on **52%** of seeds (a coin flip); Zanol
  speedup −0.01s. A faint early-window reach-rate edge on Dikow (the hardest
  roster dataset) that is within n=48 noise and gone by ~50s.
- **⚠ Sampling caution (please honour this):** a local **n=3** run showed a
  spurious **12× speedup** (0.66s vs 7.82s) that **completely vanished at n=48**.
  Time-to-optimum is high-variance; **use n≥32 per cell, wall-clock currency,
  and paired-by-seed comparison** or you will manufacture false positives. Reps
  and `candidates_evaluated` are both *unfair* currencies here (fuse replicates
  cost more; `tree_fuse` scores on an uncounted path).

**Why test it in #40 anyway** (the reason for this briefing):
1. Your corpus is **larger and harder than the mission roster** (the P1 auto-tuner
   ran 24 keys over hard datasets; `recipe-p1` memory). The roster's ceiling is
   88 tips and **no mission dataset is ≥120t**, so the `large`/`pruneReinsert`
   regime never fires here. Recombination's payoff is, by construction, confined
   to datasets where the optimum is **not** found in the first handful of restarts
   — i.e. bigger/harder matrices than I could test. The early-window reliability
   nudge that washes out at 88t may firm up there.
2. The honest open question is **reliability**, not median speed: does fuse cut
   the *variance / tail* of time-to-optimum (the unlucky-restart seeds), even if
   the median is unmoved? That needs the corpus's spread of difficulty to answer,
   and #40 is already paying for those runs.
3. It is cheap to add as one more strategy knob and costs nothing when off.

**Set expectations: I expect it to lose or tie on most classes.** The value of
including it is a clean corpus-wide verdict (and a possible win on a hard tail),
not a predicted improvement. If it shows nothing across the corpus, that closes
the thread for good.

---

## 3. How to add it to your panel

### Quick path (env flag — fine for a first screen)
Wrap the search call:
```r
withr::local_envvar(TS_FUSE_PAIRWISE = "1")   # or Sys.setenv + on.exit
MaximizeParsimony(d, ..., control = SearchControl(
  fuseInterval = 3L, poolSuboptimal = 5))     # poolSuboptimal>0 is MANDATORY
```
The flag is read per-fuse (not a process-static), so it scopes per call. Hold
`poolSuboptimal` identical between the fuse-on and fuse-off arms so you isolate
the fuse knob, **not** the pool-retention knob — note the production default is
`poolSuboptimal = 0`, so "pairwise fuse" is really a *fuse + suboptimal-pool*
package; frame the comparison as "enable the fuse+pool machinery vs plain
restarts."

### Clean path (promote to a `SearchControl` field — recommended if it earns a slot)
Mirror exactly how `fuseAcceptEqual` is plumbed:
- `R/SearchControl.R`: add `fusePairwise = FALSE` to the formals (~line 301), doc
  (~line 93), and the returned list (~line 388); add to the `"Fuse/Pool"` print
  group (line 432).
- `R/ts-driven-compat.R`: add to formals (~42) and list (~132).
- `src/ts_driven.h`: add `bool fuse_pairwise = false;` to `struct DrivenParams`
  (next to `fuse_accept_equal`, line 114).
- `src/ts_rcpp.cpp:1409`: `params.fuse_pairwise = as<bool>(ctrl["fusePairwise"]);`
- `src/ts_driven.cpp` ~983 and `src/ts_parallel.cpp` ~51: replace the
  `std::getenv("TS_FUSE_PAIRWISE")` read with `params.fuse_pairwise`.

Then a strategy preset can carry `fusePairwise = TRUE, poolSuboptimal = 5`.

---

## 4. Scoring guidance (so the verdict is trustworthy)

- **Metric:** time-to-known-optimum (anytime curve), not just final score —
  on most corpus datasets final score will tie at the optimum. Report the
  **reach-rate survival curve** (paired) and the **tail** (slowest-decile
  time-to-optimum), not only the median. Reuse `dev/benchmarks/fuse_efficiency.R`
  (`FeRun` records the `(elapsed, best-so-far)` trajectory via `progressCallback`;
  `FeAnalyze` builds the paired curve).
- **Targets:** take the optimum per dataset from
  `dev/benchmarks/headtohead_phase0.csv` (`tnt` column) or a TNT run under the
  identical `"-"→"?"` transform — **do not** retype a target table (a stale
  `.bestKnown` in `basin_diversity.R` recently held values from *different
  datasets*, e.g. Conrad2008's 1761 mislabelled as Zhu2013; see
  `dev/benchmarks/README-best-known-targets.md`).
- **Budget:** wall-clock-matched between arms; nThreads=1 (parallel makes
  wall-clock and pool dynamics nondeterministic).
- **Disposition:** land as a default only if it shows a real reliability/time
  win on some corpus class **and** survives a time-matched gate; otherwise leave
  it a default-OFF flag (or drop it). It has no quality benefit, so it should
  never displace effort that finds *better* trees — only effort that finds the
  *same* tree *sooner/more reliably*.

---

## 5. One-line summary

A validated, default-OFF recombination operator that finally makes tree-fusing
do real work; null on the small mission roster but plausibly useful on the
corpus's harder/larger datasets for *reliability/speed-to-optimum* (not quality)
— worth one slot in the panel, with a low prior and a strict time-matched gate.
