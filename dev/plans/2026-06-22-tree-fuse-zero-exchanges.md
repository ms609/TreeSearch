# tree_fuse: zero accepted exchanges in the driven search — root cause

**Date:** 2026-06-22
**Engine:** cpp-search (worktree), package 2.0.0.
**Defect:** with `verbosity=2`, every periodic fuse prints
`Fuse attempt: pool=N exchanges=0  <score> -> <score>` — `tree_fuse` fires but
recombines nothing, on Wortley2006 (37t), Eklund2004 (54t), Zanol2014 (74t), at
all pool sizes 2..10. No crash; not size-specific.

## TL;DR root cause

`exchanges` counts only **accepted** exchanges, and acceptance under the default
preset (`fuseAcceptEqual = FALSE`) requires a **strict score improvement**
(`new_score < score`, `ts_fuse.cpp:466`/`481-484`). On this dataset class a clade
swap **never strictly improves** the recipient, because the recipient is always
`pool.best().tree` (the already-TBR-optimised global best) and the donors are
within `poolSuboptimal` steps. So every trial lands in one of three rejected
buckets and `n_exchanges` stays 0.

**This is intended-but-too-strict behaviour, not a logic bug.** Split matching,
donor iteration, indexing, and `replace_subtree` are all working correctly. The
only condition never satisfied is the strict-improvement gate.

## Evidence (instrumented build, TS_FUSE_TRACE, reverted)

Temporary getenv-gated counters were added to `tree_fuse` (donors / shared splits
found / trials / accepted / rejection-bucket breakdown), built into `.agent-fuse`,
and run on the exact repro (`set.seed(3)`, Wortley2006, `fuseInterval=2`,
`poolSuboptimal=5`). Instrumentation has been fully reverted.

**Default preset (`accept_equal = FALSE`):** e.g. pool=9 attempt:
```
donors=9 shared=155 trials=155 accepted=0
| rej_worse=29  rej_equal_notopo=69  rej_equal_topo_acceptEqualOff=57
```
- `shared=155`, `trials=155` ⇒ split matching and exchange application both fire
  fine (refutes hypotheses 1, 3, 4, 5).
- `accepted=0` ⇒ no trial was a strict improvement.
- Of the 155 trials: 29 made the score **worse**; 69 were score-equal **with the
  donor's clade topology identical** to the recipient's (a literal no-op swap —
  low pool diversity); **57 were score-equal with a genuinely different clade
  topology** that `accept_equal=FALSE` discards.

**`fuseAcceptEqual = TRUE` (thorough/large preset):** same seed/data:
```
rounds=10 donors=15 shared=400 trials=215 accepted=10
| rej_worse=0 rej_equal_notopo=205 rej_equal_topo_acceptEqualOff=0
Fuse attempt: pool=3 exchanges=10  480 -> 480   (and 479 -> 479)
```
- `exchanges` is now **non-zero** (≈10/attempt × ~7 attempts ≈ the ~70 figure the
  prior audit #55 reported), confirming the gate is the whole story.
- But score is still `480 -> 480` / `479 -> 479`: **even accepting equals, fuse
  produces 0 score improvements** on this class — it shuffles equal-cost clade
  topology only. Confirms the audit's "fuse is dead weight on the >64t class".

## Where the gate lives (`src/ts_fuse.cpp`)

```
465   bool accept = false;
466   if (new_score < score) {
467     accept = true;
468   } else if (params.accept_equal && new_score == score) {
          ... // require the clade topology actually changed
478     if (changed) accept = true;
479   }
...
481   if (accept) {
483     score = new_score;
484     ++result.n_exchanges;   // <-- only path that increments
```
`accept_equal` comes from `fuseAcceptEqual`: `default`/`sprint` set it FALSE
(`R/MaximizeParsimony.R:112,128`), so the strict-improvement path at line 466 is
the only one available — and it is never taken for `recipient = pool.best()`.

Recipient is hard-wired to the pool best in both paths:
`ts_driven.cpp:962  TreeState fused = pool.best().tree;` and
`ts_parallel.cpp:51 TreeState fused = pool_.best().tree;`.

## Why strict improvement essentially never happens here

The premise of tree-fusing (Goloboff 1999) is that a *suboptimal* donor may hold a
*locally better* arrangement of a shared clade than the recipient. That requires
either (a) a recipient that is **not** already optimal in that clade, or (b)
equal-cost acceptance to migrate clade topology between basins. TreeSearch's driven
search uses (i) the **global best** as the sole recipient and (ii) **strict
improvement** by default — the intersection of the two leaves no room: the best
tree's shared clades are already locally optimal, so donor clades are either worse
(rej_worse) or score-neutral (rej_equal). With small/clean morphological matrices
the search is thorough enough that the strict-better case has measure zero.

## Relationship to the >64-tip reroot/intraFuse crash

**Distinct defects, same module, not causally linked.**
- The crash (`intraFuse` on >64t) is a *correctness/safety* problem: round-≥2 TBR
  moves tip 0 out of the root, flipped clades spuriously match a donor's
  complement, and `replace_subtree` corrupts the tree (segfault at `wps>=2`). The
  fix is per-round `reroot_at_tip0` + the `r_rest.size()!=d_rest.size()` guard
  (`ts_fuse.cpp:379, 256`) — and **both are already present in this worktree's
  `ts_fuse.cpp`**, so matching here is correct (155 valid shared splits on 37t,
  `wps=1`, no spurious matches).
- Zero-exchanges is an *acceptance-gate* property, independent of rooting. It
  reproduces at 37 tips where the reroot bug cannot fire. Fixing the reroot bug
  does nothing for zero-exchanges and vice-versa.

## Is it a bug? Verdict

**Intended-but-too-strict, AND the algorithm is structurally hamstrung for this
class.** The gate behaves as written. The real issue is design:
recipient-always-best + improvement-only ⇒ fuse cannot contribute on datasets where
the best tree is already locally optimal in every shared clade. Even the
`accept_equal` path only reshuffles equal-cost topology (0 score wins here), so
fuse currently adds cost (an extra `tree_fuse` + `score_tree` every `fuseInterval`
reps) for no quality benefit on the mission class. This corroborates audit #55
(`dev/plans/2026-06-21-search-switches-for-composition.md:142`) and the
fuse-isolation note (`dev/plans/2026-06-20-fuse-drift-isolation.md`).

## Proposed fixes (concrete; not applied)

Pick per goal. These are mutually compatible.

1. **Make fuse able to win — fuse pairwise, not best-only (the real fix).**
   In both call sites, loop the recipient over a few of the *suboptimal* pool
   trees (not just `pool.best()`), keep the best fused result, and add it to the
   pool. A suboptimal recipient has clades that a better donor *can* strictly
   improve — the configuration where Goloboff fusing actually pays. Sketch
   (`ts_driven.cpp` ~962 and `ts_parallel.cpp` ~51):
   ```cpp
   double best_fused = pool.best_score();
   TreeState best_fused_tree;
   const auto& es = pool.all();
   int n_recip = std::min<int>(es.size(), 4);      // a few worst-but-pooled
   for (int ri = 0; ri < n_recip; ++ri) {
     TreeState cand = es[ri].tree;                  // includes suboptimal trees
     FuseParams fp; fp.accept_equal = params.fuse_accept_equal; fp.max_rounds = 10;
     FuseResult fr = tree_fuse(cand, ds, pool, fp);
     double s = score_tree(cand, ds);
     if (s < best_fused) { best_fused = s; best_fused_tree = std::move(cand); }
   }
   ```
   (Cost-bounded by `n_recip`; gate behind a flag and time-match before any
   default flip — composition #40.) Self-as-donor is already a harmless no-op
   (the recipient's own splits match identically ⇒ `rej_equal_notopo`).

2. **Cheapest correctness/clarity fix — let the default gate accept equals when it
   changes topology (i.e. flip `default`/`sprint` `fuseAcceptEqual` to TRUE).**
   This makes `exchanges>0` and propagates diverse equal-cost clades into the pool
   (more donor variety for future ratchet/sectorial escapes). It does **not** by
   itself produce score wins on this class (shown above: 480→480), so it is a
   diversity lever, not a quality fix — but it removes the misleading
   `exchanges=0` and is one-line per preset.

3. **If fuse stays best-only + strict (status quo intent): stop paying for it on
   this class.** Raise `fuseInterval` (or set 0) for the auto-selected `default`
   class to reclaim the per-interval `tree_fuse`+`score_tree` overhead, since it
   provably yields nothing here. (Already half-true: `poolSuboptimal=0`→pool size
   1→fuse skipped. Make it explicit.) This is the audit #55 recommendation.

**Recommendation:** (1) is the only change that lets fuse *improve scores*; do it
behind a flag and time-match under composition #40. (2) is a safe immediate
unblock for the misleading diagnostic. (3) is the honest do-nothing if fuse is
deemed not worth it for the mission class.

## Minimal reproduction / confirmation

```r
library(TreeSearch); library(TreeTools)
e <- new.env(); data("inapplicable.phyData", package = "TreeSearch", envir = e)
m <- PhyDatToMatrix(e[["inapplicable.phyData"]][["Wortley2006"]], ambigNA = FALSE)
m[m == "-"] <- "?"; d <- MatrixToPhyDat(m); set.seed(3)
MaximizeParsimony(d, maxReplicates = 10L, targetHits = 99999L, nThreads = 1L,
  verbosity = 2L, control = SearchControl(fuseInterval = 2L, poolSuboptimal = 5))
# default preset -> "Fuse attempt: pool=N exchanges=0  s -> s"
# add fuseAcceptEqual = TRUE -> exchanges>0 but score still s -> s (no improvement)
```
To re-derive the bucket counts, re-add the `TS_FUSE_TRACE` counters around the
accept block in `tree_fuse` (`ts_fuse.cpp:460-520`) and read stderr.

## Notes on the older `test_fuse_log.txt`

The repo-root `test_fuse_log.txt` ("Pool scores: 9 10 9 / Fuse done! Score: 13")
shows a fuse producing a *worse* score than the pool best — an artefact of an
older/standalone harness predating the current accept-gate + `score_tree`
post-check at `ts_driven.cpp:968`; not representative of the current driven path.
```
