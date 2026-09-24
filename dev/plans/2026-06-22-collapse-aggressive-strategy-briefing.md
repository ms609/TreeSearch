# Briefing for #40: add AGGRESSIVE COLLAPSE (TNT `collapse 3`) as a candidate strategy

**To:** the composition agent (#40), testing recipes against the full training corpus.
**From:** the B2 / architecture-audit thread (2026-06-22; corpus-tested 2026-06-23).
**TL;DR — TESTED, NOT RECOMMENDED.** Aggressive collapse-during-search is a validated,
default-OFF, byte-identical-when-off prototype. It was hypothesised to speed up the corpus's
char-poor regime (where collapse density is genuinely high). A paired wall-clock A/B on that
exact regime (Hamilton job 17591484: 10 char-poor + 4 control datasets × 32 seeds) found **no
reliable speed win — a high-variance wash that helps some datasets and clearly hurts others
(project2144 ~1.6× slower), with no predictor of which.** Do **not** adopt it as a default or a
recipe knob. Documented here as a closed thread; the flag stays available (TS_COLLAPSE_AGGRESSIVE)
for the narrow case where someone has a specific char-poor dataset and measures a per-dataset win.

---

## 1. What it is

During TBR/SPR search, TreeSearch already skips clip/regraft candidates on *collapsed*
(zero-length) edges. The production criterion (`compute_collapsed_flags`) is **score-exact**:
it only collapses an edge when all regraft positions in the region provably produce identical
scores. On homoplasy-rich morphological data that criterion fires on **~0%** of edges (every
internal node carries downpass cost), so the machinery is inert there.

The **aggressive** criterion (`compute_collapsed_flags_aggressive`) instead collapses edges of
**minimum possible length 0** — TNT's `collapse 3` rule, an edge with *some* most-parsimonious
reconstruction that puts no change on it. This is a *heuristic neighbourhood reduction*
(Goloboff's asymmetric reachability): **scoring stays exact**, but a few improving moves may be
skipped, in exchange for a smaller, faster-to-search neighbourhood. It is the mechanism behind
TNT's documented collapse-search speedups on large/congruent matrices.

**Criterion (validated, not derived):** an internal edge (p,c) is collapsible iff
`final[p] & final[c] != 0` for every character (the up-pass MPR state sets share a state).
Validated **bit-for-bit against a brute-force MPR oracle**: 0 mismatches on all internal edges
over 120 trees, and crucially **0 false positives — it never over-collapses**
(`dev/benchmarks/b2_minlength_oracle.R`, `b2_collapsed_kernel_validate.R`).

Gated behind env flag **`TS_COLLAPSE_AGGRESSIVE=1`** (default OFF). Scoped to `tbr_search`'s
neighbourhood reduction only (plain TBR + sectorial + ratchet, all of which run through
`tbr_search`); **pool dedup keeps the exact flags**. Datasets with inapplicable characters fall
back to the conservative flags (NA soft-collapse is not yet derived), so NA data is unaffected.

**⚠ The lever only bites on NA-free data.** Because inapplicable-bearing datasets fall back to
the (inert, 0%) conservative flags, the aggressive criterion does real work **only** on
datasets with no inapplicable characters — i.e. **molecular** matrices, or morphological
matrices coded gap-as-missing (`"-"→"?"`, the standard-Fitch transform). A large *morphological*
dataset retaining inapplicables gets the inert treatment regardless of the flag. So test this on
NA-free large data; do **not** run it on an NA-heavy class and conclude "no win" — there it is
structurally a no-op, not a tested-and-failed lever.

Code (uncommitted in worktree `condescending-hawking-50a5fc`):
`src/ts_collapsed.cpp` (`compute_collapsed_flags_aggressive`), wired at `src/ts_tbr.cpp` ~1272
and ~2181, diagnostic export `ts_collapsed_flags_debug`.

---

## 2. Honest evidence so far — null on the morphological roster, by design

The mission roster (≤88-tip morphological) is **the wrong place to see a win**, and B2 is
*closed* there (see `dev/plans/2026-06-22-architecture-assumptions-audit.md`):

- **Exercised:** on random Zanol starts the aggressive flag collapses **18–22 internal edges**
  (≈28%) where the conservative flag collapses **0** — so the mechanism genuinely fires.
- **Density concentrates away from the work:** ~28% at a random start but only **~5.6% at the
  optimum**, where a parsimony descent spends the bulk of its candidate evaluations.
- **Measured work reduction on the roster ≈ 0** (the decisive speed datum): on identical
  random-start Zanol descents, `n_evaluated` ON/OFF ratio is **1.000 on 4 of 5 seeds**
  (wall-clock identical, ~0.15 s both arms); the 5th ran 1.41× *more* evals but reached a
  *better* score (1262 vs 1270) — a trajectory effect, not a speedup. So the neighbourhood
  reduction does **not** bite on this class: it fires at random starts but the descent spends
  its evaluations near the optimum where density is ~5.6%. Importantly the per-pass
  `O(n·blocks·n_states)` flag recompute did **not** inflate wall-clock (not net-slower either).
  The roster has no large NA-free dataset to exhibit the effect — so the speed premise was
  tested directly on the corpus (next bullet).
- **Corpus speed A/B (DEFINITIVE, 2026-06-23): no reliable win, a high-variance wash.** The
  mbank training corpus has a char-poor tail (14% < 1.5 chars/tip) where collapse density at
  the optimum IS high (confirmed: a char-subset of Dikow to 15 chars → 58% density). Ran the
  paired wall-clock-to-floor anytime A/B (ON vs OFF, 32 seeds, 90 s, `thorough`) on 10 char-poor
  datasets + 4 char-rich controls (Hamilton job 17591484, 896 cells). Result: **faster on ≤50%
  of cells, median speedup ≈ 0, quality tied** (dScore ≈ 0 both directions). Among the datasets
  whose search is long enough to measure, it SPLITS — project954 (83t/75ch) 14→11 s (helps),
  project3354 (78t/18ch) helps, but **project2144 (109t/124ch) 47→76 s (≈1.6× SLOWER + slightly
  worse)**. char/tip ratio does **not** predict winners (954 wins at 0.9, 2144 loses at 1.14), so
  there is no clean per-dataset selector. The asymmetric-reachability heuristic perturbs the
  trajectory in both directions; the high collapse density does **not** convert to a dependable
  wall-clock gain. Mechanistically grounded + corpus-tested ⇒ **not a reliable speed lever.**
- **Safe but no win on the roster:** wall-clock-matched runs on Zanol2014 (1261) and Dikow2009
  (1606) reach the **same optimum** with both arms (quality delta 0), and the returned-tree
  score always equals a fresh full rescore (heuristic safety confirmed).
- **Literature predicts the roster null:** the collapse family gave 50% time reduction on
  *congruent* 168-taxon data but **"no gain on incongruent data"** — morphological matrices are
  incongruent.

**Why test it in #40 anyway** (the reason for this briefing): the corpus is larger / more
molecular / more congruent than the roster, and *that* is the regime where zero-length branches
are common and collapse-search is documented to cut wall-clock substantially. This is a
**speed/efficiency** candidate (reach the same — already-reachable — optimum sooner via a
smaller neighbourhood), not a quality lever.

**Set expectations:** expect it to tie on incongruent/small classes and possibly **win on
large, NA-free, congruent classes** (molecular, or gap-as-missing morphological — see the NA
caveat in §1). The value is a clean corpus-wide verdict and a likely win on a class the roster
cannot exhibit.

---

## 3. How to add it to your panel

### Quick path (env flag — fine for a first screen)
```r
withr::local_envvar(TS_COLLAPSE_AGGRESSIVE = "1")
MaximizeParsimony(d, ...)        # works through plain TBR / sectorial / ratchet
```
Byte-identical to production when unset; no other knob needs changing.

### Clean path (promote to a `SearchControl` field — recommended if it earns a slot)
Mirror exactly how `fuseAcceptEqual` is plumbed:
- `R/SearchControl.R`: add `collapseAggressive = FALSE` to the formals, doc, the returned list,
  and a print group.
- `R/ts-driven-compat.R`: add to formals + the list passed to the driver.
- `src/ts_driven.h`: add `bool collapse_aggressive = false;` to `struct DrivenParams`.
- `src/ts_rcpp.cpp`: `params.collapse_aggressive = as<bool>(ctrl["collapseAggressive"]);`
- `src/ts_tbr.cpp` ~1269: replace the `std::getenv("TS_COLLAPSE_AGGRESSIVE")` read with the
  threaded-through param (note `tbr_search` is called by the driver/sector/ratchet — pass it on
  `TbrParams` or read from `ds`/params as those paths already do for other knobs).

Then a large-dataset preset can carry `collapseAggressive = TRUE`.

---

## 4. Scoring guidance (so the verdict is trustworthy)

- **Metric:** wall-clock **time-to-known-optimum** (anytime curve) — this is a SPEED heuristic,
  so reps and candidate counts are unfair currencies (it changes candidates-per-pass). Reuse
  `dev/benchmarks/fuse_efficiency.R` (`FeRun`/`FeAnalyze`). Report the paired survival curve and
  the tail, not just the median.
- **Quality gate is MANDATORY here** (unlike a score-exact change): aggressive collapse can skip
  improving moves (asymmetric reachability), so verify the best score reached is **not
  regressed** vs the flag-off arm on every corpus class, under a wall-clock-matched budget. If it
  reaches a worse optimum on any class, it must stay OFF for that class.
- **Targets:** optimum per dataset from `dev/benchmarks/headtohead_phase0.csv` (`tnt` column) or
  a matched TNT run — **do not** retype a target table (see
  `dev/benchmarks/README-best-known-targets.md`).
- **Budget:** wall-clock-matched between arms; nThreads=1.
- **Disposition:** land as a default (or a large-dataset preset) only if it shows a real
  wall-clock win on some corpus class **and** does not regress quality there; otherwise leave it
  a default-OFF flag.

---

## 5. One-line summary

A validated, default-OFF, never-over-collapsing TNT-`collapse 3` neighbourhood reduction; inert
on the char-rich roster *by construction*, and — when finally tested on the corpus's char-poor
regime where its density IS high — **a high-variance wash, not a reliable speed win** (helps
some datasets, ~1.6× slower on others, no selector). Quality always tied. **Closed: not a
recommended lever.** Flag retained (TS_COLLAPSE_AGGRESSIVE) for narrow per-dataset use only.
