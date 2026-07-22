# Best-known parsimony targets — provenance, regimes, and a cautionary tale

**TL;DR.** When you need a "best-known score" to gate a benchmark, read it from
`dev/benchmarks/headtohead_phase0.csv` (column `tnt`). Do **not** retype a
literal target table — that is how wrong values creep in. A target is only
meaningful paired with its **scoring regime** and its **exact dataset**; mix
either up and the number is noise.

## The canonical table

`headtohead_phase0.csv` holds, per dataset (2 seeds each), the scores TS and TNT
reach under the **apples-to-apples Fitch regime** (see below):

| dataset | tips | `ts_fitch` | `tnt` (best-known) |
|---|---|---|---|
| Wortley2006 | 37 | 481–483 | **479** |
| Eklund2004 | 54 | 440 | **440** |
| Zanol2014 | 74 | 1264–1265 | **1261** |
| Zhu2013 | 75 | 626–627 | **624** |
| Giles2015 | 78 | 671–672 | **670** |
| Dikow2009 | 88 | 1606 | **1606** |

These are corroborated by an independent TNT run (`xmult hits 10 replic 10`,
2026-06-22) which reproduced every `tnt` value exactly, and by every other
benchmark driver in this directory (`bench_iterate.R`, `bench_beam.R`,
`diag_convergence_ab.R`, `diag_collapse_sect.R`, …).

## Scoring regimes — they give *different* numbers for the same tree

The `inapplicable.phyData` datasets contain inapplicable tokens (`-`). How you
treat them changes the score, and the regimes are **not** interchangeable:

- **`"-"→"?"` (missing) — the "apples-to-apples Fitch" regime used here and by
  TNT.** Inapplicables become missing data, so the search runs *standard Fitch*.
  This is the most permissive coding ⇒ the **lowest** achievable score. It is the
  only regime in which TS and TNT scores are directly comparable.
- **Gap-as-extra-state / native-inapplicable (de-Laet–Brazeau).**
  Inapplicables are scored as structure, forcing extra steps ⇒ a **higher**
  score than the `"?"` regime for the same tree. `ts_raw` in the CSV is an
  example of a stricter column (e.g. Zanol `ts_raw`=1315 vs `ts_fitch`=1265).

**Directional rule of thumb:** any inapplicable-aware regime is `≥` the `"-"→"?"`
Fitch score. So a candidate "best-known" that is *lower* than the Fitch optimum
cannot be a re-scoring of that dataset — it belongs to a different dataset.

## The cautionary tale (2026-06-22)

`basin_diversity.R` (B1 basin-diversity harness) hardcoded a `.bestKnown` table
that had drifted to:

| dataset | wrong value | true Fitch optimum | what the wrong number actually was |
|---|---|---|---|
| Zhu2013 | 1761 | 624 | **Conrad2008**'s score (64t/360c — a *different dataset*), under a gap-aware regime |
| Giles2015 | 458 | 670 | foreign: *below* the Fitch optimum ⇒ impossible as any re-scoring of Giles2015 |
| Dikow2009 | 1075 | 1606 | foreign: likewise below the Fitch optimum |

Three independent checks pinned the diagnosis:
1. **Provenance grep.** `1761` appears in the repo only as Conrad2008's score
   (`t249_*.csv`, `t264_verify_*.csv`). It is not Zhu2013's anything.
2. **Roster sweep.** A bounded search over all 30 `inapplicable.phyData`
   datasets under `"-"→"?"` found *no* dataset scoring within ±8 of 458, 1075, or
   1761. They are not this-regime optima of anything in the package.
3. **Canonical CSV.** `headtohead_phase0.csv` — the source `basin_diversity.R`
   *claimed* to mirror — already held the correct 624 / 670 / 1606.

The wrong targets did **not** corrupt the B1 verdict (that is a *within-regime*
A/B of pairwise-fuse on vs off — target-independent), but they made the
reach-rate gate meaningless and nearly sent the analysis chasing a phantom
"TS plateaus far above optimum" gap. Only a TNT ground-truth run recovered the
true optima.

## The fix + the rule

`basin_diversity.R` now derives `.bestKnown` from `headtohead_phase0.csv` at load
(`.BdBestKnown()`), with a hardcoded fallback kept in sync with the CSV's `tnt`
column. The standing rules:

1. **Derive targets, don't retype them.** Read the CSV; if you must hardcode,
   cite the CSV row and the regime in a comment.
2. **A score without a (dataset, regime) pair is meaningless.** Never copy a bare
   number between tables.
3. **Sanity-check direction.** If a "best-known" is *below* the `"-"→"?"` Fitch
   optimum, it is from a different dataset or a stricter→looser regime mix-up —
   stop and trace it before trusting any gate built on it.
4. **When in doubt, run TNT under the identical transform.** It is ground truth
   and makes stale recorded targets irrelevant.
