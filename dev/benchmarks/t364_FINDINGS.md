# T-364/T-370 merged-fix validation: three-arm constrained-search battery

2026-07-29/30. Harness: `t364_*.R` / `t364_*.sh` in this directory.

Two questions:

1. Does the merged fix (`796a29d3`) cost anything on constrained `MaximizeParsimony()`
   across many matrices?
2. Was the old "~1.43× wall cost" ever real?

**Answers: (1) No. (2) Yes — but it belongs to complement-enforcement WITHOUT the
reroot, not to the merged fix, and not to the accidental-Wagner-restarts story.**

## Arms

| arm | commit | state |
|---|---|---|
| 1 | `bbcca1ba` | pre-fix (parent of `7685bf07`) |
| 2 | `7685bf07` | complement enforcement only, **no reroot** — T-384 exposed |
| 3 | `796a29d3` | complement + reroot at tip 0 — the merged state |

All three predate the T-384 fix deliberately: the mechanism must still be live for
arm 2 to express it. Three independent clones, 34 TUs each, three distinct `.so`
md5s. The brief said "compare against its parent state"; `796a29d3`'s literal parent
is `7685bf07`, so both `3v1` and `3v2` are reported rather than picking a reading.

## Corpus and constraints

Frozen 25-matrix `MBANK_FIXED_SAMPLE`, 21–4062 tips, **training split only** (the
loader hard-refuses a validation key). Constraints generated once against arm 3,
frozen to RDS, consumed identically by every arm — no arm can be advantaged by a
constraint it helped choose.

Shape: **clade + distant rogue** — take a real group from the matrix's own
unconstrained MP tree and add the topologically most distant tip, forcing the search
to drag a rogue in. Moderate, realistic conflict, not a maximally-interleaved group
(which on 173 tips forces an enormous penalty and a search that may not converge).

**Non-triviality: every constraint has a strictly positive score penalty** at matched
budget, +3 to +66 steps. Absence from one MP tree was *not* the acceptance test — a
split absent from one MP tree can sit in another tree of the same island at zero cost.

**Stratified on which side holds dataset tip 0**, because that is what decides
whether arm 2 can fail at all: `build_constraint()` canonicalises the mask so tip 0
is outside, so the canonical side is always the tip-0-excluding one, and the
pre-T-384 mapping needs *that* side to be a rooted clade — which fails exactly when
the root position lies inside it.

* `rogue_no0` — tip 0 in the large side ⇒ canonical side small ⇒ root rarely inside.
* `rogue_in0` — tip 0 in the small side ⇒ canonical side large ⇒ root usually inside.

For a *k*-taxon group in *n* taxa, `in0` arises with probability ≈ *k/n*: `no0` is
the common geometry, `in0` a real but minority one. **Neither is contrived — which
side holds `names(dataset)[1]` is plain matrix row order, and no user controls it.**

Groups are drawn from the reference tree's **unrooted splits, both sides**, not its
rooted clades: `MaximizeParsimony()` returns trees rooted at tip 0, so no non-root
clade ever contains tip 0 and a rooted-clade generator cannot reach the `in0`
stratum at all. The first version of the generator hit exactly this.

## Q1 — the merged fix costs nothing

423 runs, 47 (matrix × shape) pairs, 3 seeds, three arms interleaved on ONE node per
cell in an order rotated by task id. **Zero censored runs.**

| arm3 vs arm1 (MERGED vs pre-fix) | nMatrix | median ratio | slower/faster/tie | >10% slower | sign p |
|---|---|---|---|---|---|
| wall, budget 24 — all 25 matrices | 47 | **1.000** | 22 / 22 / 3 | 8 | 1 |
| wall, budget 96 — **small+medium only** | 27 | **0.967** | 10 / 17 / 0 | 4 | 0.25 |
| score, budget 24 — all 25 matrices | 47 | **+0** | 2 worse / 2 better / 43 tie | — | 1 |
| score, budget 96 — small+medium only | 27 | **+0** | 2 / 0 / 25 | — | 0.5 |

⚠ **Scope of every budget-96 figure in this document: small+medium only** (`TIERS` in
`t364_array96.sh`) — 14 matrices, 27 pairs, and the `in0` rows only 13 pairs.
Large/xlarge were not run at budget 96, so neither `0.967` nor the `2.531` below is a
corpus-wide number.

⚠ **The budget-24 `no0` comparison is budget-bound and near-vacuous.** All three arms
ran 24 / 24 / 24 replicates there (exhaustion 64% / 65% / 67%), so wall is roughly
equal by construction. `in0` at budget 24 *was* convergence-limited (16 / 24 / 17.5),
and the budget-96 run un-binds `no0` (31.5 / 29 replicates, not exhausted, arm3/arm1
0.969) — that is the informative figure for that stratum. This is a caveat on the
headline, not a refutation of it.

Median + sign count + count >10% slower, never the mean. The 8 cases >10% slower at
budget 24 are balanced by **9 cases >10% faster** (down to 0.456×), i.e. symmetric
noise rather than a concentrated regression. The 4062-tip matrix agrees across arms
to within 1%.

Compliance: arm 3 **100%** of 2634 returned trees. Arm 1 **99.96%** — one violating
tree in 2751, the rejection sampler exhausting its 100 attempts. So pre-fix the
*search* path was very nearly, but not perfectly, shielded.

**Read the 1.000 as two effects cancelling, not as "nothing happens":** arm 1 pays
rejection-sampler reshuffles (≈1/(1−p) Wagner builds per start — ~1.06× on `no0`,
~5× on `in0`), arm 2 pays blocking, arm 3 pays neither. On a corpus with a different
geometry mix the cancellation would not be exact.

## Q2 — the 1.43× is reattributed, not retracted

### Mechanism: arm 1's violating set and arm 2's complement-rooted set are the same set

21 matrices × 200 addition orders, supplied **explicitly and identically to every
arm** so the comparison is exactly paired and cannot be confounded by arms consuming
the RNG stream differently. The two distributions match in median **and** range:

| stratum | arm 1 violating | arm 2 complement-ONLY | arm 3 |
|---|---|---|---|
| `in0` | 81.2% [71.0–94.0] | **81.2% [71.0–94.0]** | 0% |
| `no0` | 6.0% [0.0–11.5] | **6.0% [0.0–11.5]** | 0% |

On the T-364 test case verbatim, the recorded 35/400 = 8.75% pre-fix violation rate
reproduces **to the digit**, and arm 2's complement-only rate on that case is *also*
35/400, on the same orders.

So an addition order that used to skip the constraint ran an effectively
*unconstrained* search — cheap, converging early — and under complement-enforcement
alone the identical order runs a fully move-blocked one. **That substitution is the
wall gap.** No accidental Wagner restarts are needed to explain it, which matters
because that attribution was wrong twice over: `AdditionTree()` has
`has_posthoc = FALSE` so it never retried, and the Wagner build is ~1–2 ms, far too
small to source a 1.43×.

**Corollary: the recorded "~1 addition order in 11" understates T-364 badly.** That
is the common-geometry figure; on `in0` constraints `AdditionTree()` violated up to
**94%** of the time.

### Magnitude, and why the budget decides it

| arm2/arm1 wall | budget 24 | budget 96 |
|---|---|---|
| `in0` (affected geometry) | 1.010 | **2.531** (11/13 slower, 11 >10%, p=0.022) |
| `no0` (common geometry) | 0.985 | 1.053 |
| pooled | 0.995 | **1.211** (19/8, 16 >10%, p=0.052) |

The ratio **grows with the budget**, which discriminates the two hypotheses: budget
exhaustion predicts growth, uniform per-unit slowness predicts invariance. Replicates
run on `in0` at budget 96: arm 1 **15**, arm 3 **17**, arm 2 **65**. Budget-24
exhaustion rate on `in0`: arm 1 50%, arm 3 47%, arm 2 **95.5%**.

`arm3/arm2` on `in0` at budget 96 = **0.382** (1/12 slower, p=0.0034) — the reroot
removes the cost entirely.

The historical 1.43× lies between the pooled 1.21× and the stratum 2.53×. An exact
match is neither expected nor needed: the T-214 battery used different matrices
(10–15 tips, 3 matrices) and an unrecorded replicate budget.

### Why it appears as score at a small budget and as wall at a large one

At budget 24 arm 2's wall is flat, but its **score is worse on 6 of 22 `in0`
matrices, 0 better (p=0.031)**, and arm 3 beats it on the same 6 (p=0.031). Given
more budget it recovers the score — and pays 2.53× to do so.

The original 1.43× was measured as *"wall to match the old score"* — holding score
fixed. This battery holds budget fixed. Same phenomenon, two ways of holding things
constant. That reconciles the old number rather than contradicting it.

## ⚠ The recorded T-384 diagnostic signature is wrong

Recorded, and flagged "worth reusing": *"every phase (TBR, XSS, ratchet) slower on a
score-IDENTICAL trajectory ⇒ moves being rejected."*

**Its sign flips with the replicate budget, on the same code and the same defect.**
Cumulative per-phase `arm2/arm1`:

| budget | verdict |
|---|---|
| 24 | arm 2 **faster** in 6 of 7 live phases (ratchet 0.79×, tbr 0.44×, wagner 0.39×, fuse 0.24×; only xss 1.26× up) |
| 96 | arm 2 **slower** in 6 of 7 (ratchet 1.60×, xss 1.95×, rss 1.67×, final_tbr 1.65×, fuse 2.10×) |

A test that reverses its verdict when a budget setting changes cannot be diagnostic.
The arithmetic reason: phase totals are summed over however many replicates ran, and
the replicate *count* is what the defect changes.

**The invariant that IS diagnostic, holding at both budgets: TBR work per replicate
roughly halves** in the blocked arm — 0.44× cumulative at budget 24, 0.498×
per-replicate at budget 96. That is what rejected moves must do: less work per
replicate, not more. Per-replicate total wall drops too (0.032 s vs 0.058/0.044).

**Correct signature: more replicates, each doing less work, same score.** A blocked
replicate never registers a hit on the best score, so the convergence rule never
trips and the search burns its entire budget. Never compare cumulative phase totals
across runs whose replicate counts differ.

## Reproducing

```
sbatch --export=ALL,ARM=1 t364_build.sh   # and ARM=2, ARM=3
sbatch --array=1-25 t364_genarray.sh                              # freeze constraints
sbatch --dependency=afterany:<gen> --array=1-3 t364_probearray.sh # probe P / discriminator
sbatch --dependency=afterany:<gen> --array=1-150 t364_array.sh    # battery, budget 24
sbatch --dependency=afterany:<batt> --array=1-150 t364_array96.sh # budget 96, small+medium
BATT_DIR=.../batt Rscript t364_analyze.R
Rscript t364_crosscheck.R   # reproduces the recorded 35/400, per arm
```

Run probe P **before** the battery: it is the gate. If arm 2's complement-only rate
is ~0 on the frozen constraints, arm 2 cannot express T-384 and question 2 is
unanswerable — and you would only discover that after spending the array. It is also
the behavioural discriminator proving the three libraries differ; three silently
identical libs would produce a clean and meaningless "no cost".

Classification is by direct edge-matrix descendant-set accumulation on tip **labels**.
Never `as.Splits` + `%in%`: S4 `%in%` on Splits silently falls through to `base::%in%`
when TreeTools is unattached and answers FALSE.
