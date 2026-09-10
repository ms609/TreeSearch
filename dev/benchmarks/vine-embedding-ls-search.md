# Vine-style embedding search: ruled out for `LeastSquaresTree()` and for sectors

**Status 2026-08-01: NEGATIVE RESULT. Do not re-investigate without new evidence.**

Prompted by Siepel, Hassett & Staklinski, *VINE: Variational inference for scalable
Bayesian reconstruction of species and cell-lineage phylogenies*
(bioRxiv `10.64898/2025.12.24.696405`). Vine embeds taxa in a continuous space,
decodes a tree by neighbour-joining, and backpropagates a likelihood gradient
through the decoder, reporting one-to-two orders of magnitude speed-up over
MCMC-based Bayesian phylogenetics.

Scripts: `dev/benchmarks/vine_ls_probe{1,4,6,7}.R`, `dev/benchmarks/vine_sector_test.R`.
All runs `nThreads = 2L`; nothing in the package was modified.

## Why this was worth testing at all

Vine's core machinery exists to make a *discrete* object (tree topology) continuously
reparameterisable. TreeSearch's parsimony objective is integer-valued and has no
branch lengths, so the chain rule terminates immediately and the method cannot be
ported. `LeastSquaresTree()` is the one place in the package with a *continuous*
objective — `RSS = Σ w_ij (d_tree(i,j) − D_ij)²` — so Vine's decoder-plus-gradient
should apply with no relaxation. It does not, for the reasons below.

## Target family (and a design trap to avoid repeating)

First attempts built the target as the average path-length matrix over an MPT set from
`MaximizeParsimony()`. **That is confounded** — the search is stochastic, so the MPT
count and therefore the problem's difficulty changes between runs (Wortley2006 gave 19
MPTs in one probe and 2 in the next, moving NJ's RSS from 33.7 to 1.5). Worse, MPTs of
a *single* dataset are all similar, so their average is nearly additive and NJ almost
solves it outright.

Use instead a reproducible, tunably hard family: backbone tree `T0` plus 11 copies each
perturbed by `m` random SPR moves; target = mean path-length matrix; conflict set by
`m`. All results below use `method = "nnls"`, the default and the Lapointe & Cucumel
convention.

## Findings

**1. The real objective has no gradient.** `LeastSquaresTree()` refits branch lengths
on each candidate, discarding the decoder's. RSS-after-refit therefore depends on the
input matrix only through the discrete topology and is piecewise constant: along a
1-D path in D-space, 241 grid points give **6 distinct values, median adjacent step
exactly 0**. The only differentiable objective is the *surrogate* — score NJ's own
branch lengths, which is what Vine actually does (241 distinct values, median step 0.4).

**2. The surrogate misranks topologies, and it is not a scale artifact.** Among
distinct immediate neighbours (RF == 2), Spearman against the true objective over five
targets: +0.10, +0.36, +0.22, +0.27, +0.17; its argmin lands in the true top decile 2
times in 5. Scale- and affine-corrected surrogates shift Spearman by **< 0.02**.
(A pooled ρ = 0.87 across perturbation scales is Simpson's paradox on ε — ignore it.)

**3. The decoder's neighbourhood is benign.** Perturbation scale is a clean dial on
topological step size: P(0 < RF ≤ 4) peaks at 0.65–0.71 around ε = 0.15–0.20,
consistently at 20/37/40 tips. The move mechanism is not the problem.

**4. Undirected D-space sampling finds the same tree, verified by identity.** The best
tree from random perturbation is the *same tree* NNI+SPR converges to (RF = 0) on every
target. An oracle control — greedy hill-climbing on the **true** objective, i.e. the
best any gradient could do — tied NNI+SPR on 3 of 5 targets and lost on 2. It never won.

**5. Honest cost of a self-terminating sampler**, measured in candidate evaluations
(both methods pay one NNLS refit per candidate), 8–10 replicate chains per cell, halting
after `k` draws without improvement:

| target | NNI+SPR cost | patience | P(hit best known) | sampler evals | ratio |
|---|---|---|---|---|---|
| n=20, m=8 | 705 | 50 / 200 / 800 | 0.30 / 0.90 / 1.00 | 62 / 322 / 1085 | 11.4× / 2.2× / **0.65×** |
| n=40, m=2 | 4 429 | 50 / 200 / 800 | 0.75 / **1.00** / 1.00 | 82 / 269 / 1152 | 54× / **16.5×** / 3.8× |
| n=40, m=8 | 7 384 | 200 / 800 / 2000 | 0.00 / 0.00 / 0.13 | 257 / 1463 / 2500 | 28.8× / 5.0× / 3.0× |

Reliability splits by **difficulty, not size**. Untuned ε ~ U(0.05, 0.30) vs tuned
ε = 0.15: P(hit) 0.75 vs 1.00 at patience 200, both 1.00 by patience 800.

**6. A pre-pass cannot cash that in, because the cost is certification.** Seeding
`LeastSquaresTree()` from the best of 150 perturbed-NJ decodes gives speed-ups of
0.69×, 0.92×, 1.03×, 1.05×, 1.07×, 1.21× over six targets — noise around 1.0. At n = 60,
730 s default vs 695 s seeded. **A better start tree buys nothing.** See the follow-up
lead below.

## Does it transfer to parsimony? No — and the reason generalizes

Everything above is a least-squares objective, where the decoder and the objective are
natively aligned (NJ *is* a distance method; the objective *is* distance fit), so the
result is partly self-fulfilling. Tested directly against the slot it would occupy:
perturbed-NJ (pNJ) vs random-addition-sequence Wagner as the start generator in
`build_ras_sector()` (`src/ts_sector.cpp:884`, called in the `ras_starts` loop at
`:1186`). pNJ competes with RAS, **not** with the TBR that follows either way.

16 sectors (8–37 tips) carved as clades from parsimony trees of Sansom2010,
Wortley2006, OLeary1999, Griswold1999; K = 5 restarts; identical fixed TBR budget
(`maxIter = 100`); equal weights. Excess = steps above the best score known for that
sector; diversity = mean pairwise normalized clustering-information distance across
restarts.

| generator | start excess | post-TBR excess | best-of-5 | P(hit sector best) | start diversity |
|---|---|---|---|---|---|
| RAS Wagner | 4.0 | **1.0** | **0.0** | **0.56** | **0.328** |
| pNJ ε = 0.15 | **2.5** | 1.5 | 1.0 | 0.44 | 0.162 |
| pNJ ε ~ U(.05,.30) | 3.0 | 1.5 | 0.5 | 0.50 | 0.201 |

Paired by sector on best-of-5: pNJ ε=0.15 beat RAS on 1, tied 7, **lost 7**;
random-ε beat RAS on 2, tied 9, lost 4.

**pNJ makes individually better start trees — 2.5 vs 4.0 steps above optimum — and
loses anyway, because its restarts are half as diverse (CID 0.162 vs 0.328).** RAS's
randomness is the point.

The mechanism corroborates itself: the *more* diverse pNJ arm (random ε, CID 0.201)
beats the less diverse one (ε = 0.15, CID 0.162) on every measure — best-of-K 0.5 vs
1.0, paired record 2/9/4 vs 1/7/7. Diversity tracks the outcome, not decoder quality.
Plain unperturbed NJ scores 2.5 excess, so ε = 0.15 barely moves away from plain NJ.

**Scope of that conclusion:** multi-start *sector solving* is diversity-limited rather
than start-quality-limited, which is reason to expect the same in other start-tree-driven
multi-start contexts (driven-search replicates). It does **not** extend to the ratchet,
whose diversity comes from perturbing the objective via character reweighting — pNJ is
not competing for that slot at all.

### Caveats
Equal weights only (implied weighting untested). Sectors modelled as the reduced dataset
on a clade's taxa, omitting the HTU composite terminals a real RSS sector carries.
Polished with the R-level `TreeSearch()` TBR, not `ts_sector.cpp`'s internal search with
its ratchet/drift options — a stronger internal search would wash out start-tree
differences further, strengthening the negative. 16 sectors, K = 5: enough for a 7-vs-1
paired split, not for a small effect.

---

## Two follow-up leads, unrelated to Vine

**(a) `LeastSquaresTree(method = "ols")` is a trap for the documented use case.** Its
winning topology, rescored under NNLS, is **worse than the untouched NJ start** on 5 of
5 targets:

| target | NJ (NNLS) | NNLS-search winner | OLS-search winner, rescored NNLS |
|---|---|---|---|
| sim n=20, m=2 | 25.884 | 25.884 | **38.333** |
| sim n=20, m=8 | 64.867 | 63.930 | **74.142** |
| sim n=40, m=2 | 273.361 | 272.221 | **360.628** |
| sim n=40, m=8 | 715.745 | 686.520 | **863.531** |
| real Longrich2010 | 0.823 | 0.823 | **1.180** |

Mathematically expected — OLS admits negative branch lengths, so the search chases
topologies that only fit well if lengths may go negative (on one target the OLS winner
reached RSS = 0.00000, an exact signed fit). Not a code defect; the function optimizes
what it is asked to. But the documented purpose is Lapointe & Cucumel average consensus,
which is NNLS/Fitch–Margoliash, so a user switching to `"ols"` for speed silently gets a
topology worse than doing nothing. Worth a sentence in `?LeastSquaresTree`.

**(b) `ts_ls_search` has no certification gate.** Finding 6 above is the LS-side
instance of the effect `dev/profiling/na-certify-gate.md` already documents for
parsimony (`exact_verify_sweep` at 97.7% of `tbr_search` wall), where
`TBRParams::certify_unrooted` is the lever and `ts_sector.cpp` already clears it at 6
sites. `src/ts_ls.{cpp,h}` contain no `certify` concept at all — it is the one search
path without the lever. At n = 40 an NNLS refit costs ≈ 10 ms and `LeastSquaresTree()`
spends ~4 400–7 400 refit-equivalents. If LS search speed ever matters, that is where
to look, not at the starting tree.

## Unexplained flake worth a separate look

`vine_sector_test.R` failed on **12 of 16 sectors** on its first run and **0 of 16** on
an identical-seed re-run, with the same sectors carved both times. Non-deterministic
under repeated same-session `MaximizeParsimony()` calls. The error text was lost to
`try(silent = TRUE)` and was not chased. Flagged only because it is a reproducibility
symptom in a package heading for release; unrelated to anything above.
