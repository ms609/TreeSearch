# Fractional (Nelson–Ladiges) quartet weighting for `QuartetConcordance()`

Branch: `frac-quart`. Goal: port Nelson & Ladiges (1992, *Syst. Biol.* 41:490–494)
fractional weighting of three-item statements into the **unrooted quartet**
concordance measure, so that logically redundant quartets are counted as the
information they actually carry (trits) rather than as raw combinatorial volume.

This file is the settled design. Everything below was derived and, where noted,
verified by brute-force enumeration / GF(2) rank (scratch scripts, this session).

---

## 1. Decisions locked in

- **Unrooted quartets only.** No polarity / no "informative state" asymmetry.
  N&L's rooted `2/n` symmetrises to the quartet weight below.
- **New arg `unit = c("quartet", "trit")`** on `QuartetConcordance()`. It
  *replaces* any `fractional` flag. `"quartet"` = current behaviour (raw counts);
  `"trit"` = N&L fractional (independent-quartet) currency.
- **Normalization = coverage by the *reported* unit** (option "b"): only a
  character-split **identical** to the tree split scores full marks; nesting and
  crossing get partial credit. See §4.
- **`weight` (TRUE/FALSE) is retained and orthogonal to `unit`.** It controls the
  *pooling amount*, not the per-pair quality. See §5.
- **`return` (edge/char/tree) retained;** it selects which unit is "reported" and
  therefore which weight normalises. See §4–5.
- This slots into the package's existing **quality × amount** decomposition
  (`ConcordanceTable()` / `QACol()`): quality = per-pair coverage ratio, amount =
  N&L fractional weight (trit content).

---

## 2. Currency: trits, quartets = 4-cycles

For `t` taxa, a binary character (no missing data) is a bipartition of sizes
`(k, t−k)`. The quartets it resolves are the **4-cycles of the complete bipartite
graph `K_{k, t−k}`** (two taxa each side = one 4-cycle). Nelson–Ladiges
entailment `(ab,cd) + (ab,ce) → (ab,de)` is exactly GF(2) cycle addition, so the
independent-quartet count is the cyclomatic number:

```
independent quartets of a split of sizes (k, t−k)  =  (k−1)(t−k−1)      # "trits"
raw quartets it resolves                           =  C(k,2)·C(t−k,2)
per-quartet fractional weight                      =  4 / (k(t−k))       # = indep/raw
```

`(k−1)(t−k−1)` is the "fractional weight" `W` of the whole character/split. It is
`0` for a constant (`k=0`) or autapomorphy (`k=1`) — so **uninformative characters
self-zero in trit currency**, no special-casing needed.

Smallest N&L case (t=5, k=2): `4/(2·3)·3 = 2` trits from 3 raw quartets. ✓

**Multistate (solved, not deferred).** A quartet is decisive only when two taxa
share one state and two share another, so every decisive quartet lives in the
`K_{n_i, n_j}` edge-block between two states `i, j`. Those blocks are
edge-disjoint, and the N&L entailment never crosses them (it forces the shared
`c,d,e` taxa into one state), so **trits add over state pairs**:

```
W_char = Σ_{i<j} (n_i−1)(n_j−1)      A = Σ_{i<j} A_ij      weight_{ij} = 4/(n_i n_j)
```

reducing to `(n−1)(t−n−1)` for binary (one pair). Each state-pair is an
independent binary sub-problem, so the binary machinery below is applied
per state-pair and summed; a binary character is the single-pair case. This is
the *same* trit currency (multistate and binary are directly comparable), not a
separate unit. Verified by GF(2) rank of the decisive **and** concordant quartet
sets on random 2/3/4-state characters (`dev/benchmarks/frac-quart/multistate_trit.R`).

---

## 3. Per-pair structure: the 2×2 table

Character `c` (bipartition, `n = n_c` taxa in one state) vs tree split `k`
(sizes `m = m_k`, `t−m`). Cross-tabulate into cells:

```
            split A    split B
state I       p          q         (n   = p+q)
state O       r          s         (t−n = r+s)
                                    (m   = p+r,  t−m = q+s)
```

**Kernel decode** (verified against `src/quartet_concordance.cpp`, binary map
state1=I, side1=A ⇒ `n1[1]=p, n0[1]=q, n1[0]=r, n0[0]=s`):

```
raw concordant  conc(k,c) = C(p,2)C(s,2) + C(q,2)C(r,2)
raw decisive    dec(k,c)  = conc + p·q·r·s              # = concordant + discordant (both-resolve)
```

**Trit (independent) quantities** — GF(2) rank of the 4-cycle sets, verified by
enumeration:

```
concordant trits   A(k,c) = (p−1)₊(s−1)₊ + (q−1)₊(r−1)₊      # ₊ = max(·,0); floors are load-bearing
```

The floors matter: without them the identity `ind_dec − ind_conc = pr + qs − 1`
pushes self-agreement above 1.

**Direction-dependence of the decisive count** (verified): the decisive set's
trit rank differs by which partition supplies the pairing, because discordant
quartets are paired differently by character vs split —
`Δ_char = (n−1)(t−n−1)`, `Δ_split = (m−1)(t−m−1)` (generic, all four cells ≥ 1).
`A` is symmetric (shared). This is *why* normalization depends on `return`.

**Three regimes by number of empty cells** (standard split-compatibility trichotomy):

| empty cells | relationship | discordant | verdict |
|---|---|---|---|
| 2 | **identical** (char split = tree split) | none | full support |
| 1 | **nested** (compatible, proper subset) | none | partial support (real!) |
| 0 | **crossing** (incompatible) | `pqrs > 0` | genuine test; conflict lowers score |

Nesting is a *passed test*, not homoplasy: no evidence against ⇒ evidence for
(cf. viviparity → Mammalia despite the platypus). It must get partial positive
credit — not exclusion.

---

## 4. The measure (option b): coverage by the reported unit

Per pair, quality = concordant trits over the **reported unit's** own weight:

```
edge  (report split k):    Q(k,c) = A(k,c) / W_k,   W_k = (m−1)(t−m−1)
char  (report char  c):    Q(k,c) = A(k,c) / W_c,   W_c = (n−1)(t−n−1)
```

`A ≤ W_reported` always (⇒ `Q ∈ [0,1]`), and **per state-pair**
`Q_ij = 1 ⟺ A_ij = W_reported ⟺ identical sub-bipartition`. So for a **binary**
character, only an identical character-split scores full marks; nesting and
crossing are partial, with conflict strictly lowering the score. Verified ladder
(t=8, edge):

| (p,q,r,s) | relationship | A | W_k | Q |
|---|---|---|---|---|
| 4,0,0,4 | identical | 9 | 9 | 1.00 |
| 3,0,2,3 | nested | 4 | 8 | 0.50 |
| 3,1,1,3 | mild crossing | 4 | 9 | 0.44 |
| 2,2,2,2 | max crossing | 2 | 9 | 0.22 |

**Multistate: `Q = 1 ⟺ the split is displayed by the character`** (§5 pools the
per-pair `Q_ij` by `M = min(W_c, W_k)`). When every state block lies wholly on
one side of the split, each split-relevant pair scores `Q_ij = 1` and the one
split-orthogonal pair has `W_k^{ij} = 0 ⇒ M = 0` and drops out, so the pooled
`Q = 1` — even though `A_tot < W_c,tot` and the multistate character is *not*
"identical" to the (binary) split. This is the correct generalisation of
"identical" (for a binary character, displayed = identical), and the intended
behaviour: a reproductive character coded `{oviparous | viviparous | ovovivip.}`
that cleanly refines Mammalia fully supports the Mammalia split. A state block
that *straddles* the split makes its pairs crossing (`Q_ij < 1`), lowering the
pooled score. Verified: `t7,t8,t9`-refining char scores `1.0` on every displayed
split, `0.667` on splits that cut a block (see `test-Concordance.R`).

Raw-currency (`unit="quartet"`) analogue: `Q = conc / (raw quartets of the
reported unit)` — same shape, `conc` and `C(m,2)C(t−m,2)` instead of `A`, `W_k`.

---

## 5. Pooling (`weight` acts on the amount, not the quality)

`W_k` is constant across characters for a fixed edge, so `weight` must act on the
**amount** (pooling weight), else it goes inert:

```
edge FQ(k) = Σ_c  M(k,c)·Q(k,c)  /  Σ_c  M(k,c)
char FQ(c) = Σ_k  M(k,c)·Q(k,c)  /  Σ_k  M(k,c)
```

- `M(k,c)` = shared information = **`min(W_c, W_k)`** (the hBest analogue used by
  `ClusteringConcordance`). Preserves "only identical = full marks" (ceiling needs
  `A = W_c = W_k`).
- `weight = FALSE` ⇒ `M ≡ 1` (per-pair-equal mean of `Q`).
- `unit = "quartet"` ⇒ `M` = raw shared-quartet count; `Q` uses raw conc.

**Quality × amount map to the existing table:** `Q` is `QACol`'s *quality* axis;
`W_c` (trit content) is the *amount* axis. The trit port reuses
`ConcordanceTable()` rather than parallelling it.

---

## 6. Build plan

1. **R prototype first** (binary, no missing). Compute cells `{p,q,r,s}` per
   (split, char) from the split logical vectors + character integer vectors;
   derive `A`, `W_c`, `W_k`, `Q`, `M`; pool per §5. Keep the existing C++ kernel
   for the raw path; the trit path can start in R since cells are a cheap crosstab.
2. **Validate**: (i) brute-force check `A` and the ladder on toy cells;
   (ii) run on `congreveLamsdellMatrices` — report edge-score movement
   `quartet → trit`, the identical/nested/crossing census per edge, and a
   quality/amount table; (iii) confirm `unit="quartet"` reproduces current
   `QuartetConcordance` bit-for-bit.
3. **Kernel** (later, perf only): the R per-state-pair path is already correct
   and validated for binary + multistate + missing data, so a C++ port of `A`
   is now a *performance* optimisation (large matrices), not a correctness
   prerequisite.
4. **Tests**: extend `tests/testthat/test-QuartetConcordance*` — `unit` switch,
   identical→1, nested partial, crossing < nested, uninformative→dropped,
   `unit="quartet"` == legacy.

## 7. Deferred / open

- **Multistate**: ~~deferred~~ **DONE.** Trits add over state pairs
  (`W_char = Σ_{i<j}(n_i−1)(n_j−1)`, `A = Σ A_ij`, per-pair weight `4/(n_i n_j)`);
  see §2. Implemented as a per-state-pair loop in `.TritConcordance()`
  (`R/Concordance.R`); GF(2)-verified; binary is the single-pair special case.
- **Missing data**: ~~deferred~~ handled by construction — the per-state-pair
  cell counts restrict to scored taxa (`scored <- !is.na(col)`), giving each
  character its own effective `t`, `W_c` and per-pair `W_k`. Validated on a
  toy with `?`/`-`/ambiguous tokens (`validate_trit.R`), but the validation
  corpus (congreveLamsdell) has no missing data, so the semantics are proven
  self-consistent (`Q∈[0,1]`, pkg==independent ref) rather than benchmarked.
- **Random-expectation baseline**: ~~deferred~~ **DONE** (`unit = "trit"`).
  New arg `normalize = FALSE / TRUE / <int>` on `QuartetConcordance()`, mirroring
  `ClusteringConcordance()`. The crossing floor (≈ 0.22 in §4) is not zero, so
  without correction a maximally conflicting character still scores positive and
  the floor is split-size dependent; `normalize` re-zeros against the
  concordance *expected by chance*, exactly as `ClusteringConcordance()` re-zeros
  MI against `miRand`.
  - **Null model (chosen): fixed-marginal randomization** — reassign each
    character's tokens across the leaves at random, holding its state counts and
    the split sizes fixed (the multivariate-hypergeometric confusion table; the
    same null as `ClusteringConcordance`'s `miRand` and `Consistency`'s
    `ExpectedLength` tip-shuffle). This is the *only* coherent null for the trit
    currency: it re-zeros each `(split, char, state-pair)` against its own
    combinatorial floor and extends unchanged to multistate and missing data.
  - **Rejected alternative: the flat Minh-style `1/3`** (three quartet
    resolutions ⇒ expected concordance `1/3`). It is a *raw-quartet* heuristic
    that ignores marginals and does **not** extend to floored trit counts (the
    trit floor is `0.22`, not `1/3`), so it cannot re-zero the trit measure
    coherently. Its only remaining use would be re-baselining the *raw* `quartet`
    unit — see the raw-unit note below.
  - **Two estimators** (mirroring `ClusteringConcordance`): `normalize = TRUE`
    computes the **exact** expected pools from the hypergeometric pmf (a small
    double sum, cached on `(n_i, n_j, M, t)`); `normalize = <int>` estimates them
    by **Monte-Carlo** tip-shuffle. Validated to agree within MC error
    (`dev/benchmarks/frac-quart/chance_baseline.R` at the cell level;
    `test-Concordance.R` end-to-end).
  - **What varies under the null.** Only `A` (and, for *multistate* pairs where a
    pair's side-A count `mA` is itself random, `wk` and `M = min(w_c, w_k)`) vary;
    for *binary* characters `mA = M`, `tP = t` so `wk` is fixed and only `A`
    varies (the clean `ClusteringConcordance` case). The exact estimator therefore
    accumulates `E[m]`, `E[m·A/w_k]`, `E[m·A/w_c]` and re-zeros the *pooled*
    edge/char ratio against the pooled expectation, with the **same** `weight`
    aggregation as the observed score.
  - **Guarantees / guards.** `.Rezero(1, z) = 1`, so a **displayed** split still
    scores exactly `1` (the ceiling invariant of §4 survives chance correction —
    only the floor is lifted); conflicting characters go **negative** ("below
    random", as in `ClusteringConcordance`). Values are returned **unclamped**
    (may fall below `−1`); clamp to `[−1, 1]` only at plot time (`QCol`/`QACol`),
    never inside the measure. The `z → 1` blow-up is guarded (→ `NA`).
    `normalize` defaults to `FALSE`, so the published measure is unchanged unless
    the user opts in.
  - **Raw `quartet` unit: DONE.** `normalize` now works for `unit = "quartet"`
    too, using the *same* fixed-marginal null. Per state-pair `conc` and `dec`
    are polynomials in the cells (no floors), so the exact expectation
    `E[conc]`, `E[dec]` is a clean sum over the same trivariate-hypergeometric
    pmf as `.ExpectedTrit` (`.ExpectedQuartet` / `.QuartetExpect`); MC via
    `.QuartetMC` reshuffles tokens and re-scores through the kernel. The pooled
    `conc/dec` ratio is re-zeroed against `E[conc]/E[dec]` (weighted) or with
    the same cell-matching as trits (unweighted). Validated exact≈MC across
    edge/char × weight × binary/multistate (`test-Concordance.R`), and
    `normalize = FALSE` is byte-identical to the published raw measure. The flat
    `1/3` baseline stays rejected (it does not adapt to marginals; the
    fixed-marginal null generalises it).
- **`M = min(W_c,W_k)`** (settled): the pooling amount is the *shared*
  information, so it must be symmetric between character and split (both
  `return`s pool by the same `M`); `M = W_k` would make the char return weight a
  pair by the *edge's* content, which is incoherent as shared information.
  `min(W_c,W_k)` also drops into `ConcordanceTable`'s `hBest·n` slot (§5).
