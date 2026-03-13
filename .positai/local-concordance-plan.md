# LocalConcordance: Entropy-based regional character agreement

## Motivation

PhyIN (Maddison 2024, PeerJ 12:e18504) identifies "garbage" regions in UCE
alignments by detecting phylogenetic incompatibility between neighbouring sites.
`LocalConcordance()` improves on this idea by replacing the binary
compatible/incompatible test with an information-theoretic measure: normalized
mutual information between site pairs, computed using clustering entropy
(as in TreeDist / ClusteringInfoDistance).

## Core idea

For each site pair (s_i, s_j) in an alignment:

1. Compute mutual information (MI) of the two characters treated as clusterings.
2. Normalize to [0, 1] where 0 = expected MI under independence (EMI) and
   1 = maximum MI (= min(H(s_i), H(s_j))).  Reuse `.ExpectedMI()` and
   `.Rezero()` from `Concordance.R`.

Then, for each site s_i, produce a score at each of several spatial scales σ:

- Weight each pair (s_i, s_j) by a Gaussian kernel f(|j − i|; σ), truncated
  at some multiple of σ for efficiency.
- Additionally weight (multiplicatively) by the joint entropy H(s_i, s_j),
  so that more informative pairs contribute more.
- The weighted average of normalized MI gives a per-site, per-scale score.

Multiple σ values (e.g. 5, 10, 25, 50, 100 sites) yield a multi-scale
"spectral" summary of local concordance.

## Implementation plan

### 1. Pairwise MI engine (`R/LocalConcordance.R`)

Write a helper that, given two integer-coded character vectors (with NA for
ambiguity), returns the normalized MI score.  This parallels the inner loop of
`ClusteringConcordance()` (lines ~200–250 of Concordance.R) but compares two
characters instead of a character and a split.

Key reusable internals:
- `entropy_int()` (imported from TreeDist) for H(s_i), H(s_j), H(s_i, s_j)
- `.ExpectedMI()` for EMI (cached, accepts contingency table marginals)
- `.Rezero()` for normalization

Consider implementing the pairwise MI computation in C++ (Rcpp) if profiling
shows the R-level loop is too slow for large alignments.  The inner computation
is simple tabulation + `entropy_int`; a C++ version could process all pairs
within the truncation window efficiently.

### 2. Alignment preprocessing

Reuse the token-cleaning logic from `ClusteringConcordance()`:
- Map the contrast matrix to unambiguous tokens; set ambiguous to NA.
- Drop autapomorphic states (singletons → NA).
- Honour the `index` attribute (so duplicate site patterns share results).

Factor this out into a shared helper if not already extracted.

### 3. Gaussian-weighted aggregation

For each site i and each scale σ:
- Determine the window: j ∈ [i − k, i + k] where k = ceiling(3σ) or similar.
- For each j in the window, look up:
  - `nmi[i, j]`: the normalized MI (from step 1).
  - `hJoint[i, j]`: joint entropy H(s_i, s_j), used as an importance weight.
  - `w[i, j]`: Gaussian weight dnorm(|j − i|, 0, σ).
- Compute: score_σ(i) = Σ_j  w[i,j] · hJoint[i,j] · nmi[i,j]
                        / Σ_j  w[i,j] · hJoint[i,j]

Return a matrix: sites × scales.

### 4. Return value and API

```r
LocalConcordance <- function(
  dataset,
  sigma = c(5, 10, 25, 50, 100),
  truncate = 3,         # truncate Gaussian at truncate * sigma
  normalize = TRUE      # EMI correction, as in ClusteringConcordance
)
```

- No `tree` argument: this is purely character-vs-character.
- `dataset`: a `phyDat` object (as used throughout TreeSearch).
- Returns a matrix (or data frame) with one row per site and one column per σ.
  Possibly also attach the raw pairwise NMI matrix as an attribute.

### 5. Visualization

A plotting function (or a `plot` method) that:

- Draws the alignment coloured by nucleotide / token (bottom panel).
- Above it, draws one heatmap row per σ value, colouring each site by its
  LocalConcordance score at that scale.
- Colour scale: diverging or sequential; low concordance (potential garbage)
  in red, high in blue/green.

This could use base graphics or ggplot2.  Base graphics is strongly preferred
as ggplot2 is not currently used in the package, and is not to the user's taste.

This should be a standalone plotting function, but may make use of
the existing `ConcordanceTable()` infrastructure (e.q. QACol())

### 6. Documentation and testing

- Roxygen documentation with `@examples` using `congreveLamsdellMatrices`.
- Unit tests (`tests/testthat/`):
  - Known MI values for hand-constructed two-character examples.
  - Gaussian weighting produces expected results for uniform-MI input.
  - Edge cases: all-ambiguous columns, single-state columns, very short
    alignments.
- A vignette or extended example showing application to a UCE-style dataset
  would strengthen the connection to PhyIN, but could come later.

### 7. Performance considerations

- For an alignment of L sites and window size k, the pairwise MI computation
  is O(L · k · n) where n = number of taxa.  For L = 5000, k = 300,
  n = 100, this is ~1.5 × 10^8 tabulations — likely needs C++ or at
  minimum vectorized R.
- Storing the full L × L NMI matrix is unnecessary; only pairs within the
  largest window are needed.  A banded storage approach (or just computing
  on the fly per σ) avoids O(L²) memory.
- The `index` attribute of phyDat means many columns share the same pattern;
  pairwise MI for pattern pairs can be computed once and looked up by index.

## File locations

| What                    | Where                                     |
|-------------------------|-------------------------------------------|
| Main implementation     | `R/LocalConcordance.R`                    |
| Plotting                | `R/LocalConcordance.R` or `R/plot_lc.R`  |
| Tests                   | `tests/testthat/test-LocalConcordance.R`  |
| C++ (if needed)         | `src/local_concordance.cpp`               |

## Open questions

- **Weighting by joint vs. marginal entropy**: The plan weights by H(s_i, s_j).
  An alternative is min(H(s_i), H(s_j)) (the normalization denominator), which
  weights by potential information rather than observed joint complexity.
  Worth trying both.
  > User comment: Yes, min(H(s_i), H(s_j)) is the natural denominator;
  > not clear that there's a strong argument for using the joint entropy,
  > but if so, let's try it?

- **Self-concordance**: What is score(i, i)?  By construction NMI = 1.
  Excluding self from the weighted average seems right (it's uninformative
  about regional agreement).  But it could serve as a check.
  > Yes, score(i, i) == 1 (always), and should be excluded.

- **Handling gaps**: Current ClusteringConcordance treats gaps as ambiguity.
  For molecular alignments, gap-rich regions are themselves a signal of
  poor alignment.  Consider an option to treat gaps as a distinct state, or
  to report gap frequency alongside concordance.
  > We treat as NA; Maddison requires >50% site occupancy but I don't know
  > that this is relevant if we're weighting by information content

- **Normalization baseline for site pairs**: `.ExpectedMI()` currently
  estimates EMI analytically from marginal counts. This is designed for
  a split (binary partition) vs. a character.  For two multi-state characters,
  the same analytical formula applies but should be verified for accuracy
  with the token counts typical of nucleotide data (4 states + gaps).
  > Agreed.
