# Concordance factors

Concordance measures the strength of support that characters in a
dataset present for each split (=edge/branch) in a tree (Minh et al.
2020; Smith 2026) .

## Usage

``` r
ClusteringConcordance(tree, dataset, return = "edge", normalize = TRUE)

MutualClusteringConcordance(tree, dataset)

QuartetConcordance(tree, dataset = NULL, weight = TRUE, return = "edge")

PhylogeneticConcordance(tree, dataset)

SharedPhylogeneticConcordance(tree, dataset)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- return:

  Character specifying the summary to return. Options are:

  - `"edge"`: average concordance of each tree split across all
    characters;

  - `"char"`: concordance of each character averaged over all splits,
    weighted by each character's information content;

  - `"tree"`: an overall tree‑level concordance score;

  - `"all"`: a full array of MI components and normalized values for
    every split–character pair.

  Matching is case‑insensitive and partial.

- normalize:

  Controls how the *expected* mutual information (the zero point of the
  scale) is determined.

  - `FALSE`: no chance correction; MI is scaled only by its maximum.

  - `TRUE`: subtract the analytical expected MI for random association.

  - `<integer>`: subtract an empirical expected MI estimated from that
    number of random trees.

  In all cases, 1 corresponds to the maximal attainable MI for the pair
  (`hBest`), and 0 corresponds to the chosen expectation.

- weight:

  Logical specifying whether to weight sites according to the number of
  quartets they are decisive for.

## Value

`ClusteringConcordance(return = "all")` returns a 3D array where each
slice corresponds to a character (site), each column to a tree split,
and each row to a different information measure. The `normalized` row
gives the normalized mutual information between each split-character
pair, scaled so that 1.0 corresponds to `hBest` (the theoretical maximum
mutual information, being the minimum of `hSplit` and `hChar`) and 0.0
corresponds to `miRand` (the expected mutual information under random
association). `hSplit` gives the entropy (information content) of each
split's bipartition; `hChar` gives the entropy of each character's state
distribution; `hJoint` gives the joint entropy of the split-character
confusion matrix; `mi` gives the raw mutual information; and `n` records
the number of informative observations. Negative normalized values
indicate observed mutual information below random expectation. `NA` is
returned when `hBest = 0` (no information potential).

`ClusteringConcordance(return = "edge")` returns a vector where each
element corresponds to a split (an edge of the tree) and gives the
normalized mutual information between that split and the character data,
averaged across all characters. When `normalize = TRUE` (default),
values are scaled relative to random expectation; when `FALSE`, raw
mutual information normalized by `hBest` is returned.

`ClusteringConcordance(return = "char")` returns a vector where each
element corresponds to a character (site) and gives the entropy-weighted
average normalized mutual information between that character and all
tree splits. Characters with higher information content receive
proportionally more weight from splits that can potentially convey more
information about them.

`ClusteringConcordance(return = "tree")` returns a single value
representing the overall concordance between the tree topology and the
character data. This averages the fit of the best-matching split for
each character. This is included for completeness, though it is not
clear that this is a useful or meaningful measure.

`MutualClusteringConcordance()` returns the mutual clustering
concordance of each character in `dataset` with `tree`. The attribute
`weighted.mean` gives the mean value, weighted by the information
content of each character.

`QuartetConcordance(return = "edge")` returns a numeric vector giving
the concordance index at each split across all sites; names specify the
number of each corresponding split in `tree`.

`QuartetConcordance(return = "char")` returns a numeric vector giving
the concordance index calculated at each site, averaged across all
splits.

`PhylogeneticConcordance()` returns a numeric vector giving the
phylogenetic information of each split in `tree`, named according to the
split's internal numbering.

`SharedPhylogeneticConcordance()` returns the shared phylogenetic
concordance of each character in `dataset` with `tree`. The attribute
`weighted.mean` gives the mean value, weighted by the information
content of each character.

## Details

`ClusteringConcordance()` measures how well each tree split reflects the
grouping structure implied by each character (Smith 2026) . Characters
and splits are treated as clusterings of taxa, and their agreement is
quantified using mutual information (MI). All reported values are scaled
so that 1 corresponds to the maximum possible mutual information for
each split–character pair (`hBest`).

The `normalize` argument specifies how the zero point is defined.

- If `normalize = FALSE`, zero corresponds to *zero* MI, without
  correcting for the positive bias that arises because MI is rarely
  exactly zero in finite samples.

- If `normalize = TRUE`, the expected MI is computed using an analytical
  approximation based on the distribution of character tokens. This is
  fast and generally accurate for large trees (~200+ taxa), but does not
  account for correlation between splits.

- If `normalize` is a positive integer `n`, the expected MI is estimated
  empirically by fitting each character to `n` uniformly random trees
  and averaging the resulting MI values. This Monte Carlo approach
  provides a more accurate baseline for small trees, for which the
  analytical approximation is biased. Monte Carlo standard errors are
  returned.

`MutualClusteringConcordance()` provides a character‑wise summary that
emphasises each character’s best‑matching split(s). It treats each
character as a simple tree and computes the mutual clustering
information between this character‑tree and the supplied phylogeny. High
values identify characters whose signal is well represented anywhere in
the tree, even if concentrated on a single edge.

`QuartetConcordance()` is the proportion of quartets (sets of four
leaves) that are decisive for a split which are also concordant with it
(the site concordance factor (Minh et al. 2020) ). For example, a
quartet with the characters `0 0 0 1` is not decisive, as all
relationships between those leaves are equally parsimonious. But a
quartet with characters `0 0 1 1` is decisive, and is concordant with
any tree that groups the first two leaves together to the exclusion of
the second.

By default, the reported value weights each site by the number of
quartets it is decisive for. This value can be interpreted as the
proportion of all decisive quartets that are concordant with a split. If
`weight = FALSE`, the reported value is the mean of the concordance
value for each site. Consider a split associated with two sites: one
that is concordant with 25% of 96 decisive quartets, and a second that
is concordant with 75% of 4 decisive quartets. If `weight = TRUE`, the
split concordance will be 24 + 3 / 96 + 4 = 27%. If `weight = FALSE`,
the split concordance will be mean(75%, 25%) = 50%.

`QuartetConcordance()` is computed exactly, using all quartets, where as
other implementations (e.g. IQ-TREE) follow Minh2020) in using a random
subsample of quartets for a faster, if potentially less accurate,
computation. Ambiguous and inapplicable tokens are treated as containing
no grouping information (i.e. `(02)` or `-` are each treated as `?`).

`PhylogeneticConcordance()` treats each character in `dataset` as a
phylogenetic hypothesis and measures the extent to which it supports the
splits of `tree`. Each character is first interpreted as a tree (or set
of trees) in which taxa sharing the same token form a clade. Only splits
for which the character contains at least four relevant taxa can
contribute information.

For each split, the function identifies which characters could
potentially support that split (i.e. those for which the induced
subtrees contain informative structure), and among these, which
characters are actually compatible with the split. The concordance value
for each split is the proportion of informative characters that support
it. A value of 1 indicates that all characters informative for that
subset of taxa support the split; a value of 0 indicates that none do.
Characters that contain only ambiguous or uninformative states for the
relevant taxa do not affect the result.

`SharedPhylogeneticConcordance()` treats each character as a simple
tree. Each token in the character corresponds to a node whose pendant
edges are the taxa with that token. The Shared Phylogenetic Concordance
for each character in `dataset` is then the Shared Phylogenetic
Information (Smith 2020) of this tree and `tree`.

## References

Minh BQ, Hahn MW, Lanfear R (2020). “New methods to calculate
concordance factors for phylogenomic datasets.” *Molecular Biology and
Evolution*, **37**(9), 2727–2733.
[doi:10.1093/molbev/msaa106](https://doi.org/10.1093/molbev/msaa106) .  
  
Smith MR (2020). “Information Theoretic Generalized Robinson-Foulds
Metrics for Comparing Phylogenetic Trees.” *Bioinformatics*, **36**(20),
5007–5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
.  
  
Smith MR (2026). “Which characters support which clades? Exploring the
distribution of phylogenetic signal using concordant information.”
*Forthcoming*.

## See also

- [`Consistency()`](https://ms609.github.io/TreeSearch/reference/Consistency.md)

Other split support functions:
[`ConcordanceTable()`](https://ms609.github.io/TreeSearch/reference/ConcordanceTable.md),
[`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md),
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md),
[`MostContradictedFreq()`](https://ms609.github.io/TreeSearch/reference/MostContradictedFreq.md),
[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples
