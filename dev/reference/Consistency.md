# Consistency and retention "indices"

`Consistency()` calculates the consistency "index" and retention index
(Farris 1989) for each character in a dataset, given a bifurcating tree.
Although there is not a straightforward interpretation of these indices,
they are sometimes taken as an indicator of the fit of a character to a
tree. Values correlate with the number of species sampled and the
distribution of taxa between character states, so are not strictly
comparable between characters in which these factors differ; and values
cannot be compared between datasets (Speed and Arbuckle 2017) .

## Usage

``` r
Consistency(dataset, tree, nRelabel = 0, compress = FALSE)
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- nRelabel:

  Integer specifying how many times to relabel leaves when computing
  MCMC estimate of null tree length for RHI calculation. Steell et
  al. (2025) recommend 1000, but suggest that 100 may suffice. If zero
  (the default), the RHI is not calculated.

- compress:

  Logical specifying whether to retain the compression of a `phyDat`
  object or to return a vector specifying to each individual character,
  decompressed using the dataset's `index` attribute.

## Value

`Consistency()` returns a matrix with named columns specifying the
consistency index (`ci`), retention index (`ri`), rescaled consistency
index (`rc`) and relative homoplasy index (`rhi`).

## Details

The **consistency "index"** (Kluge and Farris 1969) is defined as the
number of steps observed in the most parsimonious mapping of a character
to a tree, divided by the number of steps observed on the shortest
possible tree for that character. A value of one indicates that a
character's fit to the tree is optimal. Note that as the possible values
of the consistency index do not range from zero to one, it is not an
index in the mathematical sense of the term. Shortcomings of this
measure are widely documented (Archie 1989-09; Brooks et al. 1986;
Steell et al. 2025) .

The maximum length of a character (see
[`MaximumLength()`](https://ms609.github.io/TreeSearch/dev/reference/MinimumLength.md))
is the number of steps in a parsimonious reconstruction on the longest
possible tree for a character. The **retention index** is the maximum
length of a character minus the number of steps observed on a given
tree; divided by the maximum length minus the minimum length. It is
interpreted as the ratio between the observed homoplasy, and the maximum
observed homoplasy, and scales from zero (worst fit that can be
reconstructed under parsimony) to one (perfect fit).

The **rescaled consistency index** is the product of the consistency and
retention indices; it rescales the consistency index such that its range
of possible values runs from zero (least consistent) to one (perfectly
consistent).

The **relative homoplasy index** (Steell et al. 2025) is the ratio of
the observed excess tree length to the excess tree length due to chance,
taken as the median score of a character when the leaves of the given
tree are randomly shuffled.

The lengths of characters including inapplicable tokens are calculated
following Brazeau et al. (2019) , matching their default treatment in
[`TreeLength()`](https://ms609.github.io/TreeSearch/dev/reference/TreeLength.md).

## References

Archie JW (1989-09). “Homoplasy Excess Ratios: New Indices for Measuring
Levels of Homoplasy in Phylogenetic Systematics and a Critique of the
Consistency Index.” *Systematic Zoology*, **38**(3), 253.
[doi:10.2307/2992286](https://doi.org/10.2307/2992286) .  
  
Brazeau MD, Guillerme T, Smith MR (2019). “An algorithm for
morphological phylogenetic analysis with inapplicable data.” *Systematic
Biology*, **68**(4), 619–631.
[doi:10.1093/sysbio/syy083](https://doi.org/10.1093/sysbio/syy083) .  
  
Brooks DR, O'Grady RT, Wiley EO (1986). “A Measure of the Information
Content of Phylogenetic Trees, and Its Use as an Optimality Criterion.”
*Systematic Biology*, **35**(4), 571–581.
[doi:10.2307/2413116](https://doi.org/10.2307/2413116) .  
  
Farris JS (1989). “The Retention Index and the Rescaled Consistency
Index.” *Cladistics*, **5**(4), 417–419.
[doi:10.1111/j.1096-0031.1989.tb00573.x](https://doi.org/10.1111/j.1096-0031.1989.tb00573.x)
.  
  
Kluge AG, Farris JS (1969). “Quantitative Phyletics and the Evolution of
Anurans.” *Systematic Zoology*, **18**(1), 1–32.
[doi:10.1093/sysbio/18.1.1](https://doi.org/10.1093/sysbio/18.1.1) .  
  
Speed MP, Arbuckle K (2017). “Quantification Provides a Conceptual Basis
for Convergent Evolution.” *Biological Reviews*, **92**(2), 815–829.
[doi:10.1111/brv.12257](https://doi.org/10.1111/brv.12257) .  
  
Steell EM, Hsiang AY, Field DJ (2025). “Revealing Patterns of Homoplasy
in Discrete Phylogenetic Datasets with a Cross-Comparable Index.”
*Zoological Journal of the Linnean Society*, **204**(1), zlaf024.
[doi:10.1093/zoolinnean/zlaf024](https://doi.org/10.1093/zoolinnean/zlaf024)
.

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data(inapplicable.datasets)
dataset <- inapplicable.phyData[[4]]
head(Consistency(dataset, TreeTools::NJTree(dataset), nRelabel = 10))
#>             ci        ri        rc       rhi
#> [1,] 0.2500000 0.6250000 0.1562500 0.6000000
#> [2,] 0.3333333 0.3333333 0.1111111 0.6666667
#> [3,] 0.3333333 0.6000000 0.2000000 0.5000000
#> [4,] 0.2500000 0.2500000 0.0625000 1.0000000
#> [5,] 0.5000000 0.8333333 0.4166667 0.2000000
#> [6,] 0.2500000 0.2500000 0.0625000 0.7500000
```
