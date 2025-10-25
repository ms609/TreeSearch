# Calculate the parsimony score of a tree given a dataset

`TreeLength()` uses the Morphy library (Brazeau et al. 2017) to
calculate a parsimony score for a tree, handling inapplicable data
according to the algorithm of Brazeau et al. (2019) . Trees may be
scored using equal weights, implied weights (Goloboff 1993) , or profile
parsimony (Faith and Trueman 2001) .

## Usage

``` r
IWScore(tree, dataset, concavity = 10L, ...)

TreeLength(tree, dataset, concavity = Inf)

# S3 method for class 'phylo'
TreeLength(tree, dataset, concavity = Inf)

# S3 method for class 'numeric'
TreeLength(tree, dataset, concavity = Inf)

# S3 method for class 'list'
TreeLength(tree, dataset, concavity = Inf)

# S3 method for class 'multiPhylo'
TreeLength(tree, dataset, concavity = Inf)

Fitch(tree, dataset)
```

## Arguments

- tree:

  A tree of class `phylo`, a list thereof (optionally of class
  `multiPhylo`), or an integer – in which case `tree` random trees will
  be uniformly sampled.

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- concavity:

  Determines the degree to which extra steps beyond the first are
  penalized. Specify a numeric value to use implied weighting
  (Goloboff 1993) ; `concavity` specifies *k* in *k* / *e* + *k*. A
  value of 10 is recommended; TNT sets a default of 3, but this is too
  low in some circumstances (Goloboff et al. 2018; Smith 2019) . Better
  still explore the sensitivity of results under a range of concavity
  values, e.g. `k = 2 ^ (1:7)`. Specify `Inf` to weight each additional
  step equally, (which underperforms step weighting approaches (Goloboff
  et al. 2008; Goloboff et al. 2018; Goloboff and Arias 2019;
  Smith 2019) ). Specify `"profile"` to employ an approximation of
  profile parsimony (Faith and Trueman 2001) .

- ...:

  unused; allows additional parameters specified within ... to be
  received by the function without throwing an error.

## Value

`TreeLength()` returns a numeric vector containing the score for each
tree in `tree`.

## References

Brazeau MD, Guillerme T, Smith MR (2019). “An algorithm for
morphological phylogenetic analysis with inapplicable data.” *Systematic
Biology*, **68**(4), 619–631.
[doi:10.1093/sysbio/syy083](https://doi.org/10.1093/sysbio/syy083) .  
  
Brazeau MD, Smith MR, Guillerme T (2017). “MorphyLib: a library for
phylogenetic analysis of categorical trait data with inapplicability.”
[doi:10.5281/zenodo.815372](https://doi.org/10.5281/zenodo.815372) .  
  
Faith DP, Trueman JWH (2001). “Towards an inclusive philosophy for
phylogenetic inference.” *Systematic Biology*, **50**(3), 331–350.
[doi:10.1080/10635150118627](https://doi.org/10.1080/10635150118627) .  
  
Goloboff PA (1993). “Estimating character weights during tree search.”
*Cladistics*, **9**(1), 83–91.
[doi:10.1111/j.1096-0031.1993.tb00209.x](https://doi.org/10.1111/j.1096-0031.1993.tb00209.x)
.  
  
Goloboff PA, Arias JS (2019). “Likelihood approximations of implied
weights parsimony can be selected over the Mk model by the Akaike
information criterion.” *Cladistics*, **35**(6), 695–716.
[doi:10.1111/cla.12380](https://doi.org/10.1111/cla.12380) .  
  
Goloboff PA, Carpenter JM, Arias JS, Esquivel DRM (2008). “Weighting
against homoplasy improves phylogenetic analysis of morphological data
sets.” *Cladistics*, **24**(5), 758–773.
[doi:10.1111/j.1096-0031.2008.00209.x](https://doi.org/10.1111/j.1096-0031.2008.00209.x)
.  
  
Goloboff PA, Torres A, Arias JS (2018). “Weighted parsimony outperforms
other methods of phylogenetic inference under models appropriate for
morphology.” *Cladistics*, **34**(4), 407–437.
[doi:10.1111/cla.12205](https://doi.org/10.1111/cla.12205) .  
  
Smith MR (2019). “Bayesian and parsimony approaches reconstruct
informative trees from simulated morphological datasets.” *Biology
Letters*, **15**(2), 20180632.
[doi:10.1098/rsbl.2018.0632](https://doi.org/10.1098/rsbl.2018.0632) .

## See also

- Conduct tree search using
  [`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  (command line),
  [`EasyTrees()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  (graphical user interface), or
  [`TreeSearch()`](https://ms609.github.io/TreeSearch/reference/TreeSearch.md)
  (custom optimality criteria).

- See score for each character:
  [`CharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md).

Other tree scoring:
[`CharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md),
[`ExpectedLength()`](https://ms609.github.io/TreeSearch/reference/ExpectedLength.md),
[`LengthAdded()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md),
[`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md),
[`MorphyTreeLength()`](https://ms609.github.io/TreeSearch/reference/MorphyTreeLength.md),
[`TaxonInfluence()`](https://ms609.github.io/TreeSearch/reference/TaxonInfluence.md)

## Author

Martin R. Smith (using Morphy C library, by Martin Brazeau)

## Examples

``` r
data("inapplicable.datasets")
tree <- TreeTools::BalancedTree(inapplicable.phyData[[1]])
TreeLength(tree, inapplicable.phyData[[1]])
#> [1] 1117
TreeLength(tree, inapplicable.phyData[[1]], concavity = 10)
#> [1] 52.75785
TreeLength(tree, inapplicable.phyData[[1]], concavity = "profile")
#> → Inapplicable tokens treated as ambiguous for profile parsimony
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> Warning: Can handle max. 2 informative tokens. Dropping others.
#> [1] 3941.387
TreeLength(5, inapplicable.phyData[[1]])
#> [1] 1960 1927 1967 1908 1979
```
