# TreeSearch

[![Build Status](https://travis-ci.org/ms609/TreeSearch.svg?branch=master)](https://travis-ci.org/ms609/TreeSearch)
[![codecov](https://codecov.io/gh/ms609/TreeSearch/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeSearch)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeSearch)](https://cran.r-project.org/package=TreeSearch)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeSearch)](https://cran.r-project.org/package=TreeSearch)
[![DOI](https://zenodo.org/badge/98171642.svg)](https://zenodo.org/badge/latestdoi/98171642)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

This package implements maximum parsimony tree search in R, 
using the Brazeau, Guillerme & Smith (2019) approach to inapplicable data,
as implemented in the Morphy C library (Brazeau _et al_. 2017).
Implied weighting, profile parsimony and heuristic searches are implemented;
a graphical user interface is available.

# Installation

<!--Install and load the stable version from CRAN, and launch the GUI, as follows:-->
Install and load the stable version from CRAN as follows:
```r
install.packages('TreeSearch')
library('TreeSearch')
```
<!--EasyTrees()-->

The development release incorporates a major reworking of the search interface,
allowing faster and more intuitive tree search.

Install the development version and launch the GUI with:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeSearch')
library("TreeSearch")
EasyTrees()
```

# Quick start

Launch a graphical user interface (development version only)
by typing `TreeSearch::EasyTrees()` in the R console.

For more control over search settings, see [`?MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.html).

# Optimality criteria

'TreeSearch' allows the implementation of various optimality criteria, including

- Fitch parsimony with inapplicable data (Brazeau, Guillerme and Smith, 2019).
- Implied weighting (Goloboff 1993)
- The Profile Parsimony approach introduced by Faith and Trueman (2001).
- Successive Approximations weighting (Farris 1969).

It is also possible to specify bespoke optimality criteria.


# Documentation

- [Package functions](https://ms609.github.io/TreeSearch/reference) reference manual
- [Getting started](https://ms609.github.io/TreeSearch/articles/getting-started.html)
- [Loading phylogenetic data into R](https://ms609.github.io/TreeTools/articles/load-data.html)
- [Parsimony search with inapplicable data](https://ms609.github.io/TreeSearch/articles/inapplicable.html)

- [Calculating concavity profiles](https://ms609.github.io/TreeSearch/articles/profile-scores.html) for Profile Parsimony
- [Tree search with profile parsimony](https://ms609.github.io/TreeSearch/articles/profile.html)

Please note that the 'TreeSearch' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

# References

Brazeau M. D., Smith M. R., Guillerme T. (2017). 
  MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability.
  doi: [10.5281/zenodo.815372](https://doi.org/10.5281/zenodo.815372).

Brazeau, M. D., Guillerme, T. and Smith, M. R. (2019).
  An algorithm for morphological phylogenetic analysis with inapplicable data. 
  Systematic Biology, 68(4), 619-631.
  doi: [10.1093/sysbio/syy083](https://dx.doi.org/10.1093/sysbio/syy083).

Faith D. P., Trueman J. W. H. (2001).
  Towards an inclusive philosophy for phylogenetic inference.
  Systematic Biology, 50(3), 331–350. 
  doi: [10.1080/10635150118627](https://doi.org/10.1080/10635150118627).

Farris, J. S. (1969). A successive approximations approach to character weighting. 
  Systematic Biology, 18(4), 374–385.
  doi: [10.2307/2412182](https://dx.doi.org/10.2307/2412182).

Goloboff, P. A. (1993).
  Estimating character weights during tree search.
  Cladistics, 9(1), 83–91.
  doi: [10.1111/j.1096-0031.1993.tb00209.x](https://doi.org/10.1111/j.1096-0031.1993.tb00209.x).

Goloboff, P. A., Torres, A., Arias, J. S. (2018).
  Weighted parsimony outperforms other methods of phylogenetic inference under 
  models appropriate for morphology.
  Cladistics, 34(4), 407–437. doi: [10.1111/cla.12205](https://doi.org/10.1111/cla.12205).
