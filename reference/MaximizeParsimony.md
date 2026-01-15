# Find most parsimonious trees

Search for most parsimonious trees using the parsimony ratchet and TBR
rearrangements, treating inapplicable data as such using the algorithm
of Brazeau et al. (2019) .

Tree search will be conducted from a specified or
automatically-generated starting tree in order to find a tree with an
optimal parsimony score, under implied or equal weights, treating
inapplicable characters as such in order to avoid the artefacts of the
standard Fitch algorithm (see Maddison 1993; Brazeau et al. 2019) . Tree
length is calculated using the MorphyLib C library (Brazeau et al. 2017)
.

## Usage

``` r
MaximizeParsimony(
  dataset,
  tree,
  ratchIter = 7L,
  tbrIter = 2L,
  startIter = 2L,
  finalIter = 1L,
  maxHits = NTip(dataset) * 1.8,
  maxTime = 60,
  quickHits = 1/3,
  concavity = Inf,
  ratchEW = TRUE,
  tolerance = sqrt(.Machine[["double.eps"]]),
  constraint,
  verbosity = 3L
)

Resample(
  dataset,
  tree,
  method = "jack",
  proportion = 2/3,
  ratchIter = 1L,
  tbrIter = 8L,
  finalIter = 3L,
  maxHits = 12L,
  concavity = Inf,
  tolerance = sqrt(.Machine[["double.eps"]]),
  constraint,
  verbosity = 2L,
  ...
)

EasyTrees()

EasyTreesy()
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

  (optional) A bifurcating tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), containing only
  the tips listed in `dataset`, from which the search should begin. If
  unspecified, an [addition
  tree](https://ms609.github.io/TreeSearch/reference/AdditionTree.md)
  will be generated from `dataset`, respecting any supplied
  `constraint`. Edge lengths are not supported and will be deleted.

- ratchIter:

  Numeric specifying number of iterations of the parsimony ratchet
  (Nixon 1999) to conduct.

- tbrIter:

  Numeric specifying the maximum number of TBR break points on a given
  tree to evaluate before terminating the search. One "iteration"
  comprises selecting a branch to break, and evaluating each possible
  reconnection point in turn until a new tree improves the score. If a
  better score is found, then the counter is reset to zero, and tree
  search continues from the improved tree.

- startIter:

  Numeric: an initial round of tree search with `startIter` × `tbrIter`
  TBR break points is conducted in order to locate a local optimum
  before beginning ratchet searches.

- finalIter:

  Numeric: a final round of tree search will evaluate `finalIter` ×
  `tbrIter` TBR break points, in order to sample the final optimal
  neighbourhood more intensely.

- maxHits:

  Numeric specifying the maximum times that an optimal parsimony score
  may be hit before concluding a ratchet iteration or final search
  concluded.

- maxTime:

  Numeric: after `maxTime` minutes, stop tree search at the next
  opportunity.

- quickHits:

  Numeric: iterations on subsampled datasets will retain `quickHits` ×
  `maxHits` trees with the best score.

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

- ratchEW:

  Logical specifying whether to use equal weighting during ratchet
  iterations, improving search speed whilst still facilitating escape
  from local optima.

- tolerance:

  Numeric specifying degree of suboptimality to tolerate before
  rejecting a tree. The default, `sqrt(.Machine$double.eps)`, retains
  trees that may be equally parsimonious but for rounding errors.
  Setting to larger values will include trees suboptimal by up to
  `tolerance` in search results, which may improve the accuracy of the
  consensus tree (at the expense of resolution) (Smith 2019) .

- constraint:

  Either an object of class `phyDat`, in which case returned trees will
  be perfectly compatible with each character in `constraint`; or a tree
  of class `phylo`, all of whose nodes will occur in any output tree.
  See
  [`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.html)
  and
  [vignette](https://ms609.github.io/TreeSearch/articles/tree-search.html)
  for further examples.

- verbosity:

  Integer specifying level of messaging; higher values give more
  detailed commentary on search progress. Set to `0` to run silently.

- method:

  Unambiguous abbreviation of `jackknife` or `bootstrap` specifying how
  to resample characters. Note that jackknife is considered to give more
  meaningful results.

- proportion:

  Numeric between 0 and 1 specifying what proportion of characters to
  retain under jackknife resampling.

- ...:

  Additional parameters to `MaximizeParsimony()`.

## Value

`MaximizeParsimony()` returns a list of trees with class `multiPhylo`.
This lists all trees found during each search step that are within
`tolerance` of the optimal score, listed in the sequence that they were
first visited, and named according to the step in which they were first
found; it may contain more than `maxHits` elements. Note that the
default search parameters may need to be increased in order for these
trees to be the globally optimal trees; examine the messages printed
during tree search to evaluate whether the optimal score has stabilized.

The return value has the attribute `firstHit`, a named integer vector
listing the number of optimal trees visited for the first time in each
stage of the tree search. Stages are named:

- `seed`: starting trees;

- `start`: Initial TBR search;

- `ratchN`: Ratchet iteration `N`;

- `final`: Final TBR search. The first tree hit for the first time in
  ratchet iteration three is named `ratch3_1`.

`Resample()` returns a `multiPhylo` object containing a list of trees
obtained by tree search using a resampled version of `dataset`.

## Details

Tree search commences with `ratchIter` iterations of the parsimony
ratchet (Nixon 1999) , which bootstraps the input dataset in order to
escape local optima. A final round of tree bisection and reconnection
(TBR) is conducted to broaden the sampling of trees.

This function can be called using the R command line / terminal, or
through the "shiny" graphical user interface app (type `EasyTrees()` to
launch).

The optimal strategy for tree search depends in part on how close to
optimal the starting tree is, the size of the search space (which
increases super-exponentially with the number of leaves), and the
complexity of the search space (e.g. the existence of multiple local
optima).

One possible approach is to employ four phases:

1.  Rapid search for local optimum: tree score is typically easy to
    improve early in a search, because the initial tree is often far
    from optimal. When many moves are likely to be accepted, running
    several rounds of search with a low value of `maxHits` and a high
    value of `tbrIter` allows many trees to be evaluated quickly,
    hopefully moving quickly to a more promising region of tree space.

2.  Identification of local optimum: Once close to a local optimum, a
    more extensive search with a higher value of `maxHits` allows a
    region to be explored in more detail. Setting a high value of
    `tbrIter` will search a local neighbourhood more completely

3.  Search for nearby peaks: Ratchet iterations allow escape from local
    optima. Setting `ratchIter` to a high value searches the wider
    neighbourhood more extensively for other nearby peaks;
    `ratchEW = TRUE` accelerates these exploratory searches. Ratchet
    iterations can be ineffective when `maxHits` is too low for the
    search to escape its initial location.

4.  Extensive search of final optimum. As with step 2, it may be
    valuable to fully explore the optimum that is found after ratchet
    searches to be sure that the locally optimal score has been
    obtained. Setting a high value of `finalIter` performs a thorough
    search that can give confidence that further searches would not find
    better (local) trees.

A search is unlikely to have found a global optimum if:

- Tree score continues to improve on the final iteration. If a local
  optimum has not yet been reached, it is unlikely that a global optimum
  has been reached. Try increasing `maxHits`.

- Successive ratchet iterations continue to improve tree scores. If a
  recent ratchet iteration improved the score, rather than finding a
  different region of tree space with the same optimal score, it is
  likely that still better global optima remain to be found. Try
  increasing `ratchIter` (more iterations give more chance for
  improvement) and `maxHits` (to get closer to the local optimum after
  each ratchet iteration).

- Optimal areas of tree space are only visited by a single ratchet
  iteration. (See vignette: [Exploring tree
  space](https://ms609.github.io/TreeSearch/articles/tree-space.html).)
  If some areas of tree space are only found by one ratchet iteration,
  there may well be other, better areas that have not yet been visited.
  Try increasing `ratchIter`.

When continuing a tree search, it is usually best to start from an
optimal tree found during the previous iteration - there is no need to
start from scratch.

A more time consuming way of checking that a global optimum has been
reached is to repeat a search with the same parameters multiple times,
starting from a different, entirely random tree each time. If all
searches obtain the same optimal tree score despite their different
starting points, this score is likely to correspond to the global
optimum.

For detailed documentation of the "TreeSearch" package, including full
instructions for loading phylogenetic data into R and initiating and
configuring tree search, see the [package
documentation](https://ms609.github.io/TreeSearch/).

## Resampling

Note that bootstrap support is a measure of the amount of data
supporting a split, rather than the amount of confidence that should be
afforded the grouping. "Bootstrap support of 100% is not enough, the
tree must also be correct" (Phillips et al. 2004) . See discussion in
Egan (2006) ; Wagele et al. (2009) ; (Simmons and Freudenstein 2011) ;
Kumar et al. (2012) .

For a discussion of suitable search parameters in resampling estimates,
see Muller (2005) . The user should decide whether to start each
resampling from the optimal tree (which may be quicker, but result in
overestimated support values as searches get stuck in local optima close
to the optimal tree) or a random tree (which may take longer as more
rearrangements are necessary to find an optimal tree on each iteration).

For other ways to estimate clade concordance, see
[`SiteConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md).

## References

Brazeau MD, Guillerme T, Smith MR (2019). “An algorithm for
morphological phylogenetic analysis with inapplicable data.” *Systematic
Biology*, **68**(4), 619–631.
[doi:10.1093/sysbio/syy083](https://doi.org/10.1093/sysbio/syy083) .  
  
Brazeau MD, Smith MR, Guillerme T (2017). “MorphyLib: a library for
phylogenetic analysis of categorical trait data with inapplicability.”
[doi:10.5281/zenodo.815372](https://doi.org/10.5281/zenodo.815372) .  
  
Egan MG (2006). “Support versus corroboration.” *Journal of Biomedical
Informatics*, **39**(1), 72–85.
[doi:10.1016/j.jbi.2005.11.007](https://doi.org/10.1016/j.jbi.2005.11.007)
.  
  
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
  
Kumar S, Filipski AJ, Battistuzzi FU, Kosakovsky Pond SL, Tamura K
(2012). “Statistics and truth in phylogenomics.” *Molecular Biology and
Evolution*, **29**(2), 457–472.
[doi:10.1093/molbev/msr202](https://doi.org/10.1093/molbev/msr202) .  
  
Maddison WP (1993). “Missing data versus missing characters in
phylogenetic analysis.” *Systematic Biology*, **42**(4), 576–581.
[doi:10.1093/sysbio/42.4.576](https://doi.org/10.1093/sysbio/42.4.576)
.  
  
Muller KF (2005). “The efficiency of different search strategies in
estimating parsimony jackknife, bootstrap, and Bremer support.” *BMC
Evolutionary Biology*, **5**(1), 58.
[doi:10.1186/1471-2148-5-58](https://doi.org/10.1186/1471-2148-5-58) .  
  
Nixon KC (1999). “The Parsimony Ratchet, a new method for rapid
parsimony analysis.” *Cladistics*, **15**(4), 407–414. ISSN 0748-3007,
[doi:10.1111/j.1096-0031.1999.tb00277.x](https://doi.org/10.1111/j.1096-0031.1999.tb00277.x)
.  
  
Phillips MJ, Delsuc F, Penny D (2004). “Genome-scale phylogeny and the
detection of systematic biases.” *Molecular biology and evolution*,
**21**(7), 1455–8.
[doi:10.1093/molbev/msh137](https://doi.org/10.1093/molbev/msh137) .  
  
Simmons MP, Freudenstein JV (2011). “Spurious 99% bootstrap and
jackknife support for unsupported clades.” *Molecular Phylogenetics and
Evolution*, **61**(1), 177–191.
[doi:10.1016/j.ympev.2011.06.003](https://doi.org/10.1016/j.ympev.2011.06.003)
.  
  
Smith MR (2019). “Bayesian and parsimony approaches reconstruct
informative trees from simulated morphological datasets.” *Biology
Letters*, **15**(2), 20180632.
[doi:10.1098/rsbl.2018.0632](https://doi.org/10.1098/rsbl.2018.0632) .  
  
Wagele JW, Letsch H, Klussmann-Kolb A, Mayer C, Misof B, Wagele H
(2009). “Phylogenetic support values are not necessarily informative:
the case of the Serialia hypothesis (a mollusk phylogeny).” *Frontiers
in Zoology*, **6**(1), 12–29.
[doi:10.1186/1742-9994-6-12](https://doi.org/10.1186/1742-9994-6-12) .

## See also

Tree search *via* graphical user interface: `EasyTrees()`

Other split support functions:
[`ConcordanceTable()`](https://ms609.github.io/TreeSearch/reference/ConcordanceTable.md),
[`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md),
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MostContradictedFreq()`](https://ms609.github.io/TreeSearch/reference/MostContradictedFreq.md),
[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md),
[`SiteConcordance`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
## Only run examples in interactive R sessions
if (interactive()) {
  # launch "shiny" point-and-click interface
  EasyTrees()
  
  # Here too, use the "continue search" function to ensure that tree score
  # has stabilized and a global optimum has been found
}


# Load data for analysis in R
library("TreeTools")
data("inapplicable.phyData", package = "TreeSearch")
dataset <- inapplicable.phyData[["Asher2005"]]

# A very quick run for demonstration purposes
trees <- MaximizeParsimony(dataset, ratchIter = 0, startIter = 0,
                           tbrIter = 1, maxHits = 4, maxTime = 1/100,
                           concavity = 10, verbosity = 4)
#> 
#> ── BEGIN TREE SEARCH (k = 10) ──────────────────────────────────────────────────
#> → Initial score: 16.3264 
#> 
#> ── Sample local optimum ────────────────────────────────────────────────────────
#> → TBR depth 1; keeping 4 trees; k = 10
#> ℹ 2026-01-15 13:59:15: Score: 16.3264
#> ✔ 2026-01-15 13:59:16: Tree search terminated with score 15.969
names(trees)
#> [1] "final_1" "final_2" "final_3"
cons <- Consensus(trees)

# In actual use, be sure to check that the score has converged on a global
# optimum, conducting additional iterations and runs as necessary.
 
if (interactive()) {
# Jackknife resampling
nReplicates <- 10
jackTrees <- replicate(nReplicates,
  #c() ensures that each replicate returns a list of trees
  c(Resample(dataset, trees, ratchIter = 0, tbrIter = 2, startIter = 1,
             maxHits = 5, maxTime = 1 / 10,
             concavity = 10, verbosity = 0))
 )

# In a serious analysis, more replicates would be conducted, and each
# search would undergo more iterations.

# Now we must decide what to do with the multiple optimal trees from
# each replicate.

# Set graphical parameters for plotting
oPar <- par(mar = rep(0, 4), cex = 0.9)

# Take the strict consensus of all trees for each replicate
# (May underestimate support)
JackLabels(cons, lapply(jackTrees, ape::consensus))

# Take a single tree from each replicate (here, the first)
# Potentially problematic if chosen tree is not representative
JackLabels(cons, lapply(jackTrees, `[[`, 1))

# Count iteration as support if all most parsimonious trees support a split;
# as contradiction if all trees contradict it; don't include replicates where
# not all trees agree on the resolution of a split.
labels <- JackLabels(cons, jackTrees)

# How many iterations were decisive for each node?
attr(labels, "decisive")

# Show as proportion of decisive iterations
JackLabels(cons, jackTrees, showFrac = TRUE)

# Restore graphical parameters
par(oPar)
}

# Tree search with a constraint
constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
characters <- MatrixToPhyDat(matrix(
  c(0, 1, 1, 1, 0, 0,
    1, 1, 1, 0, 0, 0), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
MaximizeParsimony(characters, constraint = constraint, verbosity = 0)
#> ✔ Initialized 1 distinct constraints.
#> 1 phylogenetic tree
```
