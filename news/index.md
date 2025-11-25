# Changelog

## TreeSearch 1.7.0.9000 (development)

- Fix regression in
  [`MaximumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md).

## TreeSearch 1.7.0 (2025-08-22)

CRAN release: 2025-08-23

- [`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md)
  implements the Group Present / Contradicted measure of Goloboff et
  al. (2003).
- [`Consistency()`](https://ms609.github.io/TreeSearch/reference/Consistency.md)
  also returns the relative homoplasy index of Steell et al.  (2023).
- [`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md)
  supports multiple trees per iteration
  ([\#197](https://github.com/ms609/TreeSearch/issues/197))\[<https://github.com/ms609/TreeSearch/discussions/197>\]
- Support single-character matrices in
  [`ClusteringConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
- Fix `DoNothing(x)` to return `x` (not `NULL`)
- Remove unused `delete_rawdata()` due to implementation issues.

## TreeSearch 1.6.1 (2025-06-10)

CRAN release: 2025-06-10

- Handle invariant characters in
  [`PolEscapa()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md)
- Handle challenging root positions in
  [`PlotCharacter()`](https://ms609.github.io/TreeSearch/reference/PlotCharacter.md)
- Fix character state colours in app legend
- Tweak documentation

## TreeSearch 1.6.0 (2025-04-09)

CRAN release: 2025-04-09

### Improvements

- [`PlotCharacter()`](https://ms609.github.io/TreeSearch/reference/PlotCharacter.md)
  performs ancestral state reconstruction on consensus trees
  ([\#179](https://github.com/ms609/TreeSearch/issues/179))
- Improve support for constraints in
  [`AdditionTree()`](https://ms609.github.io/TreeSearch/reference/AdditionTree.md)
  ([\#173](https://github.com/ms609/TreeSearch/issues/173))
- Support for ordered (additive) characters via
  [`TreeTools::Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html)
- Fix SPR behaviour when move is close to root

### App improvements

- Buttons to download consensus trees in app
- Fix display of state labels in app

### Housekeeping

- Require R 4.0 (to simplify maintenance)

## TreeSearch 1.5.1 (2024-05-23)

CRAN release: 2024-05-23

- Fix calls to
  [`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.html)

## TreeSearch 1.5.0 (2024-04-03)

CRAN release: 2024-04-03

- [`MaximumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  calculates maximum possible length of characters, including with
  inapplicable tokens
- [`Consistency()`](https://ms609.github.io/TreeSearch/reference/Consistency.md)
  now returns retention index and rescaled consistency index

## TreeSearch 1.4.0 (2023-08-18)

CRAN release: 2023-08-26

### New features

- [`TaxonInfluence()`](https://ms609.github.io/TreeSearch/reference/TaxonInfluence.md)
  calculates influence of individual taxa on phylogenetic inference

### Search improvements

- Default to use equal weighting during ratchet iterations
- Support null constraints in
  [`AdditionTree()`](https://ms609.github.io/TreeSearch/reference/AdditionTree.md)

### App improvements

- Exclude taxa from search in app
- Allow search to continue when loading a new file with different taxon
  names into the app

### Housekeeping

- Update calls to
  [`DescendantEdges()`](https://ms609.github.io/TreeTools/reference/DescendantEdges.html)
  for compatibility with ‘TreeTools’ 1.10.0

## TreeSearch 1.3.2 (2023-04-27)

CRAN release: 2023-04-27

- Use
  [`PlotTools::SpectrumLegend()`](https://ms609.github.io/PlotTools/reference/SpectrumLegend.html)
  for continuous scales in app
- Restore auto-termination of `.t` files

## TreeSearch 1.3.1 (2023-03-29)

CRAN release: 2023-03-29

- [`PlotCharacter()`](https://ms609.github.io/TreeSearch/reference/PlotCharacter.md)
  now returns invisibly
- Fix missing character in Wills 2012 dataset
- Search by character text in GUI
- Call C functions using symbols

## TreeSearch 1.3.0 (2023-02-20)

CRAN release: 2023-02-20

### New features

- New function
  [`LengthAdded()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md)
  tests which characters contribute to taxon instability, per Pol &
  Escapa (2009), <doi:10.1111/j.1096-0031.2009.00258.x>
- [`WhenFirstHit()`](https://ms609.github.io/TreeSearch/reference/WhenFirstHit.md)
  recovers tree search information from tree names
- New
  [vignette](https://ms609.github.io/TreeSearch/dev/articles/tree-space.html)
  on tree space mapping
- Support `phylo` trees as constraints

### GUI improvements

- Support reading characters from Excel spreadsheets
- Allow retention of suboptimal trees
- Use K-means++ clustering

## TreeSearch 1.2.0 (2022-07-35)

CRAN release: 2022-08-10

- ‘shiny’ GUI improvements:
  - Export log of tree search commands
  - Export R scripts to reproduce figures
  - Simplify layout
  - Misc bug fixes
- New function
  [`QuartetResolution()`](https://ms609.github.io/TreeSearch/reference/QuartetResolution.md)
  evaluates how a quartet is resolved in each of a list of trees

## TreeSearch 1.1.2 (2022-05-11)

CRAN release: 2022-05-11

- Check tree order & rootedness before scoring
  ([\#133](https://github.com/ms609/TreeSearch/issues/133))
- Improve error handling
- Replace `throw` with `stop` in C++
- Remove test of elapsed times, for CRAN compliance

## TreeSearch 1.1.1 (2022-03-22)

CRAN release: 2022-03-22

- GUI allows selection of subset of trees, for easier analysis of
  Bayesian tree sets
- Miscellaneous fixes and improvements in ‘shiny’ GUI
- Test suite for ‘shiny’ GUI
- Update tests for TreeSearch 1.7

## TreeSearch 1.1.0 (2022-01-17)

CRAN release: 2022-01-17

- Improvements to ‘shiny’ GUI
- Better integration of rogue taxon exploration
- New vignette describing profile parsimony
- [`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  fully supports ambiguous applicability

## TreeSearch 1.0.1 (2021-09-27)

CRAN release: 2021-09-27

- Memory management with invalid input
- Corrections to metadata

## TreeSearch 1.0.0 (2021-09-21)

CRAN release: 2021-09-22

### New functions

- [`EasyTrees()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  ‘shiny’ graphical user interface for tree search
- [`AdditionTree()`](https://ms609.github.io/TreeSearch/reference/AdditionTree.md)
  adds each taxon in sequence to the most parsimonious place on the
  tree, generating a more parsimonious starting tree than
  neighbour-joining
- [`PlotCharacter()`](https://ms609.github.io/TreeSearch/reference/PlotCharacter.md)
  reconstructs character distributions on trees
- [`ConstrainedNJ()`](https://ms609.github.io/TreeTools/reference/ConstrainedNJ.html)
  constructs starting trees that respect a constraint
- [`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.html)
  reconciles a tree with a constraint
- [`SiteConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  calculates exact site concordance
- [`ConcordantInformation()`](https://ms609.github.io/TreeSearch/reference/ConcordantInformation.md)
  evaluates signal:noise of dataset implied by a given tree
- [`PrepareDataProfile()`](https://ms609.github.io/TreeSearch/reference/PrepareDataProfile.md)
  simplifies dataset to allow partial search when multiple applicable
  tokens are present
- [`Resample()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  conducts bootstrap and jackknife resampling
- [`Consistency()`](https://ms609.github.io/TreeSearch/reference/Consistency.md)
  calculates consistency and retention ‘indices’
- [`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  calculates minimum length of character in a dataset on any tree.

### Improvements

- [`TreeLength()`](https://ms609.github.io/TreeSearch/reference/TreeLength.md)
  supports lists of trees
- Set handling of ‘gap’ token (-) when creating Morphy object with
  `gap =`
- Label nodes with split frequencies using `JackLabels(plot = FALSE)`
- Support for topological constraints during tree search
- Remove redundant function `AsBinary()`
- Drop `nTip` parameter in
  [`RandomTreeScore()`](https://ms609.github.io/TreeSearch/reference/RandomTreeScore.md)
  (infer from `morphyObj`)
- C implementations of rearrangement functions
- Improved command line interface for search progress messaging

### Deprecations

- Remove redundant internal function `LogisticPoints()`

## TreeSearch 0.4.3 (2020-07-09)

CRAN release: 2020-07-09

- Update tests for compatibility with ‘TreeTools’ v1.1.0
- Improve memory and pointer handling

## TreeSearch 0.4.2 (2020-07-07)

CRAN release: 2020-07-07

- Update tests for compatibility with ‘TreeTools’ v1.1.0

## TreeSearch 0.4.1 (2020-06-09)

CRAN release: 2020-06-09

- Compatibility with ‘TreeTools’ v1.0.0

## TreeSearch 0.4.0 (2020-02-06)

CRAN release: 2020-02-06

### New features

- [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.html),
  complementing
  [`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.html)
- Documentation with ‘pkgdown’
- [`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md)
  helper function

### Changes

- Move tree distance measures to new package
  ‘[TreeDist](https://ms609.github.io/TreeDist/)’
- Move tree utility functions to new package
  ‘[TreeTools](https://ms609.github.io/TreeTools/)’
- Rename functions
  [`MinimumSteps()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)→[`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  and
  [`FitchSteps()`](https://ms609.github.io/TreeSearch/reference/FitchSteps.md)→[`CharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md)

### Enhancements

- Improve speed of tests (by increasing probability of false positives)
- Use `message` in place of `cat`, to allow use of
  [`suppressMessages()`](https://rdrr.io/r/base/message.html)

## TreeSearch 0.3.2 (2019-06-03)

CRAN release: 2019-06-03

- Improve text, content and build speed of vignettes

## TreeSearch 0.3.1

### New features

- `NyeTreeSimilarity()` function implements the tree similarity metric
  of Nye *et al.* (2006)
- [`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)
  function implementing the Matching Split distance of Bogdanowicz &
  Giaro (2012)

### Bug fixes

- Check whether input tree is bifurcating before attempting
  rearrangements, to avoid crashes on unsupported input

## TreeSearch 0.3.0 (2019-03-21)

CRAN release: 2019-03-21

### New features

- Implement an information theoretic tree distance measure (Smith, 2020)
- Prepare for new random number generator in R3.6.0

### Deprecations

- Function `TreeSplits()` is deprecated; use
  [`as.Splits()`](https://ms609.github.io/TreeTools/reference/Splits.html)
  instead

### Bug fixes

- Correct some mistakes in the documentation

## TreeSearch 0.2.2 (2019-01-02)

CRAN release: 2019-01-02

- Correct vignette titles

## TreeSearch 0.2.1 (2018-12-07)

CRAN release: 2018-12-07

### New features

- `CollapseNodes()` and `CollapseEdges()` allow the creation of
  polytomies
- `Tree2Splits()` lists the bipartition splits implied by a tree
  topology

### Enhancements

- [`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.html)
  now supports larger trees
- Can specify tip labels directly to
  [`ReadTntTree()`](https://ms609.github.io/TreeTools/reference/ReadTntTree.html),
  to avoid reliance on generative file

### Bug fixes

- Export missing functions

## TreeSearch 0.2.0 (2018-09-10)

CRAN release: 2018-09-11

### New features

- [`RootTree()`](https://ms609.github.io/TreeTools/reference/RootTree.html)
  allows rooting of tree on incompletely specified or single-taxon
  outgroup
- `AllTBR()` returns all trees one TBR rearrangement away
- [`TBRMoves()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  reports all possible TBR rearrangements
- [`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md)
  conducts Jackknife resampling
- [`SplitFrequency()`](https://ms609.github.io/TreeTools/reference/SplitFrequency.html)
  reports frequency of clades in a forest
- [`SupportColour()`](https://ms609.github.io/TreeTools/reference/SupportColour.html)
  allows visual marking of support values
- [`ApeTime()`](https://ms609.github.io/TreeTools/reference/ApeTime.html)
  reports the creation date of an ape-exported tree
- [`SortTree()`](https://ms609.github.io/TreeTools/reference/SortTree.html)
  flips nodes into a consistent left-right order
- `AsBinary()` supports 0

### Enhancements

- `[IW]RatchetConsensus()` renamed to `[IW]MultiRatchet()`, giving a
  better description of the function’s purpose
- Don’t warn about missing EOL when reading Nexus or TNT files
- Add new 12-colour colourblind-friendly palette
- [`FitchSteps()`](https://ms609.github.io/TreeSearch/reference/FitchSteps.md)
  now supports datasets with tips not found in tree
- Improve portability of function
  [`ReadTntTree()`](https://ms609.github.io/TreeTools/reference/ReadTntTree.html)

### Bug fixes

- `[IW]MultiRatchet()` now considers trees identical even if they’ve
  been hit a different number of times

## TreeSearch 0.1.2 (2018-03-20)

CRAN release: 2018-03-19

- Update MorphyLib library to fix C warnings
- Remove non-ASCII characters from data
- Disable slow-building and problematic vignette
- Use local copy of citation style when building vignettes

## TreeSearch 0.1.0 (2018-03-14)

CRAN release: 2018-03-14

### New features

- Helper functions to read Nexus and TNT data and trees
- Brewer palette in local data to allow easier colouring

### Enhancements

- Allow additional parameters to be passed to
  [`consensus()`](https://rdrr.io/pkg/ape/man/consensus.html) via
  [`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.html)

### Bug fixes

- C11 compliance
- `IWRatchetConsensus()` now relays concavity value to subsequent
  functions
- [`ReadCharacters()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html)
  returns labels for all characters and states if `character_num = NULL`

## TreeSearch 0.0.8

CRAN release: 2018-03-08

### New features

- Added
  [`NJTree()`](https://ms609.github.io/TreeTools/reference/NJTree.html)
  function as shortcut to generate Neighbour-Joining tree from a dataset
- Add functions to allow recovery of all trees one rearrangement from
  that input

### Efficiency gains

- Separate out
  [`NNISwap()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  functions to allow more efficient rearrangement of `edgeLists`
- \[9002\] Improve efficiency by using three-pass algorithm in place of
  four-pass precursor
- \[9004\] Bootstrap search improvements

### Bug fixes

- \[9003\] User now able to specify value of concavity constant (was
  overridden to k = 4)
- \[9003\] Bootstrap replicates now scored correctly (and without
  warning) under implied weights

## TreeSearch 0.0.7

### Inapplicable tokens:

- Integrated with this package (previously in `inapplicable`)
- Handle inapplicable data via API to Martin Brazeau’s Morphy
  Phylogenetic Library

### Profile Parsimony:

- Integrated with this package (previously in `ProfileParsimony`)
- Faster calculation of concavity profiles in C
- Persistent memoization with R.cache

## TreeSearch 0.0.6

CRAN release: 2017-11-07

- First CRAN submission
