# TreeSearch (development version)

## New features

- New function `PaintCharacters()` colours each character in a morphological
  dataset by the hue of the tree edges it most concordantly supports, using
  `ConcordanceTable()` MI weights averaged in CIELAB colour space.  Pairs with
  `TreeTools::PaintTree()` to visually map characters to clades.

- `attr(dataset, "weight")` now accepts non-integer character weights.  The
  C++ scoring engine still stores `int` weights internally; fractional
  inputs are rescaled to integer with a configurable precision (default
  0.001, controlled by `getOption("TreeSearch.fractional.scale", 1000L)`).
  Previously, fractional weights were silently truncated at the Rcpp
  boundary (e.g. `c(0.5, 1.7)` became `c(0L, 1L)`, dropping 50% / 41% of
  the respective characters' contributions).  Integer weights pass
  through unchanged.  `TreeLength()` and other scores are returned in
  units of `steps * scale` when fractional weights are present; within-
  run ranking is unaffected.

# TreeSearch 2.0.0

## Breaking changes

- Implied weighting now applies the missing-entries correction of
  Goloboff (2014) by default (`extended_iw = TRUE`).  Characters with
  many missing entries receive a reduced effective concavity, compensating
  for artificially low observed homoplasy.  Set `extended_iw = FALSE` to
  reproduce pre-2.0.0 behaviour.
- `MaximizeParsimony()` has an entirely new parameter interface.
  The previous `MaximizeParsimony()` (R-loop search using MorphyLib) has been
  renamed to `Morphy()`.
  Code that passes Morphy-style parameters (e.g. `ratchIter`, `tbrIter`,
  `maxHits`) to `MaximizeParsimony()` will be automatically forwarded to
  `Morphy()` with a deprecation warning.
  Update your code to call `Morphy()` directly, or adopt the new
  `MaximizeParsimony()` parameters.
  This compatibility shim will be removed in a future release.

## C++ search engine

`MaximizeParsimony()` is rewritten from the ground up with a native C++ search
engine, replacing the R-loop/MorphyLib backend for equal weights, implied
weights, and profile parsimony.  Typical searches are an order of magnitude
faster; inapplicable character handling (Brazeau _et al._ 2019) is built in.

### New features

- **`ScoreSpectrum()`**: Chao1-style landscape coverage estimator.  Treats
  distinct parsimony scores found across replicates as "species" and estimates
  how thoroughly the parsimony landscape has been sampled (Good-Turing sample
  coverage, Chao1 richness lower bound, unseen score-level fraction).  The
  Shiny app's confidence panel now displays the coverage estimate when
  sufficient replicates have been completed.  `MaximizeParsimony()` now
  returns a `replicate_scores` attribute containing per-replicate local-optimum
  scores for this purpose.

- **Multi-replicate driven search** pipeline: random Wagner tree → TBR →
  sectorial search (XSS, RSS, CSS) → ratchet → drift → tree fusing →
  final TBR.
- **Parallel search** via `nThreads`: replicates run on independent threads
  with a shared tree pool.
- **Timeout** via `maxSeconds`.
- **User-supplied starting tree**: when a `tree` argument is provided, the
  first replicate begins from that topology; subsequent replicates use
  random Wagner trees.
- **Adaptive strategy presets** via `strategy`: `"auto"` (default) selects
  `"sprint"`, `"default"`, or `"thorough"` based on the number of tips.
  Explicit parameters always override preset values.
- **Profile parsimony** runs natively in C++; no longer delegates to
  `Morphy()`.
- **Topological constraints** enforced natively in C++ (including during
  Wagner tree construction and sectorial search).
- **Per-phase timing** returned as a `timings` attribute on the result.
- **MPT enumeration**: after the main search converges, a TBR plateau walk
  from each pool tree discovers additional most-parsimonious topologies on the
  same and neighbouring score plateaus, up to `poolMaxSize`.
  This addresses a common complaint that the previous implementation returned
  only one tree when multiple MPTs exist.

### New parameters for `MaximizeParsimony()`

- `strategy` — `"auto"` (default), `"sprint"`, `"default"`, `"thorough"`,
  or `"none"`.
- `nThreads` — number of parallel worker threads (default 1).
- `maxSeconds` — wall-clock timeout (0 = no limit).
- `sprFirst` — run SPR before TBR in each replicate.
- `ratchetPerturbMode`, `ratchetPerturbMaxMoves`, `ratchetAdaptive` —
  configure ratchet perturbation (zero-weight, up-weight, mixed, adaptive).
- `driftCycles`, `driftAfdLimit`, `driftRfdLimit` — drift search parameters.
- `xssRounds`, `xssPartitions`, `rssRounds`, `cssRounds`, `cssPartitions`,
  `sectorMinSize`, `sectorMaxSize` — sectorial search parameters.
- `fuseInterval`, `fuseAcceptEqual` — tree fusing parameters.
- `poolMaxSize`, `poolSuboptimal` — tree pool management.
- `tbrMaxHits`, `wagnerStarts`, `tabuSize`.
- `nniFirst` — NNI warmup pass before SPR/TBR in each replicate; at
  ≥100 tips this substantially improves the Wagner starting-tree quality
  at negligible cost for small datasets.
- `postRatchetSectorial` — run a second XSS+RSS+CSS pass after ratchet
  perturbation; approximates TNT's interleaved sectorial pattern.
  Enabled by default in the `"thorough"` preset.
- `outerCycles`, `maxOuterResets` — repeat the full
  \[XSS/RSS/CSS → ratchet → NNI-perturbation → drift → TBR\] sequence
  _n_ times per replicate; budget is divided evenly.  Enabled in the
  `"thorough"` preset (`outerCycles = 2`).
- `wagnerBias`, `wagnerBiasTemp` — bias taxon addition order during Wagner
  tree construction toward taxa with more informative characters
  (Goloboff 2014), substantially improving starting-tree quality at large
  tip counts.
- `perturbStopFactor` — stop after `nTip × perturbStopFactor` consecutive
  replicates that fail to improve the best score; provides 2–7× speedup on
  converged searches at no score cost.
- `pruneReinsertCycles`, `pruneReinsertDrop`, `pruneReinsertSelection` —
  taxon pruning-reinsertion perturbation: drop a fraction of leaves, let
  the backbone re-optimise with TBR, then reinsert taxa greedily.
  Complementary to the ratchet (which perturbs character weights).
- `nniPerturbCycles`, `nniPerturbFraction` — stochastic NNI-perturbation:
  randomly apply NNI swaps to a fraction of internal branches and
  reconverge, escaping local optima without altering character weights.
- `annealCycles`, `annealPhases`, `annealTStart`, `annealTEnd`,
  `annealMovesPerPhase` — multi-cycle PCSA (simulated annealing
  perturbation) phase.
- `adaptiveLevel` — dynamically scale ratchet and drift effort per
  replicate based on the observed hit rate.
- `adaptiveStart` — Thompson-sampling bandit strategy for starting-tree
  selection; adapts over replicates to which strategies yield best scores.
- `enumTimeFraction` — fraction of `maxSeconds` reserved for the MPT
  plateau enumeration walk at the end of the search (default 10%).
- `intraFuse` — within-replicate tree fusing against pool donors after TBR
  polish; approximates TNT's within-replicate fusing pattern.
- `ratchetTaper` — gradually reduce ratchet perturbation probability as
  the pool stabilises, allowing finer local exploration late in the search.
- `consensusConstrain` — lock pool-consensus splits as topological
  constraints for subsequent replicates.
- `consensusStableReps` — stop when the strict consensus is unchanged for
  this many consecutive replicates (0 = disabled; set e.g. 3 to enable).
- `progressCallback` — R function called after each replicate (for custom
  progress reporting).

### Search output

- **Convergence summary**: when `verbosity > 0` (the default),
  `MaximizeParsimony()` now prints a one-line summary on exit reporting the
  best score, number of replicates completed, replicates since last
  improvement, number of distinct MPTs found, stop reason (time limit,
  target hits, perturbation-stop, or user interrupt), and elapsed time.
  The same information is available as named attributes on the returned
  tree list.

### Search optimizations

- **Collapsed-edge clip skipping**: TBR, SPR, and drift search skip
  clips at zero-length edges that provably cannot improve the score,
  reducing unnecessary evaluations on sparse data.
- **Conflict-guided sectorial search**: random sectorial search targets
  sectors around splits that conflict across pool trees.
- **Diversity-aware pool eviction**: when the tree pool is full, the most
  topologically similar entry is evicted to maintain diversity.
- **Cross-replicate consensus constraint tightening**: opt-in via
  `consensusConstrain = TRUE` in `SearchControl()`.
- **Consensus-stability early stopping**: when `consensusStableReps > 0` in
  `SearchControl()`, search stops when the strict consensus of best-score
  pool trees has been unchanged for that many consecutive replicates.
  Disabled by default.

### Batch resampling

- `Resample()` gains `nReplicates` and `nThreads` parameters for batch and
  parallel jackknife/bootstrap resampling via a single C++ call.
- `SuccessiveApproximations()` gains `concavity` and `constraint` parameters.

## Profile parsimony: multi-state support

- Profile parsimony now supports characters with up to 5 informative states
  (previously limited to 2).  Characters with 3--5 states use the recursive
  algorithm of Maddison & Slatkin (1991).
- New C++ function `MaddisonSlatkin()` computes the number of labelled
  histories for multi-state characters.

## Data simulation

- New function `ParsSim()` simulates morphological datasets under a parsimony
  model (equal weights, implied weights, or profile parsimony).  Each
  character starts at minimum steps; extra steps are placed one at a time,
  verified to increase the Fitch score by exactly 1.

## Scoring

- `TreeLength()` and `CharacterLength()` / `FastCharacterLength()` use the
  C++ engine for all scoring modes (equal weights, implied weights, profile
  parsimony).

## Function rename

- `TaxonInfluence()` now uses `MaximizeParsimony()` internally.
- `AdditionTree()` now uses the C++ Wagner tree engine, with native support
  for implied weights, profile parsimony, and constraints.

## Bug fixes

- Shiny: scoring error notification now shows the actual error message
  (e.g. "Trees have different numbers of edges") rather than the generic
  "Could not score all trees with dataset".
- Shiny: fix search requiring two clicks to start when trees have mixed
  topologies (polytomous/binary).  The "Search" shortcut button now appears
  only after the modal is dismissed via its own Search button, so it is never
  obscured by the modal backdrop.
- Fix output trees from `MaximizeParsimony()` having invalid preorder
  numbering (affected `DropTip()`, distance calculations, and plotting).
- Fix `fuseInterval = 0` causing a crash (division by zero).
- Fix `is_uninformative()` misclassifying ambiguous characters as
  uninformative.
- Fix `compute_fixed_steps()` undercount for all-ambiguous characters.
- Fix IW scoring with missing `min_steps` offset.
- Fix crash when dataset contains only ambiguous (`?`) tokens.

## Custom search functions

- `Ratchet()`, `MultiRatchet()`, `Jackknife()`, `MorphyBootstrap()`, and
  `TreeSearch()` are no longer deprecated.  These functions support pluggable
  `TreeScorer` and `EdgeSwapper` functions for custom scoring strategies;
  for standard parsimony, use `MaximizeParsimony()`.

## App improvements (`EasyTrees()`)

- **Async search**: the session remains responsive while a search is running.
- **Parallel search**: the search settings modal includes a thread count slider
  (when multiple cores are available).
- **Tree accumulation**: repeated "Continue search" runs accumulate trees at
  the same optimal score, with de-duplication by topology.
- **Search confidence**: after each search, the results pane shows the hit rate
  and an estimate of the replicates needed for 95% confidence.
- **Search config modal** reorganized into labelled sections (step weighting,
  parallelization, search intensity, results to keep).
- Fix `PlotCharacter()` crash on multifurcating consensus trees.
- Fix first search not appearing to update trees in memory.
- Clarified "Stop after best score found N times" slider label with help text.
- Dataset-adaptive timeout default (1–15 minutes based on dataset size).
- Internal modularization of the Shiny app into proper Shiny modules.

## Other improvements

# TreeSearch 1.8.0.9001 (2026-04-23)

- Reorder parameters in `Q[A]Col(quality, amount)`.

# TreeSearch 1.8.0.9000 (2026-02-05)

- New parameters for flexible plotting of `QALegend()`.
- `ConcordanceTable()` gains `plot` parameter.


# TreeSearch 1.8.0 (2026-01-15)

- Implements the methods of Smith (forthcoming) via `ClusteringConcordance()`,
  with visualization functions `ConcordanceTable()`, `QACol()` and `QALegend()`.
- `QuartetConcordance()` gains `return` parameter and fast C++ implementation.
- Fix regression in `MaximumLength()`.


# TreeSearch 1.7.0 (2025-08-22)

- `PresCont()` implements the Group Present / Contradicted measure of
  Goloboff et al. (2003).
- `Consistency()` also returns the relative homoplasy index of Steell et al. 
  (2023).
- `JackLabels()` supports multiple trees per iteration
  ([#197](https://github.com/ms609/TreeSearch/discussions/197))
- Support single-character matrices in `ClusteringConcordance()`
- Fix `DoNothing(x)` to return `x` (not `NULL`)
- Remove unused `delete_rawdata()` due to implementation issues.
- Port `MaximumLength()` to C++ to handle more characters, more efficiently.


# TreeSearch 1.6.1 (2025-06-10)

- Handle invariant characters in `PolEscapa()`
- Handle challenging root positions in `PlotCharacter()`
- Fix character state colours in app legend
- Tweak documentation

# TreeSearch 1.6.0 (2025-04-09)

## Improvements
- `PlotCharacter()` performs ancestral state reconstruction on consensus trees
  ([#179](https://github.com/ms609/TreeSearch/issues/179))
- Improve support for constraints in `AdditionTree()`
  ([#173](https://github.com/ms609/TreeSearch/issues/173))
- Support for ordered (additive) characters via `TreeTools::Decompose()`
- Fix SPR behaviour when move is close to root

## App improvements
- Buttons to download consensus trees in app
- Fix display of state labels in app

## Housekeeping
- Require R 4.0 (to simplify maintenance)


# TreeSearch 1.5.1 (2024-05-23)

- Fix calls to `DescendantEdges()`


# TreeSearch 1.5.0 (2024-04-03)

- `MaximumLength()` calculates maximum possible length of characters, including
  with inapplicable tokens
- `Consistency()` now returns retention index and rescaled consistency index


# TreeSearch 1.4.0 (2023-08-18)

## New features
- `TaxonInfluence()` calculates influence of individual taxa on 
  phylogenetic inference
  
## Search improvements
- Default to use equal weighting during ratchet iterations
- Support null constraints in `AdditionTree()`

## App improvements
- Exclude taxa from search in app
- Allow search to continue when loading a new file with different taxon names
  into the app
  
## Housekeeping
- Update calls to `DescendantEdges()` for compatibility with 'TreeTools' 1.10.0


# TreeSearch 1.3.2 (2023-04-27)

- Use `PlotTools::SpectrumLegend()` for continuous scales in app
- Restore auto-termination of `.t` files


# TreeSearch 1.3.1 (2023-03-29)

- `PlotCharacter()` now returns invisibly
- Fix missing character in Wills 2012 dataset
- Search by character text in GUI
- Call C functions using symbols


# TreeSearch 1.3.0 (2023-02-20)

## New features
- New function `LengthAdded()` tests which characters contribute to taxon
  instability, per Pol & Escapa (2009), doi:10.1111/j.1096-0031.2009.00258.x
- `WhenFirstHit()` recovers tree search information from tree names
- New [vignette](https://ms609.github.io/TreeSearch/dev/articles/tree-space.html) on tree space mapping
- Support `phylo` trees as constraints

## GUI improvements
- Support reading characters from Excel spreadsheets
- Allow retention of suboptimal trees
- Use K-means++ clustering


# TreeSearch 1.2.0 (2022-07-35)

- 'shiny' GUI improvements:
  - Export log of tree search commands
  - Export R scripts to reproduce figures
  - Simplify layout
  - Misc bug fixes

- New function `QuartetResolution()` evaluates how a quartet is resolved in
  each of a list of trees


# TreeSearch 1.1.2 (2022-05-11)

- Check tree order & rootedness before scoring ([#133](https://github.com/ms609/TreeSearch/issues/133))
- Improve error handling
- Replace `throw` with `stop` in C++
- Remove test of elapsed times, for CRAN compliance


# TreeSearch 1.1.1 (2022-03-22)

- GUI allows selection of subset of trees, for easier analysis of Bayesian
  tree sets
- Miscellaneous fixes and improvements in 'shiny' GUI
- Test suite for 'shiny' GUI
- Update tests for TreeSearch 1.7


# TreeSearch 1.1.0 (2022-01-17)

- Improvements to 'shiny' GUI
- Better integration of rogue taxon exploration
- New vignette describing profile parsimony
- `MinimumLength()` fully supports ambiguous applicability


# TreeSearch 1.0.1 (2021-09-27)

- Memory management with invalid input
- Corrections to metadata


# TreeSearch 1.0.0 (2021-09-21)

## New functions
 - `EasyTrees()` 'shiny' graphical user interface for tree search
 - `AdditionTree()` adds each taxon in sequence to the most parsimonious place
   on the tree, generating a more parsimonious starting tree than
   neighbour-joining
 - `PlotCharacter()` reconstructs character distributions on trees
 - `ConstrainedNJ()` constructs starting trees that respect a constraint
 - `ImposeConstraint()` reconciles a tree with a constraint
 - `SiteConcordance()` calculates exact site concordance
 - `ConcordantInformation()` evaluates signal:noise of dataset implied by a
   given tree
 - `PrepareDataProfile()` simplifies dataset to allow partial search when
   multiple applicable tokens are present
 - `Resample()` conducts bootstrap and jackknife resampling
 - `Consistency()` calculates consistency and retention 'indices'
 - `MinimumLength()` calculates minimum length of character in a dataset on any
   tree.

## Improvements
 - `TreeLength()` supports lists of trees
 - Set handling of 'gap' token (-) when creating Morphy object with `gap = `
 - Label nodes with split frequencies using `JackLabels(plot = FALSE)`
 - Support for topological constraints during tree search
 - Remove redundant function `AsBinary()`
 - Drop `nTip` parameter in `RandomTreeScore()` (infer from `morphyObj`)
 - C implementations of rearrangement functions
 - Improved command line interface for search progress messaging
 
## Deprecations
 - Remove redundant internal function `LogisticPoints()`


# TreeSearch 0.4.3 (2020-07-09)

 - Update tests for compatibility with 'TreeTools' v1.1.0
 - Improve memory and pointer handling
 
 
# TreeSearch 0.4.2 (2020-07-07)

 - Update tests for compatibility with 'TreeTools' v1.1.0


# TreeSearch 0.4.1 (2020-06-09)

 - Compatibility with 'TreeTools' v1.0.0


# TreeSearch 0.4.0 (2020-02-06)

## New features
 - `PhyDatToMatrix()`, complementing `MatrixToPhyDat()`
 - Documentation with 'pkgdown'
 - `JackLabels()` helper function
 
## Changes
 - Move tree distance measures to new package '[TreeDist](https://ms609.github.io/TreeDist/)'
 - Move tree utility functions to new package '[TreeTools](https://ms609.github.io/TreeTools/)'
 - Rename functions `MinimumSteps()`→`MinimumLength()` and 
   `FitchSteps()`→`CharacterLength()`

## Enhancements
 - Improve speed of tests (by increasing probability of false positives)
 - Use `message` in place of `cat`, to allow use of `suppressMessages()`


# TreeSearch 0.3.2 (2019-06-03)

 - Improve text, content and build speed of vignettes


# TreeSearch 0.3.1 

## New features
 - `NyeTreeSimilarity()` function implements the tree similarity metric of
   Nye _et al._ (2006)
 - `MatchingSplitDistance()` function implementing the Matching Split distance of 
   Bogdanowicz & Giaro (2012)

## Bug fixes
 - Check whether input tree is bifurcating before attempting rearrangements,
   to avoid crashes on unsupported input


# TreeSearch 0.3.0 (2019-03-21)

## New features
 - Implement an information theoretic tree distance measure (Smith, 2020)
 - Prepare for new random number generator in R3.6.0

## Deprecations
 - Function `TreeSplits()` is deprecated; use `as.Splits()` instead

## Bug fixes
 - Correct some mistakes in the documentation


# TreeSearch 0.2.2 (2019-01-02)

 - Correct vignette titles


# TreeSearch 0.2.1 (2018-12-07)

## New features
 - `CollapseNodes()` and `CollapseEdges()` allow the creation of polytomies
 - `Tree2Splits()` lists the bipartition splits implied by a tree topology

## Enhancements
 - `SplitFrequency()` now supports larger trees
 - Can specify tip labels directly to `ReadTntTree()`, to avoid reliance on
   generative file

## Bug fixes
 - Export missing functions


# TreeSearch 0.2.0 (2018-09-10)

## New features
 - `RootTree()` allows rooting of tree on incompletely specified
    or single-taxon outgroup
 - `AllTBR()` returns all trees one TBR rearrangement away
 - `TBRMoves()` reports all possible TBR rearrangements
 - `Jackknife()` conducts Jackknife resampling
 - `SplitFrequency()` reports frequency of clades in a forest
 - `SupportColour()` allows visual marking of support values
 - `ApeTime()` reports the creation date of an ape-exported tree
 - `SortTree()` flips nodes into a consistent left-right order
 - `AsBinary()` supports 0
 
## Enhancements
 - `[IW]RatchetConsensus()` renamed to `[IW]MultiRatchet()`, giving a better
     description of the function's purpose
 - Don't warn about missing EOL when reading Nexus or TNT files
 - Add new 12-colour colourblind-friendly palette
 - `FitchSteps()` now supports datasets with tips not found in tree
 - Improve portability of function `ReadTntTree()`

## Bug fixes
 - `[IW]MultiRatchet()` now considers trees identical even if they've been hit 
   a different number of times


# TreeSearch 0.1.2 (2018-03-20)

- Update MorphyLib library to fix C warnings
- Remove non-ASCII characters from data
- Disable slow-building and problematic vignette
- Use local copy of citation style when building vignettes


# TreeSearch 0.1.0 (2018-03-14)

## New features
- Helper functions to read Nexus and TNT data and trees
- Brewer palette in local data to allow easier colouring

## Enhancements
- Allow additional parameters to be passed to `consensus()` via
 `ConsensusWithout()`

## Bug fixes
- C11 compliance
- `IWRatchetConsensus()` now relays concavity value to subsequent functions
- `ReadCharacters()` returns labels for all characters and states if
  `character_num = NULL`


# TreeSearch 0.0.8

## New features
- Added `NJTree()` function as shortcut to generate Neighbour-Joining tree from 
    a dataset
- Add functions to allow recovery of all trees one rearrangement from that input

## Efficiency gains
- Separate out `NNISwap()` functions to allow more efficient rearrangement of 
  `edgeLists`
- [9002] Improve efficiency by using three-pass algorithm in place of four-pass precursor
- [9004] Bootstrap search improvements

## Bug fixes
- [9003] User now able to specify value of concavity constant
  (was overridden to k = 4)
- [9003] Bootstrap replicates now scored correctly (and without warning)
  under implied weights


# TreeSearch 0.0.7

## Inapplicable tokens:
- Integrated with this package (previously in `inapplicable`)
- Handle inapplicable data via API to Martin Brazeau's Morphy Phylogenetic Library

## Profile Parsimony:
- Integrated with this package (previously in `ProfileParsimony`)
- Faster calculation of concavity profiles in C
- Persistent memoization with R.cache


# TreeSearch 0.0.6
- First CRAN submission
