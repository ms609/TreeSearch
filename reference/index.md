# Package index

## Custom search

Functions to search for optimal phylogenetic trees.

- [`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md)
  : Jackknife resampling
- [`MorphyBootstrap()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md)
  [`Ratchet()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md)
  [`MultiRatchet()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md)
  [`RatchetConsensus()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md)
  : Parsimony Ratchet

## Tree scoring

Functions that provide tree scores, for optimization.

- [`CharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md)
  [`FastCharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md)
  : Character length
- [`ExpectedLength()`](https://ms609.github.io/TreeSearch/reference/ExpectedLength.md)
  : Expected length
- [`LengthAdded()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md)
  [`PolEscapa()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md)
  : Contribution of character to leaf instability
- [`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  [`MinimumSteps()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  [`MaximumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md)
  : Minimum and Maximum lengths possible for a character
- [`TaxonInfluence()`](https://ms609.github.io/TreeSearch/reference/TaxonInfluence.md)
  : Rank taxa by their influence on phylogenetic results
- [`IWScore()`](https://ms609.github.io/TreeSearch/reference/TreeLength.md)
  [`TreeLength()`](https://ms609.github.io/TreeSearch/reference/TreeLength.md)
  [`Fitch()`](https://ms609.github.io/TreeSearch/reference/TreeLength.md)
  : Calculate the parsimony score of a tree given a dataset
- [`ConcordantInformation()`](https://ms609.github.io/TreeSearch/reference/ConcordantInformation.md)
  [`Evaluate()`](https://ms609.github.io/TreeSearch/reference/ConcordantInformation.md)
  [`ConcordantInfo()`](https://ms609.github.io/TreeSearch/reference/ConcordantInformation.md)
  : Evaluate the concordance of information between a tree and a dataset
- [`Consistency()`](https://ms609.github.io/TreeSearch/reference/Consistency.md)
  : Consistency and retention "indices"
- [`PlotCharacter()`](https://ms609.github.io/TreeSearch/reference/PlotCharacter.md)
  : Plot the distribution of a character on a tree
- [`RandomTreeScore()`](https://ms609.github.io/TreeSearch/reference/RandomTreeScore.md)
  : Parsimony score of random postorder tree

## Tree (re)construction

Tree generation and rearrangement functions.

- [`AdditionTree()`](https://ms609.github.io/TreeSearch/reference/AdditionTree.md)
  : Addition tree

- [`RandomMorphyTree()`](https://ms609.github.io/TreeSearch/reference/RandomMorphyTree.md)
  : Random postorder tree

- [`NNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  [`cNNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  [`NNISwap()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  [`RootedNNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  [`RootedNNISwap()`](https://ms609.github.io/TreeSearch/reference/NNI.md)
  : Nearest neighbour interchange (NNI)

- [`AllSPR()`](https://ms609.github.io/TreeSearch/reference/AllSPR.md) :
  All SPR trees

- [`SPR()`](https://ms609.github.io/TreeSearch/reference/SPR.md)
  [`SPRMoves()`](https://ms609.github.io/TreeSearch/reference/SPR.md)
  [`SPRSwap()`](https://ms609.github.io/TreeSearch/reference/SPR.md)
  [`RootedSPR()`](https://ms609.github.io/TreeSearch/reference/SPR.md)
  [`RootedSPRSwap()`](https://ms609.github.io/TreeSearch/reference/SPR.md)
  : Subtree pruning and rearrangement (SPR)

- [`cSPR()`](https://ms609.github.io/TreeSearch/reference/cSPR.md) :

  [`cSPR()`](https://ms609.github.io/TreeSearch/reference/cSPR.md)
  expects a tree rooted on a single tip.

- [`TBR()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  [`TBRMoves()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  [`TBRSwap()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  [`RootedTBR()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  [`RootedTBRSwap()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
  : Tree bisection and reconnection (TBR)

- [`RearrangeEdges()`](https://ms609.github.io/TreeSearch/reference/RearrangeEdges.md)
  : Rearrange edges of a phylogenetic tree

## Tree properties

Functions that calculate support for edges.

- [`ConcordanceTable()`](https://ms609.github.io/TreeSearch/reference/ConcordanceTable.md)
  : Plot concordance table
- [`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md)
  : Label nodes with jackknife support values
- [`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md)
  : Jackknife resampling
- [`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  [`Resample()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  [`EasyTrees()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  [`EasyTreesy()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
  : Find most parsimonious trees
- [`MostContradictedFreq()`](https://ms609.github.io/TreeSearch/reference/MostContradictedFreq.md)
  : Frequency of most common contradictory split
- [`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md)
  : Group present or contradicted score
- [`ClusteringConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  [`MutualClusteringConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  [`QuartetConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  [`PhylogeneticConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  [`SharedPhylogeneticConcordance()`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
  : Concordance factors

## Profile parsimony

Functions to help calculate profile parsimony scores.

- [`Carter1()`](https://ms609.github.io/TreeSearch/reference/Carter1.md)
  [`Log2Carter1()`](https://ms609.github.io/TreeSearch/reference/Carter1.md)
  [`LogCarter1()`](https://ms609.github.io/TreeSearch/reference/Carter1.md)
  :

  Number of trees with *m* steps

- [`PrepareDataProfile()`](https://ms609.github.io/TreeSearch/reference/PrepareDataProfile.md)
  [`PrepareDataIW()`](https://ms609.github.io/TreeSearch/reference/PrepareDataProfile.md)
  : Prepare data for Profile Parsimony

- [`StepInformation()`](https://ms609.github.io/TreeSearch/reference/StepInformation.md)
  :

  Information content of a character known to contain *e* steps

- [`WithOneExtraStep()`](https://ms609.github.io/TreeSearch/reference/WithOneExtraStep.md)
  : Number of trees with one extra step

- [`profiles`](https://ms609.github.io/TreeSearch/reference/profiles.md)
  : Empirically counted profiles for small trees

## Morphy APIs

Functions that interface with the Morphy C library.

- [`GapHandler()`](https://ms609.github.io/TreeSearch/reference/GapHandler.md)
  : Read how a Morphy Object handles the inapplicable token

- [`MorphyWeights()`](https://ms609.github.io/TreeSearch/reference/MorphyWeights.md)
  [`SetMorphyWeights()`](https://ms609.github.io/TreeSearch/reference/MorphyWeights.md)
  : Set and get the character weightings associated with a Morphy
  object.

- [`PhyDat2Morphy()`](https://ms609.github.io/TreeSearch/reference/PhyDat2Morphy.md)
  :

  Initialize a Morphy object from a `phyDat` object

- [`SingleCharMorphy()`](https://ms609.github.io/TreeSearch/reference/SingleCharMorphy.md)
  : Morphy object from single character

- [`UnloadMorphy()`](https://ms609.github.io/TreeSearch/reference/UnloadMorphy.md)
  : Destroy a Morphy object

- [`is.morphyPtr()`](https://ms609.github.io/TreeSearch/reference/is.morphyPtr.md)
  : Is an object a valid Morphy object?

- [`summary(`*`<morphyPtr>`*`)`](https://ms609.github.io/TreeSearch/reference/summary.morphyPtr.md)
  : Details the attributes of a morphy object

## Datasets

Data that accompanies the package.

- [`congreveLamsdellMatrices`](https://ms609.github.io/TreeSearch/reference/congreveLamsdellMatrices.md)
  : 100 simulated data matrices
- [`inapplicable.datasets`](https://ms609.github.io/TreeSearch/reference/inapplicable.datasets.md)
  [`inapplicable.phyData`](https://ms609.github.io/TreeSearch/reference/inapplicable.datasets.md)
  [`inapplicable.trees`](https://ms609.github.io/TreeSearch/reference/inapplicable.datasets.md)
  [`inapplicable.citations`](https://ms609.github.io/TreeSearch/reference/inapplicable.datasets.md)
  : Thirty datasets with inapplicable data
- [`profiles`](https://ms609.github.io/TreeSearch/reference/profiles.md)
  : Empirically counted profiles for small trees
- [`referenceTree`](https://ms609.github.io/TreeSearch/reference/referenceTree.md)
  : Tree topology for matrix simulation

## Utility functions

Miscellaneous functions used within the package.

- [`ClusterStrings()`](https://ms609.github.io/TreeSearch/reference/ClusterStrings.md)
  : Cluster similar strings
- [`QACol()`](https://ms609.github.io/TreeSearch/reference/QACol.md)
  [`QCol()`](https://ms609.github.io/TreeSearch/reference/QACol.md)
  [`QALegend()`](https://ms609.github.io/TreeSearch/reference/QACol.md)
  : Generate colour to depict the amount and quality of observations
- [`QuartetResolution()`](https://ms609.github.io/TreeSearch/reference/QuartetResolution.md)
  : Relationship between four taxa
- [`WhenFirstHit()`](https://ms609.github.io/TreeSearch/reference/WhenFirstHit.md)
  : When was a tree topology first hit?
