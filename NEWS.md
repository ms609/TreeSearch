# To integrate into 2.0.0 notes

- Profile parsimony computes exactly for more multi-state characters, where it
  previously approximated nearly all of them.  The exact Maddison & Slatkin
  solver caches into fixed-capacity memo tables and bails out when one fills --
  a guard added to stop an unbounded probe loop -- but its reserved size was
  never matched to the feasibility gate that feeds it.  Measured against the
  worst character that gate admits, every one of them overflowed: a 3-state
  character needs up to ~28,000 memo entries against the 4,096 reserved.
  `StepInformation()` and `PrepareDataProfile(approx = "auto")` therefore fell
  back to the Monte Carlo approximation for essentially every multi-state
  character -- a documented mode, but not the one asked for.

  The 2 s wall-clock budget is unchanged, and remains what caps the wait: a
  character that cannot be solved within it still falls back to Monte Carlo.
  Only the characters that fit inside that budget are affected.

  **Information amounts for those characters will therefore change**, from a
  sampled estimate to the exact value.  Which characters those are depends on
  how fast the machine is, since the budget is what decides; pass
  `approx = "mc"` for the previous behaviour throughout.  Note that
  `approx = "exact"` waives the feasibility gate but not the budget, so it too
  can fall back on a slow machine.

  Under sanitizer builds the budget is scaled by the instrumentation's
  slowdown.  Those builds run one to two orders of magnitude slower, so a 2 s
  budget tripped on everything -- leaving the sanitizer inspecting the fallback
  rather than the algorithm it was aimed at.  There is no responsiveness to
  protect in a nightly memory check.

- `constraint` now binds the trees `MaximizeParsimony()` returns, at three
  boundaries where it did not.  A starting tree supplied through `tree` was
  never checked against the constraint; because a constrained search rejects
  every rearrangement away from a violating tree, the replicate froze on it and
  reported a score no constraint-satisfying tree could reach, which then evicted
  the compliant trees other replicates had found.  A violating start is now
  rearranged until it complies before the search begins, **with a warning**.
  Separately, a replicate's own tree entered the pool unchecked, and the final
  collapse of unsupported branches could contract the very branch that displayed
  an enforced grouping -- so under the default `collapse = TRUE` a returned tree
  could break the constraint outright.  Both paths are now checked.

  **Constrained results may therefore differ from previous versions**: scores
  can rise to the true constrained optimum, and returned trees will display the
  constrained groupings.  `MaximizeParsimony()` also warns if any replicate
  ended on a tree that could not be made to satisfy the constraint, and now
  raises an error rather than returning an unverified tree if no
  constraint-satisfying tree was found at all.
- Every part of the search now reads `constraint` the way it is documented: a
  tree is compatible with a constraint character when some edge separates the
  taxa coded `1` from those coded `0`, with `?`-coded and unmentioned taxa free
  to fall on either side.  The locked-node filter that screens individual
  rearrangements, the constrained Wagner build and the collapse pass previously
  required the `1` group to be a clade *exactly*, free taxa excluded.  That is
  strictly stronger, so the search never accepted a rearrangement that broke the
  documented constraint; but a start tree that satisfied the documented
  constraint without making either group an exact clade matched no node, every
  rearrangement was rejected, and the replicate returned its start unimproved.
  Constrained searches with `?`-coded taxa therefore reach better scores.
  The exact match also blunted the collapse protection described above: with
  free taxa it matched no branch, so the separating edge could still be
  contracted away -- the one route by which a *returned* tree could break the
  constraint.
- Random starting trees under a constraint now sample every topology the
  constraint permits.  Every tree the old generator produced was compliant, but
  it built each "together" group as an exact clade with the `?`-coded taxa held
  outside, so most compliant topologies could never be drawn at all: 15 of the
  35 on six taxa with one constraint character, and 105 of the 1155 on eight
  taxa with two.  Both are now drawn in full, and at close to equal rates.
  Constrained searches that use random starts therefore begin from the whole
  range of legal trees rather than one corner of it.
- A constraint character whose `1` or `0` group holds fewer than two taxa now
  warns and is ignored, rather than being enforced as a clade.  Every tree
  separates such a group from the rest, so the character constrains nothing
  under the documented reading.  The test is symmetric in the two groups, which
  the old one was not: `c(a = 1, b = 1, c = 0)` and `c(a = 0, b = 0, c = 1)`
  state the same constraint and are now treated the same way.  Code the taxa
  that must fall outside a group as `0`, rather than leaving them `?`, to keep
  it enforced.
- `TreeLength()`, `CharacterLength()`, `TreeScore()` and `EdgeListScore()` -- and
  so `Consistency()`, `ExpectedLength()`, `ConcordantInformation()`,
  `LengthAdded()` and `SuccessiveApproximations()`, which score trees through
  them -- now reject a
  tree that contains a polytomy, with the "`tree` must be binary" error that
  `TreeLength()` already gave for a single `phylo` tree.  Such a tree
  previously returned a number.  The scoring engine derives its node counts from
  the number of edges, which identifies a tree only if that tree is binary: a
  polytomous tree with an odd number of edges wrote past the end of the arrays
  holding its topology, and one with an even number of edges was rooted on a
  leaf and then scored from memory outside its own state buffer, so repeating
  the same call could return a different answer each time.  `MaximizeParsimony()`
  collapses the trees it returns unless `collapse = FALSE`, so scoring its output
  reached this path; search with `collapse = FALSE` to obtain trees that can be
  scored, whose lengths are the score the search reports.  Resolving a collapsed
  tree instead, with `TreeTools::MakeTreeBinary()`, does not recover that score:
  an arbitrary resolution of a polytomy need not be one of the most parsimonious
  ones.

- `inapplicable = "xform"` scores are now reported at a canonical rooting, so a
  reported score is reproducible.  The x-transformation's step matrix is
  asymmetric -- a gain costs one more than the number of secondary characters it
  brings into existence, against 1 for a loss -- which makes a tree's length
  depend on where it is rooted, unlike parsimony under the symmetric criteria.
  `MaximizeParsimony()` recorded its best score mid-search at whatever rooting
  the replicate held, while returning trees re-rooted on the first taxon, so
  `attr(result, "score")` did not match `TreeLength()` of the very tree returned
  (measured: 178 reported against 183 returned on a 36-taxon matrix), and
  re-rooting a returned tree changed its length again.  Both boundaries now
  canonicalise on the first taxon of `dataset`, so one topology has one length
  and the two agree by construction.

  **X-transformation scores may therefore differ slightly from previous
  versions**, and will not decrease: the reported value is the length of the tree
  you are handed rather than of a rooting discarded during search.  It is an
  upper bound on the rooting-free minimum, exceeding it by at most the total
  number of secondary characters across hierarchy blocks (attained exactly by
  87--98% of rootings in simulation).  This changes reporting only -- what the
  search optimises is untouched.

  `MaximizeParsimony()` now also warns when the trees it returns do not share a
  length at that common rooting, which can happen because pool membership is
  still decided on scores taken at differing rootings.  Only the x-transformation
  is affected; HSJ reporting is deliberately unchanged, since there
  rooting-invariance is a property the method requires rather than a convention
  to pick -- and it is now delivered, as the next entry describes.

- `inapplicable = "hsj"` scores no longer depend on where the tree is rooted.
  Hopkins & St John (2021) define the score as a minimum over internal-node
  labellings of a sum of *symmetric* dissimilarities across the branches of an
  *unrooted* tree, so rooting-invariance is required by the method rather than
  merely desirable.  Two defects broke it, both in how the secondary characters
  were labelled; the underlying present/absent dynamic programme was correct
  throughout.

  First, the inapplicable token was treated as an ordinary state of a secondary
  character.  Where a controlling primary codes a structure absent, its
  secondaries do not exist, so `"-"` there is not a state the character takes;
  admitting it let a node in the middle of a region where the structure *is*
  present be labelled "inapplicable", mismatching every secondary at once and
  charging that branch the full weight of the scaling parameter.  Secondaries
  are now unconstrained at tips whose primary may code the structure absent.
  Second, the remaining ambiguity was resolved by a pass whose direction was a
  property of the input rooting; that pass is now rooted canonically on the
  first taxon, inside the scoring kernel, so the labelling depends only on the
  unrooted topology.

  **HSJ scores on data with a mix of present and absent primaries may therefore
  differ from previous versions.**  Most do not: across 180 simulated
  tree--matrix pairs, 152 were unchanged, 27 fell (by up to 1) and one rose (by
  0.25).  Falls are the removal of spurious inapplicable mismatches; a rise is
  possible because a score is now taken at a fixed canonical rooting rather
  than at whichever rooting the tree happened to arrive in, and that rooting
  was sometimes the flattering one.  Scores are unchanged wherever every taxon
  shares the controlling primary's presence, and Figure 1 of the paper still
  scores 7 and 5.

  Unlike the x-transformation change above, this one alters what the search
  optimises: the criterion is now a function of the unrooted tree, so
  `MaximizeParsimony()`'s reported score matches `TreeLength()` of the trees it
  returns without any re-scoring, and every tree in a returned set shares that
  length.  A secondary character is treated as unconstrained only at taxa whose
  controlling primary cannot code the structure present; where the primary is
  ambiguous but a secondary was observed, that observation still counts.

- `inapplicable = "hsj"` scoring fixed an index-space confusion that could
  under- or over-count the controlling primary's gains and losses, and could
  score an ambiguous (`"?"`) secondary character as though it were a specific,
  conflicting state.  Internally, a tip's data value was read as an index into
  the dataset's character *states*, but it is actually an index into the
  dataset's observed *tokens* (a distinct, dataset-specific ordering that only
  sometimes coincides with state order) -- so, depending on a dataset's
  internal token ordering, an absent or inapplicable primary could be scored
  as present, a present primary as absent, and a genuinely ambiguous secondary
  character as a forced, arbitrary state.  This was independent of tree
  topology, so no search or comparison using `inapplicable = "hsj"` was
  reliable: the same dataset and tree could report different scores merely by
  virtue of the order characters happened to appear in.

  **HSJ scores may therefore differ from previous versions**, in either
  direction: scores typically rise where a genuinely present or absent
  controlling primary is now always counted, but can also fall where an
  ambiguous secondary is no longer forced into a spurious mismatch.  This is
  a correctness fix to how a tip's data is looked up; it does not touch the
  known rooting-sensitivity of HSJ scoring, which remains a separate, open
  issue.

- Zero-length-branch collapse (`collapse = TRUE`) no longer disables itself
  for an `inapplicable = "hsj"`/`"xform"` search whenever *no* hierarchy
  block actually exists in that replicate -- previously it keyed on the
  scoring mode alone.  This only affects `Resample()`, whose bootstrap and
  jackknife replicates can drop every hierarchy block from a unit while
  still passing a (now-empty) hierarchy config through; those replicates are
  ordinary Fitch data and now collapse like any other.  A replicate that
  retains any hierarchy block is unaffected.

- `inapplicable = "hsj"` scoring no longer forms a reference one element past
  the end of an internal vector.  The secondary-labelling uppass computed a
  pointer to a node's children before testing whether it had any, and for a
  childless node reached after the traversal had emitted its last child that
  pointer addressed one past the end.  No
  value was ever read through it and no score changed -- 900 of 900 HSJ and
  x-transformation lengths are bit-identical either side of the fix -- but the
  access is undefined behaviour, and any build whose standard library checks
  its own preconditions aborted on it.  That includes the container behind the
  `gcc-ASAN` workflow, which is why that workflow could not get past this
  package: it stopped on the library assertion rather than on anything the
  sanitizer itself had found.

- `MaximizeParsimony(effort = )` replaces `strategy = `, which is removed (it
  was never released).  `effort` is a **relative** offset, not an absolute
  level: `0` (the default) accepts the amount of search the dataset's size and
  character count warrant, `1` asks for one notch more, `-1` one less.  So a
  single call means "try harder than usual" whether the matrix has 20 taxa or
  200, and a user never has to know which preset it would otherwise have got.
  `effort = 0` reproduces the previous `strategy = "auto"` behaviour exactly on
  every size band.

  The rungs are the former presets — `sprint`, `default`, `thorough`, and
  `large` (which was only ever `thorough` with `maxReplicates = 500`) — so the
  ladder generalises an axis the package already had.  Beyond `large`, each
  further notch doubles BOTH the replicate budget (1000, 2000, 4000 ...) and
  the hit target, so one notch always means roughly twice the work whichever
  bound a dataset is under.  There is no policy ceiling: extra replicates cost
  wall but cannot cost reach, so the ladder stops only at rung 26, where the
  budget outgrows R's integer type.

  The rung-4 budget of 500 is measured (a 34-matrix 120--180-tip sweep found
  reach climbing from 0.68 at 96 replicates to 0.79 at 250, with the hard subset
  still climbing at 500 and no knee).  The doubling above it is an operating
  point, not a fitted constant -- nothing measures where the reach curve
  flattens, and a doubling grid over rungs 4--8 on the hard tail is what would
  replace the guess with a measurement.

  The replicate budget climbs first, because `targetHits` cannot act once that
  budget is reached -- and on hard datasets it always is.  Measured on 30
  inapplicable-bearing matrices: tripling the hit target bought 4409 extra
  replicates in total, but only 243 of them on the six matrices with anything
  left to find, and NONE on the three hardest, where the replicate cap bound
  every run of both arms.  A ladder raising the hit target first would spend its
  effort almost entirely on datasets that were already solved.

  `targetHits` is raised in step regardless, for reasons that are not reach: it
  governs when easy runs stop, so without it a notch would be inert on every
  dataset that finishes early, and under implied weights it additionally deepens
  the ratchet.  Higher rungs therefore buy confidence and distinct trees on easy
  data, and reach on hard data -- not reach uniformly.

  Anything set explicitly still wins: a `maxReplicates` or `targetHits` you
  supply is never rescaled by `effort`, and explicit `control` fields continue
  to override the rung's preset.

- Fixed: `AdditionTree(constraint = )` silently returned a constraint-violating
  tree for around one addition order in eleven.  Taxa are added to a tree seeded
  from the first three of them, which is built before the constraint is
  consulted; whenever that seed put a constrained group's taxa on both sides of
  its root, the group's ancestor was the root itself and the constraint was then
  ignored for every subsequent insertion -- without a warning, and irrecoverably,
  since an ancestor never moves back down.  A group in that position cannot be
  made monophyletic by adding leaves, so the constraint is now enforced through
  its complement, which displays the same unrooted split.  Measured on the
  previous code, a randomized `sequence` hit this in 35 of 400 seeds, and an
  explicit `sequence` whose first three taxa fall inside the constrained group
  hit it in 42 of the 120 base triples.

  A constrained Wagner tree is now also re-rooted on its first taxon before being
  returned or searched from.  The topology and score are unchanged -- the rooting
  of an addition tree is an arbitrary artefact of the order taxa were added in,
  as `?AdditionTree` notes -- but it is the rooting in which the rest of the
  constraint machinery can recognise every split the tree displays.  Without it,
  a constrained search took several times longer to reach the same score.

- Fix `TreeLength()` and `LengthAdded()` errors when scoring, under profile
  parsimony, a character with no phylogenetic information.

- New `SearchControl()` parameter `stopPatience`: stop after this many consecutive
  replicates fail to improve the best score.  Unlike `perturbStopFactor` it is a
  flat count, referring neither to the tip count nor to the number of hits, so the
  replicate at which it fires does not stretch as replicates become individually
  more expensive.  The count resets on every improvement, so a search that keeps
  improving is never cut short.  0 (the default) disables it.  As with the other
  no-improvement rules, it acts precisely only in a serial search: with
  `nThreads > 1` it is evaluated when the coordinating thread polls, so it fires
  later and less predictably.

- Implied-weights searches under `strategy = "sprint"` or `"default"` now run a
  deeper ratchet paid for by that flat patience: `sprint` takes
  `ratchetCycles = 12`, `ratchetPerturbProb = 0.25` and `stopPatience = 20`;
  `default` takes `ratchetCycles = 20` and `stopPatience = 15`.  The two knobs
  ship together because each fails on its own — the deeper ratchet improves the
  score but costs wall-clock, and stopping earlier saves wall-clock but costs
  score.  Measured over 68 training matrices (6 seeds, k = 10): on the median
  matrix `sprint` is 26% faster and `default` 18% faster, `sprint` scores better
  on 4 matrices and worse on none, and `default` better on 9 and worse on 3.
  A sweep over `stopPatience` in {10, 15, 20, 25, 30} found score and wall-clock
  both vary smoothly with the value, so these are operating points on a
  trade-off rather than tuned constants: a larger `stopPatience` buys score back
  and gives up the speed.  Not every matrix gets faster — 9 of the 44 `default`
  training matrices were more than 10% slower, being those where the patience
  rule does not fire and the deeper ratchet is not paid for.  Equal weights and
  profile parsimony are unchanged, as is `thorough`/`large`, and setting any of
  these fields yourself overrides all of it.

- Fixed: a large `targetHits` combined with a large `perturbStopFactor` stopped the
  search after two replicates and silently returned a worse tree.  The
  no-improvement rule computes `(targetHits / hits) * nTip * perturbStopFactor`,
  which overflowed for such settings; the out-of-range value became a negative
  limit, so the rule fired on the first replicate that failed to improve.  The
  limit now saturates, so settings that ask for more search get more search.
  Affected both the serial and the parallel search paths.  Note that `0`, not a
  large value, is the way to switch these rules off.

- Implied-weights searches under `strategy = "thorough"` or `"large"` now run a
  deeper parsimony ratchet (48 cycles, up from 20).  Under implied weights the
  optimum can sit in a small basin a fraction of a step below an easy-to-find
  near-optimum, and character reweighting is what crosses that gap: extra
  replicates do not substitute for it.  Over a 36-matrix grid this roughly halved
  expected time-to-optimum on the matrices that are sensitive to ratchet depth,
  and cost the rest a median 0.2 s with no change in the score reached.  Equal
  weights is unchanged: the same comparison found no benefit there.  Raising
  `targetHits` now deepens the ratchet in proportion (capped at 115 cycles),
  since no dataset property reliably predicts how much reweighting a matrix
  needs; that escalation is measured as neutral rather than beneficial, and is
  offered because a raised `targetHits` is the user's own signal that the dataset
  is hard.  Lowering `targetHits` does not make the ratchet shallower, and
  setting `ratchetCycles` yourself overrides all of this.

- `MaximizeParsimony()` now normalizes `concavity` before dispatching to the
  search engine, instead of coercing it silently later.  Previously,
  `concavity = "10"` (a numeric-coercible string) skipped the R-side
  minimum-steps calculation entirely, yet still reached the C++ engine as a
  finite value and ran in implied-weighting mode with homoplasy uncorrected
  — a silently wrong score with no error or warning.  Separately,
  `concavity = "Profile"` or `"prof"` failed the search entry's exact-match
  check (unlike scoring functions such as `TreeLength()`, which already
  matched case- and prefix-insensitively) and silently searched under equal
  weights instead of profile parsimony.  Both symptoms are now closed:
  profile-mode matching is shared with the scoring entry points, and any
  other `concavity` value is coerced with `as.numeric()` and rejected with a
  clear error if it is not a single positive number (or `Inf`).
  `concavity = 10`, `concavity = "profile"`, and `concavity = Inf` behave
  exactly as before.

- `MaximizeParsimony()` now rejects `maxReplicates < 1` with a clear error,
  rather than silently returning the random starting tree tagged with a
  bogus `attr(, "score")` of `-1` (the search loop ran zero times, so the
  pool was left empty).

- `SearchControl()` now validates `enumTimeFraction`, `nniPerturbFraction`,
  and `ratchetPerturbProb` against their documented ranges, and rejects
  `sectorMinSize > sectorMaxSize`.  Out-of-range values previously reached
  the C++ engine unchecked and produced silently degenerate search
  behaviour rather than an error.

- Fixed: the "increase `maxReplicates`" advisory warning in
  `MaximizeParsimony()` reported the character count using the
  internally-scaled integer weight (up to ~1260x the true value for
  fractional-weight datasets) rather than the dataset's actual weights,
  inflating the recommended replicate count. Cosmetic only: no change in
  search behaviour.

- `MaximizeParsimony(tree = )` now accepts a whole pool of starting trees:
  given a `multiPhylo`, replicate _i_ warm-starts from tree _i_, and any
  replicates beyond the pool build random Wagner trees as before.  Previously
  only the first tree of a `multiPhylo` was used and the rest were silently
  discarded, so resuming a search from a previous run's most-parsimonious trees
  threw away exactly the topological diversity that tree fusing exploits.  All
  supplied trees must bear the same tip labels.  One tree is consumed per
  replicate actually run, so a search that converges on `targetHits` before
  exhausting a large pool now warns rather than discarding the remainder
  silently.  Passing a single `phylo` searches exactly as it did before.

- `MaximizeParsimony()` now rejects a structurally invalid `tree` with an R
  error instead of crashing the session.  `ape::unroot()` accepts \pkg{TreeTools}'
  `order = "preorder"` attribute and then mishandles it, so unrooting a
  \pkg{TreeTools} tree returns an edge matrix containing `NA`; passing one on
  segfaulted inside the rooting code, below the level at which R can trap
  anything.  Valid trees, rooted or unrooted, are unaffected.

- `adaptiveStart = TRUE` (set by `strategy = "thorough"`) no longer credits a
  starting-tree strategy for replicates that began from a tree supplied via
  `tree = `.  Such replicates never build a start of their own, so the
  Thompson-sampling bandit was recording a `Wagner(random)` trial that never
  happened — one per supplied tree, silently biasing arm selection for the
  rest of the search.  Reseeded replicates were already excluded on the same
  grounds.  This changes which strategies later replicates sample when
  `adaptiveStart` and `tree = ` are combined; searches using either alone are
  unaffected.  `attr(, "strategy_diagnostics")$attempts` now counts only
  replicates that actually chose a strategy, so it can sum to fewer than the
  number of replicates completed, and at `verbosity = 2` the per-replicate
  `Strategy:` line is omitted for warm-started replicates rather than naming an
  arm that was never pulled.  The bandit's decay-on-improvement is deliberately
  *not* subject to this exclusion: it measures how stale the accumulated
  evidence is, not which arm ran, so it still fires when a warm-started
  replicate improves the best score.  That widens it slightly — a reseeded
  replicate that improved the score previously skipped the decay and now
  triggers it; `TS_POOL_RESEED` is off by default, so that path is dormant.

- `MaximizeParsimony(tree = )` no longer fails with "argument is of length
  zero" on some *valid* unrooted starting trees (e.g. `ape::unroot(ape::rtree())`).
  The bifurcating-tree check conflated "needs resolving" with "needs rooting":
  a valid unrooted binary tree has `nrow(edge) == 2 * NTip - 3`, which fails
  that check exactly as a genuine polytomy would, so it was passed to
  `MakeTreeBinary()`, which misread the unrooted root's legitimate degree-3
  trifurcation as a polytomy and corrupted the tree. Unrooted starts are now
  rooted (arbitrarily, on their first tip) before the bifurcating check runs.

- `MaximizeParsimony()` now contracts zero-length (unsupported) branches into
  polytomies by default (`collapse = TRUE`), deduplicating the returned trees on
  the resulting collapsed topologies, à la TNT's "collapse zero-length
  branches".  `n_topologies` counts distinct collapsed topologies, comparable
  across programs; this avoids reporting unsupported groupings and stops a single
  soft polytomy from inflating the apparent number of optimal trees by orders of
  magnitude.  Pass `collapse = FALSE` to recover fully-resolved trees (one
  arbitrary resolution per distinct topology).  Collapsed trees are returned
  rooted on the first leaf (a deterministic convention, as in TNT).  Under a
  topological `constraint`, the enforced splits are protected from collapse, so
  an enforced-but-unsupported clade stays visible while unsupported
  non-constraint branches still collapse.

- Fixed an internal inconsistency in the returned MPT set: the post-search MPT
  enumeration now deduplicates on collapsed topology, matching the main search
  loop, instead of keeping every resolved variant of one collapsed topology.
  On datasets with no zero-length branches (e.g. fully-resolved matrices) the
  returned set is unchanged; on datasets with soft polytomies it is smaller and
  internally consistent.

- `MaximizeParsimony(strategy = "thorough")` now adds 2 drift cycles and 5 Wagner
  starts, making `"thorough"` a higher-effort tier that trades wall-clock for more
  exhaustive search.  The drift cycles recover equal-score most parsimonious trees
  that sit on TBR-disconnected islands, which pure restarts rarely both seed: on a
  two-island exemplar (Zhu et al. 2013) two-island recovery rises 0.73 -> 0.95 over
  30 seeds.  A powered anytime study on 20 MorphoBank training matrices (65-120
  tips) quantifies the cost: drift is per-replicate overhead, so at a fixed budget
  `"thorough"` completes fewer replicates and reaches the optimum somewhat less
  reliably than before on the general pool -- it therefore needs a correspondingly
  larger replicate budget / more time to converge (see `maxReplicates`).  Choose
  `"thorough"` when tree quality and set-completeness matter more than turnaround.
- `strategy = "intensive"` is now a deprecated alias of `"thorough"`: the extra
  Wagner starts that once distinguished it are folded into `"thorough"`, so the two
  are identical.  Existing calls continue to work.

- **MorphyLib removed.** The Morphy Phylogenetic Library (C/C++) has
  been dropped; all parsimony scoring now runs through the native C++ kernel,
  which implements the Brazeau, Guillerme & Smith (2019) inapplicable-state
  algorithm correctly — including ambiguous-with-inapplicable tokens such as
  `{1-}`, which MorphyLib scored incorrectly.

- **`concavity` argument for `PrepareData()`.**  Implied-weights and profile
  searches with the custom-search functions no longer need a hand-written
  scorer: pass `concavity = k` (a finite constant) for implied weights,
  `concavity = "profile"` for profile parsimony, or the default `Inf` for
  equal weights, when preparing `dataset`; the default `TreeScorer`,
  [`EdgeListScore()`], honours it automatically.  (Adapted from the parallel
  T-200 work, PR #216.)

- **Custom-search scoring layer renamed** to drop the now-meaningless "Morphy"
  branding.  `PhyDat2Morphy()` → `PrepareData()`; `UnloadMorphy()` →
  `ReleaseData()`; `is.morphyPtr()` → `is.ParsimonyData()`; `SingleCharMorphy()`
  → `SingleCharData()`; `MorphyLength()` → `EdgeListScore()`;
  `MorphyTreeLength()` → `TreeScore()`; `MorphyBootstrap()` → `BootstrapTree()`;
  `RandomMorphyTree()` → `RandomPostorderTree()`.  The old names remain as
  deprecated aliases and will be removed in a future release.  `PrepareData()`
  returns a lightweight, garbage-collected `ParsimonyData` object; only the
  default `"inapplicable"` gap treatment is supported (recode data for the
  missing/extra-state treatments).

- **`TreeSearch()`, `Ratchet()` and `Jackknife()` no longer prepare `dataset`
  for you.**  Prepare it yourself -- typically with `PrepareData()` -- before
  calling these functions; a custom `TreeScorer` may take any `dataset` it
  likes (e.g. a raw `phyDat` object).  The `InitializeData` and `CleanUpData`
  arguments that formerly did this automatically are deprecated (a warning is
  issued if supplied) and will be removed in a future release: with scoring
  now handled by plain R data structures rather than external Morphy
  pointers, there is nothing left to initialize or destroy on the framework's
  behalf.  If your own `TreeScorer` holds an external resource that needs
  releasing, use your own `on.exit()`.

- `Jackknife()` and `BootstrapTree()` (formerly `MorphyBootstrap()`) now
  resample characters natively, scoring the resampled weights through the
  native kernel rather than by mutating a MorphyLib object — fixing a case
  where resampled weights could be silently ignored.

- The low-level MorphyLib bindings (`mpl_*()`), together with the
  `MorphyWeights()`, `SetMorphyWeights()`, `GapHandler()`, `MorphyErrorCheck()`,
  `GetMorphyLength()` and `C_MorphyLength()` helpers and the
  `summary.morphyPtr()` method, have been removed.

- `MaximizeParsimony()` results now carry a `candidates_evaluated` attribute:
  the number of TBR/SPR-class rearrangements examined during a single-threaded
  search (the analogue of TNT's "rearrangements examined"), for diagnosing
  search efficiency.

- New `SearchControl()` option `stallEscalateFactor` (default `1`, disabled):
  when a driven search stalls, escalate ratchet perturbation strength for
  subsequent replicates so the search adapts to a difficult dataset at runtime.

- Faster driven search: per-clip allocation churn in the TBR kernel removed
  (reusable scratch buffers and an open-addressed rerooting de-duplication
  table), and
  the debug-only topology validation no longer runs in release builds.

- Further driven-search speedup: the exact directional insertion edge-set
  computation now reuses caller-owned scratch and skips its per-clip zero-fill
  (under a write-before-read invariant, debug-asserted), saving up to ~16% wall
  on large datasets where the `O(n_node * words)` zero-fill dominated.  Search
  results are bit-identical (score and `candidates_evaluated` unchanged).

- HSJ (Hopkins & St John 2021) scoring is now invariant to the
  arbitrary ordering of a `phyDat` object's `levels`.  Both the primary
  absent/present term and the secondary-character dissimilarity term
  previously depended on the internal token ordering, so the same dataset
  could score differently under different (equivalent) `levels`.

- Hierarchical-scoring helpers renamed to the package's `BigCamelCase`
  convention (the snake_case names are removed): `recode_hierarchy()` ->
  `RecodeHierarchy()`, `hierarchy_from_names()` -> `HierarchyFromNames()`,
  `validate_hierarchy()` -> `ValidateHierarchy()`, `hierarchy_chars()` ->
  `HierarchyChars()`, `hierarchy_controlling()` -> `HierarchyControlling()`.
  The internal C++-bridge helpers (`build_tip_labels()`,
  `hierarchy_to_blocks()`, `non_hierarchy_weights()`, `hsj_absent_state()`)
  are now private and no longer exported.

- `WideSample()` now dispatches to the appropriate Max-Min diversity (MMDP)
  solver from the `Coreset` package, choosing the tier automatically
  from `length(trees)`.

- New functions `LeastSquaresTree()` and `LeastSquaresFit()` search for, and
  fit branch lengths to, the tree that best matches a target distance matrix
  under a least-squares criterion, reusing the optimised C++ rearrangement
  kernel (NNI + SPR).  Ordinary (`method = "ols"`) and non-negative
  (`method = "nnls"`) least squares are supported, with optional
  Fitch-Margoliash (`weight = "fm"`) or custom weighting.
  `LeastSquaresFit()` mirrors `phangorn::nnls.tree()` but runs in the native
  kernel.


- `attr(dataset, "weight")` now accepts non-integer character weights.  The
  C++ scoring engine still stores `int` weights internally; fractional
  inputs are rescaled to integer with a configurable precision (default
  ~0.001, controlled by `getOption("TreeSearch.fractional.scale", 1260L)`).
  Previously, fractional weights were silently truncated at the Rcpp
  boundary (e.g. `c(0.5, 1.7)` became `c(0L, 1L)`, dropping 50% / 41% of
  the respective characters' contributions).  Integer weights pass
  through unchanged.  `TreeLength()` and other scores are returned in
  units of `steps * scale` when fractional weights are present; within-
  run ranking is unaffected.

- `LengthAdded()` removes a temporary warning guard that fired on datasets
  triggering the T-302 `qmApp` scalar-unwrap fix; regression tests now cover
  both the `qmApp` (T-302) and `qm` (commit e8b318c3) scalar-unwrap paths,
  confirming all deltas are non-negative and match independent computation.

- `ExpectedLength()`'s internal cache no longer collides across trees:
  its key omitted any tree-derived component, so scoring two different trees
  against the same dataset with the same `nRelabel` could silently return one
  tree's cached result for the other, corrupting `rhi` -- a published
  statistic -- returned by `Consistency()`.

- `.SortTokens()` (used internally by `ExpectedLength()`) no longer rewrites
  a partial-ambiguity token (e.g. `(01)`) as full ambiguity when the
  dataset's contrast holds other ambiguous tokens (e.g. `?`) that are not
  present in the character being processed, another silent corruption of
  `rhi`.

- `Consistency()` now always returns a matrix, even for a dataset that
  compresses to a single character pattern; it previously returned a bare
  numeric vector in that case, breaking `[, "ci"]`-style column access.

- `Consistency()`'s documentation now states explicitly when its `ci`, `ri`,
  `rc` and `rhi` columns are `NaN` (constant, autapomorphic and
  zero-null-homoplasy characters respectively); the values themselves are
  unchanged.
  
- `QuartetResolution()` no longer errors on a tree in which the four focal
  tips form an unresolved (star) quartet -- reachable from
  `MaximizeParsimony(collapse = TRUE)` output, the default since 2026-06-24.
  Such a tree now contributes `NA` rather than raising a `vapply` error.

- `WideSample()` fixes three bugs in tree-set handling: the `effort = 1`
  (`FarFirst()`) tier returned trees in farthest-first selection order rather
  than the ascending input order its own comment described; `FarFirst()` was
  called with a mix of positional and named arguments, fragile to any future
  change to the function's argument order; and the `firstHit` attribute
  (a per-*stage* tally computed from tree names) was copied onto the
  subsetted output unchanged, so it continued to describe the pre-subset
  input rather than the trees actually returned -- `firstHit` is now dropped
  when subsetting, rather than carried over stale; call `WhenFirstHit()` on
  the result to recompute it. Other attributes (`score`, `hits_to_best`,
  etc.) are unaffected.

- `BootstrapTree()` no longer risks the classic `sample()` length-1 vector
  trap, in which a single remaining character index `k` would be sampled as
  `sample(1:k, ...)` rather than always returning `k`.

- `WhenFirstHit()`'s stage-name pattern is now anchored, so a tree or
  replicate name that merely contains a stage pattern (rather than matching
  it exactly) no longer produces a spurious, garbled stage label.

- `TaxonInfluence()`'s distance-weighted mean no longer assumes a fixed
  dimension ordering from a user-supplied `Distance` function; the returned
  matrix's shape is now checked and normalized explicitly.
- `ClusterStrings()` no longer crashes when the best clustering contains a
  singleton cluster, no longer omits the documented `silhouette` attribute
  when few unique strings are supplied, and its "no structure" branch now
  returns the documented per-element cluster-assignment vector rather than a
  bare scalar `1`.  Its internal call to `cluster::pam()` now passes the
  Levenshtein distance matrix via `as.dist()`, so it is treated as a
  dissimilarity rather than clustered on Euclidean distance between its
  rows; **silhouette scores and, in some cases, cluster assignments for the
  `pam` method may change** to more accurately reflect string similarity.

- `ParsSim()` now errors clearly, instead of silently corrupting the Fitch
  score, if asked to simulate a character with 32 or more states -- the
  internal bit-set representation of state sets overflows a 32-bit integer
  beyond that.  It also errors clearly, instead of an opaque
  `sample.int()` failure, if a tree lacks the structure to host the number
  of requested states for a character.  Simulation with `nExtraSteps > 0`
  is also faster, as the redundant saturation scan previously performed
  again on every character at return now reuses the result already
  computed during the step-placement loop.
  
- `ClusteringConcordance(normalize = TRUE)` now chance-corrects large trees,
  which it previously left uncorrected while still describing the result as
  corrected.  The expected mutual information that sets the zero point was
  accumulated by a recurrence over the hypergeometric distribution of cell
  overlaps, seeded at the smallest overlap the marginals allow.  That
  probability sinks below the smallest representable double once a character
  scores about 1080 tips, and the recurrence being multiplicative, every later
  term then stayed zero: `expected_mi()` returned exactly 0 where the seed
  vanished for every block of the character, and a silently truncated sum --
  as little as a quarter of the true value -- where it vanished for some.  The recurrence is now anchored at the
  mode of the distribution, whose probability is the largest of at most
  `N + 1` values summing to one and so is always representable.  The threshold
  is a property of the tips each character scores rather than of the tree, and
  only marginals close to even reach it: of 600 random partitions, none below
  1200 items was affected, 13 of 60 at 1500 items, and 30 of 60 at 3000.
  Values that were already correct are unchanged, to of order 1e-11 relative.

- `expected_mi()` now checks that `ni` gives exactly two block sizes, as its
  documentation always required.  A shorter vector was read past its end, and
  the arbitrary values that produced could index the log-factorial lookup
  table out of bounds and crash the session.

- `QuartetConcordance()`'s counting kernel now rejects a negative character
  state code rather than indexing its count buffers out of bounds.  State
  codes generated by the package are always positive, so no result changes.

- `TBRMoves()` now lists the complete TBR neighbourhood.  It never
  bisected the edge leading to the first-labelled tip, so every rearrangement
  that relocates that tip was missing -- around a dozen trees on a typical
  eleven-leaf tree, and enough that the output was not even a superset of
  `SPRMoves()`, which TBR contains by definition.  `TBRMoves()`
  therefore returns more trees than before, and any count derived from it will
  rise.  `MaximizeParsimony()` is unaffected: it drives a separate enumerator
  that has always swept the root edge.
- `Ratchet(stopAtScore = )` no longer returns a tree whose independently
  recomputed score disagrees with its `"score"` attribute.  Its early-exit
  paths -- meeting the target score during search, or already meeting it on
  entry -- skipped the bookkeeping that the return value depends on, so the
  *input* tree could be returned carrying the *improved* score.
  `returnAll = TRUE` no longer errors ("No trees!?") when the target score is
  met during search, and `MultiRatchet()` no longer errors when a starting
  tree already meets `stopAtScore`.

- The cache behind `ClusteringConcordance(normalize = TRUE)` keyed partitions
  on block sizes narrowed to 16 bits, so two partitions whose block sizes
  differed by a multiple of 65536 shared an entry and the second was given the
  first one's expected mutual information.  Reaching this needed a tree of at
  least 65536 tips, so no published result is affected; keys now span the full
  range of an integer, which rules the collision out rather than making it
  unlikely.

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

- `ScoreSpectrum()`: Chao1-style landscape coverage estimator.  Treats
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
- `LeastSquaresTree()` and `LeastSquaresFit()` search for, and
  fit branch lengths to, the tree that best matches a target distance matrix
  under a least-squares criterion, reusing the optimised C++ rearrangement
  kernel (NNI + SPR).  Ordinary (`method = "ols"`) and non-negative
  (`method = "nnls"`) least squares are supported, with optional
  Fitch-Margoliash (`weight = "fm"`) or custom weighting.  This provides the
  topology-search step of Lapointe & Cucumel's (1997) average consensus
  procedure; `LeastSquaresFit()` mirrors `phangorn::nnls.tree()` but runs in
  the native kernel.
- `PaintCharacters()` colours each character in a morphological
  dataset by the hue of the tree edges it most concordantly supports, using
  `ConcordanceTable()` MI weights averaged in CIELAB colour space.  Pairs with
  `TreeTools::PaintTree()` to visually map characters to clades.
- `attr(dataset, "weight")` now accepts non-integer character weights.  The
  C++ scoring engine still stores `int` weights internally; fractional
  inputs are rescaled to integer with a configurable precision (default
  ~0.001, controlled by `getOption("TreeSearch.fractional.scale", 1260L)`).
  Previously, fractional weights were silently truncated at the Rcpp
  boundary (e.g. `c(0.5, 1.7)` became `c(0L, 1L)`, dropping 50% / 41% of
  the respective characters' contributions).  Integer weights pass
  through unchanged.  `TreeLength()` and other scores are returned in
  units of `steps * scale` when fractional weights are present; within-
  run ranking is unaffected.

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

- `LengthAdded()` no longer errors on datasets whose contrast matrix contains
  zero-sum rows for tokens that are declared in the SYMBOLS list but not used
  by any taxon in the character being scored (#294).

- `LengthAdded()` no longer returns negative values when multiple rows of the
  contrast matrix satisfy the fully-ambiguous applicable condition (e.g.
  datasets with ~19 taxa and certain character structures); the first matching
  row is now used consistently (#302).

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
