# Two-island MPT regression fixture (Zhu 2013)

Files backing `dev/benchmarks/bench_thorough_two_islands.R`.

| file | what |
|------|------|
| `zhu2013_2island.nex` | Original Zhu et al. (2013) *Entelognathus* matrix, 75 taxa × 253 chars, extracted from the Nature SI. **This is the un-cleaned matrix** — it is *not* the bundled `inapplicable.phyData[["Zhu2013"]]`, which was later edited/improved and no longer exhibits the two-island pathology. |
| `zhu2013_2island_tnt.tre` | TNT MPTs (`mult = replic 100 tbr; bbreak = tbr;`, gaps recoded to missing): 224 trees at length 598. Read with `TreeTools::ReadTntTree(tipLabels = taxa)`. |
| `zhu2013_2island_island2.tre` | Golden reference: the rare **16-tree** disconnected island (Newick, full taxon labels). |
| `zhu2013_2island_main.tre` | The 256 TS-reachable main-island MPTs (Newick), for the "main island recovered" sanity check. |

## The pathology

Under gaps-as-missing Fitch, the length-598 MPTs split across **two TBR-plateau
islands that are disconnected at equal score**:

* a large **main** island every search reaches, and
* a small **island 2** of exactly **16 trees** with a narrow basin. Seeding
  TreeSearch from any one of the 16 and running a pure accept-equal TBR walk
  closes over exactly those 16 and *none* of the main island — proof they are a
  separate connected component. TS reaches island 2 only if a replicate's random
  start lands in its basin; TNT finds it, a default TS run frequently does not.

Comparison (2026-06-24): TS 256, TNT 224, **208 shared, 48 TS-only, 16 TNT-only**,
union 272. All 16 island-2 trees score exactly 598 under `TreeLength` — genuine
MPTs, not a scoring artefact.

## Why it's a regression guard

`"thorough"` exists for exactly this completeness problem, so almost every run
should recover **both** islands. The benchmark runs `thorough` over many seeds
and fails if fewer than `TS_BENCH_THRESHOLD` (default 90%) capture both. It
overrides only `poolMaxSize` and `enumTimeFraction` (so a miss reflects *seeding*,
not pool-capacity or enumeration-time artefacts); the thorough preset still
governs the search heuristics.

## Regenerating the reference

To rebuild `zhu2013_2island_ref.rds` from scratch:

1. Run TreeSearch (`inapplicable = "missing"`, `concavity = Inf`) and TNT
   (`mult = replic 100 tbr; bbreak = tbr;`, gaps→missing) on
   `zhu2013_2island.nex`; keep the unique length-598 topologies from each.
2. Canonicalise every tree as `write.tree(SortTree(RootTree(t, outgroup)))`.
3. `island2` = the TNT-only set; confirm it is one disconnected component by
   seeding TreeSearch from any one of them with ratchet/sector/drift/fuse
   disabled (`accept_equal` TBR walk only) and checking the closure contains
   exactly those trees and none of the main island.
