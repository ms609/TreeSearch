# Kernel re-slicing: ranked design to cut the ~47 ns/candidate (9-agent workflow, 2026-07-14)

Workflow `kernel-reslice-design` (Map→Design→Judge, 9 agents, 960k tok). Cost model + 4 re-slicings
developed in parallel, adversarially ranked. Baseline (real data, `kernel-speed-baseline.md`):
project970 76 ns, project2668 41 ns, project4327 18 ns per candidate. Split ≈ **17 ns per-clip
precompute + ~34 ns per-candidate score** (score = cache-cold gather+reduce over `total_words`).

## Ranked levers

| # | lever | exact? | measured/est speedup | verdict |
|---|---|---|---|---|
| 1 | **Incremental edge-set (lever-6, Scheme 1)** — patch each clip's edge_set from the whole-tree base instead of O(n·chars) from scratch | EXACT | **1.36–1.39× @482t (MEASURED, bit-identical)** | **PURSUE #1** — already BUILT on `claude/lever6-edgeset`; hardware-robust (instruction-count, not cache); caps ~1.4× (attacks only the 17 ns precompute) |
| 2 | **Block-major edge_set layout** — reorder buffer so block-0 of all nodes is contiguous (≈69 KB, L2-resident) | EXACT | est ~1.1–1.16× blended; **hardware-CONTINGENT** (may be ~5% on the authoritative EPYC 512 KB L2 if block-0 already resident) | PURSUE #2, composes with #1 (different third: gather) |
| 3 | **Char-packed early-bail** (TNT-style fields, char-sequential bail) | EXACT | unique increment only **~2–4 ns (~1.05–1.10×)** at 482t/9-state | MARGINAL — full rewrite; density-WORSE for non-power-of-2 states (9-state→16-bit field, 1.78× inflation); the only lever whose *words-touched* reduction could scale toward 1.9 ns, but not at this corpus's state counts |
| 4 | Coarsened-state (g-bucket) sound lower-bound screen on block-0 | sound screen | <5%, SPR-only | near-DEAD — vanishes on EPYC; risks the reroot MLP |

## The honest ceiling (the answer to "the kernel CAN be faster")

- **Exact re-slicing caps at ~1.4–1.5×** (51→~34 ns): incremental (~1.36×, built) + block-major
  (~1.1–1.2×, hardware-contingent), nominally additive (different thirds), bounded by the
  whole-edge-set-machinery ~1.5× ceiling (`l3b-footprint-482.md`).
- **None approach TNT's ~1.9 ns (a ~24× floor).** The per-candidate reduce+gather is O(total_words),
  and `total_words` is irreducible for uniform 5–10-state data (bit-slicing is already minimal-bits;
  char-packing rounds *up* and is denser only for power-of-2 low-state data). To touch far fewer
  words per candidate you need **incremental reoptimisation ACROSS candidates** (lever-b, judged
  infeasible on morphological data — no cross-clip locality) or a **sound bulk-reject screen** (the
  Goloboff union construct — sound only with true finals, and it cuts candidate COUNT not
  per-candidate cost). TNT's 1.9 ns per "rearrangement examined" also partly reflects cheaply
  screen-rejected candidates in its denominator (not purely a faster scored-candidate kernel), and
  TS's per-*useful*-candidate throughput already matches TNT in practice (head-to-head: TS reaches
  the optimum faster on 5/6 large datasets).

## Instrument caveat (blocking, advisor catch)

`TS_IW_TIMING` brackets the **scored-candidate body** (the ~34 ns gather), NOT the per-clip
precompute — so it shows incremental at ~0% (a FALSE KILL). Incremental's metric is **raw-TBR WALL
(min-of-runs)** or a per-clip precompute chrono. `TS_IW_TIMING`/`kern_real.R` is the right metric for
block-major / char-packed / screen (which attack the scoring body).

## Recommended build order (if pursuing the exact ~1.4–1.5×)

1. Cherry-pick lever-6 onto tip (resolve vs lever-7 monomorphized scan); gate `TS_L3B_INCREMENTAL`
   default-OFF + feature-gate to plain directional EW/IW (fallback for NA/sector/constraint/collapse).
2. Byte-identity gate: per-clip full-array oracle + C-vs-A `MaximizeParsimony` (score + n_evaluated +
   per-pass) on EW + one IW + one NA; positive control (skip Phase-A off-path sibling seeding → must
   fail oracle). Measure raw-TBR WALL @482t on Hamilton EPYC.
3. Then block-major (separate branch, composes): node-ID-order streaming writes; `_bm` scorer + x4_bm;
   gate `use_flat && !use_iw && !has_na`. A/B the composition (may negatively interact with
   incremental's patch/undo → 3 strided writes; mitigated by node-ID-order).
