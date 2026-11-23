# x4 reroot-batch software-prefetch parity (`TS_TBR_PREFETCH`) — 2026-07-16

**Status:** BUILT + byte-identity-validated (EW+IW) + local-flat. **Landed default-OFF**
(opt-in `TS_TBR_PREFETCH=1`) as a validated building block; the win/wash is settled on
Hamilton EPYC (raw-TBR + `MaximizeParsimony` wall, min-of-runs) before any default flip.

## Where it came from

Post-packing kernel-lever hunt (workflow `kernel-reroot-simd-hunt`, map→judge→adversarial-verify).
Verdict: **the per-move kernel is near its exact floor.** The per-candidate cost floor is
**memory-latency-bound on the `vroot_cache` gather** (~1.3 cyc/word measured vs ~0.3 cyc/word for
L1-resident AVX2 compute — ~4× of the cost is L2/L3 latency; matches the 0.23→0.44 ns/word growth
in `kernel-speed-baseline.md`). So the only lever that touches the true floor is **software
prefetch of the single-pass gather stream** — NOT wider SIMD (AVX2 is already under-fed: packing
cut per-block `n_states` to ~2–3 < the 4-lane width, so the reduce runs mostly scalar-tail), NOT
cache-blocking (no reuse to block on), NOT compute fusion (attacks the wrong bottleneck). All of
x8-batch / AVX-512 / multi-block-fusion / cache-block / vectorised-popcount / all three lever-6
incremental-precompute variants were **STOP** (measured washes or zero headroom).

## The lever

The dominant reroot path (77–81% of candidates) evaluates candidates through the 4-wide flat
batch kernels `fitch_indirect_cached_flat_x4` / `indirect_iw_cached_flat_x4`
(+ NA variants), each reading 4 scattered `vroot_cache[ei*tw]` rows. The **scalar** reroot path
already software-prefetches its next row (ts_tbr.cpp scalar branch); the **x4** batches issued 4
cold demand loads per batch with **no hint**, so the 4 row-heads at each batch boundary are exposed
L2/L3 latency the in-batch 4-way MLP cannot hide (the cutoff early-break is a control dependency the
OOO window can't run past).

Fix: prefetch the *next* batch's 4 `vroot_cache` row-heads during the current batch's reduce, in
both the EW x4 and IW x4 caller loops — parity with the scalar path.

- **EW x4** consumes `kept_ei` (pre-filtered ascending regraft indices) → prefetch `kept_ei[ki..ki+3]`
  after the current batch advances `ki`. Trivial.
- **IW x4** re-scans `main_edges` inline, but its skip predicate is **provably identical** to the
  `kept_ei` construction (`use_collapsed == !collapsed.empty()`, ts_tbr.cpp:1532), so `kept_ei` is
  the same candidate sequence in the same order. A parallel cursor `iw_pf` (advanced by the batch
  count) indexes the upcoming rows in `kept_ei` for the prefetch — **the candidate loop itself is
  untouched**.

Purely additive and flag-gated: with the flag off, the only added statement is a dead `iw_pf`
increment, so **flag-OFF is byte-identical to the prior committed code by construction**, and a
prefetch hint can never change a score or trajectory anyway (non-faulting, no program-visible state).

## Evidence

- **Byte-identity (must):** `ts_tbr_diagnostics` ON vs OFF on project970, both EW and IW (k=10),
  identical `score` + total `n_candidates_evaluated` (38.9M / 44.8M) + `n_clips_tried` + final edge.
  (`scratchpad/pf_validate.R`.)
- **Local A/B (Windows/MinGW — NOT decisive; workflow flagged local ≠ EPYC cache/prefetcher):**
  project970 60.0→60.2 ns/cand (flat); project510 (biggest buffer, ~1.4 MB) 57.1→57.0 (flat).
  **No regression, no local win** — consistent with the prediction that the win only appears in the
  EPYC L3-pressure regime (or is covered by EPYC's HW stride prefetcher → wash).

## Open gate (decides the default flip)

Hamilton EPYC A/B, **separate lib dir** (NOT `packbuild/lib` — the in-flight 5432 job `17887591`
uses it): raw-`tbr_search` wall (min-of-runs) + full `MaximizeParsimony` wall, prefetch ON vs OFF,
on project970 / project510 / project5432. Honest expectation: **~1.03–1.06× best case on the reroot
body, possibly a wash** if EPYC's HW prefetcher already covers the ascending-`ei` stream. If it
washes, the exact kernel is at its floor and remaining wall gains are search-policy/over-search
(Mission A / stopping), not the per-move kernel.
