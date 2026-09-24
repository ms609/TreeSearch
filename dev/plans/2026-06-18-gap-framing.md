# Framing the TS-vs-TNT wall-clock gap (2026-06-18)

Answers the directive: *before VTune, frame the right question — is TNT really
running 10× more iter/sec? Is each rearrangement of similar complexity? Or is
the gap high-order (buffering/caching → fewer operations), not line-level?*

**HEAD 25e35be7** (post Wagner directional fix 2b299e4b + `unrooted=TRUE`
default). Fresh `-O2` instrumented build `.agent-p0`. ALL prior profiling
(Rounds 1–4) is **stale** — it predates both changes. Drivers:
`dev/profiling/drivers/framing.R`, `perclip.R`.

## TL;DR — the "10×" premise is false

The matched-effort, same-machine decomposition is an identity:

```
wall_ratio(TS/TNT) = efficiency(ts_cand / tnt_rearr) × throughput(tnt_rate / ts_rate)
```

Counts are bitness-portable (local 32-bit TNT is valid for efficiency); the
throughput factor uses 32-bit TNT so it is a LOWER bound on the true per-second
gap.

| dataset | tips | gapB | efficiency | throughput (32-bit) | local wall ratio |
|---------|------|------|-----------|---------------------|------------------|
| Wortley2006 | 37 | **0** | 0.49 | 1.36 | **0.68** |
| Zanol2014   | 74 | **0** | 0.95 | 2.30 | 2.17 |
| Zhu2013     | 75 | **0** | 1.08 | 1.84 | 1.98 |
| Giles2015   | 78 | **0** | 0.42 | 1.81 | **0.75** |

1. **gapB = 0 everywhere** — TS reaches TNT's exact optimum (matched stopping;
   quality parity holds on the fresh build).
2. **Efficiency 0.42–1.08** — TS examines comparable-or-FEWER candidates than
   TNT. Full reversal from the buggy era (1.4–1.9, `headtohead_phase0.csv`).
   **Search strategy is no longer the lever** — the cost-bug + unrooted-TBR
   fixes made TS candidate-efficient.
3. **Throughput 1.36–2.30× (32-bit, same machine) — NOT 10×.** Each TS
   candidate is of *similar complexity* to a TNT rearrangement (within ~2×).
   64-bit TNT would widen this to perhaps 2–4× (needs a local 64-bit TNT to
   pin down).
4. On Wortley/Giles, TS is **faster than local 32-bit TNT** outright.

**Where did "10–16×" come from?** It conflated TNT-64bit-on-EPYC wall with
TS-run-to-a-wasteful-full-thorough-budget — not a like-for-like throughput
measurement. The like-for-like, matched-50-start, same-machine gap is
throughput-only and modest.

## The high-order finding (the buffering/caching the directive asked about)

Per-clip economy (`perclip.R`, TBR-to-convergence):

| dataset | cand/clip | ns/cand | µs/clip |
|---------|-----------|---------|---------|
| Wortley2006 | 198 | 102 | 20 |
| Zhu2013 | 590 | 46 | 28 |
| Giles2015 | 763 | 42 | 31 |

cand/clip is large (200–760) so per-clip overhead is amortized, BUT ns/cand
(42–102) materially exceeds the bare bounded scan (~18–23 ns floor) ⇒ **roughly
half the per-candidate time on the 75–78-tip sets is per-clip / per-reroot
overhead** — the per-clip state TNT maintains *incrementally* (Goloboff 1996)
but TS *recomputes every clip*:

- `compute_insertion_edge_sets` — **NEW** (Wagner fix): O(n·chars) directional
  up-pass + **3 heap allocs/clip** (`up`, `pre`, `edge_set.assign`).
- `compute_from_above` — O(subtree·chars); already Tier-1 (no per-clip alloc).
- `vroot_cache` rebuild — a per-clip **memcpy re-indexing of `edge_set_buf`**
  (ts_tbr.cpp:1487–1492); the scan could index `edge_set_buf[below]` directly.
- unrooted reroot machinery — `fitch_join_states` + `fast_hash` + dedup per
  reroot (NEW-ish: `unrooted=TRUE` default).

**This is the answer:** the residual gap is not line-level scoring (the AVX2
scan is at the compiler limit — Round 3, re-confirmed by the ~20 ns floor) and
not search strategy (efficiency at parity). It is **per-clip state
recomputation that TNT amortizes via incremental maintenance.**

## Levers (in increasing order of payoff and risk)

- **Quick, behaviour-neutral (recently-added code):** hoist
  `compute_insertion_edge_sets`'s `up`/`pre` to caller-owned buffers (emutls
  lesson — plain locals, NOT `thread_local`; see T-S3d); drop the `vroot_cache`
  memcpy by indexing `edge_set_buf` directly in the reroot scan. Verify
  score-identical + wall delta (A/B). *Pending VTune attribution of which is
  material — see Round 5 next-step.*
- **High-order (the real throughput lever, harder):** maintain per-clip state
  (`from_above` / edge sets) incrementally across clips rather than recomputing
  O(n·chars) per clip — what TNT does. Significant change; correctness-critical;
  scope separately.

## Status

Framing COMPLETE and durable (`dev/profiling/log.md` Round 5). VTune attribution
of the per-clip overhead components is the next step (read-only). Any code
change is behaviour-neutral-only and not to be committed without sign-off.
