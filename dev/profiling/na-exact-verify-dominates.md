# The native-NA wall is `exact_verify_sweep`, not per-clip overhead

**Question (2026-07-28).** The same matrix costs 74–151× more per ratchet re-opt cycle
with inapplicables kept than recoded `-`→`?`, on **identical pattern and character
counts** (Zanol2014 210 patterns both ways, Giles2015 236 both ways — only the state
count drops by one), while evaluating only ~1.9× the candidates. So the wall is not
candidate scoring. Where is it?

The standing answer (`bgs-component-map`, 2026-06) was: *"native-NA TBR wall is dominated
by PER-CLIP overhead (vroot/below_actives builds + NA rescores ≈ 98%), NOT candidate
scoring. `below_actives_cache` O(n_main×n_blocks)/clip is the obvious first suspect."*

**That is wrong, and this note replaces it.**

## Instrument

`TS_NA_TIMING` (default off, diagnostic only) accumulates disjoint brackets onto the
`DataSet`, so one R-level call yields whole-search totals:

| bracket | what |
|---|---|
| `evs` | `exact_verify_sweep` — the NA-only unrooted-TBR certification |
| `below` | `below_actives_cache` build, per clip (NA-only) |
| `vroot` | `compute_from_above` + `vroot_cache` build, per clip |
| `accept` | accept-path rescore (NA dirty passes, or `full_rescore` for a reroot) |
| `rest` | remainder = candidate scan + the down/uppasses |

The brackets are disjoint, so a **positive** `rest` is the check that nothing is
double-counted. Exposed on both `ts_tbr_search` and `ts_ratchet_search`; the latter
matters because one `DataSet` spans all ratchet cycles, which is the only way to see the
production `evs_false_cache` hit rate — a single `ts_tbr_search` call always begins with
an empty cache and so understates it.

Three regimes, because they answer different questions: **climb** (pectinate start, many
converge→improve→re-descend rounds), **certify** (re-descend the converged tree — isolates
certification cost), **ratchet** (production shape). Recoded rows report zero for
`evs`/`below` by construction — both are `has_na`-gated. That is the control, not a
finding that they are small.

## Result: one function is 89–99.6% of it

Giles2015 (78t, 51% NA), IW k=10:

| phase | total | **evs** | below | vroot | accept | rest |
|---|---|---|---|---|---|---|
| climb | 4677 ms | **89.2%** | 0.1% | 0.3% | 2.2% | 8.2% |
| certify | 3952 ms | **99.6%** | 0.0% | 0.0% | 0.1% | 0.3% |
| ratchet | 16648 ms | **96.6%** | 0.0% | 0.1% | 0.7% | 2.6% |

* `below_actives_cache` — the recorded "obvious first suspect" — is **0.001 ms per build**
  and 0.0–0.1% of wall. Refuted outright.
* `vroot` + `compute_from_above` together: 0.1–0.3%.
* Accept-path NA rescores: 0.02 ms/call, 0.7–2.2%.
* A single `certify` call costs **3936 ms**. Per executed call: 695 ms (climb),
  3936 ms (certify), 1237 ms (ratchet).

**The `evs_false_cache` hit rate is 0%** — 0 hits across 6, 1 and 13 calls. The
memoization never fires in a search that is still moving: it stores only FALSE (true
optimum) verdicts, and a search that keeps improving never revisits the same
(topology, weighting-regime) pair. Fixing the cache is therefore not the lever.

Regime penalty, whole-call wall:

| matrix | climb | certify | ratchet |
|---|---|---|---|
| Zanol2014 (74t, 64% NA) | 13.03 s vs 0.18 s (72×) | 2.57 s vs 0.02 s | 68.26 s vs 0.51 s (**134×**) |
| Giles2015 (78t, 51% NA) | 3.65 s vs 0.21 s (17×) | 2.84 s vs 0.01 s | 42.74 s vs 0.42 s (**102×**) |

## Why it costs what it costs

`exact_verify_sweep` exists because the NA indirect scan is only approximate under
Brazeau's three-pass — the clipped subtree's internal count is attachment-dependent — so
apparent convergence is not convergence. To certify, it enumerates every non-root edge's
TBR neighbourhood: clip × {identity + fragment rerootings} × regraft edges, i.e. O(n³)
candidates, each scored via `apply_tbr_move` plus a rescore.

It already uses the fast incremental (3-seed dirty) path per candidate, so this is **not**
a constant-factor problem — it is exhaustive re-enumeration. And none of it appears in
`n_candidates_evaluated`, which is why the per-move kernel campaign could measure the NA
scan as "at-limit" while 90% of the wall sat next to it, uncounted.

## What could actually be done

**(a) Certify less often — the big one.** Nothing in `tbr_search` knows whether its caller
needs a *certified* optimum. The ratchet's perturb and re-opt phases do not: the next
perturbation moves the tree anyway. Only the result actually reported needs certifying. A
`TBRParams::certify_unrooted` flag (default true, cleared by ratchet/drift/sector callers)
would remove most of those calls. **This is quality-sensitive** — uncertified convergence
returns slightly worse trees — so it must be judged on floor attainment, not on wall
(see `reach-not-escape-rate`).

**(b) Prune the sweep — blocked, and we know exactly why.** Pruning needs a sound *lower*
bound. `dev/red-team/union-of-finals-bound-proof.md` proves the union-of-finals
approximation the engine deploys is a sound **upper** bound (it over-counts), the reverse
of the in-code claim — so it cannot be used to skip candidates. The same proof's step S3
shows the union of *true MPR* finals **is** a sound lower bound, and the engine does not
compute MPR finals. So: computing MPR finals would unlock exact pruning of the function
that is 90% of the native wall. That is the principled route, and the proof already
specifies what is needed.

An upper bound is still usable in the other direction: if `upper(cand) < best`, the
candidate is a guaranteed improver and can be applied without exact scoring. That
short-circuits the *finding* case (5 of 6 climb calls, 8 of 13 ratchet calls found an
improver) but not the *certifying* case, which is the expensive one.

## It is structural, not an NA-density or size effect

Four matrices spanning 7–64% inapplicables and 23–88 tips. `exact_verify_sweep` as a share
of `tbr_search` wall:

| matrix | tips | %NA | climb | certify | ratchet | ratchet regime penalty |
|---|---|---|---|---|---|---|
| Zanol2014 | 74 | 64 | 95.6% | 99.5% | 98.1% | 134× |
| Giles2015 | 78 | 51 | 89.4% | 99.7% | 98.1% | 101× |
| Dikow2009 | 88 | 7 | 95.3% | 99.7% | 98.1% | 99× |
| Vinther2008 | 23 | — | 94.8% | 98.7% | 97.2% | 24× |

**Pooled over all 12 native rows: `evs` = 97.7% of `tbr_search` wall; per-clip scaffolding
(`below` + `vroot`) = 0.08%; accept path = 0.4%.** Minimum remainder 0.23%, so the brackets
are disjoint and nothing is double-counted.

Dikow2009 at **7% inapplicables** behaves like Zanol2014 at 64%, and 23-tip Vinther2008
behaves like 88-tip Dikow2009. The cost is not proportional to how much NA data there is —
it is the price of *certifying* at all, paid whenever `has_na` is true.

Cache hit rate pooled: **2.1%** (4 hits / 188 calls). Hits appear only in the ratchet
regime (where one `DataSet` spans cycles) and even there reach 6% at best.

Cost per executed call scales with the matrix: 31–43 ms (Vinther 23t) to 1.1–3.9 s
(74–88t). The `certify` call — proving no improver exists — is the worst case in every
row, because that is the one that cannot exit early.
