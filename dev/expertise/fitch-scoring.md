# Fitch Scoring — Design Notes & Proven Invariants

Reference for agents working on `ts_fitch.h/.cpp`, `ts_fitch_na.h`,
`ts_fitch_na_incr.h`, or the search modules that call them.

## Incremental uppass correctness (standard Fitch)

The incremental uppass (`fitch_incremental_uppass`) uses a dirty-flag
propagation scheme that does **not** explicitly revisit every node whose
prelim changed during the incremental downpass. Only nodes whose
*ancestor's final* changed are recomputed.

This looks like it could miss updates when the downpass stops before
root (prelim stabilises at some intermediate node N). Nodes between
`clip_ancestor` and N have changed prelims but their ancestors' finals
are unchanged, so the dirty-flag scheme skips them.

**This is provably correct for standard (non-NA) Fitch blocks.**

### Proof sketch

When the downpass stops at node N, `fitch(M_new, S) = fitch(M_old, S)`
where M is N's child on the downpass path and S is the sibling.

**Case 1 — both intersection-type:** `M_old ∩ S = M_new ∩ S = P`.
Then N_final ⊆ P ⊆ M_old and N_final ⊆ P ⊆ M_new. So
`uppass(N_final, M_old) = N_final ∩ M_old = N_final` and likewise for
M_new. Finals are identical.

**Case 2 — both union-type:** `M_old ∪ S = M_new ∪ S` with
`M_old ∩ S = ∅` and `M_new ∩ S = ∅`. Since the unions are equal and
both M sets are disjoint from S, `M_old = M_new`. No change.

**Case 3 — mixed types:** Intersection equals union only if both
operands are identical and the set is trivial. Not reachable in
practice (would require empty state sets).

The argument applies per-character (per bit position), so it holds
for packed 64-bit representations.

### Consequence

No code change needed. The dirty-flag scheme is an optimisation that
happens to be exact for standard Fitch, not just a heuristic.

---

## NA uppass `children_app` staleness

The NA-aware incremental uppass (`fitch_na_incremental_uppass`) has a
**theoretical staleness issue** that does NOT affect standard blocks.

The NA uppass formula at internal nodes uses:

```cpp
uint64_t children_app = 0;
for (int s = 1; s < k; ++s)
    children_app |= (tree.prelim[left + s] | tree.prelim[right + s]);
```

This `children_app` can change even when the node's own prelim is
stable, because the NA downpass aggregates children differently (using
intersection/union/strip cases) from the raw OR of children's states.

If the downpass stops at node N because N's NA-aware prelim didn't
change, but N's child M *did* change prelim, then `children_app` at N
is different from before. The dirty-flag scheme won't revisit N, so
N's `final_` for NA blocks may be stale.

### Impact

- `fitch_na_pass3_score()` uses `final_` for `ss_app` (applicability).
  A stale `ss_app` can make `divided_length` slightly wrong.
- Indirect length calculations use `final_` for virtual-root
  computation, so candidate scores can be slightly wrong.
- **Conservative**: `full_rescore()` always runs before accepting a
  move, so final results are never affected.
- Same design class as the documented `extract_divided_steps` heuristic
  (ts_tbr.cpp:39-41) which uses stale `local_cost` for NA blocks.

### If this ever needs fixing

Mark the entire rootward path from `clip_ancestor` as dirty:

```cpp
int node = clip_ancestor;
while (node != root) {
    dirty[node] = true;
    node = tree.parent[node];
}
```

This is O(depth) extra work per clip, acceptable for correctness.
Currently not worth doing because full_rescore is authoritative.

---

## upweight_mask coverage

During ratchet perturbation, `upweight_mask` doubles the contribution
of selected characters. Every function that computes EW step counts
must account for it. The pattern:

```cpp
int ns = popcount64(needs_step);
if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
extra_steps += blk.weight * ns;
```

**Sites that must have this** (all verified correct as of 2026-03-19):

| Function | File | Status |
|----------|------|--------|
| `fitch_downpass` | ts_fitch.cpp | ✓ |
| `fitch_incremental_downpass` | ts_fitch.cpp | ✓ |
| `fitch_indirect_length` | ts_fitch.cpp | ✓ |
| `fitch_indirect_length_bounded` | ts_fitch.cpp | ✓ (fixed T-096) |
| `fitch_indirect_length_cached` | ts_fitch.cpp | ✓ (fixed T-096) |
| `fitch_na_indirect_length` | ts_fitch_na_incr.h | ✓ |
| `fitch_na_indirect_length_bounded` | ts_fitch_na_incr.h | ✓ |
| `fitch_na_indirect_length_cached` | ts_fitch_na_incr.h | ✓ |
| `fitch_na_score` Pass 1 (standard blocks) | ts_fitch_na.h | ✓ |
| `fitch_na_score` Pass 3 | ts_fitch_na.h | ✓ |
| `fitch_na_pass3_score` | ts_fitch_na_incr.h | ✓ |
| `fitch_na_incremental_downpass` (standard blocks) | ts_fitch_na_incr.h | ✓ |
| `nx_cost` in TBR | ts_tbr.cpp | ✓ (fixed T-096) |
| `nx_cost` in SPR | ts_search.cpp | ✓ (fixed T-096) |
| `nx_cost` in drift | ts_drift.cpp | ✓ (fixed T-096) |
| drift RFD computation | ts_drift.cpp | ✓ (fixed T-096) |

**Does NOT need upweight_mask:**
- `extract_char_steps` / `extract_divided_steps` — these extract raw
  per-pattern step counts for IW/profile scoring, which uses
  `pattern_freq` doubling instead of `upweight_mask`.
- `fitch_downpass_node` (standalone) — callers handle weighting.
- IW indirect variants — weighting baked into `iw_delta`.
