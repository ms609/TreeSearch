# T-260: VTune TBR Per-Evaluation Overhead Analysis

**Date:** 2026-03-26
**Agent:** E
**CPU:** Intel Core i7-10700 @ 2.90 GHz (Comet Lake, 10th gen)
**Sampling:** User-mode software sampling (VTune 2025.10)
**Dataset:** Dikow2009 (88 tips, EW parsimony)
**Workload:** 50 random starts × (Wagner → NNI → 20 TBR passes) = 1000 TBR passes
**Total CPU time:** 30.96s (of which TreeSearch.dll = 23.71s = 76.6%)

## Module breakdown

| Module | CPU Time | % |
|--------|:--------:|:-:|
| TreeSearch.dll | 23.71s | 76.6% |
| ucrtbase.dll | 6.00s | 19.4% |
| R.dll | 1.10s | 3.6% |
| Other | 0.15s | 0.5% |

## Top hotspots (TreeSearch.dll + attributed ucrtbase)

### By logical category

| Category | Time | % of total | Key functions |
|----------|:----:|:----------:|---------------|
| **Full NA-aware scoring** | 9.03s | 29.2% | `fitch_na_score` (includes NNI path: 3.62s) |
| **StateSnapshot save/restore** | 4.53s | 14.6% | `save` 2.15s, `restore` 1.97s, `restore_prealloc_undo` 0.18s (memcpy in ucrtbase) |
| **Incremental scoring** | 2.28s | 7.4% | `fitch_na_indirect_length_cached` 1.02s, `fitch_na_pass3_score` 0.89s, `fitch_na_indirect_length_bounded` 0.37s |
| **Tip state reloading** | 1.62s | 5.2% | `load_tip_states` (called from `reset_states` → `full_rescore`) |
| **SIMD bit ops** | ~2.0s | 6.5% | `any_hit_reduce` 1.60s, `or_reduce` 0.21s, `any_hit_reduce3` 0.31s |
| **Buffer zeroing** | ~1.20s | 3.9% | `std::fill` in `reset_states()` — zeroes prelim, final_, down2, subtree_actives, local_cost |
| **TBR orchestration** | ~1.9s | 6.1% | `tbr_search` 1.06s, `precompute_vroot_cache` 0.46s, `fitch_join_states` 0.13s, `collect_main_edges` 0.11s, `validate_topology` 0.07s, `fast_hash` 0.06s |
| **Data setup** | ~0.9s | 2.9% | `count_state_occurrences` 0.64s, `simplify_patterns` 0.12s, `build_dataset` 0.12s |
| **Memory management** | ~0.8s | 2.6% | `malloc_base` 0.77s |
| **popcount** | ~0.43s | 1.4% | `popcount64` (multiple sites) |
| **Hash set destructor** | 0.14s | 0.4% | `unordered_set::~unordered_set` (TBR tabu set) |

### TBR-only breakdown (excluding NNI scoring)

Subtracting the NNI path (3.62s fitch_na_score + proportional overhead), the
TBR-specific budget is approximately:

| TBR phase | Time | % of TBR |
|-----------|:----:|:--------:|
| Full rescore scoring (`fitch_na_score`) | 5.41s | 28% |
| StateSnapshot save/restore | 4.53s | 23% |
| Incremental candidate screening | 2.28s | 12% |
| Buffer zeroing (`std::fill` in `reset_states`) | ~1.20s | 6% |
| Tip reloading (`load_tip_states`) | 1.60s | 8% |
| TBR orchestration | ~1.9s | 10% |
| SIMD / popcount / other | ~2.5s | 13% |
| **Total TBR** | **~19.4s** | **100%** |

## Key finding: `full_rescore` overhead

Every TBR candidate that passes incremental screening triggers:

1. `state_snap.save()` — memcpy ~190 KB (5 arrays × n_node × total_words)
2. `apply_tbr_move()` — modifies topology + states
3. `full_rescore()` = `reset_states()` + `score_tree()`
   - `reset_states()`: 5× `std::fill(0)` + `load_tip_states()`
   - `score_tree()`: `fitch_na_score()` (full 3-pass)
4. If rejected: `state_snap.restore()` — memcpy ~190 KB back

**The non-scoring overhead of a single candidate evaluation
(save + zero + load_tips + restore) totals 7.35s = 37.8% of TBR time.**

The snapshot mechanism itself (save+restore = 4.53s) is an optimization
over the alternative (re-running `full_rescore` after rejection). But the
`reset_states()` step — zeroing all arrays before the downpass overwrites
them — is likely unnecessary since the Fitch downpass will recompute all
internal node values from tips up.

## Top 3 actionable hotspots

### 1. StateSnapshot save/restore — 14.6% (4.53s)

**What:** Full-array memcpy of prelim, final_, down2, subtree_actives,
local_cost, and postorder before each candidate evaluation. Restore copies
everything back when the move is rejected.

**Why it's expensive:** At 88 tips: n_node=175, total_words≈30 → each
state array is ~42 KB. With 5 arrays + cost array + postorder, each
save/restore copies ~190 KB. At 180 tips, this doubles.

**Potential fixes:**
- **Selective save/restore**: Only save nodes affected by the TBR move
  (the clip subtree path + regraft path to root). Requires tracking dirty
  nodes in `apply_tbr_move()`.
- **Copy-on-write / versioned arrays**: Use generation counters instead
  of bulk copy.
- **Eliminate the need**: If `full_rescore()` is made cheaper (see #2),
  the restore path could simply re-run scoring instead of restoring from
  snapshot.

### 2. `reset_states()` (zero + reload tips) — 9.1% (2.82s)

**What:** `full_rescore()` calls `reset_states()` which zeroes all 5 state
arrays then copies tip data back from the dataset. This runs before every
`score_tree()`.

**Why it may be unnecessary:** The Fitch downpass computes every internal
node's `prelim` from its children's values (bottom-up), overwriting whatever
was there. The uppass similarly overwrites `final_`. The zeroing is only
needed if the scoring algorithm reads uninitialized memory — but if the
postorder traversal visits every internal node, it never does.

**Potential fix:** Replace `reset_states()` with just `load_tip_states()`.
Verify that the NA-aware passes (down2, subtree_actives) also fully
overwrite internal nodes during their traversals. If they do, save 3.9%
immediately (the std::fill cost) and reduce tip loading to only the
arrays that aren't fully recomputed.

### 3. `fitch_na_score` as authoritative rescore — 29.2% (9.03s)

**What:** The full 3-pass NA-aware Fitch algorithm is called for every
candidate that passes incremental screening. This is the authoritative
score used to accept/reject moves.

**Why it dominates:** It's the core algorithm — this is expected. But
it's called much more often than strictly necessary because incremental
scoring is only a screening heuristic.

**Potential fixes:**
- **Improve incremental accuracy**: If incremental scoring matched
  full-rescore more closely, fewer candidates would need full evaluation.
  Currently ~every clip with a viable candidate triggers full_rescore.
- **Deferred full rescore**: Accept based on incremental score, batch
  full rescores periodically (risk: score drift).
- **This is also addressed indirectly by fixes #1 and #2**: reducing
  the per-evaluation overhead means each full_rescore call is cheaper.

## Estimated impact of fixes

| Fix | Savings | Effort |
|-----|:-------:|:------:|
| Eliminate `std::fill` in `reset_states` | ~3.9% (~1.2s) | Low — verify NA invariants, remove 5 fill calls |
| Selective StateSnapshot (save/restore only dirty nodes) | ~10–12% (~3–4s) | Medium — track dirty set in apply_tbr_move |
| Reduce `load_tip_states` scope (only reload modified arrays) | ~2–3% (~0.6–0.9s) | Low — check which tip arrays are read by scoring |
| **Combined** | **~16–19%** | — |

## Raw VTune data

Results stored in `vtune-tbr-out/` (gitignored). Regenerate with:
```bash
"C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/vtune.exe" \
  -collect hotspots -result-dir vtune-tbr-out \
  -- Rscript dev/vtune-tbr-driver.R
```
