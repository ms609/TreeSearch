# Review: cpp-search (uncommitted) @ per-sector column-axis reduction (ts_sector.cpp)

## Verdict
**Block** — F1 is a guaranteed heap overflow / wrong score the instant the
reduction fires (every sector with >=1 constant column). Lens-specific
eligibility + edge-case checks are otherwise clean, conditional on F1's fix.

## Spec recap
Drop every character CONSTANT-within-{sector tips + HTU} (Fitch intersection
non-empty => 0 steps on every topology) and re-pack survivors into fewer blocks
grouped by n_states, gated to plain EW (weight 1, no upweight, no NA). Claim:
inner-sector trajectory byte-identical, only faster.

## Findings

### F1. subtree.total_words / n_blocks never re-synced after reduction — BLOCKER
**Where:** ts_sector.cpp:519-520 (set to ds.* before reduction) vs the helper's
shrink at ts_sector.cpp:401-402; consumed at ts_tree.cpp:54 (load_tip_states
memcpy n_tip*total_words), ts_fitch.cpp:88/92/103 (fitch_downpass),
ts_fitch.cpp:517 (compute_insertion_edge_sets `tw=tree.total_words`),
ts_sector.cpp:845/909 (build_ras_sector `tw=t.total_words`).
**Claim:** rd.subtree.total_words/n_blocks stay at the ORIGINAL ds values; the
state arrays are allocated with the NEW (smaller) rd.data.total_words (651,657)
but every scorer + load_tip_states stride by the stale large total_words.
**Why it matters:** load_tip_states memcpys n_sector_tips*old_tw*8 bytes from
rd.data.tip_states (now sized n_sector_tips*ntw) into prelim (sized n_node*ntw)
=> heap overflow on BOTH source and destination. Then every per-node Fitch
offset (node*old_tw) overruns the new-sized prelim/edge_set. Crash or garbage
score on essentially every sector once the flag is on.
**Evidence:** dev/red-team/reviews/cpp-search-sect-colreduce/repro-01.R (NOT run
— Windows). The desync is established by static read: only rd.data.* is updated
(grep of `total_words =` / `n_blocks =` shows 401-402 touch d=rd.data;
519-520 set subtree to ds and are never revisited). The shrink path is the
only path that reaches the scorers with mismatched strides; the no-op early
returns (333-334) leave rd.data unchanged so subtree stays consistent there.
**Suggested fix:** after the `if (kSectColReduce) reduce_sector_columns_ew(...)`
call, set rd.subtree.total_words = rd.data.total_words and rd.subtree.n_blocks =
rd.data.n_blocks (or do it inside the helper alongside d.total_words/d.n_blocks).

### F2. flag read once at static init — LOW
**Where:** ts_sector.cpp (kSectColReduce static-const lambda).
**Claim:** TS_SECT_COLREDUCE is captured at .dll load; setting it later in-process
(e.g. between MaximizeParsimony calls in one R session) is ignored.
**Why it matters:** test ergonomics only; no data race (init precedes threads).
**Suggested fix:** document "set before first sector build", or read per-call.

## Coverage notes
No new test file is included in the patch. The parallel A/B (full
MaximizeParsimony off-vs-on, dScore==0, 3 datasets x 3 seeds) is the only check.
It WILL hit F1 on run 1 on any matrix with per-sector constant columns (i.e. all
real morphological data) — it crashes (ASan) or returns dScore!=0. So F1 has low
marginal value to the team; the distinctive lens result is what holds AFTER F1.

## Eligibility + edge-case lens — per-item verdict (conditional on F1 fixed)
- upweight_mask mid-ratchet misread: NO. Gate reads the COPIED per-block masks
  (rd.data.blocks = ds.blocks at 561 carries live active_mask/upweight_mask), so
  it faithfully reflects whatever phase the live ds is in. Robust. FINE.
- all-constant sector (surv.empty): early-return, no reduction. FINE.
- surv.size() >= tot_active no-op: early-return. FINE.
- mixed n_states / n_chars==64 / single survivor: grouping + repack correct.
  ntw += nst per block is right for the transposed layout (a <=64-char block =
  exactly n_states words); (nchar==64)?~0ULL:(1<<nchar)-1 avoids shift-UB; si
  walks surv in the same sorted order blocks are emitted => alignment holds. FINE.
- partial application across sectors: each build_reduced_dataset rebuilds rd
  fresh (no cross-sector state); dropped columns contribute exactly 0 so ties
  break identically => trajectory-identical per sector. FINE (given F1).
- per-pattern array remap (pattern_index[k]=k, min_steps/precomputed_steps stale,
  ew_offset untouched): VERIFIED HARMLESS for EW. The EW indirect path
  (fitch_indirect_length_cached, fitch_indirect_cached_flat_x4,
  compute_insertion_edge_sets, fitch_downpass) reads only n_blocks,
  block_word_offset, active_mask, n_states, weight, upweight_mask — never
  pattern_index/min_steps/precomputed_steps (those are IW/profile only,
  ts_fitch.cpp:879/938/1077). ew_offset is a uniform constant added post-hoc,
  trajectory-invariant. So NO second wrong-score finding here.
- 0-step invariant for polymorphic tips + HTU-as-set: correct. constant bit c set
  iff EXISTS state s with (AND over all sector tips incl. HTU of word s) bit c set
  = "intersection of all tip state-sets non-empty" = Fitch 0 steps on every tree.
  HTU loop covered (t in [0,n_sector_tips), HTU at index n_real_tips). FINE.
- has_inapplicable / NA: gate bails per-block (305). Complete. FINE.
- IW/profile sectors: scoring_mode copied from ds; gate bails. FINE.
