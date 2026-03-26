# Agent A Progress Log

## Current Task
**Status:** IDLE

## Session: 2026-03-26 — S-RED focus 9 review

### Completed: S-RED standing task — focus area 9 (Wagner & addition trees)

Reviewed `ts_wagner.h/.cpp` (595 lines) and `ts_constraint.h/.cpp` (736+144 lines).

**Key findings:**
- No bugs found in Wagner tree construction (incremental scoring, constraint mapping, 3-taxon base case, biased addition, random constrained tree)
- Latent stale-reference issue in `impose_one_pass()` (best_node relocated when move_out_root is a direct child) — negligible severity, mitigated by retry loops and TBR enforcement
- `regraft_violates_constraint()` DFS timestamp logic verified correct
- `classify_clip_constraints()` bit masking and FORBIDDEN classification correct
- 902 constraint-related tests pass; 80/80 adversarial tests pass

No new bugs filed.

### Earlier: T-242 investigation (closed)
Confirmed T-242 was a display bug (ThreadSafePool::extract_into() resetting hits_to_best), not a search quality regression. Actual IW hit rate ~60-67%.

## Session: 2026-03-25 (evening) — Summary

### Completed: T-208 + T-211 → PR #229

Implemented `random_constrained_tree()` and fixed three `impose_constraint()`
bugs on `feature/random-constrained-tree` (worktree `TS-RCT`).

**GHA run 23557186264:** 0 FAIL, 10927 PASS on both Ubuntu and Windows.

**PR #229** created to cpp-search.

#### Commits (5)
1. `27e81942` feat: random_constrained_tree()
2. `8650522a` docs: update T-212 test comments
3. `728ec297` fix: impose_constraint() bail-out threshold and return value
4. `939ea5bf` fix: impose_constraint handles root-child moves via topology_spr()
5. `d523c99d` chore: add PCSA, reconverged, reconverges to WORDLIST

#### impose_constraint() bugs fixed
| Bug | Fix |
|-----|-----|
| Bail-out threshold n_tip/4 too aggressive | Raised to n_tip |
| Return 0 for both "no violations" and "bailed out" | Return -1 on bail-out |
| spr_clip() can't detach root children | New topology_spr() helper |

### Earlier (same session)
- Investigated random_constrained_tree history: never committed to any branch
- Commit `61fbd03d` on cpp-search: corrected T-212 test comments
