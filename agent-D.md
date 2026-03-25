# Agent D Progress Log

## Current Task: S-COORD (round 16)
**Status:** IN PROGRESS
**Started:** 2026-03-25

### Progress
- Claimed S-COORD. No unclaimed issues in issues.md.
- Verified T-213 (impose_constraint) already on cpp-search, 88 tests pass. Closed out.
- Discovered **P0 scoring regression** on cpp-search:
  - Commit `a48bfc4ad` ("Simplify 'all-[0,?]' characters in inapplicable datasets")
    incorrectly allowed simplification of characters with `?` tokens in inapplicable
    datasets. The `?` token includes the inapplicable bit → three-pass NA scoring is
    topology-dependent → simplification loses steps.
  - Impact: Vinther2008 EW pectinate: 161 vs correct 139. 59 test failures across
    `test-ts-iw.R`, `test-PolEscapa.R`, `test-SearchControl.R`, `test-ts-collapsed.R`.
  - Filed as T-218. Fix: revert to conservative bypass. PR #224, GHA 23534918452.
- T-212 GHA run 23528636505 FAILED — 59 failures are entirely due to this scoring
  regression. Should pass once PR #224 is merged.

### Coordination findings
- Agent B: S-RED focus 10 on IW/profile scoring (in progress)
- Agent C: T-214 multi-split constraint bug (in progress)
- Agent E: IDLE
- Agent G: IDLE (T-182 parked)
- T-182 PR #221: GHA PASS, ready for human merge
- T-207 PR #222: GHA PASS, ready for human merge
- T-213: DONE, cleaned up (closed out from to-do.md → completed-tasks.md)
- Pre-existing issue: test-PolEscapa.R:16 error message mismatch ("levels for 5" vs "6")
EOF 2>&1
