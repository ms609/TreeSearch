# Brief: diagnose PR #244 "checks hang at examples"

## Mission

Find why [PR #244](https://github.com/ms609/TreeSearch/pull/244) (T-302
`feature/pol-escapa-neg-delta → cpp-search`) cannot complete a full GHA
check cycle.  The user describes the symptom as "checks not finishing
because stuck at examples."  Three prior runs were cancelled at or near
the 6 h workflow timeout.  Confirm whether the issue persists on the
**latest** run, identify the precise stalling step, and propose a fix.

Report findings + recommended fix to the calling agent in under 400
words.  Do **not** push code without confirmation — diagnostic patches
are OK; structural fixes need sign-off.

## Background — what's been ruled out

1. **R_Interactive flush hang** (`b186e801`): R_FlushConsole on captured
   stdout pipe filled the buffer in R CMD check subprocesses, causing
   indefinite blocking.  Fixed in `src/ts_parallel.cpp` on cpp-search.
   PR #244 has this fix via merge from cpp-search (`c519496e`).  Touches
   the parallel progress loop only — affects tests, not examples.

2. **DEBUG_RESCORE log flood** (T-300 reverted in `b7303ee5`): the broken
   incremental rescore printed `DEBUG_RESCORE: diff=-3` on every
   accepted move, swamping stdout.  Reverted.  PR #244 has this revert
   via merge.  (Note: cpp-search HEAD now has a NEW `DEBUG_RESCORE`
   guard from `f531bbcd` + a `DEBUG_NNI_RESCORE` from `2be8228d` —
   both should produce **zero** output if the dirty-set rescore is
   correct.  If GHA shows mismatch lines, the dirty-set fix has a hole;
   report that as a separate finding.)

## What PR #244 actually changes

Only three files vs `cpp-search`:

- `NEWS.md`
- `R/PolEscapa.R` — fixes the `LengthAdded` negative-delta bug:
  `qmApp <- qmApp[[1L]]` (was a 1-element list, used as a scalar without
  unwrap); added a `#Temp` delta clamp for the integer-overflow guard.
- `tests/testthat/test-PolEscapa.R` — adds a token-6 regression test
  for the qmApp fix.

No C++ changes.  The hang, if it persists, is either in the R-level
`LengthAdded` example/test or in something inherited from cpp-search
that didn't get exercised on cpp-search's own CI.

## Recent runs (as of 2026-05-19 06:40 UTC)

```
26076620029  R-CMD-check                 in_progress (~2 h, started 04:40)
  ├─ ubuntu-24.04 (release)              pass     9m
  ├─ ubuntu-24.04 (4.1)                  pass    10m
  ├─ ubuntu-24.04 (devel)                pending  (>2 h)
  ├─ ubuntu-24.04-arm (release)          pending  (>2 h)
  ├─ macos-15-intel (release)            pass    19m
  ├─ macOS-latest (release)              pending  (>2 h)
  └─ windows-latest (release)            pass    25m
26076620046  R-CMD-check-ASAN            in_progress
  ├─ AddressSanitizer examples           pass    43m
  ├─ AddressSanitizer vignettes          pass    44m
  └─ AddressSanitizer tests              pending

26062924642, 26054234185, 26053699158    all CANCELLED (2026-05-18)
```

So **ubuntu-24.04 release + 4.1 pass in ~10 min**, but devel + arm +
macOS-latest hang past 2 h.  AddressSanitizer **examples passed** (43 m)
— so the literal "examples" stage isn't the problem on every runner.
The user's "stuck at examples" framing may be either (a) misremembered
from a prior cycle where it really did hang there, or (b) referring to
a specific platform that doesn't complete.  Verify with the actual logs.

## Investigation plan

1. **Confirm the current symptom.** For each pending job in run
   `26076620029` and `26076620046`, fetch the live log tail and
   identify the most recent output line:
   ```bash
   gh run view 76668864508 --log 2>&1 | tail -80    # ubuntu-24.04 (devel)
   gh run view 76668864476 --log 2>&1 | tail -80    # ubuntu-24.04-arm
   gh run view 76668864480 --log 2>&1 | tail -80    # macOS-latest
   gh run view 76668864520 --log 2>&1 | tail -80    # ASAN tests
   ```
   (Use `gh run view --log-failed` if the job has been killed; for
   in-progress use `gh api repos/ms609/TreeSearch/actions/jobs/<id>/logs`.)
   The last R CMD check stage emitted before the stall is the suspect.

2. **Re-pull cancelled-run last lines.** Same procedure on runs
   `26062924642` (the most recent cancellation pre-cpp-search-merge).
   If the cancelled run stalled at a different step than the current
   one, the cpp-search merge changed the failure mode — useful signal.

3. **Inspect candidates uncovered by step 1.**
   - If the hang is in `R CMD check` "Running examples..." → look at
     `man/LengthAdded.Rd` and any other Rd touched by recent commits.
     Vinther2008 `LengthAdded(trees, char)` runs on 9 trees × n_tip
     leaves.  If T-302's `qmApp <- qmApp[[1L]]` fix accidentally
     changes the iteration shape — e.g. now passes a vector where a
     scalar was expected, triggering recycling and a many-iteration
     inner loop — that would manifest as a slow example.
   - If the hang is in `tests/`  → focus on `test-PolEscapa.R` (newly
     added) and any test that calls `LengthAdded`.
   - If the hang is in `vignettes/` → already passed on ASAN, but
     check if any non-ASAN job stalls there.

4. **Local repro.** Once the suspect stage is identified, build PR #244
   locally (`R CMD INSTALL --library=/tmp/pr244 .` after `git
   checkout feature/pol-escapa-neg-delta`) and run the offending file
   directly:
   ```bash
   R --no-save -e 'library(TreeSearch, lib.loc="/tmp/pr244"); example(LengthAdded)'
   ```
   Time it.  Compare to the same invocation on `cpp-search` HEAD
   (without the T-302 changes).  If the slowdown is real, bisect the
   three files in PR #244 to localise.

5. **Special attention to the `#Temp` delta clamp.**  Mentioned in
   commit `bebae3a69`.  Find it (`git show bebae3a6 -- R/PolEscapa.R`)
   and check whether the clamp can fail to terminate a loop on certain
   inputs (e.g. NA propagating, or comparing `NA_integer_ > 0`
   evaluating to `NA` and being treated as `FALSE`).

## Tools available

- `gh` CLI for runs + jobs + logs.
- R 4.7-devel locally (Windows); GHA covers arm/devel/macOS that you
  can't repro locally.
- The package builds in ~5 min with the listed `R CMD INSTALL` flags.
- Do **not** start more GHA runs; one in-flight is enough until the
  diagnosis is clear.

## Constraints

- Token-limited regime: prefer reading one tail of one log over
  downloading entire job archives.
- Don't touch `cpp-search` or any other branch's source.  If a fix is
  needed in cpp-search (because PR #244 inherits something broken),
  surface it — don't apply unilaterally.
- If the hang turns out to be a platform-specific R-devel regression
  unrelated to PR #244, file as such and stop.

## Deliverable

- One section: "Where it hangs" (file + step + last-emitted line).
- One section: "Why" (root cause hypothesis).
- One section: "Proposed fix" (one paragraph) OR "Out of scope — file
  upstream" (if R-devel or similar).
- No code commits without confirmation.
