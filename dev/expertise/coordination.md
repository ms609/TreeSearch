# Coordination Expertise — TreeSearch

## Purpose

Review the overall state of multi-agent work. Update `coordination.md`,
propose new tasks, resolve blockers. This is the "project manager" role.

## Workflow

1. **Read all agent files** (`agent-a.md` through `agent-f.md`):
   - Who is working on what?
   - Is anyone stuck or blocked?
   - Has anyone finished a task without updating to-do.md?

2. **Read `to-do.md`**:
   - Are completed tasks moved to the Completed section?
   - Are task statuses accurate?
   - Are priorities still correct given current project state?
   - Are there enough OPEN tasks to keep all agents busy?
   - Adjust standing task priorities per the dynamic priority rule.

3. **Read `coordination.md`**:
   - Update the Agent Status table from agent files.
   - Update Known Issues if any have been resolved.
   - Add new Architecture Decisions if agents have made significant choices.

4. **Read `AGENTS.md`** (bottom sections):
   - Check for newly documented completed work.
   - Verify that documentation matches what agents report.

5. **Propose new tasks** if needed:
   - If <6 OPEN specific tasks, look at `coordination.md` strategic
     objectives and break the next one into concrete, assignable tasks.
   - If agents have reported findings (from red-team or profiling),
     ensure those are captured in to-do.md.

6. **Update all files**:
   - `coordination.md` — agent status, any new issues or decisions
   - `to-do.md` — new tasks, priority adjustments, status corrections
   - `agent-X.md` — mark your own task as complete

## Task Creation Guidelines

Good tasks are:
- **Specific**: "Profile ratchet inner loop for Zhu2013 dataset" not
  "Investigate performance"
- **Scoped**: Completable by one agent in one session (~1-2 hours)
- **Independent**: Minimal overlap with other tasks (check Blocks column)
- **Testable**: Clear success criteria (tests pass, benchmark improves, etc.)

When deriving tasks from strategic objectives:
- Break Phase 6 steps into individual tasks (T-001 through T-005 already done)
- For code quality work, group related TODOs into one task per file/module
- For documentation, one task per major section (vignettes, function docs, etc.)

## Priority Guidelines

| Priority | Criteria |
|----------|----------|
| P0 | Blocks multiple agents or causes incorrect results |
| P1 | Blocks the next strategic objective or is a correctness bug |
| P2 | Important but not blocking; performance improvements |
| P3 | Nice to have; cleanup; future-looking |

## Cross-Agent Conflict Detection

Watch for:
- Two agents modifying the same file (especially `ts_rcpp.cpp`,
  `TreeSearch-init.c`, `R/RcppExports.R`)
- Incompatible parameter changes to the same Rcpp bridge function
- One agent's optimization breaking another's assumptions

If conflicts are detected, flag them in `to-do.md` as P0 and note
which agents are affected.
