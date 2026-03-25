# Agent D Progress Log

## Current Task: T-232 (P2) — Shiny "Tips to show" input bounces back
**Status:** PARKED (GHA 23543699366)
**Started:** 2026-03-25

### Root cause
`UpdateKeepNTipsRange` reactive read `input$keepNTips` in LogMsg and
`r$oldkeepNTips` guard, creating a reactive dependency. The `observe()`
wrapper re-fired whenever the user manually changed the input, resetting
it to `nNonRogues()`.

### Fix
Wrapped `input$keepNTips` reads in `isolate()`. 3-line diff.
Committed to cpp-search: `082681da0`.

### Next
Pick up next task while GHA validates.
