# Plan: F-6 Async Search in EasyTrees Shiny App

## Problem

`StartSearch()` in `inst/Parsimony/app.R` runs `MaximizeParsimony()` synchronously
inside `withProgress()`. This blocks the entire Shiny session for the duration of
the search. The `future`/`promises` packages are already loaded but only used for
one small `MappingQuality` call.

Critically, plain `future_promise` would only unblock *other* sessions, not the
session that started the operation. `ExtendedTask` (Shiny ≥ 1.8.1) is the
official mechanism for within-session non-blocking operations.

## Approach

Use `ExtendedTask` + `future::future()` to run `MaximizeParsimony()` in a
background R process. The session stays fully responsive while the search runs.

Keep the existing `actionButton` UI unchanged (no switch to
`bslib::input_task_button`), and manage button state with `shinyjs::disable()`
/ `enable()`. This avoids touching the UI layout while getting the key benefit.

## Implementation

### 1. Bump Shiny version requirement

In `DESCRIPTION`: `shiny (>= 1.6.0)` → `shiny (>= 1.8.1)`.

### 2. Define the `ExtendedTask` at the top of the server function

```r
searchTask <- ExtendedTask$new(
  function(dataset, tree, concavity, strategy, maxReplicates,
           targetHits, maxSeconds, poolSuboptimal) {
    future::future({
      TreeSearch::MaximizeParsimony(
        dataset,
        tree = tree,
        concavity = concavity,
        strategy = strategy,
        maxReplicates = maxReplicates,
        targetHits = targetHits,
        maxSeconds = maxSeconds,
        control = TreeSearch::SearchControl(
          poolSuboptimal = poolSuboptimal
        ),
        verbosity = 0L
      )
    }, seed = TRUE)
  }
)
```

All arguments are plain R objects (phyDat, phylo, numerics, strings) that
serialize cleanly to the future worker. Namespace-qualified calls
(`TreeSearch::`) ensure the C++ code is available in the worker process.

### 3. Refactor `StartSearch()` to extract parameters then `invoke()`

The function keeps its existing synchronous preamble (start tree computation,
parameter extraction, logging). The only change: replace the
`withProgress(MaximizeParsimony(...))` block with:

```r
disable("go")
disable("modalGo")
disable("searchConfig")
r$searchNotification <- showNotification(
  paste0("Searching (", searchMaxRep, " replicates, ", searchWtType, ")..."),
  duration = NULL, type = "message", closeButton = FALSE
)
r$searchDataHash <- r$dataHash

searchTask$invoke(
  r$dataset[SearchTips()], startTree, concavity(),
  searchStrategy, searchMaxRep, searchTargetHits,
  searchMaxSeconds, searchPoolSub
)
```

The post-invoke code (`UpdateAllTrees`, button label updates, etc.)
moves to the result observer (step 4). `StartSearch()` returns immediately.

### 4. Add a result observer

```r
observe({
  newTrees <- tryCatch(
    searchTask$result(),
    error = function(e) {
      Notification(paste("Search error:", conditionMessage(e)),
                   type = "error")
      NULL
    }
  )
  # Always clean up UI
  if (!is.null(r$searchNotification)) {
    removeNotification(r$searchNotification)
    r$searchNotification <- NULL
  }
  enable("go")
  enable("modalGo")
  enable("searchConfig")

  req(newTrees)
  if (!identical(r$dataHash, r$searchDataHash)) {
    Notification("Dataset changed during search; results discarded.",
                 type = "warning")
    return()
  }

  r$sortTrees <- TRUE
  LogComment("Overwrite any previous trees with results")
  UpdateAllTrees(newTrees)
  updateSliderInput(session, "whichTree", min = 0L,
                    max = length(r[["trees"]]), value = 0L)
  updateActionButton(session, "go", "Continue")
  updateActionButton(session, "modalGo", "Continue search")
  show("displayConfig")
  Notification("Search complete", type = "message", duration = 5)
})
```

### 5. Re-entrancy protection

`ExtendedTask` naturally prevents overlapping invocations — a second
`invoke()` while one is running queues rather than overlaps. Combined with
`disable()` on the buttons, double-search is impossible.

### 6. Modal button handling

`observeEvent(input$modalGo, ...)` still calls `removeModal()` then
`StartSearch()`. No structural change needed — `StartSearch()` just ends
with `invoke()` instead of blocking.

## Files Modified

| File | Change |
|------|--------|
| `inst/Parsimony/app.R` | Refactor `StartSearch()`, add `searchTask` + result observer |
| `DESCRIPTION` | `shiny (>= 1.6.0)` → `shiny (>= 1.8.1)` |

## Not Changed

- UI layout (`fluidPage`, `actionButton`s) — unchanged
- Logging — stays synchronous, captures what was requested
- `future`/`promises` imports — already present, still used
- `plan(multisession)` — already set up

## Risks

| Risk | Mitigation |
|------|-----------|
| Data/tree serialization to worker | phyDat/phylo are plain R lists; SearchControl is a named list |
| User changes dataset mid-search | `r$searchDataHash` snapshot; discard results if hash differs |
| Future worker lacks TreeSearch | `TreeSearch::MaximizeParsimony()` triggers autoload |
| tryCatch around reactive `result()` | Establishes dependency before throwing; tryCatch catches the error |
| Shinytest breakage | `withProgress` removal may change snapshot; update expected outputs |

## Testing

- Manual: start search, interact with other controls while it runs
- Manual: verify results appear correctly on completion
- Manual: verify error notification if dataset is invalid
- Check shinytests still pass (may need snapshot update for removed progress bar)
