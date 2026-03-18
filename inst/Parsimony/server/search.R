  plan(multisession)
  
  searchTask <- ExtendedTask$new(
    function(dataset, tree, concavity, strategy, maxReplicates,
             targetHits, maxSeconds, poolSuboptimal) {
      future::future({
        args <- list(
          dataset,
          tree = tree,
          concavity = concavity,
          strategy = strategy,
          maxReplicates = maxReplicates,
          targetHits = targetHits,
          maxSeconds = maxSeconds,
          verbosity = 0L
        )
        # Only pass control when non-default, so strategy presets apply
        if (poolSuboptimal > 0) {
          args$control <- TreeSearch::SearchControl(
            poolSuboptimal = poolSuboptimal
          )
        }
        do.call(TreeSearch::MaximizeParsimony, args)
      }, seed = TRUE)
    }
  )
  
  startOpt <- options("cli.progress_show_after" = 0.1)
  
  
  LogMsg("Started server")

  concavity <- reactive({
    kExp <- if (length(input$concavity)) input$concavity else 1
    switch(weighting(),
           "on" = 10 ^ kExp,
           "off" = Inf,
           "prof" = "profile")
  })
  
  tolerance <- reactive({
    if (input$epsilon == 0) {
      sqrt(.Machine$double.eps)
    } else {
      input$epsilon
    }
  })
  
  StartSearch <- function () {
    if (!HaveData()) {
      Notification("No data loaded", type = "error")
    } else {
      startTree <- if (!AnyTrees()) {
        LogComment("Select starting tree")
        LogCode(paste0("startTree <- AdditionTree(dataset, concavity = ",
                       Enquote(concavity()), ")"))
        AdditionTree(r$dataset[SearchTips()], concavity = concavity())
      } else {
        LogComment("Select starting tree")
        treeLabels <- TipLabels(r$trees[[1]])
        if (all(SearchTips() %in% treeLabels)) {
          if (length(setdiff(treeLabels, SearchTips())) > 0) {
            if (length(r$searchWithout)) {
              LogCode(paste0(
                "searchTips <- setdiff(names(dataset), ", EnC(r$searchWithout),
                ")"),
                "startTree <- KeepTip(trees[[1]], searchTips)")
            } else {
              LogCode("startTree <- KeepTip(trees[[1]], names(dataset))")
            }
            KeepTip(r$trees[[1]], SearchTips())
          } else {
            sc <- scores()
            firstOptimal <- if (length(sc)) which.min(sc) else 1L
            LogCode(paste0("startTree <- trees[[", firstOptimal, "]]",
                           " # First tree with optimal score"))
            r$trees[[firstOptimal]]
          }
        } else {
          # Fuzzy-match labels
          matching <- TreeDist::LAPJV(adist(treeLabels, SearchTips()))$matching
          scaffold <- KeepTip(r$trees[[1]], !is.na(matching))
          scaffold[["tip.label"]] <- SearchTips()[matching[!is.na(matching)]]
          AdditionTree(r$dataset, concavity = concavity(),
                       constraint = scaffold)
        }
      }
      LogMsg("StartSearch()")
      PutData(r$dataset[SearchTips()])
      PutTree(startTree)
      LogComment("Search for optimal trees", 1)
      searchStrategy <- if (length(input$strategy)) input$strategy else "auto"
      searchMaxRep <- if (length(input$maxReplicates)) {
        as.integer(input$maxReplicates)
      } else {
        100L
      }
      searchTargetHits <- if (length(input$targetHits)) {
        as.integer(input$targetHits)
      } else {
        10L
      }
      searchMaxSeconds <- if (length(input$timeout)) {
        as.double(input$timeout) * 60
      } else {
        0
      }
      searchPoolSub <- if (length(input$epsilon) && input$epsilon > 0) {
        tolerance()
      } else {
        0
      }
      LogCode(c(
        "newTrees <- MaximizeParsimony(",
        if (length(r$searchWithout)) {
          paste0(
            "  dataset[setdiff(names(dataset), ", EnC(r$searchWithout), ")],"
          )
        } else {
          "  dataset,"
        },
        "  tree = startTree,",
        paste0("  concavity = ", Enquote(concavity()), ","),
        paste0("  strategy = \"", searchStrategy, "\","),
        paste0("  maxReplicates = ", searchMaxRep, ","),
        paste0("  targetHits = ", searchTargetHits, ","),
        if (searchMaxSeconds > 0)
          paste0("  maxSeconds = ", searchMaxSeconds, ","),
        if (searchPoolSub > 0)
          paste0("  control = SearchControl(poolSuboptimal = ", searchPoolSub, "),"),
        "  verbosity = 0",
        ")"))
      # Snapshot reactive values for the async task
      searchDataset <- r$dataset[SearchTips()]
      searchConcavity <- concavity()
      
      disable("go")
      disable("modalGo")
      disable("searchConfig")
      r$searchNotification <- showNotification(
        paste0("Searching (", searchMaxRep, " replicates, ", wtType(),
               ")\u2026"),
        duration = NULL, type = "message", closeButton = FALSE
      )
      r$searchDataHash <- r$dataHash
      output$results <- renderText(paste0(
        "Searching (", searchMaxRep, " replicates, ", wtType(), ")\u2026"
      ))
      
      searchTask$invoke(
        searchDataset, startTree, searchConcavity,
        searchStrategy, searchMaxRep, searchTargetHits,
        searchMaxSeconds, searchPoolSub
      )
    }
  }
  
  observeEvent(input$searchWithout, {
    r$searchWithout <- input$searchWithout
  }, ignoreInit = TRUE)
  
  observeEvent(input$go, StartSearch(), ignoreInit = TRUE)
  observeEvent(input$modalGo, {
    removeModal()
    StartSearch()
  }, ignoreInit = TRUE)
  
  # Handle async search completion.
  # Only searchTask$result() should be a reactive dependency;
  # isolate everything else to prevent reactive cascade re-runs.
  observe({
    newTrees <- tryCatch(
      searchTask$result(),
      validation = function(e) {
        # ExtendedTask signals validation when initial/running; not a real error
        req(FALSE)
      },
      error = function(e) {
        msg <- conditionMessage(e)
        if (nzchar(msg)) {
          Notification(paste("Search error:", msg), type = "error")
        }
        NULL
      }
    )
    isolate({
      # Clean up search-in-progress UI state
      if (!is.null(r$searchNotification)) {
        removeNotification(r$searchNotification)
        r$searchNotification <- NULL
        enable("go")
        enable("modalGo")
        enable("searchConfig")
      }

      if (is.null(newTrees)) {
        DisplayTreeScores()
        return()
      }
      if (!identical(r$dataHash, r$searchDataHash)) {
        Notification("Dataset changed during search; results discarded.",
                     type = "warning")
        DisplayTreeScores()
        return()
      }

      r$sortTrees <- TRUE

      # Accumulate trees across searches: if the new result matches the
      # current best score, merge with existing trees (dedup by topology).
      # This lets repeated "Continue search" discover additional MPTs.
      newScore    <- attr(newTrees, "score")
      newHits     <- attr(newTrees, "hits_to_best") %||% 0L
      newReps     <- attr(newTrees, "replicates")   %||% 0L
      prevCount <- length(r$allTrees)
      treesToStore <- if (
        !is.null(newScore) && !is.null(r$bestSearchScore) &&
        isTRUE(abs(newScore - r$bestSearchScore) < sqrt(.Machine$double.eps)) &&
        prevCount > 0L
      ) {
        LogComment("Same optimal score: accumulating trees across search runs")
        r$searchTotalHits <- r$searchTotalHits + newHits
        r$searchTotalReps <- r$searchTotalReps + newReps
        combined <- c(r$allTrees, newTrees)
        # Deduplicate by canonical Newick (ladderized topology string)
        nwk <- vapply(combined, function(t) {
          write.tree(ape::ladderize(t))
        }, character(1L))
        combined[!duplicated(nwk)]
      } else {
        LogComment("New or improved score: replacing trees")
        r$bestSearchScore  <- newScore
        r$searchTotalHits  <- newHits
        r$searchTotalReps  <- newReps
        newTrees
      }

      UpdateAllTrees(treesToStore)
      updateSliderInput(session, "whichTree", min = 0L,
                        max = length(r[["trees"]]), value = 0L)
      updateActionButton(session, "go", "Continue")
      updateActionButton(session, "modalGo", "Continue search")
      show("displayConfig")
      newCount <- length(r$allTrees)
      Notification(
        if (newCount > prevCount)
          paste0("Search complete — ", newCount, " trees in pool (+",
                 newCount - prevCount, " new)")
        else
          "Search complete",
        type = "message", duration = 5
      )
      r$searchCount <- r$searchCount + 1L
    })
  })
