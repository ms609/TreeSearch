# Module: Search
#
# Owns: searchTask (ExtendedTask), StartSearch(), result observer, search
# config modal, scoring, and weighting logic.
#
# Owns inputs: go, modalGo, searchConfig, strategy, maxReplicates,
#   targetHits, timeout, epsilon, searchWithout, implied.weights, concavity,
#   nThreads.
#
# Reactive args:
#   r              AppState reactiveValues
#   AnyTrees       reactive (from trees.R)
#   HaveData       reactive (from trees.R)
#   UpdateAllTrees function (from trees.R)
#   log_fns        named list: LogMsg, LogCode, LogComment
#
# Returns a list of reactives/functions consumed by other server files:
#   scores, concavity, DisplayTreeScores

# ---------------------------------------------------------------------------
# UI — returns a named list so scattered elements can be placed individually
# in ui.R (same pattern as downloads_ui).
# ---------------------------------------------------------------------------
search_ui <- function(id) {
  ns <- NS(id)
  list(
    label   = tags$label("Search", class = "control-label",
                          style = "display: block; margin-top: -15px;"),
    config  = actionButton(ns("searchConfig"), "Configure",
                           icon = Icon("gears")),
    go      = hidden(actionButton(ns("go"), "Search",
                                  icon = Icon("magnifying-glass"))),
    results = htmlOutput(ns("results"))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
search_server <- function(id, r, AnyTrees, HaveData, UpdateAllTrees, log_fns) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Unpack logging functions
    LogMsg     <- log_fns$LogMsg
    LogCode    <- log_fns$LogCode
    LogComment <- log_fns$LogComment

    ##########################################################################
    # Local helpers
    ##########################################################################

    DatasetTips <- reactive(names(r$dataset))
    SearchTips  <- reactive(setdiff(DatasetTips(), r$searchWithout))

    ##########################################################################
    # Weighting / concavity
    ##########################################################################

    weighting <- reactive(
      if (length(input$implied.weights) > 0) {
        input$implied.weights
      } else {
        "on"
      }
    )

    wtType <- reactive(switch(weighting(),
                              "on" = paste0("k = ", signif(concavity(), 3)),
                              "off" = "EW",
                              "prof" = "PP"))

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

    # Show/hide concavity slider when weighting mode changes
    observeEvent(input$implied.weights, {
      switch(input$implied.weights,
             "on" = show("concavity"),
             hide("concavity")
      )
      DisplayTreeScores()
    })

    observeEvent(input$concavity, {
      DisplayTreeScores()
    }, ignoreInit = TRUE)

    ##########################################################################
    # Scores (cached)
    ##########################################################################

    scores <- bindCache(reactive({
      if (!HaveData() || !AnyTrees()) {
        return(NULL)
      }
      PutTree(r$trees)
      PutData(r$dataset)
      LogMsg("scores(): Recalculating scores with k = ", concavity())
      withProgress(tryCatch(
        signif(TreeLength(
          RootTree(r$trees, 1),
          r$dataset,
          concavity = concavity()
        )),
        error = function (x) {
          if (HaveData() && AnyTrees()) {
            cli::cli_alert(x[[2]])
            cli::cli_alert_danger(x[[1]])
            Notification(type = "error",
                         "Could not score all trees with dataset")
          }
          NULL
       }),
       value = 0.85, message = "Scoring trees")
    }), r$treeHash, r$dataHash, concavity())

    ##########################################################################
    # DisplayTreeScores
    ##########################################################################

    DisplayTreeScores <- function () {
      LogMsg("DisplayTreeScores()")
      treeScores <- scores()
      score <- if (is.null(treeScores)) {
        "; could not be scored from dataset"
      } else if (length(unique(treeScores)) == 1) {
        paste0(", each with score ", treeScores[1], " (", wtType(), ")")
      } else {
        paste0(" with scores ", min(treeScores), " to ", max(treeScores),
               " (", wtType(), ")")
      }

      msg <- paste0(
        length(r$allTrees), " trees in memory: ",
        length(r$trees), " sampled",
        score
      )
      confText <- SearchConfidenceText(r$searchTotalHits, r$searchTotalReps)
      html <- if (!is.null(confText)) {
        paste0(msg, "<br><small style='color:#666'>", confText, "</small>")
      } else {
        msg
      }
      output$results <- renderUI(HTML(html))
      invisible(msg)
    }

    ##########################################################################
    # ExtendedTask for async search
    ##########################################################################

    searchTask <- ExtendedTask$new(
      function(dataset, tree, concavity, strategy, maxReplicates,
               targetHits, maxSeconds, poolSuboptimal, nThreads) {
        future::future({
          args <- list(
            dataset,
            tree = tree,
            concavity = concavity,
            strategy = strategy,
            maxReplicates = maxReplicates,
            targetHits = targetHits,
            maxSeconds = maxSeconds,
            nThreads = nThreads,
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

    ##########################################################################
    # StartSearch
    ##########################################################################

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
        searchNThreads <- if (length(input$nThreads)) as.integer(input$nThreads) else 1L
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
          if (searchNThreads > 1L)
            paste0("  nThreads = ", searchNThreads, "L,"),
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
                 if (searchNThreads > 1L) paste0(", ", searchNThreads, " threads") else "",
                 ")\u2026"),
          duration = NULL, type = "message", closeButton = FALSE
        )
        r$searchDataHash <- r$dataHash
        output$results <- renderUI(HTML(paste0(
          "Searching (", searchMaxRep, " replicates, ", wtType(),
          if (searchNThreads > 1L) paste0(", ", searchNThreads, " threads") else "",
          ")\u2026"
        )))

        searchTask$invoke(
          searchDataset, startTree, searchConcavity,
          searchStrategy, searchMaxRep, searchTargetHits,
          searchMaxSeconds, searchPoolSub, searchNThreads
        )
      }
    }

    ##########################################################################
    # Input observers
    ##########################################################################

    observeEvent(input$searchWithout, {
      r$searchWithout <- input$searchWithout
    }, ignoreInit = TRUE)

    observeEvent(input$go, StartSearch(), ignoreInit = TRUE)
    observeEvent(input$modalGo, {
      removeModal()
      StartSearch()
    }, ignoreInit = TRUE)

    ##########################################################################
    # Search config modal
    ##########################################################################

    observeEvent(input$searchConfig, {
      nCores <- max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)
      updateSelectInput(session, "implied.weights",
                        selected = input$implied.weights)
      updateSliderInput(session, "concavity", value = input$concavity)
      updateNumericInput(session, "epsilon", value = input$epsilon)
      updateSelectInput(session, "strategy", selected = input$strategy)
      updateSliderInput(session, "maxReplicates", value = input$maxReplicates)
      updateSliderInput(session, "targetHits", value = input$targetHits)
      updateSliderInput(session, "timeout", value = input$timeout)
      if (nCores > 1L) {
        updateSliderInput(session, "nThreads", value = input$nThreads)
      }
      showModal(modalDialog(
        easyClose = TRUE,
        fluidPage(column(6,
          # --- Step weighting ---
          tags$h6(tags$strong("Step weighting")),
          selectInput(ns("implied.weights"), "Mode",
                     list("Implied" = "on", "Profile" = "prof",
                          "Equal" = "off"), "on"),
          sliderInput(ns("concavity"), "Concavity constant", min = 0L,
                     max = 3L, pre = "10^", value = 1L),
          # --- Parallelization (only shown when > 1 core available) ---
          if (nCores > 1L) {
            tagList(
              tags$h6(tags$strong("Parallelization"),
                      style = "margin-top: 10px;"),
              sliderInput(ns("nThreads"), "Parallel search threads",
                          min = 1L, max = nCores,
                          value = if (length(input$nThreads)) input$nThreads
                                  else max(1L, floor(nCores / 2L)),
                          step = 1L)
            )
          }
        ), column(6,
          # --- Search intensity ---
          tags$h6(tags$strong("Search intensity")),
          selectInput(ns("strategy"), "Search strategy",
                     list("Auto" = "auto", "Sprint" = "sprint",
                          "Default" = "default", "Thorough" = "thorough"),
                     "auto"),
          sliderInput(ns("targetHits"), "Stop after best score found N times",
                      min = 1L, max = 50L, value = 10L, step = 1L),
          helpText("The search stops once the same best score has been",
                   "found independently N times.",
                   "A higher value gives greater confidence that the",
                   "result is the global optimum."),
          sliderInput(ns("timeout"), "Maximum run duration", min = 1,
                      max = 60, value = 5, post = "min", step = 1),
          sliderInput(ns("maxReplicates"), "Maximum replicates", min = 1L,
                      max = 500L, value = 100L, step = 1L),
          # --- Results to keep ---
          tags$h6(tags$strong("Results to keep"),
                  style = "margin-top: 10px;"),
          selectizeInput(ns("searchWithout"), "Exclude taxa", DatasetTips(),
                         r$searchWithout, multiple = TRUE),
          numericInput(ns("epsilon"), "Keep if suboptimal by \u2264", min = 0,
                      value = 0)
        )),
        title = "Tree search settings",
        footer = tagList(modalButton("Close", icon = Icon("rectangle-xmark")),
                         actionButton(ns("modalGo"), icon = Icon("magnifying-glass"),
                                      if(length(r$trees)) {
                                        "Continue search"
                                      } else {
                                        "Start search"
                                      }))
      ))
      show("go")
    })

    ##########################################################################
    # Async search result observer
    ##########################################################################

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
        newScore    <- attr(newTrees, "score")
        newHitsRaw  <- attr(newTrees, "hits_to_best")
        newRepsRaw  <- attr(newTrees, "replicates")
        newHits     <- if (is.null(newHitsRaw)) 0L else as.integer(newHitsRaw)
        newReps     <- if (is.null(newRepsRaw)) 0L else as.integer(newRepsRaw)
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
        updateActionButton(session, "go", "Continue")
        updateActionButton(session, "modalGo", "Continue search")
        shinyjs::show(selector = "#displayConfig")
        newCount <- length(r$allTrees)
        Notification(
          if (newCount > prevCount)
            paste0("Search complete \u2014 ", newCount, " trees in pool (+",
                   newCount - prevCount, " new)")
          else
            "Search complete",
          type = "message", duration = 5
        )
        r$searchCount <- r$searchCount + 1L
      })
    })

    ##########################################################################
    # Dataset change: reset search stats + update timeout default
    ##########################################################################

    observeEvent(r$dataset, {
      r$searchTotalHits <- 0L
      r$searchTotalReps <- 0L
      nTip <- length(r$dataset)
      nChar <- sum(attr(r$dataset, "weight", exact = TRUE))
      defaultTimeout <- max(1L, min(15L, ceiling(nTip * nChar / 20000L)))
      updateSliderInput(session, "timeout", value = defaultTimeout)
    })

    ##########################################################################
    # Button label management — react to tree/data state changes
    ##########################################################################

    observe({
      hasTrees <- !is.null(r$allTrees) && length(r$allTrees) > 0
      hasData  <- !is.null(r$dataset) && length(r$dataset) > 0
      if (!hasData) return()
      if (hasTrees) {
        treeTips <- r$allTrees[[1]]$tip.label
        dataTips <- names(r$dataset)
        if (length(intersect(dataTips, treeTips)) == length(r$dataset)) {
          updateActionButton(session, "go", "Continue")
        } else {
          updateActionButton(session, "go", "New search")
        }
      } else {
        updateActionButton(session, "go", "Search")
      }
    })

    ##########################################################################
    # Return values for other server files
    ##########################################################################

    list(
      scores            = scores,
      concavity         = concavity,
      DisplayTreeScores = DisplayTreeScores
    )
  })
}
