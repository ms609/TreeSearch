# Module: Search
#
# Owns: searchTask (ExtendedTask), StartSearch(), result observer, search
# config modal, scoring, and weighting logic.
#
# Owns inputs: go, modalGo, searchConfig, strategy, maxReplicates,
#   targetHits, timeout, epsilon, searchWithout, implied.weights, concavity,
#   nThreads, inapplicable, hsjAlpha.
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
    cancel  = hidden(actionButton(ns("cancel"), "Stop",
                                  icon = Icon("circle-stop"),
                                  class = "btn-danger btn-sm",
                                  style = "margin-left: 4px;")),
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

    # Adaptive note under the targetHits slider (shown inside config modal)
    output$targetHitsNote <- renderUI({
      N <- input$targetHits
      if (is.null(N) || N < 1L) return(NULL)
      # Worst-case miss probability: lim_{R->inf} (1 - N/R)^R = exp(-N)
      helpText(title = paste0(
                 "Theoretical worst-case: exp(-", N,
                 "). After searching, the results panel uses actual hit counts."
               ),
               paste0("Probability of missing best score: ",
                      FormatMissProb(exp(-N))))
    })

    ##########################################################################
    # Weighting / concavity
    ##########################################################################

    weighting <- reactive(
      if (length(input$implied.weights) > 0) {
        input$implied.weights
      } else {
        "xpiwe"
      }
    )

    wtType <- reactive(switch(weighting(),
                              "xpiwe" = paste0("k = ", signif(concavity(), 3)),
                              "on" = paste0("k = ", signif(concavity(), 3)),
                              "off" = "EW",
                              "prof" = "PP"))

    concavity <- reactive({
      kExp <- if (length(input$concavity)) input$concavity else 1
      switch(weighting(),
             "xpiwe" = 10 ^ kExp,
             "on" = 10 ^ kExp,
             "off" = Inf,
             "prof" = "profile")
    })

    # Whether to apply extended implied weighting (missing-entries correction)
    extendedIw <- reactive(identical(weighting(), "xpiwe"))

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
             "xpiwe" = , "on" = show("concavity"),
             hide("concavity")
      )
      # Weighting mode changed: old run counts no longer apply; keep trees
      r$searchTotalHits <- 0L
      r$searchTotalReps <- 0L
      r$bestSearchScore  <- NULL
      DisplayTreeScores()
    })

    observeEvent(input$concavity, {
      # Concavity constant changed: old run counts no longer apply; keep trees
      r$searchTotalHits <- 0L
      r$searchTotalReps <- 0L
      r$bestSearchScore  <- NULL
      DisplayTreeScores()
    }, ignoreInit = TRUE)

    # Show/hide hsjAlpha input when inapplicable method changes
    observeEvent(input$inapplicable, {
      if (identical(input$inapplicable, "hsj")) {
        show("hsjAlpha")
      } else {
        hide("hsjAlpha")
      }
    }, ignoreInit = TRUE)

    # Dynamic help text for hierarchy detection (shown inside config modal)
    output$hierarchyInfo <- renderUI({
      inp <- input$inapplicable
      if (is.null(inp) || identical(inp, "brazeau")) return(NULL)
      chars <- r$chars
      if (is.null(chars) || length(chars) == 0L) {
        return(helpText(
          "No character names available for hierarchy auto-detection."
        ))
      }
      h <- tryCatch(
        withCallingHandlers(
          hierarchy_from_names(chars),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) NULL
      )
      if (is.null(h)) {
        helpText(HTML(paste0(
          "No hierarchy detected. Character names must follow the convention ",
          "<code>sup_tag</code> (primary) and ",
          "<code>sub_tag[_suffix]</code> (secondary); see ",
          "<code>?hierarchy_from_names</code>."
        )))
      } else {
        n_blocks <- length(h)
        n_chars  <- length(hierarchy_chars(h))
        helpText(paste0(
          "Detected ", n_blocks, " hierarchy block(s) covering ",
          n_chars, " character(s)."
        ))
      }
    })

    ##########################################################################
    # Async profile data preparation (with progress + cancel)
    ##########################################################################

    profileDataset      <- reactiveVal(NULL)
    profileDataHash     <- reactiveVal(NULL)
    profileNotification <- reactiveVal(NULL)
    profileProgressFile <- reactiveVal(NULL)
    profileCancelFile   <- reactiveVal(NULL)

    # Inlines PrepareDataProfile() logic so the slow StepInformation loop
    # can report per-pattern progress and check a cancel file.
    # Mirrors R/data_manipulation.R::PrepareDataProfile(); keep in sync.
    profilePrepTask <- ExtendedTask$new(
      function(dataset, progressPath, cancelPath) {
        future::future({
          if ("info.amounts" %in% names(attributes(dataset))) {
            return(dataset)
          }

          at <- attributes(dataset)
          cont <- attr(dataset, "contrast")
          nTip <- length(dataset)
          index <- at[["index"]]
          allLevels <- as.character(at[["allLevels"]])

          contSums <- rowSums(cont)
          qmLevel <- which(contSums == ncol(cont))

          if (length(qmLevel) == 0) {
            attr(dataset, "contrast") <- rbind(attr(dataset, "contrast"), 1)
            attr(dataset, "allLevels") <- c(attr(dataset, "allLevels"), "{?}")
            qmLevel <- length(allLevels) + 1L
          }

          ambigs <- which(contSums > 1L & contSums < ncol(cont))
          inappLevel <- which(colnames(cont) == "-")
          if (length(inappLevel) != 0L) {
            inappLevel <- which(apply(unname(cont), 1, identical,
                                      as.double(colnames(cont) == "-")))
            dataset[] <- lapply(dataset, function(i) {
              i[i %in% inappLevel] <- qmLevel
              i
            })
          }

          if (length(ambigs) != 0L) {
            dataset[] <- lapply(dataset, function(i) {
              i[i %in% ambigs] <- qmLevel
              i
            })
          }

          nPattern <- max(index)
          mataset <- matrix(
            unlist(dataset, recursive = FALSE, use.names = FALSE), nPattern
          )
          mataset <- t(mataset)

          maxInformative <- 0L
          cappedAny <- FALSE

          for (j in seq_len(ncol(mataset))) {
            col <- mataset[, j]
            nonAmbig <- col[col != qmLevel[1]]
            if (length(nonAmbig) == 0L) next

            tab <- table(nonAmbig)
            informative <- tab > 1L
            nInf <- sum(informative)

            singletonTokens <- as.integer(names(tab[!informative]))
            if (length(singletonTokens) > 0L) {
              mataset[mataset[, j] %in% singletonTokens, j] <- qmLevel[1]
            }

            if (nInf > 5L) {
              sortedInf <- sort(tab[informative], decreasing = TRUE)
              toRemove <- as.integer(names(sortedInf)[6:length(sortedInf)])
              mataset[mataset[, j] %in% toRemove, j] <- qmLevel[1]
              nInf <- 5L
              cappedAny <- TRUE
            }

            maxInformative <- max(maxInformative, nInf)
          }

          if (maxInformative < 2L) {
            attr(dataset, "info.amounts") <- double(0)
            return(dataset[0])
          }

          AMBIG_TOKEN <- maxInformative + 1L

          for (j in seq_len(ncol(mataset))) {
            col <- mataset[, j]
            nonAmbig <- sort(unique(col[col != qmLevel[1]]))
            newCol <- rep(AMBIG_TOKEN, length(col))
            for (i in seq_along(nonAmbig)) {
              newCol[col == nonAmbig[i]] <- i
            }
            mataset[, j] <- newCol
          }

          dupCols <- duplicated(t(mataset))
          kept <- which(!dupCols)
          copies <- lapply(kept, function(i) {
            i + which(apply(
              mataset[, -seq_len(i), drop = FALSE], 2, identical, mataset[, i]
            ))
          })
          firstOccurrence <- seq_len(dim(mataset)[2])
          for (i in seq_along(copies)) {
            firstOccurrence[copies[[i]]] <- kept[i]
          }

          cipher <- seq_len(max(kept))
          cipher[kept] <- order(kept)
          index <- cipher[firstOccurrence][index]

          mataset <- mataset[, !dupCols, drop = FALSE]
          dataset[] <- lapply(
            seq_len(length(dataset)), function(i) mataset[i, ]
          )

          # --- Slow part: StepInformation per unique pattern ---
          nPatterns <- ncol(mataset)
          info <- vector("list", nPatterns)

          for (i in seq_len(nPatterns)) {
            if (file.exists(cancelPath)) return(NULL)

            info[[i]] <- TreeSearch::StepInformation(
              mataset[, i], ambiguousTokens = AMBIG_TOKEN
            )

            writeLines(paste(i, nPatterns), progressPath)
          }

          if (file.exists(cancelPath)) return(NULL)

          maxSteps <- max(vapply(
            info, function(x) max(as.integer(names(x))), integer(1)
          ))
          info <- vapply(info, function(x) {
            ret <- setNames(double(maxSteps), seq_len(maxSteps))
            x <- x[setdiff(names(x), "0")]
            if (length(x)) {
              ret[names(x)] <- max(x) - x
            }
            ret
          }, double(maxSteps))
          if (is.null(dim(info))) {
            dim(info) <- c(1L, length(info))
          }
          attr(dataset, "index") <- index
          weight <- as.integer(table(index))
          attr(dataset, "weight") <- weight
          attr(dataset, "nr") <- length(weight)
          attr(dataset, "info.amounts") <- info
          attr(dataset, "informative") <- colSums(info) > 0

          k <- maxInformative
          lvls <- as.character(seq_len(k))
          contMatrix <- rbind(diag(k), rep(1L, k))
          dimnames(contMatrix) <- list(NULL, lvls)

          attr(dataset, "levels") <- lvls
          attr(dataset, "allLevels") <- c(lvls, "?")
          attr(dataset, "contrast") <- contMatrix
          attr(dataset, "nc") <- as.integer(k)

          if (!any(attr(dataset, "bootstrap") == "info.amounts")) {
            attr(dataset, "bootstrap") <- c(
              attr(dataset, "bootstrap"), "info.amounts"
            )
          }

          dataset
        }, seed = TRUE)
      }
    )

    # Helper: start async profile data preparation. Called from
    # StartSearch() when the user requests profile scoring and data
    # hasn't been prepared yet. NOT triggered eagerly on mode change —
    # deferred until the user actually starts a search.
    startProfilePrep <- function(dataset) {
      # Cancel any in-flight prep first.
      cf <- profileCancelFile()
      if (!is.null(cf) && !file.exists(cf)) {
        file.create(cf)
      }
      status <- tryCatch(profilePrepTask$status(), error = function(e) "initial")
      if (status == "running") return(FALSE)

      profileDataset(NULL)
      LogMsg("Starting async profile data preparation")

      progPath <- tempfile("ts_profile_prog_", fileext = ".txt")
      cancPath <- tempfile("ts_profile_cancel_", fileext = ".signal")
      profileProgressFile(progPath)
      profileCancelFile(cancPath)

      nid <- showNotification("Preparing profile scores\u2026",
                              duration = NULL, type = "message")
      profileNotification(nid)
      profilePrepTask$invoke(dataset, progPath, cancPath)
      TRUE
    }

    # Poll progress file and update notification while profile prep runs
    observe({
      progFile <- profileProgressFile()
      nid <- profileNotification()
      if (is.null(progFile) || is.null(nid)) return()
      invalidateLater(500)
      progress <- tryCatch(
        readLines(progFile, warn = FALSE), error = function(e) NULL
      )
      if (is.null(progress) || length(progress) == 0L || !nzchar(progress[1])) {
        return()
      }
      parts <- strsplit(trimws(progress[1]), "\\s+")[[1]]
      if (length(parts) != 2L) return()
      current <- suppressWarnings(as.integer(parts[1]))
      total   <- suppressWarnings(as.integer(parts[2]))
      if (is.na(current) || is.na(total) || total < 1L) return()
      pct <- round(100 * current / total)
      showNotification(
        id = nid,
        paste0("Preparing profile scores\u2026 ", current, "/", total,
               " patterns (", pct, "%)"),
        duration = NULL, type = "message"
      )
    })

    # Process profile preparation result
    observe({
      result <- tryCatch(
        profilePrepTask$result(),
        shiny.silent.error = function(e) req(FALSE),
        error = function(e) {
          LogMsg("Profile data preparation failed: ", conditionMessage(e))
          NULL
        }
      )
      isolate({
        nid <- profileNotification()
        if (!is.null(nid)) {
          removeNotification(nid)
          profileNotification(NULL)
        }
        # Clean up temp files
        pf <- profileProgressFile()
        if (!is.null(pf)) {
          suppressWarnings(file.remove(pf))
          profileProgressFile(NULL)
        }
        cf <- profileCancelFile()
        if (!is.null(cf)) {
          suppressWarnings(file.remove(cf))
          profileCancelFile(NULL)
        }
        if (!is.null(result)) {
          profileDataset(result)
          profileDataHash(r$dataHash)
          Notification("Profile scores ready \u2014 click Search to start",
                       type = "message", duration = 5)
        }
        DisplayTreeScores()
      })
    })

    # Cancel profile prep if user switches away from profile mode
    observe({
      if (!identical(concavity(), "profile")) {
        nid <- profileNotification()
        if (!is.null(nid)) {
          removeNotification(nid)
          profileNotification(NULL)
        }
        cf <- profileCancelFile()
        if (!is.null(cf) && !file.exists(cf)) {
          file.create(cf)
        }
      }
    })

    ##########################################################################
    # Scores
    ##########################################################################

    scores <- reactive({
      if (!HaveData() || !AnyTrees()) {
        return(NULL)
      }
      conc <- concavity()
      ds <- if (identical(conc, "profile")) {
        pd <- profileDataset()
        if (is.null(pd)) return(NULL)
        pd
      } else {
        r$dataset
      }
      PutTree(r$trees)
      PutData(ds)
      useXpiwe <- extendedIw()
      LogMsg("scores(): Recalculating scores with k = ", conc,
             if (useXpiwe) " (extended)")
      tryCatch(
        signif(TreeLength(
          RootTree(r$trees, 1),
          ds,
          concavity = conc,
          extended_iw = useXpiwe
        )),
        error = function (x) {
          if (HaveData() && AnyTrees()) {
            cli::cli_alert(x[[2]])
            cli::cli_alert_danger(x[[1]])
            Notification(type = "error",
                         "Could not score all trees with dataset")
          }
          NULL
       })
    })

    ##########################################################################
    # DisplayTreeScores
    ##########################################################################

    DisplayTreeScores <- function () {
      # Don't overwrite "Searching..." indicator while a search is running.
      # Guard on both fields: searchNotification can be NULL if the
      # notification was dismissed externally, but searchInProgress is the
      # authoritative flag.
      if (!is.null(r$searchNotification) || isTRUE(r$searchInProgress)) return(invisible())
      LogMsg("DisplayTreeScores()")
      treeScores <- scores()
      score <- if (is.null(treeScores) && identical(concavity(), "profile") &&
                   is.null(profileDataset()) && HaveData() && AnyTrees()) {
        "; profile scores available after search"
      } else if (is.null(treeScores)) {
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
      confText <- SearchConfidenceText(r$searchTotalHits, r$searchTotalReps,
                                        r$searchCount)
      html <- if (!is.null(confText)) {
        nS <- r$searchCount
        tooltip <- paste0(
          "Estimated as exp(-K) where K = ",
          r$searchTotalHits,
          " (runs hitting best score",
          if (!is.null(nS) && nS > 1L)
            paste0(" across ", nS, " searches")
          else
            "",
          "). Assumes independent runs. ",
          "'Maximum independent runs' limits each individual search; ",
          "this tally accumulates across all continued searches. ",
          "The config dialog shows a theoretical worst-case; ",
          "this uses actual search results."
        )
        paste0(msg, "<br><small style='color:#666' title='",
               htmltools::htmlEscape(tooltip, attribute = TRUE),
               "'>", confText, "</small>")
      } else {
        msg
      }
      output$results <- renderUI(HTML(html))
      invisible(msg)
    }

    ##########################################################################
    # ExtendedTask for async search
    ##########################################################################

    # Cancel file path — created before each search, deleted on completion.
    # The C++ engine checks for this file's existence every ~200ms and stops
    # gracefully if it appears.
    cancelFile <- reactiveVal(NULL)
    # Progress file path — C++ callback writes per-replicate status here;
    # polled by an invalidateLater observer to update the notification.
    progressFile <- reactiveVal(NULL)

    searchTask <- ExtendedTask$new(
      function(dataset, tree, concavity, extendedIw, strategy,
               maxReplicates, targetHits, maxSeconds, poolSuboptimal,
               nThreads, cancelPath, progressPath,
               hierarchy, inapplicable, hsjAlpha) {
        future::future({
          on.exit({
            Sys.unsetenv("TREESEARCH_CANCEL_FILE")
            Sys.unsetenv("TREESEARCH_PROGRESS_FILE")
          })
          if (nzchar(cancelPath)) {
            Sys.setenv(TREESEARCH_CANCEL_FILE = cancelPath)
          }
          if (nzchar(progressPath)) {
            Sys.setenv(TREESEARCH_PROGRESS_FILE = progressPath)
          }
          args <- list(
            dataset,
            tree = tree,
            concavity = concavity,
            extended_iw = extendedIw,
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
          # Inapplicable handling (non-Brazeau requires hierarchy)
          if (!is.null(hierarchy) && !identical(inapplicable, "brazeau")) {
            args$hierarchy    <- hierarchy
            args$inapplicable <- inapplicable
            if (identical(inapplicable, "hsj")) {
              args$hsj_alpha <- hsjAlpha
            }
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
        return(invisible())
      }

      # Profile mode: defer search until profile data is prepared
      if (identical(concavity(), "profile") &&
          !identical(r$dataHash, profileDataHash())) {
        startProfilePrep(r$dataset)
        return(invisible())
      }

      # Read search parameters early (before any slow prep)
      searchStrategy  <- if (length(input$strategy)) input$strategy else "auto"
      searchMaxRep    <- if (length(input$maxReplicates)) {
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

      # Inapplicable handling
      searchInapplicable <- if (length(input$inapplicable)) input$inapplicable else "brazeau"
      searchHsjAlpha     <- if (length(input$hsjAlpha)) as.double(input$hsjAlpha) else 1.0
      searchHierarchy <- if (!identical(searchInapplicable, "brazeau") &&
                             !is.null(r$chars) && length(r$chars) > 0L) {
        tryCatch(
          withCallingHandlers(
            hierarchy_from_names(r$chars),
            warning = function(w) invokeRestart("muffleWarning")
          ),
          error = function(e) NULL
        )
      } else {
        NULL
      }

      # Non-Brazeau methods require a detected hierarchy; abort early
      if (!identical(searchInapplicable, "brazeau") && is.null(searchHierarchy)) {
        methodLabel <- switch(searchInapplicable,
                              hsj   = "Hopkins & St. John (HSJ)",
                              xform = "X-transformation (Goloboff)",
                              searchInapplicable)
        Notification(
          paste0(
            "The \u201c", methodLabel, "\u201d method requires a character ",
            "hierarchy. Ensure character names follow the sup_<tag> / ",
            "sub_<tag> convention (see ?hierarchy_from_names)."
          ),
          type = "error", duration = 10
        )
        return(invisible())
      }

      # Show search-in-progress indicator BEFORE tree selection (which may
      # call AdditionTree synchronously). The guard in DisplayTreeScores()
      # checks r$searchNotification to avoid overwriting this indicator.
      disable("go")
      disable("modalGo")
      disable("searchConfig")
      shinyjs::show("cancel")
      # Create unique temp file paths for cancel + progress signaling
      cancelPath <- tempfile("ts_cancel_", fileext = ".signal")
      cancelFile(cancelPath)
      progressPath <- tempfile("ts_progress_", fileext = ".txt")
      progressFile(progressPath)
      searchLabel <- paste0(
        "Searching (", searchMaxRep, " runs, ", wtType(),
        if (searchNThreads > 1L) paste0(", ", searchNThreads, " threads") else "",
        ")\u2026"
      )
      r$searchNotification <- showNotification(
        searchLabel, duration = NULL, type = "message", closeButton = FALSE
      )
      r$searchDataHash <- r$dataHash
      r$searchInProgress <- TRUE
      output$results <- renderUI(HTML(searchLabel))

      startTree <- tryCatch({
        if (!AnyTrees()) {
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
      }, error = function(e) {
        LogMsg("Starting tree error: ", conditionMessage(e), "; using fresh tree")
        LogCode(paste0("startTree <- AdditionTree(dataset, concavity = ",
                       Enquote(concavity()), ")"))
        AdditionTree(r$dataset[SearchTips()], concavity = concavity())
      })
      LogMsg("StartSearch()")
      PutData(r$dataset[SearchTips()])
      PutTree(startTree)
      LogComment("Search for optimal trees", 1)
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
        if (!searchExtendedIw && is.finite(searchConcavity))
          "  extended_iw = FALSE,",
        paste0("  strategy = \"", searchStrategy, "\","),
        paste0("  maxReplicates = ", searchMaxRep, ","),
        paste0("  targetHits = ", searchTargetHits, ","),
        if (searchMaxSeconds > 0)
          paste0("  maxSeconds = ", searchMaxSeconds, ","),
        if (searchPoolSub > 0)
          paste0("  control = SearchControl(poolSuboptimal = ", searchPoolSub, "),"),
        if (searchNThreads > 1L)
          paste0("  nThreads = ", searchNThreads, "L,"),
        if (!identical(searchInapplicable, "brazeau") && !is.null(searchHierarchy))
          paste0("  inapplicable = \"", searchInapplicable, "\","),
        if (identical(searchInapplicable, "hsj") && !is.null(searchHierarchy) &&
            searchHsjAlpha != 1.0)
          paste0("  hsj_alpha = ", searchHsjAlpha, ","),
        "  verbosity = 0",
        ")"))
      # Snapshot reactive values for the async task
      searchDataset <- r$dataset[SearchTips()]
      searchConcavity <- concavity()
      searchExtendedIw <- extendedIw()

      searchTask$invoke(
        searchDataset, startTree, searchConcavity, searchExtendedIw,
        searchStrategy, searchMaxRep, searchTargetHits,
        searchMaxSeconds, searchPoolSub, searchNThreads,
        cancelPath, progressPath,
        searchHierarchy, searchInapplicable, searchHsjAlpha
      )
    }

    ##########################################################################
    # Input observers
    ##########################################################################

    observeEvent(input$searchWithout, {
      r$searchWithout <- input$searchWithout
    }, ignoreInit = TRUE)

    observeEvent(input$go, StartSearch())
    observeEvent(input$modalGo, {
      removeModal()
      StartSearch()
    })

    # Cancel button: create the signal file so the C++ engine stops
    observeEvent(input$cancel, {
      cf <- cancelFile()
      if (!is.null(cf)) {
        file.create(cf)
        shinyjs::hide("cancel")
        # Remove search notification immediately so it doesn't linger
        if (!is.null(r$searchNotification)) {
          removeNotification(r$searchNotification)
          r$searchNotification <- NULL
        }
        output$results <- renderUI(HTML(
          "Stopping \u2014 waiting for current search phase to finish\u2026"
        ))
      }
    })

    # Poll progress file during search to update notification
    observe({
      pf <- progressFile()
      nid <- r$searchNotification
      if (is.null(pf) || is.null(nid) || !isTRUE(r$searchInProgress)) return()
      invalidateLater(500)
      if (!file.exists(pf)) return()  # C++ hasn't written first status yet
      progress <- tryCatch(
        readLines(pf, warn = FALSE),
        error = function(e) NULL
      )
      if (is.null(progress) || length(progress) == 0L ||
          !nzchar(progress[[1L]])) return()
      parts <- strsplit(progress[[1L]], " ", fixed = TRUE)[[1L]]
      if (length(parts) < 5L) return()
      rep_cur   <- parts[1L]
      rep_max   <- parts[2L]
      best      <- parts[3L]
      hits      <- parts[4L]
      target    <- parts[5L]
      msg <- paste0(
        "Searching\u2026 Rep ", rep_cur, "/", rep_max,
        " | Best: ", best,
        " | Hits: ", hits, "/", target
      )
      # Update both the results area and the toast (belt-and-suspenders: if
      # DisplayTreeScores() was called and overwrote output$results, the next
      # poll restores the progress message within 500 ms).
      output$results <- renderUI(HTML(msg))
      showNotification(msg, id = nid, duration = NULL,
                       type = "message", closeButton = FALSE)
    })

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
      # Sync inapplicable selector and show/hide hsjAlpha accordingly
      inapplicable_cur <- if (length(input$inapplicable)) input$inapplicable else "brazeau"
      updateSelectInput(session, "inapplicable", selected = inapplicable_cur)
      updateNumericInput(session, "hsjAlpha",
                         value = if (length(input$hsjAlpha)) input$hsjAlpha else 1.0)
      if (identical(inapplicable_cur, "hsj")) show("hsjAlpha") else hide("hsjAlpha")
      # Initialise all modal inputs from current values so that opening the
      # modal does not fire observeEvent(input$concavity) or
      # observeEvent(input$implied.weights), which reset the run counters.
      cur_weights   <- if (length(input$implied.weights)) input$implied.weights else "on"
      cur_concavity <- if (length(input$concavity))       input$concavity       else 1L
      cur_strategy  <- if (length(input$strategy))        input$strategy        else "auto"
      cur_maxRep    <- if (length(input$maxReplicates))   input$maxReplicates   else 100L
      cur_hits      <- if (length(input$targetHits))      input$targetHits      else 10L
      cur_timeout   <- if (length(input$timeout))         input$timeout         else 5
      cur_epsilon   <- if (length(input$epsilon))         input$epsilon         else 0
      cur_threads   <- if (length(input$nThreads))        input$nThreads        else max(1L, floor(nCores / 2L))
      showModal(modalDialog(
        easyClose = TRUE,
        fluidPage(column(6,
          selectInput(ns("implied.weights"), "Step weighting",
                     list("Implied (extended)" = "xpiwe",
                          "Implied" = "on", "Profile" = "prof",
                          "Equal" = "off"), cur_weights),
          sliderInput(ns("concavity"), "Concavity constant", min = 0L,
                     max = 3L, pre = "10^", value = cur_concavity),
          selectInput(ns("inapplicable"), "Inapplicable characters",
                      list("Brazeau et al. (default)" = "brazeau",
                           "Hopkins & St. John (HSJ)"  = "hsj",
                           "X-transformation (Goloboff)" = "xform"),
                      inapplicable_cur),
          hidden(numericInput(ns("hsjAlpha"), "HSJ \u03b1 parameter",
                              value = if (length(input$hsjAlpha)) input$hsjAlpha else 1.0,
                              min = 0, step = 0.1)),
          uiOutput(ns("hierarchyInfo")),
          if (nCores > 1L) {
            sliderInput(ns("nThreads"), "Parallel search threads",
                        min = 1L, max = nCores,
                        value = cur_threads,
                        step = 1L)
          },
          selectizeInput(ns("searchWithout"), "Exclude taxa", DatasetTips(),
                         r$searchWithout, multiple = TRUE),
          numericInput(ns("epsilon"), "Keep if suboptimal by \u2264", min = 0,
                      value = cur_epsilon)
        ), column(6,
          selectInput(ns("strategy"), "Search strategy",
                     list("Auto" = "auto", "Sprint" = "sprint",
                          "Default" = "default", "Thorough" = "thorough"),
                     cur_strategy),
          sliderInput(ns("targetHits"),
                      "Stop when N runs have hit best score",
                      min = 1L, max = 50L, value = cur_hits, step = 1L),
          uiOutput(ns("targetHitsNote")),
          sliderInput(ns("timeout"), "Maximum run duration", min = 1,
                      max = 60, value = cur_timeout, post = "min", step = 1),
          sliderInput(ns("maxReplicates"), "Maximum independent runs",
                      min = 1L, max = 500L, value = cur_maxRep, step = 1L),
          helpText("Limits each individual search. Clicking \u2018Continue\u2019",
                   "starts a fresh search; the results panel shows the",
                   "cumulative total across all continued searches.")
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
        shiny.silent.error = function(e) {
          # ExtendedTask throws shiny.silent.error (class "shiny.output.progress"
          # subclass) when status is "initial" or "running" — not a real error.
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
        # Clean up search-in-progress UI state. Gate on searchInProgress
        # (not searchNotification) because the cancel observer may have
        # already dismissed the notification.
        if (isTRUE(r$searchInProgress)) {
          if (!is.null(r$searchNotification)) {
            removeNotification(r$searchNotification)
            r$searchNotification <- NULL
          }
          enable("go")
          enable("modalGo")
          enable("searchConfig")
          shinyjs::hide("cancel")
          cf <- cancelFile()
          if (!is.null(cf)) {
            suppressWarnings(file.remove(cf))
            cancelFile(NULL)
          }
          pf <- progressFile()
          if (!is.null(pf)) {
            suppressWarnings(file.remove(pf))
            progressFile(NULL)
          }
          r$searchInProgress <- FALSE
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
        # Always refresh the display — UpdateAllTrees may short-circuit
        # when trees are unchanged, but hit/rep counts have been updated.
        DisplayTreeScores()
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
      r$bestSearchScore <- NULL
      r$searchCount <- 0L
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
      extendedIw        = extendedIw,
      DisplayTreeScores = DisplayTreeScores
    )
  })
}
