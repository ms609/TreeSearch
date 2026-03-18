  AnyTrees <- reactive({!is.null(r$trees) && length(r$trees) > 0})
  HaveData <- reactive({!is.null(r$dataset) && length(r$dataset) > 0 && inherits(r$dataset, "phyDat")})
  FetchNTree <- debounce(reactive({
    if (!is.null(r$oldNTree)) {
      if (!identical(input$nTree, r$oldNTree)) {
        r$oldNTree <- NULL
      }
    } else {
      if (UpdateNTree(input$nTree)) {
        UpdateActiveTrees()
      }
    }
  }), typingJiffy)
  
  # Return TRUE if n has changed, FALSE if not
  # Don't update active trees here: Leave this to the calling function
  UpdateNTree <- function(n) {
    if (n > length(r$allTrees)) { # nTree "max" can be beaten by typing
      r$oldNTree <- n
      n <- length(r$allTrees)
    }
    if (r$nTree == n) {
      # Return:
      FALSE
    } else {
      LogMsg("UpdateNTree(", r$nTree, " -> ", n, ")")
      r$nTree <- n
      # range <- r$treeRange[2] - r$treeRange[1]
      # if (n > range + 1L) {
      #   nTrees <- length(r$allTrees)
      #   upper <- min(nTrees, r$treeRange[1] + n - 1L)
      #   lower <- min(r$treeRange[1], upper + 1L - n)
      #   r$treeRange <- c(lower, upper)
      #   updateSliderInput(session, "treeRange", value = r$treeRange)
      # }
      if (input$nTree != n) {
        updateNumericInput(session, "nTree", value = n)
      }
      # Return:
      TRUE
    }
  }
  
  FetchTreeRange <- debounce(reactive({
    if (!is.null(r$oldTreeRange)) {
      if (!identical(input$treeRange, r$oldTreeRange)) {
        r$oldTreeRange <- NULL
      }
    } else {
      if (UpdateTreeRange(input$treeRange)) {
        UpdateActiveTrees()
      }
    }
  }), aJiffy)
  
  # Return TRUE if changed, FALSE if not
  # Don't update active trees here: Leave this to the calling function
  UpdateTreeRange <- function(range) {
    if (identical(range, r$treeRange)) {
      # Return:
      FALSE
    } else {
      LogMsg("UpdateTreeRange([", paste(r$treeRange, collapse = ", "),
             "] -> [", paste(range, collapse = ", "), "])")
      r$treeRange <- range
      span <- r$treeRange[2] - r$treeRange[1]
      if (r$nTree > span + 1L) {
        UpdateNTree(span + 1L)
      }
      
      # Return:
      TRUE
    }
  }
  
  
  UpdateActiveTrees <- reactive({
    if (r$updatingTrees) {
      LogMsg("   Skipping UpdateActiveTrees()")
      return()
    }
    r$updatingTrees <- TRUE
    on.exit(r$updatingTrees <- FALSE)
    LogMsg("UpdateActiveTrees()")
    
    nTrees <- length(r$allTrees)
    if (r$nTree == nTrees &&
        r$treeRange[1] == 1L && r$treeRange[2] == nTrees) {
      thinnedTrees <- r$allTrees
      if (!is.null(r$allTrees) && !identical(trees, thinnedTrees)) {
        LogCode("trees <- allTrees")
      }
    } else {
      thinnedTrees <- r$allTrees[
        unique(as.integer(seq.int(
          r$treeRange[1], r$treeRange[2], length.out = r$nTree)))]
      
      if (!is.null(r$allTrees) && !identical(trees, thinnedTrees)) {
        LogCode(paste0(
          "trees <- allTrees[unique(as.integer(seq.int(",
          r$treeRange[1], ", ", r$treeRange[2], ", length.out = ", r$nTree, ")))]"
        ))
      }
    }
    
    r$trees <- thinnedTrees
    r$treeHash <- rlang::hash(r$trees)
    
    DisplayTreeScores()
    
    if (AnyTrees()) {
      for (elem in c("keepNTips", "neverDrop")) {
        showElement(elem, anim = TRUE)
      }
    } else {
      for (elem in c("keepNTips", "neverDrop")) {
        hideElement(elem)
      }
    }
    
    updateSliderInput(session, "whichTree", min = 0L,
                      max = length(r[["trees"]]), value = 0L)
    UpdateKeepNTipsRange() # Updates Rogues()
    UpdateDroppedTaxaDisplay()
    if (maxProjDim() > 0) {
      updateSliderInput(inputId = "treespace-spaceDim",
                        max = max(1L, maxProjDim()),
                        value = min(maxProjDim(),
                                    input[["treespace-spaceDim"]]))
    }
    updateSelectizeInput(inputId = "neverDrop", choices = tipLabels(),
                         selected = input$neverDrop)
    UpdateOutgroupInput()
    updateSelectizeInput(inputId = "treespace-relators",
                         choices = tipLabels(),
                         selected = input[["treespace-relators"]])
  })
  
  UpdateAllTrees <- function (newTrees) {
    LogMsg("UpdateAllTrees()")
    on.exit({
      LogMsg("/UpdateAllTrees()")
    }, add = TRUE)
    
    newTrees <- c(newTrees)
    if (length(newTrees) > 1L) {
      newTrees <- RenumberTips(newTrees, newTrees[[1]]$tip.label)
    }
    if (identical(newTrees, r$newTrees)) {
      LogMsg("   <Trees unchanged; returning>")
      return()
    }
    r$newTrees <- newTrees
    
    oldNTrees <- length(r$allTrees)
    
    if (!identical(r$allTrees, newTrees)) {
      LogCode("allTrees <- newTrees")
      r$allTrees <- newTrees
    }
    nTrees <- length(newTrees)
    
    if (nTrees != oldNTrees) {
      if (!identical(input$treeRange, c(1L, nTrees))) {
        r$oldTreeRange <- input$treeRange
      }
      UpdateTreeRange(c(1L, nTrees))
      # update*Input messages are collected and sent after all the observers
      # (including outputs) have finished running.
      updateSliderInput(session, "treeRange",
                        min = 1L, max = nTrees,
                        value = r$treeRange)
    
      r$oldNTree <- input$nTree
      UpdateNTree(min(max(input$nTree, aFewTrees), nTrees))
      updateNumericInput(session, "nTree", max = nTrees,
                         value = r$nTree)
    }
    
    UpdateActiveTrees()
    if (AnyTrees()) {
      showElement("manipulateTreeset")
    } else {
      hideElement("manipulateTreeset")
    }
  }
  
