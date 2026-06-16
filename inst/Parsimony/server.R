server <- function(input, output, session) {
  
  source("server/app_state.R", local = TRUE)
  r <- AppState()
  exportTestValues(searchCount = { r$searchCount })
  
  # Async search setup (was in search.R)
  plan(multisession)
  startOpt <- options("cli.progress_show_after" = 0.1)
  
  source("server/logging.R", local = TRUE)
  LogMsg("Started server")
  
  # Forward-reference bridge for callbacks defined after the data module
  cb_ref <- new.env(parent = emptyenv())
  
  # Data module (replaces data.R + trees.R + data event bindings from events.R)
  dt <- data_server("data",
    r = r,
    parent_session = session,
    callbacks = list(
      DisplayTreeScores       = function() cb_ref$DisplayTreeScores(),
      UpdateKeepNTipsRange    = function() cb_ref$UpdateKeepNTipsRange(),
      UpdateDroppedTaxaDisplay = function() cb_ref$UpdateDroppedTaxaDisplay(),
      UpdateOutgroupInput     = function() cb_ref$UpdateOutgroupInput()
    ),
    log_fns = list(
      LogMsg     = LogMsg,
      LogComment = LogComment,
      LogCode    = LogCode,
      CacheInput = CacheInput,
      LastFile   = LastFile
    )
  )
  # Expose data module reactives for other modules
  AnyTrees            <- dt$AnyTrees
  HaveData            <- dt$HaveData
  tipLabels           <- dt$tipLabels
  nChars              <- dt$nChars
  TaxonOrder          <- dt$TaxonOrder
  DatasetMatchesTrees <- dt$DatasetMatchesTrees
  UpdateAllTrees      <- dt$UpdateAllTrees
  UpdateActiveTrees   <- dt$UpdateActiveTrees
  
  # Search module
  se <- search_server("search",
    r = r,
    AnyTrees = AnyTrees,
    HaveData = HaveData,
    UpdateAllTrees = UpdateAllTrees,
    log_fns = list(
      LogMsg     = LogMsg,
      LogCode    = LogCode,
      LogComment = LogComment
    )
  )
  scores            <- se$scores
  concavity         <- se$concavity
  DisplayTreeScores <- se$DisplayTreeScores

  # Show/hide config panels based on active plot format
  ShowConfigs <- function(visible = character(0)) {
    allConfigs <- c("whichTree", "charChooser",
                    "consConfig", "clusConfig",
                    "clusLegend", "branchLegend",
                    "spaceConfig", "treePlotConfig",
                    "mapConfig", "savePlottedTrees",
                    "droppedTips", "droppedList")
    r$visibleConfigs <- visible
    lapply(visible, show)
    lapply(setdiff(allConfigs, visible), hide)
  }

  observeEvent(input$plotFormat, {
    ShowConfigs(switch(input$plotFormat,
      "ind"   = c("whichTree", "charChooser", "treePlotConfig"),
      "cons"  = c("consConfig", "droppedTips", "savePlottedTrees",
                   "treePlotConfig", "branchLegend"),
      "clus"  = c("clusConfig", "clusLegend", "savePlottedTrees",
                   "consConfig", "droppedList", "treePlotConfig"),
      "space" = c("clusConfig", "clusLegend", "spaceConfig", "mapConfig"),
      ""))
  })

  # Clustering module
  cl <- clustering_server("clustering",
    r = r,
    distMeth = reactive(input$distMeth),
    log_fns = list(
      LogMsg      = LogMsg,
      LogCommentP = LogCommentP,
      LogCodeP    = LogCodeP,
      LogIndent   = LogIndent,
      BeginLogP   = BeginLogP,
      LogExprP    = LogExprP
    )
  )
  distances      <- cl$distances
  LogDistances   <- cl$LogDistances
  silThreshold   <- cl$silThreshold
  clusterings    <- cl$clusterings
  LogClusterings <- cl$LogClusterings

  # Treespace module
  ts <- treespace_server("treespace",
    r = r,
    clusterings = clusterings,
    silThreshold = silThreshold,
    scores = scores,
    concavity = concavity,
    distMeth = reactive(input$distMeth),
    plotFormat = reactive(input$plotFormat),
    distances = distances,
    mapLines = reactive(input$mapLines),
    LogDistances = LogDistances,
    log_fns = list(
      BeginLogP      = BeginLogP,
      LogCommentP    = LogCommentP,
      LogCodeP       = LogCodeP,
      LogIndent      = LogIndent,
      LogClusterings = LogClusterings
    )
  )
  saveDetails <- ts$saveDetails

  # Consensus module (replaces consensus.R + clustering.R + events.R bindings)
  co <- consensus_server("consensus",
    r = r,
    AnyTrees = AnyTrees,
    HaveData = HaveData,
    tipLabels = tipLabels,
    nChars = nChars,
    TaxonOrder = TaxonOrder,
    concavity = concavity,
    clusterings = clusterings,
    silThreshold = silThreshold,
    LogClusterings = LogClusterings,
    TreespacePlot    = ts$TreespacePlot,
    LogTreespacePlot = ts$LogTreespacePlot,
    dims       = ts$dims,
    nProjDim   = ts$nProjDim,
    TreeCols   = ts$TreeCols,
    treePch    = ts$treePch,
    ts_spaceCol  = ts$spaceCol,
    ts_mapLines  = ts$mapLines,
    ts_spacePch  = ts$spacePch,
    ts_relators  = ts$relators,
    plotFormat = reactive(input$plotFormat),
    plotSize  = reactive(input$plotSize),
    distMeth  = reactive(input$distMeth),
    log_fns = list(
      LogMsg      = LogMsg,
      LogComment  = LogComment,
      LogCode     = LogCode,
      LogCommentP = LogCommentP,
      LogCodeP    = LogCodeP,
      LogIndent   = LogIndent,
      BeginLogP   = BeginLogP,
      LogExprP    = LogExprP
    )
  )

  # Wire forward-reference callbacks (consensus module now defined)
  cb_ref$DisplayTreeScores       <- DisplayTreeScores
  cb_ref$UpdateKeepNTipsRange    <- co$UpdateKeepNTipsRange
  cb_ref$UpdateDroppedTaxaDisplay <- co$UpdateDroppedTaxaDisplay
  cb_ref$UpdateOutgroupInput     <- co$UpdateOutgroupInput

  # Downloads module
  downloads_server(
    "dl",
    state         = r,
    dataSource    = dt$dataSource,
    plotSize      = reactive(input$plotSize),
    cmdLogFile    = cmdLogFile,
    stashTrees    = StashTrees,
    dataFileName  = DataFileName,
    excelFileName = ExcelFileName,
    treeFileName  = TreeFileName,
    lastFile      = LastFile,
    mainPlot      = co$MainPlot,
    rCode         = co$RCode,
    saveDetails   = saveDetails
  )
  references_server("refs", weighting = se$weighting)
  
  onStop(function() {
    options(startOpt)
    if (file.exists(cmdLogFile)) {
      unlink(cmdLogFile)
    }
    # Clean cached input files from tempdir (data, tree, and excel)
    unlink(list.files(tempdir(), pattern = "^(data|tree|excel)File-",
                      full.names = TRUE))
    # T-312: also remove search/profile cancel + progress signal files; the
    # pattern above does not match them, so they otherwise leak on error /
    # interrupt / disconnect paths and accumulate across searches.
    unlink(list.files(tempdir(),
                      pattern = "^ts_(cancel|progress|profile_prog|profile_cancel)_",
                      full.names = TRUE))
    if (logging) {
      LogMsg("Session has ended")
      on.exit(close(logMsgFile))
    }
  })
}
