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
  # Expose data module reactives for source'd files + other modules
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
  
  source("server/events.R", local = TRUE)
  source("server/consensus.R", local = TRUE)
  
  # Wire forward-reference callbacks (events.R + consensus.R now defined)
  cb_ref$DisplayTreeScores       <- DisplayTreeScores
  cb_ref$UpdateKeepNTipsRange    <- UpdateKeepNTipsRange
  cb_ref$UpdateDroppedTaxaDisplay <- UpdateDroppedTaxaDisplay
  cb_ref$UpdateOutgroupInput     <- UpdateOutgroupInput
  
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
  # Expose clustering module reactives for source'd files
  distances      <- cl$distances
  LogDistances   <- cl$LogDistances
  silThreshold   <- cl$silThreshold
  clusterings    <- cl$clusterings
  LogClusterings <- cl$LogClusterings
  source("server/clustering.R", local = TRUE)
  ts <- treespace_server("treespace",
    r = r,
    clusterings = clusterings,
    silThreshold = silThreshold,
    scores = scores,
    concavity = concavity,
    distMeth = reactive(input$distMeth),
    plotFormat = reactive(input$plotFormat),
    distances = distances,
    LogDistances = LogDistances,
    log_fns = list(
      BeginLogP      = BeginLogP,
      LogCommentP    = LogCommentP,
      LogCodeP       = LogCodeP,
      LogIndent      = LogIndent,
      LogClusterings = LogClusterings
    )
  )
  # Expose treespace module reactives to source'd files (consensus, downloads)
  mapping          <- ts$mapping
  dims             <- ts$dims
  nProjDim         <- ts$nProjDim
  TreeCols         <- ts$TreeCols
  treePch          <- ts$treePch
  mstEnds          <- ts$mstEnds
  saveDetails      <- ts$saveDetails
  TreespacePlot    <- ts$TreespacePlot
  LogTreespacePlot <- ts$LogTreespacePlot
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
    mainPlot      = MainPlot,
    rCode         = RCode,
    saveDetails   = saveDetails
  )
  references_server("refs")
  
  onStop(function() {
    options(startOpt)
    if (file.exists(cmdLogFile)) {
      unlink(cmdLogFile)
    }
    # Clean cached input files from tempdir (data, tree, and excel)
    unlink(list.files(tempdir(), pattern = "^(data|tree|excel)File-",
                      full.names = TRUE))
    if (logging) {
      LogMsg("Session has ended")
      on.exit(close(logMsgFile))
    }
  })
}
