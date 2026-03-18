server <- function(input, output, session) {
  
  source("server/app_state.R", local = TRUE)
  r <- AppState()
  exportTestValues(searchCount = { r$searchCount })
  
  source("server/logging.R", local = TRUE)
  source("server/search.R", local = TRUE)
  source("server/data.R", local = TRUE)
  source("server/trees.R", local = TRUE)
  source("server/events.R", local = TRUE)
  source("server/consensus.R", local = TRUE)
  source("server/clustering.R", local = TRUE)
  ts <- treespace_server("treespace",
    r = r,
    clusterings = clusterings,
    silThreshold = silThreshold,
    scores = scores,
    concavity = concavity,
    distMeth = reactive(input$distMeth),
    plotFormat = reactive(input$plotFormat),
    log_fns = list(
      BeginLogP      = BeginLogP,
      LogCommentP    = LogCommentP,
      LogCodeP       = LogCodeP,
      LogIndent      = LogIndent,
      LogClusterings = LogClusterings
    )
  )
  # Expose module reactives to source'd files (clustering, consensus, downloads)
  distances        <- ts$distances
  mapping          <- ts$mapping
  dims             <- ts$dims
  nProjDim         <- ts$nProjDim
  TreeCols         <- ts$TreeCols
  treePch          <- ts$treePch
  mstEnds          <- ts$mstEnds
  saveDetails      <- ts$saveDetails
  TreespacePlot    <- ts$TreespacePlot
  LogTreespacePlot <- ts$LogTreespacePlot
  LogDistances     <- ts$LogDistances
  downloads_server(
    "dl",
    state         = r,
    dataSource    = reactive(input$dataSource),
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
