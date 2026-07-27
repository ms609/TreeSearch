  serverEnv <- environment()
  logIndent <- 0
  loggingOn <- TRUE
  
  cmdLogFile <- tempfile("TreeSearch-", fileext = ".R")
  Write <- function (txt, file) {
    if (serverEnv$loggingOn) {
      txt <- paste0(strrep(" ", logIndent), txt)
      con <- file(file, open = "a")
      on.exit(close(con))
      if (logging) {
        WriteLoggedCode(txt)
      }
      writeLines(txt, con)
    }
  }
  
  WriteP <- function (txt, file = NULL) {
    if (serverEnv$loggingOn) {
      txt <- paste0(strrep(" ", logIndent), txt)
      if (logging) {
        WriteLoggedCode(txt)
      }
      r$plotLog <- c(r$plotLog, as.character(txt))
    }
  }
  
  LogExpr <- function(exps, evaluate = TRUE, WriteFn = Write) {
    for (exp in exps) {
      WriteFn(as.character(exp), cmdLogFile)
      if (evaluate) {
        eval(exp)
      }
    }
  }
  
  LogExprP <- function(...) {
    LogExpr(..., WriteFn = WriteP)
  }
  
  LogIndent <- function(n) {
    serverEnv$logIndent <- serverEnv$logIndent + n
    if (serverEnv$logIndent < 0) {
      warning("Negative indent")
    }
  }
  
  systemInfo <- c(
    paste(
      "System:", Sys.info()["sysname"], Sys.info()["release"],
      Sys.info()["version"], "-",
      .Platform$OS.type, R.version$platform
    ),
    paste(
      "-", R.version$version.string
    ),
    paste("- TreeSearch", packageVersion("TreeSearch")),
    paste("- TreeTools", packageVersion("TreeTools")),
    paste("- TreeDist", packageVersion("TreeDist")),
    paste("- ape", packageVersion("ape"))
  )
  
  logCaveats <- c(
    "Before running, check that the script and any data files are in the",
    "R working directory, which can be read with getwd() and set with setwd().",
    "",
    "Please validate the code before reproducing in a manuscript, reporting",
    "any errors at https://github.com/ms609/treesearch/issues or by e-mail to",
    "the package maintainer."
  )
  
  BeginLog <- function() {
    LogComment(c(
      paste("# # TreeSearch session log:", .DateTime(), "# # #"),
      "",
      systemInfo,
      "",
      "This log was generated procedurally to facilitate the reproduction of",
      "results obtained during an interactive Shiny session.",
      "It is provided without guarantee of completeness or accuracy.",
      "In particular, code will not be logged when previously computed values",
      "are retrieved from cache.",
      "",
      logCaveats,
      "",
      "# # # # #"
    ))
    
    LogComment("Load required libraries", 2)
    LogCode(c(
      "library(\"TreeTools\", quietly = TRUE)",
      "library(\"TreeDist\")",
      "library(\"TreeSearch\")"
    ))
    
    LogComment("View recommended citations", 1)
    LogCode(c(
      "citation(\"TreeTools\")",
      "citation(\"TreeDist\")",
      "citation(\"TreeSearch\")",
      "citation(\"Rogue\")"
    ))
  }
  
  BeginLogP <- function() {
    r$plotLog <- NULL
    LogCommentP(c(
      paste("# # TreeSearch plot log:", .DateTime(), "# # #"),
      "",
      systemInfo,
      "",
      "This log was generated procedurally to facilitate the reproduction of",
      "figures obtained during an interactive Shiny session.",
      "It is provided without guarantee of completeness or accuracy.",
      "In particular, code will not be logged when previously computed values",
      "are retrieved from cache.",
      "",
      logCaveats,
      "",
      "# # # # #"
    ))
    LogCommentP("Load required libraries", 2)
    LogCodeP(c(
      "library(\"TreeTools\", quietly = TRUE)",
      "library(\"TreeDist\")",
      "library(\"TreeSearch\")"
    ))
    
    LogCommentP("View recommended citations", 1)
    LogCodeP(c(
      "citation(\"TreeTools\")",
      "citation(\"TreeDist\")",
      "citation(\"Quartet\")",
      "citation(\"TreeSearch\")",
      "citation(\"Rogue\")"
    ))
    
    LogCommentP("Check working directory", 1)
    LogCodeP("getwd() # Should match location of data / tree files",
             "setwd(\".\") # Replace . with desired/directory to change")
    
    if (HaveData()) {
      LogCommentP("Load data from file")
      LogCodeP(c(
        paste0("dataFile <- ", Enquote(DataFileName(r$dataFiles))),
        paste0("dataset <- ", r$readDataFile)
      ))
    }
    
    if (AnyTrees()) {
      LogCommentP("Load trees from file")
      LogCodeP(c(
        paste0("treeFile <- ", Enquote(TreeFileName(r$treeFiles))),
        "trees <- read.nexus(treeFile)",
        if (!identical(r$trees, r$allTrees)) {
          c(
            if (!is.null(r$thinningSeed)) {
              paste0("set.seed(", r$thinningSeed, ")")
            },
            paste0(
              "trees <- WideSample(trees[",
              r$treeRange[1], ":", r$treeRange[2],
              "], ", r$nTree, ")"
            )
          )
        }
      ))
    }
  }
  
  PauseLog <- function() {
    serverEnv$loggingOn <- FALSE
  }
  
  ResumeLog <- function() {
    serverEnv$loggingOn <- TRUE
  }
  
  LogCode <- function(..., WriteFn = Write) {
    for (line in list(...)) {
      if (!is.null(line)) {
        WriteFn(as.character(line), cmdLogFile)
      }
    }
  }
  
  LogCodeP <- function(...) {
    LogCode(..., WriteFn = WriteP)
  }
  
  LogComment <- function(exps, returns = 1, WriteFn = Write) {
    if (returns > 0) {
      WriteFn(rep("", returns), cmdLogFile)
    }
    for (exp in exps) {
      WriteFn(paste("#", exp), cmdLogFile)
    }
  }
  
  LogCommentP <- function (exps, returns = 1) {
    LogComment(exps, returns, WriteFn = WriteP)
  }
  
  r$dataFiles <- 0
  r$excelFiles <- 0
  r$treeFiles <- 0
  TwoWide <- function(n) {
    formatC(n, width = 2, flag = "0")
  }
  # tempdir() is shared by every tab of a single running app process, but
  # these counters reset to 0 per-session — so two tabs' first uploads both
  # produced "treeFile-01.txt" and silently clobbered each other. Namespace
  # by session token so concurrent tabs never collide (T-356).
  sessionTag <- gsub("[^A-Za-z0-9]", "", session$token)
  DataFileName <- function(n) if (length(n)) {
    paste0("dataFile-", sessionTag, "-", TwoWide(n), ".txt")
  }
  ExcelFileName <- function(n) if (length(n)) {
    paste0("excelFile-", sessionTag, "-", TwoWide(n), ".xlsx")
  }
  TreeFileName <- function(n) if (length(n)) {
    paste0("treeFile-", sessionTag, "-", TwoWide(n), ".txt")
  }
  LastFile <- function(type) {
    switch(pmatch(type, c("data", "excel", "tree")), 
           DataFileName(r$dataFiles),
           ExcelFileName(r$excelFiles),
           TreeFileName(r$treeFiles)
    )
  }
  CacheInput <- function(type, fileName) {
    key <- paste0(type, "Files")
    r[[key]] <- r[[key]] + 1
    file.copy(fileName, paste0(tempdir(), "/", LastFile(type)),
              overwrite = TRUE)
  }
  StashTrees <- function(trees) {
    key <- paste0("treeFiles")
    r[[key]] <- r[[key]] + 1
    write.nexus(trees, file = paste0(tempdir(), "/", LastFile("tree")))
  }
  
  BeginLog()
