  output$saveZip <- downloadHandler(
    filename = function() paste0("TreeSearch-session.zip"),
    content = function(file) {
      if (isTRUE(getOption("shiny.testmode"))) {
        file.copy(cmdLogFile, file)
      } else {
        zipDir <- tempfile("zip-")
        dir.create(zipDir)
        on.exit(unlink(zipDir))
        rFile <- paste0(zipDir, "/TreeSearch-session.R")
        file.copy(cmdLogFile, rFile, overwrite = TRUE)
        zip(file, c(
          rFile,
          if (r$dataFiles)
            paste0(tempdir(), "/", DataFileName(seq_len(r$dataFiles))),
          if (r$excelFiles)
            paste0(tempdir(), "/", ExcelFileName(seq_len(r$excelFiles))),
          if (r$treeFiles)
            paste0(tempdir(), "/", TreeFileName(seq_len(r$treeFiles)))
        ), flags = "-9Xj")
      }
    })
  
  output$savePlotZip <- downloadHandler(
    filename = function() paste0(saveDetails()$fileName, ".zip"),
    content = function(file) {
      StashTrees(r$allTrees)
      
      if (isTRUE(getOption("shiny.testmode"))) {
        rCode <- RCode()
        rCode <- sub("TreeSearch plot log: 2[\\d\\-]{9} [012][\\d:]{7}",
                     "TreeSearch plot log: <DATE-AND-TIME>", 
                     rCode, perl = TRUE)
        rCode[4] <- "# System: <SYS-INFO>"
        rCode[5:9] <- sub("^(# \\- \\w+ ).*$", "\\1<VERSION>",
                          rCode[5:9], perl = TRUE)
        rCode <- sub("dataFile <- .*$",
                     paste0("dataFile <- system.file(\"datasets/",
                            input$dataSource,
                            ".nex\", package = \"TreeSearch\") # FALSE CODE for TEST MODE"),
                     rCode,
                     perl = TRUE)
        rCode <- sub("treeFile <- .*$",
                     "treeFile <- dataFile # Test mode",
                     rCode,
                     perl = TRUE)
        writeLines(rCode, con = file)
      } else {
        tempDir <- tempfile("plot-zip-")
        dir.create(tempDir)
        on.exit(unlink(tempDir))
        rFile <- paste0(tempDir, "/", saveDetails()$fileName, ".R")
        writeLines(RCode(), con = rFile)
        
        # Create ZIP
        zip(file, c(
          rFile,
          paste0(tempdir(), "/", LastFile("data")),
          paste0(tempdir(), "/", LastFile("excel")),
          paste0(tempdir(), "/", LastFile("tree"))
        ), flags = "-r9Xj")
      }
    })
  
  output$savePng <- downloadHandler(
    filename = function() paste0(saveDetails()$fileName, ".png"),
    content = function (file) {
      png(file, width = input$plotSize, height = input$plotSize)
      MainPlot()
      dev.off()
    })
  
  output$savePdf <- downloadHandler(
    filename = function() paste0(saveDetails()$fileName, ".pdf"),
    content = function (file) {
      width <- 8
      pdf(
        file,
        title = saveDetails()$title,
        width = width,
        height = saveDetails()$asp * width
      )
      MainPlot()
      dev.off()
    })
  
  output$savePlotNwk <- downloadHandler(
    filename = "TreeSearch-consensus.nwk",
    content = function(file) {
      write.tree(r$plottedTree, file = file)
    }
  )
  
  output$savePlotNex <- downloadHandler(
    filename = "TreeSearch-consensus.nex",
    content = function(file) {
      write.nexus(r$plottedTree, file = file)
    }
  )
  
  output$saveNwk <- downloadHandler(
    filename = "TreeSearch.nwk",
    content = function(file) {
      write.tree(r$trees, file = file, tree.names = TRUE)
    }
  )
  
  output$saveNex <- downloadHandler(
    filename = "TreeSearch.nex",
    content = function(file) {
      write.nexus(r$trees, file = file)
    }
  )
