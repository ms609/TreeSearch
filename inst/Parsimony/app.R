library("methods", exclude = c('show', 'removeClass'))
library("shiny", exclude = c('runExample'))
suppressPackageStartupMessages(library("shinyjs", exclude = c('runExample')))
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeSearch")
library("TreeDist")

if (!requireNamespace('TreeDist', quietly = TRUE)) {
  install.packages('TreeDist')
  library('TreeDist')
}

palettes <- list("#91aaa7",
                 c("#497a8f", "#c1d0c0"),
                 c("#be83ae", "#2ea7af", "#fbcdcf"),
                 c("#72c5a9", "#b7c5ff", "#dbdabb", "#a28bac"),
                 c("#59c4c0", "#ea9a9a", "#7998a6", "#e9d7a9", "#9c9379"),
                 c("#e8b889", "#67c6eb", "#e5d5c2", "#938fba", "#64b69a", "#779c7b"),
                 c("#c4808f", "#5ba08f", "#f0a693", "#ccd7fe", "#cdb87e", "#c6aae2", "#d2dad8"),
                 c("#d0847f", "#63a5d7", "#d7b981", "#5a9bbb", "#9bb67e", "#dea6d5", "#91967e", "#ca7f96"),
                 c("#8b93a8", "#ccb97e", "#8e9dd7", "#57a384", "#dbb1e7", "#2da7af", "#d68986", "#75d2f9", "#e4d1f0"),
                 c("#dfcf92", "#40b3cb", "#b88a61", "#ecb2e0", "#d6dbbc", "#a28bae", "#edcfeb", "#7498ab", "#b187a0", "#8f939c"),
                 c("#a98f70", "#7be5e2", "#d295c0", "#9ae2bd", "#d3b7f1", "#eca88d", "#8993cd", "#ffc7bb", "#8695a8", "#b3e1df", "#b6878a"),
                 c("#eaa9d3", "#7ac09b", "#fdaca8", "#8ce7e4", "#eed69b", "#70a4d9", "#e8d6ba", "#589bbb", "#959672", "#d0dbd1", "#9b9282", "#d9d9c6"),
                 c("#7498ab", "#e5bd8a", "#7ed8ff", "#de8f8e", "#46bac6", "#ffc0d3", "#57ae96", "#f7cddd", "#50a098", "#b58a6d", "#add49d", "#a18da1", "#cedad9"),
                 c("#8097a4", "#d0dea9", "#a78cc3", "#aee4bf", "#bb82a8", "#5dc9c6", "#b88690", "#26a3b9", "#ad8e6f", "#a4e2ef", "#869a65", "#efcfdd", "#60a089", "#9e927b"),
                 c("#b9aae5", "#bbd69c", "#e2adde", "#77a777", "#f8abc8", "#8ee7ce", "#f2a1a5", "#81bdf1", "#f2bb91", "#b8dcfe", "#aeb276", "#f2cdef", "#e8d6b2", "#8d92b0", "#b7878d"),
                 c("#c3d79b", "#b28cc0", "#64a985", "#e3a7d4", "#2ea2aa", "#e69892", "#85c6f9", "#fbd1a0", "#7696be", "#89996c", "#ddcdff", "#719d89", "#f5cde6", "#b6e0da", "#e8d4cd", "#b5ddfa"),
                 c("#a98d83", "#84e1ff", "#bb8964", "#46b1d1", "#ffbfa5", "#6199c0", "#bbcb8f", "#bf82ab", "#85ddc4", "#eea0ba", "#c1d8ff", "#c3818b", "#c5c6ff", "#999388", "#e8cbff", "#ffb5b6", "#d2dad7"),
                 c("#faccde", "#60a987", "#c6abe4", "#6f9e77", "#c48093", "#a5e5d3", "#cc8f6f", "#499fae", "#d9dca6", "#7796b8", "#bee1ba", "#b4daff", "#919583", "#e2d3e9", "#47a19b", "#ebd4bc", "#7c9993", "#a9e3e0"),
                 c("#739e6e", "#ffbfd9", "#43b6bb", "#e8ad88", "#5e9bce", "#c2af75", "#a8e0fe", "#fad0a8", "#679e8d", "#ffc7b1", "#abe5c0", "#ac8d78", "#c5dddc", "#a48f84", "#cadfb0", "#899694", "#fdcdc1", "#d1dad5", "#dfd8c4"),
                 c("#6e9c93", "#ffb4b3", "#7ec6a2", "#eeccfe", "#cddb9d", "#8a90c5", "#dcb983", "#77bff0", "#f0ab92", "#90ddff", "#f1d3a9", "#b5c2fe", "#c1e1b7", "#7596ba", "#bce1c4", "#a88c96", "#5a9daf", "#b18b80", "#d4d6f3", "#949577"),
                 c("#e7d6bb", "#749ed5", "#f9d29d", "#67b3e2", "#d09773", "#65ccec", "#d38585", "#7fe8ef", "#cf8190", "#94e8cd", "#ae8cc1", "#b3cf95", "#cbc0fc", "#94a66c", "#eeccff", "#ada368", "#e9a6ce", "#48a297", "#ffc1df", "#799c7a", "#facbe0", "#5d9e9a", "#ffc6c1", "#619bb0", "#fccdcb", "#7197bb", "#b1e4c3", "#9390b1", "#c3e0c0", "#a98c90", "#ade3ce", "#9c927d", "#c2dafe", "#869881", "#e6d3dc", "#6e9ba4", "#bde0d0", "#8196a4", "#b2e1df", "#b9deea")
)

badToGood <- rev(c("#1AB958", "#23B956", "#2BB954", "#31B952", "#37B850", "#3CB84E", "#41B84C", "#45B74A", "#49B749", "#4DB747", "#51B645", "#54B643", "#58B641", "#5BB53F", "#5FB53D", "#62B53C", "#65B43A", "#68B438", "#6BB336", "#6DB335", "#70B333", "#73B231", "#76B230", "#78B12E", "#7BB12C", "#7DB02B", "#80B029", "#82AF28", "#85AF26", "#87AE25", "#8AAE23", "#8CAD22", "#8EAD21", "#91AC1F", "#93AC1E", "#95AB1D", "#97AB1C", "#9AAA1B", "#9CAA1A", "#9EA919", "#A0A918", "#A2A818", "#A4A717", "#A6A716", "#A8A616", "#AAA616", "#ACA515", "#AEA415", "#B0A415", "#B2A315", "#B4A315", "#B6A216", "#B8A116", "#B9A117", "#BBA017", "#BD9F18", "#BF9F18", "#C19E19", "#C29D1A", "#C49D1B", "#C69C1C", "#C79B1D", "#C99A1E", "#CB9A1F", "#CC9920", "#CE9822", "#CF9823", "#D19724", "#D29625", "#D49626", "#D59528", "#D79429", "#D8932A", "#D9932C", "#DB922D", "#DC912E", "#DD9130", "#DF9031", "#E08F33", "#E18F34", "#E28E35", "#E38D37", "#E58C38", "#E68C3A", "#E78B3B", "#E88A3D", "#E98A3E", "#EA8940", "#EB8841", "#EC8843", "#ED8744", "#EE8746", "#EE8647", "#EF8549", "#F0854A", "#F1844C", "#F2844D", "#F2834F", "#F38350", "#F48252", "#F48253", "#F58155", "#F58157", "#F68058", "#F6805A", "#F77F5B", "#F77F5D", "#F87E5E"))

Reference <- function (authors, year, title, journal = '',
                       volume = NULL, pages = NULL, doi = NULL,
                       publisher = NULL, editors = NULL) {
  nAuth <- length(authors)
  if (nAuth > 1L) {
    authors <- paste(paste0(authors[-nAuth], collapse = ', '), "&amp;", authors[nAuth])
  }
  nEd <- length(editors)
  if (nEd > 1L) {
    editors <- paste(paste0(editors[-nEd], collapse = ', '), "&amp;", editors[nEd])
  } else if (nEd < 1) {
    editors <- ''
  }
  paste0("<p class='reference'>", authors, " (", year, "). &ldquo;", title,
         "&rdquo;. ",
         if (editors != '') paste0("In: ", editors, ' (eds). ') else '',
         if (journal != "") paste0("<i>", journal, "</i>. ") else "",
         if (is.null(volume)) "" else paste0("<b>", volume, "</b>:"),
         if (is.null(publisher)) "" else paste0(publisher, '. '),
         if (is.null(pages)) "" else paste0(paste0(pages, collapse = "&ndash;"), '. '),
         if (is.null(doi)) "" else paste0(
           "doi:<a href='https://doi.org/", doi, "' title=\"CrossRef\">",
           doi, "</a>. "), 
         "</p>")
}


Brazeau2019 <- Reference(c('Brazeau, M.D.', 'Guillerme, T.', 'Smith, M.R.'), 2019,
                           title = "An algorithm for morphological phylogenetic analysis with inapplicable data",
                           journal = "Systematic Biology",
                           volume = 64,
                           pages = "619-631",
                         doi = "10.1093/sysbio/syy083")
Bien2011 <- Reference(
  c("Bien, J.", "Tibshirani, R."),
  title = "Hierarchical clustering with prototypes via minimax linkage",
  year = 2011,
  volume = 106,
  doi = "10.1198/jasa.2011.tm10183",
  pages = c(1075, 1084),
  journal = "Journal of the American Statistical Association")
Gower1966 <- Reference(  title = "Some distance properties of latent root and vector methods used in multivariate analysis",
                         author = "Gower, J.C.",
                         year = 1966,
                         volume = 53,
                         pages = c(325, 338),
                         doi = "10.2307/2333639",
                         journal = "Biometrika")
Gower1969 <- Reference(
  title = "Minimum spanning trees and single linkage cluster analysis",
  author = c("Gower, J.C.", "Ross, G.J.S."),
  year = 1969, volume = 18, pages = c(54, 64), doi = "10.2307/2346439",
  journal = "Journal of the Royal Statistical Society. Series C (Applied Statistics)")
Kaski2003 <- Reference(
  title = "Trustworthiness and metrics in visualizing similarity of gene expression",
  author = c("Kaski, S.", "Nikkil&auml;, J.", "Oja, M.", "Venna, J.",
             "T&ouml;r&ouml;nen, P.", "Castr&eacute;n, E."),
  year = 2003, volume = 4, pages = 48, doi = "10.1186/1471-2105-4-48",
  journal = "BMC Bioinformatics")
Maechler2019 <- Reference(
  title = "cluster: cluster analysis basics and extensions", year = 2019,
  author = c("Maechler, M.", "Rousseeuw, P.", "Struyf, A.", "Hubert, M.", "Hornik, K."),
  journal = "Comprehensive R Archive Network")
Morphy <- Reference(c("Brazeau, M.D.", "Smith, M.R.", "Guillerme, T."), 2017,
                    'MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability.',
                    doi = '10.5281/zenodo.815371')
Murtagh1983 <- Reference(
  title = "A survey of recent advances in hierarchical clustering algorithms",
  author = "Murtagh, F.", year = 1983, volume = 26, pages = c(354, 359),
  doi = "10.1093/comjnl/26.4.354", journal = "The Computer Journal")
Nixon1999 <- Reference(
  "Nixon, K.C.", 1999,
  journal = "Cladistics", volume = 15, pages = "407-414",
  title = "The Parsimony Ratchet, a New Method for Rapid Parsimony Analysis",
  doi = "10.1111/j.1096-0031.1999.tb00277.x")
RCoreTeam <- Reference(
  author = "R Core Team", year = 2020,
  title = "R: A language and environment for statistical computing",
  publisher = "R Foundation for Statistical Computing, Vienna, Austria")
SmithDist <- Reference('Smith, M.R.', 2020,
                       'TreeDist: distances between phylogenetic trees',
                       doi = '10.5281/zenodo.3528123', 'Comprehensive R Archive Network')
SmithQuartet <- Reference('Smith, M.R.', 2019,
                          'Quartet: comparison of phylogenetic trees using quartet and split measures',
                          'Comprehensive R Archive Network', doi = "10.5281/zenodo.2536318")
SmithSearch <- Reference('Smith, M.R.', 2018,
                         'TreeSearch: phylogenetic tree search using custom optimality criteria',
                         'Comprehensive R Archive Network', doi = "10.5281/zenodo.1042590")
Smith2020 <- Reference('Smith, M.R.', 2020,
                       'Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees',
                       'Bioinformatics', volume = 36, pages = '5007--5013',
                       doi = "10.1093/bioinformatics/btaa614")
Smith2021 <- Reference('Smith, M.R.', 2021,
                       'The importance of methodology when analyzing landscapes of phylogenetic trees',
                       'Submitted MS')
Venna2001 <- Reference(
  title = "Neighborhood preservation in nonlinear projection methods: an experimental study",
  author = c("Venna, J.", "Kaski, S."), year = 2001, pages = c(485, 491),
  journal = "Lecture Notes in Computer Science: Artificial Neural Networks&mdash;ICANN 2001",
  editors = c("Dorffner, G.", "Bischof, H.", "Hornik, K."),
  publisher = "Springer, Berlin",
  doi = "10.1007/3-540-44668-0_68")

ui <- fluidPage(theme = 'app.css',
  useShinyjs(),
  column(3,
    fluidRow(
      tags$h1("TreeSearch beta UI"),
      fileInput("datafile", "Load data", placeholder = "No data file selected"),
      radioButtons("dataFormat", "Data format", 
                   list("Nexus" = 'nex', "TNT" = 'tnt'),
                   'nex', TRUE),
      selectInput('character.weight', "Character weighting", list("Equal" = "equal"), "equal"),
      selectInput('implied.weights', "Step weighting", 
                  list("Implied" = "on", "Equal" = "off"), "on"),
      sliderInput("concavity", "Step weight concavity constant", min = 0L,
                  max = 3L, pre = '10^', value = 1L),
      sliderInput('ratchIter', "Ratchet iterations", min = 0L, max = 50L, value = 6L, step = 1L),
      # sliderInput('ratchIter', "Ratchet iterations", min = 0L, max = 50L, value = 0L, step = 1L),
      sliderInput('tbrIter', "TBR iterations", min = 1L, max = 50L, value = 6L, step = 1L),
      # sliderInput('tbrIter', "TBR iterations", min = 1L, max = 50L, value = 1L, step = 1L),
      sliderInput('maxHits', "Maximum hits", min = 0L, max = 5L, value = 2L, pre = '10^'),
      # sliderInput('maxHits', "Maximum hits", min = 0L, max = 5L, value = 0.6, pre = '10^'),
      sliderInput('finalIter', "Final iteration extra depth", min = 1L, max = 10L, value = 3L),
      # sliderInput('finalIter', "Final iteration extra depth", min = 1L, max = 10L, value = 1L),
      sliderInput('verbosity', "Notification level", min = 0L, max = 3L, value = 3L, step = 1L),
      actionButton("go", "Search"),
      textOutput("results"),
      hidden(radioButtons('plotFormat', "Display:",
                   list("Consensus tree" = 'cons',
                        "Cluster consensus trees" = 'clus',
                        "Individual tree" = 'ind',
                        "Tree space" = 'space'), 'cons')),
      hidden(sliderInput('whichTree', 'Tree to plot', min = 1L, max = 1L,
                         value = 1L, step = 1L))
      ),
  ),
  column(9,
    fluidRow(id = 'plotConfig',
      tags$span("Plot size:", id = 'plotSizeSpan'),
      sliderInput(inputId = "plotSize", label = NULL,
                                            width = '200px',
                                            min = 100, max = 2000,
                                            post = 'px', value = 600),
      tags$span("Save as: "),
      downloadButton('savePdf', 'PDF'),
      downloadButton('savePng', 'PNG'),
      downloadButton('saveNwk', 'Newick'),
      downloadButton('saveNex', 'Nexus'),
    ),
    fluidRow(
      plotOutput(outputId = "treePlot", height = "600px"),
      hidden(plotOutput('clustCons', height = "200px")),
      htmlOutput('references'),
    ),
  )
)

server <- function(input, output, session) {
  
  r <- reactiveValues()
  
  ##############################################################################
  # Load data
  ##############################################################################
  dataset <- reactive({
    fileInput <- input$datafile
    if (is.null(input)) return("No data file selected.")
    tmpFile <- fileInput$datapath
    if (is.null(tmpFile)) return ("No data file found.")
    if (input$dataFormat == 'nex') {
      ret <- ReadAsPhyDat(tmpFile)
    } else {
      ret <- ReadTntAsPhyDat(tmpFile)
    }
    if (!is.null(r$trees) && 
        length(intersect(names(ret), r$trees[[1]]$tip.label)) != 
        length(ret)) {
      r$trees <- NULL
      updateActionButton(session, "go", "New search")
    }
    ret
  })
  
  concavity <- reactive({
    if (input$implied.weights == 'on') {
      10 ^ input$concavity
    } else Inf
  })
  
  observeEvent(input$go, {
    if (!inherits(dataset(), 'phyDat')) {
      showNotification("No data loaded", type = 'error')
    } else {
      r$trees <- withProgress(c(MaximizeParsimony(dataset(),
                                   tree = if(is.null(r$trees)) 
                                     TreeTools::NJTree(dataset())
                                   else
                                     r$trees[[1]],
                                   concavity = concavity(),
                                   ratchIter = input$ratchIter,
                                   tbrIter = input$tbrIter,
                                   maxHits = ceiling(10 ^ input$maxHits),
                                   finalIter = input$finalIter,
                                   verbosity = input$verbosity,
                                   session = session)),
                              message = "Searching...")
      updateSliderInput(session, 'whichTree', min = 1L,
                        max = length(r$trees), value = 1L)
      output$results <- renderText(paste0(
        "Found ", length(r$trees), " trees with score ", 
        signif(TreeLength(r$trees[[1]], dataset(), concavity = concavity()))))
      updateActionButton(session, "go", "Continue search")
      show('plotFormat')
      PlotConsensus()
    }
  })
  
  observeEvent(input$plotFormat, PlotConsensus())
  
  RenderMainPlot <- function (x) {
    renderPlot(x, width = PlotSize(), height = PlotSize())
  }
  
  PlotConsensus <- function () (
    if (!is.null(r$trees) && inherits(r$trees, 'multiPhylo'))
    switch(input$plotFormat,
           'cons' = {
             hide('whichTree')
             output$treePlot = RenderMainPlot({
               par(mar = rep(0, 4), cex = 0.9)
               plot(consensus(r$trees))
             })
           },
           'clus' = {
             hide('whichTree')
             output$treePlot = RenderMainPlot(PlotClusterCons())
           },
           'ind' = {
             show('whichTree')
             output$treePlot = RenderMainPlot({
               par(mar = rep(0, 4), cex = 0.9)
               plot(r$trees[[input$whichTree]])
             })
           },
           'space' = {
             hide('whichTree')
             output$treePlot = RenderMainPlot(treespacePlot())
           }
    )
  )
  
  output$pcQuality <- renderPlot({
    par(mar = c(2, 0, 0, 0))
    nStop <- length(badToGood)
    
    plot(seq(0, 1, length.out = nStop), rep(0, nStop),
         pch = 15, col = badToGood,
         ylim = c(-1.5, 2.5),
         ann = FALSE, axes = FALSE)
    
    logScore <- LogScore(projQual()['TxC', dims()])
    lines(rep(logScore, 2), c(-1, 1), lty = 3)
    
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    ticks <- LogScore(tickPos)
    
    axis(1, at = ticks, labels = NA, line = 0)
    axis(1, tick = FALSE, at = ticks, labels = tickPos, line = 0)
    axis(1, line = -1, tick = FALSE,
         at = ticks[-1] - ((ticks[-1] - ticks[-length(ticks)]) / 2),
         labels = c('', 'dire', '', "ok", "gd", "excellent"))
    axis(3, at = 0.5, tick = FALSE, line = -2, 
         paste0(dims(), 'D projection quality (trustw. \ud7 contin.):'))
  })
  
  
  output$howManyDims <- renderPlot({
    par(mar = c(2.5, 2.5, 0, 0), xpd = NA, mgp = c(1.5, 0.5, 0))
    txc <- projQual()['TxC', ]
    nStop <- length(badToGood)
    
    plot(txc, type = 'n', ylim = c(min(txc, 0.5), 1),
         frame.plot = FALSE, axes = FALSE,
         xlab = 'Dimension', ylab = 'Trustw. \uD7 Contin.')
    par(xpd = FALSE)
    axis(1, 1:14)
    axis(2)
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    mids <- c(0.6, 0.75, 0.85, 0.925)
    text(rep(15, 4), mids, pos = 2, cex = 0.8,
         col = badToGood[nStop * LogScore(mids)],
         c('Essentially random', 'Dangerous', "Usable", "Good"))
    text(1, 0.975, pos = 4, "Excellent", cex = 0.8, 
         col = badToGood[LogScore(0.975) * nStop])
    for (i in tickPos[-1]) {
      abline(h = i, lty = 3, col = badToGood[LogScore(i) * nStop])
    }
    points(txc, type = 'b')
    txcNow <- txc[dims()]
    
    points(dims(), txcNow, pch = 16, col = badToGood[LogScore(txcNow) * nStop],
           cex = 1.6)
  })
  
  ##############################################################################
  # Clusterings
  ##############################################################################
  clusterings <- reactive({
    maxCluster <- min(15L, length(r$trees) - 1L)
    if (maxCluster > 1L) {
      possibleClusters <- 2:maxCluster
      
      hSil <- pamSil <- -99
      dists <- distances()
      
      nMethodsChecked <- 2
      methInc <- 1 / nMethodsChecked
      nK <- length(possibleClusters)
      kInc <- 1 / (nMethodsChecked * nK)
    
      pamClusters <- lapply(possibleClusters, function (k) {
        cluster::pam(dists, k = k)
      })
      pamSils <- vapply(pamClusters, function (pamCluster) {
        mean(cluster::silhouette(pamCluster)[, 3])
      }, double(1))
      bestPam <- which.max(pamSils)
      pamSil <- pamSils[bestPam]
      pamCluster <- pamClusters[[bestPam]]$cluster
      
      hTree <- protoclust::protoclust(dists)
      hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
      hSils <- vapply(hClusters, function (hCluster) {
        mean(cluster::silhouette(hCluster, dists)[, 3])
      }, double(1))
      bestH <- which.max(hSils)
      hSil <- hSils[bestH]
      hCluster <- hClusters[[bestH]]
    
      bestCluster <- c('none', 'pam', 'hmm')[which.max(c(0, pamSil, hSil))]
    } else {
      bestCluster <- 'none'
    }
      
    # Return:
    list(method = switch(bestCluster, pam = 'part. around medoids',
                                      hmm = 'minimax linkage',
                                      none = 'no significant clustering'),
         n = 1 + switch(bestCluster, pam = bestPam, hmm = bestH, 0),
         sil = switch(bestCluster, pam = pamSil, hmm = hSil, 0), 
         cluster = switch(bestCluster, pam = pamCluster, hmm = hCluster, 1)
    )

  })
  
  consRows <- reactive({
    cl <- clusterings()
    if (cl$sil > 0.25) ceiling(cl$n / 3) else 1L
  })
  
  consSize <- reactive({
    nLeaves() * 12 * consRows()
  })
  
  PlotClusterCons <- function () {
    cl <- clusterings()
    par(mar = c(0.2, 0, 0.2, 0), xpd = NA)
    par(cex = 0.7)
    if (cl$sil > 0.25) {
      par(mfrow = c(consRows(), ceiling(cl$n / consRows())))
      for (i in seq_len(cl$n)) {
        col <- palettes[[min(length(palettes), cl$n)]][i]
        tr <- ape::consensus(r$trees[cl$cluster == i])
        tr$edge.length <- rep(1, dim(tr$edge)[1])
        plot(tr, edge.width = 2, font = 1, cex = 1,
             edge.color = col, tip.color = col)
      }
    } else {
      tr <- ape::consensus(r$trees)
      tr$edge.length <- rep(1, dim(tr$edge)[1])
      plot(tr,edge.width = 2, font = 1, cex = 1,
           edge.color = palettes[[1]],
           tip.color = palettes[[1]])
    }
  }
  
  ##############################################################################
  # Plot settings: point style
  ##############################################################################

  pointCols <- reactive({
    cl <- clusterings()
    if (cl$sil > 0.25) {
      palettes[[min(length(palettes), cl$n)]][cl$cluster]
    } else palettes[[1]]
  })
  
  pchs <- reactive({
    cl <- clusterings()
    if (cl$sil > 0.25) {
      cl$cluster - 1
    } else 16
  })
  
  maxProjDim <- reactive({
    min(6L, length(r$trees) - 1L)
  })
  
  nProjDim <- reactive({
    dim(projection())[2]
  })
  
  dims <- debounce(reactive({
    if (mode3D()) 3L else {
      min(6, maxProjDim())
    }
  }), 400)
  
  distances <- reactive({
    if (length(r$trees) > 1L) {
      withProgress(
        message = 'Calculating distances', value = 0.99,
        TreeDist::ClusteringInfoDistance(r$trees)
      )
    } else {
      matrix(0, 0, 0)
    }
    
  })
  
  projection <- reactive({
    if (maxProjDim() > 1L) {
      withProgress(
        message = 'Projecting distances',
        value = 0.99,
        proj <- cmdscale(distances(), k = maxProjDim())
      )
      # Return:
      proj
      
    } else {
      matrix(0, 0, 0)
    }
  })
  
  mstEnds <- reactive({
    dist <- as.matrix(distances())
    nRows <- dim(dist)[1]
    withProgress(message = 'Calculating MST', {
      edges <- MSTEdges(dist)
    })
    edges
  })
  
  ##############################################################################
  # Plot tree space
  ##############################################################################
  treespacePlot <- function() {
    cl <- clusterings()
    proj <- projection()
    
    nDim <- min(dims(), nProjDim())
    plotSeq <- matrix(0, nDim, nDim)
    nPanels <- nDim * (nDim - 1L) / 2L
    plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)
    layout(t(plotSeq[-nDim, -1]))
    par(mar = rep(0.2, 4))
    withProgress(message = 'Drawing plot', {
      for (i in 2:nDim) for (j in seq_len(i - 1)) {
        incProgress(1 / nPanels)
        # Set up blank plot
        plot(proj[, j], proj[, i], ann = FALSE, axes = FALSE,
             frame.plot = nDim > 2L,
             type = 'n', asp = 1, xlim = range(proj), ylim = range(proj))
        
        # Plot MST
        apply(mstEnds(), 1, function (segment)
          lines(proj[segment, j], proj[segment, i], col = "#bbbbbb", lty = 1))
        
        # Add points
        points(proj[, j], proj[, i], pch = pchs(),
               col = paste0(pointCols(), as.hexmode(200)),
               cex = input$pt.cex)
        
        if (cl$sil > 0.25) {
          # Mark clusters
          for (clI in seq_len(cl$n)) {
            inCluster <- cl$cluster == clI
            clusterX <- proj[inCluster, j]
            clusterY <- proj[inCluster, i]
            hull <- chull(clusterX, clusterY)
            polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
                    border = palettes[[min(length(palettes), cl$n)]][clI])
            #border = '#54de25bb')
          }
        }
        if ('labelTrees' %in% input$display) {
          text(proj[, j], proj[, i], thinnedTrees())
        }
      }
    })
  }
  
  mode3D <- reactive("show3d" %in% input$display)
  PlotSize <- function () debounce(reactive(input$plotSize), 100)
  output$distPlot <- RenderMainPlot({
    if (!mode3D()) {
      if (inherits(distances(), 'dist')) {
        treespacePlot()
        output$projectionStatus <- renderText('')
      } else {
        output$projectionStatus <- renderText("No distances available.")
      }
    }
  })
  
  plotContent <- reactive({
    switch(input$plotFormat,
           'cons' = {
             par(mar = rep(0, 4), cex = 0.9)
             plot(consensus(r$trees))
           },
           'clus' = {
             PlotClusterCons()
           },
           'ind' = {
             par(mar = rep(0, 4), cex = 0.9)
             plot(r$trees[[input$whichTree]])
           },
           'space' = {
             treespacePlot()
           })
  })
  
  output$savePng <- downloadHandler(
    filename = 'TreeSearch.png',
    content = function (file) {
      png(file, width = input$plotSize, height = input$plotSize)
      plotContent()
      dev.off()
    })
  
  output$savePdf <- downloadHandler(
    filename = 'TreeSearch.pdf',
    content = function (file) {
      pdf(file, title = paste0('Tree space projection'),
          width = 10,
          height = 20)
      plotContent()
      dev.off()
    })
  
  output$saveNwk <- downloadHandler(
    filename = 'TreeSearch.nwk',
    content = function(file) {
      write.tree(r$trees, file = file)
    }
  )
  
  output$saveNex <- downloadHandler(
    filename = 'TreeSearch.nex',
    content = function(file) {
      write.nexus(r$trees, file = file)
    }
  )
  
  
  ##############################################################################
  # References
  ##############################################################################
  output$references <- renderUI({
    tags$div(style = paste0('position: relative; top: ', 
                            (input$plotSize - 600)
  #                          + if ('cons' %in% input$display) consSize() - 200 else 0
                            , 'px'),
             tagList(
               tags$h2('References for methods used'),
               tags$h3('Tree search'),
               HTML(Brazeau2019, Morphy, Nixon1999, SmithSearch),
               tags$h3('Tree space projection'),
               HTML(paste0(Gower1966, Gower1969, Kaski2003, RCoreTeam,
                           SmithDist, Smith2020, Smith2021, 
                           Venna2001)),
               ),
               tags$h3('Clustering'),
               HTML(paste0('Partitioning around medoids:', Maechler2019,
                           "Hierarchical, minimax linkage:", Bien2011,
                           Murtagh1983))
             )
  })
}

shinyApp(ui = ui, server = server)