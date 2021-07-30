library("methods", exclude = c('show', 'removeClass'))
library("cli")
suppressPackageStartupMessages({
  library("shiny", exclude = c('runExample'))
  library("shinyjs", exclude = c('runExample'))
})
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeSearch")
library("TreeDist")

if (!requireNamespace('TreeDist', quietly = TRUE)) {
  install.packages('TreeDist')
  library('TreeDist')
}

palettes <- list("#7a6c36",
                 c("#7a6c36", "#864885"),
                 c("#7a6c36", "#864885", "#427743"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020", "#c241a7"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020", "#c241a7", "#391d42"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#85c6f9", "#fbd1a0", "#7696be", "#89996c", "#ddcdff", "#719d89", "#f5cde6", "#b6e0da", "#e8d4cd", "#b5ddfa"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#bbcb8f", "#bf82ab", "#85ddc4", "#eea0ba", "#c1d8ff", "#c3818b", "#c5c6ff", "#999388", "#e8cbff", "#ffb5b6", "#d2dad7"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#cc8f6f", "#499fae", "#d9dca6", "#7796b8", "#bee1ba", "#b4daff", "#919583", "#e2d3e9", "#47a19b", "#ebd4bc", "#7c9993", "#a9e3e0"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#a8e0fe", "#fad0a8", "#679e8d", "#ffc7b1", "#abe5c0", "#ac8d78", "#c5dddc", "#a48f84", "#cadfb0", "#899694", "#fdcdc1", "#d1dad5", "#dfd8c4"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#dcb983", "#77bff0", "#f0ab92", "#90ddff", "#f1d3a9", "#b5c2fe", "#c1e1b7", "#7596ba", "#bce1c4", "#a88c96", "#5a9daf", "#b18b80", "#d4d6f3", "#949577"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#e03795", "#438f2e", "#5e2195", "#758029", "#4042b9", "#a37926", "#8364df", "#c3671f", "#444491", "#dc4c1f", "#367076", "#e2383c", "#4786b4", "#e13964", "#4c8c73", "#a53396", "#2c4422", "#b553cb", "#50381b", "#4f75d8", "#a12c1b", "#8576b8", "#bd6541", "#3a1959", "#83491f", "#2d2644", "#c45b94", "#451523", "#966883", "#782224", "#b96563", "#762254", "#95765c", "#ad355a")
)
silThreshold <- 0.4
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
Smith2021 <- Reference('Smith, M.R.', 2022,
                       'Robust analysis of phylogenetic tree space',
                       'Submitted MS')
Smith2022 <- Reference("Smith, M.R.", 2022, 
                       "Improving consensus trees by detecting rogue taxa",
                       "Submitted MS")
Stockham2002 <- Reference(
  author = c('Stockham, C.', 'Wang, L.-S.', 'Warnow, T.'), 2002,
  "Statistically based postprocessing of phylogenetic analysis by clustering",
  "Bioinformatics", 18, c('S285', 'S293'),
  doi = "10.1093/bioinformatics/18.suppl_1.S285")

Venna2001 <- Reference(
  title = "Neighborhood preservation in nonlinear projection methods: an experimental study",
  author = c("Venna, J.", "Kaski, S."), year = 2001, pages = c(485, 491),
  journal = "Lecture Notes in Computer Science: Artificial Neural Networks&mdash;ICANN 2001",
  editors = c("Dorffner, G.", "Bischof, H.", "Hornik, K."),
  publisher = "Springer, Berlin",
  doi = "10.1007/3-540-44668-0_68")

ui <- fluidPage(theme = 'app.css',
                title = 'TreeSearch',
  useShinyjs(),
  column(3,
    fluidRow(
      tags$h1("TreeSearch beta UI"),
      fileInput("datafile", "Load data",
                placeholder = "No data file selected"),
      tags$label("Search", class = "control-label", 
                 style = 'display: block; margin-top: -15px;'),
      actionButton("searchConfig", "Configure", icon = icon('cogs')),
      hidden(actionButton("go", "Search", icon = icon('search'))),
      fileInput("treeFile", label = "Load trees",
                placeholder = "No tree file selected"),
      textOutput("results"),
      hidden(radioButtons('plotFormat', "Display:",
                   list("Characters on trees" = 'ind',
                        "Consensus tree" = 'cons',
                        "Cluster consensus trees" = 'clus',
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
      htmlOutput(outputId = "plotSpacer", height = "0px"),
      hidden(plotOutput('clustCons', height = "200px")),
      hidden(tags$div(id = 'charChooser',
        tags$div(numericInput('whichChar', 'Character to map:', value = 1L,
                              min = 0L, max = 1L, step = 1L, width = 200),
                 checkboxGroupInput('mapDisplay', '', list(
                   "Align tips" = "tipsRight",
                   "Infer tips" = "updateTips"
                   )),
                 style = "float: right; width: 200px; margin-left: 2em;"),
        htmlOutput('charMapLegend'),
        htmlOutput('charNotes'),
      )),
      hidden(tags$div(id = 'consConfig',
        tags$div(style = "float: right; width: 200px; margin-left: 2em;",
          sliderInput('consP', 'Majority:', value = 1,
                      min = 0.5, max = 1, width = 200),
          sliderInput('keepTips', 'Tips to show:', value = 1L,
                      min = 1L, max = 1L, step = 1L, width = 200),
          selectizeInput('neverDrop', 'Never drop:', multiple = TRUE,
                         choices = list("NA" = "No trees loaded"))
                 ),
        htmlOutput('droppedTips')
      )),
      htmlOutput('references', style = "clear: both;"),
    ),
  )
)

server <- function(input, output, session) {
  
  r <- reactiveValues()
  
  ##############################################################################
  # Load data
  ##############################################################################
  
  datasetMatchesTrees <- reactive({
    length(intersect(names(r$dataset), r$trees[[1]]$tip.label)) == 
      length(r$dataset)
  })
  
  observeEvent(input$datafile, {
    fileInput <- input$datafile
    r$dataset <- NULL
    r$chars <- NULL
    if (is.null(input)) {
      showNotification(type = "error", "No data file selected")
      return("No data file selected.")
    }
    tmpFile <- fileInput$datapath
    if (is.null(tmpFile)) {
     showNotification(type = "error", "No data file found.")
     return ("No data file found.")
    }
    r$dataset <- tryCatch(ReadTntAsPhyDat(tmpFile),
                          error = function (e) tryCatch({
                            r$chars <- ReadCharacters(tmpFile)
                            r$charNotes <- ReadNotes(tmpFile)
                            ReadAsPhyDat(tmpFile)
                            }, error = function (e) NULL))
    
    if (is.null(r$dataset)) {
     showNotification(type = "error", "Could not read data from file")
      
      updateSliderInput(session, 'whichChar', min = 0L,
                        max = 0L, value = 0L)
     return ("Could not read data from file")
    } else {
      showNotification(type = "message", 
                       paste("Loaded", attr(r$dataset, 'nr'), "characters and",
                             length(r$dataset), "taxa"))
      
      updateSliderInput(session, 'whichChar', min = 0L,
                        max = as.integer(attr(r$dataset, 'nr')), value = 1L)
    }
    
    
    if (!is.null(r$trees)) {
      if (!datasetMatchesTrees()) {
        r$trees <- NULL
        updateActionButton(session, "go", "New search")
      } else {
        show('plotFormat')
      }
      
      output$results <- TreeSummary()
    }
    
  })
  
  observeEvent(input$searchConfig, {
    updateSelectInput(session, 'character.weight',
                      selected = input$character.weight)
    updateSelectInput(session, 'implied.weights',
                      selected = input$implied.weights)
    updateSliderInput(session, 'concavity', value = input$concavity)
    updateSliderInput(session, 'ratchIter', value = input$ratchIter)
    updateSliderInput(session, 'tbrIter', value = input$tbrIter)
    updateSliderInput(session, 'maxHits', value = input$maxHits)
    updateSliderInput(session, 'startIter', value = input$startIter)
    updateSliderInput(session, 'finalIter', value = input$finalIter)
    showModal(modalDialog(
      easyClose = TRUE,
      fluidPage(column(6,
      tagList(selectInput('character.weight', "Character weighting",
                          list("Equal" = "equal"), "equal"),
              selectInput('implied.weights', "Step weighting", 
                         list("Implied" = "on", "Profile" = "prof",
                              "Equal" = "off"), "on"),
              sliderInput("concavity", "Step weight concavity constant", min = 0L,
                         max = 3L, pre = '10^', value = 1L),
              sliderInput('ratchIter', "Ratchet iterations", min = 0L,
                          max = 50L, value = 6L, step = 1L),
              sliderInput('maxHits', "Maximum hits", min = 0L, max = 5L,
                          value = 2L, pre = '10^'),
      )), column(6, 
             tagList(
              sliderInput('tbrIter', "TBR depth", min = 1L, max = 20L,
                          value = 1L, step = 1L),
              sliderInput('startIter', "First iteration extra depth", min = 1L,
                          max = 10L, value = 3L),
              sliderInput('finalIter', "Final iteration extra depth", min = 1L,
                          max = 10L, value = 1L),
             ))
      ),
      title = "Tree search settings",
      footer = tagList(modalButton('Close', icon = icon('window-close')),
                       actionButton("modalGo", icon = icon('search'), 
                                    if(length(r$trees)) "Continue search" else "Start search"))
    ))
    show('go')
  })
  
  observeEvent(input$treeFile, {
    tmpFile <- input$treeFile$datapath
    trees <- tryCatch(read.tree(tmpFile),
                      error = function (x) tryCatch(read.nexus(tmpFile),
                                                    error = function (x) NULL))
    if (is.null(trees)) {
      showNotification("Trees not in a recognized format", type = 'error')
    } else {
      r$trees <- c(trees)
      removeModal()
      showNotification(paste("Loaded", length(r$trees), "trees"), type = "message")
      updateSliderInput(session, 'whichTree', min = 1L,
                        max = length(r$trees), value = 1L)
      updateActionButton(session, "modalGo", "Continue search")
      updateActionButton(session, "go", "Continue")
      show('plotFormat')
      
      output$results <- TreeSummary()
    }
    
  })
  
  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           'on' = show('concavity'),
           hide('concavity')
    )
  })
  
  weighting <- reactive(
    if (length(input$implied.weights) > 0) {
      input$implied.weights
    } else {
      'on'
    }
  )
  
  TreeSummary <- debounce(reactive({
    scores <- tryCatch(signif(TreeLength(r$trees, r$dataset,
                                         concavity = concavity())),
                       error = function (x) {
                         if (length(r$dataset) > 0 && 
                             length(r$trees) > 0) {
                           cli::cli_alert(x[[2]])
                           cli::cli_alert_danger(x[[1]])
                           showNotification(type = "error",
                                            "Could not score all trees with dataset")
                         }
                         NULL
                       })
    
    wtType <- switch(weighting(),
                     'on' = paste0('k = ', signif(concavity(), 3)),
                     'off' = 'EW',
                     'prof' = 'PP')
    score <- if (is.null(scores)) {
      "; could not be scored from dataset"
    } else if (length(unique(scores)) == 1) {
      paste0(", each with score ", scores[1], " (", wtType, ")")
    } else {
      paste0(" with scores ", min(scores), " to ", max(scores),
             " (", wtType, ")")
    }
    
    
    renderText(paste0(length(r$trees), " trees in memory", score))
  }), 50)
  
  observeEvent(input$implied.weights, {
    output$results <- TreeSummary()
  })
  
  observeEvent(input$concavity, {
    output$results <- TreeSummary()
  })
  
  observeEvent(r$trees, {
    tipLabels <- r$trees[[1]]$tip.label
    nTip <- length(tipLabels)
    r$trees[-1] <- RenumberTips(r$trees[-1], tipLabels)
    updateSliderInput(inputId = 'keepTips', max = nTip, value = nTip)
    updateSelectizeInput(inputId = 'neverDrop', choices = tipLabels)
  })

  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           'on' = show('concavity'),
           hide('concavity')
    )
  })

  concavity <- reactive({
    kExp <- if (length(input$concavity)) input$concavity else 1
    switch(weighting(),
           'on' = 10 ^ kExp,
           'off' = Inf,
           'prof' = 'Profile')
  })
  
  StartSearch <- function () {
    if (!inherits(r$dataset, 'phyDat')) {
      showNotification("No data loaded", type = 'error')
    } else {
      r$trees <- c(MaximizeParsimony(r$dataset,
                                     tree = if(is.null(r$trees)) 
                                       TreeTools::NJTree(r$dataset)
                                     else
                                       r$trees[[1]],
                                     
                                     concavity = concavity(),
                                     ratchIter = input$ratchIter,
                                     tbrIter = input$tbrIter,
                                     maxHits = ceiling(10 ^ input$maxHits),
                                     startIter = input$startIter,
                                     finalIter = input$finalIter,
                                     verbosity = 4L))
      
      updateSliderInput(session, 'whichTree', min = 1L,
                        max = length(r$trees), value = 1L)
      output$results <- renderText(paste0(
        "Found ", length(r$trees), " trees with score ", 
        signif(TreeLength(r$trees[[1]], r$dataset, concavity = concavity()))))
      updateActionButton(session, "go", "Continue")
      updateActionButton(session, "modalGo", "Continue search")
      show('plotFormat')
    }
  }
  
  observeEvent(input$go, StartSearch())
  observeEvent(input$modalGo, StartSearch())
  
  PlottedChar <- debounce(reactive(as.integer(input$whichChar)), 50)
  PlottedTree <- debounce(reactive(Postorder(r$trees[[input$whichTree]])), 100)
  
  RenderMainPlot <- function (x) {
    renderPlot(x, width = PlotSize(), height = PlotSize())
  }
  
  Instab <- reactive({
    TipInstability(r$trees)
  })
  
  dropSeq <- reactive(rogues()$taxon[-1])
  
  stableCol <- reactive({
    ColByStability(r$trees)
  })
  
  
  
  rogues <- reactive(
    withProgress(
      message = 'Identifying rogues', value = 0.99,
      Rogue::QuickRogue(r$trees, neverDrop = input$neverDrop,
                        fullSeq = TRUE)
    )
  )
  
  tipCols <- reactive(stableCol()) # TODO allow user to choose how to colour
  
  consP <- debounce(reactive(input$consP), 50)
  
  ConsensusPlot <- function() {
    par(mar = rep(0, 4), cex = 0.9)
    #instab <- Instab()
    #dropped <- names(instab[order(instab) > input$keepTips])
    kept <- rev(dropSeq())[seq_len(input$keepTips)]
    dropped <- setdiff(TipLabels(r$trees[[1]]), kept)
    if (length(dropped)) {
      output$droppedTips <- renderUI({tagList(
        tags$h3("Tips excluded:"),
        tags$ul(lapply(dropped, function (i) {
          tags$li(i)
          #tags$li(paste0(i, ' (', signif(instab[i]), ')'))
        }))
        )})
    } else {
      output$droppedTips <- renderUI({})
    }
    plot(ConsensusWithout(r$trees, dropped, p = consP()),
         tip.color = tipCols()[kept])
  }
  
  ShowConfigs <- function (visible = character(0)) {
    allConfigs <- c('whichTree', 'charChooser', 'consConfig')
    lapply(visible, show)
    lapply(setdiff(allConfigs, visible), hide)
  }
  
  MainPlot <- reactive({
    if (!is.null(r$trees)) {
      
    ShowConfigs(switch(input$plotFormat,
                       'cons' = 'consConfig',
                       'ind' = c('whichTree', 'charChooser'),
                       ''))
    switch(input$plotFormat,
           'cons' = {
             ConsensusPlot()
           },
           'clus' = {
             PlotClusterCons()
           },
           'ind' = {
           
             par(mar = rep(0, 4), cex = 0.9)
             n <- PlottedChar()
             if (length(n) && n > 0L) {
               UnitEdge <- function (tr) {
                 tr$edge.length <- rep_len(2, dim(tr$edge)[1])
                 tr
               }
               
               
               treeToPlot <- if('tipsRight' %in% input$mapDisplay) {
                 PlottedTree()
               } else {
                 UnitEdge(PlottedTree())
               }
               
               tryCatch({
                 pc <- PlotCharacter(treeToPlot, r$dataset, n,
                                     edge.width = 2.5,
                                     updateTips = 'updateTips' %in% input$mapDisplay)
               }, error = function (cond) {
                 cli::cli_alert_danger(cond)
                 showNotification(type = 'error',
                                  "Could not match dataset to taxa in trees")
                 plot(0, 0, type = 'n', axes = FALSE, ann = FALSE)
                 text(0, 0, 'Load dataset with\ncharacter codings\nfor taxa on tree',
                      col = '#dd6611', font = 2)
                 return()
               })
                 
               
               if (!is.null(r$chars)) {
                 output$charMapLegend <- renderUI({
                   pal <- c("#00bfc6", "#ffd46f", "#ffbcc5", "#c8a500",
                            "#ffcaf5", "#d5fb8d", "#e082b4", "#25ffd3",
                            "#a6aaff", "#e6f3cc", "#67c4ff", "#9ba75c",
                            "#60b17f")
                   
                     UCFirst <- function (str) {
                       paste0(toupper(substr(str, 1, 1)),
                              substr(str, 2, nchar(str)))
                     }
                     states <- attr(r$chars, 'state.labels')[[n]]
                     tokens <- colnames(pc)
                     appTokens <- setdiff(tokens, '-')
                     .State <- function (glyph, text = 'Error?', col = 'red') {
                       if (is.numeric(glyph)) {
                         if (glyph > length(appTokens)) return (NULL)
                         nonBlank <- states != ''
                         text <- states[nonBlank][glyph]
                         col <- pal[glyph]
                         glyph <- appTokens[glyph]
                       }
                       tags$li(style = 'margin-bottom: 2px;',
                         tags$span(glyph,
                                   style = paste("display: inline-block;",
                                                 "border: 1px solid;",
                                                 "width: 1em;",
                                                 "text-align: center;",
                                                 "line-height: 1em;",
                                                 "margin-right: 0.5em;",
                                                 "background-color:", col, ";")),
                         tags$span(UCFirst(text)))
                     }
                     tagList(
                       tags$h3(colnames(r$chars)[n]),
                       tags$ul(style = "list-style: none;",
                         .State(1), .State(2), .State(3), .State(4), .State(5),
                         .State(6), .State(7), .State(8), .State(9),
                         .State(10), .State(11), .State(12), .State(13),
                         if ('-' %in% tokens) 
                           .State("-", "Inapplicable", 'lightgrey'),
                         .State("?", "Ambiguous", 'grey')
                       )
                     )
                 })
               }
               if (!is.null(r$charNotes)) {
                 
                 output$charNotes <- renderUI({
                   charNotes <- r$charNotes[[n]]
                   description <- charNotes[[1]]
                   notes <- charNotes[[2]]
                   states <- attr(r$chars, 'state.labels')[[n]]
                   tokens <- colnames(pc)
                   
                   tagList(
                     if (length(description) > 0) {
                       tags$div(id = 'char-description',
                                lapply(strsplit(description, '\n')[[1]], tags$p))
                     },
                     tags$ul(class = 'state-notes', {
                       PrintNote <- function(note) {
                         taxa <- names(note)[note]
                         tags$li(class = 'state-note',
                                 tags$span(class = 'state-note-label',
                                           paste(gsub('_', ' ', fixed = TRUE,
                                                      taxa), collapse = ', ')),
                                 tags$span(class = 'state-note-detail',
                                           notes[taxa[1]]))
                       }
                       
                       DuplicateOf <- function(x) {
                         duplicates <- duplicated(x)
                         masters <- x[!duplicates]
                         vapply(masters, function(d) x == d, logical(length(x)))
                       }
                       if (length(notes) == 1) {
                         onlyOne <- TRUE
                         names(onlyOne) <- names(notes)
                         PrintNote(onlyOne)
                       }
                       else {
                         notes <- notes[order(names(notes))]
                         duplicates <- DuplicateOf(toupper(notes))
                         apply(duplicates, 2, PrintNote)
                       }
                     }),
                     if (!states[[1]] %in% c("", "''")
                         && any(tokens == '-')) {
                       tags$p("Brazeau _et al._ (2019) advise that neomorphic characters should not contain inapplicable tokens (-).")
                     }
                   )
                })
               }
             } else {
               output$charMapLegend <- renderUI({})
               output$charNotes <- renderUI({})
               plot(PlottedTree(), tip.color = tipCols())
             }
           },
           'space' = {
             treespacePlot()
           }
    )}
  })
  
  PlotSize <- function () debounce(reactive(input$plotSize), 100)
  
  output$treePlot <- renderPlot(MainPlot(),
                                width = PlotSize(),
                                height = PlotSize())
  
  output$pcQuality <- renderPlot({
    par(mar = c(2, 0, 0, 0))
    nStop <- length(badToGood)
    
    plot(seq.int(from = 0, to = 1, length.out = nStop), numeric(nStop),
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
    text(rep.int(15, 4), mids, pos = 2, cex = 0.8,
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
    
      bestCluster <- c('none', 'pam', 'hmm')[which.max(c(silThreshold, pamSil, hSil))]
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
    par(cex = 0.75)
    if (cl$sil > silThreshold) {
      par(mfrow = c(consRows(), ceiling(cl$n / consRows())))
      for (i in seq_len(cl$n)) {
        col <- palettes[[min(length(palettes), cl$n)]][i]
        tr <- ape::consensus(r$trees[cl$cluster == i])
        tr$edge.length <- rep.int(1, dim(tr$edge)[1])
        plot(tr, edge.width = 2, font = 1, cex = 1,
             edge.color = col, tip.color = col)
      }
    } else {
      tr <- ape::consensus(r$trees)
      tr$edge.length <- rep.int(1, dim(tr$edge)[1])
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
    min(5L, length(r$trees) - 1L)
  })
  
  nProjDim <- reactive({
    dim(projection())[2]
  })
  
  dims <- debounce(reactive({
    if (mode3D()) 3L else {
      min(5L, maxProjDim())
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
  
  plotContent <- reactive({
    switch(input$plotFormat,
           'cons' = {
             par(mar = rep(0, 4), cex = 0.9)
             ConsensusPlot()
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
  output$plotSpacer <- renderUI({
    tags$div(style = paste0('margin-bottom: ', 
                            (input$plotSize - 600)
                            , 'px'),
             " ")
  })
  
  output$references <- renderUI({
    tagList(
     tags$h2('References for methods used'),
     tags$h3('Tree search'),
     HTML(Brazeau2019, Morphy, Nixon1999, SmithSearch),
     tags$h3('Tree space projection'),
     HTML(paste0(Gower1966, Gower1969, Kaski2003, RCoreTeam,
                 SmithDist, Smith2020, Smith2021, 
                 Venna2001)),
     tags$h3('Clustering'),
     HTML(paste("Cluster consensus trees:", Stockham2002)),
     HTML(paste0('Partitioning around medoids:', Maechler2019,
                 "Hierarchical, minimax linkage:", Bien2011,
                 Murtagh1983)),
     tags$h3("Rogue taxa"),
     HTML(paste0(Smith2022))
    )
  })

}

shinyApp(ui = ui, server = server)