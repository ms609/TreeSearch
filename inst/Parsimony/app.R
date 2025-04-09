# options("TreeSearch.logging" = TRUE) # Log function entry and exit
# options("TreeSearch.write.code" = TRUE) # Show code as it is written to log
logging <- isTRUE(getOption("TreeSearch.logging"))
options(shiny.maxRequestSize = 1024 ^ 3) # Allow max 1 GB files


library("methods", exclude = c("show", "removeClass"))
library("cli")
library("TreeSearch") # load now: inapplicable.datasets required within ui
.DateTime <- function() { # Copy, because not exported
  format(Sys.time(), "%Y-%m-%d %T")
}

suppressPackageStartupMessages({
  library("shiny", exclude = c("runExample"))
  library("shinyjs", exclude = c("runExample"))
})


if (logging) {
  logMsgFile <- file("log.lg", open = "w+")
  LogMsg <- function (...) {
    message(.DateTime(), ": ", ...)
    writeLines(.DateTime(), con = logMsgFile)
    writeLines(paste0("  ", ...), con = logMsgFile)
  }
  Put <- function (..., file) {
    dput(..., file = file)
    writeLines(gsub("<pointer: [^.]+>", "NULL", readLines(file)),
               file)
  }
  PutTree <- function (...) {
    Put(..., file = "tree.lg")
  }
  PutData <- function (...) {
    Put(..., file = "dataset.lg")
  }
} else {
  PutData <- PutTree <- LogMsg <- function (...) {}
}

WriteLoggedCode <- if (isTRUE(getOption("TreeSearch.write.code"))) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    function(txt) {
      for (line in txt) cat(if (substr(trimws(line), 0, 1) == "#") {
        crayon::green("  ", line, "\n")
      } else {
        crayon::yellow("  ", line, "\n")
      })
    }
  } else {
    function(txt) message("       ", txt)
  }
} else {
  function(txt) {}
}

Notification <- function (...) {
  if (!isTRUE(getOption("shiny.testmode"))) {
    showNotification(...)
  }
}

Icon <- function(...) icon(..., class = "fas")

aJiffy <- 42 # ms, default debounce period for input sliders etc
typingJiffy <- 2.5 * aJiffy # slightly slower if might be typing
aFewTrees <- 48L # Too many and rogues / tree space are slowed
NO_OUTGROUP <- "! TREESEARCH_no outgroup specified ."

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

ErrorPlot <- function (...) {
  plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
  text(0, 0, paste0(..., collapse = "\n"),
       col = "#dd6611", font = 2)
}

badToGood <- rev(c("#1AB958", "#23B956", "#2BB954", "#31B952", "#37B850", "#3CB84E", "#41B84C", "#45B74A", "#49B749", "#4DB747", "#51B645", "#54B643", "#58B641", "#5BB53F", "#5FB53D", "#62B53C", "#65B43A", "#68B438", "#6BB336", "#6DB335", "#70B333", "#73B231", "#76B230", "#78B12E", "#7BB12C", "#7DB02B", "#80B029", "#82AF28", "#85AF26", "#87AE25", "#8AAE23", "#8CAD22", "#8EAD21", "#91AC1F", "#93AC1E", "#95AB1D", "#97AB1C", "#9AAA1B", "#9CAA1A", "#9EA919", "#A0A918", "#A2A818", "#A4A717", "#A6A716", "#A8A616", "#AAA616", "#ACA515", "#AEA415", "#B0A415", "#B2A315", "#B4A315", "#B6A216", "#B8A116", "#B9A117", "#BBA017", "#BD9F18", "#BF9F18", "#C19E19", "#C29D1A", "#C49D1B", "#C69C1C", "#C79B1D", "#C99A1E", "#CB9A1F", "#CC9920", "#CE9822", "#CF9823", "#D19724", "#D29625", "#D49626", "#D59528", "#D79429", "#D8932A", "#D9932C", "#DB922D", "#DC912E", "#DD9130", "#DF9031", "#E08F33", "#E18F34", "#E28E35", "#E38D37", "#E58C38", "#E68C3A", "#E78B3B", "#E88A3D", "#E98A3E", "#EA8940", "#EB8841", "#EC8843", "#ED8744", "#EE8746", "#EE8647", "#EF8549", "#F0854A", "#F1844C", "#F2844D", "#F2834F", "#F38350", "#F48252", "#F48253", "#F58155", "#F58157", "#F68058", "#F6805A", "#F77F5B", "#F77F5D", "#F87E5E"))

Reference <- function (authors, year, title, journal = "",
                       volume = NULL, pages = NULL, doi = NULL,
                       publisher = NULL, editors = NULL) {
  nAuth <- length(authors)
  if (nAuth > 1L) {
    authors <- paste(paste0(authors[-nAuth], collapse = ", "), "&amp;", authors[nAuth])
  }
  nEd <- length(editors)
  if (nEd > 1L) {
    editors <- paste(paste0(editors[-nEd], collapse = ", "), "&amp;", editors[nEd])
  } else if (nEd < 1) {
    editors <- ""
  }
  paste0("<p class=\"reference\">", authors, " (", year, "). &ldquo;", title,
         "&rdquo;. ",
         if (editors != "") paste0("In: ", editors, " (eds). ") else "",
         if (journal != "") paste0("<i>", journal, "</i> ") else "",
         if (is.null(volume)) "" else paste0("<b>", volume, "</b>:"),
         if (is.null(publisher)) "" else paste0(publisher, ". "),
         if (is.null(pages)) "" else paste0(paste0(pages, collapse = "&ndash;"), ". "),
         if (is.null(doi)) "" else paste0(
           "doi:<a href=\"https://doi.org/", doi, "\" title=\"CrossRef\">",
           doi, "</a>. "), 
         "</p>")
}


Arthur2007 <- Reference(
  c("Arthur, D.", "Vassilvitskii, S"),
  title = "k-means++: the advantages of careful seeding",
  year = 2007,
  journal = "Proceedings of the Eighteenth Annual ACM-SIAM Symposium on Discrete Algorithms",
  pages = c(1027, 1035)
)
Brazeau2019 <- Reference(c("Brazeau, M.D.", "Guillerme, T.", "Smith, M.R."), 2019,
                           title = "An algorithm for morphological phylogenetic analysis with inapplicable data",
                           journal = "Systematic Biology",
                           volume = 64,
                           pages = c(619, 631),
                         doi = "10.1093/sysbio/syy083")
Bien2011 <- Reference(
  c("Bien, J.", "Tibshirani, R."),
  title = "Hierarchical clustering with prototypes via minimax linkage",
  year = 2011,
  volume = 106,
  doi = "10.1198/jasa.2011.tm10183",
  pages = c(1075, 1084),
  journal = "Journal of the American Statistical Association")
Gower1966 <- Reference(title = "Some distance properties of latent root and vector methods used in multivariate analysis",
                       authors = "Gower, J.C.",
                       year = 1966,
                       volume = 53,
                       pages = c(325, 338),
                       doi = "10.2307/2333639",
                       journal = "Biometrika")
Gower1969 <- Reference(
  title = "Minimum spanning trees and single linkage cluster analysis",
  authors = c("Gower, J.C.", "Ross, G.J.S."),
  year = 1969, volume = 18, pages = c(54, 64), doi = "10.2307/2346439",
  journal = "Journal of the Royal Statistical Society Series C (Applied Statistics)")
Hartigan1979 <- Reference(
  title = "Algorithm AS 136: a <i>K</i>-means clustering algorithm",
  authors = c("Hartigan, J.A.", "Wong, M.A."),
  journal = "Journal of the Royal Statistical Society Series C (Applied Statistics)",
  year = 1979, volume = 28, pages = c(100, 108),
  doi = "10.2307/2346830")
Kaski2003 <- Reference(
  title = "Trustworthiness and metrics in visualizing similarity of gene expression",
  authors = c("Kaski, S.", "Nikkil&auml;, J.", "Oja, M.", "Venna, J.",
             "T&ouml;r&ouml;nen, P.", "Castr&eacute;n, E."),
  year = 2003, volume = 4, pages = 48, doi = "10.1186/1471-2105-4-48",
  journal = "BMC Bioinformatics")
Klopfstein2019 <- Reference(
  title = "Illustrating phylogenetic placement of fossils using RoguePlots: An example from ichneumonid parasitoid wasps (Hymenoptera, Ichneumonidae) and an extensive morphological matrix.",
  authors = c("Klopfstein, S.", "Spasojevic, T."), year = 2019, 
  journal = "PLoS ONE", volume = 14, pages = "e0212942",
  doi = "10.1371/journal.pone.0212942"
)
Maechler2019 <- Reference(
  title = "cluster: cluster analysis basics and extensions", year = 2022,
  authors = c("Maechler, M.", "Rousseeuw, P.", "Struyf, A.", "Hubert, M.", "Hornik, K."),
  journal = "Comprehensive R Archive Network")
Morphy <- Reference(
  c("Brazeau, M.D.", "Smith, M.R.", "Guillerme, T."), 2017,
  "MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability",
  doi = "10.5281/zenodo.815371")
Murtagh1983 <- Reference(
  title = "A survey of recent advances in hierarchical clustering algorithms",
  authors = "Murtagh, F.", year = 1983, volume = 26, pages = c(354, 359),
  doi = "10.1093/comjnl/26.4.354", journal = "The Computer Journal")
Nixon1999 <- Reference(
  "Nixon, K.C.", 1999,
  journal = "Cladistics", volume = 15, pages = c(407, 414),
  title = "The Parsimony Ratchet, a new method for rapid parsimony analysis",
  doi = "10.1111/j.1096-0031.1999.tb00277.x")
Pol2009 <- Reference(
  title = "Unstable taxa in cladistic analysis: identification and the assessment of relevant characters",
  authors = c("Pol, D.", "Escapa, I.H."),
  journal = "Cladistics", 2009, 25, pages = c(515, 527), 
  doi = "10.1111/j.1096-0031.2009.00258.x")
RCoreTeam <- Reference(
  authors = "R Core Team", year = 2020,
  title = "R: A language and environment for statistical computing",
  publisher = "R Foundation for Statistical Computing, Vienna, Austria")
Rousseeuw1987 <- Reference(
  title = "Silhouettes: a graphical aid to the interpretation and validation of cluster analysis",
  author = "Rousseeuw, P.J.", year = 1987,
  journal = "Journal of Computational and Applied Mathematics",
  volume = 20, pages = c(53, 65), doi = "10.1016/0377-0427(87)90125-7"
)
SmithDist <- Reference(
  "Smith, M.R.", "2020a", "TreeDist: distances between phylogenetic trees",
  doi = "10.5281/zenodo.3528123", "Comprehensive R Archive Network")
SmithQuartet <- Reference(
  "Smith, M.R.", 2019,
  "Quartet: comparison of phylogenetic trees using quartet and split measures",
  "Comprehensive R Archive Network", doi = "10.5281/zenodo.2536318")
SmithSearch <- Reference(
  "Smith, M.R.", 2023, "TreeSearch: morphological phylogenetic analysis in R",
  "R Journal", volume = 14, pages = c(305, 315),
  doi = "10.32614/RJ-2023-019")
Smith2020 <- Reference(
  "Smith, M.R.", "2020b",
  "Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees",
  "Bioinformatics", volume = 36, pages = c("5007", "5013"),
  doi = "10.1093/bioinformatics/btaa614")
SmithSpace <- Reference(
  "Smith, M.R.", "2022a", "Robust analysis of phylogenetic tree space",
  "Systematic Biology", 71, pages = c("1255", "1270"),
  doi = "10.1093/sysbio/syab100")
SmithRogue <- Reference(
  "Smith, M.R.", "2022b",
  "Using information theory to detect rogue taxa and improve consensus trees",
  "Systematic Biology", 71, pages = c("1088", "1094"),
  doi = "10.1093/sysbio/syab099")
Stockham2002 <- Reference(
  authors = c("Stockham, C.", "Wang, L.-S.", "Warnow, T."), 2002,
  "Statistically based postprocessing of phylogenetic analysis by clustering",
  "Bioinformatics", 18, c("S285", "S293"),
  doi = "10.1093/bioinformatics/18.suppl_1.S285")

Venna2001 <- Reference(
  title = "Neighborhood preservation in nonlinear projection methods: an experimental study",
  authors = c("Venna, J.", "Kaski, S."), year = 2001, pages = c(485, 491),
  journal = "Lecture Notes in Computer Science: Artificial Neural Networks&mdash;ICANN 2001",
  editors = c("Dorffner, G.", "Bischof, H.", "Hornik, K."),
  publisher = "Springer, Berlin",
  doi = "10.1007/3-540-44668-0_68")














ui <- fluidPage(
  theme = "app.css",
  title = "TreeSearch",
  
  if (isTRUE(getOption("shiny.testmode"))) {
    tags$head(
      tags$style(HTML("#shiny-notification-panel {visibility: hidden;}")
      )
    )
  },
  useShinyjs(),
  column(3,
    fluidRow(
      tags$h1("TreeSearch", style = "margin-top: 0.4em;"),
      selectInput(
        "dataSource",
        "Dataset",
        c("< Load from file below >" = "file",
          "Agnarsson 2004" = "Agnarsson2004",
          "Sun et al. 2018" = "Sun2018",
          "Wills et al. 2012" = "Wills2012",
          if (logging) setNames(names(inapplicable.datasets),
                                names(inapplicable.datasets))
        )
      ),
      fileInput("dataFile",
                tags$span(
                  tags$i(class="fas fa-solid fa-table"),
                  tags$span("Load data from file")
                  ),
                placeholder = "No data file selected"),
      hidden(tags$span(id = "readxl.options",
        selectInput("readxl.sheet", "Excel sheet to read:", "Sheet 1", "Sheet 1"),
        tags$span("First character row & column:"),
        numericInput("readxlSkip",
                     label = NULL,
                     min = 2L, value = 2L, step = 1L),
        numericInput("readxlSkipCols",
                     label = NULL,
                     min = 2L, value = 2L, step = 1L),
        htmlOutput("readxl.chars", style = "clear: both;"),
        htmlOutput("readxl.taxa", style = "clear: both; margin-bottom: 1em;")
      )),
      tags$label("Search", class = "control-label", 
                 style = "display: block; margin-top: -15px;"),
      actionButton("searchConfig", "Configure", icon = Icon("gears")),
      hidden(actionButton("go", "Search", icon = Icon("magnifying-glass"))),
      downloadButton("saveZip", "Save log", icon = Icon("download")),
      fileInput("treeFile",
                label = tags$span(
                  tags$i(class="fas fa-solid fa-tree"),
                  tags$span("Load trees")
                ),
                placeholder = "No tree file selected"),
      textOutput("results"),
      hidden(tags$div(id = "manipulateTreeset",
        numericInput("nTree",
                     label = HTML("Sample <i>n</i> trees from range:"),
                     min = 1L, value = 1L, step = 1L),
        sliderInput("treeRange", label = "", min = 1L, max = 1L,
                    step = 1L, value = c(1, 1)),
        tags$label("Save chosen trees:", class = "control-label"),
        tags$div(style = "display: inline-block",
          downloadButton("saveNwk", "Newick", icon = Icon("download")),
          downloadButton("saveNex", "Nexus", icon = Icon("download"))
        )
      )),
      hidden(
        tags$div(id = "displayConfig",
                 radioButtons("plotFormat", "Display:",
                   list("Characters on trees" = "ind",
                        "Consensus tree" = "cons",
                        "Cluster consensus trees" = "clus",
                        "Tree space" = "space"),
                   # "ind"),
                   "cons"),
                 hidden(sliderInput("whichTree", "Tree to plot", value = 0L,
                                    min = 0L, max = 1L, step = 1L)),
                 hidden(tags$div(id = "treePlotConfig",
                   selectizeInput("outgroup", "Root on:", multiple = TRUE,
                                  choices = list()),
                   selectizeInput(
                     "concordance",
                     "Split support:",
                     choices = list(
                       "None" = "none",
                       "% trees containing" = "p",
                       "Quartet concordance" = "qc",
                       "Clustering concordance" = "clc",
                       "Phylogenetic concordance" = "phc",
                       "Mutual Clustering conc." = "mcc",
                       "Shared Phylog. conc." = "spc"
                     ))
                 )),
                 hidden(tags$div(id = "mapConfig",
                   checkboxGroupInput("mapLines", "Connect:",
                                      choices = list(
                                        "Cluster convex hulls" = "hull",
                                        "Minimum spanning tree" = "mst",
                                        "Trees in sequence" = "seq"
                                      ), selected = c("hull", "mst"))
                 ))
        )
      ),
    ),
  ),
  column(9,
    fluidRow(id = "plotConfig",
      tags$div(id = "plotSizer", 
               tags$span("Plot size:", id = "plotSizeSpan"),
               sliderInput(inputId = "plotSize",
                           label = NULL, width = "200px",
                           min = 100, max = 2000,
                           post = "px", value = 600),
      ),
      tags$div(id = "saveAs", 
               tags$span("Save\ua0plot: "),
               downloadButton("savePlotZip", "R script", icon = Icon("download")),
               downloadButton("savePdf", "PDF", icon = Icon("download")),
               downloadButton("savePng", "PNG", icon = Icon("download"))
      ),
      tags$div(id = "savePlottedTrees",
               downloadButton("savePlotNwk", "Newick", icon = Icon("download")),
               downloadButton("savePlotNex", "Nexus", icon = Icon("download"))
      )
    ),
    fluidRow(
      plotOutput(outputId = "treePlot", height = "600px"),
      hidden(plotOutput("clustCons", height = "200px")),
      hidden(tags$div(id = "charChooser",
        tags$div(
          numericInput("plottedChar", "Character to map:", value = 1L,
                       min = 0L, max = 1L, step = 1L, width = 200),
          selectizeInput("searchChar", "Search characters:", multiple = FALSE,
                         choices = list()),
          checkboxGroupInput("mapDisplay", "", list(
            "Align tips" = "tipsRight",
            "Infer tips" = "updateTips"
          )),
          style = "float: right; width: 200px; margin-left: 2em;"),
        htmlOutput("charMapLegend"),
        htmlOutput("charNotes"),
      )),
      hidden(tags$div(id = "consConfig",
        tags$div(style = "float: right; width: 200px; margin-left: 2em;",
          sliderInput("consP", "Majority:", value = 1,
                      min = 0.5, max = 1, width = 200),
          numericInput("keepNTips", "Tips to show:", value = 0L,
                       min = 3L, max = 2L, step = 1L, width = 200),
          selectizeInput("neverDrop", "Never drop:", multiple = TRUE,
                         choices = c())
                 ),
        tags$div(id = "consLegend",
                 tags$span(id = "instabLegend",
                          tagList(
                            tags$span(class = "legendLeft", "Stable"),
                            tags$span(class = "infernoScale legendBar", "\ua0"),
                            tags$span(class = "legendRight", "Unstable"),
                          )),
                 htmlOutput("branchLegend", inline = TRUE)),
        tags$div(id = "droppedTips",
          selectInput("excludedTip", "Show excluded tip", choices = list())),
        tags$div(id = "droppedList", style = "float: left;"),
      )),
      hidden(tags$div(id = "clusLegend",
                      htmlOutput("instabLegend2", inline = TRUE)
      )),
      hidden(tags$div(id = "clusConfig",
                      style = "float: right; width: 200px; margin-left: 2em;",
          sliderInput("clThresh", "Cluster threshold:", value = 0.5,
                      min = 0, max = 1, width = 200),
          selectInput("distMeth", "Distance method:", selected = "cid",
                      choices = list("Clustering Information" = "cid",
                                     "Phylogenetic information" = "pid",
                                     "Matching split info" = "msid",
                                     "Robinson-Foulds (fast, iffy)" = "rf",
                                     "Quartet (slower)" = "qd"),
                      width = 200)
      )),
      hidden(tags$div(
        id = "spaceConfig",
        tags$div(id = "spaceLegend",
                 style = "float: left;",
                 plotOutput(outputId = "pcQuality",
                            height = "72px", width = "240px"),
                 htmlOutput("stressLegend", inline = TRUE)
        ),
        tags$div(
          style = "float: right; width: 200px; margin-left: 2em;",
          sliderInput("spaceDim", "Dimensions:", value = 5,
                      min = 1, max = 12, step = 1, width = 200),
          selectInput("spaceCol", "Colour trees by:",
                      list("Cluster membership" = "clust",
                           "Parsimony score" = "score",
                           "When first found" = "firstHit")),
          selectInput("spacePch", "Plotting symbols:",
                      selected = "relat",
                      list("Cluster membership" = "clust",
                           "Relationships" = "relat",
                           "Tree name" = "name")),
          selectizeInput("relators", "Show relationship between:",
                         choices = list(), multiple = TRUE),
        ),
      )),
      htmlOutput("references", style = "clear: both;"),
    ),
  )
)


















Enquote <- function(x, ...) {
  if (mode(x) == "character") {
    paste0("\"", x, "\"")
  } else {
    signif(x, ...)
  }
}

EnC <- function(...) {
  if (length(...) == 1) {
    Enquote(...)
  } else {
    paste0("c(", paste(sapply(..., Enquote), collapse = ", "), ")")
  }
}

server <- function(input, output, session) {
  
  r <- reactiveValues(
    dataFileVisible = TRUE,
    ignoreTreeRange = TRUE,
    ignoreNTree = TRUE,
    nTree = 0L,
    oldOutgroup = NO_OUTGROUP,
    sortTrees = FALSE, # May be arranged nicely in input files
    treeRange = c(1L, 1L),
    updatingTrees = FALSE # TODO DELETE?
  )
  
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
          paste0(
            "trees <- trees[unique(as.integer(seq.int(",
            r$treeRange[1], ", ", r$treeRange[2],
            ", length.out = ", r$nTree, ")))]"
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
  DataFileName <- function(n) if (length(n)) {
    paste0("dataFile-", TwoWide(n), ".txt")
  }
  ExcelFileName <- function(n) if (length(n)) {
    paste0("excelFile-", TwoWide(n), ".xlsx")
  }
  TreeFileName <- function(n) if (length(n)) {
    paste0("treeFile-", TwoWide(n), ".txt")
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
  
  if (!requireNamespace("TreeDist", quietly = TRUE)) {
    install.packages("TreeDist")
  }
  
  library("TreeTools", quietly = TRUE)
  library("TreeDist")
  library("TreeSearch")
  
  BeginLog()
  
  library("future")
  library("promises")
  plan(multisession)
  
  startOpt <- options("cli.progress_show_after" = 0.1)
  
  
  LogMsg("Started server")
  
  
  ##############################################################################
  # Load data
  ##############################################################################
  
  tipLabels <- reactive({r$trees[[1]][["tip.label"]]})
  
  nChars <- reactive({
    if (HaveData()) {
      as.integer(length(attr(r$dataset, "index")))
    } else {
      0L
    }
  })
  
  TaxonOrder <- reactive({
    if (HaveData()) {
      names(r$dataset)
    } else {
      tipLabels()
    }
  })
  
  DatasetMatchesTrees <- reactive({
    length(intersect(names(r$dataset), tipLabels())) == length(r$dataset)
  })
  
  UpdateData <- reactive({
    source <- input$dataSource
    if (source == "file") {
      if (!r$dataFileVisible) {
        showElement("dataFile")
        r$dataFileVisible <- TRUE
        runjs("console.log($('#dataFile-label'))")
        runjs(paste0(
          "$('#dataFile-label').parent()",
          ".css({'outline': 'dashed #428bca 20px', ", 
          "'width': '100%'})", 
          ".animate({'outline-width': '0px'}, 'slow');"))
        return()
      }
      
      fileInput <- input$dataFile
      r$dataset <- NULL
      r$chars <- NULL
      if (is.null(fileInput)) {
        # How can this be?
        Notification(type = "error", "No data file selected")
        return("No data file selected.")
      }
      dataFile <- fileInput$datapath
      if (is.null(dataFile)) {
        Notification(type = "error", "No data file found.")
        return ("No data file specified.")
      }
      
      LogMsg("UpdateData(): from file")
      r$sortTrees <- FALSE # Trees loaded from dataset may be in sequence
      r$readDataFile <- NULL
      
      if (length(grep("\\.xlsx?$", dataFile))) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          install.packages("readxl")
        }
        showElement("readxl.options", anim = TRUE)
        
        r$dataset <- tryCatch({
          sheets <- readxl::excel_sheets(dataFile)
          updateSelectInput(session,
                            inputId = "readxl.sheet",
                            choices = setNames(sheets, sheets),
                            selected = if (input$readxl.sheet %in% sheets) {
                                input$readxl.sheet
                              } else {
                                sheets[1]
                              })

          tibble <- readxl::read_excel(
            path = dataFile,
            sheet = match(input$readxl.sheet, sheets, nomatch = 1L),
            skip = max(0L, input$readxlSkip - 2L),
            .name_repair = "minimal",
            col_types = "text"
          )
          
          firstCol <- input$readxlSkipCols - 1L
          chars <- colnames(tibble)[-seq_len(firstCol)]
          taxNames <- gsub(" ", "_", trimws(unlist(tibble[, firstCol])))
          output$readxl.taxa <- renderUI(HTML(paste(
            "<em>Taxon names</em>:",
            paste(taxNames[1:3], collapse = ", "),
            "...\n")))
          output$readxl.chars <- renderUI(HTML(paste(
            "<em>Character names</em>:",
            # not r$chars, which may be modified before output updated
            paste(chars[1:3], collapse = ", "),
            "..."
          )))
          r$chars <- chars
          
          dat <- as.matrix(tibble[, -seq_len(firstCol)])
          rownames(dat) <- taxNames
          dat <- MatrixToPhyDat(dat)
          if (attr(dat, "nr") == 0) {
            stop("No characters loaded; throw error")
          }
          
          # Lines that could cause an error must come before log
          
          LogComment("Load data from spreadsheet", 2)
          if (r$excelFiles == 0 ||
              tools::md5sum(dataFile) != 
              tools::md5sum(paste0(tempdir(), "/", LastFile("excel")))) {
            CacheInput("excel", dataFile)
          }
          LogCode(c(
            paste0("dataFile <- \"", LastFile("excel"), "\""),
            "excelSheet <- readxl::read_excel(",
            "  path = dataFile,",
            paste0("  sheet = ", match(input$readxl.sheet, sheets, 1L), ","),
            paste0("  skip = ", max(0L, input$readxlSkip - 2L), ","),
            "  .name_repair = \"minimal\",",
            "  col_types = \"text\"",
            ")",
            paste0("dat <- as.matrix(excelSheet[, -seq_len(", firstCol, ")])"),
            paste0("rownames(dat) <- unlist(excelSheet[, ", firstCol, "])"),
            "dataset <- MatrixToPhyDat(dat)"
          ))
          
          # Return:
          dat
        }, error = function(e) {
          NULL
        })
      } else {
        hideElement("readxl.options")
      }
      
      if (is.null(r$dataset)) suppressWarnings({
        r$dataset <- tryCatch({
          r$readDataFile <- "ReadTntAsPhyDat(dataFile)"
          
          # Return:
          ReadTntAsPhyDat(dataFile)
        }, error = function(e) tryCatch({
          r$chars <- tryCatch(
            ReadCharacters(dataFile),
            error = function(e) {
              Notification(type = "error", "Error reading characters from file")
              # Return:
              NULL
            })
          
          r$charNotes <- tryCatch(
            ReadNotes(dataFile),
            error = function(e) {
              Notification(type = "error", "Error reading character notes")
              # Return:
              NULL
            })
          
          r$readDataFile <- "ReadAsPhyDat(dataFile)"
          
          # Return:
          ReadAsPhyDat(dataFile)
        }, error = function(e) {
          r$readDataFile <- NULL
          # Return:
          NULL
        }))
        
        if (!is.null(r$dataset)) {
          LogComment("Load data from file", 2)
          CacheInput("data", dataFile)
          LogCode(c(
            paste0("dataFile <- \"", LastFile("data"), "\""),
            paste0("dataset <- ", r$readDataFile)
          ))
        }
      })
    } else {
      LogMsg("UpdateData(): from package")
      
      r$sortTrees <- TRUE # Nicer plots
      
      r$dataFileVisible <- FALSE
      hideElement("dataFile")
      
      dataFile <- system.file(paste0("datasets/", source, ".nex"),
                              package = "TreeSearch")
      CacheInput("data", dataFile)
      r$chars <- ReadCharacters(dataFile)
      r$charNotes <- ReadNotes(dataFile)
      r$readDataFile <- "ReadAsPhyDat(dataFile)"
      r$dataset <- ReadAsPhyDat(dataFile)
      LogComment("Load dataset file from TreeSearch package")
      LogCode(c(
        paste0("dataFile <- system.file(\"datasets/", source,
               ".nex\", package = \"TreeSearch\")"),
        "dataset <- ReadAsPhyDat(dataFile)"
      ))
    }
    if (is.null(r$dataset)) {
      Notification(type = "error", "Could not read data from file")
      
      updateNumericInput(session, "plottedChar", min = 0L,
                         max = 0L, value = 0L)
      updateSelectizeInput(session, "searchChar", choices = NULL)
      return ("Could not read data from file")
    } else {
      Notification(type = "message", 
                       paste("Loaded", nChars(), "characters and",
                             length(r$dataset), "taxa"))
      
      updateNumericInput(session, "plottedChar", min = 0L,
                         max = nChars(), value = 1L)
      updateSelectizeInput(session, "searchChar",
                           choices = paste0(seq_len(nChars()), ": ", 
                                            colnames(r$chars)),
                           selected = "",
                           server = TRUE)
    }
    
    tryCatch({
      dataFileTrees <- read.nexus(dataFile)
      LogComment("Read trees from dataset file")
      LogCode("newTrees <- read.nexus(dataFile)")
      UpdateAllTrees(dataFileTrees)
      CacheInput("tree", dataFile)
      r$readTreeFile <- "read.nexus(treeFile)"
    }, error = function (e) NULL)
    if (!AnyTrees() || !DatasetMatchesTrees()) {
      updateActionButton(session, "go", "New search")
    } else {
      show("displayConfig")
    }
    
    DisplayTreeScores()
  })
  
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
      updateSliderInput(inputId = "spaceDim", max = max(1L, maxProjDim()),
                        value = min(maxProjDim(), input$spaceDim))
    }
    updateSelectizeInput(inputId = "neverDrop", choices = tipLabels(),
                         selected = input$neverDrop)
    UpdateOutgroupInput()
    updateSelectizeInput(inputId = "relators", choices = tipLabels(),
                         selected = input$relators)
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
  
  ##############################################################################
  # Event listeners
  ##############################################################################
  
  observeEvent(input$dataSource, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$dataFile, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxl.sheet, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxlSkip, UpdateData(), ignoreInit = TRUE)
  observeEvent(input$readxlSkipCols, UpdateData(), ignoreInit = TRUE)
  
  observeEvent(r$dataset, {
    r$dataHash <- rlang::hash(r$dataset)
  })
  observeEvent(input$plotSize, {
    px <- paste0("'", input$plotSize, "px'")
    runjs(paste0("$('#treePlot').css({height: ", px, ", width: ", px, "});"))
  })
  
  observeEvent(input$searchConfig, {
    #updateSelectInput(session, "character.weight",
    #                  selected = input$character.weight)
    updateSelectInput(session, "implied.weights",
                      selected = input$implied.weights)
    updateSliderInput(session, "concavity", value = input$concavity)
    updateNumericInput(session, "epsilon", value = input$epsilon)
    updateSliderInput(session, "ratchIter", value = input$ratchIter)
    updateSliderInput(session, "tbrIter", value = input$tbrIter)
    updateSliderInput(session, "maxHits", value = input$maxHits)
    updateSliderInput(session, "startIter", value = input$startIter)
    updateSliderInput(session, "finalIter", value = input$finalIter)
    showModal(modalDialog(
      easyClose = TRUE,
      fluidPage(column(6,
      tagList(
        #selectInput("character.weight", "Character weighting",
        #                  list("Equal" = "equal"), "equal"),
              selectInput("implied.weights", "Step weighting", 
                         list("Implied" = "on", "Profile" = "prof",
                              "Equal" = "off"), "on"),
              sliderInput("concavity", "Step weight concavity constant", min = 0L,
                         max = 3L, pre = "10^", value = 1L),
              numericInput("epsilon", "Keep if suboptimal by \u2264", min = 0,
                          value = 0),
              sliderInput("ratchIter", "Ratchet iterations", min = 0L,
                          max = 50L, value = 6L, step = 1L),
              sliderInput("timeout", "Maximum run duration", min = 1,
                          max = 600, value = 30, post = "min", step = 1),
      )), column(6, 
             tagList(
              sliderInput("maxHits", "Maximum hits", min = 0L, max = 5L,
                         value = 2L, pre = "10^"),
              sliderInput("tbrIter", "TBR depth", min = 1L, max = 20L,
                          value = 1L, step = 1L),
              sliderInput("startIter", "First iteration extra depth", min = 1L,
                          max = 10L, value = 3L, pre = "\ud7"),
              sliderInput("finalIter", "Final iteration extra depth", min = 1L,
                          max = 10L, value = 1L, pre = "\ud7"),
              selectizeInput("searchWithout", "Exclude taxa", DatasetTips(),
                             r$searchWithout, multiple = TRUE)
             ))
      ),
      title = "Tree search settings",
      footer = tagList(modalButton("Close", icon = Icon("rectangle-xmark")),
                       actionButton("modalGo", icon = Icon("magnifying-glass"), 
                                    if(length(r$trees)) {
                                      "Continue search" 
                                    } else {
                                      "Start search"
                                    }))
    ))
    show("go")
  })
  
  observeEvent(input$treeFile, {
    tmpFile <- input$treeFile$datapath
    newTrees <- tryCatch({
        r$readTreeFile <- "read.tree(treeFile)"
        LogMsg("Trying read.tree()")
        read.tree(tmpFile)
      },
      error = function (x) tryCatch({
        r$readTreeFile <- "read.nexus(treeFile)"
          LogMsg("Trying read.nexus()")
          read.nexus(tmpFile)
        },
        error = function (err) tryCatch(
          {
            if (grepl("NA/NaN argument", err)) {
              LogMsg("Terminating tree block")
              # Unterminated tree block, perhaps because a search is ongoing
              withEnd <- tempfile()
              on.exit(unlink(withEnd))
              writeLines(c(readLines(tmpFile), "\nEND;"), withEnd)
              read.nexus(withEnd)
            } else {
              stop("Next handler, please")
            }
          },
          error = function (x) tryCatch({
              r$readTreeFile <- "ReadTntTree(treeFile)"
              ReadTntTree(tmpFile)
            }, warning = function (x) tryCatch({
              Notification(as.character(x), type = "warning")
              tryLabels <- TipLabels(r$dataset)
              if (length(tryLabels) > 2) {
                Notification("Inferring tip labels from dataset",
                                 type = "warning")
                r$readTreeFile <- 
                  "ReadTntTree(treeFile, tipLabels = TipLabels(dataset))"
                ReadTntTree(tmpFile, tipLabels = tryLabels)
              } else {
                NULL
              }
            }, error = NULL
            )
          )
        )
      )
    )
    if (is.null(newTrees)) {
      Notification("Trees not in a recognized format", type = "error")
    } else {
      LogComment("Load tree from file", 2)
      CacheInput("tree", tmpFile)
      LogCode(paste0("treeFile <- \"", LastFile("tree"), "\""))
      LogCode(paste0("newTrees <- ", r$readTreeFile))
      
      UpdateAllTrees(newTrees) # updates r$trees
      
      removeModal()
      Notification(paste("Loaded", length(r$trees), "trees"), type = "message")
      updateActionButton(session, "modalGo", "Continue search")
      updateActionButton(session, "go", "Continue")
      show("displayConfig")
    }
    
  })
  
  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           "on" = show("concavity"),
           hide("concavity")
    )
    DisplayTreeScores()
  })
  
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
    output$results <- renderText(msg)
    msg
  }
  
  observeEvent(input$concavity, {
    DisplayTreeScores()
  }, ignoreInit = TRUE)
  
  TipsInTree <- reactive({
    if (AnyTrees()) {
      length(r$trees[[1]]$tip.label)
    } else {
      0L
    }
  })
  
  UpdateKeepNTipsRange <- reactive({
    if (AnyTrees() && "consConfig" %in% r$visibleConfigs) {
      nTip <- TipsInTree()
      LogMsg("UpdateKeepNTipsRange(", input$keepNTips, " -> ", nTip, ")")
      r$keepNTips <- nNonRogues()
      if (r$keepNTips != input$keepNTips) {
        r$oldkeepNTips <- input$keepNTips
      }
      updateNumericInput(inputId = "keepNTips",
                         label = paste0("Tips to show (/", nTip, "):"),
                         min = max(3L, length(input$neverDrop)),
                         max = nTip,
                         value = nNonRogues())
    }
  })
  
  UpdateExcludedTipsInput <- reactive({
    if (AnyTrees() && "consConfig" %in% r$visibleConfigs) {
      LogMsg("UpdateExcludedTipsInput()")
      dropList <- dropSeq()[seq_along(DroppedTips())]
      updateSelectInput(inputId = "excludedTip",
                        choices = dropList,
                        selected = if(input$excludedTip %in% DroppedTips())
                          input$excludedTip else dropSeq()[1])
      html("droppedList",
           paste0("<label class=\"control-label\">Dropped tips:</label>", 
                  "<ul>", 
                  paste0("<li style=\"color: ", TipCols()[dropList], "\">",
                         dropList, "</li>", collapse = "\r\n"),
                  "</ul>"))
    }
  })
  
  UpdateDroppedTaxaDisplay <- reactive({
    LogMsg("UpdateDroppedTaxaDisplay()")
    if ("consConfig" %in% r$visibleConfigs) {
      if (length(DroppedTips())) {
        UpdateExcludedTipsInput()
        if ("droppedTips" %in% r$visibleConfigs) {
          show("droppedTips")
        }
        if ("droppedList" %in% r$visibleConfigs) {
          show("droppedList")
        }
      } else {
        hide("droppedTips")
        hide("droppedList")
      }
    }
  })
  
  observeEvent(r$visibleConfigs, {
    UpdateDroppedTaxaDisplay()
  })
  
  UpdateOutgroupInput <- reactive({
    if (AnyTrees() && "treePlotConfig" %in% r$visibleConfigs) {
      LogMsg("UpdateOutgroupInput()")
      r$outgroup <- intersect(r$outgroup, KeptTips())
      if (length(r$outgroup) == 0) {
        r$outgroup <- if (HaveData()) {
          intersect(names(r$dataset), KeptTips())[1]
        } else {
          KeptTips()[1]
        }
      }
      
      if (!identical(sort(r$outgroup), sort(input$outgroup))) {
        r$oldOutgroup <- if (is.null(input$outgroup)) {
          NO_OUTGROUP
        } else {
          input$outgroup
        }
      }
      
      updateSelectizeInput(
        inputId = "outgroup",
        selected = r$outgroup,
        choices = KeptTips()
        )
    }
  })
  
  observeEvent(input$implied.weights, {
    switch(input$implied.weights,
           "on" = show("concavity"),
           hide("concavity")
    )
  })
  
  ShowConfigs <- function (visible = character(0)) {
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
                       "ind" = c("whichTree", "charChooser",
                                 "treePlotConfig"),
                       "cons" = c("consConfig", "droppedTips",
                                  "savePlottedTrees",
                                  "treePlotConfig", "branchLegend"),
                       "clus" = c("clusConfig", "clusLegend",
                                  "savePlottedTrees",
                                  "consConfig", "droppedList",
                                  "treePlotConfig"),
                       "space" = c("clusConfig", "clusLegend",
                                   "spaceConfig", "mapConfig"),
                       ""))
  })
  
  
  output$branchLegend <- renderUI({
    if (!AnyTrees()) {
      return()
    }
    LogMsg("renderUI(branchLegend)")
    on.exit(LogMsg("/renderUI(branchLegend)"))
    kept <- KeptTips()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      consTrees <- lapply(r$trees, DropTip, setdiff(dropped, input$excludedTip))
      plotted <- TreeTools::RoguePlot(
        trees = consTrees,
        tip = input$excludedTip,
        p = consP(),
        plot = FALSE
      )
      tagList(
        tags$span(class = "legendLeft", "1 tree"),
        tags$span(id = "blackToGreen", class = "legendBar", "\ua0"),
        tags$span(class = "legendRight",
                  paste(max(c(plotted$onEdge, plotted$atNode)), "trees")),
      )
    }
  })
  
  concavity <- reactive({
    kExp <- if (length(input$concavity)) input$concavity else 1
    switch(weighting(),
           "on" = 10 ^ kExp,
           "off" = Inf,
           "prof" = "Profile")
  })
  
  tolerance <- reactive({
    if (input$epsilon == 0) {
      sqrt(.Machine$double.eps)
    } else {
      input$epsilon
    }
  })
  
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
            firstOptimal <- which.min(scores())
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
      LogCode(c(
        "newTrees <- MaximizeParsimony(",
        if (length(r$searchWithout)) {
          paste0(
            "  dataset[setdiff(names(dataset), ", EnC(r$searchWithout), ")]"
          )
        } else {
          "  dataset,"
        },
        "  tree = startTree,",
        paste0("  concavity = ", Enquote(concavity()), ","),
        paste0("  ratchIter = ", input$ratchIter, ","), 
        paste0("  tbrIter = ", input$tbrIter, ","), 
        paste0("  maxHits = ", ceiling(10 ^ input$maxHits), ","), 
        paste0("  maxTime = ", input$timeout, ","),
        paste0("  startIter = ", input$startIter, ","),
        paste0("  finalIter = ", input$finalIter, ","),
        if (input$epsilon > 0) paste0("  tolerance = ", tolerance(), ","),
        "  verbosity = 4",
        ")"))
      newTrees <- withProgress(
        MaximizeParsimony(r$dataset[SearchTips()],
                          tree = startTree,
                          concavity = concavity(),
                          ratchIter = input$ratchIter,
                          tbrIter = input$tbrIter,
                          maxHits = ceiling(10 ^ input$maxHits),
                          maxTime = input$timeout,
                          startIter = input$startIter,
                          finalIter = input$finalIter,
                          tolerance = tolerance(),
                          verbosity = 4L),
        value = 0.85, message = "Finding MPT",
        detail = paste0(ceiling(10^input$maxHits), " hits; ", wtType())
      )
      r$sortTrees <- TRUE # No meaning in order; display nicely
      LogComment("Overwrite any previous trees with results")
      LogCode(c(
        "if (inherits(newTrees, \"phylo\")) {",
        "  trees <- list(newTrees)",
        "  attr(trees, \"firstHit\") <- attr(newTrees, \"firstHit\")",
        "  attr(trees[[1]], \"firstHit\") <- NULL",
        "}"
      ))
      UpdateAllTrees(newTrees)
      if (inherits(newTrees, "phylo")) {
        attr(r$trees, "firstHit") <- attr(newTrees, "firstHit")
        attr(r$trees[[1]], "firstHit") <- NULL
      }
      
      updateSliderInput(session, "whichTree", min = 0L,
                        max = length(r[["trees"]]), value = 0L)
      
      updateActionButton(session, "go", "Continue")
      updateActionButton(session, "modalGo", "Continue search")
      show("displayConfig")
    }
  }
  
  observeEvent(input$searchWithout, {
    r$searchWithout <- input$searchWithout
  }, ignoreInit = TRUE)
  
  observeEvent(input$go, StartSearch(), ignoreInit = TRUE)
  observeEvent(input$modalGo, {
    removeModal()
    StartSearch()
  }, ignoreInit = TRUE)
  
  UserRoot <- function(tree) {
    outgroupTips <- intersect(r$outgroup, tree$tip.label)
    if (length(outgroupTips)) {
      # DELETE? tr <- deparse(substitute(tree))
      RootTree(tree, outgroupTips)
    } else {
      tree
    }
  }
  
  LogUserRoot <- function(tree = "cons", dropped = character(0)) {
    outgroupTips <- setdiff(r$outgroup, dropped)
    if (length(outgroupTips)) {
      LogCommentP("Root tree")
      LogCodeP(paste0(tree, " <- RootTree(", tree, ", ", EnC(outgroupTips), ")"))
    }
  }
  
  PlottedChar <- debounce(reactive({
    typed <- max(0L, as.integer(input$plottedChar), na.rm = TRUE)
    if (nChars() > 0 && typed > nChars()) {
      Notification(type = "warning",
                   paste("Dataset contains", nChars(), "characters.")
      )
      updateNumericInput(session, "plottedChar", value = nChars())
    }
    min(typed, nChars())
  }), aJiffy)
  
  observeEvent(PlottedChar(), {
    if (PlottedChar() > 0) {
      showElement("mapDisplay")
    } else {
      hideElement("mapDisplay")
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$searchChar, {
    searchResult <- as.numeric(strsplit(input$searchChar, ": ")[[1]][1])
    if (!is.na(searchResult)) {
      updateNumericInput(session, "plottedChar", value = searchResult)
    }
  })
  
  whichTree <- debounce(reactive(input$whichTree), aJiffy)
  
  PlottedTree <- reactive({
    if (length(r$trees) > 0L) {
      plottedTree <- if (whichTree() > 0) {
        r$trees[[whichTree()]]
      } else {
        Consensus(r$trees, p = 1)
      }
      plottedTree <- UserRoot(plottedTree)
      plottedTree <- SortEdges(plottedTree)
      if (!("tipsRight" %in% input$mapDisplay)) {
        plottedTree$edge.length <- rep_len(2, dim(plottedTree[["edge"]])[[1]])
      }
      plottedTree
    }
  })
  LogPlottedTree <- function() {
    if (whichTree() > 0) {
      LogCodeP(paste0("plottedTree <- trees[[", whichTree(), "]]"))
    } else {
      LogCodeP("plottedTree <- Consensus(trees, p = 1)")
    }
    LogUserRoot("plottedTree")
    if (!("tipsRight" %in% input$mapDisplay)) {
      LogCommentP("Set uniform edge length", 0)
      LogCodeP(
        "plottedTree$edge.length <- rep.int(2, nrow(plottedTree$edge))"
      )
    }
    LogSortEdges("plottedTree")
  }
  
  Instab <- reactive({
    TipInstability(r$trees)
  })
  
  dropSeq <- reactive({
    LogMsg("dropSeq()")
    Rogues()$taxon[-1]
  })
  
  stableCol <- reactive({
    Rogue::ColByStability(r$trees)
  })
  
  Rogues <- bindCache(reactive({
    if (AnyTrees() && inherits(r$trees, "multiPhylo")) {
      LogComment("Check for rogue taxa", 2)
      LogComment(paste0(
        "Use RogueTaxa() in place of QuickRogue() for a more complete ",
        "analysis"))
      LogCode(c(
        "rogues <- Rogue::QuickRogue(",
        "  trees,",
        if (length(input$neverDrop)) paste0(
          "  neverDrop = ", EnC(input$neverDrop), ","
        ),
        "  fullSeq = TRUE,",
        paste0("  p = ", Enquote(consP())),
        ")",
        "print(rogues) # Detailed results of rogue analysis",
        "print(rogues$taxon[-1]) # Sequence of taxa to drop"
      ))
      withProgress(
        message = "Identifying rogues", value = 0.99,
        rogues <- Rogue::QuickRogue(r$trees, neverDrop = input$neverDrop,
                          fullSeq = TRUE, p = consP())
      )
      # TODO delete once Rogue 2.1.2 released -- return QuickRogue above.
      rogues[!rogues$taxon %in% input$neverDrop, ]
    } else {
      data.frame(num = 0, taxNum = NA_integer_, taxon = NA_character_,
                 rawImprovement = NA_real_, IC = 0)
    }
  }), r$treeHash, input$neverDrop, consP())
  
  unitEdge <- reactive({
    TRUE
  })
  
  SortEdges <- function (tr, force = FALSE) {
    if (force || r$sortTrees) {
      # Return:
      SortTree(tr, order = TaxonOrder())
    } else {
      # Return:
      tr
    }
  }
  LogSortEdges <- function(tr) (
    if (r$sortTrees) {
      LogCommentP("Rotate nodes, to display clades in order of size", 0)
      LogCodeP(paste0(
        tr, " <- SortTree(", tr, ", order = ", 
        if (HaveData()) {
          "names(dataset)"
        } else {
          "trees[[1]]$tip.label"
        },
        ")"
      ))
    }
  )
  
  nNonRogues <- reactive({
    LogMsg("nNonRogues()")
    on.exit(LogMsg("nNonRogues: ", nrow(Rogues()) - which.max(Rogues()$IC)))
    nrow(Rogues()) - which.max(Rogues()$IC)
  })
  
  TipCols <- reactive(stableCol()) # TODO allow user to choose how to colour
  
  TipColLegend <- function() {
    PlotTools::SpectrumLegend(
      "bottomleft", horiz = TRUE, inset = 0.01, bty = "n", xpd = NA,
      palette = hcl.colors(131, "inferno")[1:101],
      legend = c("Stable", "Unstable"),
      title = "Leaf stability",
      title.font = 2
    )
  }
  
  consP <- debounce(reactive(signif(input$consP)), 50)
  observeEvent(consP(), {
    if (AnyTrees()) {
      LogMsg("Observed consP()")
      UpdateKeepNTipsRange()
      UpdateDroppedTaxaDisplay()
      r$concordance <- list()
    }
  }, ignoreInit = TRUE)
  
  concordance <- bindCache(reactive({
    LogMsg("concordance()")
    # Return:
    switch(input$concordance,
          "p" = SplitFrequency(r$plottedTree, r$trees) / length(r$trees),
          "qc" = QuartetConcordance(r$plottedTree, r$dataset),
          "mcc" = MutualClusteringConcordance(r$plottedTree, r$dataset),
          "spc" = SharedPhylogeneticConcordance(r$plottedTree, r$dataset),
          "clc" = ClusteringConcordance(r$plottedTree, r$dataset),
          "phc" = PhylogeneticConcordance(r$plottedTree, r$dataset),
          NULL
    )
  }), r$plottedTree, r$treeHash, r$dataHash, input$concordance)
  
  LabelConcordance <- \() {
    LogMsg("LabelConcordance()")
    if (input$concordance != "none" &&
        !is.null(r$plottedTree)) {
      LabelSplits(r$plottedTree, signif(concordance(), 3),
                  col = SupportColor(concordance()),
                  frame = "none", pos = 3L)
    }
  }
  
  LogConcordance <- function(plottedTree = "plottedTree") {
    if (input$concordance != "none") {
      LogCommentP("Calculate split concordance", 1)
      concCode <- switch(
        input$concordance,
        "p"   = paste0("SplitFrequency(", plottedTree,
                       ", trees) / length(trees)"),
        "qc"  = paste0("QuartetConcordance(", plottedTree, ", dataset)"),
        "clc" = paste0("ClusteringConcordance(", plottedTree, ", dataset)"),
        "phc" = paste0("PhylogeneticConcordance(", plottedTree, ", dataset)"),
        "mcc" = paste0("MutualClusteringConcordance(", plottedTree,
                       ", dataset)"),
        "spc" = paste0("SharedPhylogeneticConcordance(", plottedTree,
                       ", dataset)"),
        NULL
      )
      LogCodeP(paste0("concordance <- ", concCode))
      LogCommentP("Annotate splits by concordance", 1)
      LogCodeP("LabelSplits(",
              paste0("  tree = ", plottedTree, ","),
              "  labels = signif(concordance, 3),",
              "  col = SupportColor(concordance),",
              "  frame = \"none\",",
              "  pos = 3",
              ")")
    }
  }
  
  observeEvent(input$keepNTips, {
    if (!is.null(r$oldkeepNTips)) {
      if (!identical(input$keepNTips, r$oldkeepNTips)) {
        r$oldkeepNTips <- NULL
      }
    } else {
      LogMsg("Observed input$keepNTips -> ", EnC(input$keepNTips))
      r$keepNTips <- max(length(input$neverDrop), 3L,
                         min(input$keepNTips, TipsInTree()))
      UpdateOutgroupInput()
      UpdateDroppedTaxaDisplay()
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$neverDrop, {
    LogMsg("Observed input$neverDrop -> ", EnC(input$neverDrop))
    UpdateKeepNTipsRange()
    UpdateOutgroupInput()
    UpdateDroppedTaxaDisplay()
  }, ignoreInit = TRUE)
  
  observeEvent(input$outgroup, {
    if (!is.null(r$oldOutgroup)) {
      if (!identical(input$outgroup, r$oldOutgroup)) {
        r$oldOutgroup <- NULL
      }
    } else {
      LogMsg("Observed input$outgroup -> ", EnC(input$outgroup))
      r$outgroup <- input$outgroup
    }
  }, ignoreInit = TRUE)
  
  DatasetTips <- reactive(names(r$dataset))
  SearchTips <- reactive(setdiff(DatasetTips(), r$searchWithout))
  
  KeptTips <- reactive({
    LogMsg("KeptTips()")
    n <- r$keepNTips
    maxN <- length(tipLabels())
    if (is.na(n) || is.null(n)) {
      n <- maxN
    }
    if (n < 3L) {
      n <- 3L 
    }
    nNeverDrop <- length(input$neverDrop)
    if (n < nNeverDrop) {
      n <- nNeverDrop
    }
    nFromDropSeq <- n - nNeverDrop
    
    # Return:
    if (nFromDropSeq > length(dropSeq())) {
      c(input$neverDrop, dropSeq())
    } else {
      c(input$neverDrop, rev(dropSeq())[seq_len(nFromDropSeq)])
    }
  })
  
  DroppedTips <- reactive({
    LogMsg("DroppedTips()")
    if (length(KeptTips()) > 1) {
      setdiff(tipLabels(), KeptTips())
    } else {
      character(0)
    }
  })
  
  ConsensusPlot <- function() {
    LogMsg("ConsensusPlot()")
    on.exit(LogMsg("/ConsensusPlot()"))
    
    par(mar = rep(0, 4), cex = 0.9)
    kept <- KeptTips()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      
      if (length(setdiff(dropped, input$excludedTip))) {
        consTrees <- lapply(r$trees, DropTip,
                            setdiff(dropped, input$excludedTip))
      } else {
        consTrees <- r$trees
      }
      
      plotted <- TreeTools::RoguePlot(
        consTrees,
        input$excludedTip,
        p = consP(),
        edgeLength = 1,
        outgroupTips = r$outgroup,
        tip.color = TipCols()[intersect(consTrees[[1]]$tip.label, kept)]
      )
      r$plottedTree <- plotted$cons
      
      LabelConcordance()
    } else {
      without <- intersect(dropped, tipLabels()) # `dropped` might be outdated
      if (length(without)) {
      } else {
      }
      cons <- ConsensusWithout(r$trees, without, p = consP())
      cons <- UserRoot(cons)
      
      if (unitEdge()) {
        cons$edge.length <- rep.int(1, dim(cons$edge)[1])
      }
      cons <- SortEdges(cons)
      
      r$plottedTree <- cons
      plot(r$plottedTree, tip.color = TipCols()[intersect(cons$tip.label, kept)])
      LabelConcordance()
    }
  }
  
  LogConsensusPlot <- function() {
    BeginLogP()
    LogPar()
    dropped <- DroppedTips()
    
    if (length(dropped) &&
        length(input$excludedTip) &&
        nchar(input$excludedTip) &&
        input$excludedTip %in% tipLabels()) {
      
      LogCommentP("Prepare reduced consensus tree", 1)
      if (length(setdiff(dropped, input$excludedTip))) {
        LogCodeP(paste0("exclude <- ",
                       EnC(setdiff(dropped, input$excludedTip))))
        LogCodeP("consTrees <- lapply(trees, DropTip, exclude)")
        LogCodeP("labels <- setdiff(consTrees[[1]]$tip.label, exclude)")
      } else {
        LogCodeP("consTrees <- trees",
                "labels <- consTrees[[1]]$tip.label")
      }
      
      LogCommentP(paste0(
        "Colour tip labels according to their original 'instability' ",
        "(Smith 2022)")
      )
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)",
        paste0(
          "tipCols <- tipCols[setdiff(labels, ",
          Enquote(input$excludedTip), ")]"
        )
      )
      LogCommentP(paste0(
        "Plot the reduced consensus tree, showing position of ",
        gsub("_", " ", input$excludedTip, fixed = TRUE))
      )
      LogCodeP("plotted <- RoguePlot(",
              "  trees = consTrees,",
              paste0("  tip = ", Enquote(input$excludedTip), ","),
              paste0("  p = ", consP(), ","),
              "  edgeLength = 1,",
              if(length(r$outgroup)) {
                  paste0("  outgroupTips = ", EnC(r$outgroup), ",")
              },
              "  tip.color = tipCols",
              ")")
      
      LogCommentP("Store tree to plot concordance")
      LogCodeP("plottedTree <- plotted$cons")
      
      LogConcordance()
    } else {
      without <- intersect(dropped, tipLabels()) # `dropped` might be outdated
      LogCommentP("Calculate consensus tree")
      if (length(without)) {
        LogCodeP(
          "cons <- ConsensusWithout(",
          "  trees,",
          paste0("  ", EnC(without), ","),
          paste0("  p = ", consP()),
          ")")
      } else {
        LogCodeP(paste0(
          "cons <- Consensus(trees, p = ", consP(), ")"
        ))
      }
      LogUserRoot(dropped = without)
      if (unitEdge()) {
        LogCodeP("cons$edge.length <- rep.int(1L, nrow(cons$edge))")
      }
      LogSortEdges("cons")
      LogCommentP("Plot consensus tree")
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)[cons$tip.label]",
        "plot(cons, tip.color = tipCols)")
      LogConcordance("cons")
    }
  }
  
  PolEscVal <- reactive({
    LengthAdded(r$trees,
                r$dataset[tipLabels(), PlottedChar()],
                concavity())
  })
  
  CharacterwisePlot <- function() {
    par(mar = rep(0, 4), cex = 0.9)
    n <- PlottedChar()
    if (whichTree() > 0) {
      LogMsg("Plotting PlottedTree(", whichTree(), ", ", n, ")")
    }
    r$plottedTree <- PlottedTree()
    if (length(n) && n > 0L) {
      pc <- tryCatch({
        extraLen <- PolEscVal()
        roguishness <- if (max(extraLen) == 0) {
          "black"
        } else {
          hcl.colors(256, "inferno")[
            (192 * extraLen[r$plottedTree$tip.label] / max(extraLen)) + 1
          ]
        }
        PlotCharacter(
          if (whichTree() > 0) r$plottedTree else lapply(r$trees, UserRoot),
          r$dataset,
          n,
          edge.width = 2.5,
          updateTips = "updateTips" %in% input$mapDisplay,
          tip.color = roguishness,
          Display = function(tr) {
            tr <- UserRoot(tr)
            if (unitEdge()) {
              tr$edge.length <- rep.int(1, dim(tr$edge)[[1]])
            }
            SortEdges(tr)
          }
        )
        if (max(extraLen) > 0) {
          PlotTools::SpectrumLegend(
            "bottomleft", bty = "n",
            palette = hcl.colors(256, "inferno")[1:193],
            title = "Mean tree score\nimpact",
            title.font = 2,
            y.intersp = 1.42,
            legend = c(signif(4:1 * max(extraLen) / 4, 3), "No impact")
          )
        }
      },
      error = function (cond) {
        message("ASHOGAIH")
        message(cond)
        cli::cli_alert_danger(cond)
        Notification(type = "error",
                     "Could not match dataset to taxa in trees")
        ErrorPlot("Load dataset with\n", "character codings\n",
                  "for taxa on tree")
        return()
      }
      )
      
      LabelConcordance()
    } else {
      plot(r$plottedTree, tip.color = TipCols()[r$plottedTree$tip.label])
      TipColLegend()
    }
  }
  
  LogPar <- function() {
    LogCommentP("Set up plotting area")
    LogCodeP(c(
      "par(",
      "  mar = c(0, 0, 0, 0), # Zero margins",
      "  cex = 0.9            # Smaller font size",
      ")"
    ))
  }
  
  LogCharacterwisePlot <- function() {
    BeginLogP()
    LogPar()
    n <- PlottedChar()
    if (whichTree() > 0) {
      LogComment(paste("Select tree", whichTree(), "from tree set"))
    }
    LogPlottedTree()
    if (length(n) && n > 0L) {
      if (whichTree() > 0) {
        LogCommentP(paste("Map character", n, "onto tree", whichTree()))
      } else {
        LogCommentP(paste("Map character", n, "onto consensus tree"))
      }
      LogCodeP(
        "PlotCharacter(",
        if (whichTree() > 0) "  tree = plottedTree," else 
          paste0("  tree = RootTree(trees, ", EnC(r$outgroup), "),"),
        "  dataset = dataset,",
        paste0("  char = ", n, ","),
        paste0("  updateTips = ", "updateTips" %in% input$mapDisplay, ","),
        "  Display = function(tr) {",
        paste0("    tr <- RootTree(tr, ", EnC(r$outgroup), ")"),
        "    tr$edge.length <- rep.int(2, nrow(tr$edge))",
        "    SortTree(tr)",
        "  },",
        "  edge.width = 2.5",
        ")"
      )
      LogConcordance()
    } else {
      LogCommentP("Plot single tree")
      LogCodeP(
        "tipCols <- Rogue::ColByStability(trees)[plottedTree$tip.label]",
        "plot(plottedTree, tip.color = tipCols)"
      )
    }
  }
  
  MainPlot <- function() {
    if (AnyTrees()) {
      LogMsg("MainPlot()")
      switch(
        input$plotFormat,
        "cons" = {
          ConsensusPlot()
        },
        "clus" = {
          PlotClusterCons()
        },
        "ind" = {
          CharacterwisePlot()
        },
        "space" = {
          TreespacePlot()
        }
      ) # end switch
    }
  }
  ReactiveMainPlot <- reactive({MainPlot()})
  
  output$treePlot <- renderCachedPlot(
    ReactiveMainPlot(),
    cacheKeyExpr = { # Must be identical to RCode below
      switch(
        input$plotFormat,
        
        "clus" = list(r$treeHash, input$plotFormat,
                      r$keepNTips, input$excludedTip,
                      consP(),
                      input$neverDrop, r$outgroup,
                      input$distMeth,
                      input$concordance,
                      silThreshold(),
                      input$consP, input$concordance),
        "cons" = list(r$treeHash, input$plotFormat,
                      r$keepNTips, input$excludedTip,
                      consP(),
                      input$neverDrop, r$outgroup,
                      input$concordance),
        "ind" = list(PlottedChar(),
                     whichTree(),
                     input$concordance,
                     r$outgroup,
                     concavity(),
                     input$mapDisplay,
                     r$dataHash, r$treeHash), 
        "space" = list(r$treeHash, input$plotFormat,
                       min(dims(), nProjDim()),
                       TreeCols(),
                       treePch(),
                       input$distMeth,
                       input$spaceCol,
                       input$mapLines,
                       concavity(),
                       input$spacePch,
                       if (input$spacePch == "relat") input$relators,
                       silThreshold(),
                       input$display)
      )
    },
    sizePolicy = function(x) rep(input$plotSize, 2)
  )
  
  RCode <- bindCache(reactive({
    switch(
      input$plotFormat,
      "cons" = {
        LogConsensusPlot()
      },
      "clus" = {
        LogPlotClusterCons()
      },
      "ind" = {
        LogCharacterwisePlot()
      },
      "space" = {
        LogTreespacePlot()
      }
    )
    
    # Return:
    r$plotLog
  }),  # Must be identical to output$treePlot above
    switch(
      input$plotFormat,
      
      "clus" = list(r$treeHash, input$plotFormat,
                    r$keepNTips, input$excludedTip,
                    consP(),
                    input$neverDrop, r$outgroup,
                    input$distMeth,
                    input$concordance,
                    silThreshold(),
                    input$consP, input$concordance),
      "cons" = list(r$treeHash, input$plotFormat,
                    r$keepNTips, input$excludedTip,
                    consP(),
                    input$neverDrop, r$outgroup,
                    input$concordance),
      "ind" = list(PlottedChar(),
                   whichTree(),
                   input$concordance,
                   r$outgroup,
                   input$mapDisplay,
                   r$dataHash, r$treeHash), 
      "space" = list(r$treeHash, input$plotFormat,
                     min(dims(), nProjDim()),
                     TreeCols(),
                     treePch(),
                     input$distMeth,
                     input$spaceCol,
                     input$mapLines,
                     concavity(),
                     input$spacePch,
                     if (input$spacePch == "relat") input$relators,
                     silThreshold(),
                     input$display)
    )
  )
  
  UCFirst <- function (str) {
    paste0(toupper(substr(str, 1, 1)),
           substr(str, 2, nchar(str)))
  }
  
  nonAmbigContrast <- reactive({
    cont <- attr(r$dataset, "contrast")
    applic <- cont[, setdiff(colnames(cont), "-")]
    cont[rowSums(applic) == dim(applic)[[2]], ] <- 0
    
    # Return:
    cont
  })
  
  plottedTokens <- reactive({
    n <- PlottedChar()
    # `phyDat[,]` returns a new phyDat object with a different "contrast"
    # Hence we manually extract the compressed character tokens:
    phyColumn <- vapply(r$dataset, `[[`, integer(1),
                        attr(r$dataset, "index")[[n]], USE.NAMES = FALSE)
    tokens <- colSums(nonAmbigContrast()[phyColumn, ]) > 0L
    names(tokens[tokens])
  })
  
  output$charMapLegend <- bindCache(
    renderUI({
      n <- PlottedChar()
      if (length(n) && n > 0L && !is.null(r$chars)) {
      
        pal <- c("#00bfc6", "#ffd46f", "#ffbcc5", "#c8a500",
                 "#ffcaf5", "#d5fb8d", "#e082b4", "#25ffd3",
                 "#a6aaff", "#e6f3cc", "#67c4ff", "#9ba75c",
                 "#60b17f")
        
        states <- attr(r$chars, "state.labels")[[n]]
        tokens <- plottedTokens()
        appTokens <- setdiff(tokens, "-")
        .State <- function (glyph, text = "Error?", col = "red") {
          if (is.numeric(glyph)) {
            if (glyph > length(appTokens)) {
              return(NULL)
            }
            nonBlank <- states != ""
            text <- states[nonBlank][glyph]
            col <- pal[glyph]
            glyph <- appTokens[glyph]
          }
          
          tags$li(style = "margin-bottom: 2px;",
                  tags$span(glyph,
                            style = paste("display: inline-block;",
                                          "border: 1px solid;",
                                          "width: 1em;",
                                          "text-align: center;",
                                          "line-height: 1em;",
                                          "margin-right: 0.5em;",
                                          "background-color:", col, ";")
                  ),
                  tags$span(UCFirst(text)))
        }
        
        tagList(
          tags$h3(colnames(r$chars)[n]),
          tags$ul(style = "list-style: none;",
                  .State(1), .State(2), .State(3), .State(4), .State(5),
                  .State(6), .State(7), .State(8), .State(9),
                  .State(10), .State(11), .State(12), .State(13),
                  if ("-" %in% tokens) 
                    .State("-", "Inapplicable", "lightgrey"),
                  .State("?", "Ambiguous", "grey")
          )
        )
      }
    }),
    PlottedChar(),
    r$chars,
    r$dataset
  )
    
  
  output$charNotes <- bindCache(
    renderUI({
      n <- PlottedChar()
      if (length(n) && n > 0L
          && is.list(r$charNotes) && is.list(r$charNotes[[1]])
          && length(r$charNotes) >= n) {
      
        charNotes <- r$charNotes[[n]]
        description <- charNotes[[1]]
        notes <- charNotes[[2]]
        states <- attr(r$chars, "state.labels")[[n]]
        tokens <- plottedTokens()
        
        tagList(
          if (length(description) > 0) {
            tags$div(id = "char-description",
                     lapply(strsplit(description, "\n")[[1]], tags$p))
          },
          if (!is.null(notes)) tags$ul(class = "state-notes", {
            PrintNote <- function(note) {
              taxa <- names(note)[note]
              tags$li(class = "state-note",
                      tags$span(class = "state-note-label",
                                paste(gsub("_", " ", fixed = TRUE,
                                           taxa), collapse = ", ")),
                      tags$span(class = "state-note-detail",
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
            } else {
              notes <- notes[order(names(notes))]
              duplicates <- DuplicateOf(toupper(notes))
              apply(duplicates, 2, PrintNote)
            }
          }),
          if (!states[[1]] %in% c("", "''")
              && any(tokens == "-")) {
            tags$p(tags$em("Brazeau et al. (2019) advise that neomorphic (0/1) characters should not contain inapplicable tokens (-)."))
          }
        )
      }
    }),
    PlottedChar(),
    r$dataset,
    r$chars,
    r$charNotes
  )
  
  LogScore <- function (x) {
    (-(log10(1 - pmin(1, x) + 1e-2))) / 2
  }
  
  QualityPlot <- function (quality) {
    par(mar = c(2, 0, 0, 0))
    nStop <- length(badToGood) + 1L
    
    # LogMsg("QualityPlot()")
    plot(NULL, xlim = c(0, 1), ylim = c(-1.5, 2.5),
         ann = FALSE, axes = FALSE)
    x <- seq.int(from = 0, to = 1, length.out = nStop)
    segments(x[-nStop], numeric(nStop), x[-1], lwd = 5, col = badToGood)
    
    trust <- quality[["Trustworthiness"]]
    cont <- quality[["Continuity"]]
    txc <- quality[["sqrtTxC"]]
    
    if (trust > 1) {
      LogMsg("Preternaturally high Trustworthiness: ", trust)
    }
    if (cont > 1) {
      LogMsg("Preternaturally high Continuity: ", cont)
    }
    LogMsg(trust * nStop)
    segments(LogScore(txc), -1, y1 = 1, lty = 3)
    text(LogScore(trust), 1, "T", col = badToGood[LogScore(trust) * nStop])
    text(LogScore(cont), -1, "C", col = badToGood[LogScore(cont) * nStop])
    
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    ticks <- LogScore(tickPos)
    
    axis(1, at = ticks, labels = NA, line = 0)
    axis(1, tick = FALSE, at = ticks, labels = tickPos, line = 0)
    axis(1, line = -1, tick = FALSE,
         at = ticks[-1] - ((ticks[-1] - ticks[-length(ticks)]) / 2),
         labels = c("", "dire", "", "ok", "gd", "excellent"))
    axis(3, at = 0.5, tick = FALSE, line = -2, 
         paste0(dims(), "D mapping quality (trustw. / contin.):"))
  }
  
  output$pcQuality <- renderCachedPlot({
    if (length(r$trees) < 3) {
      return()
    }
    dstnc <- distances()
    mppng <- mapping()
    mppng <- mapping()[, seq_len(min(dim(mppng)[2], dims()))]
    neighbs <- min(10L, length(r$trees) / 2)
    future_promise(
      TreeDist::MappingQuality(dstnc, dist(mppng), neighbs),
      seed = NULL) %...>% QualityPlot
  }, cacheKeyExpr = {
    list(r$treeHash, input$distMeth, dims())
  },
    sizePolicy = function (dims) dims
  )
  
  
  output$howManyDims <- renderPlot({
    par(mar = c(2.5, 2.5, 0, 0), xpd = NA, mgp = c(1.5, 0.5, 0))
    txc <- projQual()["TxC", ]
    nStop <- length(badToGood)
    
    plot(txc, type = "n", ylim = c(min(txc, 0.5), 1),
         frame.plot = FALSE, axes = FALSE,
         xlab = "Dimension", ylab = "Trustw. \uD7 Contin.")
    par(xpd = FALSE)
    axis(1, 1:14)
    axis(2)
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    mids <- c(0.6, 0.75, 0.85, 0.925)
    text(rep.int(15, 4), mids, pos = 2, cex = 0.8,
         col = badToGood[nStop * LogScore(mids)],
         c("Essentially random", "Dangerous", "Usable", "Good"))
    text(1, 0.975, pos = 4, "Excellent", cex = 0.8, 
         col = badToGood[LogScore(0.975) * nStop])
    for (i in tickPos[-1]) {
      abline(h = i, lty = 3, col = badToGood[LogScore(i) * nStop])
    }
    points(txc, type = "b")
    txcNow <- txc[dims()]
    
    points(dims(), txcNow, pch = 16, col = badToGood[LogScore(txcNow) * nStop],
           cex = 1.6)
  })
  
  observeEvent(input$clThresh, {
    classes <- c("meaningless", "weak", "good", "strong")
    liveClass <- classes[as.integer(cut(input$clThresh, c(0, 0.25, 0.5, 0.7, 1),
                                        include.lowest = TRUE, right = FALSE))]
    addClass("clThresh-label", liveClass)
    removeClass("clThresh-label", setdiff(classes, liveClass))
  })
  silThreshold <- debounce(reactive({
    input$clThresh
  }), 50)
  
  ##############################################################################
  # Clusterings
  ##############################################################################
  clusterings <- bindCache(reactive({
    ## CAUTION: Update LogClusterings() to reflect any changes made
    ## to this function 
    LogMsg("clusterings()")
    maxCluster <- min(15L, length(r$trees) - 1L)
    if (maxCluster > 1L) {
      possibleClusters <- 2:maxCluster
      
      hSil <- pamSil <- -99
      dists <- distances()
      
      nMethodsChecked <- 3L
      cli::cli_progress_bar("Computing clusterings", "K-means",
                            total = nMethodsChecked)
      
      nK <- length(possibleClusters)
    
      kClusters <- lapply(possibleClusters,
                          function (k) TreeDist::KMeansPP(dists, k))
      kSils <- vapply(kClusters, function (kCluster) {
        mean(cluster::silhouette(kCluster$cluster, dists)[, 3])
      }, double(1))
      bestK <- which.max(kSils)
      kSil <- kSils[bestK]
      kCluster <- kClusters[[bestK]]$cluster
      
      cli::cli_progress_update(1, status = "PAM")
      pamClusters <- lapply(possibleClusters, function (k) {
        cluster::pam(dists, k = k)
      })
      pamSils <- vapply(pamClusters, function (pamCluster) {
        mean(cluster::silhouette(pamCluster)[, 3])
      }, double(1))
      bestPam <- which.max(pamSils)
      pamSil <- pamSils[bestPam]
      pamCluster <- pamClusters[[bestPam]]$cluster
      
      cli::cli_progress_update(1, status = "Hierarchical")
      hTree <- protoclust::protoclust(dists)
      hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
      hSils <- vapply(hClusters, function (hCluster) {
        mean(cluster::silhouette(hCluster, dists)[, 3])
      }, double(1))
      bestH <- which.max(hSils)
      hSil <- hSils[bestH]
      hCluster <- hClusters[[bestH]]
      cli::cli_progress_update(1, status = "Done")
      
      bestCluster <- c("none", "pam", "hmm", "kmn")[
        which.max(c(silThreshold(), pamSil, hSil, kSil))]
    } else {
      bestCluster <- "none"
    }
     
    LogMsg("Best clustering: ", bestCluster, 
        "; sil: ", signif(switch(bestCluster, pam = pamSil, hmm = hSil, kmn = kSil, 0)))
    # Return:
    list(method = switch(bestCluster, pam = "part. around medoids",
                                      hmm = "minimax linkage",
                                      kmn = "k-means",
                                      none = "no significant clustering"),
         n = 1 + switch(bestCluster, pam = bestPam, hmm = bestH, kmn = bestK, 0),
         sil = switch(bestCluster, pam = pamSil, hmm = hSil, kmn = kSil, 0), 
         cluster = switch(bestCluster, pam = pamCluster, hmm = hCluster, kmn = kCluster, 1)
    )

  }), r$treeHash, silThreshold(), input$distMeth)
  
  LogClusterings <- function() {
    maxCluster <- min(15L, length(r$trees) - 1L)
    if (maxCluster > 1L) {
      possibleClusters <- paste(2, maxCluster, sep = ":")
      
      hSil <- pamSil <- -99
      LogDistances()
      dists <- distances()
      
      LogCommentP("Compute clusters of trees", 2)
      nK <- length(possibleClusters)
      LogCommentP("Try K-means++ clustering (Arthur & Vassilvitskii 2007):")
      LogCodeP(
        paste0(
          "kClusters <- lapply(", possibleClusters, ", ",
          "function (k) KMeansPP(dists, k)", ")"
        ),
        "kSils <- vapply(kClusters, function (kCluster) {",
        "  mean(cluster::silhouette(kCluster$cluster, dists)[, 3])",
        "}, double(1))",
        "bestK <- which.max(kSils)",
        "kSil <- kSils[bestK] # Best silhouette coefficient",
        "kCluster <- kClusters[[bestK]]$cluster # Best solution"
      )
      
      LogCommentP("Try partitioning around medoids (Maechler et al. 2019):")
      LogCodeP(
        paste0(
          "pamClusters <- lapply(", possibleClusters, ", ",
          "function (k) cluster::pam(dists, k = k)", ")"
        ),
        "pamSils <- vapply(pamClusters, function (pamCluster) {",
        "  mean(cluster::silhouette(pamCluster)[, 3])",
        "}, double(1))",
        "bestPam <- which.max(pamSils)",
        "pamSil <- pamSils[bestPam] # Best silhouette coefficient",
        "pamCluster <- pamClusters[[bestPam]]$cluster # Best solution"
      )
      
      
      LogCommentP(
        paste("Try hierarchical clustering with minimax linkage",
              "(Bien & Tibshirani 2011):")
      )
      LogCodeP(
        "hTree <- protoclust::protoclust(dists)",
        paste0(
          "hClusters <- lapply(", possibleClusters, ", ", 
          "function (k) cutree(hTree, k = k)", ")"
        ),
        "hSils <- vapply(hClusters, function (hCluster) {",
        "  mean(cluster::silhouette(hCluster, dists)[, 3])",
        "}, double(1))",
        "bestH <- which.max(hSils)",
        "hSil <- hSils[bestH] # Best silhouette coefficient",
        "hCluster <- hClusters[[bestH]] # Best solution"
      )
      
      LogCommentP("Set threshold for recognizing meaningful clustering")
      LogCommentP("no support < 0.25 < weak < 0.5 < good < 0.7 < strong", 0)
      LogCodeP(paste0("threshold <- ", silThreshold()))
      
      LogCommentP("Compare silhouette coefficients of each method")
      LogCodeP(
        "bestMethodId <- which.max(c(threshold, pamSil, hSil, kSil))",
        "bestCluster <- c(\"none\", \"pam\", \"hmm\", \"kmn\")[bestMethodId]"
      )
      if (clusterings()$n == 1) {
        LogCommentP("No significant clustering was found.")
        LogCodeP("clustering <- 1 # Assign all trees to single cluster")
      } else {
        LogCommentP(paste0("Best clustering was ", clusterings()$method, ":"))
        LogCommentP(paste0("Silhouette coefficient = ",
                          signif(clusterings()$sil)), 0)
        LogCommentP(paste0("Store the cluster to which each tree is ",
                          "optimally assigned:"))
        LogCodeP(paste0(
          "clustering <- switch(bestCluster, pam = pamCluster, hmm = hCluster,",
          " kmn = kCluster, 1)"),
          paste0("nClusters <- length(unique(clustering))"),
          paste0(
          "clusterCol <- ",
          EnC(palettes[[min(length(palettes), clusterings()$n)]]),
          " # Arbitrarily"
          )
        )
      }
    } else {
      LogCommentP("Not enough trees for clustering analysis")
      LogCodeP("bestCluster <- \"none\"")
      LogCodeP("nClusters <- 1")
    }
  }
  
  PlotClusterCons <- function() {
    LogMsg("PlotClusterCons()")
    on.exit(LogMsg("/PlotClusterCons()"))
    
    cl <- clusterings()
    
    kept <- KeptTips()
    dropped <- if (length(kept) > 1) {
      setdiff(TipLabels(r$trees[[1]]), kept)
    } else {
      character(0)
    }
    par(mar = c(0.2, 0, 0.2, 0), xpd = NA)
    if (cl$sil > silThreshold()) {
      nRow <- ceiling(cl$n / 3)
      r$plottedTree <- vector("list", cl$n)
      par(mfrow = c(nRow, ceiling(cl$n / nRow)))

      for (i in seq_len(cl$n)) {
        col <- palettes[[min(length(palettes), cl$n)]][i]
        PutTree(r$trees)
        PutData(cl$cluster)
        
        cons <- ConsensusWithout(r$trees[cl$cluster == i], dropped, p = consP())
        cons <- UserRoot(cons)
        if (unitEdge()) {
          cons$edge.length <- rep.int(1, dim(cons$edge)[1])
        }
        cons <- SortEdges(cons)
        r$plottedTree[[i]] <- cons
        plot(cons, edge.width = 2, font = 3, cex = 0.83,
             edge.color = col, tip.color = TipCols()[cons$tip.label])
        legend("topright", paste0("Cluster ", i), pch = 15, col = col,
               pt.cex = 1.5, bty = "n")
        LabelConcordance()
      }
    } else {
      PutTree(r$trees)
      cons <- ConsensusWithout(r$trees, dropped, p = consP())
      cons <- UserRoot(cons)
      if (unitEdge()) {
        cons$edge.length <- rep.int(1, dim(cons$edge)[1])
      }
      cons <- SortEdges(cons)
      r$plottedTree <- cons
      plot(cons, edge.width = 2, font = 3, cex = 0.83,
           edge.color = palettes[[1]], tip.color = TipCols()[cons$tip.label])
      LabelConcordance()
      legend("topright", "No clustering", pch = 16, col = palettes[[1]],
             bty = "n")
    }
  }
  
  LogPlotClusterCons <- function() {
    LogMsg("PlotClusterCons()")
    on.exit(LogMsg("/PlotClusterCons()"))
    
    BeginLogP()
    
    cl <- clusterings()
    LogClusterings()
    
    kept <- KeptTips()
    dropped <- if (length(kept) > 1) {
      setdiff(TipLabels(r$trees[[1]]), kept)
    } else {
      character(0)
    }
    if (cl$sil > silThreshold()) {
      nRow <- ceiling(cl$n / 3)
      LogCommentP("Plot consensus of each tree cluster", 2)
      LogCodeP(paste0(
        "par(mfrow = c(", nRow, ", ",
        ceiling(cl$n / nRow), "))",
        " # Plotting area layout"
      ))
      LogCodeP(
        paste0(
          "tipCols <- Rogue::ColByStability(trees)", 
          " # Colour tips by stability"
        )
      )
      LogCommentP("Plot each consensus tree in turn:", 1)
      LogCodeP(paste0("for (i in seq_len(", cl$n, ")) {"))
      LogIndent(+2)
      LogCodeP(
        "clusterTrees <- trees[clustering == i]",
        "cons <- ConsensusWithout(", 
        "  trees = clusterTrees,",
        paste0("  tip = ", EnC(dropped), ","),
        paste0("  p = ", consP()),
        ")"
      )
      LogUserRoot(dropped = dropped)
      if (unitEdge()) {
        LogExprP("cons$edge.length <- rep.int(1, nrow(cons$edge))")
      }
      LogSortEdges("cons")
      LogCodeP("plot(",
              "  cons,",
              "  edge.width = 2,             # Widen lines",
              "  font = 3,                   # Italicize labels",
              "  cex = 0.83,                 # Shrink tip font size",
              "  edge.color = clusterCol[i], # Colour tree",
              "  tip.color = tipCols[cons$tip.label]",
              ")")
      LogCodeP("legend(", 
              "  \"bottomright\",",
              "  paste(\"Cluster\", i),",
              "  pch = 15,            # Filled circle icon",
              "  pt.cex = 1.5,        # Increase icon size",
              "  col = clusterCol[i],",
              "  bty = \"n\"            # Don't plot legend in box",
              ")")
      LogConcordance("cons")
      LogIndent(-2)
      LogCodeP("}")
    } else {
      LogCommentP("No clustering structure: Plot consensus tree")
      LogCodeP(
        if (length(dropped)) {
          c("cons <- ConsensusWithout(",
            "  trees = trees,",
            paste0("  tip = ", EnC(dropped), ","),
            paste0("  p = ", consP()),
            ")"
          )
        } else {
          paste0("cons <- Consensus(trees, p = ", consP(), ")")
        }
      )
      LogUserRoot("cons", dropped = dropped)
      if (unitEdge()) {
        LogCommentP("Set unit edge length", 0)
        LogCodeP("cons$edge.length <- rep.int(1, nrow(cons$edge))")
      }
      LogSortEdges("cons")
      LogCodeP("plottedTree <- cons # Store for future reference")
      
      LogCodeP("tipCols <- Rogue::ColByStability(trees)[cons$tip.label]")
      LogCommentP("Plot consensus tree")
      LogCodeP(
        "plot(",
        "  cons,",
        "  edge.width = 2, # Widen lines",
        "  font = 3,       # Italicize labels",
        "  cex = 0.83,     # Shrink tip font size",
        "  tip.color = tipCols",
        ")"
      )
      LogConcordance()
    }
  }
  
  ##############################################################################
  # Plot settings: point style
  ##############################################################################

  spaceCex <- reactive(1.7)
  spaceLwd <- reactive(2)
  
  FirstHit <- reactive({
    r$trees <- WhenFirstHit(r$trees)
    
    # Return:
    attr(r$trees, "firstHit")
  })
  
  LogFirstHit <- function() {
    LogCodeP("whenHit <- gsub(\"(seed|start|ratch\\\\d+|final)_\\\\d+\", \"\\\\1\",
              names(trees), perl = TRUE)")
    LogCodeP("attr(trees, \"firstHit\") <- table(whenHit)[unique(whenHit)]")
  }
  
  FirstHitCols <- reactive({
    if (is.null(FirstHit())) {
      palettes[[1]]
    } else {
      hcl.colors(length(FirstHit()), "viridis")
    }
  })
  
  LogFirstHitCols <- reactive({
    if (is.null(FirstHit())) {
      paste0(palettes[[1]], " # Arbitrarily")
    } else {
      "hcl.colors(length(firstHit), \"viridis\")"
    }
  })
  
  TreeCols <- reactive({
    switch(
      input$spaceCol,
      "clust" = {
        cl <- clusterings()
        if (cl$sil > silThreshold()) {
          palettes[[min(length(palettes), cl$n)]][cl$cluster]
        } else {
          palettes[[1]]
        }
      }, "score" = {
        if (is.null(scores()) || length(unique(scores())) == 1L) {
          palettes[[1]]
        } else {
          norm <- scores() - min(scores())
          norm <- (length(badToGood) - 1L) * norm / max(norm)
          rev(badToGood)[1 + norm]
        }
      }, "firstHit" = {
        if (is.null(FirstHit())) {
          Notification("Data not available; were trees loaded from file?",
                           type = "warning")
          palettes[[1]]
        } else {
          rep(FirstHitCols(), FirstHit())
        }
      },
      "black"
    )
  })
  
  LogTreeCols <- reactive({
    beige <- paste0("treeCols <- ", Enquote(palettes[[1]]), " # Arbitrarily")
    switch(
      input$spaceCol,
      "clust" = {
        cl <- clusterings()
        if (cl$sil > silThreshold()) {
          paste0("treeCols <- ", 
                 EnC(palettes[[min(length(palettes), cl$n)]]), 
                 "[clustering]")
        } else {
          beige
        }
      }, "score" = {
        if (is.null(scores()) || length(unique(scores())) == 1L) {
          beige
        } else {
          c(paste0("scores <- TreeLength(trees, dataset, concavity = ",
                   Enquote(concavity()), ")"),
            "normalized <- scores - min(scores)",
            "normalized <- 107 * normalized / max(normalized)",
            "goodToBad <- hcl.colors(108, \"Temps\")",
            "treeCols <- goodToBad[1 + normalized]"
          )
        }
      }, "firstHit" = {
        if (is.null(FirstHit())) {
          beige
        } else {
          c("trees <- WhenFirstHit(trees)",
            "firstHit <- attr(trees, \"firstHit\")",
            paste0("treeCols <- rep(", LogFirstHitCols(), ", firstHit))")
          )
        }
      },
      "treeCols <- black"
    )
  })
  
  treeNameClustering <- reactive({
    ClusterStrings(names(r$trees))
  })
 
  treePch <- reactive({
    switch(
      input$spacePch,
      "clust" = {
        cl <- clusterings()
        if (cl$sil > silThreshold()) {
          cl$cluster - 1
        } else {
          16 # Filled circle
        }
      }, "relat" = {
        quartet <- input$relators
        if (length(quartet) == 4) {
          QuartetResolution(r$trees, input$relators)
        } else {
          Notification("Select four taxa to show relationships")
          0
        }
      }, "name" = {
        if (is.null(names(r$trees))) {
          Notification("Trees lack names", type = "warning")
          16
        } else {
          indices <- treeNameClustering()
          # Match pch from BGS2019 Fig. 9 for pre-loaded datasets.
          # Embarrassingly, in BGS19 I plotted ambigAbsent instead of ambiguous.
          # Sadly, Systematic Biology will not allow a correction.
          c(1, 3, 4, 2, seq_len(max(indices))[-(1:4)])[indices]
        }
      }, 0)
  })
  
  LogTreePch <- function() {
    switch(
      input$spacePch,
      "clust" = {
        cl <- clusterings()
        if (cl$sil > silThreshold()) {
          "cl$cluster - 1"
        } else {
          "16 # No clustering structure: Use filled circle"
        }
      }, "relat" = {
        quartet <- input$relators
        if (length(quartet) == 4) {
          paste0("QuartetResolution(trees, ", EnC(input$relators), ")")
        } else {
          "0 # Square"
        }
      }, "name" = {
        if (is.null(names(r$trees))) {
          "16 # Filled circle"
        } else {
          "ClusterStrings(names(trees))"
        }
      }, "0 # Square")
  }
  
  maxProjDim <- reactive({
    min(12, max(0L, length(r$trees) - 1L))
  })
  
  nProjDim <- reactive({
    dim(mapping())[2]
  })
  
  dims <- debounce(reactive({
    if (mode3D()) 3L else {
      min(input$spaceDim, maxProjDim())
    }
  }), 400)
  
  Quartet <- function (...) {
    if (!requireNamespace("Quartet", quietly = TRUE)) {
      Notification("Installing required package \"Quartet\"",
                       type = "warning", duration = 20)
      install.packages("Quartet")
    }
    as.dist(Quartet::QuartetDivergence(
      Quartet::ManyToManyQuartetAgreement(...), similarity = FALSE))
  }
  
  distances <- bindCache(reactive({
    ## CAUTION: LogDistances() must be updated to reflect any changes to
    ## this code
    LogMsg("distances(): ", input$distMeth)
    if (length(r$trees) > 1L) {
      Dist <- switch(input$distMeth,
                     "cid" = TreeDist::ClusteringInfoDistance,
                     "pid" = TreeDist::PhylogeneticInfoDistance,
                     "msid" = TreeDist::MatchingSplitInfoDistance,
                     "rf" = TreeDist::RobinsonFoulds,
                     "qd" = Quartet)
      withProgress(
        message = "Initializing distances...", value = 0.99,
        Dist(r$trees)
      )
    } else {
      matrix(0, 0, 0)
    }
    
  }), input$distMeth, r$treeHash)
  
  LogDistances <- function() {
    LogCommentP("Compute tree distances")
    LogCodeP(switch(
        input$distMeth,
        "cid" = "dists <- TreeDist::ClusteringInfoDistance(trees)",
        "pid" = "dists <- TreeDist::PhylogeneticInfoDistance(trees)",
        "msid" = "dists <- TreeDist::MatchingSplitInfoDistance(trees)",
        "rf" = "dists <- TreeDist::RobinsonFoulds(trees)",
        "qd" = c("dists <- as.dist(Quartet::QuartetDivergence(",
          "  Quartet::ManyToManyQuartetAgreement(trees),",
          "  similarity = FALSE)", ")")
    ))
  }
  
  mapping <- bindCache(reactive({
    LogMsg("mapping()")
    if (maxProjDim() > 1L) {
      withProgress(
        message = "Mapping trees",
        value = 0.99,
        tryCatch(cmdscale(distances(), k = maxProjDim()),
                 warning = function (e) {
                   nDim <- as.integer(substr(e$message, 6, 7))
                   updateSliderInput(inputId = "spaceDim",
                                     value = min(nDim, input$spaceDim),
                                     max = nDim)
                   message("Max dimensions available for mapping: ", nDim, ".")
                   cmdscale(distances(), k = nDim)
                 })
      )
    } else {
      matrix(0, 0, 0)
    }
  }), r$treeHash, input$distMeth, maxProjDim())
  
  LogMapping <- function() {
    k <- dim(mapping())[2]
    if (!is.null(k) && k > 0) {
      LogCommentP(paste0(
        "Generate first ", k, " dimensions of tree space using PCoA"
      ))
      LogCodeP(paste0("map <- cmdscale(dists, k = ", k, ")"))
    }
  }
  
  mstEnds <- bindCache(reactive({
    dist <- as.matrix(distances())
    withProgress(message = "Calculating MST", {
      edges <- MSTEdges(dist)
    })
    edges
  }), input$distMeth, r$treeHash)
  
  ##############################################################################
  # Plot tree space
  ##############################################################################
  # CAUTION: Remember to update accompanying logging function below.
  TreespacePlot <- function() {
    if (length(r$trees) < 3) {
      return(ErrorPlot("Need at least\nthree trees to\nmap tree space"))
    }
    
    cl <- clusterings()
    map <- mapping()
    
    nDim <- min(dims(), nProjDim())
    if (nDim < 2) {
      if (dim(map)[2] == 1L) {
        map <- cbind(map, 0)
      } else {
        map[, 2] <- 0
      }
      nDim <- 2L
      nPanels <- 1L
    } else {
      plotSeq <- matrix(0, nDim, nDim)
      nPanels <- nDim * (nDim - 1L) / 2L
      plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)
      if (nDim > 2) {
        plotSeq[nDim - 1, 2] <- max(plotSeq) + 1L
      }
      layout(t(plotSeq[-nDim, -1]))
    }
    
    par(mar = rep(0.2, 4))
    withProgress(message = "Drawing plot", {
      for (i in 2:nDim) for (j in seq_len(i - 1)) {
        incProgress(1 / nPanels)
        # Set up blank plot
        plot(map[, j], map[, i], ann = FALSE, axes = FALSE,
             frame.plot = nDim > 2L,
             type = "n", asp = 1, xlim = range(map), ylim = range(map))
        
        # Connect sequential trees
        if ("seq" %in% input$mapLines) {
          lines(map[, j], map[, i], col = "#ffcc33", lty = 2)
        }
        
        # Plot MST
        if ("mst" %in% input$mapLines) {
          segments(map[mstEnds()[, 1], j], map[mstEnds()[, 1], i],
                   map[mstEnds()[, 2], j], map[mstEnds()[, 2], i],
                   col = "#bbbbbb", lty = 1)
        }
        
        
        # Add points
        points(map[, j], map[, i], pch = treePch(),
               col = paste0(TreeCols(), as.hexmode(200)),
               cex = spaceCex(),
               lwd = spaceLwd()
               )#input$pt.cex)
        
        if (cl$sil > silThreshold() && "hull" %in% input$mapLines) {
          # Mark clusters
          for (clI in seq_len(cl$n)) {
            inCluster <- cl$cluster == clI
            clusterX <- map[inCluster, j]
            clusterY <- map[inCluster, i]
            hull <- chull(clusterX, clusterY)
            polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
                    border = palettes[[min(length(palettes), cl$n)]][clI])
          }
        }
        if ("labelTrees" %in% input$display) {
          text(map[, j], map[, i], names(r$trees))
        }
      }
      if (nDim > 2) {
        plot.new()
      }
      if (input$spacePch == "relat") {
        if (length(input$relators) == 4L) {
          legend(
            "topright",
            bty = "n",
            pch = 1:3,
            xpd = NA,
            pt.cex = spaceCex(),
            pt.lwd = spaceLwd(),
            gsub("_", " ", fixed = TRUE,
                 paste(input$relators[2:4], "&", input$relators[[1]]))
          )
        }
      } else if (input$spacePch == "name") {
        clstr <- treeNameClustering()
        clusters <- unique(clstr)
        if (length(clusters) > 1L) {
          legend(bty = "n", "topright", xpd = NA,
                 pch = c(1, 3, 4, 2, seq_len(max(clstr))[-(1:4)])[clusters],
                 paste0("~ ", attr(clstr, "med"), " (", table(clstr), ")"))
        }
      }
      if (input$spaceCol == "firstHit" && length(FirstHit())) {
        legend(bty = "n", "topleft", pch = 16, col = FirstHitCols(),
               pt.cex = spaceCex(),
               names(FirstHit()), title = "Iteration first hit")
      } else if (input$spaceCol == "score") {
        legendRes <- length(badToGood)
        leg <- rep(NA, legendRes)
        leg[c(legendRes, 1)] <- signif(range(scores()))
        legend("bottomright", bty = "n", border = NA,
               legend = leg, fill = rev(badToGood),
               y.intersp = 0.04, cex = 1.1)
      }
    })
  }
  
  LogTreespacePlot <- function() {
    BeginLogP()
    
    LogClusterings()
    LogMapping()
    
    map <- mapping()
    nDim <- min(dims(), nProjDim())
    if (nDim < 2) {
      LogCommentP("Prepare 1D map", 0)
      if (dim(map)[2] == 1L) {
        LogCodeP("map <- cbind(map, 0)")
      } else {
        LogCodeP("map[, 2] <- 0")
      }
      nDim <- 2L
      nPanels <- 1L
    } else {
      LogCommentP("Prepare plot layout")
      
      LogCodeP(c(
        paste0("nDim <- ", nDim, " # Number of dimensions to plot"),
        "nPanels <- nDim * (nDim - 1L) / 2L # Lower-left triangle",
        "plotSeq <- matrix(0, nDim, nDim)",
        "plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)",
        if (nDim > 2) {
           "plotSeq[nDim - 1, 2] <- max(plotSeq) + 1L"
        },
        "layout(t(plotSeq[-nDim, -1]))"
      ))
    }
    
    LogCommentP("Set plot margins", 0)
    LogCodeP("par(mar = rep(0.2, 4))")
    
    LogCommentP("Set up tree plotting symbols")
    LogCodeP(paste0("treePch <- ", LogTreePch()),
             LogTreeCols(),
             "treeCols <- paste0(treeCols, as.hexmode(200)) # Semitransparent"
    )
    
    LogCodeP("for (i in 2:nDim) for (j in seq_len(i - 1)) {")
    LogIndent(+2)
    LogCommentP("Set up blank plot")
    LogCodeP("plot(",
            "  x = map[, j],",
            "  y = map[, i],",
            "  ann = FALSE,        # No annotations",
            "  axes = FALSE,       # No axes",
            paste0("  frame.plot = ", 
                   if(nDim > 2L) {
                     "TRUE,  # Border around plot"
                   } else {
                     "FALSE, # No border around plot"  
                   }),
            "  type = \"n\",         # Don't plot any points yet",
            "  asp = 1,            # Fix aspect ratio to avoid distortion",
            "  xlim = range(map),  # Constant X range for all dimensions",
            "  ylim = range(map)   # Constant Y range for all dimensions",
            ")")
    
    if ("seq" %in% input$mapLines) {
      LogCommentP("Connect trees in sequence")
      LogCodeP("lines(",
               "  x = map[, j],",
               "  y = map[, i],",
               "  col = \"#ffcc33\", # Orange",
               "  lty = 2 # dashed",
               ")")
    }
    
    if ("mst" %in% input$mapLines) {
      LogCommentP("Plot minimum spanning tree (Gower 1969)")
      LogCodeP(
        "mst <- MSTEdges(as.matrix(dists))",
        "segments(",
        "  x0 = map[mst[, 1], j],",
        "  y0 = map[mst[, 1], i],",
        "  x1 = map[mst[, 2], j],",
        "  y1 = map[mst[, 2], i],",
        "  col = \"#bbbbbb\", # Light grey",
        "  lty = 1          # Solid lines",
        ")"
      )
    }
    
    LogCommentP("Add points")
    LogCodeP(
      "points(",
      "  x = map[, j],",
      "  y = map[, i],",
      "  pch = treePch,",
      "  col = treeCols,",
      paste0("  cex = ", spaceCex(), ", # Point size"),
      paste0("  lwd = ", spaceLwd(), " # Line width"),
      ")"
    )
    
    cl <- clusterings()
    if (cl$sil > silThreshold() && "hull" %in% input$mapLines) {
      LogCommentP("Mark clusters")
      LogCodeP("for (clI in seq_len(nClusters)) {")
      LogIndent(+2)
      LogCodeP(
        "inCluster <- clustering == clI",
        "clusterX <- map[inCluster, j]",
        "clusterY <- map[inCluster, i]",
        "hull <- chull(clusterX, clusterY)",
        "polygon(",
        "  x = clusterX[hull],",
        "  y = clusterY[hull],",
        "  lty = 1, # Solid line style",
        "  lwd = 2, # Wider line width",
        "  border = clusterCol[clI]",
        ")")
      LogIndent(-2)
      LogCodeP("}")
    }
    if ("labelTrees" %in% input$display) {
      #TODO input$display doesn't exist. If useful, implement below too.
      LogCodeP("text(map[, j], map[, i], trees)")
    }
    
    LogIndent(-2)
    LogCodeP("}")
    
    if (nDim > 2) {
      LogCodeP("plot.new() # Use new panel to plot legends")
    }
    
    if (input$spacePch == "relat") {
      if (length(input$relators) == 4L) {
        LogCommentP("Add legend for plotting symbols")
        LogCodeP(
          "legend(",
          "  \"topright\",",
          "  bty = \"n\", # No legend border box",
          "  pch = 1:3, # Legend symbols",
          "  xpd = NA, # Display overflowing text",
          paste0("  pt.cex = ", spaceCex(), ", # Point size"),
          paste0("  pt.lwd = ", spaceLwd(), ", # Line width"),
          paste0("  ",
                 EnC(gsub("_", " ", fixed = TRUE,
                          paste(input$relators[2:4], "&", input$relators[[1]])))
          ), ")"
        )
      }
    } else if (input$spacePch == "name") {
      clstr <- treeNameClustering()
      clusters <- unique(clstr)
      if (length(clusters) > 1L) {
        LogCommentP("Add legend for plotting symbols")
        LogCodeP(
          "nameClusters <- ClusterStrings(names(trees))",
          "uniqueClusters <- unique(nameClusters)",
          "legend(",
          "  \"topright\",",
          "  bty = \"n\", # No legend border box",
          "  xpd = NA, # Display overflowing text",
          paste0(
            "  pch = ",
            EnC(c(1, 3, 4, 2, seq_len(max(clstr))[-(1:4)])[clusters]),
            ", # Legend symbols"
          ), paste0("  ", 
                    EnC(paste0("~ ", attr(clstr, "med"),
                                   " (", table(clstr), ")"))
          ),
          ")")
      }
    }
    if (input$spaceCol == "firstHit" && length(FirstHit())) {
      LogCommentP("Record when trees first hit")
      LogFirstHit()
      
      LogCommentP("Add legend for symbol colours")
      LogCodeP(
        "legend(",
        "  \"topleft\",",
        "  bty = \"n\", # No legend border box",
        "  pch = 16, # Circle symbol",
        "  xpd = NA, # Display overflowing text",
        paste0("  col = ", LogFirstHitCols(), ","),
        paste0("  pt.cex = ", spaceCex(), ", # Point size"),
        paste0("  ", EnC(names(FirstHit())), ","),
        "  title = \"Iteration first hit\"",
        ")"
      )
    } else if (input$spaceCol == "score") {
      LogCommentP("Add legend for symbol colours")
      LogCodeP(
        "goodToBad <- hcl.colors(108, \"Temps\")",
        "leg <- rep_len(NA, 108)",
        paste0("leg[c(1, 108)] <- ", EnC(rev(signif(range(scores()))))),
        "legend(",
        "  \"bottomright\",",
        "  legend = leg,",
        "  bty = \"n\", # No legend border box",
        "  border = NA, # No border around plot icons",
        "  xpd = NA, # Display overflowing text",
        "  fill = goodToBad,",
        "  y.intersp = 0.04, # Compress squares to make gradient scale",
        "  cex = 1.1 # Increase font and icon size slightly",
        ")"
      )
    }
  }
  
  mode3D <- reactive("show3d" %in% input$display)
  
  saveDetails <- reactive({
    switch(input$plotFormat,
           "cons" = {
             list(
               fileName = "ConsensusTrees",
               title = "Consensus tree - TreeSearch",
               asp = 2L
             )
           },
           "clus" = {
             list(
               fileName = "ClusterCons",
               title = "Cluster Consensus trees - TreeSearch",
               asp = 1.6
             )
           },
           "ind" = {
             list(
               fileName = "OptimalTree",
               title = "Optimal tree - TreeSearch",
               asp = 2L
             )
           },
           "space" = {
             list(
               fileName = "TreeSpace",
               title = "Tree space - TreeSearch",
               asp = 1L
             )
           })
  })
  
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
  
  ##############################################################################
  # References
  ##############################################################################
  
  output$references <- renderUI({
    tagList(
     tags$h2("References for methods used"),
     tags$h3("Tree search"),
     HTML(Brazeau2019, Morphy, Nixon1999, SmithSearch),
     tags$h3("Tree space mapping"),
     HTML(paste0(Gower1966, Gower1969, Kaski2003, RCoreTeam,
                 SmithDist, Smith2020, SmithSpace, 
                 Venna2001)),
     tags$h3("Clustering"),
     HTML(paste("Cluster consensus trees:", Stockham2002)),
     HTML(paste0(
       "k-means++:", Arthur2007, Hartigan1979, 
       "Partitioning around medoids:", Maechler2019,
       "Hierarchical, minimax linkage:", Bien2011, Murtagh1983,
       "Clustering evaluation:", Rousseeuw1987
       )),
     tags$h3("Rogue taxa"),
     HTML(paste("Detection:", SmithRogue)),
     HTML(paste("Plotting:", Klopfstein2019)),
     HTML(paste("Character analysis:", Pol2009)),
    )
  })

  onStop(function() {
    options(startOpt)
    if (file.exists(cmdLogFile)) {
      unlink(cmdLogFile)
    }
    unlink(DataFileName("*"))
    unlink(TreeFileName("*"))
    if (logging) {
      LogMsg("Session has ended")
      on.exit(close(logMsgFile))
    }
  })
}


shinyApp(ui = ui, server = server)
