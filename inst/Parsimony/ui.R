fluidPage(
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
      data_ui_elems$data_source,
      data_ui_elems$data_file,
      data_ui_elems$readxl_options,
      se_ui$label,
      se_ui$config,
      se_ui$go,
      dl_ui$save_zip,
      data_ui_elems$tree_file,
      se_ui$results,
      hidden(tags$div(id = "manipulateTreeset",
        data_ui_elems$nTree_input,
        data_ui_elems$treeRange_input,
        tags$label("Save chosen trees:", class = "control-label"),
        tags$div(style = "display: inline-block",
          dl_ui$save_nwk,
          dl_ui$save_nex
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
               dl_ui$save_plot_zip,
               dl_ui$save_pdf,
               dl_ui$save_png
      ),
      tags$div(id = "savePlottedTrees",
               dl_ui$save_plot_nwk,
               dl_ui$save_plot_nex
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
          clustering_ui("clustering"),
          selectInput("distMeth", "Distance method:", selected = "cid",
                      choices = list("Clustering Information" = "cid",
                                     "Phylogenetic information" = "pid",
                                     "Matching split info" = "msid",
                                     "Robinson-Foulds (fast, iffy)" = "rf",
                                     "Quartet (slower)" = "qd"),
                      width = 200)
      )),
      hidden(treespace_ui("treespace")),
      references_ui("refs"),
    ),
  )
)
