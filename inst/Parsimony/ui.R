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
      tags$div(
        style = "display: flex; align-items: center; gap: 8px; margin-top: 0.4em;",
        tags$span(
          style = "flex-shrink: 0; line-height: 0;",
          HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32" width="32" height="32" aria-label="TreeSearch logo">
            <circle cx="12" cy="12" r="9" fill="white" stroke="#111" stroke-width="2.5"/>
            <line x1="12" y1="19" x2="12" y2="15" stroke="#111" stroke-width="1.5" stroke-linecap="round"/>
            <line x1="12" y1="15" x2="7" y2="10" stroke="#111" stroke-width="1.5" stroke-linecap="round"/>
            <line x1="12" y1="15" x2="17" y2="10" stroke="#111" stroke-width="1.5" stroke-linecap="round"/>
            <line x1="17" y1="10" x2="15" y2="7" stroke="#111" stroke-width="1.2" stroke-linecap="round"/>
            <line x1="17" y1="10" x2="19" y2="7" stroke="#111" stroke-width="1.2" stroke-linecap="round"/>
            <line x1="18.5" y1="18.5" x2="28" y2="28" stroke="#111" stroke-width="2.5" stroke-linecap="round"/>
          </svg>')
        ),
        tags$h1("TreeSearch", style = "margin: 0;")
      ),
      data_ui_elems$data_source,
      data_ui_elems$data_file,
      data_ui_elems$readxl_options,
      se_ui$label,
      se_ui$config,
      se_ui$go,
      se_ui$cancel,
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
                   "cons"),
                 hidden(tags$div(id = "whichTree",
                   co_ui$which_tree
                 )),
                 hidden(tags$div(id = "treePlotConfig",
                   co_ui$tree_plot_config
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
      co_ui$tree_plot,
      hidden(tags$div(id = "charChooser", co_ui$char_chooser)),
      hidden(tags$div(id = "consConfig", co_ui$cons_config)),
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
