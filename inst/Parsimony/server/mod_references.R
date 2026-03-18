# Module: References panel
#
# Renders the references section. No reactive inputs or state reads/writes.
# Citation variables are defined in global.R and passed via `cites`.

references_ui <- function(id) {
  ns <- NS(id)
  htmlOutput(ns("references"), style = "clear: both;")
}

#' @param id Module namespace id.
#' @param cites Named list of citation HTML strings. Defaults to looking up
#'   each variable in the calling environment (i.e. global.R when run as app).
references_server <- function(id, cites = NULL) {
  # If no cites list supplied, collect from the caller's environment so the
  # app's global.R assignments are found automatically.
  if (is.null(cites)) {
    e <- parent.frame()
    get_cite <- function(nm) get(nm, envir = e, inherits = TRUE)
    cites <- list(
      Brazeau2019    = get_cite("Brazeau2019"),
      Morphy         = get_cite("Morphy"),
      Nixon1999      = get_cite("Nixon1999"),
      SmithSearch    = get_cite("SmithSearch"),
      Gower1966      = get_cite("Gower1966"),
      Gower1969      = get_cite("Gower1969"),
      Kaski2003      = get_cite("Kaski2003"),
      RCoreTeam      = get_cite("RCoreTeam"),
      SmithDist      = get_cite("SmithDist"),
      Smith2020      = get_cite("Smith2020"),
      SmithSpace     = get_cite("SmithSpace"),
      Venna2001      = get_cite("Venna2001"),
      Stockham2002   = get_cite("Stockham2002"),
      Arthur2007     = get_cite("Arthur2007"),
      Hartigan1979   = get_cite("Hartigan1979"),
      Maechler2019   = get_cite("Maechler2019"),
      Bien2011       = get_cite("Bien2011"),
      Murtagh1983    = get_cite("Murtagh1983"),
      Rousseeuw1987  = get_cite("Rousseeuw1987"),
      SmithRogue     = get_cite("SmithRogue"),
      Klopfstein2019 = get_cite("Klopfstein2019"),
      Pol2009        = get_cite("Pol2009")
    )
  }

  moduleServer(id, function(input, output, session) {
    output$references <- renderUI({
      tagList(
        tags$h2("References for methods used"),
        tags$h3("Tree search"),
        HTML(cites$Brazeau2019, cites$Morphy, cites$Nixon1999, cites$SmithSearch),
        tags$h3("Tree space mapping"),
        HTML(paste0(cites$Gower1966, cites$Gower1969, cites$Kaski2003,
                    cites$RCoreTeam, cites$SmithDist, cites$Smith2020,
                    cites$SmithSpace, cites$Venna2001)),
        tags$h3("Clustering"),
        HTML(paste("Cluster consensus trees:", cites$Stockham2002)),
        HTML(paste0(
          "k-means++:", cites$Arthur2007, cites$Hartigan1979,
          "Partitioning around medoids:", cites$Maechler2019,
          "Hierarchical, minimax linkage:", cites$Bien2011, cites$Murtagh1983,
          "Clustering evaluation:", cites$Rousseeuw1987
        )),
        tags$h3("Rogue taxa"),
        HTML(paste("Detection:", cites$SmithRogue)),
        HTML(paste("Plotting:", cites$Klopfstein2019)),
        HTML(paste("Character analysis:", cites$Pol2009)),
      )
    })
  })
}
