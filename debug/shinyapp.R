
input <- list(
  dataSource = "Sun2018",
  mapDisplay = "",
  plotFormat = "ind"
)
dataFile <- system.file(paste0("datasets/", input$dataSource, ".nex"),
                        package = "TreeSearch")
r <- list(
  chars = ReadCharacters(dataFile),
  charNotes = ReadNotes(dataFile),
  dataset = ReadAsPhyDat(dataFile),
  trees = read.nexus(dataFile),
  null = NULL
)

conc <- function() 10
n=16

rogueCont <- PolEscapa(r$trees,
                       r$dataset[r$trees[[1]][["tip.label"]], n],
                       10)
roguishness <- if (max(rogueCont) == 0) {
  "black"
} else {
  hcl.colors(256, "inferno")[
    (192 * rogueCont[r$trees[[1]][["tip.label"]]] / max(rogueCont)) + 1
  ]
}
par(mar = rep(0, 4), cex=0.9)
r$plottedTree = r$trees[[2]]
PlotCharacter(r$plottedTree, r$dataset, n,
              edge.width = 2.5,
              updateTips = "updateTips" %in% input$mapDisplay,
              tip.color = roguishness)

distances <- ClusteringInfoDist(r$trees)
mapped <- cmdscale(distances, k = 6)


dataset <- source("inst/Parsimony/dataset.lg")$value
