input <- list(
  dataSource = "Agnarsson2004",
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

r$plottedTree = r$trees[[1]]

distances <- ClusteringInfoDist(r$trees)
mapped <- cmdscale(distances, k = 6)
