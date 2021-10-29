file.pattern <- ".*?(\\w+\\d{4})\\.nex$"
files <- list.files('data-raw', file.pattern, full.names = TRUE)
inapplicable.datasets <- lapply(files, ape::read.nexus.data)
names(inapplicable.datasets) <- inapplicable.names <- gsub(file.pattern, "\\1", files)
inapplicable.phyData <- setNames(lapply(files, TreeTools::ReadAsPhyDat),
                                 inapplicable.names)
inapplicable.trees <- setNames(lapply(files, ape::read.nexus),
                               inapplicable.names)

usethis::use_data(inapplicable.datasets,
                  inapplicable.phyData,
                  inapplicable.trees,
                  overwrite = TRUE)
