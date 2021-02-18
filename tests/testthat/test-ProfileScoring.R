test_that("Profile scoring (TEMPORARY)", {
  data('congreveLamsdellMatrices')
  dataset <- congreveLamsdellMatrices[[42]]
  dataset <- PrepareDataProfile(dataset)
  tree <- NJTree(dataset)
  edge <- Preorder(tree)$edge
  at <- attributes(dataset)
  profiles <- attr(dataset, 'info.amounts')
  charSeq <- seq_along(dataset[[1]]) - 1L

  characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                            byTaxon = FALSE, concatenate = FALSE)
  startWeights <- at$weight
  morphyObjects <- lapply(characters, SingleCharMorphy)
  on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
          add = TRUE)
    
  expect_equal(ProfileScore(tree, dataset),
               morphy_profile(edge, morphyObjects, startWeights, 
                              charSeq, profiles, Inf))
  
  MaximizeParsimony(dataset, tree, concavity = 'profile')
})