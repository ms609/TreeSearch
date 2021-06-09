test_that("PlotCharacter()", {
  
  Character <- function (str, plot = FALSE) {
    tree <- ape::read.tree(text = 
     "((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));")
    dataset <- TreeTools::StringToPhyDat(str, tips = tree)
    PlotCharacter(tree, dataset,
                  edge.width = 2, plot = plot)
  }
  
  expect_equal(structure(c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, 
                           TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, 
                           TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, 
                           FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, 
                           FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, 
                           TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                           TRUE, FALSE, FALSE, FALSE, TRUE, TRUE), .Dim = c(23L, 5L), .Dimnames = list(
                             NULL, c("-", "0", "1", "2", "3"))),
               Character("23--1??--032"))
})
