test_that("PlotCharacter()", {
  
  Character <- function (str, plot = FALSE, ...) {
    tree <- ape::read.tree(text = 
     "((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));")
    dataset <- TreeTools::StringToPhyDat(str, tips = tree)
    PlotCharacter(tree, dataset,
                  edge.width = 3, plot = plot, ...)
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
               Character("23--1??--032", updateTips = TRUE))
  
  skip_if_not_installed('vdiffr')
  
  Test <- if (interactive()) {
    function (str) invisible(Character(str, plot = TRUE))
  } else {
    function (str) {
      vdiffr::expect_doppelganger(
        paste0('PlotChar_',
               gsub('?', 'Q',
                    gsub('(', 'd',
                         gsub(')', 'b',
                              gsub('-', 'I', str,
                                   fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)),
        function() Character(str, plot = TRUE))
    }
  }
  Test("23--1??--032")
  Test("23--1??(-0)-(01)32")
  Test("23??1????032")
  Test("11--????--11")
  Test("000011????00")
  Test("????????????")
  Test("-------?????")
  Test("------------")
  Test("1234(45)AACGTTT")
  
  # From TGuillerme testing suite:
  Test("11-------100")
  Test("1100----1100")
  Test("000011110000")
  Test("1---1111---1")
  Test("----1111---1")
  Test("01----010101")
  Test("01---1010101")
  Test("1??--??--100")
  Test("21--3??--032")
  Test("11--1??--111")
  Test("11--1000001-")
  Test("01------0101")
  Test("110--?---100")
  Test("210--100--21")
  Test("????----1???")
  Test("23--1----032")
  Test("1----1----1-")
  Test("-1-1-1--1-1-")
  
  Test("--------0101")
  Test("10101-----01")
  Test("011--?--0011")
  Test("110--??--100")
  Test("21--1----012")
  Test("11----111111")
  Test("210210------")
  Test("----1111----")
  Test("230--??1--32")
  Test("023--??1--32")
  Test("023-???1--32")
  Test("23--1?1--023")
  Test("----1010----")
  Test("------11---1")
  Test("10----11---1")
  Test("320--??3--21")
})


test_that("Out-of-sequence works", {
  skip_if_not_installed('vdiffr')
  vdiffr::expect_doppelganger('PlotChar_out-of-sequence',
                              function () {
  PlotCharacter(ape::read.tree(text = '(a, (b, (c, d)));'),
                TreeTools::StringToPhyDat('1342', tips = c('a', 'c', 'd', 'b'))
                )
                              })
})