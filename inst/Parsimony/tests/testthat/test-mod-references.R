library(shiny)

# Source the module under test (relative to tests/testthat/)
source("../../server/mod_references.R")

# Minimal citation stubs — actual HTML is defined in global.R; stubs are
# sufficient here since we only need to verify rendering logic.
stub_cites <- list(
  Brazeau2019    = "<p>Brazeau2019</p>",
  Goloboff1993   = "<p>Goloboff1993</p>",
  Goloboff1999   = "<p>Goloboff1999</p>",
  Goloboff2014   = "<p>Goloboff2014</p>",
  Nixon1999      = "<p>Nixon1999</p>",
  SmithSearch    = "<p>SmithSearch</p>",
  Gower1966      = "<p>Gower1966</p>",
  Gower1969      = "<p>Gower1969</p>",
  Kaski2003      = "<p>Kaski2003</p>",
  RCoreTeam      = "<p>RCoreTeam</p>",
  SmithDist      = "<p>SmithDist</p>",
  Smith2020      = "<p>Smith2020</p>",
  SmithSpace     = "<p>SmithSpace</p>",
  Venna2001      = "<p>Venna2001</p>",
  Stockham2002   = "<p>Stockham2002</p>",
  Arthur2007     = "<p>Arthur2007</p>",
  Hartigan1979   = "<p>Hartigan1979</p>",
  Maechler2019   = "<p>Maechler2019</p>",
  Bien2011       = "<p>Bien2011</p>",
  Murtagh1983    = "<p>Murtagh1983</p>",
  Rousseeuw1987  = "<p>Rousseeuw1987</p>",
  SmithRogue     = "<p>SmithRogue</p>",
  Klopfstein2019 = "<p>Klopfstein2019</p>",
  Pol2009        = "<p>Pol2009</p>"
)

test_that("references_server renders section headings", {
  shiny::testServer(references_server, args = list(cites = stub_cites), {
    result <- output$references
    expect_false(is.null(result))
    rendered <- paste(as.character(result), collapse = "")
    expect_true(grepl("Tree search", rendered, fixed = TRUE))
    expect_true(grepl("Clustering",  rendered, fixed = TRUE))
    expect_true(grepl("Rogue taxa",  rendered, fixed = TRUE))
  })
})

test_that("EW mode shows standing refs but not IW/XPIWE refs", {
  wt <- reactiveVal("off")
  shiny::testServer(references_server,
    args = list(weighting = wt, cites = stub_cites), {
    rendered <- paste(as.character(output$references), collapse = "")
    expect_true(grepl("SmithSearch",   rendered, fixed = TRUE))
    expect_true(grepl("Goloboff1999",  rendered, fixed = TRUE))
    expect_true(grepl("Nixon1999",     rendered, fixed = TRUE))
    expect_true(grepl("Brazeau2019",   rendered, fixed = TRUE))
    expect_false(grepl("Goloboff1993", rendered, fixed = TRUE))
    expect_false(grepl("Goloboff2014", rendered, fixed = TRUE))
  })
})

test_that("IW mode adds Goloboff 1993 but not Goloboff 2014", {
  wt <- reactiveVal("on")
  shiny::testServer(references_server,
    args = list(weighting = wt, cites = stub_cites), {
    rendered <- paste(as.character(output$references), collapse = "")
    expect_true(grepl("Goloboff1993",  rendered, fixed = TRUE))
    expect_false(grepl("Goloboff2014", rendered, fixed = TRUE))
  })
})

test_that("XPIWE mode adds both Goloboff 1993 and 2014", {
  wt <- reactiveVal("xpiwe")
  shiny::testServer(references_server,
    args = list(weighting = wt, cites = stub_cites), {
    rendered <- paste(as.character(output$references), collapse = "")
    expect_true(grepl("Goloboff1993", rendered, fixed = TRUE))
    expect_true(grepl("Goloboff2014", rendered, fixed = TRUE))
  })
})
