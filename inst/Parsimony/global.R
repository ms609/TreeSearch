
aJiffy <- 42 # ms, default debounce period for input sliders etc
typingJiffy <- 2.5 * aJiffy # slightly slower if might be typing
aFewTrees <- 48L # Too many and rogues / tree space are slowed

palettes <- list("#7a6c36",
                 c("#7a6c36", "#864885"),
                 c("#7a6c36", "#864885", "#427743"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020", "#c241a7"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#824eca", "#b3622a", "#452580", "#417f81", "#ca4172", "#6171ca", "#364020", "#c241a7", "#391d42"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#85c6f9", "#fbd1a0", "#7696be", "#89996c", "#ddcdff", "#719d89", "#f5cde6", "#b6e0da", "#e8d4cd", "#b5ddfa"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#bbcb8f", "#bf82ab", "#85ddc4", "#eea0ba", "#c1d8ff", "#c3818b", "#c5c6ff", "#999388", "#e8cbff", "#ffb5b6", "#d2dad7"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#cc8f6f", "#499fae", "#d9dca6", "#7796b8", "#bee1ba", "#b4daff", "#919583", "#e2d3e9", "#47a19b", "#ebd4bc", "#7c9993", "#a9e3e0"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#a8e0fe", "#fad0a8", "#679e8d", "#ffc7b1", "#abe5c0", "#ac8d78", "#c5dddc", "#a48f84", "#cadfb0", "#899694", "#fdcdc1", "#d1dad5", "#dfd8c4"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#dcb983", "#77bff0", "#f0ab92", "#90ddff", "#f1d3a9", "#b5c2fe", "#c1e1b7", "#7596ba", "#bce1c4", "#a88c96", "#5a9daf", "#b18b80", "#d4d6f3", "#949577"),
                 c("#7a6c36", "#864885", "#427743", "#4c5c86", "#cb4745", "#73383b", "#e03795", "#438f2e", "#5e2195", "#758029", "#4042b9", "#a37926", "#8364df", "#c3671f", "#444491", "#dc4c1f", "#367076", "#e2383c", "#4786b4", "#e13964", "#4c8c73", "#a53396", "#2c4422", "#b553cb", "#50381b", "#4f75d8", "#a12c1b", "#8576b8", "#bd6541", "#3a1959", "#83491f", "#2d2644", "#c45b94", "#451523", "#966883", "#782224", "#b96563", "#762254", "#95765c", "#ad355a")
)

ErrorPlot <- function (...) {
  plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
  text(0, 0, paste0(..., collapse = "\n"),
       col = "#dd6611", font = 2)
}

badToGood <- rev(c("#1AB958", "#23B956", "#2BB954", "#31B952", "#37B850", "#3CB84E", "#41B84C", "#45B74A", "#49B749", "#4DB747", "#51B645", "#54B643", "#58B641", "#5BB53F", "#5FB53D", "#62B53C", "#65B43A", "#68B438", "#6BB336", "#6DB335", "#70B333", "#73B231", "#76B230", "#78B12E", "#7BB12C", "#7DB02B", "#80B029", "#82AF28", "#85AF26", "#87AE25", "#8AAE23", "#8CAD22", "#8EAD21", "#91AC1F", "#93AC1E", "#95AB1D", "#97AB1C", "#9AAA1B", "#9CAA1A", "#9EA919", "#A0A918", "#A2A818", "#A4A717", "#A6A716", "#A8A616", "#AAA616", "#ACA515", "#AEA415", "#B0A415", "#B2A315", "#B4A315", "#B6A216", "#B8A116", "#B9A117", "#BBA017", "#BD9F18", "#BF9F18", "#C19E19", "#C29D1A", "#C49D1B", "#C69C1C", "#C79B1D", "#C99A1E", "#CB9A1F", "#CC9920", "#CE9822", "#CF9823", "#D19724", "#D29625", "#D49626", "#D59528", "#D79429", "#D8932A", "#D9932C", "#DB922D", "#DC912E", "#DD9130", "#DF9031", "#E08F33", "#E18F34", "#E28E35", "#E38D37", "#E58C38", "#E68C3A", "#E78B3B", "#E88A3D", "#E98A3E", "#EA8940", "#EB8841", "#EC8843", "#ED8744", "#EE8746", "#EE8647", "#EF8549", "#F0854A", "#F1844C", "#F2844D", "#F2834F", "#F38350", "#F48252", "#F48253", "#F58155", "#F58157", "#F68058", "#F6805A", "#F77F5B", "#F77F5D", "#F87E5E"))

Reference <- function (authors, year, title, journal = "",
                       volume = NULL, pages = NULL, doi = NULL,
                       publisher = NULL, editors = NULL) {
  nAuth <- length(authors)
  if (nAuth > 1L) {
    authors <- paste(paste0(authors[-nAuth], collapse = ", "), "&amp;", authors[nAuth])
  }
  nEd <- length(editors)
  if (nEd > 1L) {
    editors <- paste(paste0(editors[-nEd], collapse = ", "), "&amp;", editors[nEd])
  } else if (nEd < 1) {
    editors <- ""
  }
  paste0("<p class=\"reference\">", authors, " (", year, "). &ldquo;", title,
         "&rdquo;. ",
         if (editors != "") paste0("In: ", editors, " (eds). ") else "",
         if (journal != "") paste0("<i>", journal, "</i>. ") else "",
         if (is.null(volume)) "" else paste0("<b>", volume, "</b>:"),
         if (is.null(publisher)) "" else paste0(publisher, ". "),
         if (is.null(pages)) "" else paste0(paste0(pages, collapse = "&ndash;"), ". "),
         if (is.null(doi)) "" else paste0(
           "doi:<a href=\"https://doi.org/", doi, "\" title=\"CrossRef\">",
           doi, "</a>. "), 
         "</p>")
}


Brazeau2019 <- Reference(c("Brazeau, M.D.", "Guillerme, T.", "Smith, M.R."), 2019,
                         title = "An algorithm for morphological phylogenetic analysis with inapplicable data",
                         journal = "Systematic Biology",
                         volume = 64,
                         pages = "619-631",
                         doi = "10.1093/sysbio/syy083")
Bien2011 <- Reference(
  c("Bien, J.", "Tibshirani, R."),
  title = "Hierarchical clustering with prototypes via minimax linkage",
  year = 2011,
  volume = 106,
  doi = "10.1198/jasa.2011.tm10183",
  pages = c(1075, 1084),
  journal = "Journal of the American Statistical Association")
Gower1966 <- Reference(title = "Some distance properties of latent root and vector methods used in multivariate analysis",
                       authors = "Gower, J.C.",
                       year = 1966,
                       volume = 53,
                       pages = c(325, 338),
                       doi = "10.2307/2333639",
                       journal = "Biometrika")
Gower1969 <- Reference(
  title = "Minimum spanning trees and single linkage cluster analysis",
  authors = c("Gower, J.C.", "Ross, G.J.S."),
  year = 1969, volume = 18, pages = c(54, 64), doi = "10.2307/2346439",
  journal = "Journal of the Royal Statistical Society. Series C (Applied Statistics)")
Hartigan1979 <- Reference(
  title = "Algorithm AS 136: a <i>K</i>-means clustering algorithm",
  authors = c("Hartigan, J.A.", "Wong, M.A."),
  journal = "Journal of the Royal Statistical Society. Series C (Applied Statistics)",
  year = 1979, volume = 28, pages = c(100, 108),
  doi = "10.2307/2346830")
Kaski2003 <- Reference(
  title = "Trustworthiness and metrics in visualizing similarity of gene expression",
  authors = c("Kaski, S.", "Nikkil&auml;, J.", "Oja, M.", "Venna, J.",
              "T&ouml;r&ouml;nen, P.", "Castr&eacute;n, E."),
  year = 2003, volume = 4, pages = 48, doi = "10.1186/1471-2105-4-48",
  journal = "BMC Bioinformatics")
Klopfstein2019 <- Reference(
  title = "Illustrating phylogenetic placement of fossils using RoguePlots: An example from ichneumonid parasitoid wasps (Hymenoptera, Ichneumonidae) and an extensive morphological matrix.",
  authors = c("Klopfstein, S.", "Spasojevic, T."), year = 2019, 
  journal = "PLoS ONE", volume = 14, pages = "e0212942",
  doi = "10.1371/journal.pone.0212942"
)
Maechler2019 <- Reference(
  title = "cluster: cluster analysis basics and extensions", year = 2019,
  authors = c("Maechler, M.", "Rousseeuw, P.", "Struyf, A.", "Hubert, M.", "Hornik, K."),
  journal = "Comprehensive R Archive Network")
Morphy <- Reference(c("Brazeau, M.D.", "Smith, M.R.", "Guillerme, T."), 2017,
                    "MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability.",
                    doi = "10.5281/zenodo.815371")
Murtagh1983 <- Reference(
  title = "A survey of recent advances in hierarchical clustering algorithms",
  authors = "Murtagh, F.", year = 1983, volume = 26, pages = c(354, 359),
  doi = "10.1093/comjnl/26.4.354", journal = "The Computer Journal")
Nixon1999 <- Reference(
  "Nixon, K.C.", 1999,
  journal = "Cladistics", volume = 15, pages = "407-414",
  title = "The Parsimony Ratchet, a new method for rapid parsimony analysis",
  doi = "10.1111/j.1096-0031.1999.tb00277.x")
RCoreTeam <- Reference(
  authors = "R Core Team", year = 2020,
  title = "R: A language and environment for statistical computing",
  publisher = "R Foundation for Statistical Computing, Vienna, Austria")
SmithDist <- Reference("Smith, M.R.", 2020,
                       "TreeDist: distances between phylogenetic trees",
                       doi = "10.5281/zenodo.3528123", "Comprehensive R Archive Network")
SmithQuartet <- Reference("Smith, M.R.", 2019,
                          "Quartet: comparison of phylogenetic trees using quartet and split measures",
                          "Comprehensive R Archive Network", doi = "10.5281/zenodo.2536318")
SmithSearch <- Reference("Smith, M.R.", 2018,
                         "TreeSearch: phylogenetic tree search using custom optimality criteria",
                         "Comprehensive R Archive Network", doi = "10.5281/zenodo.1042590")
Smith2020 <- Reference("Smith, M.R.", 2020,
                       "Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees",
                       "Bioinformatics", volume = 36, pages = "5007--5013",
                       doi = "10.1093/bioinformatics/btaa614")
SmithSpace <- Reference("Smith, M.R.", "2022b",
                        "Robust analysis of phylogenetic tree space",
                        "Systematic Biology", pages = "syab100",
                        doi = "10.1093/sysbio/syab100")
SmithRogue <- Reference("Smith, M.R.", "2022a", 
                        "Using information theory to detect rogue taxa and improve consensus trees",
                        "Systematic Biology", pages = "syab099",
                        doi = "10.1093/sysbio/syab099")
Stockham2002 <- Reference(
  authors = c("Stockham, C.", "Wang, L.-S.", "Warnow, T."), 2002,
  "Statistically based postprocessing of phylogenetic analysis by clustering",
  "Bioinformatics", 18, c("S285", "S293"),
  doi = "10.1093/bioinformatics/18.suppl_1.S285")

Venna2001 <- Reference(
  title = "Neighborhood preservation in nonlinear projection methods: an experimental study",
  authors = c("Venna, J.", "Kaski, S."), year = 2001, pages = c(485, 491),
  journal = "Lecture Notes in Computer Science: Artificial Neural Networks&mdash;ICANN 2001",
  editors = c("Dorffner, G.", "Bischof, H.", "Hornik, K."),
  publisher = "Springer, Berlin",
  doi = "10.1007/3-540-44668-0_68")
