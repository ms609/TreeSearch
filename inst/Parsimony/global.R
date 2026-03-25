# options("TreeSearch.logging" = TRUE) # Log function entry and exit
# options("TreeSearch.write.code" = TRUE) # Show code as it is written to log
logging <- isTRUE(getOption("TreeSearch.logging"))
options(shiny.maxRequestSize = 1024 ^ 3) # Allow max 1 GB files

# Development: prepend .agent-shiny library so library("TreeSearch") finds
# the pre-built v2.0.0 install, preventing pkgload from intercepting and
# attempting a debug recompile (which fails when src/*.o files are stale).
local({
  shiny_lib <- normalizePath(
    file.path(dirname(dirname(getwd())), ".agent-shiny"),
    mustWork = FALSE
  )
  if (dir.exists(shiny_lib)) {
    .libPaths(c(shiny_lib, .libPaths()))
  }
})

library("methods", exclude = c("show", "removeClass"))
library("cli")
library("TreeSearch") # load now: inapplicable.datasets required within ui
.DateTime <- function() { # Copy, because not exported
  format(Sys.time(), "%Y-%m-%d %T")
}

local({
  needed <- c("cluster", "future", "PlotTools", "promises",
              "protoclust", "Rogue", "shinyjs")
  miss <- needed[!vapply(needed, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(miss)) {
    message("Installing packages required by EasyTrees(): ",
            paste(miss, collapse = ", "))
    utils::install.packages(miss)
  }
})

suppressPackageStartupMessages({
  library("shiny", exclude = c("runExample"))
  library("shinyjs", exclude = c("runExample"))
})
library("TreeTools", quietly = TRUE)
library("TreeDist", quietly = TRUE)
library("future")
library("promises")


if (logging) {
  logMsgFile <- file("log.lg", open = "w+")
  LogMsg <- function (...) {
    message(.DateTime(), ": ", ...)
    writeLines(.DateTime(), con = logMsgFile)
    writeLines(paste0("  ", ...), con = logMsgFile)
  }
  Put <- function (..., file) {
    dput(..., file = file)
    writeLines(gsub("<pointer: [^.]+>", "NULL", readLines(file)),
               file)
  }
  PutTree <- function (...) {
    Put(..., file = "tree.lg")
  }
  PutData <- function (...) {
    Put(..., file = "dataset.lg")
  }
} else {
  PutData <- PutTree <- LogMsg <- function (...) {}
}

WriteLoggedCode <- if (isTRUE(getOption("TreeSearch.write.code"))) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    function(txt) {
      for (line in txt) cat(if (substr(trimws(line), 0, 1) == "#") {
        crayon::green("  ", line, "\n")
      } else {
        crayon::yellow("  ", line, "\n")
      })
    }
  } else {
    function(txt) message("       ", txt)
  }
} else {
  function(txt) {}
}

Notification <- function (...) {
  if (!isTRUE(getOption("shiny.testmode"))) {
    showNotification(...)
  }
}

Icon <- function(...) icon(..., class = "fas")

aJiffy <- 42 # ms, default debounce period for input sliders etc
typingJiffy <- 2.5 * aJiffy # slightly slower if might be typing
aFewTrees <- 48L # Too many and rogues / tree space are slowed
NO_OUTGROUP <- "! TREESEARCH_no outgroup specified ."

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
         if (journal != "") paste0("<i>", journal, "</i> ") else "",
         if (is.null(volume)) "" else paste0("<b>", volume, "</b>:"),
         if (is.null(publisher)) "" else paste0(publisher, ". "),
         if (is.null(pages)) "" else paste0(paste0(pages, collapse = "&ndash;"), ". "),
         if (is.null(doi)) "" else paste0(
           "doi:<a href=\"https://doi.org/", doi, "\" title=\"CrossRef\">",
           doi, "</a>. "), 
         "</p>")
}


Arthur2007 <- Reference(
  c("Arthur, D.", "Vassilvitskii, S"),
  title = "k-means++: the advantages of careful seeding",
  year = 2007,
  journal = "Proceedings of the Eighteenth Annual ACM-SIAM Symposium on Discrete Algorithms",
  pages = c(1027, 1035)
)
Brazeau2019 <- Reference(c("Brazeau, M.D.", "Guillerme, T.", "Smith, M.R."), 2019,
                           title = "An algorithm for morphological phylogenetic analysis with inapplicable data",
                           journal = "Systematic Biology",
                           volume = 64,
                           pages = c(619, 631),
                         doi = "10.1093/sysbio/syy083")
Goloboff2014 <- Reference("Goloboff, P.A.", 2014,
                           "Extended implied weighting",
                           "Cladistics", volume = 30,
                           pages = c(260, 272),
                           doi = "10.1111/cla.12047")
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
  journal = "Journal of the Royal Statistical Society Series C (Applied Statistics)")
Hartigan1979 <- Reference(
  title = "Algorithm AS 136: a <i>K</i>-means clustering algorithm",
  authors = c("Hartigan, J.A.", "Wong, M.A."),
  journal = "Journal of the Royal Statistical Society Series C (Applied Statistics)",
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
  title = "cluster: cluster analysis basics and extensions", year = 2022,
  authors = c("Maechler, M.", "Rousseeuw, P.", "Struyf, A.", "Hubert, M.", "Hornik, K."),
  journal = "Comprehensive R Archive Network")
Morphy <- Reference(
  c("Brazeau, M.D.", "Smith, M.R.", "Guillerme, T."), 2017,
  "MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability",
  doi = "10.5281/zenodo.815371")
Murtagh1983 <- Reference(
  title = "A survey of recent advances in hierarchical clustering algorithms",
  authors = "Murtagh, F.", year = 1983, volume = 26, pages = c(354, 359),
  doi = "10.1093/comjnl/26.4.354", journal = "The Computer Journal")
Nixon1999 <- Reference(
  "Nixon, K.C.", 1999,
  journal = "Cladistics", volume = 15, pages = c(407, 414),
  title = "The Parsimony Ratchet, a new method for rapid parsimony analysis",
  doi = "10.1111/j.1096-0031.1999.tb00277.x")
Pol2009 <- Reference(
  title = "Unstable taxa in cladistic analysis: identification and the assessment of relevant characters",
  authors = c("Pol, D.", "Escapa, I.H."),
  journal = "Cladistics", 2009, 25, pages = c(515, 527), 
  doi = "10.1111/j.1096-0031.2009.00258.x")
RCoreTeam <- Reference(
  authors = "R Core Team", year = 2020,
  title = "R: A language and environment for statistical computing",
  publisher = "R Foundation for Statistical Computing, Vienna, Austria")
Rousseeuw1987 <- Reference(
  title = "Silhouettes: a graphical aid to the interpretation and validation of cluster analysis",
  author = "Rousseeuw, P.J.", year = 1987,
  journal = "Journal of Computational and Applied Mathematics",
  volume = 20, pages = c(53, 65), doi = "10.1016/0377-0427(87)90125-7"
)
SmithDist <- Reference(
  "Smith, M.R.", "2020a", "TreeDist: distances between phylogenetic trees",
  doi = "10.5281/zenodo.3528123", "Comprehensive R Archive Network")
SmithQuartet <- Reference(
  "Smith, M.R.", 2019,
  "Quartet: comparison of phylogenetic trees using quartet and split measures",
  "Comprehensive R Archive Network", doi = "10.5281/zenodo.2536318")
SmithSearch <- Reference(
  "Smith, M.R.", 2023, "TreeSearch: morphological phylogenetic analysis in R",
  "R Journal", volume = 14, pages = c(305, 315),
  doi = "10.32614/RJ-2023-019")
Smith2020 <- Reference(
  "Smith, M.R.", "2020b",
  "Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees",
  "Bioinformatics", volume = 36, pages = c("5007", "5013"),
  doi = "10.1093/bioinformatics/btaa614")
SmithSpace <- Reference(
  "Smith, M.R.", "2022a", "Robust analysis of phylogenetic tree space",
  "Systematic Biology", 71, pages = c("1255", "1270"),
  doi = "10.1093/sysbio/syab100")
SmithRogue <- Reference(
  "Smith, M.R.", "2022b",
  "Using information theory to detect rogue taxa and improve consensus trees",
  "Systematic Biology", 71, pages = c("1088", "1094"),
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















Enquote <- function(x, ...) {
  if (mode(x) == "character") {
    paste0("\"", x, "\"")
  } else {
    signif(x, ...)
  }
}

#' Confidence text for post-search results display.
#'
#' Given K hits to best score in R total runs, returns a plain-text
#' summary: "K of R runs hit best score. Probability that a better tree
#' exists: ~X%".
#'
#' @param K integer. Cumulative hits to best score.
#' @param R integer. Cumulative runs completed.
#' @return character(1) or NULL if no search data.
FormatMissProb <- function(prob) {
  pct <- prob * 100
  if (pct >= 1) paste0("~", round(pct), "%")
  else if (pct >= 0.1) "<1%"
  else if (pct >= 0.01) "<0.1%"
  else "<0.01%"
}

SearchConfidenceText <- function(K, R, nSearches = 1L,
                                 nTopologies = NULL,
                                 lastImprovedRep = NULL,
                                 stopReason = NULL) {
  if (is.null(K) || is.null(R) || R <= 0L || K <= 0L) return(NULL)
  K <- min(K, R)

  # Tightened binomial bound: (1 - K/R)^R is tighter than exp(-K) when K < R.
  # Falls back to exp(-K) when K == R, since (1 - 1)^R = 0 is overconfident.
  prob_miss <- if (K < R) (1 - K / R) ^ R else exp(-K)

  runs_label <- if (!is.null(nSearches) && nSearches > 1L) {
    paste0("total runs across ", nSearches, " searches")
  } else {
    "runs"
  }

  # Only warn when a single topology limits the independence assumption
  topo_note <- if (!is.null(nTopologies) && nTopologies == 1L) {
    " [single topology \u2014 limited independence]"
  } else {
    ""
  }

  # Trajectory info
  trajectory_note <- if (!is.null(lastImprovedRep) && R > 1L) {
    paste0(" Last improvement: replicate ", lastImprovedRep, ".")
  } else {
    ""
  }

  # Landscape ruggedness flag
  rugged_note <- if (K / R < 0.3 && R >= 5L) {
    paste0(" Hit rate low (", round(100 * K / R),
           "%) \u2014 more replicates may help.")
  } else {
    ""
  }

  # Nudge for small K == R
  small_sample_note <- if (K == R && R <= 5L) {
    paste0(" \u2014 increase \u2018Stop when N runs hit best\u2019 for a ",
           "tighter estimate")
  } else {
    ""
  }

  stop_note <- if (identical(stopReason, "consensus")) {
    " Search stopped: consensus tree unchanged across recent replicates."
  } else if (identical(stopReason, "timeout")) {
    " Search stopped: time limit reached."
  } else {
    ""
  }

  paste0(K, " of ", R, " ", runs_label, " hit best score. ",
         "Probability that a better score exists: ",
         FormatMissProb(prob_miss),
         topo_note, trajectory_note, rugged_note, small_sample_note,
         stop_note)
}

EnC <- function(...) {
  if (length(...) == 1) {
    Enquote(...)
  } else {
    paste0("c(", paste(sapply(..., Enquote), collapse = ", "), ")")
  }
}

# Shiny modules — sourced here so ui.R can call xxx_ui() at definition time
source("server/mod_references.R")
source("server/mod_downloads.R")
dl_ui <- downloads_ui("dl")
source("server/mod_search.R")
se_ui <- search_ui("search")
source("server/mod_data.R")
source("server/mod_clustering.R")
source("server/mod_treespace.R")
source("server/mod_consensus.R")
data_ui_elems <- data_ui("data")
co_ui <- consensus_ui("consensus")
