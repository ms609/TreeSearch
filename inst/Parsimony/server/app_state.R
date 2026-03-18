# AppState: Centralized reactive state for the TreeSearch Shiny app.
#
# All reactive values used by server modules are defined here with explicit
# initial values and domain grouping. This replaces the ad-hoc
# reactiveValues() call in server.R.
#
# Usage in server.R:
#   r <- AppState()
#
# Modules access fields via r$fieldName. See field documentation below.

AppState <- function() {
  reactiveValues(

    # -- Data domain --
    # Primary dataset and metadata loaded from file
    dataset       = NULL,          # phyDat object (or NULL before load)
    chars         = NULL,          # character matrix from ReadCharacters()
    charNotes     = NULL,          # character notes from ReadNotes()
    dataHash      = NULL,          # rlang::hash() of dataset (change trigger)
    dataFileVisible = TRUE,        # whether file-upload UI is shown
    readDataFile  = NULL,          # string: R expression used to read data file

    # -- File tracking (logging) --
    # Counters for unique file uploads per session (used by logging)
    dataFiles     = 0,             # count of data file uploads
    excelFiles    = 0,             # count of Excel file uploads
    treeFiles     = 0,             # count of tree file uploads

    # -- Tree domain --
    # Trees loaded from files or produced by search
    allTrees      = NULL,          # multiPhylo: full tree set (unsorted/unfiltered)
    trees         = NULL,          # multiPhylo: active subset (after range/thin)
    treeHash      = NULL,          # rlang::hash() of trees (change trigger)
    newTrees      = NULL,          # multiPhylo: trees from most recent search
    sortTrees     = FALSE,         # logical: sort trees by score before display
    readTreeFile  = NULL,          # string: R expression used to read tree file

    # -- Tree subsetting state --
    nTree         = 0L,            # integer: current max trees to display
    treeRange     = c(1L, 1L),     # integer[2]: active range of tree indices
    updatingTrees = FALSE,         # reentrancy guard for UpdateActiveTrees()

    # -- "Old" values for change detection --
    # These track previous input values so observers can detect real changes
    # vs reactive re-fires. Will be removed when modules handle own state.
    oldNTree      = NULL,          # previous value of input$nTree
    oldTreeRange  = NULL,          # previous value of input$treeRange
    oldOutgroup   = NO_OUTGROUP,   # previous value of input$outgroup
    oldkeepNTips  = NULL,          # previous value of input$keepNTips

    # -- Search domain --
    searchCount        = 0L,       # integer: how many searches have been run
    searchDataHash     = NULL,     # hash of dataset at search time
    searchNotification = NULL,     # Shiny notification ID (for dismissal)
    searchWithout      = NULL,     # character: taxa excluded from search
    bestSearchScore    = NULL,     # numeric: best score seen across all searches (for accumulation)
    searchTotalHits    = 0L,       # integer: cumulative hits_to_best across runs at current best score
    searchTotalReps    = 0L,       # integer: cumulative replicates run across runs at current best score

    # -- Consensus / plotting domain --
    outgroup      = NULL,          # character: selected outgroup taxa
    keepNTips     = NULL,          # integer: tips retained in consensus
    plottedTree   = NULL,          # phylo or list: tree(s) currently plotted
    concordance   = list(),        # list: cached concordance results
    plotLog       = NULL,          # character vector: R code log for plot

    # -- Cluster domain --
    # (r$cluster is not a state field; clustering.R uses local variables)

    # -- UI state --
    visibleConfigs = NULL          # character: which config panels are visible
  )
}
