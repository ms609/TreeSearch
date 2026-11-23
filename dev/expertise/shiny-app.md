# Shiny App Expertise — TreeSearch

## Purpose

This document provides best practices and troubleshooting guidance for developing and maintaining the TreeSearch Shiny interactive application (`inst/Parsimony/app.R`). The app provides a user-friendly interface for phylogenetic tree search with real-time feedback, logging, and publication-ready visualization.

## App Architecture

### High-level Structure

```
app.R (3683 lines)
├── UI (lines 264-471)
│   ├── Left sidebar (3-column)
│   │   ├── Data loading (file, package datasets)
│   │   ├── Search controls (configure, start, save log)
│   │   ├── Tree loading and sampling
│   │   └── Display configuration (format, outgroup, etc.)
│   └── Main panel (9-column)
│       ├── Plot area with dynamic sizing
│       ├── Plot controls (size, export, concordance, clustering)
│       └── Tree/space visualization panels (conditional display)
│
├── Server (lines 506-3683)
│   ├── Logging infrastructure (Write, LogCode, LogComment, etc.)
│   ├── Data loading (UpdateData, Excel/TNT/PhyDat parsers)
│   ├── Tree management (UpdateAllTrees, UpdateActiveTrees, filtering)
│   ├── Search execution (StartSearch, MaximizeParsimony dispatch)
│   ├── Display rendering (consensus, clustering, tree space visualization)
│   ├── User interactions (observeEvent handlers, reactive computations)
│   └── Export functionality (Newick, Nexus, PDF, PNG, R script logging)
│
└── Supporting Elements
    ├── Palettes (56+ color schemes for taxa)
    ├── References (formatted bibliography)
    ├── Helper functions (Enquote, EnC, Icon, ErrorPlot)
    └── Notification system (Notification function wrapping showNotification)
```

### Key Reactive Values (lines 508-517)

- `r$dataFiles`, `r$excelFiles`, `r$treeFiles` — file counters for temp caching
- `r$dataset` — loaded phyDat object
- `r$allTrees`, `r$trees` — all vs. displayed tree subset
- `r$outgroup` — selected outgroup taxa for rooting
- `r$searchWithout` — taxa to exclude from search
- `r$sortTrees` — whether to reorder edges by clade size (for display)
- `r$plotLog`, `r$cmdLogFile` — logging outputs for export

### Data Flow

1. **Data load** → `UpdateData()` (line 797)
   - Detects file type (Excel, TNT, PhyDat)
   - Caches to temp directory
   - Logs code for reproducibility
   - Attempts to load trees from same file

2. **Search** → `StartSearch()` (line 1566)
   - Builds or uses existing starting tree
   - Dispatches to `MaximizeParsimony()` (C++ engine)
   - Logs search code with all parameters
   - Updates tree display

3. **Display** → Reactive plot rendering (lines 1731+)
   - User selects plot format (individual trees, consensus, clustering, tree space)
   - Conditional UI elements show/hide based on selection
   - Plots render via R base graphics (not ggplot2)

## Critical Functions by Purpose

### Data Loading

| Function | Lines | Role |
|----------|-------|------|
| `UpdateData()` | 797 | Main dispatcher; handles file/package sources |
| Excel parsing | 830-903 | readxl-based with skip/column controls |
| TNT/PhyDat parsing | 908-949 | Tries multiple formats; caches successfully read files |
| `CacheInput()` | 739 | Copies file to temp for reproducibility |
| Character extraction | 961 | Reads character names/notes for display |

### Tree Management

| Function | Lines | Role |
|----------|-------|------|
| `UpdateAllTrees()` | 1145 | Replace all trees; renumber tips consistently |
| `UpdateActiveTrees()` | 1086 | Thin to user-selected range and count |
| `UpdateTreeRange()` | 1067 | Sync range slider with data structures |
| `UpdateNTree()` | 1026 | Update tree count; validate against range |
| `FetchNTree()`, `FetchTreeRange()` | 1012, 1053 | Debounced reactive accessors |

### Search & Scoring

| Function | Lines | Role |
|----------|-------|------|
| `StartSearch()` | 1566 | Build starting tree, dispatch MaximizeParsimony, log code |
| `scores()` | 1344 | Cached TreeLength() call on active trees |
| `DisplayTreeScores()` | 1369 | Update results text; show score range and weighting |
| `concavity()` | 1550 | Parse IW exponent or profile mode from input |
| `weighting()` | 1332 | Map UI "on"/"off"/"prof" to concavity values |

### Rogue Taxon Detection

| Function | Lines | Role |
|----------|-------|------|
| `Rogues()` | 1775 | Cached Rogue::QuickRogue() call |
| `nNonRogues()` | 1834 | Rogue count at selected p-value |
| `KeptTips()`, `DroppedTips()` | 1949, 1973 | Filter tree tips by rogue analysis |
| `UpdateKeepNTipsRange()` | 1402 | Validate user input; sync with rogue count |

### Visualization

| Function | Lines | Role |
|----------|-------|------|
| `PlottedTree()` | 1731 | Consensus or individual tree, rooted/sorted |
| `concordance()` | 1862 | Calculate split support (multiple measures) |
| `LabelConcordance()` | 1876 | Annotate tree with support values |
| `ConsensusPlot()` | 1982 | Render consensus with rogue drop sequence |
| `TipCols()` | 1840 | Color tips by stability (Rogue::ColByStability) |

### Logging & Export

| Function | Lines | Role |
|----------|-------|------|
| `BeginLog()` | 590 | Initialize search log with system info |
| `LogCode()`, `LogComment()` | 692, 704 | Append to R script log |
| `Write()` | 524 | Append to temp log file with indentation |
| `StashTrees()` | 745 | Save trees to Nexus in temp for export |

## Best Practices

### 1. Reactive Programming Patterns

**Use `reactive()` for derived values, `bindCache()` for expensive calls:**
```r
# Simple derived value
weighting <- reactive(switch(input$implied.weights, "on" = Inf, ...))

# Cached function (re-run only if dependencies change)
scores <- bindCache(reactive({ TreeLength(r$trees, ...) }),
                     r$treeHash, r$dataHash, concavity())
```

**Avoid:**
- Direct `input$*` reads in observers (use reactive() wrapper)
- Computing the same expensive value multiple times
- Calling `reactive()` inside `observe()`/`observeEvent()`

### 2. File Handling

**Always cache input files to temp directory for reproducibility:**
```r
CacheInput("data", fileName)  # Copies to tempdir() + DataFileName(counter)
LogCode(paste0("dataFile <- \"", LastFile("data"), "\""))
```

**Supported formats (auto-detect by extension):**
- `.xlsx` / `.xls` — Excel (readxl + configurable skip/columns)
- `.nex` — Nexus (read.nexus)
- `.tre` / `.txt` — TNT or Newick (ReadTntTree or read.tree/read.nexus)
- Any phyDat-compatible text format (ReadAsPhyDat)

### 3. Logging Code Reproducibility

**Every significant user action must log equivalent R code:**
```r
LogCode(c(
  "newTrees <- MaximizeParsimony(",
  "  dataset,",
  "  concavity = 10,",
  "  maxReplicates = 100",
  ")"
))
```

**Use `EnC()` to quote parameters safely:**
```r
# EnC(c("a", "b")) → "c(\"a\", \"b\")"
# EnC("profile") → "\"profile\""
# EnC(10) → "10"
```

**Indentation via `LogIndent()` for nested scopes:**
```r
LogIndent(2)  # Indent +2 spaces
LogCode("for (tree in trees) {")
LogIndent(2)
LogCode("  tree <- Consensus(tree, p = 0.5)")
LogIndent(-2)
LogCode("}")
LogIndent(-2)
```

### 4. Observing User Input

**Use debounce for high-frequency inputs (sliders, text boxes):**
```r
PlottedChar <- debounce(reactive({ as.integer(input$plottedChar) }), aJiffy)
```

**Use `ignoreInit = TRUE` to skip initialization:**
```r
observeEvent(input$searchConfig, { ... }, ignoreInit = TRUE)
```

**Cache tree hashes to detect changes (avoid spurious recalculations):**
```r
observeEvent(r$dataset, {
  r$dataHash <- rlang::hash(r$dataset)
})
r$trees <- thinnedTrees
r$treeHash <- rlang::hash(r$trees)
```

### 5. Conditional UI & Show/Hide Elements

**Use bslib-style id-based show/hide (not class-based):**
```r
# Define in UI with hidden(...) wrapper
hidden(tags$div(id = "displayConfig", ...))

# Toggle in server
show("displayConfig", anim = TRUE)    # With fade-in animation
hide("displayConfig")                  # Fade-out
showElement("displayConfig")           # JavaScript show() without animation
hideElement("displayConfig")
```

**Manage multiple related configs via `ShowConfigs()`:**
```r
observeEvent(input$plotFormat, {
  ShowConfigs(switch(input$plotFormat,
    "ind"   = c("whichTree", "charChooser", "treePlotConfig"),
    "cons"  = c("consConfig", "branchLegend", "savePlottedTrees"),
    "clus"  = c("clusConfig", "clusLegend", "savePlottedTrees"),
    ""      # Default: hide all
  ))
})
```

### 6. Modal Dialogs for Configuration

**Example: Search configuration modal (line 1220):**
```r
observeEvent(input$searchConfig, {
  # Pre-populate with current values
  updateSelectInput(session, "concavity", selected = input$concavity)
  
  showModal(modalDialog(
    fluidPage(column(6, ...), column(6, ...)),
    title = "Tree search settings",
    footer = tagList(
      modalButton("Close", icon = Icon("rectangle-xmark")),
      actionButton("modalGo", "Start search", icon = Icon("magnifying-glass"))
    ),
    easyClose = TRUE
  ))
})

observeEvent(input$modalGo, {
  removeModal()
  StartSearch()
})
```

## Common Issues & Troubleshooting

### Issue 1: File Upload Not Working

**Symptom:** User selects file, nothing happens.

**Checks:**
- File size < `shiny.maxRequestSize` (default 5MB; app sets 1GB at line 4)
- File extension recognized (Excel, TNT, Nexus, text)
- `readxl` installed for Excel files (auto-install at line 831)
- Check browser console for error messages
- If TNT format: tip labels must be inferrable (will try 4 caterpillar orderings)

### Issue 2: Search Hangs or No Results

**Symptom:** Click "Search", progress bar shows, but never completes.

**Checks:**
- Dataset is valid phyDat (not NULL, has tips)
- Tree space not empty or trivial (≥4 tips recommended)
- Replicates/timeout reasonable (maxReplicates ≥ 1, timeout > search time)
- Check `maxSeconds` timeout — if 0, no timeout; if very small, search aborts early
- Parallel mode (nThreads > 1) is non-deterministic; may find different trees

**Debugging:**
```r
# In console:
ds <- ReadAsPhyDat("data.nex")
attr(ds, "nr")  # Check character count
length(ds)      # Check taxon count
tree <- AdditionTree(ds)  # Should complete quickly
```

### Issue 3: Trees Don't Display / Blank Plot

**Symptom:** Plot area is empty; no error message.

**Checks:**
- Trees loaded? (r$trees length > 0)
- Dataset loaded? (needed for consensus/character display)
- Display format selected? (default "cons" should show something)
- Outgroup valid? (must be in tree tips)
- Rogue-dropping valid? (can't drop all tips)

**Debugging:**
```r
# In console:
length(app_env$r$trees)                    # Should be > 0
app_env$AnyTrees()                         # Should be TRUE
app_env$Consensus(app_env$r$trees, p=1)   # Should render
```

### Issue 4: Logging Code Mismatch

**Symptom:** Exported R script doesn't reproduce results.

**Checks:**
- File paths in log correct? (should use temp files like "dataFile-00.txt")
- Parameters logged correctly? (check `Enquote()` results)
- Library calls present? (BeginLog should include all imports)
- Character encoding OK? (use system-appropriate paths)

**Prevention:**
- Always use `LogCode()` immediately after performing an action
- Test exported script manually in a fresh R session
- Check `tempdir()` for actual cached files

### Issue 5: Rogue Analysis Crashes or Misses Taxa

**Symptom:** `Rogues()` returns NULL, or taxa don't appear in drop sequence.

**Checks:**
- Dataset properly loaded (not NULL)
- Trees properly loaded (at least 1 tree, tip labels match)
- `p` parameter reasonable (0.5 to 1.0; default 1.0 = strict majority rule)
- Run `Rogue::QuickRogue()` manually to test:
  ```r
  rogues <- Rogue::QuickRogue(r$trees, neverDrop = input$neverDrop, 
                              fullSeq = TRUE, p = consP())
  ```

### Issue 6: Memory Leak or Slowdown Over Time

**Symptom:** App slows down after many searches; process memory grows.

**Checks:**
- File caching in `tempdir()` consuming space? (e.g., 1000 searches → 1000s of cached files)
- Large tree objects retained? (clear old results before new search)
- Image caches building up? (plots rendered reactively, may leak if observer not cleaned up)

**Prevention:**
- Periodically clear `tempdir()` (not auto-cleared by default)
- Use `on.exit()` to clean up temporary objects:
  ```r
  observeEvent(input$clearCache, {
    do.call(file.remove, list(dir(tempdir(), full.names=TRUE)))
    Notification("Cache cleared", type="message")
  })
  ```

## Integration with C++ Engine

### Strategy Presets (line 1231)

- **"auto"** — Auto-selects based on dataset size (sprint ≤30, default 31-60, thorough 61+)
- **"sprint"** — 3 ratchet cycles, no drift; minimal sectorial
- **"default"** — 5 ratchet, 2 drift; XSS+RSS+CSS
- **"thorough"** — 20 ratchet, 12 drift; intensive sectorial; adaptive ratchet

### Weighting Mode (line 1224)

- **"on"** (Implied) — IW with concavity exponent (k = 10^exponent)
- **"off"** (Equal) — EW (all characters weight 1)
- **"prof"** (Profile) — Profile parsimony (info-theoretic weighting)

## Testing Checklist

Before deploying app updates:

- [ ] Data loads: Excel (with skip/columns), TNT, Nexus, generic text
- [ ] Search runs: EW, IW, profile; small (4 tips), medium (25), large (75+)
- [ ] Logging: exported R script runs in fresh session, reproduces trees
- [ ] Display: individual, consensus, clustering, tree space all render
- [ ] Rogue analysis: correctly identifies and drops unstable taxa
- [ ] Outgroup: rooting works; must be in tree and dataset
- [ ] Export: PDF, PNG, Newick, Nexus files valid
- [ ] Performance: 50+ searches don't slow app significantly
- [ ] Parallel: nThreads=2 works; results reasonable (non-deterministic)
- [ ] Edge cases: 3-tip tree, single-character dataset, all inapplicable, empty pool

## Performance Tips

1. **Limit active tree display** — reduce `whichTree` max range if >100 trees
2. **Cache tree hashes** — avoid re-scoring unchanged trees
3. **Use bounded indirect** — ensure TBR/drift/SPR use `_bounded` variants
4. **Debounce slider inputs** — high-frequency slider updates (default aJiffy ≈ 42ms)
5. **Profile big plots** — use `system.time({ ... })` for consensus/space rendering

## References

- **app.R**: Main application file (3683 lines)
- **Related packages**: shiny, shinyjs, bslib, TreeTools, TreeSearch, Rogue, TreeDist
- **C++ search**: MaximizeParsimony() documented in `R/MaximizeParsimony.R`
- **Logging infrastructure**: BeginLog, LogCode, Write functions (lines 590-715)
