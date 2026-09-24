# Development library path for Shiny app testing.
#
# When shinytest2 launches this app inside the package source tree, pkgload
# intercepts library("TreeSearch") and tries to recompile from source
# (debug build). That fails when src/*.o files are stale or locked by AV.
#
# Fix: prepend .agent-shiny to .libPaths() so library("TreeSearch") finds
# the pre-built 2.0.0 installation. The version match prevents pkgload from
# deciding a recompile is needed.
#
# Only active when .agent-shiny exists (i.e. in the dev workspace).
local({
  shiny_lib <- normalizePath(
    file.path(dirname(dirname(getwd())), ".agent-shiny"),
    mustWork = FALSE
  )
  if (dir.exists(shiny_lib)) {
    .libPaths(c(shiny_lib, .libPaths()))
  }
})
