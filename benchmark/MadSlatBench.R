library("TreeSearch", quietly = TRUE)
states <- c(2, 3, 0, 4, 0, 0, 0, 2)

bench::mark(
  single = {MaddisonSlatkin_clear_cache(); MaddisonSlatkin(3L, states)},
  multi  = {MaddisonSlatkin_clear_cache(); MaddisonSlatkin_steps(0L, 4L, states)},
  max_iterations = 12,
  check = FALSE
)

bench::mark(
  single = {MaddisonSlatkin_clear_cache(); MaddisonSlatkin(3L, states)},
  multi  = {MaddisonSlatkin_clear_cache(); MaddisonSlatkin_steps(0L, 4L, states)},
  check = FALSE,
  min_time = 20
)
