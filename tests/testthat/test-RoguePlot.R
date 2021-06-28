test_that('Rogue plot successful', {
   trees <- list(read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                 read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                 read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                 read.tree(text = '(a1, (b1, (c1, (d1, (rogue, e1)))));'))
   expect_equal(c(0, 0, 3, 0, 1, 0, 0, 0),
                RoguePlot(trees, 'rogue', plot = FALSE))
   
   
})