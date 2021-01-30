message("\n\n\n #### ", Sys.time(), ": RUN COMPLETE ####\n")

pdf(file = paste0('results-', nTip, '-', nRep, '.pdf'),
    title = paste0("Results: ", nTip, " leaves, ", nRep, " replicates"),
    width = 9, height = 4.5)
par(mfrow = c(1, 2), mgp = c(1.9, 0.7, 0), mar = rep(3, 4))
# par(mfrow = c(2, 3), mgp = c(1.9, 0.7, 0), mar = rep(3, 4))
#boxplot(t(log(results['betterThanGen', , ])),
#        ylab = '(log) Trees better than generative', notch = TRUE, 
#        frame = FALSE, col = cols)
boxplot(log10(t(results['betterThanGen', , ] + results['equalToGen', , ])),
        notch = TRUE,
        main = paste0(nTip, ' leaves, ', nRep, ' replicates'),
        ylab = 'Log10 trees equal or better than generative', 
        frame = FALSE, col = cols)
boxplot(t(results['cidFromGen', , ]), notch = TRUE,
        ylab = 'Mean CID, best to generative',
        main = Sys.time(),
        frame = FALSE, col = cols)
# boxplot(t(results['qdFromGen', , ]),
#         ylab = 'QD from generative',
#         frame = FALSE, col = cols)
# boxplot(t(results['tbrFromGen', , ]), notch = TRUE,
#         ylab = 'TBR from generative',
#         frame = FALSE, col = cols)
# plot(type = 'n', 0, 0, main = paste0(nTip, ' tips, ', nRep, ' replicates'),
#      frame = FALSE, axes = FALSE, xlab = '', ylab = '')
dev.off()
