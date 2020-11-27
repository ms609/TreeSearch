message("\n\n\n #### ", Sys.time(), ": RUN COMPLETE ####\n")

results <- vapply(seq_len(nRep), CompareMethods, matrix(0, 5, 5))
write.csv(results, file = paste0('results-', nTip, '-', nRep, '.csv'))

pdf(file = paste0('results-', nTip, '-', nRep, '.pdf'),
    title = paste0("Results: ", nTip, " tips, ", nRep, " replications"),
    width = 8, height = 6)
par(mfrow = c(2, 2), mgp = c(1.9, 0.7, 0), mar = rep(3, 4))
# par(mfrow = c(2, 3), mgp = c(1.9, 0.7, 0), mar = rep(3, 4))
boxplot(t(results['betterThanGen', , ]),
        ylab = 'Trees better than generative', notch = TRUE, 
        frame = FALSE, col = cols)
boxplot(t(results['betterThanGen', , ] + results['equalToGen', , ]),
        ylab = 'Trees equal or better than generative', 
        frame = FALSE, col = cols)
boxplot(t(results['cidFromGen', , ]), notch = TRUE,
        ylab = 'CID from generative',
        frame = FALSE, col = cols)
# boxplot(t(results['qdFromGen', , ]),
#         ylab = 'QD from generative',
#         frame = FALSE, col = cols)
# boxplot(t(results['tbrFromGen', , ]), notch = TRUE,
#         ylab = 'TBR from generative',
#         frame = FALSE, col = cols)
plot(type = 'n', 0, 0, main = paste0(nTip, ' tips, ', nRep, ' replicates'),
     frame = FALSE, axes = FALSE, xlab = '', ylab = '')
dev.off()
