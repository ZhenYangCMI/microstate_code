

rm (list=ls(all= TRUE))

library(aplpack)
library(R.matlab)
seed <-4
state <-4
base <- "/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/"
fileName <-file.path(base,sprintf("duration_seed%dState%d_2ses.mat",seed,state))
matData <- readMat(fileName)
varName <- sprintf("seed%dState%d", seed, state)
duration <-matData[[varName]]

figDir <-"/Users/zhenyang/Desktop/Zhen/figs/5_3_13"
figName <-file.path(figDir,sprintf("bagplot_seed%dState%d", seed, state))
png(filename=figName, width=480)
bagplot(duration[,1],duration[,2], xlab='session1', ylab="session2", main=varName, pch = 16, cex = 1)
dev.off()
