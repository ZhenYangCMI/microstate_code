

rm (list=ls(all= TRUE))

library(aplpack)
library(R.matlab)
seed <- 2
stateT <-2
stateT1 <-4
base <- "/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate"
figDir <-"/Users/zhenyang/Desktop/Zhen/figs/5_13_13"

#plot raw data
fileName <-file.path(base,sprintf("TM_seed%dState%dToState%d.mat",seed,stateT, stateT1))
matData <- readMat(fileName)
varName <- "sesCmb"
TM <-matData[[varName]]

figName <-file.path(figDir,sprintf("TM_seed%dState%dToState%d",seed,stateT,stateT1))

png(filename=figName, width=480)
bagplot(TM[,1],TM[,2], xlab='session1', ylab="session2", main=varName, pch = 16, cex = 1)
dev.off()

