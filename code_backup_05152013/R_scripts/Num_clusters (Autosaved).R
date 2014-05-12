
# Summarize the number of clusters
rm(list = ls(all = TRUE))

sessionType <- c('session1','session2')
sessionNum <- 1:length(sessionType)
winType <- c("Full","Partial")
winTypeNum <- 1:length(winType)

for (i in sessionNum) {
	for (j in winTypeNum) {
base <- "/Users/zhen.yang/Documents/Zhen_CMI/microstate/results/lambdaOptimalPerSub"
fileName <- sprintf("/Users/zhen.yang/Documents/Zhen_CMI/microstate/results/lambdaOptimalPerSub/winAllSeeds_%sCorLasso_%s.txt", winType[j], sessionType[i])	
indxFile <- read.table(fileName)
table(indxFile)
}
}