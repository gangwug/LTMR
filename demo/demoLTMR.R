######==========================run the demo using the example data
rm(list=ls())
library(LTMR)
######load the example data
load("../data/mClockBenchD.rda")
load("../data/exampleD.rda")
######LTMprep: normalize the data
outPrepD = LTMprep(prepD = exampleD, quantNorm = TRUE, uniStyle = "mad", removeLowQuant = 0.1, bluntLowQuant = 0.025, bluntHighQuant = 0.975)
######LTMcut: filter the data for a quick test
outCutD = LTMcut(cutD = outPrepD, exemptGenes = mClockBenchD$geneSym, minExp = 10, minFold = 1.5)
######LTMheat: prepare the input data for latter analysis / takes ~half an hour per run / parallel computing will use less time when set 'nCores' larger than 1
outHeatA = LTMheat(heatD = outCutD, benchD = mClockBenchD, cvGenes = mClockBenchD$geneSym, qnum = 4, nCores = 1, releaseNote = TRUE)
outHeatA = outHeatA$screen
outHeatB = LTMheat(heatD = outCutD, benchD = mClockBenchD, cvGenes = mClockBenchD$geneSym, qnum = 5, nCores = 1, releaseNote = TRUE)
outHeatB = outHeatB$screen
######LTMcook: do the correlation for each quantile group
outCookA = LTMcook(cookD = outHeatA, corMethod = "pearson")
outCookB = LTMcook(cookD = outHeatB, corMethod = "pearson")
######LTMdish: get the LTM output results
LTMdish(dishL = list("Q4" = outCookA, "Q5" = outCookB), targetMeasures = c("zmantel", "zncv"), fileName = "example_LTMdish.csv", outDir = "./" )
#####take a look at the output file
