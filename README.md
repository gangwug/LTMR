## Introduction about LTMR package
LTM is an in silico screen to infer genetic influences on circadian clock function. LTM uses natural variation in gene expression data and directly links gene expression variation to clock strength independent of longitudinal data.

## Installation
Use **devtools** to install this version from Github:

  ```r
# install 'devtools' in R(>3.0.2)
# install.packages("devtools")
# install LTMR
devtools::install_github('gangwug/LTMR')
```

## Usage
```r
library(LTMR)
###### follow the below script to run the example data
######load the example data
head(mClockBenchD)
head(exampleD)
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
```
## Advice for LTMR users
1. LTM screened genes in the complex tissues (eg., tumors) may pick other stronger factors (eg., fractions of non-cancer cells) than clock-coupled pathways. Considering multiple factors may influence the clock strength in human population datasets, non-biological factors (eg., batch effects) may also affect the interpretation of LTM screening results. It is important to know the dataset before applying LTM. 
2. The genes identified by LTM are based on correlations. The strong correlation is not the direct evidence of causal effect. Therefore, the LTM predicted clock-coupled pathways need experimental validation and/or literature evidence. 
3. LTM can not screen genes correlated with the period and phase variation of the circadian clock. 
4. Ranking samples by the gene expression are strongly influenced by the circadian phase for strong cyclers, which makes it difficult to accurately quantify the clock strength when separating samples by circadian phase. So LTM may be biased to non-cycling genes.
5. The intesntiy values are suggested for array data. For RNA-seq data, we suggest to use read-normalized values (e.g., TPM, FPKM or RPKM; TPM is preferred) instead of the raw read counts. 

## More information

Wu G, Ruben MD, Francey LJ,  Lee Y, Anafi RC, Hogenesch JB. An in silico genome-wide screen for circadian clock strength in human samples. bioRxiv, 2022, https://www.biorxiv.org/content/10.1101/2022.05.10.491250v2.
