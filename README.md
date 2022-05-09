# LTMR
An in silico genome-wide screen for circadian clock strength

## Introduction about nCV package
Robust oscillation of clock genes is a core feature of the circadian system. Relative amplitude ([rAMP](https://journals.sagepub.com/doi/10.26599/BSA.2020.9050005)) measures the robustness of clock gene oscillations, but only works for longitudinal samples. We lack a method for estimating robust oscillations from human samples without labeled time. Therefore, we developed the normalized coefficient of variation ([nCV](https://www.biorxiv.org/content/10.1101/2021.07.28.454045v1.full)) method to address this challenge, and implement it into an R package. 

The nCV package has two functions: nCVnet and nCVgene. nCVnet can test whether there is a functional clock network in population scale data. nCVgene can evaluate the robustness of clock genes in population scale data. 

## Installation
Use **devtools** to install this version from Github:

  ```r
# install 'devtools' in R(>3.0.2)
# install.packages("devtools")
# install nCV
devtools::install_github('gangwug/nCV')
```

## Usage
```r
library(nCV)
# see detail introduction and associated application examples
?nCVnet
?nCVgene
```
## Advice for nCV users

### Transcriptome data format

The intesntiy values are suggested for array data. For RNA-seq data, we suggest to use read-normalized values (e.g., TPM, FPKM or RPKM; TPM is preferred) instead of the raw read counts. 

### Requirements of gene expression in the data for nCV calculation

Before applying nCV package, please make sure 1) there is not obvious batch effects (e.g., biopsy or sequenciing center bias) in the population scale data, 2) remove genes with zero variance (i.e., same expression value for all the samples), 3) majority clock genes are deteced by the array or RNA-seq. 

### The number of samples in the population scale data

The test using mouse liver data suggests 1) the linear relationship between nCV and nrAMP keeps even only containing samples collected in a 12 hour window, 2) the linear correlation value becomes larger if containing samples across the whole circle. Since majority human samples do not have sample collection time, including a large number of samples is the only way of helping guarantee the nCV calculation with samples across the circle. Our current experience suggests that the minimal number of samples suitable for nCV calculation is 50. For human population data analysis, more samples are always preferred. 

### First run nCVnet, then run nCVgene

The nCV only indicate clock robustness if nCVnet returns a significant p-value. 

### The list of target clock genes for nCV calculation

The 'nCVgene' function can calculate the nCV of all genes in the input dataset. However, 'nCVgene' will only output the nCV for the targeted output genes set by the parameter of 'cgenes'. The suggested 10 clock genes (ARNTL, CLOCK, NPAS2, NR1D1, CIART, DBP, PER1, CRY2, PER2, and CRY1) vary in molecular clock phase and robustness level. The nCV may also indicate the oscillation robustness of other ubiquious cyclers (e.g., Usp2, Tsc22d3). 

### Considering the relative large individual difference in human samples caused by many factors (e.g., genetic background, age, sex, nutrition etc.), we suggest to reach the consistent result by testing multiple datasets, especially for estimating the robustness of a single clock gene. For example, we reach the conclusion that the robustness of BMAL1 and PER2 are decreased in tumor samples compared to non-tumor samples based on 8 different datasets. Without testing multiple thyroid carcinoma datasets, we can not confidently say that the robustness of NR1D1 is decreased in tumor samples (possibly true) from patients with thyroid carcinoma yet. 

## More information

Wu G, Francey LJ, Ruben MD, Hogenesch JB. Normalized coefficient of variation (nCV): a method to evaluate circadian clock robustness in population scale data. Bioinformatics, 2021, https://academic.oup.com/bioinformatics/article/37/23/4581/6415819.
