######==========================================================================================================
##sub-function of mantel test; make sure inputD and benchCor with the same gene order
mantelF <- function(inputD, benchCor, permutationValue = 1, digitsNum = 3)  {
  ##the first column of inputD is gene identification
  corA <- stats::cor(t(inputD[,-1]), method = "spearman")
  simA <- ape::mantel.test(corA, benchCor, nperm = permutationValue)  ##calculate the similarity between matrix with 'ape' package
  return(round(simA$z.stat, digitsNum))
}
######==========================================================================================================
##subfunction of running mantel test for each gene; zindexL is one gene with multiple sub-groups
runMantelF <- function(zindexL, expD, benchCor, gorder, cvGeneV) {
  benchCor = as.matrix(benchCor[gorder,gorder])
  zstatD = lapply(zindexL, function(zD) {
    ##make sure no sd = 0 for any gene before running mantelF
    mantChekA = apply(expD[gorder,zD$index], 1, sd)
    if (length(which(mantChekA == 0))) {
      zstat = NA
      zcvD = data.frame(nCV = NA)
    } else {
      ##calculate the zstat
      zstat = mantelF(inputD = expD[gorder,c(1, zD$index)], benchCor = benchCor)
      ##calculate cv for all genes and get the nCV for target genes
      zsdv = apply(expD[,zD$index], 1, sd)
      zmeanv = apply(expD[,zD$index], 1, mean)
      zcvD = data.frame("geneSym" = unlist(expD[,1]), "cva" = zsdv / zmeanv, stringsAsFactors = FALSE) %>%
        dplyr::mutate(nCV = cva / mean(cva)) %>%
        dplyr::filter(geneSym %in% cvGeneV) %>%
        dplyr::select(geneSym, nCV)
      # ##get the nCV for target genes using nCV package (meet the aim of nCV test, but not efficient enough for LTM test)
      # zcvD = nCV::nCVgene(inputD = expD[, c(1, zD$index)], cgenes = cvGeneV)
    }
    zout = data.frame(gname = unique(zD$gname), tag = unique(zD$tag), zmean = unique(zD$meanv),
                      zmantel = zstat, zncv = mean(zcvD$nCV), stringsAsFactors = FALSE) %>%
      bind_cols(tidyr::spread(zcvD, key = "geneSym", value = "nCV"))
    return(zout)
  }) %>% bind_rows()
  return(zstatD)
}

######==========================================================================================================
######get the correlation value
topTabF <- function(screen_CorD = NULL, target_Measures = NULL, only_ConsistentTrend = FALSE) {
  colnames(screen_CorD)[1] = "geneSym"
  topD = screen_CorD %>%
    dplyr::filter( measure %in% target_Measures ) %>%
    base::split(.$geneSym) %>%
    purrr::map( function(zD) {
      ###only keep consistent positive or negative trend
      if (only_ConsistentTrend) {
        zoutD = zD %>% dplyr::mutate(trendNum = length(unique(sign(corv))),
                                     uniTrend = paste(unique(sign(corv)), collapse = "//"),
                                     #weightCor = prod(corv),
                                     LTM_pre = mean(corv),
                                     #fisherPva = fisher.method(p.value),
                                     ##the max as the representative p-vlaue of nCV and mantel test
                                     maxPva = max(p.value) ) %>%
          #dplyr::filter(trendNum == 1, abs(weightCor) >= weit_CorCut) %>%
          dplyr::filter(trendNum == 1) %>%
          #dplyr::arrange(uniTrend, weightCor) %>%
          dplyr::select(-trendNum)
      }
      ###keep all
      if (!only_ConsistentTrend) {
        zoutD = zD %>% dplyr::mutate(uniTrend = paste(sign(corv), collapse = "//"),
                                     LTM_pre = mean(corv),
                                     #fisherPva = fisher.method(p.value),
                                     maxPva = max(p.value) ) #%>%
          #dplyr::filter(abs(weightCor) >= weit_CorCut) %>%
          #dplyr::arrange(uniTrend, LTM_pre) %>%
          #dplyr::select(-trendNum)
      }
      return(zoutD)
    } )  %>% dplyr::bind_rows() %>%
    dplyr::select(-p.value) %>%
    tidyr::spread(key = measure, value = corv) %>%
    dplyr::arrange(uniTrend, LTM_pre)
  return(topD[,c("geneSym", "uniTrend", "LTM_pre", "maxPva", target_Measures)])
}
######==========================================================================================================
######combine p-values by performing Fisher's method (source code in MetaCycle)
fisher.method <- function(pvals, zero.sub=1e-100){
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  pvals <- pvals[!is.na(pvals)]
  S = -2*sum(log(pvals))
  num.p = length(pvals)
  fisher.pva <- 1 - pchisq(S, df=2*num.p)
  return(fisher.pva)
}
######==========================================================================================================
tcgaTN <- function(inputD, outputTN = "tumor") {
  ###the first column should be the gene symbol or other annotation name, and other columns should be expression values
  colnames(inputD)[1] <- "Gene.Symbol"
  ##get the tumor samples(outputTN = "tumor")/normal samples(outputTN = "normal")/tumor-normal matched( outputTN = "paired")/all samples( outputTN = "all")
  if (tolower(outputTN) != "all") {
    sampleD = data.frame(sampleID = colnames(inputD)[-1], stringsAsFactors = FALSE) %>%
      dplyr::mutate( subjectID = gsub("(TCGA\\.\\S+\\.\\S+)\\.\\d\\d\\w\\.\\d\\d\\w\\.\\S+", "\\1", sampleID),
                     group = gsub("TCGA\\.\\S+\\.\\S+\\.(\\d\\d)\\w\\.\\d\\d\\w\\.\\S+", "\\1", sampleID) ) %>%
      dplyr::mutate( group = as.numeric(group) )
    ##tumors with 01-09; normal with 10-19; control with 20-29; take normal and control as the same group
    if (tolower(outputTN) == "tumor") {
      sampleD = dplyr::filter(sampleD, group < 10)
    }  else if (tolower(outputTN) == "normal") {
      sampleD = dplyr::filter(sampleD, group >= 10)
    } else if (tolower(outputTN) == "paired") {
      sampleD = sampleD %>% dplyr::mutate( tag = case_when(group < 10 ~ -1,
                                                           group >=10 ~ 1) ) %>%
        group_by(subjectID) %>%
        dplyr::mutate( numv = length(sampleID), flag = sum(tag) ) %>%
        dplyr::filter( numv == 2, flag == 0) %>%
        as.data.frame()
    } else {
      sampleD = NULL
      cat("Error: wrong parameter values for outputTN!!!\n")
    }
    inputD = inputD[,c("Gene.Symbol", sampleD$sampleID)]
  }
  return(inputD)
}
######==========================================================================================================
######transform the gene symbols to ensembl ID
topEnsemF <- function(inDF, species = "Hs") {
  colnames(inDF)[1] = "gname"
  topD = indDF %>%
    dplyr::distinct() %>%
    dplyr::mutate( geneSym = limma::alias2SymbolTable(gname, species = species) )
  ###transform gene symbol to ENSEMBL_GENE_ID, which is much efficient for running DAVID on line
  #listDatasets(useEnsembl(biomart="ensembl"))
  if (species == "Hs") {
    ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    #listAttributes(ensembl)
    ### minor revision here, replace 'topD$gname' with 'topD$geneSym'
    ensemD = biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters = 'hgnc_symbol',
                            values = topD$geneSym, mart = ensembl, uniqueRows = TRUE) %>%
      dplyr::group_by(hgnc_symbol) %>%
      dplyr::mutate(uniGeneID = paste(ensembl_gene_id, collapse = "_")) %>%
      as.data.frame() %>%
      dplyr::select(-ensembl_gene_id) %>%
      dplyr::distinct() %>%
    ### minor revision here, replace 'c("hgnc_symbol" = "gname")' with 'c("hgnc_symbol" = "geneSym")'
      dplyr::right_join(topD, by = c("hgnc_symbol" = "geneSym") ) %>%
      dplyr::filter(!is.na(uniGeneID))
  }
  if (species == "Mm") {
    ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
    ensemD = biomaRt::getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), filters = 'mgi_symbol',
                   values = topD$gname, mart = ensembl, uniqueRows = TRUE) %>%
      dplyr::group_by(mgi_symbol) %>%
      dplyr::mutate(uniGeneID = paste(ensembl_gene_id, collapse = "_")) %>%
      as.data.frame() %>%
      dplyr::select(-ensembl_gene_id) %>%
      dplyr::distinct() %>%
      dplyr::right_join(topD, by = c("mgi_symbol" = "gname") ) %>%
      dplyr::filter(!is.na(uniGeneID))
  }
  return(ensemD)
}
