#' Prepare the input data for LTM analysis
#'
#' This function can do quantile normalization, get the representative profile for rows with duplicated gene symbols,
#' filter out genes with low expression values, and blunt the outliers of the given expression matrix.
#'
#' @param prepD  a data frame. The first column is the gene symbol, and other columns are samples. One row per gene.
#' @param quantNorm logical. Set \code{TRUE} will do quantile normalization across samples.
#' @param uniStyle a character string. Select the way of getting the representative profile for rows with duplicated gene symbols.
#' It must be \code{"mad"}(default), \code{"sum"}, or \code{"none"}, represent 'select the row with max mean absolute deviation', 'sum up duplicate rows'
#' or 'there is no duplicate gene name, and skip the step of getting the representative profile'.
#' @param removeLowQuant a numeric value. Filter out genes based on their median expression value. If the median expression value locate below the \code{0.25}(default) quantile,
#' this gene will be filtered out. Set it to \code{0} for keeping all the input genes.
#' @param bluntLowQuant a numeric value. Blunt the outlier value across samples gene by gene. If an expression value is below the \code{0.025}(default) quantile,
#' this outlier value will be blunted.
#' @param bluntHighQuant a numeric value. Blunt the outlier value across samples gene by gene. If an expression value is above the \code{0.975}(default) quantile,
#' this outlier value will be blunted.
#' @param digitsNum an integer value. The integer indicates the number of decimal places to be used in the prepared data frame.
#' @param outFileSymbol a character string. A prefix exists in the name of output file.
#' Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store output file.
#' @return a data fame prepared for LTM analysis
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom dplyr mutate group_by filter select summarise_if bind_rows
#' @importFrom purrr map

LTMprep <- function(prepD, quantNorm = TRUE, uniStyle = "mad", removeLowQuant = 0.1, bluntLowQuant = 0.025, bluntHighQuant = 0.975, digitsNum = 4, outFileSymbol = NULL, outDir = NULL) {
  ##the first column should be the gene symbol or other annotation name, and other columns should be expression values
  inputD = as.data.frame(prepD, stringsAsFactors = FALSE)
  colnames(inputD)[1] = "geneSym"
  ##quantile normalization
  if ( quantNorm ) {
    normD = preprocessCore::normalize.quantiles(as.matrix(log2(inputD[,-1] + 10**(-1*digitsNum))))   ##set 1e-10 is not good, reset this as a parameter, related it to digitsNum
    normD = as.data.frame(2^normD) %>%
      dplyr::mutate(geneSym = inputD$geneSym)
    colnames(normD) = c(colnames(inputD)[-1], "geneSym")
    normD = normD[,colnames(inputD)]
    inputD = normD
  }
  ##remove redundant symbols
  if (tolower(uniStyle) == "sum") {
    inputD = dplyr::group_by(inputD, geneSym)  %>%
      dplyr::mutate( sdv = apply(inputD[,-1], 1, sd) ) %>%
      dplyr::filter(geneSym != "---" & geneSym != "" & !is.na(geneSym), sdv != 0) %>%
      dplyr::select(-sdv) %>%
      dplyr::mutate( rowNum = length(geneSym) ) %>%
      as.data.frame()
    uniD = dplyr::filter(inputD, rowNum == 1)
    multiD = dplyr::filter(inputD, rowNum > 1) %>%
      split(.$geneSym) %>%
      purrr::map( function(zD) {
        zoutD = zD %>%
          dplyr::summarise_if(is.numeric, sum, na.rm = TRUE) %>%
          dplyr::mutate(geneSym = unique(zD$geneSym))
        return(zoutD)
      } ) %>% bind_rows() %>%
      as.data.frame()
    outD = dplyr::bind_rows(uniD, multiD[,colnames(uniD)]) %>%
      dplyr::select(-rowNum)
  }  else if (tolower(uniStyle) == "mad") {
    outD =  dplyr::mutate(inputD, sdv=apply(inputD[,-1], 1, sd), madv=apply(inputD[,-1], 1, mad), index=1:nrow(inputD) ) %>%
      dplyr::filter(geneSym != "---" & geneSym != "" & !is.na(geneSym), sdv != 0) %>%
      dplyr::group_by(geneSym) %>%
      dplyr::filter(madv == max(madv))  %>%
      dplyr::filter(index == min(index)) %>%
      dplyr::select(-sdv, -madv, -index)
  }  else if (tolower(uniStyle) == "none")  {
    outD = inputD
  }  else  {
    cat("Error:please set uniStyle to 'mad', 'sum' or 'none'.\n")
  }
  ##remove values in the low quantile
  if ( (removeLowQuant > 0) & (removeLowQuant < 1) ) {
    medianv = apply(outD[,-1], 1, median)
    low_cut = quantile(medianv, probs = c(0, removeLowQuant, 1))
    row_index = which(medianv > as.numeric(low_cut[2]) )
    outD = outD[row_index,]
  }
  ##blunt low outliers
  if ((bluntLowQuant > 0) & (bluntLowQuant < 1) )  {
    outB = apply(outD[,-1], 1, function(z)  {
      lowv = as.numeric( quantile(z, probs = c(0, bluntLowQuant, 1))[2] )
      z[z < lowv] = lowv
      return(round(z, digitsNum))
    } )
    ##transform the data
    outB = as.data.frame(t(outB))
    colnames(outB) = colnames(outD)[-1]
    outB = dplyr::mutate(outB, geneSym = outD$geneSym)
    outD = outB[,colnames(outD)]
  }
  ##blunt high outliers
  if ( (bluntHighQuant < 1) & (bluntHighQuant > 0) )  {
    outB = apply(outD[,-1], 1, function(z)  {
      highv = as.numeric( quantile(z, probs = c(0, bluntHighQuant, 1))[2] )
      z[z > highv] = highv
      return(round(z, digitsNum))
    } )
    ##transform the data
    outB = as.data.frame(t(outB))
    colnames(outB) = colnames(outD)[-1]
    outB = dplyr::mutate(outB, geneSym = outD$geneSym)
    outD = outB[,colnames(outD)]
  }
  ##check and remove those rows with constant values again
  outD = dplyr::mutate(outD, sdv = apply(outD[,-1], 1, sd))  %>%
    dplyr::filter(sdv != 0) %>%
    dplyr::select(-sdv)
  ##output the file
  if (length(outFileSymbol)) {
    if (length(outDir)) {
      write.csv( outD, file = paste0(outDir, "/", outFileSymbol, ".LTMprep.csv"), row.names = FALSE )
    }  else  {
      write.csv( outD, file = paste0(outFileSymbol, ".LTMprep.csv"), row.names = FALSE )
    }
  }
  if (!length(outFileSymbol)) {
    return( outD )
  }
}


#' Prepare the input data for a quick run of LTMheat
#'
#' This function will select genes with high expression values and large expression variations across samples for a quick run of LTMheat.
#'
#' @param cutD  a data frame. The first column is the gene symbol, and other columns are samples. One row per gene.
#' @param exemptGenes a character vector. A list of clock genes that will be kept in the output data frame without getting into the filtration process.
#' @param minExp a numeric value. The minimal median expression value in quantile group 1.
#' @param minFold a numeric value. The minimal fold change between quantile group 4 and quantile group 1. Median value in each quantile group is used.
#' @param outFileSymbol a character string. A prefix exists in the name of output file.
#' Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store output file.
#' @return a data fame prepared for a quick run of LTM
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom dplyr filter

LTMcut <- function(cutD = NULL, exemptGenes = NULL, minExp = -1e10, minFold = 1, outFileSymbol = NULL, outDir = NULL) {
  ##a data frame; the first column is gene name, other columns are numeric values
  colnames(cutD)[1] = "geneSym"
  testD = dplyr::filter(cutD, !(geneSym %in% exemptGenes))
  exemptD = dplyr::filter(cutD, geneSym %in% exemptGenes)
  outindex = apply(testD[,-1], 1, function(inputv) {
    ##cut-off: Grp[Q75-Q100] / Grp[Q0-Q25] >= minfold; Grp[Q0-Q25] > minexp; use the middle of each group
    qexpv = quantile(inputv, probs = c(0.125, 0.875), na.rm = TRUE)
    if ( (qexpv[1] >= minExp) & (qexpv[2] / qexpv[1] >= minFold) ) {
      return(1)
    } else {
      return(0)
    }
  } )
  outD = bind_rows(testD[which(outindex > 0),], exemptD)
  ##output the file
  if (length(outFileSymbol)) {
    if (length(outDir)) {
      write.csv( outD, file = paste0(outDir, "/", outFileSymbol, ".LTMcut.csv"), row.names = FALSE )
    }  else  {
      write.csv( outD, file = paste0(outFileSymbol, ".LTMcut.csv"), row.names = FALSE )
    }
  }
  if (!length(outFileSymbol)) {
    return( outD )
  }
}


#' Prepare the screen and sample data used for running LTMcook
#'
#' This function will separate the samples into given number of quantile groups, calculate the mantel's zstat value,
#' nCV values of given clock gene list, gene by gene. It will output a list with screen and sample data frame used for running LTMcook.
#'
#' @param heatD  a data frame. The first column is the gene symbol, and other columns are samples. One row per gene.
#' @param benchD a data frame. The expression correlation values of paired clock genes.
#' @param cvGenes a vector string. The list of clock genes representing the oscillation robustness of circadian clock.
#' @param outFileSymbol a character string. A prefix exists in the name of output file.
#' Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store output file.
#' @param qnum a numeric value. The number of quantile groups when separate samples based on expression level gene by gene.
#' @param nCores an integer. Number of cores to use. For parallel computing, set \code{nCores} a integer larger than 1.
#' @param logical. If TRUE, time used during the analysis will be released on the screen.
#' @return a list with screen and sample data frame.
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom dplyr mutate group_by filter select inner_join bind_rows
#' @importFrom parallel mclapply

LTMheat <- function(heatD, benchD, cvGenes, outFileSymbol = NULL, outDir = NULL, qnum = 4, nCores = 1, releaseNote = TRUE) {
  ###start time
  run_start=proc.time()
  ###check the bench data frame
  colnames(benchD)[1] = "geneSym"
  if ( nrow(benchD) != length(unique(benchD$geneSym)) ) {
    stop( c("Please make sure that there is no duplciate names in the first column of 'benchD'.\n") );
  }
  rownames(benchD) = benchD$geneSym
  ##the output data and sample data
  dataD = as.data.frame(heatD)
  colnames(dataD)[1] = "geneSym"
  if ( nrow(dataD) != length(unique(dataD$geneSym)) ) {
    stop( c("Please make sure that there is no duplciate names in the first column of each data frame of the 'heatD' data frame list.\n") );
  }
  rownames(dataD) = dataD$geneSym
  cnames = colnames(dataD)
  commonD = dplyr::inner_join(data.frame(geneSym = benchD$geneSym, stringsAsFactors = FALSE), dataD, by = "geneSym")
  ##run gene one by one
  indexL = apply(dataD, 1, function(za) {
    zout = data.frame(expv = as.numeric(za[-1]), index = 2:length(za)) %>%
      dplyr::mutate(tag = paste0("Q", dplyr::ntile(expv, qnum)), gname = unlist(za[1]) ) %>%
      dplyr::group_by(tag) %>%
      dplyr::mutate(meanv = mean(expv)) %>%
      dplyr::select(-expv) %>%
      as.data.frame() %>%
      split(.$tag)
    return(zout)
  } )
  ##the mantel test, cv calculation
  if (nCores > 1)  {
    ##the mantel test, cv calculation, parallel version
    outputD = parallel::mclapply(indexL, runMantelF, expD = dataD, benchCor = benchD,
                                 gorder = commonD$geneSym, cvGeneV = cvGenes, mc.cores = nCores) %>% dplyr::bind_rows()
  } else {
    outputD = lapply(indexL, runMantelF, expD = dataD, benchCor = benchD,
                     gorder = commonD$geneSym, cvGeneV = cvGenes) %>% dplyr::bind_rows()
  }
  ###get the sample list
  sampleD = lapply(indexL, function(zinL) {
    zoutD = lapply(zinL, function(zD) {
      zindex = zD$index
      zsampleD = zD %>% dplyr::mutate(sampleID = cnames[zindex]) %>%
        dplyr::select(-index, -meanv)
    }) %>% dplyr::bind_rows()
    return(zoutD)
  } ) %>% dplyr::bind_rows()
  ##output the files
  if (length(outFileSymbol)) {
    if (length(outDir)) {
      write.csv( outputD, file = paste0(outDir, "/", outFileSymbol, ".Q", qnum, ".screen.csv"), row.names = FALSE )
      write.csv( sampleD, file = paste0(outDir, "/", outFileSymbol, ".Q", qnum, ".sample.csv"), row.names = FALSE )
    }  else  {
      write.csv( outputD, file = paste0(outFileSymbol, ".Q", qnum, ".screen.csv"), row.names = FALSE )
      write.csv( sampleD, file = paste0(outFileSymbol, ".Q", qnum, ".sample.csv"), row.names = FALSE )
    }
  }
  if (releaseNote)  {
    cat(c("The LTMheat analysis is done.\n"));
    print( c("Time used:", (proc.time()-run_start)) );
    cat("\n\n");
  }
  if (!length(outFileSymbol)) {
    return( list("screen" = outputD, "sample" = sampleD) )
  }
}


#' Prepare the data for running LTMdish
#'
#' This function calculates the correlation value between expression level and measures of clock strength, gene by gene.
#' The input screen data is from LTMheat.
#'
#' @param cookD  a data frame. The screen data from LTMheat.
#' @param corMethod a character string. Select the correlation coefficient is to be used for
#' the correlation test between expression level and measures of clock strength.
#' It must be \code{"pearson"}(default), \code{"kendall"}, or \code{"spearman"}.
#' @param outFileSymbol a character string. A prefix exists in the name of output file.
#' Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store output file.
#' @return a data fame
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tidyr gather
#' @importFrom purrr map

LTMcook <- function(cookD, corMethod = "pearson", outFileSymbol = NULL, outDir = NULL) {
  if ( (length(corMethod) == 1) & (corMethod %in% c("pearson", "kendall", "spearman")) ) {
    ##get the correlation data frame
    screenD = as.data.frame(cookD)
    colnames(screenD)[1] = "geneSym"
    screenCorD = screenD %>%
      split(.$geneSym) %>% purrr::map( function(zD) {
        zDC = zD %>% tidyr::gather(key = "measure", value = "zva", -geneSym, -tag, -zmean) %>%
          split(.$measure) %>% purrr::map( function(zaD) {
            zaout = data.frame(geneSym = unique(zaD$geneSym),
                               measure = unique(zaD$measure),
                               corv = cor(zaD$zmean, zaD$zva, method = corMethod),
                               p.value = cor.test(zaD$zmean, zaD$zva, method = corMethod)["p.value"], stringsAsFactors = FALSE)
            return(zaout)
          } ) %>% dplyr::bind_rows()
      } )  %>% dplyr::bind_rows()
    ##output the correlation files
    if (length(outFileSymbol)) {
      if (length(outDir)) {
        write.csv( screenCorD, file = paste0(outDir, "/", outFileSymbol, ".screenCorv.csv"), row.names = FALSE )
      }  else  {
        write.csv( screenCorD, file = paste0(outFileSymbol, ".screenCorv.csv"), row.names = FALSE )
      }
    } else {
      return(screenCorD)
    }
  } else {
    cat("Please set 'corMethod' as 'pearson', 'kendall', or 'spearman'.\n")
  }
}


#' Generate the LTM analysis result for the input expression dataset
#'
#' This function calculates the LTM_abs value for each gene, which indicates the correlation between gene expression and clock strength.
#'
#' @param dishL  a list of data frame. Each data frame is from LTMcook output.
#' @param targetMeasures a vector of character string. The target measures of clock strength used to calculate LTM_abs value for each gene.
#' @param onlyConsistentTrend logical. Set \code{TRUE} will only keep those genes show consistent trend across targetMeasures.
#' @param fileName a character string. The name of output file. Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store output file.
#' @return a data fame
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom dplyr mutate group_by select bind_rows  arrange
#' @importFrom purrr map

LTMdish <- function( dishL, targetMeasures = c("zmantel", "zncv"), onlyConsistentTrend = FALSE, fileName = NULL, outDir = NULL)  {
  ###get the LTM_pre value
  outD = purrr::map( dishL, function( screenCorD ) {
    sigCorD = topTabF(screen_CorD = screenCorD, target_Measures = targetMeasures,
                      only_ConsistentTrend = onlyConsistentTrend)
    return(sigCorD)
  }) %>% dplyr::bind_rows()
  ###get the LTM_ori and LTM_abs value
  outD = outD %>%
    dplyr::group_by( geneSym ) %>%
    dplyr::summarise( LTM_ori = mean(LTM_pre),
                   #LTM_pvalue = max(p.adjust(maxPva, method = adjMethod)),
                   LTM_pvalue = max(maxPva) ) %>%
    dplyr::mutate( LTM_abs = abs(LTM_ori) ) %>%
    as.data.frame() %>%
    dplyr::mutate( LTM_BH.Q = p.adjust(LTM_pvalue, method = "BH") ) %>%
    #dplyr::select( geneSym, LTM_ori, LTM_abs, LTM_pvalue, LTM_BH.Q)
    dplyr::select( geneSym, LTM_ori, LTM_abs) %>%
    dplyr::arrange( -LTM_abs )
  ##output the file or return a data frame
  if (length(fileName)) {
    if (length(outDir)) {
      write.csv( outD, file = paste0(outDir, "/", fileName), row.names = FALSE )
    }  else  {
      write.csv( outD, file = paste0(fileName), row.names = FALSE )
    }
  } else {
    return(outD)
  }
}


#' Integrate the LTM analysis results of multiple datasets.
#'
#' This function get the meta values by integrating LTM analysis results of multiple datasets, e.g., multiple skin datasets.
#'
#' @param metaL  a list of data frame. Each data frame is from LTMdish output.
#' @param fileName a character string. The name of output file. Set it \code{NULL}(default) if prefer to return a data frame.
#' @param outDir a character string. The name of directory used to store the output file.
#' @return a data fame
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export
#' @importFrom dplyr mutate group_by filter select summarise bind_rows arrange
#' @importFrom purrr map2

LTMmeta <- function(metaL, fileName = NULL, outDir = NULL) {
  tags = paste("meta", 1:length(metaL))
  metaD = purrr::map2(metaL, tags, function(ztepD = .x, tag = .y) {
    colnames(ztepD)[1] = "geneSym"
    ztepD = dplyr::mutate(ztepD, grp = tag)
    return(ztepD)
  } ) %>% dplyr::bind_rows()
  ###integrate result
  combD = metaD %>%
    dplyr::filter( !is.na(geneSym) ) %>%
    dplyr::group_by(geneSym) %>%
    dplyr::summarise( meta_LTM_ori = sum(LTM_ori)/length(tags) ) %>%
    as.data.frame() %>%
    dplyr::mutate( meta_LTM_abs = abs(meta_LTM_ori) ) %>%
    dplyr::select(geneSym, meta_LTM_ori, meta_LTM_abs) %>%
    dplyr::arrange( -meta_LTM_abs )
  # ##get the integrated p-value
  # if (tolower(combinePvalue) == "fisher") {
  #   combD = metaD %>% group_by(geneSym) %>%
  #     dplyr::summarise( meta_LTM_ori = sum(LTM_ori)/length(tags),
  #                       meta_LTM_pvalue = fisher.method(LTM_pvalue) ) %>%
  #     as.data.frame()
  # } else if (tolower(combinePvalue) == "max")  {
  #   combD = metaD %>% group_by(geneSym) %>%
  #     dplyr::summarise( meta_LTM_ori = sum(LTM_ori)/length(tags),
  #                       meta_LTM_pvalue = max(LTM_pvalue) ) %>%
  #     as.data.frame()
  # } else if (tolower(combinePvalue) == "min") {
  #   combD = metaD %>% group_by(geneSym) %>%
  #     dplyr::summarise( meta_LTM_ori = sum(LTM_ori)/length(tags),
  #                       meta_LTM_pvalue = min(LTM_pvalue) ) %>%
  #     as.data.frame()
  # } else if (tolower(combinePvalue) == "none") {
  #   combD = metaD %>% group_by(geneSym) %>%
  #     dplyr::summarise( meta_LTM_ori = sum(LTM_ori)/length(tags),
  #                       meta_LTM_pvalue = paste(LTM_pvalue, collapse = "///") ) %>%
  #     as.data.frame()
  # } else {
  #   stop( c("Please set a value to 'combinePvalue' as 'fisher', 'max', 'min', or 'none'.\n") );
  # }
  # combD = combD %>%
  #   dplyr::mutate( meta_LTM_abs = abs(meta_LTM_ori),
  #                  meta_LTM_BH.Q = p.adjust(meta_LTM_pvalue, method = "BH") ) %>%
  #   dplyr::select(geneSym, meta_LTM_ori, meta_LTM_abs, meta_LTM_pvalue, meta_LTM_BH.Q) %>%
  #   dplyr::arrange( -meta_LTM_abs )
  ##output the file or return a data frame
  if (length(fileName)) {
    if (length(outDir)) {
      write.csv( combD, file = paste0(outDir, "/", fileName), row.names = FALSE )
    }  else  {
      write.csv( combD, file = paste0(fileName), row.names = FALSE )
    }
  } else {
    return(combD)
  }
}


#' Delete unnecessary files
#'
#' This function helps delete intermediate files generated during LTM analysis.
#'
#' @param washDir a character string. The name of directory used to store intermediate files during LTM analysis.
#' @param washKey a character string. The shared key word among targeted files under deletion.
#' @param checkWashList logical. Set \code{TRUE} will return a list of files.
#' @return a data frame of deleted files.
#' @examples
#' ## please refer to the webpage of LTMR package
#' @export

######LTMwash(delete or compress temporary files)
LTMwash <- function(washDir, washKey = NULL, checkWashList = FALSE) {
  delFiles = grep(washKey, dir(washDir), value = TRUE)
  delFiles = paste0(washDir, "/", delFiles)
  if (length(delFiles) & (washWay == "delete") ) {
    zflag = file.remove(delFiles)
  }
  if (checkWashList) {
    return(cbind(delFiles, zflag))
  }
}
