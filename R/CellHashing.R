#' @include LabKeySettings.R
#' @import Seurat
#' @import Rlabkey




#' @title ProcessCiteSeqCount
#'
#' @description A description
#' @param bFile, A
#' @param doRowFilter, A
#' @return A
#' @keywords CITE-Seq,
#' @export
#' @importFrom knitr kable
#' @importFrom knitr kable
ProcessCiteSeqCount <- function(bFile=NA, doRowFilter = T) {
  if (is.na(bFile)){
    stop("No file set: change bFile")
  }

  if (!file.exists(bFile)){
    stop(paste0("File does not exist: ", bFile))
  }

  if (dir.exists(bFile)) {
    #CITE-seq-Count 1.4.2 and higher creates a folder
    bData <- Read10X(bFile, gene.column=1)
    bData <- bData[which(!(rownames(bData) %in% c('unmapped'))), , drop = F]
    bData <- as.matrix(bData)
  } else {
    # older versions create a CSV file
    bData <- utils::read.table(bFile, sep = ',', header = T, row.names = 1)
    bData <- bData[which(!(rownames(bData) %in% c('no_match', 'total_reads'))),]
  }

  print(paste0('Initial barcodes in HTO data: ', ncol(bData)))
  bData <- DoCellFiltering(bData)

  if (doRowFilter) {
    bData <- DoRowFiltering(bData)

    # repeat colsum filter.  because we potentially dropped rows, some cells might now drop out
    bData <- DoCellFiltering(bData)
  } else {
    print('Row filtering will not be performed')

    toDrop <- rowSums(bData) == 0
    if (sum(toDrop) > 0){
      print(paste0('HTOs dropped due to zero counts in all cells: ', sum(toDrop)))
      print(paste(rownames(bData)[toDrop], collapse = ', '))
      bData <- bData[which(rowSums(bData) > 0), , drop = FALSE]
      print(paste0('HTOs after filter: ', nrow(bData)))
      print(paste(rownames(bData), collapse = ', '))
    }
  }

  #repeat summary:
  if (nrow(bData) == 0) {
    print('No HTOs remaining')
  } else {
    rowSummary <- GenerateByRowSummary(bData)
    print(kable(rowSummary, caption = 'HTO Summary After Filter', row.names = F))
  }

  return(bData)
}

#' @title DoRowFiltering
#'
#' @description A description
#' @return A modified Seurat object.
DoRowFiltering <- function(bData, minRowSum = 5,
                           minRowMax = 20,
                           minRowMean = 0.2,
                           minMeanNonZeroCount = 2){

  #thresholdRowSum <- InferThresholds(rowSums(bData), dataLabel = 'Row Sums')
  #print(thresholdRowSum)
  #TODO: consider using this
  #minRowSum <- thresholdRowSum$ElbowThreshold

  #thresholdRowMax <- InferThresholds(rowSummary$max, dataLabel = 'Row Max')
  #print(thresholdRowMax)
  #TODO: consider using this
  #minRowMax <- thresholdRowMax$ElbowThreshold

  thresholdRowMean <- InferThresholds(rowMeans(bData), dataLabel = 'Row Mean')
  print(thresholdRowMean)

  #rowsum
  toDrop <- sum(rowSums(bData) < minRowSum)
  if (toDrop > 0){
    print(paste0('HTOs dropped due to low total counts (', minRowSum, '): ', toDrop))
    print(paste(rownames(bData)[which(rowSums(bData) < minRowSum)], collapse = ', '))
    bData <- bData[which(rowSums(bData) >= minRowSum), , drop = FALSE]
    print(paste0('HTOs after filter: ', nrow(bData)))
    print(paste(rownames(bData), collapse = ', '))
  }

  #summarise
  rowSummary <- GenerateByRowSummary(bData)
  print(kable(rowSummary, caption = 'HTO Summary', row.names = F))

  #rowmean
  toDrop <- sum(rowMeans(bData) < minRowMean)
  if (toDrop > 0){
    print(paste0('HTOs dropped due to low row means (', minRowMean, '): ', toDrop))
    print(paste(rownames(bData)[which(rowMeans(bData) < minRowMean)], collapse = ', '))
    bData <- bData[which(rowMeans(bData) >= minRowMean), , drop = FALSE]
    print(paste0('HTOs after filter: ', nrow(bData)))
    print(paste(rownames(bData), collapse = ', '))
  }

  #summarise
  rowSummary <- GenerateByRowSummary(bData)
  print(kable(rowSummary, caption = 'HTO Summary', row.names = F))

  #Drop HTOs with zero strong cells:
  barcodeMatrix <- as.matrix(bData)
  rowMaxes <- apply(barcodeMatrix, 1, max)
  toDrop <- rowMaxes < minRowMax
  if (sum(toDrop) > 0){
    print(paste0('HTOs dropped due to low max counts (', minRowMax, '): ', sum(toDrop)))
    print(paste(rownames(bData)[toDrop], collapse = ', '))
    bData <- bData[!toDrop, ,  drop = FALSE ]
    print(paste0('HTOs after filter: ', nrow(bData)))
    print(paste(rownames(bData), collapse = ', '))
  }

  #Now filter HTO by mean count among non-zero cells:
  barcodeMatrix <- as.matrix(bData)
  meanNonZero <- (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix))
  meanNonZeroRatio <- meanNonZero / rowMeans(barcodeMatrix)
  print(kable(data.frame(HTO = rownames(barcodeMatrix), MeanCountOfNonZeroCells = meanNonZero, RatioOfMeanToNonZeroMean = meanNonZeroRatio), row.names = F))

  toDrop <- meanNonZero < minMeanNonZeroCount
  if (sum(toDrop) > 0){
    print(paste0('HTOs dropped due to insufficient mean non-zero count (', minMeanNonZeroCount, '): ', sum(toDrop)))
    print(paste(rownames(bData)[toDrop], collapse = ', '))
    bData <- bData[!toDrop, , drop = FALSE ]
    print(paste0('HTOs after filter: ', nrow(bData)))
    print(paste(rownames(bData), collapse = ', '))
  }

  return(bData)
}


#' @importFrom naturalsort naturalfactor
GenerateByRowSummary <- function(barcodeData) {
  if (nrow(barcodeData) == 0) {
    print('No rows in data, cannot generate summary')
    return(data.frame(HTO = character()))
  }

  barcodeMatrix <- as.matrix(barcodeData)
  df <- data.frame(HTO = naturalfactor(rownames(barcodeData)), min = apply(barcodeMatrix, 1, min), max = apply(barcodeMatrix, 1, max), mean = apply(barcodeMatrix, 1, mean), logmean = log(apply(barcodeMatrix, 1, mean) + 1), nonzero = apply(barcodeMatrix, 1, function(x){
    sum(x > 0)
  }), mean_nonzero = (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix)), total_gt1 = apply(barcodeMatrix, 1, function(x){
    sum(x > 1)
  }), mean_gt1 = apply(barcodeMatrix, 1, function(x){
    mean(sapply(x, function(y){if (y > 1) y else NA}), na.rm = T)
  }))

  df <- df[order(df$HTO),]

  return(df)
}

#' @title A Title
#'
#' @description A description
#' @return A modified Seurat object.
InferThresholds <- function(data, dataLabel, minQuant = 0.05, plotFigs = T, FindElbowMinY = NA) {
  print(paste0('Inferring thresholds for: ', dataLabel))
  ret <- list()
  if (length(data) == 0) {
    print('Unable to infer thresholds, data was empty')
    return(ret)
  }

  #this method requires one to set the quantile of outliers desired.
  ld <- log2(data + 1)

  if (plotFigs){

    boxplot(ld, main = paste0(dataLabel, " Threshold Based on Quantile"), ylab = paste0(dataLabel, " (log2)"))
    abline(h=quantile(ld, c(minQuant)), col="red")
    abline(h=quantile(ld, c(1-minQuant)), col="red")
  }

  ret$QuantileThreshold <- unname(exp(quantile(ld, c(minQuant))))
  print(paste0('Threshold inferred by quantile: ', ret$QuantileThreshold))

  #Find elbow:
  if (length(data) > 30){
    tempDF.plot <- as.data.frame(cbind(x=1:30,
                                       y = unlist(lapply(1:30, function(xN){
                                         sum(data < xN)
                                       })))); NoOfSteps = 30
  } else {
    #if it is less than 3 columns wide, then just choose 1/2 so not to break the pipeline
    tempDF.plot <- as.data.frame(cbind(x=1:round(length(data)/2),
                                       y = unlist(lapply(1:round(length(data)/2), function(xN){
                                         sum(data < xN)
                                       })))); NoOfSteps = round(length(data)/2)
  }

  if (nrow(tempDF.plot) > 1) {
    tempElbow <- FindElbow(y=(tempDF.plot$y), plot=T, min.y = FindElbowMinY, ignore.concavity = T)

    #since the FindElbow is ordering y decendingly, and we did 1:30
    tempElbow <- NoOfSteps - tempElbow
    if (plotFigs){
      plot(tempDF.plot, main = paste0(dataLabel, " Threshold Based on Elbow"), xlab = dataLabel)
      abline(v=tempElbow, col="red")
    }

    ret$ElbowThreshold <- tempElbow
    print(paste0('Threshold inferred by elbow: ', ret$ElbowThreshold))

    remove(tempElbow, tempDF.plot, NoOfSteps)
  } else {
    print('Too few points, unable to infer by elbow')
  }
  return(ret)
}

#' @title A Title
#'
#' @description A description
#' @return A modified Seurat object.
DoCellFiltering <- function(bData, minQuant = 0.05, maxValueForColSumFilter = 5){
  thresholdsSum <- InferThresholds(colSums(bData), dataLabel = "Column Sums", minQuant = minQuant, FindElbowMinY = 5000)
  minColSum <- thresholdsSum$ElbowThreshold
  minColSum <- min(minColSum, maxValueForColSumFilter)

  #colsum filter
  toDrop <- sum(colSums(bData) < minColSum)
  if (toDrop > 0){
    print(paste0('cells dropped due to low total counts per column (', minColSum, '): ', toDrop))
    bData <- bData[,which(colSums(bData) >= minColSum), drop = FALSE]
    print(paste0('Final cell barcodes: ', ncol(bData)))
  }

  #colmax filter
  #barcodeMatrix <- as.matrix(bData)
  #cm <- apply(barcodeMatrix, 2, max)
  #thresholdsMax <- InferThresholds(cm, dataLabel = "Column Max", minQuant = minQuant, FindElbowMinY = 5000)
  #minColMax <- thresholdsMax$ElbowThreshold

  #toDrop <- sum(cm < minColMax)
  #if (toDrop > 0){
  #  print(paste0('cells dropped due to low max counts per column (', minColMax,'): ', toDrop))
  #  bData <- bData[,(cm >= minColMax), drop = FALSE]
  #  print(paste0('Final cell barcodes: ', ncol(bData)))
  #}

  return(bData)
}

utils::globalVariables(
  names = c('Barcode1'),
  package = 'OOSAP',
  add = TRUE
)

#' @title GenerateQcPlots
#'
#' @description Generate QC plots for HTO/barcode data
#' @return A modified Seurat object.
#' @importFrom knitr kable
#' @export
#' @import ggplot2
GenerateQcPlots <- function(barcodeData){
  print('Generating QC Plots')

  #Plot counts/cell:
  countsPerCell <- Matrix::colSums(barcodeData)
  countsPerCell <- sort(countsPerCell)
  countAbove <-unlist(lapply(countsPerCell, function(x){
    sum(countsPerCell >= x)
  }))
  plot(log10(countAbove), log10(countsPerCell), pch=20, ylab = "log10(Reads/Cell)", xlab = "log10(Total Cells)")

  topBarcodes <- sort(tail(countsPerCell, n = 20), decreasing = T)

  print(kable(data.frame(CellBarcode = names(topBarcodes), Count = topBarcodes), row.names = F))

  #boxplot per HTO:
  barcodeMatrix <- as.matrix(barcodeData)
  melted <- setNames(reshape2::melt(barcodeMatrix), c('HTO', 'CellBarcode', 'Count'))
  print(ggplot(melted, aes(x = HTO, y = Count)) +
          geom_boxplot() +
          xlab('HTO') +
          ylab('Count') +
          ggtitle('Counts By HTO') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  melted$Count <- melted$Count + 0.5
  print(ggplot(melted, aes(x = HTO, y = Count)) +
          geom_boxplot() +
          xlab('HTO') +
          scale_y_continuous(trans='log10') +
          ylab('Count') +
          ggtitle('Counts By HTO (log)') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  #normalize columns, print top barcode fraction:
  normalizedBarcodes <- sweep(barcodeMatrix,2,colSums(barcodeMatrix),`/`)
  topValue <- apply(normalizedBarcodes,2,function(x){
    max(x)
  })

  df <- data.frame(Barcode1 = topValue)
  print(ggplot(df, aes(x = Barcode1)) +
          geom_histogram(binwidth = 0.05) +
          xlab('Top Barcode Fraction') +
          ylab('Count')
  )

  print(paste0('Total cells where top barcode is >0.75 of counts: ', length(topValue > 0.75)))

}

utils::globalVariables(
  names = c('p_val_adj', 'avg_logFC', 'cluster'),
  package = 'OOSAP',
  add = TRUE
)


utils::globalVariables(
  names = c('HTO_Classification', 'HTO_classification.global', 'HTO', 'Count'),
  package = 'OOSAP',
  add = TRUE
)

#' @title GenerateCellHashCallsSeurat
#'
#' @description Generates final cell hashing calls using Seurat3 HTODemux
#' @return A data table of results
#' @import data.table
GenerateCellHashCallsSeurat <- function(barcodeData, positive.quantile = 0.99, attemptRecovery = F, minCellsForRecovery = 500) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')

  tryCatch({
    seuratObj <- DoHtoDemux(seuratObj, positive.quantile = positive.quantile)
    dt <- data.table(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$hash.ID, HTO_classification.all = seuratObj$HTO_classification, HTO_classification.global = seuratObj$HTO_classification.global, key = c('Barcode'), stringsAsFactors = F)

    #attempt recovery:
    if (attemptRecovery && length(seuratObj$HTO_classification.global == 'Negative') > minCellsForRecovery) {
      print('Attempting recovery of negatives')
      seuratObj2 <- subset(seuratObj, subset = HTO_classification.global == 'Negative')
      seuratObj2 <- CreateSeuratObject(seuratObj2@assays$HTO@counts, assay = 'HTO')
      print(paste0('Initial negative cells : ', ncol(seuratObj2)))

      tryCatch({
        seuratObj2 <- DoHtoDemux(seuratObj2, positive.quantile = positive.quantile, label = 'Seurat HTODemux (2nd Round)')
        seuratObj2 <- subset(seuratObj2, subset = HTO_classification.global == 'Singlet')
        print(paste0('Total cells rescued by second HTODemux call: ', ncol(seuratObj2)))
        if (ncol(seuratObj2) > 0) {
          dt2 <- data.table(Barcode = as.factor(colnames(seuratObj2)), HTO_classification = seuratObj2$hash.ID, HTO_classification.all = seuratObj2$HTO_classification, HTO_classification.global = seuratObj2$HTO_classification.global, key = c('Barcode'), stringsAsFactors = F)
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification <- dt2$HTO_classification
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification.all <- dt2$HTO_classification.all
          dt[dt$Barcode %in% dt2$Barcode]$HTO_classification.global <- dt2$HTO_classification.global
        }
      }, error = function(e){
        print('Error generating second round of seurat calls, using first round')
        return(dt)
      })
    }

    return(dt)
  }, error = function(e){
    print(e)
    print('Error generating seurat calls, aborting')
    return(NA)
  })
}

utils::globalVariables(
  names = c('sortOrder'),
  package = 'OOSAP',
  add = TRUE
)

#' @title AppendCellHashing
#'
#' @description Appends cell hashing calls to a seurat object
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
#' @importFrom dplyr arrange
AppendCellHashing <- function(seuratObj, barcodeCallFile, barcodePrefix = NULL) {
  initialCells <- ncol(seuratObj)
  print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))

  if (!file.exists(barcodeCallFile)) stop("Barcode File Not found")

  barcodeCallTable <- utils::read.table(barcodeCallFile, sep = '\t', header = T)
  if (!is.null(barcodePrefix)) {
    barcodeCallTable$CellBarcode <- paste0(barcodePrefix, '_', barcodeCallTable$CellBarcode)

    initialCells <- sum(seuratObj$BarcodePrefix == barcodePrefix)
    print(paste0('Initial cell barcodes in GEX data for prefix: ', initialCells))
  }

  #Hack until we figure this out upstream
  #TODO: Find discordant duplicates add as second col or convert to dublets or some 99 err
  barcodeCallTable <- unique(barcodeCallTable)


  barcodeCallTable <- barcodeCallTable[barcodeCallTable$HTO != 'Negative',]
  if (nrow(barcodeCallTable)==0) stop("Something is wrong, table became 0 rows")

  print(paste0('Non-negative cell barcodes in HTO calls: ', nrow(barcodeCallTable)))

  #dup <- barcodeCallTable[barcodeCallTable$CellBarcode %in% barcodeCallTable$CellBarcode[duplicated(barcodeCallTable$CellBarcode)],]
  #write.table(dup, file='duplicates.txt', sep = '\t', quote = F, row.names = F)

  # Dont overwrite in case we already added data for another dataset
  if (!('HTO' %in% names(seuratObj@meta.data))) {
    print('Adding HTO columns to seurat object')
    seuratObj$HTO <- c(NA)
    seuratObj$HTO_Classification <- c(NA)
  }

  datasetSelect <- seuratObj$BarcodePrefix == barcodePrefix
  df <- data.frame(CellBarcode = colnames(seuratObj)[datasetSelect])
  df$sortOrder = 1:nrow(df)
  df <- merge(df, barcodeCallTable, all.x = T, all.y = F, by = c('CellBarcode'))
  df <- dplyr::arrange(df, sortOrder)
  df <- df[names(df) != 'sortOrder']

  if (sum(datasetSelect) != nrow(df)) {
    stop('Length of data select and df do not match!')
  }

  df$HTO <- as.character(df$HTO)
  df$HTO_Classification <- as.character(df$HTO_Classification)

  df$HTO[is.na(df$HTO)] <- 'ND'
  df$HTO_Classification[is.na(df$HTO_Classification)] <- 'ND'
  print(paste0('Total HTOs added: ', nrow(df[!is.na(df$HTO),])))

  # Check barcodes match before merge
  if (sum(df$cellBarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix]) > 0) {
    stop(paste0('Seurat and HTO barcodes do not match after merge, differences: ', sum(df$cellBarcode != colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix])))
  }

  HTO <- as.character(seuratObj$HTO)
  HTO_Classification <- as.character(seuratObj$HTO_Classification)

  HTO[datasetSelect] <- df$HTO
  HTO_Classification[datasetSelect] <- df$HTO_Classification

  seuratObj$HTO <- as.factor(HTO)
  seuratObj$HTO_Classification <- as.factor(HTO_Classification)

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @param seurObj, A Seurat object.
#' @importFrom cluster clara
#' @importFrom Matrix t
#' @return A modified Seurat object.
DebugDemux <- function(seuratObj, assay = 'HTO', reportKmeans = FALSE) {
  print('Debugging information for Seurat HTODemux:')
  data <- GetAssayData(object = seuratObj, assay = assay)
  ncenters <- (nrow(x = data) + 1)

  init.clusters <- clara(
    x = t(x = data),
    k = ncenters,
    samples = 100
  )
  Idents(object = seuratObj, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering

  # Calculate tSNE embeddings with a distance matrix
  perplexity <- .InferPerplexity(seuratObj, 100)
  seuratObj[['hto_tsne']] <- RunTSNE(dist(t(data)), assay = assay, perplexity = perplexity)
  P <- DimPlot(seuratObj, reduction = 'hto_tsne', label = TRUE)
  P <- P + ggtitle('Clusters: clara')
  print(P)

  average.expression <- AverageExpression(
    object = seuratObj,
    assays = c(assay),
    verbose = FALSE
  )[[assay]]

  print(knitr::kable(average.expression, label = 'clara'))

  if (reportKmeans) {
    print('kmeans:')
    init.clusters <- kmeans(
      x = t(x = data),
      centers = ncenters,
      nstart = 100
    )
    Idents(object = seuratObj, cells = names(x = init.clusters$cluster), drop = TRUE) <- init.clusters$cluster

    # Calculate tSNE embeddings with a distance matrix
    P <- DimPlot(seuratObj, label = TRUE)
    P <- P + ggtitle('Clusters: kmeans')
    print(P)

    average.expression <- AverageExpression(
      object = seuratObj,
      assays = c(assay),
      verbose = FALSE
    )[[assay]]

    print(knitr::kable(average.expression, label = 'kmeans'))
  }
}

#' @title A Title
#'
#' @description A description
#' @param seurObj, A Seurat object.
#' @return A modified Seurat object.
DoHtoDemux <- function(seuratObj, positive.quantile = 0.99, label = 'Seurat HTODemux', plotDist = FALSE) {
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", verbose = FALSE)

  DebugDemux(seuratObj)

  seuratObj <- HTODemux2(seuratObj, positive.quantile =  positive.quantile, plotDist = plotDist)

  HtoSummary(seuratObj, label = label, htoClassificationField = 'hash.ID', globalClassificationField = 'HTO_classification.global')

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @return A modified Seurat object.
#' @import data.table
GenerateCellHashCallsMultiSeq <- function(barcodeData) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')

  tryCatch({
    seuratObj <- DoMULTIseqDemux(seuratObj)

    return(data.table(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$MULTI_ID, HTO_classification.global = seuratObj$MULTI_classification.global, key = c('Barcode')))
  }, error = function(e){
    print(e)
    print('Error generating multiseq calls, aborting')
    return(NA)
  })
}

#' @title GenerateCellHashingCalls
#'
#' @description A description
#' @return A data table of results.
#' @export
#' @import data.table
GenerateCellHashingCalls <- function(barcodeData, positive.quantile = 0.99, attemptRecovery = FALSE, useSeurat = TRUE, useMultiSeq = TRUE, outFile = 'combinedHtoCalls.txt', allCallsOutFile = NA) {
  sc <- NA
  if (useSeurat) {
    sc <- GenerateCellHashCallsSeurat(barcodeData, positive.quantile = positive.quantile, attemptRecovery = attemptRecovery )
  }

  mc <- NA
  if (useMultiSeq) {
    mc <- GenerateCellHashCallsMultiSeq(barcodeData)
  }

  dt <- ProcessEnsemblHtoCalls(mc, sc, barcodeData, outFile = outFile, allCallsOutFile = allCallsOutFile)

  return(dt)
}

#' @title A Title
#'
#' @description A description
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
DoMULTIseqDemux <- function(seuratObj, assay = 'HTO', autoThresh = TRUE, quantile = NULL, maxiter = 20, qrange = seq(from = 0.2, to = 0.95, by = 0.05)) {

  ## Normalize Data: Log2 Transform, mean-center
  counts <- GetAssayData(seuratObj, assay = assay, slot = 'counts')
  log2Scaled <- as.data.frame(log2(counts))
  for (i in 1:ncol(counts)) {
    ind <- which(is.finite(log2Scaled[,i]) == FALSE)
    log2Scaled[ind,i] <- 0
    log2Scaled[,i] <- log2Scaled[,i]-mean(log2Scaled[,i])
  }
  seuratObjMS <- CreateSeuratObject(counts, assay = 'MultiSeq')
  seuratObjMS[['MultiSeq']]@data <- as.matrix(log2Scaled)

  seuratObjMS <- MULTIseqDemux(seuratObjMS, assay = "MultiSeq", quantile = quantile, verbose = TRUE, autoThresh = autoThresh, maxiter = maxiter, qrange = qrange)

  seuratObj$MULTI_ID <- as.character(seuratObjMS$MULTI_ID)
  seuratObj$MULTI_classification.global <- as.character(seuratObjMS$MULTI_ID)
  seuratObj$MULTI_classification.global[!(seuratObjMS$MULTI_ID %in% c('Negative', 'Doublet'))] <- 'Singlet'
  seuratObj$MULTI_classification.global <- as.factor(seuratObj$MULTI_classification.global)

  HtoSummary(seuratObj, label = 'MULTI-SEQ', htoClassificationField = 'MULTI_ID', globalClassificationField = 'MULTI_classification.global')

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @param seurObj, A Seurat object.
#' @return A modified Seurat object.
HtoSummary <- function(seuratObj, htoClassificationField, globalClassificationField, label, doHeatmap = T, doTSNE = T, assay = 'HTO') {
  #report outcome
  print(table(seuratObj[[htoClassificationField]]))
  print(table(seuratObj[[globalClassificationField]]))

  # Group cells based on the max HTO signal
  seuratObj_hashtag <- seuratObj
  Idents(seuratObj_hashtag) <- globalClassificationField
  htos <- rownames(GetAssayData(seuratObj_hashtag, assay = assay))
  for (hto in naturalsort::naturalsort(htos)){
    print(VlnPlot(seuratObj_hashtag, features = c(hto), assay = assay, ncol = 1, log = T) + ggtitle(paste0(label, ": ", hto, " (log)")))
    print(VlnPlot(seuratObj_hashtag, features = c(hto), assay = assay, ncol = 1, log = F) + ggtitle(paste0(label, ": ", hto)))
  }

  if (doTSNE) {
    perplexity <- .InferPerplexity(seuratObj, 100)
    seuratObj[['hto_tsne']] <- RunTSNE(dist(Matrix::t(GetAssayData(seuratObj, slot = "data", assay = assay))), assay = assay, perplexity = perplexity)
    print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = htoClassificationField, label = TRUE) + ggtitle(label))
    print(DimPlot(seuratObj, reduction = 'hto_tsne', group.by = globalClassificationField, label = TRUE) + ggtitle(label))
  }

  if (doHeatmap) {
    print(HTOHeatmap(seuratObj, assay = assay, classification = htoClassificationField, global.classification = globalClassificationField, ncells = min(3000, ncol(seuratObj)), singlet.names = NULL) + ggtitle(label))
  }
}

utils::globalVariables(
  names = c('HTO_classification.global.MultiSeq', 'HTO_classification.global.Seurat', 'n', 'HTO_classification.MultiSeq', 'HTO_classification.Seurat'),
  package = 'OOSAP',
  add = TRUE
)

#' @title ProcessEnsemblHtoCalls
#'
#' @description A description
#' @return A modified Seurat object.
#' @export
#' @import data.table
#' @import ggplot2
#' @param mc Multiseq calls dataframe
#' @param sc Seurat calls dataframe
#' @param barcodeData The barcode count matrix
#' @param outFile The output TSV file
#' @param allCallsOutFile If provided, a more detailed output will be written here
#' @importFrom dplyr %>% group_by summarise
ProcessEnsemblHtoCalls <- function(mc, sc, barcodeData,
                                   outFile = 'combinedHtoCalls.txt',
                                   allCallsOutFile = NA) {

  if (all(is.na(sc)) && all(is.na(mc))){
    print('MULTI-seq and Seurat failed to produce calls, aborting')
    return()
  }

  if (all(is.na(sc))){
    print('No calls for seurat found')
    dt <- data.table(CellBarcode = mc$Barcode, HTO = mc$HTO_classification, HTO_Classification = mc$HTO_classification.global, key = 'CellBarcode', Seurat = c(F), MultiSeq = c(T))
    dt <- PrintFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  if (all(is.na(mc))){
    print('No calls for MULTI-seq found')
    dt <- data.table(CellBarcode = sc$Barcode, HTO = sc$HTO_classification, HTO_Classification = sc$HTO_classification.global, key = 'CellBarcode', Seurat = c(T), MultiSeq = c(F))
    dt <- PrintFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  mc$Barcode <- as.character(mc$Barcode)
  sc$Barcode <- as.character(sc$Barcode)
  merged <- merge(mc, sc, all = T, by = 'Barcode', suffixes = c('.MultiSeq', '.Seurat'))

  merged$HTO_classification.MultiSeq <- as.character(merged$HTO_classification.MultiSeq)
  merged$HTO_classification.MultiSeq[is.na(merged$HTO_classification.MultiSeq)] <- 'Negative'
  merged$HTO_classification.MultiSeq <- naturalsort::naturalfactor(merged$HTO_classification.MultiSeq)

  merged$HTO_classification.Seurat <- as.character(merged$HTO_classification.Seurat)
  merged$HTO_classification.Seurat[is.na(merged$HTO_classification.Seurat)] <- 'Negative'
  merged$HTO_classification.Seurat <- naturalsort::naturalfactor(merged$HTO_classification.Seurat)

  merged$HTO_classification.global.MultiSeq <- as.character(merged$HTO_classification.global.MultiSeq)
  merged$HTO_classification.global.MultiSeq[is.na(merged$HTO_classification.global.MultiSeq)] <- 'Negative'
  merged$HTO_classification.global.MultiSeq <- naturalsort::naturalfactor(merged$HTO_classification.global.MultiSeq)

  merged$HTO_classification.global.Seurat <- as.character(merged$HTO_classification.global.Seurat)
  merged$HTO_classification.global.Seurat[is.na(merged$HTO_classification.global.Seurat)] <- 'Negative'
  merged$HTO_classification.global.Seurat <- naturalsort::naturalfactor(merged$HTO_classification.global.Seurat)

  merged$HasSeuratCall <- merged$HTO_classification.Seurat != 'Negative'
  merged$HasMultiSeqCall <- merged$HTO_classification.MultiSeq != 'Negative'

  #dont count situations where one side is negative as discordant
  merged$Concordant <- as.character(merged$HTO_classification.MultiSeq) == as.character(merged$HTO_classification.Seurat)
  merged$Concordant[!merged$Concordant & (merged$HTO_classification.MultiSeq == 'Negative' | merged$HTO_classification.Seurat == 'Negative')] <- TRUE

  merged$GlobalConcordant <- as.character(merged$HTO_classification.global.MultiSeq) == as.character(merged$HTO_classification.global.Seurat)
  merged$GlobalConcordant[!merged$GlobalConcordant & (merged$HTO_classification.global.MultiSeq == 'Negative' | merged$HTO_classification.global.Seurat == 'Negative')] <- TRUE

  print(paste0('Total concordant: ', nrow(merged[merged$Concordant])))
  print(paste0('Total discordant (HTO call): ', nrow(merged[!merged$Concordant])))
  print(paste0('Total discordant (HTO classification): ', nrow(merged[!merged$GlobalConcordant])))

  tbl <- data.frame(table(Concordant = merged$Concordant))

  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(names(tbl)))), "Set1"))
  colorValues <- getPalette(length(names(tbl)))

  print(ggplot(tbl, aes(x="", y=Freq, fill=Concordant)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
		coord_polar("y", start=0) +
    scale_fill_manual(values = colorValues) +
		theme_minimal() +
		theme(
			axis.text.x=element_blank(),
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			panel.grid  = element_blank()
		) +
    ggtitle('Seurat/MultiSeq Concordance')
	)

  discord <- merged[!merged$GlobalConcordant]
  discord <- discord %>% group_by(HTO_classification.global.MultiSeq, HTO_classification.global.Seurat) %>% summarise(Count = dplyr::n())

  print(qplot(x=HTO_classification.global.MultiSeq, y=HTO_classification.global.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By Global Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  discord <- merged[!merged$Concordant]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = dplyr::n())

  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By HTO Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  # These calls should be identical, except for possibly negatives from one method that are non-negative in the other
  # For the time being, accept those as correct.
  merged$FinalCall <- as.character(merged$HTO_classification.MultiSeq)
  merged$FinalCall[merged$HTO_classification.MultiSeq == 'Negative'] <- as.character(merged$HTO_classification.Seurat[merged$HTO_classification.MultiSeq == 'Negative'])
  merged$FinalCall[!merged$Concordant] <- 'Discordant'
  merged$FinalCall <- naturalsort::naturalfactor(merged$FinalCall)

  merged$FinalClassification <- as.character(merged$HTO_classification.global.MultiSeq)
  merged$FinalClassification[merged$HTO_classification.global.MultiSeq == 'Negative'] <- as.character(merged$HTO_classification.global.Seurat[merged$HTO_classification.global.MultiSeq == 'Negative'])
  merged$FinalClassification[!merged$GlobalConcordant] <- 'Discordant'
  merged$FinalClassification[!merged$Concordant] <- 'Discordant'
  merged$FinalClassification <- as.factor(merged$FinalClassification)

  df <- data.frame(
    TotalSinglet = c(sum(merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet')),
    ConcordantSinglet = c(sum(merged$Concordant & merged$HTO_classification.global.Seurat == 'Singlet'), sum(merged$Concordant & merged$HTO_classification.global.MultiSeq == 'Singlet'), sum(merged$FinalClassification == 'Singlet'))
  )
  rownames(df) <- c('Seurat', 'MultiSeq', 'Final')
  df <- t(df)
  print(knitr::kable(df))

  if (!is.na(allCallsOutFile) && nrow(merged) > 0) {
    write.table(merged, file = allCallsOutFile, row.names = F, sep = '\t', quote = F)
  }

  if (nrow(merged) > 0){
    dt <- data.frame(CellBarcode = merged$Barcode, HTO = merged$FinalCall, HTO_Classification = merged$FinalClassification, key = 'CellBarcode', Seurat = merged$HasSeuratCall, MultiSeq = merged$HasMultiSeqCall)
    dt <- PrintFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)

  } else {
    print('No rows, not saving ')
  }
}


utils::globalVariables(
  names = c('HTO_Classification', 'TotalCounts'),
  package = 'OOSAP',
  add = TRUE
)

#' @title PrintFinalSummary
#'
#' @return A modified Seurat object.
#' @importFrom naturalsort naturalfactor
#' @importFrom knitr kable
#' @importFrom data.table melt
#' @param The data table with calls
#' @param The barcode counts matrix
#' @import ggplot2
PrintFinalSummary <- function(dt, barcodeData){
  #Append raw counts:
  bc <- t(barcodeData)
  x <- reshape2::melt(bc)
  names(x) <- c('CellBarcode', 'HTO', 'Count')

  merged <- merge(dt, x, by = c('CellBarcode', 'HTO'), all.x = T, all.y = F)

  bc <- as.data.frame(bc)
  bc$CellBarcode <- rownames(bc)
  merged <- merge(merged, bc, by = c('CellBarcode'), all.x = T, all.y = F)

  merged$HTO <- as.character(merged$HTO)
  merged$HTO[is.na(merged$HTO)] <- c('Negative')
  merged$HTO <- as.factor(merged$HTO)

  merged$HTO_Classification <- as.character(merged$HTO_Classification)
  merged$HTO_Classification[is.na(merged$HTO_Classification)] <- c('Negative')
  merged$HTO_Classification <- as.factor(merged$HTO_Classification)

  #summarise reads by type:
  barcodeMatrix <- as.matrix(barcodeData)
  cs <- colSums(barcodeMatrix)
  cs <- cs[merged$CellBarcode]
  merged$TotalCounts <- cs

  htoNames <- simplifyHtoNames(as.character(merged$HTO))

  merged$HTO <- naturalfactor(as.character(htoNames))

  t <- table(SeuratCall = merged$Seurat, MultiSeqCall = merged$MultiSeq)

  colnames(t)[colnames(t) == T] <- c('MultiSeq Call')
  colnames(t)[colnames(t) == F] <- c('MultiSeq No Call')

  rownames(t)[rownames(t) == T] <- c('Seurat Call')
  rownames(t)[rownames(t) == F] <- c('Seurat No Call')

  print(kable(t))

  print(ggplot(merged, aes(x = HTO)) +
          geom_bar(stat = 'count') +
          xlab('HTO') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  tbl <- table(HTO = merged$HTO)
  df <- data.frame(tbl)
  df$Pct <- round((df$Freq / sum(df$Freq)) * 100, 2)
  print(kable(df))

  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(names(tbl)))), "Set1"))
  colorValues <- getPalette(length(names(tbl)))

  print(ggplot(df, aes(x = '', y=Freq, fill=HTO)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = colorValues) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    ) +
    ggtitle('HTO')
  )

  print(ggplot(merged, aes(x = HTO)) +
          geom_bar(stat = 'count') +
          xlab('HTO') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  tbl <- table(HTO_Classification = merged$HTO_Classification)
  df <- data.frame(tbl)
  df$Pct <- round((df$Freq / sum(df$Freq)) * 100, 2)
  print(kable(df))

  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(names(tbl)))), "Set1"))
  colorValues <- getPalette(length(names(tbl)))

  print(ggplot(df, aes(x = '', y=Freq, fill=HTO_Classification)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = colorValues) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank()
    ) +
    ggtitle('HTO Classification')
  )

  print(ggplot(merged, aes(x = HTO_Classification, y = TotalCounts)) +
          geom_boxplot()  +
          xlab('HTO Classification') +
          ylab('Counts Per Cell (log10)') +
          ggtitle('Counts By Call Type') +
          scale_y_continuous(trans='log10') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  return(merged)
}



simplifyHtoNames <- function(v) {
  return(sapply(v, function(x){
    x <- unlist(strsplit(x, '-'))
    if (length(x) > 1) {
      x <- x[-(length(x))]
    }

    paste0(x, collapse = "-")
  }))
}



#' @title GenerateSummaryForExpectedBarcodes
#'
#' @description A description
#' @export
GenerateSummaryForExpectedBarcodes <- function(dt, whitelistFile, outputFile, barcodeData) {
  categoryName <- "Cell Hashing Concordance"

  whitelist <- utils::read.table(whitelistFile, sep = '\t', header = F)
  names(whitelist) <- c('CellBarcode')
  df <- data.frame(Category = categoryName, MetricName = "InputBarcodes", Value = length(whitelist$CellBarcode))


  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCounts", Value = sum(barcodeData)))

  #Any called cell:
  calledCellBarcodes <- dt$CellBarcode[dt$HTO_Classification != 'Negative' & dt$HTO_Classification != 'Discordant']
  calledIntersect <- intersect(whitelist$CellBarcode, calledCellBarcodes)

  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalled", Value = length(calledCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalledOverlapping", Value = length(calledIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputCalled", Value = (length(calledIntersect) / length(whitelist$CellBarcode))))

  totalCalledNotInInput <- length(calledCellBarcodes) - length(calledIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCalledNotInInput", Value = totalCalledNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionCalledNotInInput", Value = (totalCalledNotInInput / length(calledCellBarcodes))))

  #Singlets:
  singletCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Singlet']
  singletIntersect <- intersect(whitelist$CellBarcode, singletCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSinglet", Value = length(singletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSingletOverlapping", Value = length(singletIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputSinglet", Value = (length(singletIntersect) / length(whitelist$CellBarcode))))

  totalSingletNotInInput <- length(singletCellBarcodes) - length(singletIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalSingletNotInInput", Value = totalSingletNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionSingletNotInInput", Value = (totalSingletNotInInput / length(singletCellBarcodes))))

  #Doublets:
  doubletCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Doublet']
  doubletIntersect <- intersect(whitelist$CellBarcode, doubletCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoublet", Value = length(doubletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoubletOverlapping", Value = length(doubletIntersect)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputDoublet", Value = (length(doubletIntersect) / length(whitelist$CellBarcode))))

  totalDoubletNotInInput <- length(doubletCellBarcodes) - length(doubletIntersect)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDoubletNotInInput", Value = totalDoubletNotInInput))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionDoubletNotInInput", Value = (totalDoubletNotInInput / length(doubletCellBarcodes))))

  #Discordant:
  discordantCellBarcodes <- dt$CellBarcode[dt$HTO_Classification == 'Discordant']
  discordantIntersect <- intersect(whitelist$CellBarcode, discordantCellBarcodes)
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalDiscordant", Value = length(doubletCellBarcodes)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "FractionOfInputDiscordant", Value = (length(discordantIntersect) / length(whitelist$CellBarcode))))


  #By caller:
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "SeuratNonNegative", Value = sum(dt$HasSeuratCall & dt$HTO_Classification != 'Discordant')))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "MultiSeqNonNegative", Value = sum(dt$HasMultiSeqCall & dt$HTO_Classification != 'Discordant')))

  df$Value[is.na(df$Value)] <- 0

  write.table(df, file = outputFile, quote = F, row.names = F, sep = '\t')
}


#' @title DownloadAndAppendCellHashing
#'
#' @description A description
#' @param seurObject, A Seurat object.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendCellHashing <- function(seuratObject, outPath = '.'){
  if (is.null(seuratObject[['BarcodePrefix']])){
    stop('Seurat object lacks BarcodePrefix column')
  }

  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    print(paste0('Possibly adding cell hashing data for prefix: ', barcodePrefix))

    cellHashingId <- FindMatchedCellHashing(barcodePrefix)
    if (is.null(cellHashingId)){
      print(paste0('Cell hashing not used for prefix: ', barcodePrefix, ', skipping'))
      next(seuratObject)
    } else if (is.na(cellHashingId)){
      stop(paste0('Unable to find cellHashing calls table file for prefix: ', barcodePrefix))
    }

    callsFile <- file.path(outPath, paste0(barcodePrefix, '_cellHashingCalls.csv'))
    DownloadOutputFile(outputFileId = cellHashingId, outFile = callsFile, overwrite = T)
    if (!file.exists(callsFile)){
      stop(paste0('Unable to download calls table for prefix: ', barcodePrefix))
    }

    seuratObject <- AppendCellHashing(seuratObj = seuratObject, barcodeCallFile = callsFile, barcodePrefix = barcodePrefix)
  }

  return(seuratObject)
}


#' @title FindMatchedCellHashing
#'
#' @description A description
#' @return A modified Seurat object.
#' @import Rlabkey
FindMatchedCellHashing <- function(loupeDataId){
  #Note: the seurat object gets associated with the GEX readset, so look based on this:
  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="readset,library_id",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) == 0) {
    print(paste0("Loupe File ID: ", loupeDataId, " not found"))
    return(NA)
  }

  readset <- unique(rows[['readset']])

  if (is.na(readset) || is.null(readset)) {
    print("readset is NA/NULL")
    return(NA)
  }

  libraryId <- unique(rows[['library_id']])

  #determine whether we expect cell hashing to be used:
  cDNAs <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="tcrdb",
    queryName="cdnas",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("readsetId", "EQUALS", readset)),
    colSelect="rowid,readsetid,hashingreadsetid",
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(cDNAs) == 0) {
    stop(paste0('No cDNA records found for GEX readset: ', readset))
  } else if (sum(!is.na(cDNAs$hashingreadsetid)) == 0) {
    print(paste0('The cDNA library does not use cell hashing, aborting'))
    return(NULL)
  }

  rowsB <- suppressWarnings(labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    colSort="-rowid",
    colSelect="rowid",
    colFilter=makeFilter(c("readset", "EQUAL", readset),
                         c("category", "EQUAL", "Seurat Cell Hashing Calls"),
                         c("library_id", "EQUAL", libraryId)),
    containerFilter=NULL,
    colNameOpt="rname"
  ))

  ret <- NULL
  if (nrow(rowsB) == 0){
    print(paste0("Output of type 'Seurat Cell Hashing Calls' not found.  Readset: ", readset, ", libraryId: ", libraryId))
  } else {
    ret <- rowsB[1]
  }

  if (all(is.null(ret))){
    print("Trying to find output of type: '10x GEX Cell Hashing Calls'")
    rowsB <- suppressWarnings(labkey.selectRows(
      baseUrl=lkBaseUrl,
      folderPath=lkDefaultFolder,
      schemaName="sequenceanalysis",
      queryName="outputfiles",
      colSort="-rowid",
      colSelect="rowid,",
      colFilter=makeFilter(c("readset", "EQUAL", readset),
                           c("category", "EQUAL", "10x GEX Cell Hashing Calls"),
                           c("library_ld", "EQUAL", libraryId)),
      containerFilter=NULL,
      colNameOpt="rname"
    ))

    if (nrow(rowsB) == 0){
      print("Not found")
    } else {
      ret <- rowsB[1]
      print("Found output")
    }
  }

  if (all(is.null(ret))) {
    return(NA)
  } else {
    if (length(ret) > 0){
      print('More than one matching file found, using most recent')
    }

    return(ret$rowid[1])
  }
}


#' @title DownloadOutputFile
#'
#' @import Rlabkey
DownloadOutputFile <- function(outputFileId, outFile, overwrite = T) {
  if (is.na(outputFileId)) {
    stop('Output file ID cannot be NA')
  }

  if (file.exists(outFile) & !overwrite) {
    print(paste0("File exists, will not overwrite: ", outFile))
    return(outFile)
	}

  #There should be a file named all_contig_annotations.csv in the same directory as the VLoupe file
  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
    colFilter=makeFilter(c("rowid", "EQUAL", outputFileId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    stop(paste0('More than one matching file found, this should not occur.  RowId: ', outputFileId))
  }

  wb <- rows[['workbook_workbookid']]
  if (is.na(wb) || is.null(wb)){
    wb <- ''
  }

  remotePath <- rows[['dataid_webdavurlrelative']]

  success <- labkey.webdav.get(
    baseUrl=lkBaseUrl,
    folderPath=paste0(lkDefaultFolder,wb),
    remoteFilePath = remotePath,
    overwrite = overwrite,
    localFilePath = outFile
  )

  if (!success || !file.exists(outFile)) {
    stop(paste0('labkey.webdav.get failed for file: ', remotePath))
  }

  return(outFile)
}
