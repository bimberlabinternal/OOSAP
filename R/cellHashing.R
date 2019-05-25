library(data.table)


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
processCiteSeqCount <- function(bFile=NA, doRowFilter = T) {
  if (is.na(bFile)){
    stop("No file set: change bFile")
  }

  if (dir.exists(bFile)) {
    #CITE-seq-Count 1.4.2 and higher creates a folder
    bData <- Read10X(bFile, gene.column=1)
    bData <- bData[which(!(rownames(bData) %in% c('unmapped'))), , drop = F]
    bData <- as.matrix(bData)
  } else {
    # older versions create a CSV file
    bData <- read.table(bFile, sep = ',', header = T, row.names = 1)
    bData <- bData[which(!(rownames(bData) %in% c('no_match', 'total_reads'))),]
  }

  print(paste0('Initial barcodes in HTO data: ', ncol(bData)))
  bData <- doCellFiltering(bData)

  if (doRowFilter) {
    bData <- doRowFiltering(bData)

    # repeat colsum filter.  because we potentially dropped rows, some cells might now drop out
    bData <- doCellFiltering(bData)
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
    rowSummary <- generateByRowSummary(bData)
    print(kable(rowSummary, caption = 'HTO Summary After Filter', row.names = F))
  }

  return(bData)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
doRowFiltering <- function(bData, minRowSum = 5,
                           minRowMax = 20,
                           minRowMean = 0.2,
                           minMeanNonZeroCount = 2){

  #thresholdRowSum <- inferThresholds(rowSums(bData), dataLabel = 'Row Sums')
  #print(thresholdRowSum)
  #TODO: consider using this
  #minRowSum <- thresholdRowSum$ElbowThreshold

  #thresholdRowMax <- inferThresholds(rowSummary$max, dataLabel = 'Row Max')
  #print(thresholdRowMax)
  #TODO: consider using this
  #minRowMax <- thresholdRowMax$ElbowThreshold

  thresholdRowMean <- inferThresholds(rowMeans(bData), dataLabel = 'Row Mean')
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

  #summarize
  rowSummary <- generateByRowSummary(bData)
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

  #summarize
  rowSummary <- generateByRowSummary(bData)
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

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
generateByRowSummary <- function(barcodeData) {
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
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
inferThresholds <- function(data, dataLabel, minQuant = 0.05, plotFigs = T, findElbowMinY = NA) {
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
    tempElbow <- findElbow(y=(tempDF.plot$y), plot=F, min.y = findElbowMinY, ignore.concavity = T)

    #since the findElbow is ordering y decendingly, and we did 1:30
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
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
doCellFiltering <- function(bData, minQuant = 0.05, maxValueForColSumFilter = 5){
  thresholdsSum <- inferThresholds(colSums(bData), dataLabel = "Column Sums", minQuant = minQuant, findElbowMinY = 5000)
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
  #thresholdsMax <- inferThresholds(cm, dataLabel = "Column Max", minQuant = minQuant, findElbowMinY = 5000)
  #minColMax <- thresholdsMax$ElbowThreshold

  #toDrop <- sum(cm < minColMax)
  #if (toDrop > 0){
  #  print(paste0('cells dropped due to low max counts per column (', minColMax,'): ', toDrop))
  #  bData <- bData[,(cm >= minColMax), drop = FALSE]
  #  print(paste0('Final cell barcodes: ', ncol(bData)))
  #}

  return(bData)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
generateQcPlots <- function(barcodeData){
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
  melted <- setNames(melt(barcodeMatrix), c('HTO', 'CellBarcode', 'Count'))
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

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
generateCellHashCallsSeurat <- function(barcodeData) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')

  tryCatch({
    seuratObj <- doHtoDemux(seuratObj)

    return(data.table(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$hash.ID, HTO_classification.all = seuratObj$HTO_classification, HTO_classification.global = seuratObj$HTO_classification.global, key = c('Barcode')))
  }, error = function(e){
    print(e)
    print('Error generating seurat calls, aborting')
    return(NA)
  })
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
appendCellHashing <- function(seuratObj, barcodeCallFile, barcodePrefix = NULL) {
  initialCells <- ncol(seuratObj)
  print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))

  if(!file.exists(barcodeCallFile)) stop("Barcode File Not found")

  barcodeCallTable <- read.table(barcodeCallFile, sep = '\t', header = T)
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
  df <- arrange(df, sortOrder)
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
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
debugDemux <- function(seuratObj) {
  print('Debugging information for Seurat HTODemux:')
  data <- GetAssayData(object = seuratObj, assay = 'HTO')
  ncenters <- (nrow(x = data) + 1)

  init.clusters <- clara(
    x = t(x = GetAssayData(object = seuratObj, assay = 'HTO')),
    k = ncenters,
    samples = 100
  )
  #identify positive and negative signals for all HTO
  Idents(object = seuratObj, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering

  average.expression <- AverageExpression(
    object = seuratObj,
    assay = "HTO",
    verbose = FALSE
  )[["HTO"]]

  print(average.expression)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
doHtoDemux <- function(seuratObj) {
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", display.progress = FALSE)

  debugDemux(seuratObj)

  seuratObj <- HTODemux2(seuratObj, positive.quantile =  0.99)

  htoSummary(seuratObj, field1 = 'HTO_classification.global', field2 = 'hash.ID')

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
generateCellHashCallsMultiSeq <- function(barcodeData) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')

  tryCatch({
    seuratObj <- doMULTIseqDemux(seuratObj)

    return(data.table(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$MULTI_ID, HTO_classification.global = seuratObj$MULTI_classification.global, key = c('Barcode')))
  }, error = function(e){
    print(e)
    print('Error generating multiseq calls, aborting')
    return(NA)
  })
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
doMULTIseqDemux <- function(seuratObj) {
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", display.progress = FALSE)

  seuratObj <- MULTIseqDemux(seuratObj, assay = "HTO", quantile = 0.7, autoThresh = T, maxiter = 50, qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

  seuratObj$MULTI_classification.global <- as.character(seuratObj$MULTI_ID)
  seuratObj$MULTI_classification.global[!(seuratObj$MULTI_ID %in% c('Negative', 'Doublet'))] <- 'Singlet'
  seuratObj$MULTI_classification.global <- as.factor(seuratObj$MULTI_classification.global)

  htoSummary(seuratObj, field1 = 'MULTI_classification.global', field2 = 'MULTI_ID', doHeatmap = F)

  return(seuratObj)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
htoSummary <- function(seuratObj, field1, field2, doHeatmap = T) {
  #report outcome
  print(table(seuratObj[[field1]]))
  print(table(seuratObj[[field2]]))

  # Group cells based on the max HTO signal
  seuratObj_hashtag <- seuratObj
  Idents(seuratObj_hashtag) <- field2
  htos <- rownames(GetAssayData(seuratObj_hashtag,assay = "HTO"))
  for (hto in htos){
    print(RidgePlot(seuratObj_hashtag, features = c(hto), assay = 'HTO', ncol = 1))
  }

  if (doHeatmap) {
    print(HTOHeatmap(seuratObj, assay = "HTO", classification = field2, global.classification = field1, ncells = min(3000, ncol(seuratObj)), singlet.names = NULL))
  }
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
processEnsemblHtoCalls <- function(mc, sc, barcodeData,
                                   outFile = 'combinedHtoCalls.txt',
                                   allCallsOutFile = NA) {

  if (all(is.na(sc)) && all(is.na(mc))){
    print('MULTI-seq and Seurat failed to produce calls, aborting')
    return()
  }

  if (all(is.na(sc))){
    print('No calls for seurat found')
    dt <- data.table(CellBarcode = mc$Barcode, HTO = mc$HTO_classification, HTO_Classification = mc$HTO_classification.global, key = 'CellBarcode', Seurat = c(F), MultiSeq = c(T))
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  if (all(is.na(mc))){
    print('No calls for MULTI-seq found')
    dt <- data.table(CellBarcode = sc$Barcode, HTO = sc$HTO_classification, HTO_Classification = sc$HTO_classification.global, key = 'CellBarcode', Seurat = c(T), MultiSeq = c(F))
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)
  }

  mc$Barcode <- as.character(mc$Barcode)
  sc$Barcode <- as.character(sc$Barcode)
  merged <- merge(mc, sc, all = T, by = 'Barcode', suffixes = c('.MultiSeq', '.Seurat'))

  merged$HTO_classification.MultiSeq[is.na(merged$HTO_classification.MultiSeq)] <- 'Negative'
  merged$HTO_classification.Seurat[is.na(merged$HTO_classification.Seurat)] <- 'Negative'

  merged$HTO_classification.global.MultiSeq[is.na(merged$HTO_classification.global.MultiSeq)] <- 'Negative'
  merged$HTO_classification.global.Seurat[is.na(merged$HTO_classification.global.Seurat)] <- 'Negative'

  merged$Concordant <- as.character(merged$HTO_classification.MultiSeq) == as.character(merged$HTO_classification.Seurat)
  merged$ConcordantNoNeg <- !(!merged$Concordant & merged$HTO_classification.MultiSeq != 'Negative' & merged$HTO_classification.Seurat != 'Negative')
  merged$GlobalConcordant <- as.character(merged$HTO_classification.global.MultiSeq) == as.character(merged$HTO_classification.global.Seurat)
  merged$HasSeuratCall <- !is.na(merged$HTO_classification.Seurat) & merged$HTO_classification.Seurat != 'Negative'
  merged$HasMultiSeqCall <- !is.na(merged$HTO_classification.MultiSeq) & merged$HTO_classification.MultiSeq != 'Negative'

  merged$Seurat <- merged$HTO_classification.Seurat != 'Negative'
  merged$MultiSeq <- merged$HTO_classification.MultiSeq != 'Negative'

  print(paste0('Total concordant: ', nrow(merged[merged$Concordant])))
  print(paste0('Total discordant: ', nrow(merged[!merged$Concordant])))
  print(paste0('Total discordant, ignoring negatives: ', nrow(merged[!merged$ConcordantNoNeg])))
  print(paste0('Total discordant global calls: ', nrow(merged[!merged$GlobalConcordant])))

  discord <- merged[!merged$GlobalConcordant]
  discord <- discord %>% group_by(HTO_classification.global.MultiSeq, HTO_classification.global.Seurat) %>% summarise(Count = n())

  print(qplot(x=HTO_classification.global.MultiSeq, y=HTO_classification.global.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By Global Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  discord <- merged[!merged$Concordant]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = n())

  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By HTO Call') + ylab('Seurat') + xlab('MULTI-seq')
  )

  discord <- merged[!merged$ConcordantNoNeg]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = n())

  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By HTO Call, Ignoring Negatives') + ylab('Seurat') + xlab('MULTI-seq')
  )
  ret <-merged[merged$ConcordantNoNeg,]

  # These calls should be identical, except for possibly negatives from one method that are non-negative in the other
  # For the time being, accept those as correct.
  ret$FinalCall <- ret$HTO_classification.MultiSeq
  ret$FinalCall[ret$HTO_classification.MultiSeq == 'Negative'] <- ret$HTO_classification.Seurat[ret$HTO_classification.MultiSeq == 'Negative']

  ret$FinalClassification <- ret$HTO_classification.global.MultiSeq
  ret$FinalClassification[ret$HTO_classification.global.MultiSeq == 'Negative'] <- ret$HTO_classification.global.Seurat[ret$HTO_classification.global.MultiSeq == 'Negative']

  if (!is.na(allCallsOutFile) && nrow(merged) > 0) {
    write.table(merged, file = allCallsOutFile, row.names = F, sep = '\t', quote = F)
  }

  if (nrow(ret) > 0){
    dt <- data.table(CellBarcode = ret$Barcode, HTO = ret$FinalCall, HTO_Classification = ret$FinalClassification, key = 'CellBarcode', Seurat = ret$HasSeuratCall, MultiSeq = ret$HasMultiSeqCall)
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)

    return(dt)

  } else {
    print('No rows, not saving ')
  }
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
printFinalSummary <- function(dt, barcodeData){
  #Append raw counts:
  bc <- t(barcodeData)
  x <- melt(bc)
  names(x) <- c('CellBarcode', 'HTO', 'Count')

  merged <- merge(dt, x, by = c('CellBarcode', 'HTO'), all.x = T, all.y = F)

  bc <- as.data.frame(bc)
  bc$CellBarcode <- rownames(bc)
  merged <- merge(merged, bc, by = c('CellBarcode'), all.x = T, all.y = F)

  merged$HTO[is.na(merged$HTO)] <- c('Negative')
  merged$HTO_Classification[is.na(merged$HTO_Classification)] <- c('Negative')

  #summarize reads by type:
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

  print(kable(table(Classification = merged$HTO)))

  print(ggplot(merged, aes(x = HTO_Classification)) +
          geom_bar(stat = 'count') +
          xlab('Classification') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )

  print(kable(table(Classification = merged$HTO_Classification)))

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

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
simplifyHtoNames <- function(v) {
  return(sapply(v, function(x){
    x <- unlist(strsplit(x, '-'))
    if (length(x) > 1) {
      x <- x[-(length(x))]
    }

    paste0(x, collapse = "-")
  }))
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
HTODemux2 <- function(

  object,
  assay = "HTO",
  positive.quantile = 0.99,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  verbose = TRUE
) {
  # This is a hack around Seurat's method.  If this is improved, shift to use that:

  #initial clustering
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = GetAssayData(object = object, assay = assay)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
        x = t(x = GetAssayData(object = object, assay = assay)),
        k = ncenters,
        samples = nsamples
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    },
    stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
  )
  #average hto signals per cluster
  #work around so we don't average all the RNA levels which takes time
  average.expression <- AverageExpression(
    object = object,
    assay = assay,
    verbose = FALSE
  )[[assay]]

  #TODO: checking for any HTO negative in all clusters:

  #if (sum(average.expression == 0) > 0) {
  #  stop("Cells with zero counts exist as a cluster.")
  #}

  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]

    minNonZero <- which.min(x = average.expression[iter,average.expression[iter, ] > 0])
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[[minNonZero]]
    )]
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    }
  }
  # now assign cells to HTO based on discretized values
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
    }
  )])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
  )
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )
  colnames(x = classification.metadata) <- paste(
    assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
    sep = '_'
  )
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  # object@meta.data$hash.ID <- Idents(object)
  object$hash.ID <- Idents(object = object)
  return(object)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
generateSummaryForExpectedBarcodes <- function(dt, whitelistFile, outputFile, barcodeData) {
  categoryName <- "Cell Hashing Concordance"

  whitelist <- read.table(whitelistFile, sep = '\t', header = F)
  names(whitelist) <- c('CellBarcode')
  df <- data.frame(Category = categoryName, MetricName = "InputBarcodes", Value = length(whitelist$CellBarcode))


  df <- rbind(df, data.frame(Category = categoryName, MetricName = "TotalCounts", Value = sum(barcodeData)))

  #Any called cell:
  calledCellBarcodes <- dt$CellBarcode[dt$HTO_Classification != 'Negative']
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

  #By caller:
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "SeuratNonNegative", Value = sum(dt$Seurat)))
  df <- rbind(df, data.frame(Category = categoryName, MetricName = "MultiSeqNonNegative", Value = sum(dt$MultiSeq)))

  df$Value[is.na(df$Value)] <- 0

  write.table(df, file = outputFile, quote = F, row.names = F, sep = '\t')
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
downloadAndAppendCellHashing <- function(seuratObject, outPath = '.'){
  if (is.null(seuratObject[['BarcodePrefix']])){
    stop('Seurat object lacks BarcodePrefix column')
  }

  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    print(paste0('Adding cell hashing data for prefix: ', barcodePrefix))

    cellHashingId <- findMatchedCellHashing(barcodePrefix)
    if (is.na(cellHashingId)){
      stop(paste0('Unable to find cellHashing calls table file for prefix: ', barcodePrefix))
    }

    callsFile <- file.path(outPath, paste0(barcodePrefix, '_cellHashingCalls.csv'))
    downloadOutputFile(outputFileId = cellHashingId, outFile = callsFile, overwrite = T)
    if (!file.exists(callsFile)){
      stop(paste0('Unable to download calls table for prefix: ', barcodePrefix))
    }

    seuratObject <- appendCellHashing(seuratObj = seuratObject, barcodeCallFile = callsFile, barcodePrefix = barcodePrefix)
  }

  return(seuratObject)
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
findMatchedCellHashing <- function(loupeDataId){
  #Note: the seurat object gets associated with the GEX readset, so look based on this:
  rows <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="readset",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    return(NA)
  }

  readset <- rows[['readset']]
  if (is.na(readset) || is.null(readset)) {
    return(NA)
  }

  rows <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,",
    colFilter=makeFilter(c("readset", "EQUAL", readset), c("category", "EQUAL", "Seurat Cell Hashing Calls")),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) > 1){
    rows <- rows[1]
  }

  return(rows[['rowid']])
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
downloadOutputFile <- function(outputFileId, outFile, overwrite = T) {
  #There should be a file named all_contig_annotations.csv in the same directory as the VLoupe file
  rows <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
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
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/",wb),
    remoteFilePath = remotePath,
    overwrite = overwrite,
    localFilePath = outFile
  )

  if (!success | !file.exists(outFile)) {
    stop(paste0('failed to download file: ', remotePath))
  }

  return(outFile)
}
