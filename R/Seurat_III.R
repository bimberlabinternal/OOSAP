#' @import Seurat
#' @import Rlabkey

Rlabkey::labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")

utils::globalVariables(
  names = c('nCount_RNA', 'nFeature_RNA', 'p.mito'),
  package = 'OOSAP',
  add = TRUE
)


#' @title Read and Filter 10X files.
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords ReadAndFilter10X
#' @export
#' @importFrom Seurat Read10X
ReadAndFilter10xData <- function(dataDir, datasetName) {
  if (!file.exists(dataDir)){
    stop(paste0("File does not exist: ", dataDir))
  }

  if (!dir.exists(dataDir)){
    stop(paste0("File is not a directory: ", dataDir))
  }

  seuratRawData <- Read10X(data.dir = dataDir)
  seuratRawData <- PerformEmptyDropletFiltering(seuratRawData)

  seuratObj <- CreateSeuratObj(seuratRawData, project = datasetName)
  PrintQcPlots(seuratObj)

  return(seuratObj)
}



#' @title Create a Seurat 3 object
#'
#' @description Create Seurat Object from Read10X().
#' @param seuratData, A Seurat input data from Read10X().
#' @param project, Sets the project name for the Seurat object.
#' @param minFeatures, Include cells where at least this many features are detected.
#' @param minCells, Include features detected in at least this many cells.
#' @return A Seurat object with p.mito calculated.
#' @keywords CreateSeuratObj
#' @importFrom Matrix colSums
CreateSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25, minCells = 0, MitoGenesPattern = "^MT-"){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = minCells, min.features = minFeatures, project = project)

  mito.features <- grep(pattern = MitoGenesPattern, x = rownames(x = seuratObj), value = TRUE)
  p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['p.mito']] <- p.mito

  return(seuratObj)
}



#' @title A Title
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object
#' @importFrom Matrix colSums
PrintQcPlots <- function(seuratObj) {
  print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "p.mito"), ncol = 3))
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "p.mito"))
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

  #10x-like plot
  nUMI <- Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))
  nUMI <- sort(nUMI)

  countAbove <-unlist(lapply(nUMI, function(x){
    sum(nUMI >= x)
  }))

  plot(log(countAbove), log(nUMI), pch=20, ylab = "UMI/Cell", xlab = "# Cells")
}



#' @importFrom DropletUtils barcodeRanks
PerformEmptyDropletFiltering <- function(seuratRawData, fdrThreshold=0.01, emptyDropNIters=10000) {
  br.out <- DropletUtils::barcodeRanks(seuratRawData)

  # Making a plot.
  plot(br.out$rank, br.out$total+1, log="xy", xlab="Rank", ylab="Total")

  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=br.out$knee, col="dodgerblue", lty=2)
  abline(h=br.out$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection")
  )

  e.out <- PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters, fdrThreshold = fdrThreshold)

  toPlot <- e.out[is.finite(e.out$LogProb),]
  if (nrow(toPlot) > 0) {
    plot(toPlot$Total, -toPlot$LogProb, col=ifelse(toPlot$is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability")
  } else {
    print('Probabilities all -Inf, unable to plot')
  }

  if (nrow(toPlot) != nrow(e.out)) {
    print(paste0('Total rows with non-finite probabilities: ', (nrow(e.out) - nrow(toPlot))))
  }

  passingCells <- rownames(e.out)[e.out$is.cell]

  return(seuratRawData[,passingCells])
}



#' @title PerformEmptyDrops
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @importFrom DropletUtils emptyDrops
PerformEmptyDrops <- function(seuratRawData, emptyDropNIters, fdrThreshold=0.01){
  print(paste0('Performing emptyDrops with ', emptyDropNIters, ' iterations'))

  e.out <- DropletUtils::emptyDrops(seuratRawData, niters = emptyDropNIters)

  print(paste0('Input cells: ', nrow(e.out)))
  e.out <- e.out[!is.na(e.out$LogProb),]
  e.out$is.cell <- e.out$FDR <= fdrThreshold
  print(paste0('Cells passing FDR: ', sum(e.out$is.cell, na.rm=TRUE)))
  print(paste0('Cells failing FDR: ', sum(!e.out$is.cell, na.rm=TRUE)))

  #If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
  print(table(Limited=e.out$Limited, Significant=e.out$is.cell))
  totalLimited <- sum(e.out$Limited[e.out$Limited == T] & e.out$Significant == F)
  if (totalLimited == 0){
    return(e.out)
  } else {
    return(PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters * 2, fdrThreshold = fdrThreshold))
  }
}



HasStepRun <- function(seuratObj, name) {
  return(!is.null(seuratObj@misc[[paste0(name, 'Run')]]))
}



MarkStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)
  }

  return(seuratObj)
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
MergeSeuratObjs <- function(seuratObjs, metadata=NULL, alignData = T, MaxCCAspaceDim = 20, MaxPCs2Weight = 20, projectName = NULL, PreProcSeur = F, useAllFeatures = F, nVariableFeatures = 2000, includeCellCycleGenes = T){
  nameList <- ifelse(is.null(metadata), yes = names(seuratObjs), no = names(metadata))

  for (exptNum in nameList) {
    print(paste0('adding dataset: ', exptNum))
    prefix <- paste0(exptNum)
    so <- seuratObjs[[exptNum]]

    if (!('BarcodePrefix' %in% names(so@meta.data))) {
      so <- RenameCells(object = so, add.cell.id = prefix)
      so[['BarcodePrefix']] <- c(prefix)
    } else {
      print('Barcode prefix already added')
    }

    if (alignData && length(seuratObjs) > 1) {
      if (!HasStepRun(so, 'NormalizeData')) {
        print('Normalizing')
        so <- NormalizeData(object = so, verbose = F)
      } else {
        print('Normalization performed')
      }

      if (!HasStepRun(so, 'FindVariableFeatures')) {
        print('Finding variable features')
        so <- FindVariableFeatures(object = so, verbose = F, selection.method = "vst", nfeatures = nVariableFeatures)
      } else {
        print('FindVariableFeatures performed')
      }
    }

    seuratObjs[[exptNum]] <- so
  }

  seuratObj <- NULL
  if (alignData && length(seuratObjs) > 1) {
    # dims here means : Which dimensions to use from the CCA to specify the neighbor search space
    anchors <- FindIntegrationAnchors(object.list = seuratObjs, dims = 1:MaxCCAspaceDim, scale = !PreProcSeur, verbose = F)

    #always run using intersection of all features
    features <- NULL
    if (useAllFeatures) {
      features <- rownames(seuratObjs[[1]])
      if (length(seuratObjs) > 1) {
        for (i in 2:length(seuratObjs)) {
          features <- c(features, rownames(seuratObjs[[i]]))
        }

      }

      print(paste0('Total features in common: ', length(features)))
    } else if (includeCellCycleGenes) {
        features <- unique(c(cc.genes, g2m.genes.orig))
        if (length(seuratObjs) > 1) {
            for (i in 2:length(seuratObjs)) {
                features <- intersect(features, rownames(seuratObjs[[i]]))
            }
        }

        #Merge with the default set Seurat will use
        features <- unique(c(features, slot(object = anchors, name = "anchor.features")))
    }

    # dims here means : #Number of PCs to use in the weighting procedure
    seuratObj <- IntegrateData(anchorset = anchors, dims = 1:MaxPCs2Weight, verbose = F, features.to.integrate = features, new.assay.name = "Integrated")
    DefaultAssay(seuratObj) <- "Integrated"

    # This will prevent repeating this step downstream
    seuratObj <- MarkStepRun(seuratObj, 'NormalizeData')
  }
  else {
    for (exptNum in nameList) {
      if (is.null(seuratObj)) {
        seuratObj <- seuratObjs[[exptNum]]
      } else {
        seuratObj <- merge(x = seuratObj,
                           y = seuratObjs[[exptNum]],
                           project = projectName,
                           do.normalize = F)
      }
    }

    print('after merge')
    print(seuratObj)
  }

  return(seuratObj)
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
ProcessSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
                           nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                           nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf, forceReCalc = F,
                           variableGeneTable = NULL, variableFeatureSelectionMethod = 'vst', nVariableFeatures = 2000, printDefaultPlots = T,
                           npcs = 50){

  if (doCellFilter & (forceReCalc | !HasStepRun(seuratObj, 'FilterCells'))) {
    print("Filtering Cells...")
    seuratObj@misc$OriginalCells <- length(colnames(x = seuratObj))
    seuratObj <- subset(x = seuratObj, subset = nCount_RNA > nGene.low & nCount_RNA < nGene.high)
    seuratObj <- subset(x = seuratObj, subset = nFeature_RNA > nUMI.low & nFeature_RNA < nUMI.high)
    seuratObj <- subset(x = seuratObj, subset = p.mito > pMito.low & p.mito < pMito.high)

    print(paste0('Initial cells: ', seuratObj@misc$OriginalCells, ', after filter: ', length(colnames(x = seuratObj))))

    seuratObj <- MarkStepRun(seuratObj, 'FilterCells')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'NormalizeData')) {
    seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)
    seuratObj <- MarkStepRun(seuratObj, 'NormalizeData', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'FindVariableFeatures')) {
    seuratObj <- FindVariableFeatures(object = seuratObj, mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), verbose = F, selection.method = variableFeatureSelectionMethod, nVariableFeatures = NULL)
    seuratObj <- MarkStepRun(seuratObj, 'FindVariableFeatures', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ScaleData')) {
    seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "percent.mito"), display.progress = F, verbose = F)
    seuratObj <- MarkStepRun(seuratObj, 'ScaleData')
  }

  if (doCellCycle & (forceReCalc | !HasStepRun(seuratObj, 'CellCycle'))) {
    seuratObj <- RemoveCellCycle(seuratObj)
    seuratObj <- MarkStepRun(seuratObj, 'CellCycle', saveFile)
  }

  vg <- VariableFeatures(object = seuratObj)
  if (forceReCalc | !HasStepRun(seuratObj, 'RunPCA')) {
    seuratObj <- RunPCA(object = seuratObj, features = vg, verbose = F, npcs = npcs)
    seuratObj <- MarkStepRun(seuratObj, 'RunPCA')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ProjectDim')) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- MarkStepRun(seuratObj, 'ProjectDim')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
    seuratObj <- MarkStepRun(seuratObj, 'JackStraw', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- MarkStepRun(seuratObj, 'ScoreJackStraw')
  }

  print(paste0('Total variable genes: ', length(vg)))
  if (!is.null(variableGeneTable)){
    write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
  }

  if (printDefaultPlots){
    print(VizDimLoadings(object = seuratObj, dims = 1:2))
    print(DimPlot(object = seuratObj))

    print(DimHeatmap(object = seuratObj, dims = 1, cells = 500, balanced = TRUE, fast = F) + NoLegend())
    print(DimHeatmap(object = seuratObj, dims = 1:20, cells = 500, balanced = TRUE, fast = F) + NoLegend())

    print(JackStrawPlot(object = seuratObj, dims = 1:20))
    print(ElbowPlot(object = seuratObj))
  }

  return(seuratObj)
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendTcrClonotypes <- function(seuratObject, outPath = '.', dropExisting = T, metaFeat = NULL){
  #metaFeat = NULL is the default behavior colname(seuratObject)
  #         = "metaFeat" means use the specified metadata column instead
  
  if (is.null(seuratObject[['BarcodePrefix']])){
    stop('Seurat object lacks BarcodePrefix column')
  }

  i <- 0
  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    i <- i + 1
    print(paste0('Adding TCR clonotypes for prefix: ', barcodePrefix))

    vloupeId <- FindMatchedVloupe(barcodePrefix)
    if (is.na(vloupeId)){
      stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
    }

    clonotypeFile <- file.path(outPath, paste0(barcodePrefix, '_clonotypes.csv'))
    DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = T)
    if (!file.exists(clonotypeFile)){
      stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
    }

    doDropExisting <- i == 1 && dropExisting
    seuratObject <- AppendTcrClonotypes(seuratObject, clonotypeFile, barcodePrefix = barcodePrefix, dropExisting = doDropExisting, metaFeat = metaFeat)
  }

  return(seuratObject)
}


AppendTcrClonotypes <- function(seuratObject = NA, clonotypeFile = NA, barcodePrefix = NULL, dropExisting = F, metaFeat = NULL){
  tcr <- ProcessAndAggregateTcrClonotypes(clonotypeFile)

  if (!is.null(barcodePrefix)){
    tcr$barcode <- as.character(tcr$barcode)
    tcr$barcode <- paste0(barcodePrefix, '_', tcr$barcode)
    tcr$barcode <- as.factor(tcr$barcode)
  }

  origRows <- nrow(tcr)

  datasetSelect <- seuratObject$BarcodePrefix == barcodePrefix
  #metaFeat = NULL is the default behavior colname(seuratObject)
  #         = "metaFeat" means use the specified metadata column instead
  if(is.null(metaFeat)) gexBarcodes <- colnames(seuratObject)[datasetSelect] else {
    gexBarcodes <- as.vector(seuratObject@meta.data[,metaFeat])
    #gexBarcodes <- rownames(seuratObject@meta.data)[datasetSelect]
  }
  # This change is because, one cannot change the cell names of a Seurat object, but the metadata DF rownames cange be changed. 
  # gexBarcodes <- rownames(seuratObject@meta.data)[datasetSelect]
  
  tcr <- tcr[tcr$barcode %in% gexBarcodes,]
  pct <- nrow(tcr) / origRows * 100

  print(paste0('Barcodes with clonotypes: ', origRows, ', intersecting with GEX data: ', nrow(tcr), " (", pct, "%)"))

  merged <- merge(data.frame(barcode = gexBarcodes, sortOrder = 1:length(gexBarcodes)), tcr, by = c('barcode'), all.x = T)
  rownames(merged) <- merged$barcode
  merged <- dplyr::arrange(merged, sortOrder)
  merged <- merged[colnames(merged) != 'sortOrder']

  # Check barcodes match before merge
  if (sum(merged$barcode != gexBarcodes) > 0) {
    #stop(paste0('Seurat and TCR barcodes do not match after merge, total different: ', sum(merged$barcode != colnames(seuratObject)[datasetSelect])))
    stop(paste0('Seurat and TCR barcodes do not match after merge, total different: ', sum(merged$barcode != gexBarcodes )))
    
  }

  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    toAdd <- as.character(merged[[colName]])
    names(toAdd) <- merged[['barcode']]

    if ((colName %in% names(seuratObject@meta.data)) && dropExisting) {
      seuratObject@meta.data[colName] <- NULL
    }

    if (!(colName %in% names(seuratObject@meta.data))) {
      toUpdate <- rep(NA, ncol(seuratObject))
    } else {
      toUpdate <- unlist(seuratObject[[colName]])
    }

    names(toUpdate) <- colnames(seuratObject)
    toUpdate[datasetSelect] <- toAdd
    seuratObject[[colName]] <- as.factor(toUpdate)
  }

  return(seuratObject)
}


#' @import Rlabkey
FindMatchedVloupe <- function(loupeDataId) {
  rows <- labkey.selectRows(
    folderPath="/Labs/Bimber/",
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="readset/cdna/enrichedReadsetId",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeDataId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    return(NA)
  }

  tcrReadsetId <- rows[['readset_cdna_enrichedreadsetid']]
  if (is.na(tcrReadsetId) || is.null(tcrReadsetId)) {
    return(NA)
  }

  rows <- labkey.selectRows(
    folderPath="/Labs/Bimber/",
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,",
    colFilter=makeFilter(c("readset", "EQUAL", tcrReadsetId), c("category", "EQUAL", "10x VLoupe")),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) > 1){
    rows <- rows[1]
  }

  return(rows[['rowid']])
}


#' @import Rlabkey
DownloadCellRangerClonotypes <- function(vLoupeId, outFile, overwrite = T) {
  #There should be a file named all_contig_annotations.csv in the same directory as the VLoupe file
  rows <- labkey.selectRows(
    folderPath="/Labs/Bimber/",
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
    colFilter=makeFilter(c("rowid", "EQUAL", vLoupeId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    return(NA)
  }

  wb <- rows[['workbook_workbookid']]
  if (is.na(wb) || is.null(wb)){
    wb <- ''
  }

  remotePath <- rows[['dataid_webdavurlrelative']]
  remotePath <- paste0(dirname(remotePath), '/all_contig_annotations.csv')

  success <- labkey.webdav.get(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/",wb),
    remoteFilePath = remotePath,
    overwrite = overwrite,
    localFilePath = outFile
  )

  if (!success | !file.exists(outFile)) {
    return(NA)
  }

  return(outFile)
}



#' @import Rlabkey
#' @importFrom dplyr %>% coalesce group_by summarise
#' @importFrom naturalsort naturalsort
ProcessAndAggregateTcrClonotypes <- function(clonotypeFile){
  tcr <- utils::read.table(clonotypeFile, header=T, sep = ',')
  tcr <- tcr[tcr$cdr3 != 'None',]

  # drop cellranger '-1' suffix
  tcr$barcode <- gsub("-1", "", tcr$barcode)

  #Download named clonotypes and merge:
  # Add clone names:
  labelDf <- labkey.selectRows(
    folderPath="/Labs/Bimber/",
    schemaName="tcrdb",
    queryName="clones",
    showHidden=TRUE,
    colSelect=c('clonename','chain','cdr3','animals', 'displayname', 'vgene'),
    containerFilter=NULL,
    colNameOpt='rname'
  )

  labelDf$LabelCol <- coalesce(as.character(labelDf$displayname), as.character(labelDf$clonename))

  labelDf <- labelDf %>%
    group_by(chain, cdr3) %>%
    summarise(CloneName = paste0(sort(unique(LabelCol)), collapse = ","))

  tcr <- merge(tcr, labelDf, by.x = c('chain', 'cdr3'), by.y = c('chain', 'cdr3'), all.x = TRUE, all.y = FALSE)

  # Many TRDV genes can be used as either alpha or delta TCRs.  10x classifies and TRDV/TRAJ/TRAC clones as 'Multi'.  Re-classify these:
  tcr$chain[tcr$chain == 'Multi' & grepl(pattern = 'TRD', x = tcr$v_gene) & grepl(pattern = 'TRAJ', x = tcr$j_gene) & grepl(pattern = 'TRAC', x = tcr$c_gene)] <- c('TRA')

  # Add chain-specific columns:
  tcr$ChainCDR3s <- paste0(tcr$chain, ':', tcr$cdr3)
  for (l in c('TRA', 'TRB', 'TRD', 'TRG')){
    tcr[[l]] <- c(NA)
    tcr[[l]][tcr$chain == l] <- as.character(tcr$cdr3[tcr$chain == l])

    target <- paste0(l, 'V')
    tcr[[target]] <- c(NA)
    tcr[[target]][tcr$chain == l] <- as.character(tcr$v_gene[tcr$chain == l])
  }

  # Summarise, grouping by barcode
  tcr <- tcr %>% group_by(barcode) %>% summarise(
    ChainCDR3s = paste0(sort(unique(ChainCDR3s)), collapse = ","),
    CDR3s = paste0(sort(unique(cdr3)), collapse = ","),
    TRA = paste0(sort(unique(as.character(TRA))), collapse = ","),
    TRB = paste0(sort(unique(as.character(TRB))), collapse = ","),
    TRD = paste0(sort(unique(as.character(TRD))), collapse = ","),
    TRG = paste0(sort(unique(as.character(TRG))), collapse = ","),
    TRAV = paste0(sort(unique(as.character(TRAV))), collapse = ","),
    TRBV = paste0(sort(unique(as.character(TRBV))), collapse = ","),
    TRDV = paste0(sort(unique(as.character(TRDV))), collapse = ","),
    TRGV = paste0(sort(unique(as.character(TRGV))), collapse = ","),

    CloneNames = paste0(sort(unique(CloneName)), collapse = ",")  #this is imprecise b/c we count a hit if we match any chain, but this is probably what we often want
  )

  # Note: we should attempt to produce a more specfic call, assuming we have data from multiple chains
  # The intent of this was to allow a A- or B-only hit to produce a call, but if we have both A/B, take their intersect.

  tcr$CloneNames <- sapply(strsplit(as.character(tcr$CloneNames), ",", fixed = TRUE), function(x) paste0(naturalsort(unique(x)), collapse = ","))

  tcr$barcode <- as.factor(tcr$barcode)
  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    v <- tcr[[colName]]
    v <- as.character(v)
    v[v == ''] <- NA

    tcr[[colName]] <- as.factor(v)
  }

  return(tcr)
}



#' @title RemoveCellCycle
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @importFrom cowplot plot_grid
RemoveCellCycle <- function(seuratObj, runPCAonVariableGenes = F) {
  print("Performing cell cycle cleaning ...")

  # Cell cycle genes were obtained from the Seurat example (See regev_lab_cell_cycle_genes.txt)
  # and stored using use_data(internal = T) (https://github.com/r-lib/usethis and use_data)
  # cc.genes
  # g2m.genes.orig
  if (length(cc.genes) != 97) {
    stop('Something went wrong loading cc.genes list')
  }

  if (length(g2m.genes.orig) != 200) {
    stop('Something went wrong loading g2m.genes list')
  }

  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- cc.genes[1:43]
  g2m.genes <- unique(c(g2m.genes.orig, cc.genes[44:97]))

  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]

  if (length(g2m.genes) < 20 || length(s.genes) < 20) {
    print(paste0("Warning, the number of g2m (", length(g2m.genes), ") and/or s genes (", length(s.genes), ") in your data is low"))
  }

  #proceeds <20 but warns, but <5 is fishy and cant use
  if (length(g2m.genes) < 5 || length(s.genes) < 5) {
    print("Error, the number of g2m and/or s genes < 5")
    #break()
    seuratObj <- MarkStepRun(seuratObj, 'FAIL_RemoveCellCycle')
    return(seuratObj) # for pipeline not breaking,... but
  }

  print("Running PCA with cell cycle genes")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE, verbose = F)
  print(DimPlot(object = seuratObj, reduction = "pca"))

  #store values to append later
  SeuratObjsCCPCA <- as.data.frame(seuratObj@reductions$pca@cell.embeddings)
  colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")

  # Eisa: do we still need a custom version of this?
  seuratObj <- CellCycleScoring(object = seuratObj,
                                       s.features = s.genes,
                                       g2m.features = g2m.genes,
                                       set.ident = TRUE)

  print(cowplot::plot_grid(plotlist = list(DimPlot(object = seuratObj, reduction = "pca", dims = c(1, 2)),
                                           DimPlot(object = seuratObj, reduction = "pca", dims = c(2, 3)),
                                           DimPlot(object = seuratObj, reduction = "pca", dims = c(3, 4)),
                                           DimPlot(object = seuratObj, reduction = "pca", dims = c(4, 5)) ))
  )

  print("Regressing out S and G2M score ...")
  seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = F, verbose = F)

  #Note: normally this is run upstream, which supports more options for how variable genes are defined.
  if (runPCAonVariableGenes) {
    print("Running PCA with variable genes ...")
    seuratObj <- RunPCA(object = seuratObj, pc.genes = VariableFeatures(object = seuratObj), do.print = F, verbose = F)
  }

  for (colName in colnames(SeuratObjsCCPCA)) {
    seuratObj[[colName]] <- SeuratObjsCCPCA[colnames(seuratObj),colName]
  }

  return(seuratObj)
}



#' @title FindClustersAndDimRedux
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
FindClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, saveFile = NULL, forceReCalc = F) {
  if (is.null(dimsToUse)) {
    elbow <- FindSeuratElbow(seuratObj)
    print(paste0('Inferred elbow: ', elbow))

    dimsToUse <- 1:elbow
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'FindNeighbors')) {
    seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
    seuratObj <- MarkStepRun(seuratObj, 'FindNeighbors')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'FindClusters')) {
    for (resolution in c(0.2, 0.4, 0.8, 1.2, 0.6)){
      seuratObj <- FindClusters(object = seuratObj, resolution = resolution)
      seuratObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = seuratObj, verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'FindClusters', saveFile)
    }
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'RunTSNE')) {
    # To avoid Rtsne 'perplexity is too large for the number of samples' error.
    # See check here: https://github.com/jkrijthe/Rtsne/blob/ebed20f612712987fd160386132c17289169b4d8/R/utils.R
    # See also: https://github.com/satijalab/seurat/issues/167
    # Infer the max allowable based on their formula
    perplexity <- 30
    if (ncol(seuratObj) - 1 < 3 * perplexity) {
      print(paste0('Perplexity is too large for the number of samples: ', ncol(seuratObj)))
      perplexityNew <- floor((ncol(seuratObj) - 1) / 3)
      print(paste0('lowering from ', perplexity, ' to: ', perplexityNew))
      perplexity <- perplexityNew
    }

    seuratObj <- RunTSNE(object = seuratObj, dims.use = dimsToUse, check_duplicates = FALSE, perplexity = perplexity)
    seuratObj <- MarkStepRun(seuratObj, 'RunTSNE', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'RunUMAP')) {
    seuratObj <- RunUMAP(seuratObj,
                           dims = dimsToUse,
                           n.neighbors = 40L,
                           min.dist = 0.2,
                           metric = "correlation",
                           seed.use = 1234)
    seuratObj <- MarkStepRun(seuratObj, 'RunUMAP', saveFile)
  }

  for (reduction in c('tsne', 'umap')){
    plot1 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.2", label = TRUE) + ggtitle('Resolution: 0.2')
    plot2 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.4", label = TRUE) + ggtitle('Resolution: 0.4')
    plot3 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.6", label = TRUE) + ggtitle('Resolution: 0.6')
    plot4 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.8", label = TRUE) + ggtitle('Resolution: 0.8')
    plot5 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_1.2", label = TRUE) + ggtitle('Resolution: 1.2')

    print(CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), legend = 'none'))
  }

  return(seuratObj)
}


#' @title FindMarkers
#' @importFrom dplyr %>% coalesce group_by summarise filter top_n
#' @export
FindMarkers <- function(seuratObj, resolutionToUse, outFile, saveFileMarkers = NULL,
                        testsToUse = c('wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'),
                        numGenesToSave = 20, onlypos = F) {

  Idents(seuratObj) <- seuratObj[[paste0('ClusterNames_',resolutionToUse)]]

  if (file.exists(saveFileMarkers)) {
    print('resuming from file')
    seuratObj.markers <- readRDS(saveFileMarkers)
  } else {
    seuratObj.markers <- NA
    for (test in testsToUse) {
      print(paste0('Running using test: ', test))
      tMarkers <- FindAllMarkers(object = seuratObj, only.pos = onlypos, min.pct = 0.25, logfc.threshold = 0.25, verbose = F, test.use = test)
      if (nrow(tMarkers) == 0) {
        print('No genes returned, skipping')
      } else {
        tMarkers$test <- c(test)
        tMarkers$cluster <- as.character(tMarkers$cluster)

        toBind <- data.frame(test = tMarkers$test,
                             cluster = tMarkers$cluster,
                             gene = tMarkers$gene,
                             avg_logFC = tMarkers$avg_logFC,
                             pct.1 = tMarkers$pct.1,
                             pct.2 = tMarkers$pct.2,
                             p_val_adj = tMarkers$p_val_adj
        )
        if (all(is.na(seuratObj.markers))) {
          seuratObj.markers <- toBind
        } else {
          seuratObj.markers <- rbind(seuratObj.markers, toBind)
        }
      }
    }

    if (!('cluster' %in% names(seuratObj.markers))) {
      warning('cluster column not found!')
    } else {
      seuratObj.markers$cluster <- as.factor(seuratObj.markers$cluster)
    }

    saveRDS(seuratObj.markers, file = saveFileMarkers)
  }

  if (nrow(seuratObj.markers) == 0) {
    print('No significant markers were found')
    return()
  } else {
    toWrite <- seuratObj.markers %>% filter(p_val_adj < 0.001) %>% filter(avg_logFC > 0.5) %>% group_by(cluster, test) %>% top_n(numGenesToSave, avg_logFC)
    write.table(toWrite, file = outFile, sep = '\t', row.names = F, quote = F)

    if (nrow(toWrite) == 0) {
      print('No significant markers were found')
    } else {
      print(DimPlot(object = seuratObj, reduction = 'tsne'))

      topGene <- seuratObj.markers %>% filter(p_val_adj < 0.001) %>% filter(avg_logFC > 0.5) %>% group_by(cluster, test) %>% top_n(20, avg_logFC)
      print(DoHeatmap(object = seuratObj, features = unique(topGene$gene)))
    }
  }
}


#' @title FindSeuratElbow
#' @import ggplot2
FindSeuratElbow <- function(object, ndims = 25, reduction = "pca", print.plot = T, min.y = 1.3) {
  data.use <- Stdev(object = object, reduction = reduction)

  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use),
            " reductions")
    ndims <- length(x = data.use)
  }

  #1 sd = 1.3
  elbowX <- try(FindElbow(data.use[1:ndims], plot = T, ignore.concavity = F, min.y = min.y))
  if (class(elbowX)=="try-error" || elbowX[1]==2) {
    if (is.null(ndims)){
      elbowX = 2
    } else {
      elbowX = ndims
    }
  }

  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = "dims", y = "stdev")) +
    labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])),
         y = "Standard Deviation") + theme_bw() + geom_vline(xintercept = elbowX) + ggtitle("Elbow Identification")

  if (print.plot) {
    print(plot)
  }

  return(elbowX)
}

#' @title FindElbow
#' @importFrom stats coef
FindElbow <- function(y, plot = FALSE, ignore.concavity = FALSE, min.y = NA, min.x = NA) {

  # minor modification to debug specic scenarios when fail to find elbow
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.

  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  #' @examples
  #' \dontrun{
  #' tmp <- FindElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
  #' 	plot = TRUE) # wandering
  #' tmp <- FindElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
  #' 	plot = TRUE) # late rise
  #' tmp <- FindElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
  #' 	plot = TRUE) # gradual, no obvious break
  #'
  #' # Not the usual way to choose the number of PCs:
  #' library("chemometrics")
  #' data(glass)
  #' pca <- prcomp(glass)
  #' eigensum <- sum(pca$sdev * pca$sdev)
  #' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
  #' cs <- cumsum(vv)
  #' tmp <- FindElbow(vv, plot = TRUE)
  #' tmp <- FindElbow(cs, plot = TRUE)
  #' }

  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }

  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if (any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if (any(u < 0.00001) || any(u > 1)) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      #TODO: giving warning b/c length of ix/iy can be >1.  maybe if needs any() or all()??
      if (length(ix) > 1 || length(iy) > 1) {
        warning(paste0('length GT 1: ', length(ix), '/', length(iy)))
        print('ix:')
        print(ix)
        print('iy:')
        print(iy)
        print(x1)
        print(y1)
        print(x2)
        print(y2)
      }

      if (ix > iy)  ans <- iy
      else ans <- ix
    } else {
      ## Intersecting point is on the line, use the formula
      ix <- x1 + u * (x2 - x1)
      iy <- y1 + u * (y2 - y1)
      ans <- lineMagnitude(px, py, ix, iy)
    }
    ans
  }

  # End of helper functions by PB

  ### Now for the actual FindElbow function!

  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).


  y <- sort(y, decreasing = T)

  # Add an index to argument values for easy plotting
  DF <- data.frame(x = 1:length(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]

  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.

  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if (ignore.concavity) concave <- TRUE

  if (!concave) {
    stop("Your curve doesn't appear to be concave")
  }

  # Calculate the orthogonal distances
  if (is.na(min.x)){
    if (!is.na(min.y)){
      if (!length(which(DF$y<=min.y))<1){
        min.x = min(DF[which(DF$y<=min.y), ]$x)
      } else {
        print("min.y greater than smallest y")
        min.x = 2
      }
    } else {
      print("min.x and min.y are NA")
      min.x = 2
    }
  }

  use     <- min.x:(nrow(DF)-1)
  elbowd  <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- rep(NA, nrow(DF))
  DF$dist[use]  <- elbowd # c(NA, elbowd, NA) # first & last points don't have a distance

  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
         main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
             DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")
    points(DF$x[edm], DF$y[edm], pch = 20)
  }

  if (is.na(which.max(DF$dist))) {
    #if all fails return 2
    return(2)
  } else {
    return(which.max(DF$dist))
  }
}


#' @title WriteSummaryMetrics
#' @export
#' @param SeurObj, A Seurat object.
WriteSummaryMetrics <- function(seuratObj, file) {
  df <- data.frame(Category = "Seurat", MetricName = "TotalCells", Value = ncol(seuratObj))
  df <- rbind(df, data.frame(Category = "Seurat", MetricName = "TotalFeatures", Value = nrow(seuratObj)))

  write.table(df, file = file, quote = F, row.names = F, sep = '\t')
}



#' @title WriteCellBarcodes
#' @description Writes a table of cell barcodes to the provided file
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
WriteCellBarcodes <- function(seuratObj, file) {
  df <- data.frame(CellBarcode = colnames(seuratObj))

  write.table(df, file = file, quote = F, row.names = F, sep = ',', col.names = F)
}



#' @title SaveDimRedux
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
#' @import data.table
SaveDimRedux <- function(seuratObj, reductions=c("pca", "tsne", "umap"),
                         file=NA, maxPCAcomps=10, returnResults=F){

  if (is.na(file)){
    stop("file is required")
  }

  tempDT <- data.table(cbind(1:ncol(seuratObj), colnames(seuratObj)))
  rownames(tempDT) <- colnames(seuratObj)
  colnames(tempDT) <- c("nID", "cID")

  if ("tsne" %in% reductions) {
    if (is.null(seuratObj@reductions$tsne)) print("tsne slot NULL") else {

    }
    tsneDT <- data.table(seuratObj@reductions$tsne@cell.embeddings)
    tsneDT$cID <- rownames(seuratObj@reductions$tsne@cell.embeddings)
    tempDT <- merge(tempDT, tsneDT, by="cID")
  }

  if ("pca" %in% reductions) {
    if (is.null(seuratObj@reductions$pca)) print("pca slot NULL") else {
      pcaDT <- data.table(seuratObj@reductions$pca@cell.embeddings[,1:maxPCAcomps])
      pcaDT$cID <- rownames(seuratObj@reductions$pca@cell.embeddings)
      tempDT <- merge(tempDT, pcaDT, by="cID")
    }

  }

  if ("umap" %in% reductions) {
    if (is.null(seuratObj@reductions$umap)) print("umap slot NULL") else {
      umapDT <- data.table(seuratObj@reductions$umap@cell.embeddings)
      umapDT$cID <- rownames(seuratObj@reductions$umap@cell.embeddings)
      tempDT <- merge(tempDT, umapDT, by="cID")
    }
  }

  print("saving DimRedux")

  write.csv(tempDT, file = file, row.names=TRUE)

  if (returnResults) {
    return(tempDT)
  }
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
AddTitleToMultiPlot <- function(plotGrid, title, relHeights = c(0.1, 1)) {
  return(cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title), plotGrid, ncol = 1, rel_heights = relHeights))
}



#' @title A Title
#'
#' @import ggplot2
ggUMAP <- function(object,
                          colFac = NULL,
                          col_vector=NULL,
                          ptSize=0.15, ptAlpha=0.5,
                          add2title="",
                          legendTitle="",
                          cells.use = NULL){

  #updated March/22/19 - EM :: use.cell is needed to subset cells
  #updated April/03/19 - EM :: debug


  datat.temp <- as.data.frame(object@reductions$umap@cell.embeddings)
  colnames(datat.temp) <- c("UMAP1", "UMAP2")

  if (!is.null(cells.use)){
    datat.temp <- datat.temp[cells.use, ]
  }


  if (!is.null(colFac)){

    if (!is.factor(colFac)) colFac <- factor(colFac)

    if (is.null(col_vector)) col_vector = gg_color_hue(length(levels(colFac)))


    myColors <- col_vector[1:length(levels(colFac))]
    names(myColors) <- levels(colFac)

    colScale <- scale_colour_manual(name = ifelse(legendTitle=="", "grp", legendTitle), values = myColors)


    datat.temp$colFac <- colFac

    ggplot(datat.temp, aes(UMAP1, UMAP2 , col=colFac)) +
      geom_point(size=ptSize, alpha = ptAlpha) +
      theme(legend.position = "bottom") +
      ggtitle(paste("Expression of factor on UMAP", add2title, sep="")) +
      theme_bw() + colScale + guides(colour = guide_legend(override.aes = list(size=8, alpha=1)))

  } else {

    ggplot(datat.temp, aes(UMAP1, UMAP2)) +
      geom_point(size=ptSize, alpha = ptAlpha) +
      theme(legend.position = "bottom") +
      ggtitle(paste("UMAP 2D Space", add2title, sep="")) +
      theme_bw()


  }
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
MakeSerObjs_10XFolders <- function(counts.path = NULL,
                                          min.cells = 0,
                                          min.genes = 0,
                                          ProjName="10X",
                                          save.path = NULL,
                                          returnList=F, path.exclude="raw", string.exclude=NULL){

  #this function searches for 10X folder.
  #this is an update from the MakeSerObjs_10XFolders() which was done in Seurat II



  # counts.path = "../../../Bimber/Expts/214/10X/bin/all"
  # min.cells = 0
  # min.genes = 0
  # ProjName="10X"
  # save.path = NULL
  # path.exclude="filtered"
  # string.exclude = c("/raw_gene_bc_matrices", "/cellRanger-3204293", " 10x Loupe File", "/TempOut")
  # returnList = F

  if (returnList) TempLS <- list()

  if (!is.null(counts.path)){

    if (is.null(save.path)) save.path <- counts.path

    exp.dirs <- list.files(counts.path, recursive = T, full.names = T, pattern = ".mtx")
    exp.dirs <- exp.dirs[!grepl(".gz", exp.dirs)]
    exp.dirs <- gsub("/matrix.mtx", "", exp.dirs)
    if (!path.exclude=="") exp.dirs <- exp.dirs[!grepl(path.exclude, exp.dirs)]

    FileNames2Save <- exp.dirs
    if (is.null(string.exclude)) string.exclude <- paste(counts.path, "/", sep="")

    if (!is.null(string.exclude)){
      string.exclude <- c(string.exclude, paste(counts.path, "/", sep=""))
      for (iN in 1:length(string.exclude)){
        FileNames2Save <- gsub(string.exclude[iN], "", FileNames2Save)
      }
    }

    print("Found files... here some..")
    head(exp.dirs)

    if (returnList) TempLS$exp.dirs <- exp.dirs

    for (xN in 1:length(exp.dirs)){
      # xN = 1
      # use the list.files below to make sure the expected files are there....
      #list.files(exp.dirs[xN], full.names = T, recursive = T)

      if (returnList) TempLS$SeuratObjs <- list()

      print(exp.dirs[xN])
      print(paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))

      if (!file.exists(paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))){

        SeuratObjs <- ReadAndFilter10xData(dataDir=exp.dirs[xN],
                                                  datasetName=paste(ProjName, FileNames2Save[xN], sep="-"))

        # print("Reading in 10X folder...")
        # Seurat10X  <- Read10X(data.dir = exp.dirs[xN])
        #
        # print("Converting to Seurat Obj....")
        # rownames(Seurat10X)
        #
        # SeuratObjs <- CreateSeuratObj(seuratData = Seurat10X,
        #                                      minCells = min.cells,   #genes expressed in >= 5 cells
        #                                      minFeatures = min.genes, #Keep all cells with at least 200 detected genes
        #                                      project = paste(ProjName, FileNames2Save[xN], sep="-"))



        if (returnList) TempLS$SeuratObjs[[basename(exp.dirs[xN])]] <- SeuratObjs

        print("saving...")
        saveRDS(SeuratObjs,
                paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))

      } else {
        print(paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))
        print("already Exists...")
      }



    }; remove(xN)

    if (returnList) return(TempLS)
  }

}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
PreProcess_SerObjs <- function(SerObj.path = NULL, SerObjRDSKey="SeuratObj.rds",
                                      ProjName="10X",
                                      save.path = NULL, save.fig.path = NULL,
                                      returnList = F, returnDimReduxOnly=T,
                                      save.fig = T,
                                      nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                                      nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf,
                                      fvg.x.low.cutoff = 0.01, fvg.x.high.cutoff = 4.5, fvg.y.cutoff = 1.5,
                                      KeepGene.LS =NULL,
                                      nDimPCA=NA, findPCAElbow = T,
                                      RemoveCellCycle=F,
                                      path2CCfiles="./data/CellCycle", doUMAP = T){







  if (returnList) TempLS <- list()

  if (!is.null(SerObj.path)){

    if (is.null(save.path)) {
      save.path <- paste(SerObj.path, "/SerProc", sep="")
      if (!dir.exists(save.path)) dir.create(save.path, recursive = T)
    }
    if (is.null(save.fig.path)) {
      save.fig.path <- paste(SerObj.path, "/SerProcFigs", sep="")
      if (!dir.exists(save.fig.path)) dir.create(save.fig.path, recursive = T)
    }



    all_RDS  <- list.files(SerObj.path, full.names = T, pattern = ".rds")
    SeurObj_RDS <-  all_RDS[grep(SerObjRDSKey, all_RDS)]

    print("Found files... examples:")
    print(head(SeurObj_RDS))

    if (returnList){
      TempLS$SeuratObjs <- list()
      TempLS$DimReduxComps <- list()
    }



    for (xN in 1:length(SeurObj_RDS)){
      # xN=1
      bName <- gsub("\\.", "", gsub( "-", "", gsub( "_", "-", gsub("SeuratObj.rds", "", basename(SeurObj_RDS[xN])))))


      if (!file.exists(paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))){



        print("Reading in...")
        print(SeurObj_RDS[xN])


        SeuratObjs <- readRDS(SeurObj_RDS[xN])

        ### call fx here....
        #SeuratObjs <- ProcessSeurat1(SeuratObjs,
        #                                    dispersion.cutoff = c(fvg.y.cutoff, Inf),
        #                                    mean.cutoff = c(fvg.x.low.cutoff, fvg.x.high.cutoff),
        #                                    saveFile = NULL, doCellFilter=T,
        #                                    RemoveCellCycle = F,
        #                                    nUMI.high = nUMI.high, nGene.high = nGene.high, pMito.high = pMito.high,
        #                                    nUMI.low = nUMI.low, nGene.low = nGene.low, pMito.low = pMito.low)

        # nUMI.high = nUMI.high, nGene.high = nGene.high, pMito.high = pMito.high,
        # nUMI.low = nUMI.low, nGene.low = nGene.low, pMito.low = pMito.low

        # saveRDS(SeuratObjs, "./tempSeuratObjs.rds")
        # SeuratObjs <- readRDS("./tempSeuratObjs.rds")




        # #Genes that dont map to a specific name
        # noGeneSYM <- rownames(SeuratObjs)[grepl(ENSMB.tag, rownames(SeuratObjs))]
        #
        # length(noGeneSYM)

        # write.table(noGeneSYM,
        #             "./10X/Rhesus_ENSMMUG.csv",
        #             sep=", ", , row.names = F, quote = F,
        #             col.names = F)









        if (is.na(nDimPCA) || findPCAElbow) {
          nDimPCA <- FindSeuratElbow(SeuratObjs)
        }

        print(VizDimLoadings(object = SeuratObjs, dims = 1:4))
        print(DimHeatmap(object = SeuratObjs, dims = 1:nDimPCA, cells = 200, balanced = TRUE))

        # print(JackStrawPlot(object = seuratObj, dims = 1:20))

        SeuratObjs <- FindClustersAndDimRedux(seuratObj = SeuratObjs, dimsToUse=1:nDimPCA)

        print("saving ...")
        saveRDS(SeuratObjs,
                paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))

        if (returnList) TempLS$SeuratObjs[[bName]] <- SeuratObjs

      } else {
        if (returnList) {

          #to return a full list of data, need to read what is already processed and saved...
          TempLS$SeuratObjs[[bName]] <- readRDS(
            paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))

          TempLS$DimReduxComps[[bName]] <- SaveDimRedux(seuratObj= TempLS$SeuratObjs[[bName]],
                                                               file = paste(save.path, "/", bName, "_DimReduxComps.csv", sep=""), returnResults=T)

        }

      }




    }

    if (returnList) {
      if (returnDimReduxOnly) {
        TempLS$SeuratObjs <- NULL
      }
      return(TempLS)
    }

  }
}




#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
WilcoxDETest <- function(
  object,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  object <- object[, rownames(x = group.info), drop = FALSE]
  my.sapply <- ifelse(
    test = verbose && PlanThreads() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- my.sapply(
    X = 1:nrow(x = object),
    FUN = function(x) {
      return(wilcox.test(object[x, ] ~ group.info[, "group"], ...)$p.value)
    }
  )
  return(data.frame(p_val, row.names = rownames(x = object)))
}

#' @export
GetXYDataFromPlot <- function(plot, cellNames) {
  xynames <- Seurat:::GetXYAesthetics(plot = plot)

  plot.data <- plot$data[cellNames, ]
  names(plot.data)[names(plot.data) == xynames$x] <- 'x'
  names(plot.data)[names(plot.data) == xynames$y] <- 'y'

  return(plot.data)
}

#' @export
#' @import ggplot2
AddClonesToPlot <- function(seuratObject, plot, colorShapes = "black") {
  cellNames <- colnames(seuratObject)[!is.na(seuratObject$CloneNames)]
  plot.data <- GetXYDataFromPlot(plot, cellNames)
  plot.data$CloneName <- seuratObject$CloneNames[!is.na(seuratObject$CloneNames)]

  plot <- plot + geom_point(
  mapping = aes_string(x = 'x', y = 'y', shape = 'CloneName'),
  color = colorShapes,
  data = plot.data

  )

  return(plot)
}

#' @export
FilterCloneNames <- function(seuratObject, minValue) {
  ct <- table(seuratObject$CloneNames)
  ct <- ct[ct < minValue]

  seuratObject$CloneNames[seuratObject$CloneNames %in% names(ct)] <- NA

  return(seuratObject)
}



#' @title AvgCellExprs
#'
#' @param SeurObj, A Seurat object.
#' @return A data.frame of avg expression per var
#' @importFrom Matrix rowMeans
#' @export
AvgCellExprs <- function(seuratObj, varName = "ClusterNames_0.2", Genes){
  # seuratObj = SERObjLS$LymphAxLN
  # Genes = MarkersOfInterest
  
  AvlLevels <- factor(as.character(FetchData(seuratObj, varName)[,1]))
  
  ClustLS <- list()
  
  for(lev in levels(AvlLevels)){
    print(lev)
    ClustLS[[lev]] <- as.data.frame(Matrix::rowMeans(GetAssayData(object = seuratObj, 
                                                                  features = Genes)[Genes, colnames(seuratObj)[which(AvlLevels==lev)]  ]))
  }
  
  ClustDF <- as.data.frame(ClustLS)
  colnames(ClustDF) <- paste0("clus", levels(AvlLevels))
  
  return(ClustDF)
}