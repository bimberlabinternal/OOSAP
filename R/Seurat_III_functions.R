
Rlabkey::labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")


#' @title Read and Filter 10X files
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using performEmptyDropletFiltering and returns a Seurat object.
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords ReadAndFilter10X
#' @export
#' @examples
#' readAndFilter10xData(dataDir="./10Xdata", datasetName = "my10X")
readAndFilter10xData <- function(dataDir, datasetName) {
  seuratRawData <- Read10X(data.dir = dataDir)
  seuratRawData <- performEmptyDropletFiltering(seuratRawData)

  seuratObj <- createSeuratObj(seuratRawData, project = datasetName)
  printQcPlots(seuratObj)

  return(seuratObj)
}



#' @title Create a Seurat object
#'
#' @description Create Seurat Object from Read10X().
#' @param seuratData, A Seurat input data from Read10X().
#' @param project, Sets the project name for the Seurat object.
#' @param minFeatures, Include cells where at least this many features are detected.
#' @param minCells, Include features detected in at least this many cells.
#' @return A Seurat object with p.mito calculated.
#' @keywords CreateSeuratObject
#' @export
#' @examples
createSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25, minCells = 0, MitoGenesPattern = "^MT-"){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = minCells, min.features = minFeatures, project = project)

  mito.features <- grep(pattern = MitoGenesPattern, x = rownames(x = seuratObj), value = TRUE)
  p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['p.mito']] <- p.mito

  return(seuratObj)
}



#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
printQcPlots <- function(seuratObj) {
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





#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
performEmptyDropletFiltering <- function(seuratRawData, fdrThreshold=0.01, emptyDropNIters=10000) {
  br.out <- barcodeRanks(seuratRawData)

  # Making a plot.
  plot(br.out$rank, br.out$total+1, log="xy", xlab="Rank", ylab="Total")

  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=br.out$knee, col="dodgerblue", lty=2)
  abline(h=br.out$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))

  e.out <- performEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters, fdrThreshold = fdrThreshold)

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



#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
performEmptyDrops <- function(seuratRawData, emptyDropNIters, fdrThreshold=0.01){
  print(paste0('Performing emptyDrops with ', emptyDropNIters, ' iterations'))

  e.out <- emptyDrops(seuratRawData, niters = emptyDropNIters)

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
    return(performEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters * 2, fdrThreshold = fdrThreshold))
  }
}



#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
hasStepRun <- function(seuratObj, name) {
  return(!is.null(seuratObj@misc[[paste0(name, 'Run')]]))
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
markStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)
  }

  return(seuratObj)
}




#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
mergeSeuratObjs <- function(seuratObjs, data, alignData = T, MaxCCAspaceDim = 20, MaxPCs2Weight = 20){
  for (exptNum in names(data)) {
    print(paste0('adding expt: ', exptNum))
    prefix <- paste0(exptNum)
    so <- seuratObjs[[exptNum]]

    if (!('BarcodePrefix' %in% names(so@meta.data))) {
      so <- RenameCells(object = so, add.cell.id = prefix)
      so[['BarcodePrefix']] <- c(prefix)
    } else {
      print('Barcode prefix already added')
    }

    if (alignData && length(seuratObjs) > 1) {
      if (!hasStepRun(so, 'NormalizeData')) {
        print('Normalizing')
        so <- NormalizeData(object = so, verbose = F)
      } else {
        print('Normalization performed')
      }

      if (!hasStepRun(so, 'FindVariableFeatures')) {
        print('Finding variable features')
        so <- FindVariableFeatures(object = so, verbose = F, selection.method = "vst", nfeatures = 2000)
      } else {
        print('FindVariableFeatures performed')
      }
    }

    seuratObjs[[exptNum]] <- so
  }

  seuratObj <- NULL
  if (alignData && length(seuratObjs) > 1) {
    # dims here means : Which dimensions to use from the CCA to specify the neighbor search space
    anchors <- FindIntegrationAnchors(object.list = seuratObjs, dims = 1:MaxCCAspaceDim, scale = T, verbose = T)
    # dims here means : #Number of PCs to use in the weighting procedure
    seuratObj <- IntegrateData(anchorset = anchors, dims = 1:MaxPCs2Weight)
    DefaultAssay(seuratObj) <- "integrated"
  }
  else {
    for (exptNum in names(data)) {
      if (is.null(seuratObj)) {
        seuratObj <- seuratObjs[[exptNum]]
      } else {
        seuratObj <- merge(x = seuratObj,
                           y = seuratObjs[[exptNum]],
                           project = outPrefix,
                           do.normalize = F)
      }
    }

    print('after merge')
    print(seuratObj)
  }

  return(seuratObj)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
processSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
                           nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                           nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf, forceReCalc = F,
                           variableGeneTable = NULL, variableFeatureSelectionMethod = 'vst'){

  if (doCellFilter & (forceReCalc | !hasStepRun(seuratObj, 'FilterCells'))) {
    print("Filtering Cells...")
    seuratObj@misc$OriginalCells <- length(colnames(x = seuratObj))
    seuratObj <- subset(x = seuratObj, subset = nCount_RNA > nGene.low & nCount_RNA < nGene.high)
    seuratObj <- subset(x = seuratObj, subset = nFeature_RNA > nUMI.low & nFeature_RNA < nUMI.high)
    seuratObj <- subset(x = seuratObj, subset = p.mito > pMito.low & p.mito < pMito.high)

    print(paste0('Initial cells: ', seuratObj@misc$OriginalCells, ', after filter: ', length(colnames(x = seuratObj))))

    seuratObj <- markStepRun(seuratObj, 'FilterCells')
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'NormalizeData')) {
    seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)
    seuratObj <- markStepRun(seuratObj, 'NormalizeData', saveFile)
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'FindVariableFeatures')) {
    seuratObj <- FindVariableFeatures(object = seuratObj, mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), verbose = F, selection.method = variableFeatureSelectionMethod)
    seuratObj <- markStepRun(seuratObj, 'FindVariableFeatures', saveFile)
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'ScaleData')) {
    seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "percent.mito"), display.progress = F, verbose = F)
    seuratObj <- markStepRun(seuratObj, 'ScaleData')
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'RunPCA')) {
    vg <- VariableFeatures(object = seuratObj)
    print(paste0('Total variable genes: ', length(vg)))
    if (!is.null(variableGeneTable)){
      write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
    }

    seuratObj <- RunPCA(object = seuratObj, features = vg, verbose = F)
    seuratObj <- markStepRun(seuratObj, 'RunPCA')
  }

  if (doCellCycle & (forceReCalc | !hasStepRun(seuratObj, 'CellCycle'))) {
    seuratObj <- removeCellCycle(seuratObj)
    seuratObj <- markStepRun(seuratObj, 'CellCycle', saveFile)
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'ProjectDim')) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- markStepRun(seuratObj, 'ProjectDim')
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
    seuratObj <- markStepRun(seuratObj, 'JackStraw', saveFile)
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- markStepRun(seuratObj, 'ScoreJackStraw')
  }

  print(paste0('Variable genes: ', length(x = VariableFeatures(object = seuratObj))))

  print(VizDimLoadings(object = seuratObj, dims = 1:2))
  print(DimPlot(object = seuratObj))

  print(DimHeatmap(object = seuratObj, dims = 1, cells = 500, balanced = TRUE, fast = F) + NoLegend())
  print(DimHeatmap(object = seuratObj, dims = 1:20, cells = 500, balanced = TRUE, fast = F) + NoLegend())

  print(JackStrawPlot(object = seuratObj, dims = 1:20))
  print(ElbowPlot(object = seuratObj))

  return(seuratObj)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
downloadAndAppendTcrClonotypes <- function(seuratObject, outPath = '.'){
  if (is.null(seuratObject[['BarcodePrefix']])){
    stop('Seurat object lacks BarcodePrefix column')
  }

  for (barcodePrefix in unique(unique(unlist(seuratObject[['BarcodePrefix']])))) {
    print(paste0('Adding TCR clonotypes for prefix: ', barcodePrefix))

    vloupeId <- findMatchedVloupe(barcodePrefix)
    if (is.na(vloupeId)){
      stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
    }

    clonotypeFile <- file.path(outPath, paste0(barcodePrefix, '_clonotypes.csv'))
    downloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = T)
    if (!file.exists(clonotypeFile)){
      stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
    }

    seuratObject <- appendTcrClonotypes(seuratObject, clonotypeFile, barcodePrefix)
  }

  return(seuratObject)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
appendTcrClonotypes <- function(seuratObject = NA, clonotypeFile = NA, barcodePrefix = NULL){
  tcr <- processAndAggregateTcrClonotypes(clonotypeFile)

  if (!is.null(barcodePrefix)){
    tcr$barcode <- as.character(tcr$barcode)
    tcr$barcode <- paste0(barcodePrefix, '_', tcr$barcode)
    tcr$barcode <- as.factor(tcr$barcode)
  }

  origRows <- nrow(tcr)

  datasetSelect <- seuratObject$BarcodePrefix == barcodePrefix
  gexBarcodes <- colnames(seuratObject)[datasetSelect]

  tcr <- tcr[tcr$barcode %in% gexBarcodes,]
  pct <- nrow(tcr) / origRows * 100

  print(paste0('Barcodes with clonotypes: ', origRows, ', intersecting with GEX data: ', nrow(tcr), " (", pct, "%)"))

  merged <- merge(data.frame(barcode = gexBarcodes, sortOrder = 1:length(gexBarcodes)), tcr, by = c('barcode'), all.x = T)
  rownames(merged) <- merged$barcode
  merged <- arrange(merged, sortOrder)
  merged <- merged[colnames(merged) != 'sortOrder']

  # Check barcodes match before merge
  if (sum(merged$barcode != colnames(seuratObject)[datasetSelect]) > 0) {
    stop(paste0('Seurat and TCR barcodes do not match after merge, total different: ', sum(merged$barcode != colnames(seuratObject)[datasetSelect])))
  }

  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    toAdd <- as.character(merged[[colName]])
    names(toAdd) <- merged[['barcode']]

    if (!(colName %in% names(seuratObject@meta.data))) {
      seuratObject[[colName]] <- c(NA)
    }

    toUpdate <- as.character(seuratObject[[colName]])
    toUpdate[datasetSelect] <- toAdd
    seuratObject[[colName]] <- as.factor(toUpdate)
  }

  return(seuratObject)
}

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
findMatchedVloupe <- function(loupeDataId) {
  rows <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
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
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
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


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
downloadCellRangerClonotypes <- function(vLoupeId, outFile, overwrite = T) {
  #There should be a file named all_contig_annotations.csv in the same directory as the VLoupe file
  rows <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu",
    folderPath=paste0("/Labs/Bimber/"),
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


#' FunctionName.
#'
#' @param clonotypeFile, Expects all_contig_annotations.csv from cellranger vdj
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
processAndAggregateTcrClonotypes <- function(clonotypeFile){
  tcr <- read.table(clonotypeFile, header=T, sep = ',')
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

  labelDf$LabelCol <- coalesce(labelDf$displayname, labelDf$clonename)

  labelDf <- labelDf %>%
    group_by(chain, cdr3) %>%
    summarize(CloneName = paste(sort(unique(LabelCol)), collapse = ","))

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

  # Summarize, grouping by barcode
  tcr <- tcr %>% group_by(barcode) %>% summarize(
    ChainCDR3s = paste(sort(unique(ChainCDR3s)), collapse = ","),
    CDR3s = paste(sort(unique(cdr3)), collapse = ","),
    TRA = paste(sort(unique(as.character(TRA))), collapse = ","),
    TRB = paste(sort(unique(as.character(TRB))), collapse = ","),
    TRD = paste(sort(unique(as.character(TRD))), collapse = ","),
    TRG = paste(sort(unique(as.character(TRG))), collapse = ","),
    TRAV = paste(sort(unique(as.character(TRAV))), collapse = ","),
    TRBV = paste(sort(unique(as.character(TRBV))), collapse = ","),
    TRDV = paste(sort(unique(as.character(TRDV))), collapse = ","),
    TRGV = paste(sort(unique(as.character(TRGV))), collapse = ","),

    CloneNames = paste(sort(unique(CloneName)), collapse = ",")  #this is imprecise b/c we count a hit if we match any chain, but this is probably what we often want
  )

  # Note: we should attempt to produce a more specfic call, assuming we have data from multiple chains
  # The intent of this was to allow a A- or B-only hit to produce a call, but if we have both A/B, take their intersect.

  tcr$CloneNames <- sapply(strsplit(as.character(tcr$CloneNames), ",", fixed = TRUE), function(x) paste(naturalsort(unique(x)), collapse = ","))
  tcr$CloneNames[tcr$CloneNames == 'NA'] <- NA

  tcr$barcode <- as.factor(tcr$barcode)
  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    v <- tcr[[colName]]
    v <- as.character(v)
    v[v == ''] <- NA

    tcr[[colName]] <- as.factor(v)
  }

  return(tcr)
}

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
removeCellCycle <- function(seuratObj) {
  print("Performing cell cycle cleaning ...")

  con <- url("https://raw.githubusercontent.com/bbimber/rnaseq/master/data/cellCycle/regev_lab_cell_cycle_genes.txt")
  cc.genes <- readLines(con = con, warn = F)
  close(con)
  if (length(cc.genes) != 97) {
    stop('Something went wrong downloading gene list')
  }

  con <- url("https://raw.githubusercontent.com/bbimber/rnaseq/master/data/cellCycle/G2M.txt")
  g2m.genes <- readLines(con =  con, warn = F)
  close(con)

  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- cc.genes[1:43]
  g2m.genes <- unique(c(g2m.genes, cc.genes[44:97]))

  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]

  if (length(g2m.genes) < 20 || length(s.genes) < 20) {
    print("Warning, the number of g2m and/or s genes in your data has low coverage")
  }

  #proceeds <20 but warns, but <5 is fishy and cant use
  if (length(g2m.genes) < 5 || length(s.genes) < 5) {
    print("Error, the number of g2m and/or s genes < 5")
    #break()
    seuratObj <- markStepRun(seuratObj, 'FAIL_removeCellCycle')
    return(seuratObj) # for pipeline not breaking,... but
  }

  print("Running PCA with cell cycle genes")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE, verbose = F)
  print(DimPlot(object = seuratObj, reduction = "pca"))

  #store values to append later
  SeuratObjsCCPCA <- as.data.frame(seuratObj@reductions$pca@cell.embeddings)
  colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")

  seuratObj <- CellCycleScoring_SERIII(object = seuratObj,
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

  print("Running PCA with variable genes ...")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = VariableFeatures(object = seuratObj), do.print = F, verbose = F)

  for (colName in colnames(SeuratObjsCCPCA)) {
    seuratObj[[colName]] <- SeuratObjsCCPCA[colnames(seuratObj),colName]
  }

  return(seuratObj)
}

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
findClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, saveFile = NULL, forceReCalc = F) {
  if (is.null(dimsToUse)) {
    elbow <- findSeuratElbow(seuratObj)
    print(paste0('Inferred elbow: ', elbow))

    dimsToUse <- 1:elbow
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'FindNeighbors')) {
    seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
    seuratObj <- markStepRun(seuratObj, 'FindNeighbors')
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'FindClusters')) {
    for (resolution in c(0.2, 0.4, 0.8, 1.2, 0.6)){
      seuratObj <- FindClusters(object = seuratObj, resolution = resolution)
      seuratObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = seuratObj, verbose = F)
      seuratObj <- markStepRun(seuratObj, 'FindClusters', saveFile)
    }
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'RunTSNE')) {
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
    seuratObj <- markStepRun(seuratObj, 'RunTSNE', saveFile)
  }

  if (forceReCalc | !hasStepRun(seuratObj, 'RunUMAP')) {
    seuratObj <- RunUMAP(seuratObj,
                         dims = dimsToUse,
                         n.neighbors = 40L,
                         min.dist = 0.2,
                         metric = "correlation",
                         seed.use = 1234)
    seuratObj <- markStepRun(seuratObj, 'RunUMAP', saveFile)
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

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
findMarkers <- function(seuratObj, resolutionToUse, outFile, saveFileMarkers = NULL,
                        testsToUse = c('wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'),
                        numGenesToSave = 20) {

  Idents(seuratObj) <- seuratObj[[paste0('ClusterNames_',resolutionToUse)]]

  if (file.exists(saveFileMarkers)) {
    print('resuming from file')
    seuratObj.markers <- readRDS(saveFileMarkers)
  } else {
    seuratObj.markers <- NA
    for (test in testsToUse) {
      print(paste0('Running using test: ', test))
      tMarkers <- FindAllMarkers(object = seuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F, test.use = test)
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
      print('cluster column not found!')
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

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
createExampleData <- function(nRow = 100, nCol = 10){
  my.counts <- matrix(rpois(1000, lambda=5), ncol=nCol, nrow=nRow)
  my.counts <- as(my.counts, "dgCMatrix")
  cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))

  ngenes <- nrow(my.counts)
  gene.ids <- paste0("ENSG0000", seq_len(ngenes))
  gene.symb <- paste0("GENE", seq_len(ngenes))

  # Writing this to file:
  tmpdir <- tempfile()
  write10xCounts(tmpdir, my.counts, gene.id=gene.ids,
                 gene.symbol=gene.symb, barcodes=cell.ids)
  return(tmpdir)
}

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
CellCycleScoring_SERIII <- function (object, s.features, g2m.features, set.ident = FALSE) {
  enrich.name <- 'Cell Cycle'
  genes.list <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
  object <- AddModuleScore_SERIII(
    object = object,
    genes.list = genes.list,
    enrich.name = enrich.name,
    ctrl.size = min(vapply(X = genes.list, FUN = length, FUN.VALUE = numeric(1)))
  )
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object@meta.data))
  cc.scores <- object@meta.data[, cc.columns]

  gc(verbose = FALSE)
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        return(c(first, second)[which(x = scores == max(scores))])
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]

  object$S.Score <- cc.scores$S.Score
  object$G2M.Score <- cc.scores$G2M.Score
  object$Phase <- cc.scores$Phase

  if (set.ident) {
    object$old.or.idents <- Idents(object = object)
    Idents(object = object) <- cc.scores$Phase
  }
  return(object)
}

#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
AddModuleScore_SERIII <- function(
  #using Seurat v2, AddModuleScore and converting to v3
  #this worked in v2 but their new version 3 breaks so this replaces it

  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster",
  random.seed = 1
) {
  set.seed(seed = random.seed)
  genes.old <- genes.list
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(object)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (!all(Seurat:::LengthCheck(values = genes.list))) {
    warning(paste(
      'Could not find enough genes in the object from the following gene lists:',
      paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
      'Attempting to match case...'
    ))
    genes.list <- lapply(
      X = genes.old,
      FUN = CaseMatch, match = rownames(x = object@data)
    )
  }
  if (!all(Seurat:::LengthCheck(values = genes.list))) {
    stop(paste(
      'The following gene lists do not have enough genes present in the object:',
      paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
      'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(object)
  }
  data.avg <- Matrix::rowMeans(object@assays[[DefaultAssay(object)]]@data[genes.pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object@assays[[DefaultAssay(object)]]@data)
  )

  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object@assays[[DefaultAssay(object)]]@data[genes.use, ])
  }
  genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object@assays[[DefaultAssay(object)]]@data)
  )
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- object@assays[[DefaultAssay(object)]]@data[genes.use, , drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
  rownames(x = genes.scores.use) <- colnames(x = object@assays[[DefaultAssay(object)]]@data)

  for (colName in colnames(genes.scores.use)) {
    object[[colName]] <- genes.scores.use[colnames(object),colName]
  }

  gc(verbose = FALSE)

  return(object)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
findSeuratElbow <- function(seuratObj, ndims = 25, reduction = "pca", print.plot = T){
  data.use <- Stdev(object = seuratObj, reduction = reduction)

  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use),
            " reductions")
    ndims <- length(x = data.use)
  }

  #TODO: what if this cant infer the elbow?  I think it will return a data.frame.
  #1 sd = 1.3
  #the findElbow should not fail, internally if it does it should return 2. But for whatever reason unpredicted it fails, this will capture it.
  elbowX <- try(findElbow(data.use[1:ndims], plot = T, ignore.concavity = F, min.y = 1.3))
  if (class(elbowX)=="try-error" || elbowX[1]==2) {
    if (is.null(ndims)){
      elbowX = 2
    } else {
      elbowX = ndims
    }
  }

  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = "dims", y = "stdev")) +
    labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = seuratObj[[reduction]])),
         y = "Standard Deviation") + theme_bw() + geom_vline(xintercept = elbowX) + ggtitle("Elbow Identification")

  if (print.plot) {
    print(plot)
  }

  return(elbowX)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
findElbow <- function(y, plot = FALSE, ignore.concavity = FALSE, min.y = NA, min.x = NA) {

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
  #' tmp <- findElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
  #' 	plot = TRUE) # wandering
  #' tmp <- findElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
  #' 	plot = TRUE) # late rise
  #' tmp <- findElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
  #' 	plot = TRUE) # gradual, no obvious break
  #'
  #' # Not the usual way to choose the number of PCs:
  #' library("chemometrics")
  #' data(glass)
  #' pca <- prcomp(glass)
  #' eigensum <- sum(pca$sdev * pca$sdev)
  #' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
  #' cs <- cumsum(vv)
  #' tmp <- findElbow(vv, plot = TRUE)
  #' tmp <- findElbow(cs, plot = TRUE)
  #'

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
    print(ix)
    print(iy)
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if (any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if (any((u < 0.00001) || (u > 1))) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      #TODO: giving warning.  needs any() or all()??
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

  ### Now for the actual findElbow function!

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
    print("Your curve doesn't appear to be concave")
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


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
writeSummaryMetrics <- function(seuratObj, file) {
  df <- data.frame(Category = "Seurat", MetricName = "TotalCells", Value = ncol(seuratObj))
  df <- rbind(df, data.frame(Category = "Seurat", MetricName = "TotalFeatures", Value = nrow(seuratObj)))

  write.table(df, file = file, quote = F, row.names = F, sep = '\t')
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
writeCellBarcodes <- function(seuratObj, file) {
  df <- data.frame(CellBarcode = colnames(seuratObj))

  write.table(df, file = file, quote = F, row.names = F, sep = ',', col.names = F)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
saveDimRedux <- function(serObj, reductions=c("pca", "tsne", "umap"),
                         file=NA, maxPCAcomps=10){

  if (is.na(file)){
    stop("file is required")
  }

  tempDT <- data.table(cbind(1:ncol(serObj), colnames(serObj)))
  rownames(tempDT) <- colnames(serObj)
  colnames(tempDT) <- c("nID", "cID")

  if("tsne" %in% reductions) {
    if(is.null(serObj@reductions$tsne)) print("tsne slot NULL") else {

    }
    tsneDT <- data.table(serObj@reductions$tsne@cell.embeddings)
    tsneDT$cID <- rownames(serObj@reductions$tsne@cell.embeddings)
    tempDT <- merge(tempDT, tsneDT, by="cID")
  }

  if("pca" %in% reductions) {
    if(is.null(serObj@reductions$pca)) print("pca slot NULL") else {
      pcaDT <- data.table(serObj@reductions$pca@cell.embeddings[,1:maxPCAcomps])
      pcaDT$cID <- rownames(serObj@reductions$pca@cell.embeddings)
      tempDT <- merge(tempDT, pcaDT, by="cID")
    }

  }

  if("umap" %in% reductions) {
    if(is.null(serObj@reductions$umap)) print("umap slot NULL") else {
      umapDT <- data.table(serObj@reductions$umap@cell.embeddings)
      umapDT$cID <- rownames(serObj@reductions$umap@cell.embeddings)
      tempDT <- merge(tempDT, umapDT, by="cID")
    }

  }

  print("saving DimRedux")

  write.csv(tempDT, file = file, row.names=TRUE)
}


#' FunctionName.
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)
addTitleToMultiPlot <- function(plotGrid, title, relHeights = c(0.1, 1)) {
  return(plot_grid(ggdraw() + draw_label(title), plotGrid, ncol = 1, rel_heights = relHeights))
}
