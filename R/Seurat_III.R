#' @include LabKeySettings.R
#' @include Seurat_III_Fixes.R
#' @import Seurat
#' @import Rlabkey

utils::globalVariables(
  names = c('nCount_RNA', 'nFeature_RNA', 'p.mito'),
  package = 'OOSAP',
  add = TRUE
)


#' @title Read and Filter 10X files.
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
#' @param dataDir The directory holding raw count data
#' @param datasetName A name to use when creating the Seurat object
#' @param emptyDropNIters The number of iterations to use with PerformEmptyDrops()
#' @return A Seurat object.
#' @keywords ReadAndFilter10X
#' @export
#' @importFrom Seurat Read10X
ReadAndFilter10xData <- function(dataDir, datasetName, emptyDropNIters=10000) {
  if (!file.exists(dataDir)){
    stop(paste0("File does not exist: ", dataDir))
  }

  if (!dir.exists(dataDir)){
    stop(paste0("File is not a directory: ", dataDir))
  }

  seuratRawData <- Read10X(data.dir = dataDir)

  #Cannot have underscores in feature names, Seurat will replace with hyphen anyway.  Perform upfront to avoid warning
  if (sum(grepl(x = rownames(seuratRawData), pattern = '_')) > 0) {
    print('Replacing underscores with hyphens in feature names')
    rownames(seuratRawData) <- gsub(x = rownames(seuratRawData), pattern = '_', replacement = '-')
  }

  seuratRawData <- PerformEmptyDropletFiltering(seuratRawData, emptyDropNIters=emptyDropNIters)

  seuratObj <- CreateSeuratObj(seuratRawData, project = datasetName)
  PrintQcPlots(seuratObj)

  return(seuratObj)
}


#' @title Create a Seurat 3 object
#'
#' @description Create Seurat object from count data (usually from Read10X()). This also sets pct.mito.
#' @param seuratData, Seurat input data, usually from Read10X().
#' @param project, Sets the project name for the Seurat object.
#' @param minFeatures, Include cells where at least this many features are detected.
#' @param minCells, Include features detected in at least this many cells.
#' @param mitoGenesPattern The expression to use when identfying mitochondrial genes
#' @return A Seurat object with p.mito calculated.
#' @keywords CreateSeuratObj
#' @importFrom Matrix colSums
CreateSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25, minCells = 0, mitoGenesPattern = "^MT-"){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = minCells, min.features = minFeatures, project = project)

  mito.features <- grep(pattern = mitoGenesPattern, x = rownames(x = seuratObj), value = TRUE)
  p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['p.mito']] <- p.mito

  return(seuratObj)
}


#' @title PrintQcPlots
#'
#' @description Prints a number of QC plots from a seurat object, generally following the basic Seurat vignette
#' @param seuratObj A Seurat object.
#' @return NULL
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


#' @title PerformEmptyDropletFiltering
#'
#' @param seuratRawData Raw data
#' @param fdrThreshold FDR threshold, passed directly to PerformEmptyDrops()
#' @param emptyDropNIters Number of iterations, passed directly to PerformEmptyDrops()
#' @return Plot
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
#' @param seuratRawData Raw data
#' @param fdrThreshold FDR threshold, passed directly to PerformEmptyDrops()
#' @param emptyDropNIters Number of iterations, passed directly to PerformEmptyDrops()
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


#' @title HasStepRun
#' @description An internal method to mark steps as complete, to avoid recalculating
#' @param seuratObj The seurat object
#' @param name The name of the step to mark complete
#' @return A modified Seurat object.
HasStepRun <- function(seuratObj, name, forceReCalc = F, printResult = T) {
  ret <- !is.null(seuratObj@misc[[paste0(name, 'Run')]])
  if (ret && printResult) {
    if (forceReCalc) {
      print(paste0('Step already run, but forceReCalc=TRUE, so will repeat: ', name))
    } else {
      print(paste0('Step already run, will not repeat: ', name))
    }
  }

  return(ret)
}


#' @title MarkStepRun
#' @description An internal method to determine is a step has been marked as complete
#' @param seuratObj The seurat object
#' @param name The name of the step to test
#' @param saveFile If provided, the seurat object will be saved as RDS to this location
# @return A modified Seurat object.
MarkStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)
  }

  return(seuratObj)
}


#' @title MergeSeuratObjs
#' @description Merges a list of Seurat objects, using Seurat::IntegrateData()
#' @param seuratObjs A list of seurat objects, optionally named (in which case these will be used as dataset names)
#' @param metadata A list of metadata.  If provided, the names of this list will be used as dataset names
#' @param alignData If true, data will be aligned using Seurat::IntegrateData()
#' @param maxCCAspaceDim The number of dims to use with FindIntegrationAnchors()
#' @param maxPCs2Weight The number of dims to use with IntegrateData()
#' @param projectName The project name when creating the final seuratObj
#' @param doScaleEach If true, scale=T will be passed to FindIntegrationAnchors(); its default behavior to find anchors. 
#' @param useAllFeatures If true, the resulting object will contain all features, as opposed to just VariableGenes (not recommended)
#' @param nVariableFeatures The number of VariableFeatures to identify
#' @param includeCellCycleGenes If true, the cell cycles genes will always be included with IntegrateData(), as opposed to just VariableGenes
#' @return A modified Seurat object.
#' @export
#' @importFrom methods slot
MergeSeuratObjs <- function(seuratObjs, metadata=NULL, alignData = T, maxCCAspaceDim = 20, maxPCs2Weight = 20,
projectName = NULL, doScaleEach = T, useAllFeatures = F, nVariableFeatures = 2000,
includeCellCycleGenes = T){
  nameList <- NULL
  if (is.null(metadata)){
    nameList <- names(seuratObjs)
  } else {
    nameList <- names(metadata)
  }

  for (exptNum in nameList) {
    print(paste0('adding dataset: ', exptNum))
    prefix <- paste0(exptNum)
    so <- seuratObjs[[exptNum]]

    if (!('BarcodePrefix' %in% names(so@meta.data))) {
      print(paste0('Adding barcode prefix: ', prefix))
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
      
      if (HasStepRun(so, 'ScaleData') & doScaleEach) {
         warning("doScaleEach = T and this object is scaled; adding to time complexity")
      } else {
        print('ScaleData not prev. performed')
      }
      
      

      print(LabelPoints(plot = VariableFeaturePlot(so), points = head(VariableFeatures(so), 20), repel = TRUE, xnudge = 0, ynudge = 0))
    }

    seuratObjs[[exptNum]] <- so
  }

  seuratObj <- NULL
  if (alignData && length(seuratObjs) > 1) {
    CheckDuplicatedCellNames(seuratObjs)

    # dims here means : Which dimensions to use from the CCA to specify the neighbor search space
    anchors <- FindIntegrationAnchors(object.list = seuratObjs, dims = 1:maxCCAspaceDim, scale = doScaleEach, verbose = F)

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
    seuratObj <- IntegrateData(anchorset = anchors, dims = 1:maxPCs2Weight, verbose = F, features.to.integrate = features, new.assay.name = "Integrated")
    Seurat::DefaultAssay(seuratObj) <- "Integrated"

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


#' @title CheckDuplicatedCellNames
#'
#' @description A variant on Seurat CheckDuplicatedCellNames that reports the top duplicated barcodes.  Mainly for debugging
#' @param object.list A lists of Seurat objects
#' @param stop Boolean that determines if stop() or print() is used when duplicates are found
#' @return NULL
CheckDuplicatedCellNames <- function(object.list, stop = TRUE){
  cell.names <- unlist(
  x = sapply(
  X = 1:length(x = object.list),
  FUN = function(x) Cells(object.list[[x]])
  )
  )

  dups <- duplicated(x = cell.names)
  if (any(dups)) {
    if (stop){
      stop(paste0('There were duplicated cell names: ',  paste0(head(unique(cell.names[dups])), collapse = ',')))
    } else {
      print('There were duplicated cell names')
      print(head(unique(cell.names[dups])))
    }
  }
}


#' @title Run the primary seurat processing steps.
#'
#' @description This is the primary entry point for processing scRNAseq data with Seurat
#' @param seuratObj, A Seurat object.
#' @param saveFile If provided, the seuratObj will be saved here as it is processed, providing some ability to resume if there is a failure
#' @param doCellCycle If true, CellCycle genes will be regressed
#' @param doCellFilter If true, basic filtering will be performed using nCount_RNA, nFeature_RNA, and pMito
#' @param uUMI.high If doCellFilter=T, cells with nUMI above this value will be filtered
#' @param uUMI.low If doCellFilter=T, cells with nUMI below this value will be filtered
#' @param nGene.high If doCellFilter=T, cells with nGene above this value will be filtered
#' @param nGene.low If doCellFilter=T, cells with nGene below this value will be filtered
#' @param pMito.high If doCellFilter=T, cells with percent mito above this value will be filtered
#' @param pMito.low If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param forceReCalc If true, all steps will be repeated even if already marked as complete
#' @param variableGeneTable If provided, a table of variable genes will be written to this file
#' @param variableFeatureSelectionMethod The selection method to be passed to FindVariableFeatures()
#' @param useSCTransform If true, SCTransform will be used in place of the standard Seurat workflow (NormalizeData, ScaleData, FindVariableFeatures)
#' @param nVariableFeatures The number of variable features to find
#' @param printDefaultPlots If true, the default set of QC plots will be printed
#' @param npcs Number of PCs to use for RunPCA()
#' @param ccPcaResultFile If provided, the PCA results from cell cycle regression will be written here
#' @return A modified Seurat object.
#' @export
ProcessSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
                            nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                            nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf, forceReCalc = F,
                            variableGeneTable = NULL, variableFeatureSelectionMethod = 'vst', 
                            nVariableFeatures = 2000, printDefaultPlots = T,
                            npcs = 50, ccPcaResultFile = NULL, useSCTransform = F, 
                            mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), 
                            spikeGenes = NULL){

  if (doCellFilter & (forceReCalc | !HasStepRun(seuratObj, 'FilterCells', forceReCalc = forceReCalc))) {
    print("Filtering Cells...")
    seuratObj@misc$OriginalCells <- length(colnames(x = seuratObj))

    #See: https://github.com/satijalab/seurat/issues/1053#issuecomment-454512002
    expr <- Seurat::FetchData(object = seuratObj, vars = 'nCount_RNA')
    seuratObj <- seuratObj[, which(x = expr > nGene.low & expr < nGene.high)]

    expr <- Seurat::FetchData(object = seuratObj, vars = 'nFeature_RNA')
    seuratObj <- seuratObj[, which(x = expr > nUMI.low & expr < nUMI.high)]

    expr <- Seurat::FetchData(object = seuratObj, vars = 'p.mito')
    seuratObj <- seuratObj[, which(x = expr > pMito.low & expr < pMito.high)]

    print(paste0('Initial cells: ', seuratObj@misc$OriginalCells, ', after filter: ', length(colnames(x = seuratObj))))

    seuratObj <- MarkStepRun(seuratObj, 'FilterCells')
  }

  if (!useSCTransform) {
    if (forceReCalc | !HasStepRun(seuratObj, 'NormalizeData', forceReCalc = forceReCalc)) {
      seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'NormalizeData', saveFile)
    }

    if (forceReCalc | !HasStepRun(seuratObj, 'FindVariableFeatures', forceReCalc = forceReCalc)) {
      seuratObj <- FindVariableFeatures(object = seuratObj, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff , verbose = F, selection.method = variableFeatureSelectionMethod, nfeatures = nVariableFeatures)
      seuratObj <- MarkStepRun(seuratObj, 'FindVariableFeatures', saveFile)
    }

    if (forceReCalc | !HasStepRun(seuratObj, 'ScaleData', forceReCalc = forceReCalc)) {
      seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "p.mito"), verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'ScaleData')
    }

    if (doCellCycle & (forceReCalc | !HasStepRun(seuratObj, 'CellCycle', forceReCalc = forceReCalc))) {
      seuratObj <- RemoveCellCycle(seuratObj, pcaResultFile = ccPcaResultFile, useSCTransform = F)
      seuratObj <- MarkStepRun(seuratObj, 'CellCycle', saveFile)
    }
  } else {
    print('Using SCTransform')
    seuratObj <- SCTransform(seuratObj, vars.to.regress = c("nCount_RNA", "p.mito"), verbose = FALSE, return.only.var.genes = F)
    
    if (doCellCycle & (forceReCalc | !HasStepRun(seuratObj, 'CellCycle', forceReCalc = forceReCalc))) {
      seuratObj <- RemoveCellCycle(seuratObj, pcaResultFile = ccPcaResultFile, useSCTransform = T)
      seuratObj <- MarkStepRun(seuratObj, 'CellCycle', saveFile)
    }
  }
  
  if(!is.null(spikeGenes)){
    VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), spikeGenes))
  }

  vg <- VariableFeatures(object = seuratObj)
  
  if (forceReCalc | !HasStepRun(seuratObj, 'RunPCA', forceReCalc = forceReCalc)) {
    seuratObj <- RunPCA(object = seuratObj, features = vg, verbose = F, npcs = npcs)
    seuratObj <- MarkStepRun(seuratObj, 'RunPCA')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ProjectDim', forceReCalc = forceReCalc)) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- MarkStepRun(seuratObj, 'ProjectDim')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'JackStraw', forceReCalc = forceReCalc)) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
    seuratObj <- MarkStepRun(seuratObj, 'JackStraw', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ScoreJackStraw', forceReCalc = forceReCalc)) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- MarkStepRun(seuratObj, 'ScoreJackStraw')
  }

  print(paste0('Total variable genes: ', length(vg)))
  if (!is.null(variableGeneTable)){
    write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
  }

  if (printDefaultPlots){
    print(VizDimLoadings(object = seuratObj, dims = 1:2))
    print(LabelPoints(plot = VariableFeaturePlot(seuratObj), points = head(VariableFeatures(seuratObj), 20), repel = TRUE, xnudge = 0, ynudge = 0))

    print(DimPlot(object = seuratObj))
    if (doCellCycle) {
      print(cowplot::plot_grid(plotlist = list(
        DimPlot(object = seuratObj, reduction = "pca", dims = c(1, 2), group.by = 'Phase'),
        DimPlot(object = seuratObj, reduction = "pca", dims = c(2, 3), group.by = 'Phase'),
        DimPlot(object = seuratObj, reduction = "pca", dims = c(3, 4), group.by = 'Phase'),
        DimPlot(object = seuratObj, reduction = "pca", dims = c(4, 5), group.by = 'Phase')
      )))
    }

    print(DimHeatmap(object = seuratObj, dims = 1, cells = 500, balanced = TRUE, fast = F) + NoLegend())
    print(DimHeatmap(object = seuratObj, dims = 1:20, cells = 500, balanced = TRUE, fast = F) + NoLegend())

    print(JackStrawPlot(object = seuratObj, dims = 1:20))
    print(ElbowPlot(object = seuratObj))
  }

  return(seuratObj)
}


#' @title RemoveCellCycle
#' @param seuratObj The seurat object
#' @param pcaResultFile If provided, cell cycle PCA results will be written here
#' @param useSCTransform If true, SCTransform will be used in place of the standard Seurat workflow (NormalizeData, ScaleData, FindVariableFeatures)
#' @param do.scale If true, scale residuals to have unit variance; SCTransform default = F; ScaldData default = T
#' @param do.center If true,center residuals to have mean zero; SCTransform default = T; ScaldData default = T
#' @return A modified Seurat object.
#' @importFrom cowplot plot_grid
RemoveCellCycle <- function(seuratObj, pcaResultFile = NULL, 
                            useSCTransform = F, do.scale = T, do.center = T) {
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
  seuratObj <- RunPCA(object = seuratObj, features = c(s.genes, g2m.genes), do.print = FALSE, verbose = F)
  print(DimPlot(object = seuratObj, reduction = "pca"))

  #store values to append later
  SeuratObjsCCPCA <- as.data.frame(seuratObj@reductions$pca@cell.embeddings)
  colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")

  seuratObj <- CellCycleScoring(object = seuratObj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

  print(cowplot::plot_grid(plotlist = list(
    DimPlot(object = seuratObj, reduction = "pca", dims = c(1, 2)),
    DimPlot(object = seuratObj, reduction = "pca", dims = c(2, 3)),
    DimPlot(object = seuratObj, reduction = "pca", dims = c(3, 4)),
    DimPlot(object = seuratObj, reduction = "pca", dims = c(4, 5))
  )))

  print(table(seuratObj$Phase))

  print("Regressing out S and G2M score ...")

  if (!useSCTransform) {
    seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), 
                           verbose = F, features = rownames(x = seuratObj),
                           do.scale = do.scale, do.center = do.center)
  } else {
    seuratObj <- SCTransform(seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE, 
                             return.only.var.genes = F, do.scale = do.scale, do.center = do.center)
  }

  if (!is.null(pcaResultFile)) {
    SeuratObjsCCPCA$CellBarcode <- colnames(seuratObj)
    write.table(SeuratObjsCCPCA, file = pcaResultFile, sep = '\t', row.names = F, quote = F)
  }
  
  return(seuratObj)
}


#' @title FindClustersAndDimRedux
#' @param seuratObj, A Seurat object.
#' @param dimsToUse The number of dims to use.  If null, this will be inferred using FindSeuratElbow()
#' @param minDimsToUse The minimum numbers of dims to use.  If dimsToUse is provided, this will override.
#' @param saveFile If provided, the seurat object will be saved as RDS to this location
#' @param forceReCalc If true, all steps will be performed even if already marked complete
#' @param umap.method The UMAP method, either uwot or umap-learn, passed directly to Seurat::RunUMAP
#' @return A modified Seurat object.
#' @export
FindClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, saveFile = NULL, forceReCalc = F, minDimsToUse = NULL, umap.method = 'umap-learn',
                                   UMAP_NumNeib = 40L, UMAP_MinDist = 0.2, UMAP_Seed = 1234, UMAP.NumEpoc = 500) {
  if (is.null(dimsToUse)) {
    dimMax <- FindSeuratElbow(seuratObj)
    print(paste0('Inferred elbow: ', dimMax))

    dimsToUse <- 1:dimMax

    if (!is.null(minDimsToUse)) {
      print(paste0('Min dims to use: ', minDimsToUse))
      dimMax <- max(minDimsToUse, dimMax)
    }

    print(paste0('Selected dimsToUse: 1:', dimMax))
    dimsToUse <- 1:dimMax
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'FindNeighbors', forceReCalc = forceReCalc)) {
    seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
    seuratObj <- MarkStepRun(seuratObj, 'FindNeighbors')
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'FindClusters', forceReCalc = forceReCalc)) {
    for (resolution in c(0.2, 0.4, 0.8, 1.2, 0.6)){
      seuratObj <- FindClusters(object = seuratObj, resolution = resolution)
      seuratObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = seuratObj, verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'FindClusters', saveFile)
    }
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'RunTSNE', forceReCalc = forceReCalc)) {
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

  if (forceReCalc | !HasStepRun(seuratObj, 'RunUMAP', forceReCalc = forceReCalc)) {
    seuratObj <- RunUMAP(seuratObj,
                           dims = dimsToUse,
                           n.neighbors = UMAP_NumNeib,
                           min.dist = UMAP_MinDist,
                           metric = "correlation",
                           umap.method = umap.method,
                           seed.use = UMAP_Seed, n.epochs = UMAP.NumEpoc)
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


#' @title Find_Markers
#' @param seuratObj A seurat object
#' @param resolutionToUse The resolution, generally computed during FindClustersAndDimRedux()
#' @param outFile A file where a table of markers will be saved
#' @param saveFileMarkers A file where the dataframe of markers will be saved as RDS
#' @param testsToUse A vector of tests to be used.  Each will be used to run FindAllMarkers() and the results merged
#' @param numGenesToSave The number of top markers per cluster to save
#' @param onlyPos If true, only positive markers will be saved
#' @importFrom dplyr %>% coalesce group_by summarise filter top_n
#' @export
Find_Markers <- function(seuratObj, resolutionToUse, outFile, saveFileMarkers = NULL,
testsToUse = c('wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'),
numGenesToSave = 20, onlyPos = F) {

  Idents(seuratObj) <- seuratObj[[paste0('ClusterNames_',resolutionToUse)]]

  if (!is.null(saveFileMarkers) && file.exists(saveFileMarkers)) {
    print('resuming from file')
    seuratObj.markers <- readRDS(saveFileMarkers)
  } else {
    seuratObj.markers <- NA
    tMarkers <- NA
    for (test in testsToUse) {
      print(paste0('Running using test: ', test))
      tryCatch({
        tMarkers <- FindAllMarkers(object = seuratObj, only.pos = onlyPos, min.pct = 0.25, logfc.threshold = 0.25, verbose = F, test.use = test)
        if (nrow(tMarkers) == 0) {
          print('No genes returned, skipping')
        } else {
          tMarkers$test <- c(test)
          tMarkers$cluster <- as.character(tMarkers$cluster)
          if (test == 'roc') {
            toBind <- data.frame(
              test = tMarkers$test,
              cluster = tMarkers$cluster,
              gene = tMarkers$gene,
              pct.1 = tMarkers$pct.1,
              pct.2 = tMarkers$pct.2,
              avg_logFC = NA,
              p_val_adj = NA,
              myAUC = tMarkers$myAUC,
              power = tMarkers$power,
              avg_diff = tMarkers$avg_diff
            )
          } else {
            toBind <- data.frame(
              test = tMarkers$test,
              cluster = tMarkers$cluster,
              gene = tMarkers$gene,
              pct.1 = tMarkers$pct.1,
              pct.2 = tMarkers$pct.2,
              avg_logFC = tMarkers$avg_logFC,
              p_val_adj = tMarkers$p_val_adj,
              myAUC = NA,
              power = NA,
              avg_diff = NA
            )
          }

          if (all(is.na(seuratObj.markers))) {
            seuratObj.markers <- toBind
          } else {
            seuratObj.markers <- rbind(seuratObj.markers, toBind)
          }
        }
      }, error = function(e){
        print(paste0('Error running test: ', test))
        print(e)
        print(str(tMarkers))
        print(str(seuratObj.markers))
      })
    }

    if (is.na(seuratObj.markers)) {
      print('All tests failed, no markers returned')
      return()
    }

    if (!('cluster' %in% names(seuratObj.markers))) {
      warning('cluster column not found!')
    } else {
      seuratObj.markers$cluster <- as.factor(seuratObj.markers$cluster)
    }

    if (!is.null(saveFileMarkers)){
      saveRDS(seuratObj.markers, file = saveFileMarkers)
    }
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
#' @param object A seurat object
#' @param ndims The number of dims, passed to FindElbow()
#' @param reduction The reduction, such as 'pca'
#' @param print.plot If true, summary plots will be printed
#' @param min.y The min.y, passed directly to FindElbow()
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
    if (any(lineMag < 0.00000001)) {
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)

    ans <- c()
    for (i in 1:length(u)) {
      if (u[i] < 0.00001 || u[i] > 1) {
        print('closest point does not fall within the line segment, take the shorter distance to an endpoint')
        ix <- lineMagnitude(px[i], py[i], x1[i], y1[i])
        iy <- lineMagnitude(px[i], py[i], x2[i], y2[i])
        ans <- c(ans, min(ix, iy))
      } else {
        ## Intersecting point is on the line, use the formula
        ix <- x1 + u * (x2 - x1)
        iy <- y1 + u * (y2 - y1)
        ans <- c(ans, lineMagnitude(px, py, ix, iy)[i])
      }
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
    print('Max not found, defaulting to 2')
    return(2)
  } else {
    return(which.max(DF$dist))
  }
}


#' @title WriteSummaryMetrics
#' @export
#' @param seuratObj, A Seurat object.
#' @param file The file where metrics will be written
WriteSummaryMetrics <- function(seuratObj, file) {
  df <- data.frame(Category = "Seurat", MetricName = "TotalCells", Value = ncol(seuratObj))
  df <- rbind(df, data.frame(Category = "Seurat", MetricName = "TotalFeatures", Value = nrow(seuratObj)))

  write.table(df, file = file, quote = F, row.names = F, sep = '\t')
}


#' @title WriteCellBarcodes
#' @description Writes a table of cell barcodes to the provided file
#' @return A modified Seurat object.
#' @param seuratObj The seurat object
#' @param The file to which barcodes will be written
#' @export
WriteCellBarcodes <- function(seuratObj, file) {
  df <- data.frame(CellBarcode = colnames(seuratObj))

  write.table(df, file = file, quote = F, row.names = F, sep = ',', col.names = F)
}


#' @title SaveDimRedux
#' @return A modified Seurat object.
#' @param seuratObj The seurat object
#' @param reductions Reductions to save
#' @param file The output file
#' @param maxPCAcomps The maximum number of PCs to save
#' @param returnResults Whether to return results
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


#' @title AddTitleToMultiPlot
#'
#' @description Add overall title to a plot_grid, as as the results generated by FeaturePlot() from multiple genes
#' @param plotGrid The list of plots
#' @param title The title for this plot
#' @param relHeights The relative heights, passed to plot_grid
#' @export
AddTitleToMultiPlot <- function(plotGrid, title, relHeights = c(0.1, 1)) {
  return(cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title), plotGrid, ncol = 1, rel_heights = relHeights))
}


utils::globalVariables(
names = c('UMAP1', 'UMAP2'),
package = 'OOSAP',
add = TRUE
)


#' @title ggUMAP
#' @description Generated a ggplot object from the UMAP data of a seurat object
#' @param object A seurat object
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


#' @title GetXYDataFromPlot
#' @description Get XY data from a seurat plot
#' @param plot The plot object
#' @param cellNames The set of cells to export
#' @export
GetXYDataFromPlot <- function(plot, cellNames) {
  xynames <- Seurat:::GetXYAesthetics(plot = plot)

  plot.data <- plot$data[cellNames, ]
  names(plot.data)[names(plot.data) == xynames$x] <- 'x'
  names(plot.data)[names(plot.data) == xynames$y] <- 'y'

  return(plot.data)
}


#' @title AddClonesToPlot
#' @description Can be used to highlight a set of cells from a seurat plot, such as overlaying specific clonotypes
#' @param seuratObj The seurat object
#' @param The plot object, such as the result from FeaturePlot()
#' @param colorShapes A string or vector that is passed to geom_point()
#' @export
#' @import ggplot2
AddClonesToPlot <- function(seuratObj, plot, colorShapes = "black") {
  cellNames <- colnames(seuratObj)[!is.na(seuratObj$CloneNames)]
  plot.data <- GetXYDataFromPlot(plot, cellNames)
  plot.data$CloneName <- seuratObj$CloneNames[!is.na(seuratObj$CloneNames)]

  plot <- plot + geom_point(
  mapping = aes_string(x = 'x', y = 'y', shape = 'CloneName'),
  color = colorShapes,
  data = plot.data
  )

  return(plot)
}


#' @title FilterCloneNames
#' @description Filter Clone Names
#' @param seuratObj The seurat object
#' @param minValue Filters clones not present in at least this many cells
#' @export
FilterCloneNames <- function(seuratObj, minValue) {
  ct <- table(seuratObj$CloneNames)
  ct <- ct[ct < minValue]

  seuratObj$CloneNames[seuratObj$CloneNames %in% names(ct)] <- NA

  return(seuratObj)
}


#' @title AvgCellExprs
#'
#' @param seuratObj, A Seurat object.
#' @param varName The resolution/ident to use
#' @param genes A vector of genes to include
#' @param slot The slot to use, passed to GetAssayData
#' @return A data.frame of avg expression per var
#' @importFrom Matrix rowMeans
#' @export
AvgCellExprs <- function(seuratObj, varName = "ClusterNames_0.2", genes, slot = "scale.data"){
  #Slot : Specific information to pull (i.e. counts, data, scale.data, ...)

  AvlLevels <- factor(as.character(FetchData(seuratObj, varName)[,1]))

  ClustLS <- list()

  for(lev in levels(AvlLevels)){
    print(lev)
    ClustLS[[lev]] <- as.data.frame(Matrix::rowMeans(GetAssayData(object = seuratObj, features = genes, slot = slot)[genes, colnames(seuratObj)[which(AvlLevels==lev)]  ]))
  }

  ClustDF <- as.data.frame(ClustLS)
  colnames(ClustDF) <- paste0("clus", levels(AvlLevels))

  return(ClustDF)
}


#' @title AddModuleScoreAvg
#'
#' @description Instead of the status quo score of Seurat which is 1 score 1 gene, this function takes a list of genes and computes per set, the average of the individual scores.
#' @param object, A Seurat object.
#' @param genes.list, Gene list to obtain a score for
#' @param genes.pool, Gene list to base as the pool; NULL = all.
#' @param n.bin, number of bins to evaluate score across; default 25.
#' @param seed.use, seed use.
#' @param ctrl.size, control gene set size.
#' @param enrich.name, A name for the assesment
#' @param random.seed, random seed
#' @return A modified Seurat object.
#' @keywords Seurat, average, gene
#' @export
#' @importFrom Hmisc cut2
#' @importFrom Matrix colMeans rowMeans
AddModuleScoreAvg <- function(
#May-2019 version

#this is a modified version of the AddModuleScore
#returnScore = F/T controls the output.
#if T, just the scores are returned,
#if F, the scores are put in the Seurat obj and the Seurat SeurObj is returned.
#Also this FX is modified to work for Seurat V3
SeurObj,
genes.list = NULL,
genes.pool = NULL,
n.bin = 25,
seed.use = 1,
ctrl.size = 100,
enrich.name = "Cluster",
random.seed = 1, returnScore = F) {
  set.seed(seed = random.seed)
  genes.old <- genes.list


  if (is.null(x = genes.list)) {
    stop("Missing input gene list")
  }

  genes.list <- lapply(
    X = genes.list,
    FUN = function(x) {
      return(intersect(x = x, y = rownames(SeurObj)))
    }
  )

  cluster.length <- length(x = genes.list)

  if (!all(Seurat:::LengthCheck(values = genes.list))) {
    warning(paste(
    'Could not find enough genes in the SeurObj from the following gene lists:',
    paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
    'Attempting to match case...'
    ))

    genes.list <- lapply(
    X = genes.old,
    FUN = CaseMatch, match = rownames(SeurObj)
    )
  }

  if (!all(Seurat:::LengthCheck(values = genes.list))) {
    stop(paste(
    'The following gene lists do not have enough genes present in the SeurObj:',
    paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
    'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(SeurObj)
  }
  data.avg <- Matrix::rowMeans(Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj))[genes.pool, ])
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
    ncol = ncol(x = Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj)))
  )

  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj))[genes.use, ])
  }
  genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj)))
  )

  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj))[genes.use, , drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
  rownames(x = genes.scores.use) <- colnames(x = Seurat::GetAssayData(SeurObj, assay = Seurat::DefaultAssay(SeurObj)))

  for (colName in colnames(genes.scores.use)) {
    SeurObj[[colName]] <- genes.scores.use[colnames(SeurObj), colName]
  }

  gc(verbose = FALSE)

  if(!returnScore){
    return(SeurObj)
  } else {
    SeurObj@meta.data$cID <- rownames(SeurObj@meta.data)
    return(SeurObj@meta.data[, c("cID", colnames(genes.scores.use))] )
  }
}


#' @title PlotImmuneMarkers
#'
#' @description Generate a set of Seurat FeaturePlots for common immune cell markers
#' @param seuratObj A seurat object
#' @param reduction The reduction to use
#' @export
#' @import Seurat
PlotImmuneMarkers <- function(seuratObj, reduction = 'tsne') {
  PlotMarkerSet(seuratObj, reduction, 'CD8/CD4 Markers', c('CD8A', 'CD8B', 'CD4', 'IL7R'))

  #Eff v. Mem:
  #IL7R = CD127
  #IL2RA = CD25
  #PTPRC = CD45
  #SELL = CD62-L / CD-197
  PlotMarkerSeries(seuratObj, reduction, c('CCR7', 'SELL', 'GZMB', 'CCR5', 'IL2RA', 'PTPRC', 'IL7R', 'CTLA4'), 'Effector vs. Memory')

  #CD8 Activation
  PlotMarkerSeries(seuratObj, reduction, c('CCL4', 'IFNG', 'CD69', 'TNF', 'NFKBID', 'LTB', 'TNFRSF9', 'CCL4L2'), 'CD8 Activation Markers')

  PlotMarkerSeries(seuratObj, reduction, c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'), 'Cytotoxicity')

  PlotMarkerSet(seuratObj, reduction, 'B-cell Markers', c('MS4A1', 'CD79A', 'CD74', 'DRA'))

  PlotMarkerSet(seuratObj, reduction, 'Monocyte', c('LYZ', 'CST3', 'S100A6', 'VIM'))

  PlotMarkerSet(seuratObj, reduction, 'Transcription Factors', c('TBX21', 'GATA3', 'RORC', 'FOXP3'))

  #LILR/KIR:
  PlotMarkerSeries(seuratObj, reduction, c('LILRA5','LILRA6','LILRB4','LILRB5','KIR2DL4','KIR3DX1', 'MAMU-KIR'), 'KIR/LILR')

  PlotMarkerSeries(seuratObj, reduction, c('FCGR1A','FCGR2B','FCGR3'), 'FCGR')

  #Cytokines
  cytokines <- c('IL1A','IL1B','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL2','IL2RA','IL2RB','IL2RG','IL3','IL3RA','IL4','IL4I1','IL4R','IL5','IL5RA','IL6','IL6R','IL6ST','IL7','IL7R','IL9','IL10','IL10RA','IL11','IL12A','IL12B','IL12RB1','IL12RB2','IL13','IL13RA2','IL15','IL15Ra','IL16','IL17A','IL17B','IL17C','IL17D','IL17F','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL18BP','IL18R1','IL18RAP','IL19','IL20','IL20RA','IL20RB','IL21','IL21R','IL22','IL22RA2','IL23A','IL24','IL25','IL26','IL27','IL27RA','IL31','IL31RA','IL33','IL34','IL36A','IL36B','IL36G','IL37','ILDR1','ILDR2','ILF2','ILF3','ILK','ILKAP','ILVBL')
  PlotMarkerSeries(seuratObj, reduction, cytokines, 'Cytokines/Receptors')

  klrs <- c('KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRF2', 'KLRG1', 'KLRG2')
  PlotMarkerSeries(seuratObj, reduction, klrs, 'KLRs')

  #chemokines
  chemokines <- c('CCL1','CCL2','CCL4','CCL4L2','CCL5','CCL8','CCL11','CCL13','CCL14','CCL16','CCL17','CCL18','CCL19','CCL20','CCL21','CCL22','CCL23','CCL24','CCL25','CCL26','CCL27','CCL28')
  chemokines <- c(chemokines, c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CCRL2'))
  chemokines <- c(chemokines, c('CXCL1','CXCL3','CXCL9','CXCL10','CXCL13','CXCL14','CXCL16','CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','XCR1'))

  PlotMarkerSeries(seuratObj, reduction, chemokines, 'Chemokines/Receptors')
}


#' @title PlotMarkerSeries
#'
#' @description Generate a set of Seurat FeaturePlots for the provided markers
#' @param seuratObj The seurat object
#' @param reduction The reduction to use
#' @param features The features to plot
#' @param title The title for the plot
#' @param setSize The number of markers per plot
#' @import Seurat
PlotMarkerSeries <- function(seuratObj, reduction, features, title, setSize = 4) {
  featuresToPlot <- unique(intersect(features, row.names(seuratObj)))
  featuresToPlot <- RemoveUnchangedOrZero(seuratObj = seuratObj, reduction = reduction, features = featuresToPlot)
  steps <- ceiling(length(featuresToPlot) / setSize) - 1

  for (i in 0:steps) {
    start <- (i * 4) + 1
    end <- min((start + 3), length(featuresToPlot))
    genes <- featuresToPlot[start:end]


    PlotMarkerSet(seuratObj, reduction, title, genes)
  }
}


#' @title RemoveUnchangedOrZero
#' @description Starting with a set of features, remove any that are unchanged or zeros across all cells
#' @param seuratObj The seurat object
#' @param reduction The reduction to use
#' @param features The features to evalute
#' @return The updated set of features
#' @import Seurat
RemoveUnchangedOrZero <- function(seuratObj, reduction, features) {
  ret <- c()
  #Remove zeros or unchanged:
  for (feature in features) {
    dims <- paste0(Key(object = seuratObj[[reduction]]), c(1,2))
    data <- FetchData(object = seuratObj, vars = c(dims, features), cells = colnames(x = seuratObj))
    if (!all(data[, feature] == data[, feature][1])) {
      ret <- c(ret, feature)
    }
  }

  return(ret)
}


#' @title PlotMarkerSeries
#'
#' @description Generate a labeled FeaturePlot for the provided markers
#' @param seuratObj The seurat object
#' @param reduction The reduction to use
#' @param title The title for the plot
#' @param features The features to plot
#' @import Seurat
PlotMarkerSet <- function(seuratObj, reduction, title, features) {
  featuresToPlot <- intersect(features, row.names(seuratObj))
  featuresToPlot <- RemoveUnchangedOrZero(seuratObj, reduction, featuresToPlot)

  if (length(features) != length(featuresToPlot)){
    missing <- features[!(features %in% featuresToPlot)]
    print(paste0('The following features were requested, but not present: ', paste0(missing, collapse = ',')))
  }

  if (length(featuresToPlot) == 0){
    print('None of the requested features were present, skipping')
    return()
  }

  print(AddTitleToMultiPlot(FeaturePlot(seuratObj, features = featuresToPlot, reduction = reduction), title))
}
