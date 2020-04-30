#' @include LabKeySettings.R
#' @include Seurat_III_Fixes.R
#' @include Utils.R
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
#' @param storeGeneIds If true, a map to translate geneId and name (by default rownames will use gene name)
#' @return A Seurat object.
#' @keywords ReadAndFilter10X
#' @export
#' @importFrom Seurat Read10X
ReadAndFilter10xData <- function(dataDir, datasetName, emptyDropNIters=10000, storeGeneIds=TRUE) {
  if (!file.exists(dataDir)){
    stop(paste0("File does not exist: ", dataDir))
  }

  if (!dir.exists(dataDir)){
    stop(paste0("File is not a directory: ", dataDir))
  }

  seuratRawData <- Read10X(data.dir = dataDir, strip.suffix = TRUE)

  #Cannot have underscores in feature names, Seurat will replace with hyphen anyway.  Perform upfront to avoid warning
  if (sum(grepl(x = rownames(seuratRawData), pattern = '_')) > 0) {
    print('Replacing underscores with hyphens in feature names')
    rownames(seuratRawData) <- gsub(x = rownames(seuratRawData), pattern = '_', replacement = '-')
  }

  seuratRawData <- PerformEmptyDropletFiltering(seuratRawData, emptyDropNIters=emptyDropNIters)

  seuratObj <- CreateSeuratObj(seuratRawData, project = datasetName)
  PrintQcPlots(seuratObj)

  if (storeGeneIds) {
    #store IDs in assay metadata
    geneIds <- rownames(Read10X(data.dir = dataDir, gene.column = 1, strip.suffix = TRUE))
    names(geneIds) <- rownames(seuratObj)
    assayName <- DefaultAssay(seuratObj)
    seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds, col.name = 'GeneId')
  }

  return(seuratObj)
}


#' @title Retrieve the gene IDs from a seuratObj created using ReadAndFilter10xData.
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
#' @param seuratObj The seurat object
#' @param datasetName A name to use when creating the Seurat object
#' @param throwIfGenesNotFound If true and any of the requested gene names are not found, an error will be thrown.  Otherwise, the result will contain NAs
#' @return A named vector of the gene IDs
#' @export
GetGeneIds <- function(seuratObj, geneNames, throwIfGenesNotFound = TRUE) {
	ret <- NULL

  featureMeta <- GetAssay(seuratObj)@meta.features
  if ('GeneId' %in% colnames(featureMeta)) {
    ret <- featureMeta$GeneId
    names(ret) <- rownames(seuratObj)
    ret <- ret[geneNames]
  }
  #NOTE: in previous versions we stored geneIDs here:
	else if ('geneIds' %in% names(seuratObj@misc)) {
    ret <- seuratObj@misc$geneIds[geneNames]
  }

  if (is.null(ret)) {
  	stop('Expected gene IDs to be stored under GetAssay(seuratObj)@meta.features or seuratObj@misc$geneIds')
  }

  if (throwIfGenesNotFound & sum(is.na(ret)) > 0) {
    notFound <- paste0(geneNames[is.na(ret)], collapse = ',')
    stop(paste0('Gene names not found: ', notFound))
  } else {
    names(ret) <- geneNames
  }

  return(ret)
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
#' @return A modified Seurat object.
MarkStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)
  }

  #Write to logfile to keep record in case this crashes
  .WriteLogMsg(paste0('Step complete: ', name))

  return(seuratObj)
}



#' @title doMergeCCA
#' @description An internal method to do CCA merging
#' @param seuratObjs A list of seurat objects
#' @param nameList A list of names
#' @param maxCCAspaceDim The number of dims to use with FindIntegrationAnchors()
#' @param maxPCs2Weight The number of dims to use with IntegrateData()
#' @param assay NULL for default assay or specify
#' @param useAllFeatures If true, the resulting object will contain all features, as opposed to just VariableGenes (not recommended)
#' @param nVariableFeatures The number of VariableFeatures to identify
#' @param includeCellCycleGenes If true, the cell cycles genes will always be included with IntegrateData(), as opposed to just VariableGenes
#' @param spike.genes If NULL ignored, but a vector of unique genes
#' @param spikeImmuneGenes If CCA is used, this will spike a pre-determined set of immune-related genes into the merged object
#' @param normalization.method Normalization method
#' @param plotFigs If true, QC figures will be generated
#' @return A modified Seurat object.
doMergeCCA <- function(seuratObjs, nameList, 
                       maxCCAspaceDim, 
                       assay = NULL, 
                       normalization.method = "LogNormalize", 
                       useAllFeatures = F, includeCellCycleGenes = T, nVariableFeatures = NULL,
                       plotFigs = T, spike.genes = NULL, spikeImmuneGenes=T, maxPCs2Weight=NULL){


  if (length(seuratObjs) == 1) {
		print("Only one file in list. no need to merge. returning the single object")
		return(seuratObjs[[1]])
	}

  #Pre process each object in list of objects
	for (exptNum in nameList) {
		print(paste0('adding dataset: ', exptNum))
		so <- seuratObjs[[exptNum]]
		if (!HasStepRun(so, 'NormalizeData')) {
			print('Normalizing')
			so <- NormalizeData(object = so, verbose = F, normalization.method = normalization.method)
		} else {
			print('Normalization performed')
		}

		if (!HasStepRun(so, 'FindVariableFeatures')) {
			print('Finding variable features')
			so <- FindVariableFeatures(object = so, verbose = F, selection.method = "vst", nfeatures = nVariableFeatures)
		} else {
			print('FindVariableFeatures performed')
		}

		# No scaling needed at this step.

		if (plotFigs) {
			print(LabelPoints(plot = VariableFeaturePlot(so), points = head(VariableFeatures(so), 20), repel = TRUE, xnudge = 0, ynudge = 0))
		}

		seuratObjs[[exptNum]] <- so
  }
  
  seuratObj <- NULL
	CheckDuplicatedCellNames(seuratObjs)
    
	print("Performing FindIntegrationAnchors to find anchors...")
    
	# dims here means : Which dimensions to use from the CCA to specify the neighbor search space
	anchors <- FindIntegrationAnchors(object.list = seuratObjs, dims = 1:maxCCAspaceDim, verbose = F, assay = assay)

	allFeatures <- rownames(seuratObjs[[1]])
    
	#if all of the objects have the same genes, then above suffice
	#if objects have different genes, their intersect is what we want, i.e., genes represented across samples
	for (i in 2:length(seuratObjs)) {
		allFeatures <- intersect(allFeatures, rownames(seuratObjs[[i]]))
	}

	#run using intersection of all features
	if (useAllFeatures) {
		features <- allFeatures
		print(paste0('Total features in common: ', length(features)))
	} else {
		#use the variable genes, across all samples, get the intersect
		features <- VariableFeatures(seuratObjs[[1]])

		for (i in 2:length(seuratObjs)) {
			features <- intersect(features, VariableFeatures(seuratObjs[[i]]))
		}

		#now add the anchor.features, should be the same as variable features above, unless pipeline is changed
		features <- unique(c(features, slot(object = anchors, name = "anchor.features")))

		if (includeCellCycleGenes) features <- unique(c(features, .GetSPhaseGenes(), .GetG2MGenes()))
		if (!is.null(spike.genes)) features <- unique(c(features, spike.genes))
    if (spikeImmuneGenes) {
      SGS.LS <- Phenotyping_GeneList()
      immuneSpikeGenes <- unique(c(
        SGS.LS$WBC,
        SGS.LS$Erythrocyte,
        SGS.LS$Lymphoid,
        SGS.LS$TCellCanonical, SGS.LS$TCellSecondary,
        SGS.LS$CD8Canonical, SGS.LS$CD8Subphenos1, SGS.LS$MAIT,
        SGS.LS$CD4Canonical, SGS.LS$CD4Subphenos1, SGS.LS$CD4Subphenos2,
        SGS.LS$BCellCanonical, SGS.LS$BCellSecondary, SGS.LS$BCellCanonicalV3,
        SGS.LS$NKCanonical,
        SGS.LS$HA_ImpGenes, SGS.LS$HA_CART_gini, SGS.LS$HighlyActivated3, SGS.LS$LessActivated,
        SGS.LS$Myeloid,
        SGS.LS$Monocytes, SGS.LS$MonocytesCD34p, SGS.LS$MonocytesFCGR3A,
        SGS.LS$Vilani_Mono1_classical_CD14high_CD16neg,
        SGS.LS$Vilani_Mono2_nonclassical_CD14posCD16high,
        SGS.LS$Vilani_Mono3_undef1,
        SGS.LS$Vilani_Mono4_undef2,
        SGS.LS$Macrophage
      ))

      features <- unique(c(features, immuneSpikeGenes))
	  }
  }

	#canonicals or other spikes may not be in the set of genes, so errors are given. Needs to remove genes not in data.
	features <- features[features %in% allFeatures]

	print(paste0("Number of final feature kept: ", length(features)))

	print("Starting IntegrateData to combine using CCA")

	# dims here means : #Number of PCs to use in the weighting procedure
	seuratObj <- IntegrateData(anchorset = anchors, dims = 1:maxPCs2Weight,
														 verbose = F,
														 features.to.integrate = features,
														 new.assay.name = "Integrated")

	DefaultAssay(seuratObj) <- "Integrated"

	# This will prevent repeating this step downstream
	seuratObj <- MarkStepRun(seuratObj, 'NormalizeData')

	return(seuratObj)
}


#' @title doMergeSimple
#' @description An internal method to do simple merging of seurat objects from a list of them.
#' @param seuratObjs A list of seurat objects, optionally named (in which case these will be used as dataset names). Also can use SplitObject(, split.by =)
#' @param nameList A list of names from MergeSeuratObjs()
#' @param projectName The projectName to pass to Seurat
doMergeSimple <- function(seuratObjs, nameList, projectName){
  seuratObj <- NULL

  for (exptNum in nameList) {
    print(exptNum)
    if (is.null(seuratObj)) {
      seuratObj <- seuratObjs[[exptNum]]
    } else {
			assayName <- DefaultAssay(seuratObj)
      geneIds1 <- GetAssay(seuratObj)@meta.features$GeneId
			names(geneIds1) <- rownames(seuratObj)
      geneIds2 <- GetAssay(seuratObjs[[exptNum]])@meta.features$GeneId
			names(geneIds2) <- rownames(seuratObjs[[exptNum]])

      if (any(rownames(seuratObj) != rownames(seuratObjs[[exptNum]]))) {
        stop('Gene names are not equal!')
      }

      seuratObj <- merge(x = seuratObj,
                         y = seuratObjs[[exptNum]],
                         project = projectName)

      if (any(is.na(geneIds1)) & !any(is.na(geneIds2))) {
				seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds2, col.name = 'GeneId')
      } else if (!any(is.na(geneIds1)) & any(is.na(geneIds2))) {
				seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
      } else if (!any(is.na(geneIds1)) & !any(is.na(geneIds2))) {
        if (any(geneIds1 != geneIds2)) {
          stop('Gene IDs did not match between seurat objects!')
        }

        seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
      } else {
        seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
			}
    }
  }
  
  return(seuratObj)
}



#' @title MergeSeuratObjs
#' @description Merges a list of Seurat objects, using Seurat::IntegrateData()
#' @param seuratObjs A list of seurat objects, optionally named (in which case these will be used as dataset names). Also can use SplitObject(, split.by =)
#' @param metadata A list of metadata.  If provided, the names of this list will be used as dataset names
#' @param method A string either simple or cca
#' @param maxCCAspaceDim The number of dims to use with FindIntegrationAnchors()
#' @param maxPCs2Weight The number of dims to use with IntegrateData()
#' @param projectName The project name when creating the final seuratObj
#' @param useAllFeatures If true, the resulting object will contain all features, as opposed to just VariableGenes (not recommended due to complexity). Also having all genes, means possibly more unwanted noise. Consider Spiking interesting/canonical genes as oppose to all. 
#' @param nVariableFeatures The number of VariableFeatures to identify
#' @param includeCellCycleGenes If true, the cell cycles genes will always be included with IntegrateData(), as opposed to just VariableGenes
#' @param spike.genes If CCA is used, this list of list of genes will be spiked into the merged object, beyond the default VariableGenes()
#' @param assay The assay to use
#' @param normalization.method Normalization method
#' @param spikeImmuneGenes If CCA is used, this will spike a pre-determined set of immune-related genes into the merged object
#' @return A modified Seurat object.
#' @export
#' @importFrom methods slot
MergeSeuratObjs <- function(seuratObjs, metadata=NULL, 
                            method = c("simple", "cca"),
                            maxCCAspaceDim = 20, maxPCs2Weight = 20, 
                            useAllFeatures = F, nVariableFeatures = 2000,
                            includeCellCycleGenes = T, assay = NULL,
                            normalization.method = "LogNormalize",
                            spike.genes = NULL, spikeImmuneGenes = T){


  method <- match.arg(method)

  print(paste0("Starting merge.  Method: ", method))

  nameList <- NULL
  if (is.null(metadata)){
    nameList <- names(seuratObjs)
  } else {
    nameList <- names(metadata)
  }

  #ensure barcodes unique:
  for (exptNum in nameList) {
    print(paste0('adding dataset: ', exptNum))
    prefix <- paste0(exptNum)
    so <- seuratObjs[[exptNum]]
    if (!('BarcodePrefix' %in% names(so@meta.data))) {
      print(paste0('Adding barcode prefix: ', prefix))
      so <- RenameCells(object = so, add.cell.id = prefix)
      so[['BarcodePrefix']] <- c(prefix)
      seuratObjs[[exptNum]] <- so
    } else {
      print('Barcode prefix already added')
    }
  }
  
  if (method == "cca"){
    seuratObj <- doMergeCCA(seuratObjs = seuratObjs, 
                            nameList = nameList, 
                            plotFigs = T, 
                            maxCCAspaceDim = maxCCAspaceDim, 
                            assay = assay, 
                            normalization.method = normalization.method, 
                            useAllFeatures = useAllFeatures,
                            includeCellCycleGenes = includeCellCycleGenes,
                            nVariableFeatures = nVariableFeatures, 
                            spike.genes = spike.genes,
                            spikeImmuneGenes = spikeImmuneGenes,
                            maxPCs2Weight=maxPCs2Weight)
  } else if (method == "simple") {
    seuratObj <- doMergeSimple(seuratObjs = seuratObjs, 
                                           nameList = nameList, 
                                           projectName = "simplemerge")
  }

  #store method used
  seuratObj@misc['MergeMethod'] <- method

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
#' @param nCount_RNA.high If doCellFilter=T, cells with nCount_RNA above this value will be filtered
#' @param nCount_RNA.low If doCellFilter=T, cells with nCount_RNA below this value will be filtered
#' @param nFeature.high If doCellFilter=T, cells with nFeature above this value will be filtered
#' @param nFeature.low If doCellFilter=T, cells with nFeature below this value will be filtered
#' @param pMito.high If doCellFilter=T, cells with percent mito above this value will be filtered
#' @param pMito.low If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param forceReCalc If true, all steps will be repeated even if already marked as complete
#' @param variableGeneTable If provided, a table of variable genes will be written to this file
#' @param variableFeatureSelectionMethod The selection method to be passed to FindVariableFeatures()
#' @param useSCTransform If true, SCTransform will be used in place of the standard Seurat workflow (NormalizeData, ScaleData, FindVariableFeatures)
#' @param nVariableFeatures The number of variable features to find
#' @param dispersion.cutoff Passed directly to FindVariableFeatures
#' @param mean.cutoff Passed directly to FindVariableFeatures
#' @param spikeGenes If provided these will be appended to the set of VariableFeatures
#' @param printDefaultPlots If true, the default set of QC plots will be printed
#' @param npcs Number of PCs to use for RunPCA()
#' @param ccPcaResultFile If provided, the PCA results from cell cycle regression will be written here
#' @return A modified Seurat object.
#' @export
ProcessSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
                            nCount_RNA.high = 20000, nFeature.high = 3000, pMito.high = 0.15,
                            nCount_RNA.low = 0.99, nFeature.low = 200, pMito.low = -Inf, forceReCalc = F,
                            variableGeneTable = NULL, variableFeatureSelectionMethod = 'vst', 
                            nVariableFeatures = 2000, printDefaultPlots = T,
                            npcs = 50, ccPcaResultFile = NULL, useSCTransform = F, 
                            mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), 
                            spikeGenes = NULL){

	if (!forceReCalc && HasStepRun(seuratObj, 'ProcessSeurat1', forceReCalc = forceReCalc)) {
    if (printDefaultPlots){
      .PrintSeuratPlots(seuratObj, doCellCycle)
    }

		return(seuratObj)
	}

  if (doCellFilter & (forceReCalc | !HasStepRun(seuratObj, 'FilterCells', forceReCalc = forceReCalc))) {
    seuratObj <- .DoCellFilter(seuratObj = seuratObj,
			nCount_RNA.high = nCount_RNA.high,
			nFeature.high = nFeature.high,
			pMito.high = pMito.high,
			nCount_RNA.low = nCount_RNA.low,
			nFeature.low = nFeature.low,
			pMito.low = pMito.low
    )

    seuratObj <- MarkStepRun(seuratObj, 'FilterCells', saveFile)
  }

  if (!useSCTransform) {
    if (forceReCalc | !HasStepRun(seuratObj, 'NormalizeData', forceReCalc = forceReCalc)) {
      print('Normalizing data:')
      seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'NormalizeData', saveFile)
    }

    if (forceReCalc | !HasStepRun(seuratObj, 'FindVariableFeatures', forceReCalc = forceReCalc)) {
      print('Find variable features:')
      seuratObj <- FindVariableFeatures(object = seuratObj, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff , verbose = F, selection.method = variableFeatureSelectionMethod, nfeatures = nVariableFeatures)
      seuratObj <- MarkStepRun(seuratObj, 'FindVariableFeatures', saveFile)
    }

    if (forceReCalc | !HasStepRun(seuratObj, 'ScaleData', forceReCalc = forceReCalc)) {
      print('Scale data:')
      seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "p.mito"), verbose = F)
      seuratObj <- MarkStepRun(seuratObj, 'ScaleData', saveFile)
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
  
  if (!is.null(spikeGenes)){
    VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), spikeGenes))
  }

  vg <- VariableFeatures(object = seuratObj)
  
  if (forceReCalc | !HasStepRun(seuratObj, 'RunPCA', forceReCalc = forceReCalc)) {
    seuratObj <- RunPCA(object = seuratObj, features = vg, verbose = F, npcs = npcs)
    seuratObj <- MarkStepRun(seuratObj, 'RunPCA', saveFile)
  }

  if (forceReCalc | !HasStepRun(seuratObj, 'ProjectDim', forceReCalc = forceReCalc)) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- MarkStepRun(seuratObj, 'ProjectDim', saveFile)
  }

  #Verify data exists.  This appears to get reset, possibly by DimRedux steps
  runJackStraw <- forceReCalc | !HasStepRun(seuratObj, 'JackStraw', forceReCalc = forceReCalc)
  if (!useSCTransform) {
    if (!runJackStraw && length(seuratObj@reductions$pca@jackstraw$empirical.p.values) == 0) {
      warning('JackStraw marked as complete, but seurat object lacks data')
      runJackStraw <- TRUE
    }
  }

  if (!useSCTransform && runJackStraw) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- MarkStepRun(seuratObj, 'JackStraw', saveFile)
  }

  print(paste0('Total variable genes: ', length(vg)))
  if (!is.null(variableGeneTable)){
    write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
  }

  if (printDefaultPlots){
    .PrintSeuratPlots(seuratObj, doCellCycle)
  }

  seuratObj <- MarkStepRun(seuratObj, 'ProcessSeurat1', saveFile = saveFile)

  print(seuratObj)

  return(seuratObj)
}


.DoCellFilter <- function(seuratObj, nCount_RNA.high = 20000, nFeature.high = 3000, pMito.high = 0.15, nCount_RNA.low = 0.99, nFeature.low = 200, pMito.low = -Inf) {
  print("Filtering Cells...")
  seuratObj@misc$OriginalCells <- length(colnames(x = seuratObj))
  print(paste0('Initial cells: ', length(colnames(x = seuratObj))))

  P1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "p.mito")
  P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_hline(aes(yintercept=pMito.low), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_hline(aes(yintercept=pMito.high), color="blue", linetype="dashed", size=1)
  print(P1)

  P1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_hline(aes(yintercept=nFeature.low), color="blue", linetype="dashed", size=1)
  P1 <- P1 + geom_hline(aes(yintercept=nFeature.high), color="blue", linetype="dashed", size=1)
  print(P1)

  #See: https://github.com/satijalab/seurat/issues/1053#issuecomment-454512002
  expr <- Seurat::FetchData(object = seuratObj, vars = 'nCount_RNA')
  seuratObj <- seuratObj[, which(x = expr > nCount_RNA.low & expr < nCount_RNA.high)]
  print(paste0('After nCount_RNA filter: ', length(colnames(x = seuratObj))))

  expr <- Seurat::FetchData(object = seuratObj, vars = 'nFeature_RNA')
  seuratObj <- seuratObj[, which(x = expr > nFeature.low & expr < nFeature.high)]
  print(paste0('After nFeature_RNA filter: ', length(colnames(x = seuratObj))))

  expr <- Seurat::FetchData(object = seuratObj, vars = 'p.mito')
  if (!all(is.na(expr)) && max(expr) != 0) {
		seuratObj <- seuratObj[, which(x = expr > pMito.low & expr < pMito.high)]
		print(paste0('After p.mito filter: ', length(colnames(x = seuratObj))))
  } else {
    print('Either p.mito was NA or all values were 0')
  }

  print(paste0('Final: ', length(colnames(x = seuratObj))))
  
  return(seuratObj)
}


.PrintSeuratPlots <- function(seuratObj, doCellCycle) {
  print(VizDimLoadings(object = seuratObj, dims = 1:2))
  print(LabelPoints(plot = VariableFeaturePlot(seuratObj), points = head(VariableFeatures(seuratObj), 20), repel = TRUE, xnudge = 0, ynudge = 0))

  print(DimPlot(object = seuratObj))
  if (doCellCycle && ('Phase' %in% names(seuratObj@meta.data))) {
    print(cowplot::plot_grid(plotlist = list(
      DimPlot(object = seuratObj, reduction = "pca", dims = c(1, 2), group.by = 'Phase'),
      DimPlot(object = seuratObj, reduction = "pca", dims = c(2, 3), group.by = 'Phase'),
      DimPlot(object = seuratObj, reduction = "pca", dims = c(3, 4), group.by = 'Phase'),
      DimPlot(object = seuratObj, reduction = "pca", dims = c(4, 5), group.by = 'Phase')
    )))
  }

  print(DimHeatmap(object = seuratObj, dims = 1, cells = 500, balanced = TRUE, fast = F) + NoLegend())
  Try2Prit <-try(print(DimHeatmap(object = seuratObj, dims = 1:20, cells = 500, balanced = TRUE, fast = F) + NoLegend()), silent = T)
  if(class(Try2Prit) == "try-error") try(print(DimHeatmap(object = seuratObj, dims = 1:6, cells = 500, balanced = TRUE, fast = F) + NoLegend()), silent = T)

  if (length(seuratObj@reductions$pca@jackstraw$empirical.p.values) == 0) {
    print('Unable to display JackStrawPlot, data not available')
  } else {
    print(JackStrawPlot(object = seuratObj, dims = 1:20))
  }

  print(ElbowPlot(object = seuratObj))
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

  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- .GetSPhaseGenes()
  g2m.genes <- .GetG2MGenes()

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
#' @param maxTsneIter The value of max_iter to provide to RunTSNE.  Increasing can help large datasets.
#' @param forceReCalc If true, all steps will be performed even if already marked complete
#' @param umap.method The UMAP method, either uwot or umap-learn, passed directly to Seurat::RunUMAP
#' @return A modified Seurat object.
#' @export
FindClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, saveFile = NULL, forceReCalc = F, minDimsToUse = NULL, umap.method = 'umap-learn',
                                   UMAP_NumNeib = 40L, UMAP_MinDist = 0.2, UMAP_Seed = 1234, UMAP.NumEpoc = 500, maxTsneIter = 2000) {
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
    perplexity <- .InferPerplexityFromSeuratObj(seuratObj)
    seuratObj <- RunTSNE(object = seuratObj, dims.use = dimsToUse, check_duplicates = FALSE, perplexity = perplexity, max_iter = maxTsneIter)
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
#' @importFrom dplyr %>% coalesce group_by summarise filter top_n select everything
#' @import DESeq2
#' @import MAST
#' @importFrom knitr kable
#' @export
Find_Markers <- function(seuratObj, resolutionToUse, outFile = NULL, saveFileMarkers = NULL,
testsToUse = c('wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'),
numGenesToSave = 20, onlyPos = F, includeNameTranslation = T) {

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

          print(paste0('Total genes: ', nrow(toBind)))

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

    if (all(is.na(seuratObj.markers))) {
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
    if (nrow(toWrite) == 0) {
      print('No significant markers were found')
    } else {
      if (includeNameTranslation) {
        toWrite$translatedName <- AliasGeneNames(toWrite$gene)
        toWrite <- toWrite %>% select(gene, translatedName, everything())
      }

      if (!is.null(outFile)) {
        write.table(toWrite, file = outFile, sep = '\t', row.names = F, quote = F)
      }

      print(DimPlot(object = seuratObj, reduction = 'tsne'))

      topGene <- toWrite %>% group_by(cluster, test) %>% top_n(20, avg_logFC)
      print(DoHeatmap(object = seuratObj, features = unique(as.character(topGene$gene))))
  		print(kable(topGene))
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
#' @param file The output file
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

  for (lev in levels(AvlLevels)){
    print(lev)
    ClustLS[[lev]] <- as.data.frame(Matrix::rowMeans(GetAssayData(object = seuratObj, features = genes, slot = slot)[genes, colnames(seuratObj)[which(AvlLevels==lev)]  ]))
  }

  ClustDF <- as.data.frame(ClustLS)
  colnames(ClustDF) <- paste0("clus", levels(AvlLevels))

  return(ClustDF)
}


#' @title PlotImmuneMarkers
#'
#' @description Generate a set of Seurat FeaturePlots for common immune cell markers
#' @param seuratObj A seurat object
#' @param reduction The reduction to use
#' @export
#' @import Seurat
PlotImmuneMarkers <- function(seuratObj, reduction = 'tsne') {
  #ENSMMUG00000003532=CD8b
  PlotMarkerSet(seuratObj, reduction, 'CD8/CD4 Markers', c('CD8A', 'ENSMMUG00000003532', 'CD4', 'IL7R'))

  #Eff v. Mem:
  #IL7R = CD127
  #IL2RA = CD25
  #PTPRC = CD45
  #SELL = CD62-L / CD-197
  PlotMarkerSeries(seuratObj, reduction, c('CCR7', 'SELL', 'GZMB', 'CCR5', 'IL2RA', 'PTPRC', 'IL7R', 'CTLA4'), 'Effector vs. Memory')

  #CD8 Activation
  PlotMarkerSeries(seuratObj, reduction, c('IFNG', 'CD69', 'TNF', 'NFKBID', 'LTB', 'TNFRSF9', 'CCL4L1', 'NR4A3', 'TNFSF14', 'CD82', 'PIGT', 'IRF8', 'IRF4', 'RGCC', 'PD1', 'PDCD1', 'TNFAIP3', 'ENSMMUG00000013779', 'ENSMMUG00000060218', 'ENSMMUG00000008111'), 'CD8 Activation Markers')

  PlotMarkerSeries(seuratObj, reduction, c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'), 'Cytotoxicity')

  PlotMarkerSet(seuratObj, reduction, 'B-cell Markers', c('MS4A1', 'CD79A', 'CD74', 'DRA'))

  PlotMarkerSet(seuratObj, reduction, 'Monocyte', c('LYZ', 'CST3', 'S100A6', 'VIM'))

  PlotMarkerSet(seuratObj, reduction, 'Transcription Factors', c('TBX21', 'GATA3', 'RORC', 'FOXP3', 'BCL6', 'EOMES', 'TOX'))

  PlotMarkerSet(seuratObj, reduction, 'Inhibitory Markers', c('TIGIT', 'CTLA4', 'BTLA', 'PDCD1', 'CD274'))

  #DAP10/12
  PlotMarkerSet(seuratObj, reduction, 'Signaling', c('HCST', 'TYROBP', 'SYK', 'ZAP70'))

  #LILR/KIR:
  PlotMarkerSeries(seuratObj, reduction, c('LILRA5','LILRA6','LILRB4','LILRB5','KIR2DL4','KIR3DX1', 'MAMU-KIR', 'KIR2DL4', 'KIR3DL2'), 'KIR/LILR')

  PlotMarkerSeries(seuratObj, reduction, c('FCGR1A','FCGR2B','FCGR3'), 'FCGR')

  #Cytokines
  cytokines <- c('IL1A','IL1B','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL2','IL2RA','IL2RB','IL2RG','IL3','IL3RA','IL4','IL4I1','IL4R','IL5','IL5RA','IL6','IL6R','IL6ST','IL7','IL7R','IL9','IL10','IL10RA','IL11','IL12A','IL12B','IL12RB1','IL12RB2','IL13','IL13RA2','IL15','IL15Ra','IL16','IL17A','IL17B','IL17C','IL17D','IL17F','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL18BP','IL18R1','IL18RAP','IL19','IL20','IL20RA','IL20RB','IL21','IL21R','IL22','IL22RA2','IL23A','IL24','IL25','IL26','IL27','IL27RA','IL31','IL31RA','IL33','IL34','IL36A','IL36B','IL36G','IL37','ILDR1','ILDR2','ILF2','ILF3','ILK','ILKAP','ILVBL')
  PlotMarkerSeries(seuratObj, reduction, cytokines, 'Cytokines/Receptors')

  klrs <- c('KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRF2', 'KLRG1', 'KLRG2', 'ENSMMUG00000050862')
  PlotMarkerSeries(seuratObj, reduction, klrs, 'KLRs')

  PlotMarkerSet(seuratObj, reduction, 'Resident Memory', c('ITGAE', 'CD69', 'CXCR6'))

  #chemokines
  chemokines <- c('CCL1','CCL11','CCL13','CCL16','CCL17','CCL18','CCL19','CCL2','CCL20','CCL21','CCL22','CCL24','CCL25','CCL26','CCL27','CCL28','CCL5','CCL7','CCL8')
  chemokines <- c(chemokines, c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CCRL2'))
  chemokines <- c(chemokines, c('CXCL1','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17','CXCL5','CXCL6','CXCL8','CXCL9','CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','XCR1'))

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

  print(AddTitleToMultiPlot(FeaturePlot(seuratObj, features = featuresToPlot, reduction = reduction, min.cutoff = 'q05', max.cutoff = 'q95'), title))
}

