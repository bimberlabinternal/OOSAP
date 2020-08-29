#' @include LabKeySettings.R
#' @include Utils.R
#' @import Seurat
#' @import Rlabkey

#' @title DownloadAndAppendCiteSeq
#'
#' @description Downloads matching Cite-seq counts using barcodePrefix on the seurat object
#' @param seuratObj, A Seurat object.
#' @param featureLabelTable An optional TSV file listing the markers and metadata.  It must contain a header row with at least the columns (lowercase): tagname, sequence. If it contains the column markername, this will be used to replace the rownames of the matrix.
#' @param minRowSum  If provided, any ADTs with a rowSum less than this value will be dropped.
#' @param adtWhitelist An optional character vector of tag names, which must exactly match the rownames of the ADT count matrix.  If provided, data will be limited to these tags.
#' @param assayName The name of the assay to store the ADT data.
#' @return A modified Seurat object.
#' @export
DownloadAndAppendCiteSeq <- function(seuratObj, outPath = '.', assayName = 'ADT', featureLabelTable = NULL, minRowSum = 5, adtWhitelist = NULL){
	if (is.null(seuratObj[['BarcodePrefix']])){
		stop('Seurat object lacks BarcodePrefix column')
	}

	for (barcodePrefix in unique(unique(unlist(seuratObj[['BarcodePrefix']])))) {
		print(paste0('Possibly adding cell hashing data for prefix: ', barcodePrefix))

		citeseqId <- .FindMatchedCiteSeq(barcodePrefix)
		if (is.null(citeseqId)){
			print(paste0('CITE-seq not used for prefix: ', barcodePrefix, ', skipping'))
			next
		} else if (is.na(citeseqId)){
			stop(paste0('Unable to find CITE-seq matrix file for prefix: ', barcodePrefix))
		}

		countDir <- file.path(outPath, paste0(barcodePrefix, '_citeseqCounts'))
		if (!dir.exists(countDir)) {
			dir.create(countDir)
		}

		countDir <- .DownloadCiteSeqDir(outputFileId = citeseqId, localBaseDir = countDir)
		if (!dir.exists(countDir)){
			stop(paste0('Unable to download calls table for prefix: ', barcodePrefix, ', expected file: ', countDir))
		}

		seuratObj <- AppendCiteSeq(seuratObj = seuratObj, countMatrixDir = countDir, barcodePrefix = barcodePrefix, assayName = assayName, featureLabelTable = featureLabelTable, adtWhitelist = adtWhitelist, minRowSum = minRowSum, skipNormalize = T)

		seuratObj <- NormalizeData(seuratObj, assay = assayName, normalization.method = "CLR")
		seuratObj <- ScaleData(seuratObj, assay = assayName)
	}

	.PlotCiteSeqCountData(seuratObj, assayName = assayName)

	return(seuratObj)
}


.DownloadCiteSeqDir <- function(outputFileId, localBaseDir = './', overwriteFiles = T, mergeFolders = F) {
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

	#The database tracks the matrix, but we want the whole folder
	remotePath <- rows[['dataid_webdavurlrelative']]
	remotePath <- gsub(x = remotePath, pattern = 'matrix.mtx.gz', replacement = '')

	success <- labkey.webdav.downloadFolder(
		baseUrl=lkBaseUrl,
		folderPath=paste0(lkDefaultFolder,wb),
		remoteFilePath = remotePath,
		overwriteFiles = overwriteFiles,
		mergeFolders = mergeFolders,
		localBaseDir = localBaseDir
	)

	if (!success || !dir.exists(localBaseDir) || !file.exists(paste0(localBaseDir, '/umi_count/matrix.mtx.gz'))) {
		stop(paste0('labkey.webdav.downloadFolder failed, expected: ', localBaseDir))
	}

	return(paste0(localBaseDir, '/umi_count'))
}

#' @import Rlabkey
.FindMatchedCiteSeq <- function(loupeDataId){
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
		colSelect="rowid,readsetid,citeseqreadsetid",
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(cDNAs) == 0) {
		stop(paste0('No cDNA records found for GEX readset: ', readset))
	} else if (sum(!is.na(cDNAs$citeseqreadsetid)) == 0) {
		print(paste0('The cDNA library does not use citeseq, aborting'))
		return(NULL)
	}

	rows <- suppressWarnings(labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=lkDefaultFolder,
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		colSort="-rowid",
		colSelect="rowid",
		colFilter=makeFilter(c("readset", "EQUAL", readset),
		c("category", "EQUAL", "Seurat CITE-Seq Count Matrix"),
		c("library_id", "EQUAL", libraryId)),
		containerFilter=NULL,
		colNameOpt="rname"
	))

	ret <- NULL
	if (nrow(rows) == 0){
		print(paste0("Output of type 'Seurat CITE-Seq Count Matrix' not found.  Readset: ", readset, ", libraryId: ", libraryId))
	} else {
		ret <- rows[1]
	}

	if (all(is.null(ret))){
		print("Trying to find output of type: 'CITE-Seq Count Matrix' using GEX readset")

		rows <- suppressWarnings(labkey.selectRows(
			baseUrl=lkBaseUrl,
			folderPath=lkDefaultFolder,
			schemaName="sequenceanalysis",
			queryName="outputfiles",
			colSort="-rowid",
			colSelect="rowid,",
			colFilter=makeFilter(c("readset", "EQUAL", readset),
				c("category", "EQUAL", "CITE-Seq Count Matrix"),
				c("library_ld", "EQUAL", libraryId)
			),
			containerFilter=NULL,
			colNameOpt="rname"
		))

		if (nrow(rows) == 0){
			print("Not found")
		} else {
			ret <- rows[1]
			print("Found output")
		}
	}

	if (all(is.null(ret))){
		print("Trying to find output of type: 'CITE-Seq Count Matrix', using cite-seq readset")
		citeseqReadset <- unique(cDNAs$citeseqreadsetid)

		rows <- suppressWarnings(labkey.selectRows(
			baseUrl=lkBaseUrl,
			folderPath=lkDefaultFolder,
			schemaName="sequenceanalysis",
			queryName="outputfiles",
			colSort="-rowid",
			colSelect="rowid,",
			colFilter=makeFilter(c("readset", "EQUAL", citeseqReadset),
				c("category", "EQUAL", "CITE-Seq Count Matrix"),
				c("library_ld", "EQUAL", libraryId)
			),
			containerFilter=NULL,
			colNameOpt="rname"
		))

		if (nrow(rows) == 0){
			print("Not found")
		} else {
			ret <- rows[1]
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

#' @importFrom dplyr arrange
AppendCiteSeq <- function(seuratObj, countMatrixDir, barcodePrefix = NULL, assayName = 'ADT', minRowSum = 5, featureLabelTable = NULL, adtWhitelist = NULL, skipNormalize = F) {
	gexCells <- colnames(seuratObj)
	print(paste0('Initial cell barcodes in GEX data: ', length(gexCells)))

	if (!dir.exists(countMatrixDir))
	stop("Count matrix not found")

	bData <- Read10X(countMatrixDir, gene.column=1, strip.suffix = TRUE)
	bData <- bData[which(!(rownames(bData) %in% c('unmapped'))), , drop = F]
	if (!is.null(barcodePrefix)) {
		colnames(bData) <- paste0(barcodePrefix, '_', colnames(bData))

		gexCells <- gexCells[seuratObj$BarcodePrefix == barcodePrefix]
		print(paste0('Initial cell barcodes in GEX data for prefix: ', length(gexCells)))
	}

	print(paste0('Initial cells in cite-seq matrix: ', ncol(bData)))
	bData <- bData[,which(colnames(bData) %in% gexCells), drop = F]
	print(paste0('Intersect with GEX data: ', ncol(bData)))
	if (ncol(bData) == 0) {
		print(head(gexCells))
		print(head(colnames(bData)))
		stop('No cells are shared')
	}
	bData <- as.sparse(bData)

	# Dont filter if a whitelist is used.
	toDrop <- rowSums(bData) < minRowSum
	if (is.null(adtWhitelist) && sum(toDrop) > 0){
		print(paste0('ADTs dropped due to low counts across cells: ', sum(toDrop)))
		print(paste0(rownames(bData)[toDrop], collapse = ','))
		bData <- bData[!toDrop, , drop = FALSE]
		print(paste0('ADTs after filter: ', nrow(bData)))
		if (nrow(bData) == 0) {
			stop('No ADTs remain after filter')
		}
	}

	# Apply this prior to rename:
	if (!is.null(adtWhitelist)) {
		print('Filtering HTOs based on whitelist')
		bData <- bData[rownames(bData) %in% adtWhitelist, , drop = FALSE]
		print(paste0('ADTs after filter: ', nrow(bData)))
		if (nrow(bData) == 0) {
			stop('No rows left after filter!')
		}
	}

	#rename features based on DB
	ft <- NULL
	shouldAddMetadata <- F
	if (!is.null(featureLabelTable)) {
		ft <- read.table(featureLabelTable, sep = '\t', header = T, stringsAsFactors = F)
		if (!('tagname' %in% colnames(ft))) {
			stop('Table must contain the column tagname')
		}

		rownames(ft) <- ft$featureName

		if ('tagname' %in% colnames(ft)){
			print('Renaming ADTs')
			ft$tagname <- paste0(ft$tagname, '-', ft$sequence)
			newRows <- data.frame(tagname = rownames(bData), sortorder = 1:nrow(bData), stringsAsFactors = F)
			newRows <- merge(newRows , ft, by = 'tagname', all.x = T, all.y = F)
			newRows <- newRows %>% arrange(sortorder)
			newRows <- newRows[names(newRows) != 'sortorder']

			newRows$markername <- dplyr::coalesce(newRows$markername, newRows$tagname)
			newRows$markernameUnique <- make.unique(newRows$markername)

			print(paste0('Total renamed: ', sum(newRows$markernameUnique != newRows$tagname)))

			rownames(bData) <- newRows$markernameUnique
			ft <- newRows
			shouldAddMetadata <- T
		}
	}

	#Note: we need to support merging with any existing data and account for barcodePrefix:
	assayData <- NULL
	if (assayName %in% names(seuratObj@assays)) {
		print('existing assay found, merging counts')
		assayData <- GetAssayData(seuratObj, assay = assayName, slot = 'counts')

		# Add any new ADTs from this dataset, if needed:
		missing <- rownames(bData)[!(rownames(bData) %in% rownames(assayData))]
		if (length(missing) > 0) {
			missingMat <- matrix(rep(0, ncol(assayData) * length(missing)), ncol = ncol(assayData), nrow = length(missing))
			rownames(missingMat) <- missing
			print(paste0('total ADT rows added to assay: ', length(missing)))
			assayData <- as.sparse(rbind(assayData, missingMat))
		}

		assayData[rownames(bData),colnames(bData), drop = F] <- bData
		assayData <- Seurat::CreateAssayObject(counts = assayData)

		metaToAdd <- GetAssay(seuratObj, assay = assayName)@meta.features
		if (ncol(metaToAdd) > 0) {
			assayData <- Seurat::AddMetaData(object = assayData, metadata = metaToAdd)
		}

		if (sum(colnames(seuratObj) != colnames(assayData)) > 0) {
			stop('The columns of the ADT matrix do not equal the seurat object')
		}

	} else {
		print('no pre-existing assay, creating new')
		#now add empty cells for those lacking ADTs:
		#Append blank cells for any in GEX but missing in ADT:
		missing <- colnames(seuratObj)[!(colnames(seuratObj) %in% colnames(bData))]
		if (length(missing) > 0) {
			missingMat <- matrix(rep(0, nrow(bData) * length(missing)), nrow = nrow(bData), ncol = length(missing))
			colnames(missingMat) <- missing
			print(paste0('total cells lacking data: ', length(missing)))
			bData <- cbind(bData, missingMat)
		}
		bData <- as.sparse(bData[,colnames(seuratObj), drop = F])
		if (sum(colnames(seuratObj) != colnames(bData)) > 0) {
			stop('The columns of the ADT matrix do not equal the seurat object')
		}
		assayData <- Seurat::CreateAssayObject(counts = bData)
	}

	seuratObj[[assayName]] <- assayData

	if (shouldAddMetadata) {
		seuratObj[[assayName]] <- .AddFeatureMetadata(ft, seuratObj[[assayName]])
	}

	if (!skipNormalize) {
		seuratObj <- NormalizeData(seuratObj, assay = assayName, normalization.method = "CLR")
		seuratObj <- ScaleData(seuratObj, assay = assayName)
	}

	return(seuratObj)
}

.AddFeatureMetadata <- function(ft, assayData) {
	print('Adding feature metadata')
	for (colname in names(ft)) {
		print(paste0('adding column: ', colname))
		toAdd <- ft[[colname]]
		names(toAdd) <- ft$markernameUnique
		select <- names(toAdd) %in% rownames(assayData)
		toAdd <- toAdd[select]
		if (length(toAdd) != nrow(ft)) {
			print(paste0('Not all markers in feature table were present in found table, missing: ', paste0(ft$markernameUnique[!select])))
		}

		if (length(toAdd) > 0) {
			assayData <- Seurat::AddMetaData(object = assayData, metadata = toAdd, col.name = colname)
		}
	}

	return(assayData)
}


.PlotCiteSeqCountData <- function(seuratObj, assayName = 'ADT') {
	assayData <- GetAssayData(seuratObj, slot = "counts", assay = assayName)

	featuresToPlot <- rownames(assayData)
	setSize <- 4
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	P <- NULL
	for (i in 0:steps) {
		start <- (i * setSize) + 1
		end <- min((start + setSize - 1), length(featuresToPlot))
		features <- featuresToPlot[start:end]

		P1 <- RidgePlot(seuratObj, assay = assayName, features = features)
		if (is.null(P)) {
			P <- P1
		} else {
			P <- P | P1
		}
	}

	print(P)

	# Also total per ADT
	countsPerAdt <- rowSums(as.matrix(assayData))
	countsPerAdt <- data.frame(Marker = names(countsPerAdt), TotalCount = countsPerAdt)

	P1 <- ggplot(countsPerAdt, aes(x = TotalCount)) +
		geom_density() +
		xlab('Total Count/ADT') +
		ylab('Density') +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P1)

	P2 <- ggplot(countsPerAdt, aes(x = Marker, y = TotalCount)) +
		geom_bar(stat = 'identity') +
		xlab('Marker') +
		ylab('Total Count') +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P2)
}

ProcessCiteSeqData <- function(seuratObj, assayName = 'ADT'){
	origAssay <- DefaultAssay(seuratObj)
	DefaultAssay(seuratObj) <- assayName
	print(paste0('Processing ADT data, features: ', paste0(rownames(seuratObj), collapse = ',')))

	#PCA:
	seuratObj <- RunPCA(seuratObj, features = rownames(seuratObj), reduction.name = "pca_adt", reduction.key = "pcaadt_", verbose = FALSE)
	print(DimPlot(seuratObj, reduction = "pca_adt"))

	adt.data <- GetAssayData(seuratObj, slot = "data")
	adt.dist <- dist(t(adt.data))

	# Before we recluster the data on ADT levels, we'll stash the original cluster IDs for later
	seuratObj[["origClusterID"]] <- Idents(seuratObj)

	# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
	seuratObj[["tsne_adt"]] <- RunTSNE(adt.dist, assay = assayName, reduction.key = "adtTSNE_")
	seuratObj[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
	seuratObj <- FindClusters(seuratObj, resolution = 0.2, graph.name = "adt_snn")

	clustering.table <- table(Idents(seuratObj), seuratObj$origClusterID)
	print(clustering.table)

	#Restore original state:
	Idents(seuratObj) <- seuratObj[["origClusterID"]]
	DefaultAssay(seuratObj) <- origAssay

	#Compare new/old:
	tsne_orig <- DimPlot(seuratObj, reduction = "tsne_adt", group.by = "origClusterID", combine = FALSE) + NoLegend()
	tsne_orig <- tsne_orig + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
	tsne_orig <- LabelClusters(plot = tsne_orig, id = "origClusterID", size = 4)

	tsne_adt <- DimPlot(seuratObj, reduction = "tsne_adt", pt.size = 0.5, combine = FALSE) + NoLegend()
	tsne_adt <- tsne_adt + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
	tsne_adt <- LabelClusters(plot = tsne_adt, id = "ident", size = 4)

	# Note: for this comparison, both the RNA and protein clustering are visualized on a tSNE generated using the ADT distance matrix.
	print(patchwork::wrap_plots(list(tsne_orig, tsne_adt), ncol = 2))

	return(seuratObj)
}