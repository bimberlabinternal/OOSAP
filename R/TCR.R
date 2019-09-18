#' @include LabKeySettings.R
#' @include Utils.R



#' @title CalculateTCRFreqForActivatedCells
#' @description For the supplied cDNA Rows, this will query their gene expression and TCR readsets, identifying a) existing seurat objects, b) vloupe files. If both are found,
#  it will download them, build a whitelist of highly activated cells (using SGS), calculate TCR frequencies, and return a dataframe
#' @return A dataframe of results
#' @param cDndIds A vector of cDNA rowIDs
#' @param geneSetName The gene set name to use for SGS
#' @param positivityThreshold The threshold to use for calling cells as positive
#' @param outPrefix A string that will be prepended to all saved files
#' @param invert If TRUE, those cells NOT positive for the gene set will be summarized, instead of positive cells
#' @param doCleanup If TRUE, any downloaded files will be deleted on completion
#' @export
#' @import Seurat
#' @importFrom Biostrings readDNAStringSet
#' @importFrom dplyr %>% group_by select n summarize
CalculateTCRFreqForActivatedCells <- function(cDndIds, geneSetName = 'HighlyActivated', positivityThreshold = 0.5, outPrefix = './', invert = FALSE, doCleanup = FALSE) {
	rows <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=lkDefaultFolder,
		schemaName="tcrdb",
		queryName="cdnas",
		viewName="",
		colSort="-rowid",
		colFilter = makeFilter(c("rowid", "IN", paste0(cDndIds, collapse = ";"))),
		colSelect="rowid,readsetid,enrichedreadsetid,hashingreadsetid",
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(rows) != length(cDndIds)) {
		print(paste0('Not all requested cDNAs found.  Row IDs found: ', paste0(unique(rows$rowid), collapse = ',')))
		return(NA)
	}

	# Identify, download seuratObj, created from the appropriate readsetId:
	seuratRows <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=lkDefaultFolder,
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		colSort="-rowid",
		colSelect="rowid",
		colFilter=makeFilter(
			c("readset", "IN", paste0(unique(rows$readsetid), collapse = ';')),
			c("category", "EQUAL", "Seurat Data")
		),
		containerFilter=NULL,
		colNameOpt="rname"
	)

	if (nrow(rows) != length(cDndIds)) {
		print(paste0('Not all requested cDNAs found.  Row IDs found: ', paste0(unique(rows$rowid), collapse = ',')))
		return(NA)
	}

	downloadedFiles <- c()
	ret <- NA
	for (idx in 1:nrow(seuratRows)) {
		row <- seuratRows[idx,,drop = F]

		f <- paste0(outPrefix, row[['rowid']], '.seurat.rds')
		DownloadOutputFile(row[['rowid']], f, overwrite = F)
		downloadedFiles <- c(downloadedFiles, f)

		# For each, apply metadata, TCR clones
		seuratObj <- readRDS(f)

		seuratObj <- DownloadAndAppendCellHashing(seuratObject = seuratObj)
		seuratObj <- QueryAndApplyCdnaMetadata(seuratObj)
		seuratObj <- ClassifySGSAndApply(seuratObj = seuratObj, geneSetName = 'Positive', geneList = OOSAP::Phenotyping_GeneList()[[geneSetName]], positivityThreshold = positivityThreshold)
		if (invert) {
			print('Selecting cells without the provided signature')
			barcodeWhitelist <- colnames(seuratObj)[!seuratObj$Positive.Call]
		} else {
			barcodeWhitelist <- colnames(seuratObj)[seuratObj$Positive.Call]
		}

		i <- 0
		for (barcodePrefix in unique(unique(unlist(seuratObj[['BarcodePrefix']])))) {
			i <- i + 1

			cDNA <- unique(seuratObj$cDNA_ID[seuratObj$BarcodePrefix == barcodePrefix])
			if (length(cDNA) > 1) {
				stop(paste0('More than one cDNA ID found for barcodePrefix: ', barcodePrefix))
			}

			vloupeId <- .FindMatchedVloupe(barcodePrefix)
			if (is.na(vloupeId)){
				stop(paste0('Unable to find VLoupe file for loupe file: ', barcodePrefix))
			}

			#TCR libraryId
			tcrLibRows <- labkey.selectRows(
				baseUrl=lkBaseUrl,
				folderPath=lkDefaultFolder,
				schemaName="sequenceanalysis",
				queryName="outputfiles",
				colSort="-rowid",
				colSelect="library_id,analysis_id,readset",
				colFilter=makeFilter(c("rowid", "EQUALS", vloupeId)),
				containerFilter=NULL,
				colNameOpt="rname"
			)
			libraryId <- tcrLibRows$library_id[1]
			analysisId <- tcrLibRows$analysis_id[1]
			tcrReadset <- tcrLibRows$readset[1]

			bamRows <- labkey.selectRows(
				baseUrl=lkBaseUrl,
				folderPath=lkDefaultFolder,
				schemaName="sequenceanalysis",
				queryName="outputfiles",
				colSort="-rowid",
				colSelect="dataid",
				colFilter=makeFilter(
					c("library_id", "EQUALS", libraryId),
					c("readset", "EQUALS", tcrReadset),
					c("category", "EQUALS", "Alignment")
				),
				containerFilter=NULL,
				colNameOpt="rname"
			)
			if (nrow(bamRows) == 0) {
				stop('Unable to find alignment')
			}

			alignmentId <- bamRows[1][['dataid']]

			#All clonotypes
			clonotypeFile <- file.path(outPrefix, paste0(barcodePrefix, '_all_contig_annotations.csv'))
			.DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = clonotypeFile, overwrite = T)
			if (!file.exists(clonotypeFile)){
				stop(paste0('Unable to download clonotype file for prefix: ', barcodePrefix))
			}
			downloadedFiles <- c(downloadedFiles, clonotypeFile)

			#FASTA:
			fastaFile <- file.path(outPrefix, paste0(barcodePrefix, '_consensus.fasta'))
			.DownloadCellRangerClonotypes(vLoupeId = vloupeId, outFile = fastaFile, overwrite = T, fileName = 'consensus.fasta')
			if (!file.exists(fastaFile)){
				stop(paste0('Unable to download clonotype FASTA for prefix: ', barcodePrefix))
			}
			downloadedFiles <- c(downloadedFiles, fastaFile)

			tcrData <- .ProcessTcrClonotypes(clonotypeFile)
			if (!is.null(barcodePrefix)){
				tcrData$barcode <- as.character(tcrData$barcode)
				tcrData$barcode <- paste0(barcodePrefix, '_', tcrData$barcode)
				tcrData$barcode <- as.factor(tcrData$barcode)
			}

			retain <- intersect(barcodeWhitelist, tcrData$barcode)

			print(paste0('initial barcodes with TCR call: ', nrow(tcrData)))
			pct1 <- round(length(retain) / length(unique(tcrData$barcode)), 2)
			pct2 <- round(length(retain) / length(barcodeWhitelist), 2)

			print(paste0('overlapping with activated cells: ', length(retain), ' (', pct1, ' of TCR calls, ',pct2,' of positive cells)'))
			tcrData <- tcrData[tcrData$barcode %in% retain,]
			if (nrow(tcrData) == 0) {
				print(paste0('no rows for prefix: ', barcodePrefix))
				next
			}

			# Group
			tcrData$barcodePrefix <- barcodePrefix
			tcrData <- tcrData %>% group_by(barcodePrefix, chain, cdr3, v_gene, d_gene, j_gene, c_gene, raw_clonotype_id, raw_consensus_id) %>% summarize(count = dplyr::n())
			names(tcrData) <- c('barcodePrefix', 'locus', 'cdr3', 'vHit', 'dHit', 'jHit', 'cHit', 'cloneId', 'consensus_id', 'count')

			tcrData <- tcrData %>% group_by(barcodePrefix) %>% mutate(totalCells = dplyr::n())
			tcrData$fraction = tcrData$count / tcrData$totalCells

			#summarize metadata
			meta <- data.frame(
				barcodePrefix = barcodePrefix,
				SubjectId = as.character(seuratObj$SubjectId[seuratObj$BarcodePrefix == barcodePrefix]),
				Stim = as.character(seuratObj$Stim[seuratObj$BarcodePrefix == barcodePrefix]),
				population = as.character(seuratObj$Population[seuratObj$BarcodePrefix == barcodePrefix]),
				date = as.character(seuratObj$SampleDate[seuratObj$BarcodePrefix == barcodePrefix]),
				cdna = as.character(seuratObj$cDNA_ID[seuratObj$BarcodePrefix == barcodePrefix]),
				libraryId = c(libraryId),
				alignmentId = c(alignmentId),
				analysisId = c(analysisId)
			)

			meta <- unique(meta)
			meta$SampleName <- paste0(meta$SubjectId, '_', meta$Stim)
			tcrData <- merge(tcrData, meta, by = c('barcodePrefix'), all.x = T)

			# Merge sequence:
			fastaData <- Biostrings::readDNAStringSet(fastaFile)
			seqDf <- data.frame(consensus_id = names(fastaData), sequence = paste(fastaData))

			tcrData <- merge(tcrData, seqDf, by = c('consensus_id'), all.x = T)
			tcrData <- tcrData[!(names(tcrData) %in% c('consensus_id', 'totalCells', 'Stim', 'barcodePrefix'))]
			tcrData[tcrData == 'None'] <- NA

			if (all(is.na(ret))) {
				ret <- tcrData
			} else {
				ret <- rbind(ret, tcrData)
			}
		}

		if (doCleanup) {
			print('Cleaning up downloaded files')
			for (f in downloadedFiles) {
				unlink(f)
			}
		}
		return(ret)
	}
}

#' @title CalculateTCRFreqForActivatedCellsAndImport
#' @description For the supplied cDNA Rows, this will call CalculateTCRFreqForActivatedCells(), and also import the results into the provided assay
#' @return Returns the object representation of the experiment batch (from Rlabkey::labkey.experiment.saveBatch)
#' @param cDndIds A vector of cDNA rowIDs
#' @param geneSetName The gene set name to use for SGS
#' @param positivityThreshold The threshold to use for calling cells as positive
#' @param outPrefix A string that will be prepended to all saved files
#' @param workbook The target workbook to save results
#' @param assayName The name of the target assay
#' @param populationNameSuffix This string will be added to the end of the population field on the saved assay rows, appended to the associated cDNA population
#' @param invert If TRUE, those cells NOT positive for the gene set will be summarized, instead of positive cells
#' @param doCleanup If TRUE, any downloaded files will be deleted on completion
#' @export
#' @import Seurat
#' @importFrom dplyr %>% group_by select n summarize
CalculateTCRFreqForActivatedCellsAndImport <- function(cDndIds, workbook = NULL, geneSetName = 'HighlyActivated', positivityThreshold = 0.5, outPrefix = './', assayName = 'TCRdb', populationNameSuffix = '-HA', invert = FALSE, doCleanup = FALSE) {
	folder <- lkDefaultFolder
	if (!is.null(workbook)) {
		folder <- paste0(folder, workbook)
	}

	resultDataFrame <- CalculateTCRFreqForActivatedCells(cDndIds, geneSetName = geneSetName, positivityThreshold = positivityThreshold, outPrefix = outPrefix, invert = invert)
	resultDataFrame$calculatedPopulation <- paste0(as.character(resultDataFrame$population), populationNameSuffix)

	analysisIds <- unique(resultDataFrame$analysisId)
	print(paste0('Total analyses: ', length(analysisIds)))

	runList <- list()
	for (analysisId in analysisIds) {
		df <- resultDataFrame[resultDataFrame$analysisId == analysisId,]
		print(paste0('Preparing run for analysis Id: ',analysisId,', total rows: ', nrow(df)))
		run <- labkey.experiment.createRun(list(name = paste0('AnalysisId: ', analysisId), properties = list(assayName = '10x', analysisId = analysisId)), dataRows = df)
		runList <- append(runList, list(run))
	}

	print(paste0('total runs: ', length(runList)))

	labkey.experiment.saveBatch(
		baseUrl=lkBaseUrl,
		folderPath=folder,
		assayConfig=list(assayName=assayName, providerName='TCRdb'),
		runList = runList
	)
}