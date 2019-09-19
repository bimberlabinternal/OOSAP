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
	print(paste0('Total cDNA records: ', length(cDndIds)))
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

	gexReadsets <- unique(rows$readsetid)
	print(paste0('total GEX readsets: ', length(gexReadsets)))

	# Identify, download seuratObj, created from the appropriate readsetId:
	seuratRows <- labkey.selectRows(
		baseUrl=lkBaseUrl,
		folderPath=lkDefaultFolder,
		schemaName="sequenceanalysis",
		queryName="outputfiles",
		colSort="-rowid",
		colSelect="rowid,readset",
		colFilter=makeFilter(
			c("readset", "IN", paste0(gexReadsets, collapse = ';')),
			c("category", "EQUAL", "Seurat Data")
		),
		containerFilter=NULL,
		colNameOpt="rname"
	)

	# Possible to have a duplicate:
	seuratRows <- unique(seuratRows)

	if (length(unique(seuratRows$readset)) != length(gexReadsets)) {
		missing <- gexReadsets[!(gexReadsets %in% unique(seuratRows$readset))]
		print(paste0('Not all requested cDNAs have seurat objects.  Readsets missing: ', paste0(unique(missing), collapse = ',')))

		return(NA)
	}

	downloadedFiles <- c()
	ret <- NA
	for (gexReadset in gexReadsets) {
		print(paste0('processing readset: ', gexReadset))
		row <- seuratRows[seuratRows$readset == gexReadset,,drop = F]
		if (nrow(row) > 1) {
			print('More than one seurat row found, using the most recent')
			row <- row[1,,drop = F]
		}

		if (is.na(row[['rowid']])) {
			warning(paste0('Error: RowID was NA, skipping'))
			print(row)
			next
		}

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
			barcodeWhitelist <- colnames(seuratObj)[!seuratObj$Positive.Call & !is.na(seuratObj$cDNA_ID)]
		} else {
			barcodeWhitelist <- colnames(seuratObj)[seuratObj$Positive.Call & !is.na(seuratObj$cDNA_ID)]
		}

		i <- 0
		for (barcodePrefix in unique(unlist(seuratObj[['BarcodePrefix']]))) {
			i <- i + 1

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
			} else if (nrow(bamRows) > 1) {
				print('More than one BAM found')
			}

			alignmentId <- bamRows$dataid[1]

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

			#summarize metadata
			meta <- data.frame(
				barcode = colnames(seuratObj)[seuratObj$BarcodePrefix == barcodePrefix],
				SubjectId = as.character(seuratObj$SubjectId[seuratObj$BarcodePrefix == barcodePrefix]),
				Stim = as.character(seuratObj$Stim[seuratObj$BarcodePrefix == barcodePrefix]),
				population = as.character(seuratObj$Population[seuratObj$BarcodePrefix == barcodePrefix]),
				date = as.character(seuratObj$SampleDate[seuratObj$BarcodePrefix == barcodePrefix]),
				cdna = as.character(seuratObj$cDNA_ID[seuratObj$BarcodePrefix == barcodePrefix]),
				libraryId = c(libraryId),
				alignmentId = c(alignmentId),
				analysisId = c(analysisId)
			)
			meta$SampleName <- paste0(meta$SubjectId, '_', meta$Stim)
			tcrData <- merge(tcrData, meta, by = c('barcode'), all.x = T)

			if (sum(is.na(tcrData$SubjectId)) > 0) {
				f <- paste0(outPrefix, 'temp.txt')
				write.table(tcrData, file = f, sep = '\t', quote = F, row.names = F)
				stop(paste0('Missing subject IDs!  See ', f, ' for table of results'))
			}

			# Group
			tcrData <- tcrData %>% group_by(SampleName, SubjectId, population, date, cdna, libraryId, alignmentId, analysisId, chain, cdr3, v_gene, d_gene, j_gene, c_gene, raw_clonotype_id, raw_consensus_id) %>% summarize(count = dplyr::n())
			names(tcrData) <- c('SampleName', 'SubjectId', 'population', 'date', 'cdna', 'libraryId', 'alignmentId', 'analysisId', 'locus', 'cdr3', 'vHit', 'dHit', 'jHit', 'cHit', 'cloneId', 'consensus_id', 'count')

			tcrData <- tcrData %>% group_by(cdna) %>% mutate(totalCells = dplyr::n())
			tcrData$fraction = tcrData$count / tcrData$totalCells

			# Merge sequence:
			fastaData <- Biostrings::readDNAStringSet(fastaFile)
			seqDf <- data.frame(consensus_id = names(fastaData), sequence = paste(fastaData))

			tcrData <- merge(tcrData, seqDf, by = c('consensus_id'), all.x = T)
			tcrData <- tcrData[!(names(tcrData) %in% c('consensus_id', 'totalCells', 'Stim'))]
			tcrData[tcrData == 'None'] <- NA

			if (all(is.na(ret))) {
				ret <- tcrData
			} else {
				ret <- rbind(ret, tcrData)
			}
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
#' @param minCells If provided, data will only be imported if at least this many cells exist for the cDNA library
#' @export
#' @import Seurat
#' @importFrom dplyr %>% group_by select n summarize
CalculateTCRFreqForActivatedCellsAndImport <- function(cDndIds, workbook = NULL, geneSetName = 'HighlyActivated', positivityThreshold = 0.5, outPrefix = './', assayName = 'TCRdb', populationNameSuffix = '-HA', invert = FALSE, doCleanup = FALSE, minCells = NULL) {
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
		print(paste0('Preparing run for analysis Id: ',analysisId,', total rows: ', nrow(df), ', total cells: ', sum(df$count)))
		if (!is.null(minCells)) {
			totals <- df %>% group_by(cdna) %>% summarize(total = sum(count))
			toKeep <- unique(totals$cdna[totals$total >= minCells])
			if (length(toKeep) != length(unique(df$cdna))) {
				print('The following cDNA libraries will be skipped due to low cells:')
				print(paste0(totals$cdna[totals$total < minCells], ': ', totals$total[totals$total < minCells]))

				df <- df[df$cdna %in% toKeep,]
			}
		}

		if (sum(is.na(df$SubjectId)) > 0) {
			f <- paste0(outPrefix, 'temp.txt')
			write.table(df, file = f, sep = '\t', quote = F, row.names = F)
			stop(paste0('Missing subject IDs!  See ', f, ' for table of results'))
		}

		if (nrow(df) > 0) {
			print('run passed, adding')
			run <- labkey.experiment.createRun(list(name = paste0('AnalysisId: ', analysisId, populationNameSuffix), properties = list(assayName = '10x', analysisId = analysisId)), dataRows = df)
			runList <- append(runList, list(run))
		} else {
			print('no rows passed, skipping')
		}
	}

	print(paste0('total runs: ', length(runList)))
	if (length(runList) > 0) {
		labkey.experiment.saveBatch(
			baseUrl=lkBaseUrl,
			folderPath=folder,
			assayConfig=list(assayName=assayName, providerName='TCRdb'),
			runList = runList
		)
	} else {
		print('no runs to save')
	}
}