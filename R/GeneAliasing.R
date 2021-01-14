
# Perform the actual query against STRINGdb
.QuerySTRINGdb <- function(inputIds, speciesId, score_threshold = 0, stringDBVersion = "11"){
	print('Querying STRINGdb')
  string_db <- STRINGdb::STRINGdb$new(version = stringDBVersion,
                            species = speciesId, 
                            score_threshold = score_threshold, 
                            input_directory = "")
  
  ## map inputIds to stringID mapIds
  inputIds.map <- string_db$map(my_data_frame = data.frame(InputTerm = as.character(inputIds)),
                                my_data_frame_id_col_names = "InputTerm",
                                takeFirst = T, 
                                removeUnmappedRows = FALSE)
	## get all species aliases
  stringdb.alias <- string_db$get_aliases()
  stringdb.alias <- stringdb.alias %>%
		dplyr::group_by(STRING_id) %>%
		dplyr::summarize(STRING.aliases = paste0(sort(unique(alias)), collapse = ','))

  stringdb.alias <- merge(inputIds.map, stringdb.alias, by = c("STRING_id"), all.x = F)

	stringdb.alias <- stringdb.alias %>%
		dplyr::group_by(InputTerm) %>%
		dplyr::summarize(STRING_id = paste0(sort(unique(STRING_id)), collapse = ','), STRING.aliases = paste0(sort(unique(STRING.aliases)), collapse = ','))

	print(paste0('Found ', sum(!is.na(stringdb.alias$STRING_id)), ' of ', length(inputIds)))

  return(stringdb.alias)
}


# Perform the actual query against DAVID
.QueryDAVID <- function(inputIds, email, idType){
	print(paste0('Querying DAVID by: ', idType))
  david <- RDAVIDWebService::DAVIDWebService$new(email = email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  
  david.addList <- RDAVIDWebService::addList(object = david,
                           inputIds = inputIds, 
                           listName = 'MyListName',
                           idType = idType, 
                           listType = 'Gene')
	print(paste0('Found ', (100*david.addList$inDavid), '% of ', length(inputIds)))

  davidReport <- david$getGeneListReport()
  davidResult.df <- data.frame(InputIds = as.character(davidReport$ID), DAVID.GeneName = as.character(davidReport$Name), DAVID.Id = as.character(davidReport$ID), stringsAsFactors = FALSE)

	# Extract gene symbol from name string (i.e. "microRNA 214(MIR214)")
  davidResult.df$DAVID.Symbol <- stringr::str_match(davidResult.df$DAVID.GeneName, '\\(.*?\\)$')
  davidResult.df$DAVID.Symbol <- sapply(davidResult.df$DAVID.Symbol, function(x){
		return(unlist(strsplit(x, '[()]'))[2])
	})

  return(davidResult.df)
}


# Basic argument checking
.CheckGeneInputs <- function(ensemblIds, geneSymbols){
  if (.IsEmpty(ensemblIds) && .IsEmpty(geneSymbols)) {
    stop('Must provide either ensemblIds or geneSymbols')
  }
  else if (!.IsEmpty(ensemblIds) & !.IsEmpty(geneSymbols)) {
    if (length(ensemblIds) != length(geneSymbols)) {
      stop('EnsemblIds and geneSymbol must be of equal length.')
    }
  }
}


# Utility function to test if the input is NA or NULL
.IsEmpty <- function(x) {
	return(all(is.na(x)) || all(is.null(x)))  
}


#' @title TranslateToEnsembl
#' @param ensemblIds A vector of ensembl IDs, passed to the biomaRt::getBM() to query against ensembl_gene_id
#' @param geneSymbols A vector of gene symbols, passed to the biomaRt::getBM() to query against hgnc_symbol
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly to biomaRt::useEnsembl
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr %>% group_by summarize mutate coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export
TranslateToEnsembl <- function(ensemblIds = NULL, geneSymbols = NULL, dataset = "mmulatta_gene_ensembl", ensemblVersion = NULL, ensemblMirror = NULL, biomart = "ensembl"){
  .CheckGeneInputs(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
  if (is.null(ensemblIds)) {
    ensemblIds <- NA
  }

  if (is.null(geneSymbols)) {
    geneSymbols <- NA
  }
  
	combinedEnsembl <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	combinedEnsembl$Order <- 1:nrow(combinedEnsembl)

	ensemblById <- NA
	if (!.IsEmpty(ensemblIds)) {
		ensemblById <- .QueryEnsembl(inputIds = ensemblIds,
                             queryField = "ensembl_gene_id",
                             dataset = dataset, 
                             ensemblVersion = ensemblVersion,
														 ensemblMirror = ensemblMirror,
                             biomart = biomart)

		ensemblById$EnsemblId <- ensemblById$ensembl_gene_id
		ensemblById <- merge(combinedEnsembl, ensemblById, by = 'EnsemblId', all.x = T)
		ensemblById <- dplyr::arrange(ensemblById, Order)
  }

	ensemblBySymbol1 <- NA
	ensemblBySymbol2 <- NA
  if (!.IsEmpty(geneSymbols)) {
		ensemblBySymbol1 <- .QueryEnsembl(inputIds = geneSymbols,
                             queryField = "external_gene_name",
                             dataset = dataset, 
                             ensemblVersion = ensemblVersion,
														 ensemblMirror = ensemblMirror,
                             biomart = biomart)
		ensemblBySymbol1$GeneSymbol <- ensemblBySymbol1$external_gene_name
		ensemblBySymbol1 <- merge(combinedEnsembl, ensemblBySymbol1, by = 'GeneSymbol', all.x = T)
		ensemblBySymbol1 <- dplyr::arrange(ensemblBySymbol1, Order)

		ensemblBySymbol2 <- .QueryEnsembl(inputIds = geneSymbols,
														queryField = "hgnc_symbol",
														dataset = dataset,
														ensemblVersion = ensemblVersion,
														ensemblMirror = ensemblMirror,
														biomart = biomart)
		ensemblBySymbol2$GeneSymbol <- ensemblBySymbol2$hgnc_symbol
		ensemblBySymbol2 <- merge(combinedEnsembl, ensemblBySymbol2, by = 'GeneSymbol', all.x = T)
		ensemblBySymbol2 <- dplyr::arrange(ensemblBySymbol2, Order)
	}

	#Concat in preferential order, based on EnsemblId.  We can assume any resolved hit from Ensembl will have an Ensembl ID
	ret <- data.frame(EnsemblId = combinedEnsembl$EnsemblId, GeneSymbol = combinedEnsembl$GeneSymbol, ensembl_gene_id = NA, hgnc_symbol = NA, external_gene_name = NA, Order = 1:nrow(combinedEnsembl), stringsAsFactors=FALSE)
	ret <- .ConcatPreferentially(fieldToTest = 'ensembl_gene_id', datasets = list(ensemblById, ensemblBySymbol1, ensemblBySymbol2), baseDf = ret)
	ret <- ret[names(ret) != 'Order']

	return(ret)
}

.ConcatPreferentially <- function(fieldToTest, datasets, baseDf){
	datasets <- datasets[!is.na(datasets)]
	if (!('Order' %in% names(baseDf))) {
		baseDf$Order <- 1:nrow(baseDf)
	}

	ret <- baseDf[FALSE, TRUE]  #zero rows

	for (dataset in datasets) {
		toAppend <- dataset[!is.na(dataset[[fieldToTest]]) & !(dataset[[fieldToTest]] %in% ret[[fieldToTest]]),]
		ret <- rbind(ret, toAppend[names(ret)])
	}

	#now ensure all input terms represented, in order:
	ret <- rbind(ret, baseDf[!(baseDf$Order %in% ret$Order),])
	ret <- dplyr::arrange(ret, Order)

  return(ret)
}


# Translate a list of gene IDs to Ensembl IDs, based on a target field
.QueryEnsembl <- function(inputIds, queryField, biomart, dataset, ensemblVersion, ensemblMirror){
	print(paste0('Querying Ensembl using: ', queryField))
	ensembl <- biomaRt::useEnsembl(biomart = biomart,
		dataset = dataset,
		version = ensemblVersion,
		mirror = ensemblMirror
	)

	ensemblResults <- NULL
	tryCatch(expr = {
		ensemblResults <- biomaRt::getBM(
			attributes = c('ensembl_gene_id', 'hgnc_symbol', 'external_gene_name'),
			filters = c(queryField),
			values = inputIds,
			mart = ensembl
		)
	}, error = function(e){
		print(e)
		stop(paste0('Error querying ensembl, using genes: ', paste0(inputIds, collapse = ';')))
	})

	#Drop duplicates.  Note: could consider group_concat on variables?
	ensemblResults <- ensemblResults %>% group_by_at(queryField) %>% dplyr::mutate(total = dplyr::n())
	ensemblResults <- ensemblResults[ensemblResults$total == 1,]
	ensemblResults <- ensemblResults[names(ensemblResults) != 'total']

	print(paste0('Found ', sum(!is.na(ensemblResults$ensembl_gene_id)), ' of ', length(inputIds)))

	return(ensemblResults)
}

  
#' @title TranslateToStringDb
#' @param ensemblIds A vector of ensembl IDs
#' @param geneSymbols A vector of gene symbols
#' @param speciesId Species ID. see Stringdb reference for list of avialable species
#' @param replaceUnmatched Logical. If TRUE, removes NA's and replaces with inputs
#' @import STRINGdb 
#' @importFrom dplyr %>% group_by summarize mutate coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

TranslateToStringDb <- function(ensemblIds = NULL, geneSymbols = NULL, speciesId = 9606){
	.CheckGeneInputs(ensemblIds, geneSymbols)
  if (is.null(ensemblIds)) {
    ensemblIds <- NA
  }
  
  if (is.null(geneSymbols)) {
    geneSymbols <- NA
  }
  
  queryIds <- character()
	if (!.IsEmpty(ensemblIds)) {
    queryIds <- c(queryIds, ensemblIds)
  }

	if (!.IsEmpty(geneSymbols)) {
    queryIds <- c(queryIds, geneSymbols)
  }
  
  queryIds <- unique(queryIds)
  stringRes <- .QuerySTRINGdb(inputIds = queryIds, speciesId = speciesId)

	#Preferentially accept results by EnsemblId
	resultsBase <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	resultsBase$Order <- 1:nrow(resultsBase)

	resultsById <- merge(resultsBase, stringRes, by.x = 'EnsemblId', all.x = T, by.y = 'InputTerm')
	resultsById <- resultsById[!is.na(resultsById$STRING_id),]

	resultsBySymbol <- merge(resultsBase, stringRes, by.x = 'GeneSymbol', all.x = T, by.y = 'InputTerm')
	resultsBySymbol <- resultsBySymbol[!(resultsBySymbol$Order %in% resultsById$Order),]

	results <- rbind(resultsById, resultsBySymbol)
	results <- rbind(results, resultsBase[!(resultsBase$Order %in% results$Order),])
	results <- dplyr::arrange(results, Order)
	results <- results[names(results) != 'Order']

	return(results)
}


#' @title TranslateToDAVID
#' @param ensemblIds A vector of ensembl IDs
#' @param geneSymbols A vector of gene symbols
#' @param email email account registered with DAVIDWebService 
#' @rawNamespace import(RDAVIDWebService, except = c('counts', 'cluster'))
#' @importFrom dplyr %>% group_by summarize mutate coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

TranslateToDAVID <- function(ensemblIds = NULL, geneSymbols = NULL, email){
	.CheckGeneInputs(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
  if (is.null(ensemblIds)) {
    ensemblIds <- NA
  }
  
  if (is.null(geneSymbols)) {
    geneSymbols <- NA
  }

	results <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	results$Order <- 1:nrow(results)

	davidResById <- NA
	if (!.IsEmpty(ensemblIds)) {
    davidResById <- .QueryDAVID(inputIds = ensemblIds, email = email, idType = 'ENSEMBL_GENE_ID')
		davidResById <- merge(results, davidResById, by.x = 'EnsemblId', all.x = F, by.y = 'InputIds')
  }

	davidResBySymbol <- NA
	if (!.IsEmpty(geneSymbols)) {
		# OFFICIAL_GENE_SYMBOL is not available in RDAVIDWebService for now
    #davidResBySymbol <- .QueryDAVID(inputIds = geneSymbols, email = email, idType = 'OFFICIAL_GENE_SYMBOL')
		#davidResBySymbol <- merge(results, davidResById, by.x = 'GeneSymbol', all.x = F, by.y = 'InputIds')
  }

	ret <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, DAVID.Id = NA, DAVID.GeneName = NA, DAVID.Symbol = NA, stringsAsFactors=FALSE)
	ret$Order <- 1:nrow(ret)
	ret <- .ConcatPreferentially(fieldToTest = 'DAVID.Id', datasets = list(davidResById, davidResBySymbol), baseDf = ret)
	ret <- ret[names(ret) != 'Order']

  return(ret)
}


#' @title aliasTable
#' @param ensemblIds A vector of ensembl IDs
#' @param geneSymbols A vector of gene symbols
#' @param ensemblAttributes A vector of ensembl attributes
#' @param ensemblDataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly to biomaRt::useEnsembl
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param davidEmail email account registered with DAVIDWebService 
#' @param stringSpeciesId species ID. see Stringdb for list of available species
#' @param aliasPriorityOrder vector containig priority order of alias database. Must be UPPERCASE. Current databases: ENSEMBL, STRING, DAVID
#' @importFrom biomaRt useEnsembl getBM
#' @import STRINGdb
#' @rawNamespace import(RDAVIDWebService, except = c('counts', 'cluster'))
#' @importFrom dplyr %>% group_by summarize mutate coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

TranslateGeneNames <- function(ensemblIds = NULL, geneSymbols = NULL, davidEmail,
                       ensemblDataset = "mmulatta_gene_ensembl", ensemblVersion = NULL, ensemblMirror = NULL, biomart = 'ensembl',
                       stringSpeciesId = 9606,
											 useEnsembl = TRUE, useSTRINGdb = TRUE, useDAVID = TRUE){

	.CheckGeneInputs(ensemblIds, geneSymbols)
	if (is.null(ensemblIds)) {
		ensemblIds <- NA
	}

	if (is.null(geneSymbols)) {
		geneSymbols <- NA
	}

	inputDf <- data.frame(EnsemblId = ensemblIds, GeneSymbol = geneSymbols, stringsAsFactors = FALSE)
	inputDf$Order <- 1:nrow(inputDf)

	if (useEnsembl) {
	  ret.ensembl <- TranslateToEnsembl(ensemblIds = ensemblIds,
                              geneSymbols = geneSymbols, 
                              dataset = ensemblDataset,
                              ensemblVersion = ensemblVersion,
															ensemblMirror = ensemblMirror,
                              biomart = biomart)

		if (nrow(inputDf) != nrow(ret.ensembl)) {
			stop('Rows not equal for ensembl result')
		}

		ret.ensembl$Order <- 1:nrow(ret.ensembl)
		ret.ensembl <- ret.ensembl[!(names(ret.ensembl) %in% c('EnsemblId', 'GeneSymbol'))]
		inputDf <- merge(inputDf, ret.ensembl, by = 'Order', all.x = TRUE)
	}

	if (useSTRINGdb) {
  	ret.string <- TranslateToStringDb(ensemblIds = ensemblIds,
                              geneSymbols = geneSymbols,
                              speciesId = stringSpeciesId)

		if (nrow(inputDf) != nrow(ret.string)) {
			stop('Rows not equal for STRINGdb result')
		}

		ret.string$Order <- 1:nrow(ret.string)
		ret.string <- ret.string[!(names(ret.string) %in% c('EnsemblId', 'GeneSymbol'))]
		inputDf <- merge(inputDf, ret.string, by = 'Order', all.x = TRUE)
	}

	if (useDAVID) {
  	ret.david <- TranslateToDAVID(ensemblIds = ensemblIds, geneSymbols = geneSymbols, email = davidEmail)
		if (nrow(inputDf) != nrow(ret.david)) {
			stop('Rows not equal for DAVID result')
		}

		ret.david$Order <- 1:nrow(ret.david)
		ret.david <- ret.david[!(names(ret.david) %in% c('EnsemblId', 'GeneSymbol'))]
		inputDf <- merge(inputDf, ret.david, by = 'Order', all.x = TRUE)
	}

	inputDf <- dplyr::arrange(inputDf, Order)
	inputDf <- inputDf[names(inputDf) != 'Order']

  return(inputDf)
}


