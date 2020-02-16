####################################################################################################################
mapENSEMBL <- function(inputIds, biomart, dataset, version, mirror, attributes, filters){
  ensembl <- useEnsembl(biomart = biomart, 
                        dataset = dataset, 
                        version = version, 
                        mirror = mirror)
  ensemblResults <- getBM(attributes = attributes,
                          filters = filters, 
                          values = inputIds, 
                          mart = ensembl)
  return(ensemblResults)
}

####################################################################################################################
mapSTRINGdb <- function(inputIds, speciesId, score_threshold = 0, version = "10"){
  string_db <- STRINGdb$new(version = version, 
                            species = speciesId, 
                            score_threshold = score_threshold, 
                            input_directory = "")
  
  ## map inputIds to stringID mapIds
  inputIds.df <- data.frame(InputIds = as.character(inputIds))
  inputIds.map <- string_db$map(my_data_frame = inputIds.df, 
                                my_data_frame_id_col_names = "InputIds", 
                                takeFirst = T, 
                                removeUnmappedRows = FALSE)
  
  ## get all species aliases
  stringdb.alias <- string_db$get_aliases()
  stringdb.alias <- stringdb.alias %>% 
    group_by(STRING_id) %>% 
    summarize(STRING.aliases = paste0(alias, collapse = ','))
  
  ## megre alias table with mapped inputIds to get inputIds' aliases
  stringdbResults <- merge(stringdb.alias, inputIds.map, by = c("STRING_id"), all.x = F)
  ## get 1st alias(gene name) from list of aliases 
  stringdbResults$STRING.symbol <- sapply(strsplit(stringdbResults$STRING.aliases, ','), head, 1)
  
  return(stringdbResults)
}

####################################################################################################################
mapDAVID <- function(inputIds, email, idType, listType = "Gene", listName = "MyListName"){
  david <- DAVIDWebService$new(email = email, 
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  
  david.addList <- addList(object = david, 
                           inputIds = inputIds, 
                           listName = listName, 
                           idType = idType, 
                           listType = listType)
  
  davidReport <- david$getGeneListReport()
  
  davidResult.df <- data.frame(InputIds = as.character(davidReport$ID), DAVID.gene_name = as.character(davidReport$Name))
  
  davidResult.df$DAVID.symbol <- str_match(davidResult.df$DAVID.gene_name, '\\(.*?\\)$')
  davidResult.df$DAVID.symbol <- lapply(strsplit(davidResult.df$DAVID.symbol, '[()]'), tail, n = 1L)
  return(davidResult.df)
}

####################################################################################################################
inputOrder <- function(input, mergeWith, inputLen){
  ##add newly translated genes to main gene table 
  if (is.null(input)){
    input.df <- data.frame(InputIds = NA, SortOrder = as.numeric(1:inputLen))
  } else{
    input.df <- data.frame(InputIds = as.character(input), SortOrder = as.numeric(1:inputLen))
  }
  
  merged <- merge(input.df, mergeWith, by = c("InputIds"), all.x = T)
  merged <- merged[order(merged$SortOrder),]
  
}

####################################################################################################################
coalesceAlias <- function(ensInput, symbInput, toCoalesce_x, toCoalesce_y, replaceUnmatched = FALSE){
  ##to coalesce symbols outputted from ensemblId input and geneSymbols input
  merged <- merge(ensInput, symbInput, by = c("SortOrder"))
  
  merged <- mutate_all(merged, as.character)
  if (nrow(merged) > 0){
    merged[merged == "" | merged == "NA"] <- NA
  }
  
  if (replaceUnmatched == TRUE){
    coalesced <- coalesce(merged[[toCoalesce_x]], merged[[toCoalesce_y]], merged[["InputIds.y"]], merged[["InputIds.x"]])
  } else{
    coalesced <- coalesce(merged[[toCoalesce_x]], merged[[toCoalesce_y]])
  }
  
  coalesced <- data.frame(SortOrder = as.character(merged[["SortOrder"]]), 
                          InputIds.EnsemblId = as.character(merged[["InputIds.x"]]), 
                          InputIds.GeneSymbol = as.character(merged[["InputIds.y"]]), 
                          Coalesced = as.vector(coalesced))
  
  return(coalesced)
}

####################################################################################################################
checkInputLen <- function(ensemblIds, geneSymbols){
  lenEns <- length(ensemblIds)
  lenSym <- length(geneSymbols)
  
  inputLen <- max(lenEns, lenSym)
  
  if (!lenEns == lenSym){
    if (!is.null(ensemblIds) & !is.null(geneSymbols)){
      stop("Unequal length of inputs. If ensemblIds and geneSymbol is non-null, both must be off equal length. ")
    }
  }
  return(inputLen)
}

####################################################################################################################
#' @title aliasENSEMBL
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param geneSymbols A vector of gene symbols, passed to the getBM() filters argument
#' @param attributes A vector of ensembl attributes
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @param version Passed directly to biomaRt::useEnsembl
#' @param mirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param replaceUnmatched Logical. If TRUE, removes NA's and replaces with inputs
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

aliasENSEMBL <- function(ensemblIds = NULL, geneSymbols = NULL, 
                         attributes = c('hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'),  
                         dataset = "mmulatta_gene_ensembl", version = NULL, mirror = "uswest", biomart = "ensembl", 
                         replaceUnmatched = TRUE){
  
  #check ensemblIds, geneSymbols for NA/NULL.
  # if both non-null, make sure they are the same length
  inputLen <- checkInputLen(ensemblIds = ensemblIds, geneSymbols = geneSymbols)

  if (!is.null(ensemblIds)){
    ensemblRes <- mapENSEMBL(inputIds = ensemblIds, 
                             attributes = c('hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'), 
                             filters = c("ensembl_gene_id"), 
                             dataset = dataset, 
                             version = version, 
                             mirror = mirror, 
                             biomart = biomart)
    colnames(ensemblRes)[colnames(ensemblRes) == "ensembl_gene_id"] <- "InputIds"
    ensemblResEnsId <- inputOrder(input = ensemblIds, mergeWith = ensemblRes, inputLen = inputLen)
  } else{
    ##create empty dataframe of length(input)
    ensemblResEnsId <- data.frame(InputIds = NA, 
                                  SortOrder = as.numeric(1:length(geneSymbols)), 
                                  hgnc_symbol = NA, 
                                  external_gene_name = NA)
  }
  
  if (!is.null(geneSymbols)){
    ## better to use filters = c("hgnc_symbol") than 'external_gene_names'  
    ##because of numerous ensemblIds matches with generic names like
    ##5S_rRNA, Metazoa_SRP, Y_RNA, U1-6. 
    ensemblRes <- mapENSEMBL(inputIds = geneSymbols, 
                             attributes = c('hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'), 
                             filters = c("hgnc_symbol"), 
                             dataset = dataset, 
                             version = version, 
                             mirror = mirror, 
                             biomart = biomart)
    
    colnames(ensemblRes)[colnames(ensemblRes) == "hgnc_symbol"] <- "InputIds"
    ensemblResSymbId <- inputOrder(input = geneSymbols, mergeWith = ensemblRes, inputLen = inputLen)
  } else{
    ##create empty dataframe of length(input)
    ensemblResSymbId <- data.frame(InputIds = NA, 
                                   SortOrder = as.numeric(1:length(ensemblIds)), 
                                   ensembl_gene_id = NA, 
                                   external_gene_name = NA)
  }
  
  #reconcile/consensus between ensembl Id and symbol
  if (replaceUnmatched == T){
    coalesceEnsembl <- coalesceAlias(ensInput = ensemblResEnsId, symbInput = ensemblResSymbId, 
                                    toCoalesce_x = "external_gene_name.x", toCoalesce_y = "external_gene_name.y",
                                    replaceUnmatched = TRUE)
  } else{
    coalesceEnsembl <- coalesceAlias(ensInput = ensemblResEnsId, symbInput = ensemblResSymbId, 
                                    toCoalesce_x = "external_gene_name.x", toCoalesce_y = "external_gene_name.y")
  }
  
  
  colnames(coalesceEnsembl)[colnames(coalesceEnsembl) == "Coalesced"] <- "Coalesced.ENSEMBL"
  #write.table(coalesceEnsembl, file = "ensemblCoalesce.csv", sep = ",", row.names = F)
  
  return(coalesceEnsembl)
}

####################################################################################################################
#' @title aliasSTRINGdb
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param geneSymbols A vector of gene symbols, passed to the getBM() filters argument
#' @param speciesId Species ID. see Stringdb reference for list of avialable species
#' @param replaceUnmatched Logical. If TRUE, removes NA's and replaces with inputs
#' @import STRINGdb 
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

aliasSTRINGdb <- function(ensemblIds = NULL, geneSymbols = NULL, 
                          speciesId = 9606,  
                          replaceUnmatched = TRUE){
  #check ensemblIds, geneSymbols for NA/NULL.
  # if both non-null, make sure they are the same length
  inputLen <- checkInputLen(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
  
  queryIds<- character()
  if (!is.null(ensemblIds)){
    queryIds <- c(queryIds, ensemblIds)
    
    # stringRes <- mapSTRINGdb(inputIds = geneSymbols, speciesId = speciesId, score_threshold = score_threshold)
    # stringResSymbId <- inputOrder(input = geneSymbols, mergeWith = stringRes)
  } else{
    ##create empty dataframe of length(input)
    stringResEnsId <- data.frame(InputIds = NA, 
                                 SortOrder = as.numeric(1:length(geneSymbols)), 
                                 STRING_id = NA, 
                                 STRING.aliases = NA, 
                                 STRING.symbol = NA)
  }
  
  if (!is.null(geneSymbols)){
    queryIds <- c(queryIds, geneSymbols)
    # stringRes <- mapSTRINGdb(inputIds = geneSymbols, speciesId = speciesId, score_threshold = score_threshold)
    # stringResSymbId <- inputOrder(input = geneSymbols, mergeWith = stringRes)
  } else{
    ##create empty dataframe of length(input)
    stringResSymbId <- data.frame(InputIds = NA, 
                                  SortOrder = as.numeric(1:length(ensemblIds)), 
                                  STRING_id = NA, 
                                  STRING.aliases = NA, 
                                  STRING.symbol = NA)
  }
  
  queryIds <- unique(queryIds)
  stringRes <- mapSTRINGdb(inputIds = queryIds, 
                           speciesId = speciesId)
  stringResEnsId <- inputOrder(input = ensemblIds, mergeWith = stringRes, inputLen = inputLen)
  stringResSymbId <- inputOrder(input = geneSymbols, mergeWith = stringRes, inputLen = inputLen)
  

  #reconcile/consensus between ensembl Id and symbol
  if (replaceUnmatched == TRUE){
    coalesceString <- coalesceAlias(ensInput = stringResEnsId, symbInput = stringResSymbId, 
                                   toCoalesce_x = "STRING.symbol.x", toCoalesce_y = "STRING.symbol.y",
                                   replaceUnmatched = TRUE)
  } else{
    coalesceString <- coalesceAlias(ensInput = stringResEnsId, symbInput = stringResSymbId, 
                                   toCoalesce_x = "STRING.symbol.x", toCoalesce_y = "STRING.symbol.y")
  }
  
  colnames(coalesceString)[colnames(coalesceString) == "Coalesced"] <- "Coalesced.STRINGDB"
  #write.table(coalesceString, file = "stringCoalesce.csv", sep = ",", row.names = F)
  
  return(coalesceString)
}

####################################################################################################################
#' @title aliasDAVID
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param geneSymbols A vector of gene symbols, passed to the getBM() filters argument
#' @param email email account registered with DAVIDWebService 
#' @param replaceUnmatched Logical. If TRUE, removes NA's and replaces with inputs
#' @import RDAVIDWebService 
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

aliasDAVID <- function(ensemblIds = NULL, geneSymbols = NULL, email, 
                       replaceUnmatched = TRUE){
  ##geneSymbols = NULL for now as idType = 'OFFICIAL_GENE_SYMBOL' is not available in RDAVIDWebService for now
  
  #check ensemblIds, geneSymbols for NA/NULL.
  # if both non-null, make sure they are the same length
  inputLen <- checkInputLen(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
  
  if (!is.null(ensemblIds)){
    davidRes <- mapDAVID(inputIds = ensemblIds,
                         email = email,
                         idType = 'ENSEMBL_GENE_ID')
    davidResEnsId <- inputOrder(input = ensemblIds, mergeWith = davidRes, inputLen = inputLen)
  } else{
    ##create empty dataframe of length(input)
    davidResEnsId <- data.frame(InputIds = NA, 
                                SortOrder = as.numeric(1:length(geneSymbols)), 
                                DAVID.gene_name = NA, 
                                DAVID.symbol = NA)
  }
  
  ##idType = 'OFFICIAL_GENE_SYMBOL' is not available in RDAVIDWebService for now
  if (!is.null(geneSymbols)){
    davidRes <- mapDAVID(inputIds = geneSymbols, 
                         email = email,
                         idType = 'ENSEMBL_GENE_ID')
    davidResSymbId <- inputOrder(input = geneSymbols, mergeWith = davidRes, inputLen = inputLen)
  } else{
    ##create empty dataframe of length(input)
    davidResSymbId <- data.frame(InputIds = NA, 
                                 SortOrder = as.numeric(1:length(ensemblIds)), 
                                 DAVID.gene_name = NA, 
                                 DAVID.symbol = NA)
  }
  
  #reconcile/consensus between ensembl Id and symbol
  if (replaceUnmatched == TRUE){
    coalesceDavid <- coalesceAlias(ensInput = davidResEnsId, symbInput = davidResSymbId, 
                                   toCoalesce_x = "DAVID.symbol.x", toCoalesce_y = "DAVID.symbol.y",
                                   replaceUnmatched = TRUE)
  } else{
    coalesceDavid <- coalesceAlias(ensInput = davidResEnsId, symbInput = davidResSymbId, 
                                   toCoalesce_x = "DAVID.symbol.x", toCoalesce_y = "DAVID.symbol.y")
  }
  
  colnames(coalesceDavid)[colnames(coalesceDavid) == "Coalesced"] <- "Coalesced.DAVID"
  #write.table(coalesceDavid, file = "davidCoalesce.csv", sep = ",", row.names = F)
  
  return(coalesceDavid)
}

####################################################################################################################
#' @title aliasTable
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param geneSymbols A vector of gene symbols, passed to the getBM() filters argument
#' @param ensemblAttributes A vector of ensembl attributes
#' @param ensemblDataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly to biomaRt::useEnsembl
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param davidEmail email account registered with DAVIDWebService 
#' @param stringSpeciesId species ID. see Stringdb for list of avialable species
#' @param aliasPriorityOrder vector containig priority order of alias database. Must be UPPERCASE. Current databases: ENSEMBL, STRING, DAVID
#' @importFrom biomaRt useEnsembl getBM
#' @import STRINGdb
#' @import RDAVIDWebService 
#' @importFrom dplyr %>% group_by summarize mutate_all coalesce
#' @importFrom BiocGenerics order
#' @importFrom stringr str_match
#' @importFrom stats setNames
#' @export

aliasTable <- function(ensemblIds = NULL, geneSymbols = NULL, davidEmail,
                       ensemblAttributes = c('hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'), 
                       ensemblDataset = "mmulatta_gene_ensembl", ensemblVersion = NULL, ensemblMirror = "uswest", biomart = 'ensembl', 
                       stringSpeciesId = 9606,
                       aliasPriorityOrder = c("ENSEMBL", "STRING", "DAVID")){
  
  ret.ensembl <- aliasENSEMBL(ensemblIds = ensemblIds, 
                              geneSymbols = geneSymbols, 
                              dataset = ensemblDataset, 
                              attributes = ensemblAttributes, 
                              version = ensemblVersion, 
                              mirror = ensemblMirror, 
                              biomart = biomart, 
                              replaceUnmatched = FALSE)
  
  ##Aadd aliasSTRINGdb to aliasTable after bioconductor update of STRING_db (when mmlulatta is available)
  #https://bioc.ism.ac.jp/packages/3.1/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
  ret.string <- aliasSTRINGdb(ensemblIds = ensemblIds, 
                              geneSymbols = geneSymbols,
                              speciesId = stringSpeciesId, 
                              replaceUnmatched = FALSE)
  
  ret.david <- aliasDAVID(ensemblIds = ensemblIds, 
                          geneSymbols = geneSymbols,
                          email = davidEmail, 
                          replaceUnmatched = FALSE)
  
  
  ##merge outputs of ensembl, string, david
  merged <- merge(ret.ensembl, ret.string, by = c("SortOrder", "InputIds.EnsemblId", "InputIds.GeneSymbol"), all.x = T)
  final <- merge(merged, ret.david, by = c("SortOrder", "InputIds.EnsemblId", "InputIds.GeneSymbol"), all.x = T)
  final$SortOrder <- as.numeric(as.character(final$SortOrder))
  final <- final[order(as.numeric(final[["SortOrder"]])),]
  final <- final[ , !(names(final) %in% "SortOrder")]
  
  ##add consensus column
  final <- mutate_all(final, as.character)
  final[final == "" | final == "NA"] <- NA
  
  ##order according to preference for coalesce
  ##TODO: add check for aliasPriorityOrder
  consensus <- coalesce(final[[grep(aliasPriorityOrder[1], colnames(final))]],
                        final[[grep(aliasPriorityOrder[2], colnames(final))]],
                        final[[grep(aliasPriorityOrder[3], colnames(final))]],
                        final[["InputIds.GeneSymbol"]], 
                        final[["InputIds.EnsemblId"]])
                        
  final <- data.frame(final, Consensus = as.vector(consensus))
  
  ##TODO: fix unique alias list
  # final$All_Aliases <- apply(final, 1, function(x)unique(x))
  # final$All_Aliases <- paste(sapply(final$All_Aliases, paste, collapse=', '))
  
  return(final)
}


