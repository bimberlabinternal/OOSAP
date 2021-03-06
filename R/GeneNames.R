
#' @title QueryEnsemblSymbolAndHumanHomologs
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param ensemblFilters A vector of ensembl IDs, passed to the getBM() filters argument.
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @param version Passed directly to biomaRt::useEnsembl
#' @param extraAttrs A vector of attributes that biomaRt::getBM() should get beyond default set. see listAttributes.
#' @param ensemblMirror Passed directly to biomaRt::useEnsembl.  Defaults to uswest since the default (NULL) was given SSL errors in docker
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr %>% group_by summarise
#' @export
QueryEnsemblSymbolAndHumanHomologs <- function(ensemblIds, biomart = "ensembl", dataset = "mmulatta_gene_ensembl", ensemblFilters = c('ensembl_gene_id'), version = NULL, ensemblMirror = NULL, extraAttrs = NULL) {
    ensembl = biomaRt::useEnsembl(biomart=biomart, dataset=dataset, version = version, mirror = ensemblMirror)
    homologAttrs <- c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name')
    if (!is.null(extraAttrs)) {
        homologAttrs <- unique(c(homologAttrs, extraAttrs))
    }

    homologs <- NULL
    tryCatch(expr = {
        ensemblIdsToQuery <- ensemblIds[grepl(ensemblIds, pattern = '^ENS')]
        if (length(ensemblIdsToQuery) > 0) {
            homologs <- biomaRt::getBM(attributes=homologAttrs, filters = ensemblFilters, values = ensemblIds, mart = ensembl)
            homologs <- homologs %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::summarise(
                ensembl_transcript_id = paste(sort(unique(ensembl_transcript_id)), collapse=','),
                external_gene_name = paste(sort(unique(external_gene_name)), collapse=','),
                hsapiens_homolog_ensembl_gene = paste(sort(unique(hsapiens_homolog_ensembl_gene)), collapse=','),
                hsapiens_homolog_associated_gene_name = paste(sort(unique(hsapiens_homolog_associated_gene_name)), collapse=',')
            )

            if (nrow(homologs) > 0 && sum(homologs == '') > 0) {
                homologs[homologs == ''] <- NA
            }
        } else {
            print('None of the input IDs started with ENS, skipping emsembl query')
        }
    }, error = function(e){
        print(paste0('Error: ', e))
        stop(paste0('Error querying ensembl, using genes: ', paste0(ensemblIds, collapse = ';')))
    })

    ret <- data.frame(ensembl_gene_id = ensemblIds, SortOrder = 1:length(ensemblIds))
    ret$Label <- as.character(ret$ensembl_gene_id)
    if (!all(is.null(homologs))) {
        ret <- merge(ret, homologs, by = c('ensembl_gene_id'), all.x = T)
        ret$Label[!is.na(ret$hsapiens_homolog_associated_gene_name)] <- paste0(ret$hsapiens_homolog_associated_gene_name[!is.na(ret$hsapiens_homolog_associated_gene_name)], '(hs)')
        ret$Label[!is.na(ret$external_gene_name)] <- ret$external_gene_name[!is.na(ret$external_gene_name)]
    } else {
        ret$ensembl_transcript_id <- NA
        ret$external_gene_name <- NA
        ret$hsapiens_homolog_ensembl_gene <- NA
        ret$hsapiens_homolog_associated_gene_name <- NA
    }

		ret <- ret[order(ret$SortOrder),]
    ret <- ret[names(ret) != 'SortOrder']

    return(ret)
}

#Based on: https://www.genenames.org/data/genegroup/#!/group/471
RenameGenesUsingCD <- function(geneSymbols) {
    df <- data.frame(GeneSymbol = geneSymbols, SortOrder = 1:length(geneSymbols), stringsAsFactors = F)
    ret <- merge(df, OOSAP::cdGenes, all.x = T, by = c('GeneSymbol'))
    ret <- ret[order(ret$SortOrder),]

    ret$Label <- as.character(ret$GeneSymbol)
    sel <- !is.na(ret$PreviousSymbols) & !is.null(ret$PreviousSymbols) & ret$PreviousSymbols != ''
    ret$Label[sel] <- paste0(ret$Label[sel], ' (', as.character(ret$PreviousSymbols[sel]), ')')

    return(ret$Label)
}


#' @title AliasGeneNames
#' @description This function accepts a vector of gene names. It will first query Ensembl to attempt to translate ID to gene symbol. Any names not translated will be compared for human homologs. Finally, these symbols will be compared to CD genes, and the CD name appended (if found).
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param ensemblFilters A vector of ensembl IDs, passed to the getBM() filters argument.
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @param ensemblVersion Passed directly as version to biomaRt::useEnsembl
#' @importFrom biomaRt useEnsembl getBM
#' @export
AliasGeneNames <- function(ensemblIds, biomart = "ensembl", dataset = "mmulatta_gene_ensembl", ensemblFilters = c('ensembl_gene_id'), ensemblVersion = NULL) {
    ret <- QueryEnsemblSymbolAndHumanHomologs(ensemblIds, biomart = biomart, dataset = dataset, ensemblFilters = ensemblFilters, version = ensemblVersion)

    return(RenameGenesUsingCD(ret$Label))
}