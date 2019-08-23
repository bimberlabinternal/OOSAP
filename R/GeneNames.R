
#' @title QueryEnsemblSymbolAndHumanHomologs
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param ensemblFilters A vector of ensembl IDs, passed to the getBM() filters argument.
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @importFrom biomaRt useEnsembl getBM
#' @export
QueryEnsemblSymbolAndHumanHomologs <- function(ensemblIds, biomart = "ensembl", dataset = "mmulatta_gene_ensembl", ensemblFilters = c('ensembl_gene_id')) {
    ensembl = useEnsembl(biomart=biomart, dataset=dataset)
    homologAttrs <- c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name')
    homologs <- getBM(attributes=homologAttrs, filters = ensemblFilters, values = ensemblIds, mart = ensembl)

    ret <- data.frame(ensembl_gene_id = ensemblIds, SortOrder = 1:length(ensemblIds))
    ret <- merge(ret, homologs, by = c('ensembl_gene_id'), all.x = T)
    ret <- ret[order(ret$SortOrder),]

    ret$Label <- as.character(ret$ensembl_gene_id)
    ret$Label[!is.na(ret$hsapiens_homolog_associated_gene_name)] <- ret$hsapiens_homolog_associated_gene_name[!is.na(ret$hsapiens_homolog_associated_gene_name)]
    ret$Label[!is.na(ret$external_gene_name)] <- ret$external_gene_name[!is.na(ret$external_gene_name)]

    return(ret)
}

#Based on: https://www.genenames.org/data/genegroup/#!/group/471
RenameGenesUsingCD <- function(geneSymbols) {
    df <- data.frame(GeneSymbol = geneSymbols, SortOrder = 1:length(geneSymbols), stringsAsFactors = F)
    ret <- merge(df, cdGenes, all.x = T, by = c('GeneSymbol'))
    ret <- ret[order(ret$SortOrder),]

    ret$Label <- as.character(ret$GeneSymbol)
    ret$Label[!is.na(ret$PreviousSymbols)] <- paste0(ret$Label[!is.na(ret$PreviousSymbols)], ' (', as.character(ret$PreviousSymbols[!is.na(ret$PreviousSymbols)]), ')')

    return(ret$Label)
}


#' @title AliasGeneNames
#' @description This function accepts a vector of gene names. It will first query Ensembl to attempt to translate ID to gene symbol. Any names not translated will be compared for human homologs. Finally, these symbols will be compared to CD genes, and the CD name appended (if found).
#' @param ensemblIds A vector of ensembl IDs, passed to the getBM() filters argument
#' @param ensemblFilters A vector of ensembl IDs, passed to the getBM() filters argument.
#' @param biomart Passed directly to biomaRt::useEnsembl
#' @param dataset Passed directly to biomaRt::useEnsembl
#' @importFrom biomaRt useEnsembl getBM
#' @export
AliasGeneNames <- function(ensemblIds, biomart = "ensembl", dataset = "mmulatta_gene_ensembl", ensemblFilters = c('ensembl_gene_id')) {
    ret <- QueryEnsemblSymbolAndHumanHomologs(ensemblIds, biomart = biomart, dataset = dataset, ensemblFilters = ensemblFilters)

    return(RenameGenesUsingCD(ret$Label))
}