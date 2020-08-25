
#' @title QueryEnsemblSymbolAndHumanHomologs
#' @param so A Seurat Object
#' @param group.by a parameter in the Seurat Objects's meta to sample againts factor-levels
#' @param downSample If True random sampling is performed to return equal number of cells across group.by. if sampleNmax = NULL this is set by min cells found across all group.by factor levels
#' @param sampleNmax maximum number of cells to sample from, if greater than any given group.by factor level since relace = F the maximum number avaiable for that level is taken. leave as NULL if equal number of samples across all levels is desired determined by the lowest number of cells across group.by factor-levels. 
#' @param slot Passed directly to Seurat::AverageExpression slot parameter. 
#' @param genes A vector of genes desired, leave as NULL for all genes. 
cellsToSamples <- function(so=NULL, group.by=NULL, downSample=F, sampleNmax = NULL, slot="data", genes = NULL){
  # so = ComboSeuratObj_QC; group.by = "ExpID"
  if(is.null(so)) stop("No Seurat Object provided")
  if(is.null(group.by)) stop("group.by is null")
  if(!group.by %in% colnames(so@meta.data)) stop("group.by not in Seurat Object")
  if(is.null(genes)) genes = rownames(so)
  
  so$CellId_2_samp <- colnames(so)
  
  if(is.null(sampleNmax)) sampleNmax = min(table(so[[group.by]]))
  print(paste0("smallest sample size is: ", min(table(so[[group.by]]))))
  
  GroupNames <- table(so[[group.by]]) %>% unlist() %>% names()
  
  selectedCellIdsPerGroup <- lapply(GroupNames, function(GN){
    if(downSample){
      rowIds = rownames(so@meta.data)[so@meta.data[,group.by]==GN]
      sampleN = min(c(sampleNmax, length(rowIds)))
      rowIds[sample(1:length(rowIds), sampleN, replace = F)]
    } else {
      rownames(so@meta.data)[so@meta.data[,group.by]==GN]
    }
    
  })
  so <- subset(so, CellId_2_samp %in% unlist(selectedCellIdsPerGroup))
  
  AvgExprDGE <- OOSAP::AvgCellExprs(so, varName = group.by, slot=slot, genes = genes)
  
  
  return(AvgExprDGE)
  
  
}