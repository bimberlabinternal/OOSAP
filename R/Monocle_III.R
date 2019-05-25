

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
Ser2Monicle_MakeNProcess <- function(SeurObj = NULL, retunMon = T, PCAnDim = 20,
                                     doUMAP = T, min_dist=.3, n_neighbors = 40,
                                     dotSNE = F,
                                     doClust = T, ClusLouvRes = 0.00005, louvain_iter = 3,
                                     KeepTopNgenes = 3000){

  SeurObj.counts   <- SeurObj@assays$RNA@counts
  SeurObj.meta.data <- SeurObj@meta.data

  gene_ann <- data.frame(gene_short_name = row.names(SeurObj.counts), row.names = row.names(SeurObj.counts))

  pd <- new("AnnotatedDataFrame",data=SeurObj.meta.data)
  fd <- new("AnnotatedDataFrame",data=gene_ann)

  SeurObj_cds <- newCellDataSet(SeurObj.counts,
                                phenoData = pd, featureData =fd,
                                expressionFamily = negbinomial.size(),
                                lowerDetectionLimit=1)


  SeurObj_cds <- detectGenes(SeurObj_cds, min_expr = 0.1)


  # plot(density(SeurObj_cds$num_genes_expressed))



  SeurObj_cds <- estimateSizeFactors(SeurObj_cds)
  SeurObj_cds <- estimateDispersions(SeurObj_cds)

  disp_table = dispersionTable(SeurObj_cds)

  disp_table = disp_table %>% mutate(excess_disp =
                                       (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))

  top_subset_genes = as.character(head(disp_table, KeepTopNgenes)$gene_id)

  SeurObj_cds = setOrderingFilter(SeurObj_cds, top_subset_genes)


  pData(SeurObj_cds)$Total_mRNAs <- Matrix::colSums(exprs(SeurObj_cds))


  upper_bound <- 10^(mean(log10(pData(SeurObj_cds)$Total_mRNAs)) +
                       2*sd(log10(pData(SeurObj_cds)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(SeurObj_cds)$Total_mRNAs)) -
                       2.5*sd(log10(pData(SeurObj_cds)$Total_mRNAs)))




  SeurObj_cds <- SeurObj_cds[,pData(SeurObj_cds)$Total_mRNAs > lower_bound &
                               pData(SeurObj_cds)$Total_mRNAs < upper_bound]

  SeurObj_cds <- detectGenes(SeurObj_cds, min_expr = 0.1)



  SeurObj_cds <- preprocessCDS(SeurObj_cds,
                               method = 'PCA',
                               norm_method = 'log',
                               num_dim = PCAnDim,
                               verbose = T)

  if(doUMAP) SeurObj_cds <- reduceDimension(SeurObj_cds, max_components = 2,
                                            reduction_method = 'UMAP',
                                            metric="correlation",
                                            min_dist = min_dist,
                                            n_neighbors = n_neighbors,
                                            verbose = T)

  # plot_pc_variance_explained(SeurObj_cds, return_all = F) # norm_method='log'

  if(dotSNE) SeurObj_cds <- reduceDimension(SeurObj_cds, max_components = 2, num_dim = 15,
                                            reduction_method = 'tSNE', verbose = T)

  if(doClust) SeurObj_cds <- clusterCells(SeurObj_cds,
                                          method = 'louvain',
                                          res = ClusLouvRes,
                                          louvain_iter = louvain_iter,
                                          verbose = T)





  return(SeurObj_cds)


}

