
#' @title RecreateSeuratObjFromSDAmatrix
#' @param cellnames cell names
#' @param recomposedMat recomposed matrix
#' @param metaDF metadata datframe
#' @param projectName Title of project
#' @param dimsToUse The number of dims to use
#' @param spikeGenes list of genes to be added beyond the default VariableGenes()
#' @param nfeatures The number of VariableFeatures to identify
#' @param tSNE_perplexity tSNE perplexity. Passed directly to Seurat::RunTSNE()
#' @param UMAP_MinDist UMAP min. distance. Passed directly to Seurat::RunUMAP()
#' @param UMAP_NumNeib UMAP number of neighboring points. Passed directly to Seurat::RunUMAP()
#' @param UMAP.NumEpoc UMAP min. distance. Passed directly to Seurat::RunUMAP()
#' @import Seurat
#' @export
RecreateSeuratObjFromSDAmatrix <- function(cellnames = NULL, recomposedMat = NULL, metaDF, projectName = "SDAproject", dimsToUse = 15, spikeGenes = NULL, nfeatures = 500, 
                                           tSNE_perplexity = 0.3, UMAP_MinDist = 0.3, UMAP_NumNeib = 35L, UMAP.NumEpoc = 500) {
  
  #confirm that gene are rows and columns are cells
  if(any(grepl("CD1", rownames(recomposedMat)))) {
    print("genes in rownames, good")
  } else { recomposedMat <- t(recomposedMat) }
  
  recomposedMat <- recomposedMat[, cellnames]
  
  recomposedMat <- as(recomposedMat, "dgCMatrix")
  # recomposedMat <- recomposedMat + 10 #abs(round(min(recomposedMat)))+1
  metaDF <-  as.data.frame(metaDF[which(metaDF$Barcode %in% cellnames), ])
  
  rownames(metaDF) <- metaDF$Barcode
  #genes should be in rows and cells in columns. 
  
  # with imputed full matrix, create seurat object 
  seuratObj <- CreateSeuratObject(counts = recomposedMat, min.features = 10, project = projectName, meta.data = metaDF)
  
  #TODO: if the imputed is pre-normalized is this redundant?
  seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  tempData <- as.data.frame(seuratObj@assays$RNA@data)
  tempData[is.na(tempData)] <- 0
  tempData <- as(as.matrix(tempData), "dgCMatrix")
  seuratObj@assays$RNA@data <- tempData
  
  DefaultAssay(seuratObj) <- "RNA"
  seuratObj <- FindVariableFeatures(object = seuratObj, selection.method="vst", nfeatures = nfeatures)
  if(!is.null(spikeGenes)) {
    seuratObj@assays$RNA@var.features <- unique(c(spikeGenes, seuratObj@assays$RNA@var.features))
  }
  
  # (is.na(SeurComboObj_SDA@assays$RNA@counts))
  #NAs to 0, since some cells become NA by normalization step and mess things downstream for DE
  seuratObj <- ScaleData(object = seuratObj, min.cells.to.block = 2000, features = rownames(seuratObj))
  seuratObj <- RunPCA(object = seuratObj, npcs = 40)
  ElbowPlot(seuratObj)
  
  seuratObj <- FindNeighbors(object = seuratObj, dims = 1:dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = 0.1)
  seuratObj <- FindClusters(object = seuratObj, resolution = 0.6)
  seuratObj <- FindClusters(object = seuratObj, resolution = 1.2)
  seuratObj <- RunTSNE(object = seuratObj, dims = 1:dimsToUse, perplexity = tSNE_perplexity, max_iter = 5000, num_threads = 10)
  seuratObj <- RunUMAP(object = seuratObj, dims = 1:dimsToUse, n.neighbors = UMAP_NumNeib, min.dist = UMAP_MinDist, 
                       metric = "correlation", seed.use = 11358, umap.method = "umap-learn", n.epochs = UMAP.NumEpoc)
  
  DimPlot(seuratObj)
  
  return(seuratObj)
}




#' @title RecomposeSDAmatrix
#' @param SDAseuratObj Final seurat object from SDA processing
#' @param compKeepSet Three (3) options, "all", "qc", or "batch" default to "batch" which is comps removed by ShinySDA.
#' @param metaDF metadata datframe
#' @param spikeGenes list of genes to be added beyond the default VariableGenes()
#' @param nfeatures The number of VariableFeatures to identify
#' @param tSNE_perplexity tSNE perplexity. Passed directly to Seurat::RunTSNE()
#' @param UMAP_MinDist UMAP min. distance. Passed directly to Seurat::RunUMAP()
#' @param UMAP_NumNeib UMAP number of neighboring points. Passed directly to Seurat::RunUMAP()
#' @param UMAP.NumEpoc UMAP min. distance. Passed directly to Seurat::RunUMAP()
#' @param saveFile If provided, the seurat object will be saved as RDS to this location
#' @import Seurat
#' @export
RecomposeSDAmatrix <- function(SDAseuratObj = NULL, compKeepSet = "batch", 
                               metaDF = NULL, spikeGenes = NULL, nfeatures = 2000, 
                               tSNE_perplexity = 100, UMAP_MinDist = 0.5, UMAP_NumNeib = 60L, saveFile = NULL) {
  
  
  if(! compKeepSet %in% c("all", "qc", "batch")) {
    warning("compKeepSet was not correctly parameterized: choose all, qc, or batch default to batch")
    compKeepSet = "batch"
  }
  
  SDAseuratObj$Barcode <- rownames(SDAseuratObj@meta.data)
  
  #combine SDAseuratObj metadata and raw seurat metadata 
  if(!is.null(metaDF)){
    if(!"Barcode" %in% colnames(metaDF)) 
      metaDF$Barcode <- rownames(metaDF)
    if("Barcode" %in% colnames(metaDF)){
      metaDF <- merge(SDAseuratObj@meta.data, metaDF, by="Barcode", all = T, suffixes = c("SDA", "Seurat"))
      if(nrow(metaDF) == nrow(SDAseuratObj@meta.data)){
        rownames(metaDF) <- metaDF$Barcode
        metaDF <- metaDF[rownames(SDAseuratObj@meta.data), ]
      } else {
        metaDF <- SDAseuratObj@meta.data
      }
    }
  } else {
    metaDF <- SDAseuratObj@meta.data
  }
  
  #keep only passing componenets or all componenets
  if(compKeepSet == "all") {
    compsToKeep <- 1:ncol(SDAseuratObj@reductions$SDA@cell.embeddings)
  } else {
    if(compKeepSet == "qc"){
      compsToKeep <- setdiff(1:ncol(SDAseuratObj@reductions$SDA@cell.embeddings), as.numeric(SDAseuratObj@misc$SDA_processing_results$QC_components))
    } else {
      if(compKeepSet == "batch"){
        compsToKeep <- setdiff(1:ncol(SDAseuratObj@reductions$SDA@cell.embeddings), as.numeric(SDAseuratObj@misc$SDA_processing_results$Remove_comps))
      }
    }
  }
    
    
  
  #do dot product of cell embeddings and gene loadings to get approximate/imputed raw seurat matrix
  print("Performing dot product...")
  recomposedMat <-  t(SDAseuratObj@reductions$SDA@cell.embeddings[,compsToKeep] %*% t(SDAseuratObj@reductions$SDA@feature.loadings)[compsToKeep,])
  print("Finished dot product.")
  
  if(saveFile=="" | is.null(saveFile)){ saveFile = "TempName" }
  
  #create seurat object with imputed matrix and run standard seuart processes
  print("Starting RecreateSeuratObjFromSDAmatrix...")
  imputedSeuratObj <- RecreateSeuratObjFromSDAmatrix(cellnames = colnames(recomposedMat), recomposedMat = recomposedMat,  metaDF = metaDF, 
                                                     projectName = gsub(pattern = ".rds", replacement = "", (saveFile)), spikeGenes = spikeGenes, 
                                                     tSNE_perplexity = tSNE_perplexity, UMAP_MinDist = UMAP_MinDist, UMAP_NumNeib = UMAP_NumNeib, nfeatures = nfeatures)
  return(imputedSeuratObj)
}


