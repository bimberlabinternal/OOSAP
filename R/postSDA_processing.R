
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
RecreateSeuratObjFromSDAmatrix <- function (cellnames = NULL, recomposedMat = NULL, metaDF, projectName = "SDAproject", 
                                             minDimsToUse = 15, npcs = 50, 
                                             spikeGenes = NULL, nVariableFeatures = 2000, variableFeatureSelectionMethod = "vst",
                                             tSNE_perplexity = 30, tSNE_maxIter = 2000,
                                             UMAP_NumNeib = 40L, UMAP_MinDist = 0.2, UMAP_Seed = 1234, UMAP.NumEpoc = 500,
                                             mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf)) 
{
  if (any(grepl("CD1", rownames(recomposedMat)))) {
    print("Confirmed that genes are in rows.")
  }
  else {
    recomposedMat <- t(recomposedMat)
  }
  recomposedMat <- recomposedMat[, cellnames]
  recomposedMat <- as(recomposedMat, "dgCMatrix")
  metaDF <- as.data.frame(metaDF[which(metaDF$Barcode %in% 
                                         cellnames), ])
  rownames(metaDF) <- metaDF$Barcode
  seuratObj <- CreateSeuratObject(counts = recomposedMat, 
                                  min.features = 10, project = projectName, meta.data = metaDF)
  
  print("Normalizing data:")
  seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", 
                             #scale.factor = 10000,
                             verbose = T)
  
  #drop NaN values in matrix (otherwise DE test errors)    
  data <- as.data.frame(seuratObj@assays$RNA@data)
  data[is.na(data)] <- 0
  data <- as(as.matrix(data), "dgCMatrix")
  seuratObj@assays$RNA@data <- data
  
  DefaultAssay(seuratObj) <- "RNA"
  
  print("Find variable features:")
  seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = variableFeatureSelectionMethod, nfeatures = nVariableFeatures,
                                    mean.cutoff=mean.cutoff, dispersion.cutoff=dispersion.cutoff, verbose = F)
  
  if (!is.null(spikeGenes)) {
    VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), spikeGenes))
  }
  
  seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), 
                         #min.cells.to.block = 2000, 
                         #vars.to.regress = c("nCount_RNA", "p.mito"), #slow computation: https://github.com/satijalab/seurat/issues/618
                         verbose = F)
  
  seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), verbose = F, npcs = npcs)
  
  seuratObj <- ProjectDim(object = seuratObj)
  
  ElbowPlot(seuratObj)
  
  print("Running FindClustersAndDimRedux:")
  seuratObj <- OOSAP::FindClustersAndDimRedux(seuratObj = seuratObj, forceReCalc = FALSE, minDimsToUse = 15, 
                                       maxTsneIter = tSNE_maxIter,
                                       umap.method = "uwot", 
                                       UMAP_NumNeib=UMAP_NumNeib, UMAP_MinDist=UMAP_MinDist, UMAP_Seed=UMAP_Seed, UMAP.NumEpoc=UMAP.NumEpoc)
  
  DimPlot(seuratObj)
  
  return(seuratObj)
}



#' @title RecomposeSDAmatrix
#' @param SDAseuratObj Final seurat object from SDA processing
#' @param WhichComp Three (3) options, "all", "qc", or "batch" default to "batch" which is comps removed by ShinySDA.
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
RecomposeSDAmatrix <- function(SDAseuratObj = NULL, WhichComp = "batch", metaDF = NULL, projectName = 'imputedSeuratObj',
                                spikeGenes = NULL, nVariableFeatures = 2000,  variableFeatureSelectionMethod = "vst",
                                tSNE_perplexity = 30, tSNE_maxIter = 2000,
                                UMAP_NumNeib = 40L, UMAP_MinDist = 0.2, UMAP_Seed = 1234, UMAP.NumEpoc = 500
) 
{
  if(! WhichComp %in% c("all", "qc", "batch")) {
    warning("WhichComp was not correctly parameterized: choose all, qc, or batch default to batch")
    WhichComp = "batch"
  }
  
  SDAseuratObj$Barcode <- rownames(SDAseuratObj@meta.data)
  if (!is.null(metaDF)) {
    if (!"Barcode" %in% colnames(metaDF)) 
      metaDF$Barcode <- rownames(metaDF)
    if ("Barcode" %in% colnames(metaDF)) {
      metaDF <- merge(SDAseuratObj@meta.data, metaDF, 
                      by = "Barcode", all = T, suffixes = c(".SDA", 
                                                            ".Seurat"))
      if (nrow(metaDF) == nrow(SDAseuratObj@meta.data)) {
        rownames(metaDF) <- metaDF$Barcode
        metaDF <- metaDF[rownames(SDAseuratObj@meta.data),]
      }
      else {
        metaDF <- SDAseuratObj@meta.data
      }
    }
  } else {
    metaDF <- SDAseuratObj@meta.data
  }
  
  if(WhichComp == "all") {
    cpsToKeep <- SDAseuratObj@misc$SDA_processing_results$QC_components
  } else {
    if(WhichComp == "qc"){
      compsToKeep <- SDAseuratObj@misc$SDA_processing_results$QC_components
    } else {
      if(WhichComp == "batch"){
        compsToKeep <- as.numeric(SDAseuratObj@misc$SDA_processing_results$Remove_comps)
      }
    }
  }
  
  print("Performing dot product...")
  recomposedMat <- t(SDAseuratObj@reductions$SDA@cell.embeddings[,cpsToKeep] %*% t(SDAseuratObj@reductions$SDA@feature.loadings)[cpsToKeep,])
  print("Finished dot product.")
  
  print("Starting RecreateSeuratObjFromSDAmatrix...")
  imputedSeuratObj <- RecreateSeuratObjFromSDAmatrix(cellnames = colnames(recomposedMat), 
                                                      recomposedMat=recomposedMat, metaDF=metaDF, projectName=projectName,
                                                      spikeGenes = NULL, nVariableFeatures=nVariableFeatures, 
                                                      tSNE_perplexity=tSNE_perplexity, tSNE_maxIter=tSNE_maxIter,
                                                      UMAP_NumNeib=UMAP_NumNeib, UMAP_MinDist=UMAP_MinDist, UMAP_Seed=UMAP_Seed, UMAP.NumEpoc=UMAP.NumEpoc)
  
  return(imputedSeuratObj)
}


