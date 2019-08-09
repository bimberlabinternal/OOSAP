


#' @title generate SingleR object
#'
#' @description Compute SingleR classification on a Seurat object
#' @param SeurObjPath, path to Seurat object
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @return SingleR object or Seurat object
#' @keywords SingleR Classification
#' @export
#' @import Seurat
#' @import SingleR
#' @importFrom SingleR hpca
generateSingleR <- function(SeurObj = NULL, Annotation = Annotation, ProjName = ProjName, MinGenes = MinGenes, Species = Species,
                            ClusteringName = ClusteringName, NumCores = NumCores, NormGeneLength = F, VarGeneMethod = "de",
                            FineTune = T, DoSig = T, MainTypes=T, RedFS = T){
  
  Counts.Mat <- as.matrix(SeurObj@assays$RNA@counts)

  #Attempt to force-loading of the annotations.  for some reason this must be defined in the global scope?  could have something to do with lazy data loading?
  hpca <<- SingleR::hpca

  # if(length(grep("MTOR", rownames(SeurObj)))>0 | length(grep("CD1", rownames(SeurObj)))>0 | length(grep("ENS", rownames(SeurObj)))>0){
  #   print("genes in rownames, good!")
  # } else {
  #   Counts.Mat <- t(Counts.Mat)
  # }
  singler = CreateSinglerObject(counts = Counts.Mat,
                                annot = Annotation,
                                project.name = ProjName,
                                min.genes = MinGenes,
                                species = Species, 
                                citation = "",
                                ref.list = list(), 
                                normalize.gene.length = NormGeneLength, 
                                variable.genes = "de", # or sd
                                fine.tune = FineTune, 
                                do.signatures = DoSig, 
                                clusters = ClusteringName, 
                                do.main.types = MainTypes, 
                                reduce.file.size = RedFS, 
                                numCores = NumCores)

  #cleanup global variable
  rm(hpca)

  return(singler)
  
  
}


#' @title SingleR my Seurat Object
#'
#' @description Compute SingleR classification on a Seurat object from path or direcly
#' @param SeurObjPath, path to Seurat object
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @return SingleR object or Seurat object
#' @keywords SingleR Classification
#' @export
#' @import Seurat
#' @import SingleR
SingleRmySerObj <- function(SeurObjPath = NULL, SeurObj = NULL, PlotFigs = T, 
                            Annotation = NULL, ProjName = NULL , 
                            MinGenes = 500, Species = "Human", NumCores = 4, 
                            ClusteringName = NULL, SavePath = NULL, ReturnSeurObj = F){
  

  #Annotation can be any Idents from Metadata
  #ClusteringName input cluster id for each of the cells with at least min.genes, if NULL uses SingleR clusterings.
  
  if(is.null(SeurObj) & is.null(SeurObjPath)) stop("Seurat object not def, SeurObjPath or SeurObj needed.")
  
  if((!is.null(SavePath)) & (!dir.exists(dirname(SavePath))) ) stop("Save path does not exist")
  
    
   
    
    
  
  if(!is.null(SeurObjPath)) SeurObj <- readRDS(SeurObjPath)
  
  dim(SeurObj)
  
  if(PlotFigs) DimPlot(SeurObj) + theme_bw()
  
  if(is.null(ProjName)) ProjName = SeurObj@project.name
  
  
  # SeuratMeta <- as.data.frame(SeurObj@meta.data)
  # if(!is.null(Annotation))
  
  
  singler <- generateSingleR(SeurObj = SeurObj, Annotation = Annotation, ProjName = ProjName, MinGenes = MinGenes, Species = Species, 
                  ClusteringName = ClusteringName, NumCores = NumCores, NormGeneLength = F, VarGeneMethod = "de",
                  FineTune = T, DoSig = T, MainTypes=T, RedFS = T)
  
  
  
  
  
  if(!is.null(SavePath)) saveRDS(singler, SavePath)
  
  
  if(ReturnSeurObj){
    SeurObj$SingleR_Labels1 <- "unk"
    SeurObj@meta.data[names(singler$singler[[1]]$SingleR.single$labels[,1]),]$SingleR_Labels1 <- singler$singler[[1]]$SingleR.single$labels[,1]

    SeurObj$SingleR_Labels2 <- "unk"
    SeurObj@meta.data[names(singler$singler[[2]]$SingleR.single$labels[,1]),]$SingleR_Labels2 <- singler$singler[[2]]$SingleR.single$labels[,1]

    SeurObj$SingleR_LabelsOther <- "unk"
    SeurObj@meta.data[names(singler$other),]$SingleR_LabelsOther <- singler$other

    return(SeurObj)
  } else {
    return(singler)
  }
  
  

  
}





#' @title Save SingleR to SeurObj
#'
#' @description Add a pre-Computed SingleR classification into a Seurat object from path or direcly
#' @param SeurObjPath, path to Seurat object
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @param SingleRPath path to SingleR object
#' @return modified Seurat object
#' @keywords SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
SaveSingleRtoSeurObj <- function(SeurObjPath = NULL, SeurObj = NULL, SingleRPath = NULL){
  
  if(is.null(SeurObj) & is.null(SeurObjPath)) stop("Seurat object not def, SeurObjPath or SeurObj needed.")
  if(!is.null(SeurObjPath)) SeurObj <- readRDS(SeurObjPath)
  
  if(!is.null(SingleRPath)) singler <- readRDS(SingleRPath) else stop("SingleR path NULL")
  
  SeurObj$SingleR_Labels1 <- "unk"
  SeurObj@meta.data[names(singler$singler[[1]]$SingleR.single$labels[,1]),]$SingleR_Labels1 <- singler$singler[[1]]$SingleR.single$labels[,1]
  
  SeurObj$SingleR_Labels2 <- "unk"
  SeurObj@meta.data[names(singler$singler[[2]]$SingleR.single$labels[,1]),]$SingleR_Labels2 <- singler$singler[[2]]$SingleR.single$labels[,1]
  
  SeurObj$SingleR_LabelsOther <- "unk"
  SeurObj@meta.data[names(singler$other),]$SingleR_LabelsOther <- singler$other
  
  return(SeurObj)
  

}


#' @title DimPlot SingleR Class Lables
#' @description Dimplot Seurobject with SingleR class lables
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @keywords Dimplot SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
DimPlot_SingleRClassLabs <- function(SeurObj){
  cowplot::plot_grid(
  DimPlot(SeurObj, group.by = "SingleR_Labels1") + theme_bw() + ggtitle("SingleR_Labels1") +theme(legend.position="bottom"),
  DimPlot(SeurObj, group.by = "SingleR_Labels2") + theme_bw() + ggtitle("SingleR_Labels2") +theme(legend.position="bottom"),
  DimPlot(SeurObj, group.by = "SingleR_LabelsOther") + theme_bw() + ggtitle("SingleR_LabelsOther") +theme(legend.position="bottom"), 
  ncol = 1)
  
}


#' @title Tabulate SingleR Class Lables
#' @description Tabulate Seurobject with SingleR class lables
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @keywords Tabulate SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
Tabulate_SingleRClassLabs <- function(SeurObj){
  cowplot::plot_grid(
    ggplot(melt(table(SeurObj$SingleR_Labels1)), aes(x=Var1, y = value, fill=Var1))  + 
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR predicted classification  labels #1:: TestisII \n Total Contribution") + 
      ylab("Number of cells"),
    
    
    ggplot(melt(table(SeurObj$SingleR_Labels2)), aes(x=Var1, y = value, fill=Var1))  + 
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR predicted classification labels #2:: TestisII \n Total Contribution") + 
      ylab("Number of cells"),
    
    
    ggplot(melt(table(SeurObj$SingleR_LabelsOther)), aes(x=Var1, y = value, fill=Var1))  + 
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR predicted classification labels other:: TestisII \n Total Contribution") + 
      ylab("Number of cells"),
    ncol = 1)
  
}




