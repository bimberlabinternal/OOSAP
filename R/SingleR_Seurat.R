

#' @title SingleR my Seurat Object
#'
#' @description Compute SingleR classification on a Seurat object from path or direcly
#' @param SeurObjPath, path to Seurat object
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @return SingleR object or Seurat object
#' @keywords SingleR Classification
#' @export
#' @import Seurat SingleR
SingleRmySerObj <- function(SeurObjPath = NULL, SeurObj = NULL, PlotFigs = T, 
                            Annotation = NULL, ProjName = NULL , 
                            MinGenes = 500, Species = "Human", NumCores = 4, 
                            ClusteringName = NULL, SavePath = NULL, ReturnSeurObj = F){
  
  # SeurObjPath = "./data/296-6-GEX.seurat.rds"; SeurObj = NULL; PlotFigs = T
  # Annotation = NULL; ProjName = NULL; MinGenes = 500; Species = "Human"; NumCores = 4
  
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
  
  Counts.Mat <- as.matrix(SeurObj@assays$RNA@counts)
  
  
  if(length(grep("MTOR", rownames(SeurObj)))>0 | length(grep("CD1", rownames(SeurObj)))>0 | length(grep("ENS", rownames(SeurObj)))>0){
    print("genes in rownames, good!")
  } else {
    Counts.Mat <- t(Counts.Mat)
  }
  
  singler = CreateSinglerObject(counts = Counts.Mat,
                                annot = Annotation,
                                project.name = ProjName,
                                min.genes = MinGenes,
                                species = Species, 
                                citation = "",
                                ref.list = list(), 
                                normalize.gene.length = F, 
                                variable.genes = "de", # or sd
                                fine.tune = T, 
                                do.signatures = T, 
                                clusters = ClusteringName, 
                                do.main.types = T, 
                                reduce.file.size = T, 
                                numCores = NumCores)
  
  if(!is.null(SavePath)) saveRDS(singler, SavePath)
  
  
  if(ReturnSeurObj){
    
    SeurObj$SingleR_Laels1 <- "unk"
    SeurObj@meta.data[names(singler$singler[[1]]$SingleR.single$labels[,1]),]$SingleR_Laels1 <- singler$singler[[1]]$SingleR.single$labels[,1]
    
    table(SeurObj$SingleR_Laels1)
    
    
    
    
    SeurObj$SingleR_Laels2 <- "unk"
    SeurObj@meta.data[names(singler$singler[[2]]$SingleR.single$labels[,1]),]$SingleR_Laels2 <- singler$singler[[2]]$SingleR.single$labels[,1]
    
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
#' @import Seurat SingleR
SaveSingleRtoSeurObj <- function(SeurObjPath = NULL, SeurObj = NULL, SingleRPath = NULL){
  
  if(is.null(SeurObj) & is.null(SeurObjPath)) stop("Seurat object not def, SeurObjPath or SeurObj needed.")
  if(!is.null(SeurObjPath)) SeurObj <- readRDS(SeurObjPath)
  
  if(!is.null(SingleRPath)) singler <- readRDS(SingleRPath) else stop("SingleR path NULL")
  
  SeurObj$SingleR_Laels1 <- "unk"
  SeurObj@meta.data[names(singler$singler[[1]]$SingleR.single$labels[,1]),]$SingleR_Laels1 <- singler$singler[[1]]$SingleR.single$labels[,1]
  
  table(SeurObj$SingleR_Laels1)
  
  
  
  
  SeurObj$SingleR_Laels2 <- "unk"
  SeurObj@meta.data[names(singler$singler[[2]]$SingleR.single$labels[,1]),]$SingleR_Laels2 <- singler$singler[[2]]$SingleR.single$labels[,1]
  
 

  return(SeurObj)
  

}


#' @title DimPlot SingleR Class Lables
#' @description Dimplot Seurobject with SingleR class lables
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @keywords Dimplot SingleR Classification 
#' @export
#' @import Seurat SingleR cowplot
DimPlot_SingleRClassLabs <- function(SeurObj){
  cowplot::plot_grid(
  DimPlot(SeurObj, group.by = "SingleR_Laels1") + theme_bw() + ggtitle("SingleR_Laels1") +theme(legend.position="bottom"),
  DimPlot(SeurObj, group.by = "SingleR_Laels2") + theme_bw() + ggtitle("SingleR_Laels2") +theme(legend.position="bottom"), 
  ncol = 1)
  
}


#' @title Tabulate SingleR Class Lables
#' @description Tabulate Seurobject with SingleR class lables
#' @param SeurObj a Seurat object, but if path given, path is prioritized. 
#' @keywords Tabulate SingleR Classification 
#' @export
#' @import Seurat SingleR cowplot
Tabulate_SingleRClassLabs <- function(SeurObj){
  cowplot::plot_grid(
    ggplot(melt(table(SeurObj$SingleR_Laels1)), aes(x=Var1, y = value, fill=Var1))  + 
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR predicted classification  labels #1:: TestisII \n Total Contribution") + 
      ylab("Number of cells"),
    
    
    ggplot(melt(table(SeurObj$SingleR_Laels2)), aes(x=Var1, y = value, fill=Var1))  + 
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR predicted classification labels #2:: TestisII \n Total Contribution") + 
      ylab("Number of cells"),
    ncol = 1)
  
}




