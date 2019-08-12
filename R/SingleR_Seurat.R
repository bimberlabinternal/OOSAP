


#' @title generate SingleR object
#'
#' @description Compute SingleR classification on a Seurat object
#' @param seuratObj A Seurat object
#' @param dataset The dataset (see singleR docs) to use as a reference
#' @param assay The assay in the seurat object to use
#' @return SingleR object
#' @keywords SingleR Classification
#' @import Seurat
#' @import SingleR
GenerateSingleR <- function(seuratObj = NULL, dataset = 'hpca', assay = NULL){
    data <- GetAssayData(object = seuratObj, slot = 'counts', assay = assay)
    ref <- SingleR::getReferenceDataset(dataset)
    pred <- SingleR(test=data, ref=ref$data, labels=ref$types)

    return(pred)
}


#' @title SingleR my Seurat Object
#'
#' @description Compute SingleR classification on a Seurat object from path or direcly
#' @param seuratObj A Seurat object
#' @param dataset The dataset (see singleR docs) to use as a reference
#' @param assay The assay in the seurat object to use
#' @param singleRSavePath If provided, the singleR object will be saved here
#' @return The modified Seurat object
#' @keywords SingleR Classification
#' @export
#' @import Seurat
#' @import SingleR
SingleRmySerObj <- function(seuratObj = NULL, dataset = 'hpca', assay = NULL, singleRSavePath = NULL) {
    if (is.null(seuratObj)){
        stop("Seurat object is required")
    }
  
    if (!is.null(singleRSavePath)) {
        if (!dir.exists(dirname(singleRSavePath))) {
            stop("Save path does not exist")
        }
    }

    singler <- GenerateSingleR(seuratObj = seuratObj, dataset = dataset)

    if (!is.null(singleRSavePath)) {
        saveRDS(singler, singleRSavePath)
    }

    seuratObj <- SaveSingleRtoSeuratObj(seuratObj, singler)

    return(seuratObj)
}


#' @title Save SingleR to SeurObject
#'
#' @description Add a pre-Computed SingleR classification into a Seurat object from path or direcly
#' @return modified Seurat object
#' @keywords SingleR Classification
#' @import Seurat
SaveSingleRtoSeuratObj <- function(seuratObj, singler) {
  seuratObj$SingleR_Labels1 <- "unk"
  seuratObj@meta.data[names(singler$singler[[1]]$SingleR.single$labels[,1]),]$SingleR_Labels1 <- singler$singler[[1]]$SingleR.single$labels[,1]
  
  seuratObj$SingleR_Labels2 <- "unk"
  seuratObj@meta.data[names(singler$singler[[2]]$SingleR.single$labels[,1]),]$SingleR_Labels2 <- singler$singler[[2]]$SingleR.single$labels[,1]
  
  seuratObj$SingleR_LabelsOther <- "unk"
  seuratObj@meta.data[names(singler$other),]$SingleR_LabelsOther <- singler$other
  
  return(seuratObj)
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
      ncol = 1
  )
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




