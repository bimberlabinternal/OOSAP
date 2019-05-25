
#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
Garnett_Classification_Seurat <- function(SerObj,
                                          marker_file_path = "./data/Garnett/pbmc_classification.txt",
                                          reutrnMonObj=F, returnTrainedClassifier=T){
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  library(monocle)
  library(garnett)
  library(org.Hs.eg.db)

  tempStorageLS <- list()

  tempStorageLS$MonComboObj <- importCDS(SerObj, import_all = TRUE)


  tempStorageLS$MonComboObj <- estimateSizeFactors(tempStorageLS$MonComboObj)

  head(pData(tempStorageLS$MonComboObj))



  marker_check <- check_markers(tempStorageLS$MonComboObj,
                                marker_file_path,
                                db=org.Hs.eg.db,
                                cds_gene_id_type = "SYMBOL",
                                marker_file_gene_id_type = "SYMBOL")

  tempStorageLS$Garnett_classifier <- train_cell_classifier(cds = tempStorageLS$MonComboObj,
                                                            marker_file = marker_file_path,
                                                            db=org.Hs.eg.db,
                                                            cds_gene_id_type = "SYMBOL",
                                                            num_unknown = 50,
                                                            marker_file_gene_id_type = "SYMBOL")



  tempStorageLS$MonComboObj <- classify_cells(tempStorageLS$MonComboObj,
                                              tempStorageLS$Garnett_classifier,
                                              db = org.Hs.eg.db,
                                              cluster_extend = TRUE,
                                              cds_gene_id_type = "SYMBOL")
  tempStorageLS$return <- list()
  tempStorageLS$return$Meta.Data <- pData(tempStorageLS$MonComboObj)
  if(reutrnMonObj) tempStorageLS$return$MonObj <- tempStorageLS$MonComboObj else tempStorageLS$return$MonObj <- "notSaved"
  if(returnTrainedClassifier) tempStorageLS$return$Garnett_classifier <- tempStorageLS$Garnett_classifier else tempStorageLS$return$Garnett_classifier <- "notSaved"
  return(tempStorageLS$return)

}

