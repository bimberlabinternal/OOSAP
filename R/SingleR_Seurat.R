#' @title generate SingleR object
#'
#' @description Compute SingleR classification on a Seurat object
#' @param seuratObj A Seurat object
#' @param dataset The dataset (see singleR docs) to use as a reference
#' @param assay The assay in the seurat object to use
#' @param resultTableFile If provided, a table of results will be saved here
#' @return The modified seurat object
#' @keywords SingleR Classification
#' @import Seurat
#' @import SingleR
#' @export
#' @importFrom scater logNormCounts
RunSingleR <- function(seuratObj = NULL, dataset = 'hpca', assay = NULL, resultTableFile = NULL){
    if (is.null(seuratObj)){
        stop("Seurat object is required")
    }

    sce <- Seurat::as.SingleCellExperiment(seuratObj)

    if (dataset == 'hpca'){
        ref <- SingleR::HumanPrimaryCellAtlasData()
    } else {
        stop('hpca is currently the only supported reference dataset')
    }

    #Subset genes:
    genesPresent <- intersect(rownames(seuratObj), rownames(ref))
    ref <- ref[genesPresent,]
    seuratObj <- seuratObj[genesPresent,]

    #Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(seuratObj)
    sce <- scater:::logNormCounts(sce)

    refAssay <- 'logcounts'
    if (!('logcounts' %in% names(SummarizedExperiment::assays(ref)))) {
        refAssay <- 'normcounts'
    }
    pred.results <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main, method = 'single', assay.type.ref = refAssay)
    pred.results$labels[is.na(pred.results$labels)] <- 'Unknown'

    if (sum(colnames(seuratObj) != rownames(pred.results)) > 0) {
        stop('Cell barcodes did not match for all results')
    }

    seuratObj[['SingleR_Labels']] <- pred.results$labels

    pred.results <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, method = 'single', assay.type.ref = refAssay)
    pred.results$labels[is.na(pred.results$labels)] <- 'Unknown'

    if (sum(colnames(seuratObj) != rownames(pred.results)) > 0) {
        stop('Cell barcodes did not match for all results')
    }

    seuratObj[['SingleR_Labels_Fine']] <- pred.results$labels

    #sanity check:
    if (length(colnames(seuratObj)) != length(rownames(pred.results))) {
        stop('SingleR did not produce results for all cells')
    }

    if (!is.null(resultTableFile)){
        write.table(file = resultTableFile, data.frame(CellBarcodes = rownames(pred.results), SingleR_Labels = seuratObj$SingleR_Labels, SingleR_Labels_Fine = seuratObj$SingleR_Labels_Fine), sep = '\t', row.names = F, quote = F)
    }

    return(seuratObj)
}


#' @title DimPlot SingleR Class Lables
#' @description Dimplot Seurobject with SingleR class lables
#' @param seuratObject a Seurat object, but if path given, path is prioritized.
#' @param plotIndividually If true, two separate plots will be printed.  Otherwise a single plot wil be printed with one above the other
#' @keywords Dimplot SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
DimPlot_SingleRClassLabs <- function(seuratObject, plotIndividually = F){
    plots <- list(
      DimPlot(seuratObject, group.by = "SingleR_Labels") + theme_bw() + ggtitle("SingleR Predicted Classification") + theme(legend.position="bottom"),
      DimPlot(seuratObject, group.by = "SingleR_Labels_Fine") + theme_bw() + ggtitle("SingleR Predicted Classification (Fine)") + theme(legend.position="bottom")
    )

    if (plotIndividually){
        plot(plots[[1]])
        plot(plots[[2]])
    } else {
        cowplot::plot_grid(plots[[1]], plots[[2]], ncol = 1)
    }
}


#' @title Tabulate SingleR Class Lables
#' @description Tabulate Seurobject with SingleR class lables
#' @param seuratObject a Seurat object, but if path given, path is prioritized.
#' @param plotIndividually If true, two separate plots will be printed.  Otherwise a single plot wil be printed with one above the other
#' @keywords Tabulate SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
Tabulate_SingleRClassLabs <- function(seuratObject, plotIndividually = F){
  plots <- list(
    ggplot(melt(table(seuratObject$SingleR_Labels)), aes(x=Var1, y = value, fill=Var1))  +
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR Predicted Classification:") +
      ylab("Number of cells"),

    ggplot(melt(table(seuratObject$SingleR_Labels_Fine)), aes(x=Var1, y = value, fill=Var1))  +
      geom_bar(stat="identity", position="dodge", width = 0.7) + 
      # scale_fill_manual(values=col_vector) +
      theme_bw() + 
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle("SingleR Predicted Classification (Fine):") +
      ylab("Number of cells")
    )

    if (plotIndividually){
        plot(plots[[1]])
        plot(plots[[2]])
    } else {
        cowplot::plot_grid(plots[[1]], plots[[2]], ncol = 1)
    }
}




