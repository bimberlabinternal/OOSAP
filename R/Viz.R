# Functions with primary goal or only output is a plot such as customized viz fx are here.

#' @importFrom grDevices colorRampPalette colors hcl recordPlot
#' @importFrom graphics abline boxplot legend lines plot points segments


#' @title PlotAvgExpr
#'
#' @description Plots average expression of each cluster/group for all associated cells
#' @param x, numbers.
#' @return histo_numers
#' @export
PlotAvgExpr <- function(GenesNames2Show, X_avg, Y_avg, features=NULL, Xlab="Xlab", Ylab="Ylab", Title = "Title", HighColor = "dodgerblue"){
  
  X_avg$gene <- rownames(X_avg)
  Y_avg$gene <- rownames(Y_avg)
  
  if(is.null(features)) features = rownames(Y_avg)
  
  avg.combo.cells <- merge(X_avg[features,], Y_avg[features,], by = "gene")
  
  colnames(avg.combo.cells) <- c("gene", "X", "Y")
  
  
  avg.combo.cells$gene3 <- avg.combo.cells$gene
  
  avg.combo.cells$gene2 <- ifelse(avg.combo.cells$gene %in% GenesNames2Show, "show", "hide")
  
  avg.combo.cells[which(avg.combo.cells$gene2=="show"),]$gene3 <- avg.combo.cells[which(avg.combo.cells$gene2=="show"),]$gene
  
  avg.combo.cells[which(avg.combo.cells$gene2=="hide"),]$gene3 <- NA
  
  
  p1 <- ggplot(avg.combo.cells, aes(X, Y)) + geom_point() + 
    geom_text(aes(label=gene3), size=3, colour=HighColor,
              vjust=0, hjust=-0.1) +
    ggtitle(Title) + xlab(Xlab) + ylab(Ylab) + 
    theme_bw() +
    geom_point(data=subset(avg.combo.cells, gene2 == "show"), aes(x=X, y=Y), colour="dodgerblue", size=2)
  
  print(p1)
  
  
}

#' @title PlotMyTable
#'
#' @description Plots average expression of each cluster/group for all associated cells
#' @param MyTable, A 2-way table e.g. MyTable = table(x, y)
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotMyTable <- function(MyTable, Title="", legend.position="bottom", PlotCombo = F){
  tempMeltDF <- melt(MyTable)
  
  
  p1 <- (ggplot(tempMeltDF) +
          geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
          theme_bw()  + scale_fill_manual(values=col_vector) +
          theme(legend.position=legend.position,
                legend.direction="horizontal",
                legend.title = element_blank()) +
          ggtitle(paste0(Title, "\n Total Contribution")) + ylab("Total No."))
  
  
  p2<- (ggplot(tempMeltDF) +
          geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
          theme_bw()  + scale_fill_manual(values=col_vector) +
          theme(legend.position=legend.position,
                legend.direction="horizontal",
                legend.title = element_blank())+
          scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                             labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
          ggtitle(paste0(Title, "\n Relative % Contribution")) + ylab("(%)"))
  
  if(PlotCombo) print(cowplot::plot_grid(p1, p2, ncol = 1)) else {
    print(p1)
    print(p2)
  }
  
}

#' @title PlotBin2D
#'
#' @description Plots a set of plots with geom_bin2d() split by factor levels
#' @param DF, A dataframe with at least 3 columns, 2 for 2D scatter and one for the levels Named V1, V2, and Fac.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotBin2D <- function(DF, bins = 75, ncol = 3){
  
  if(!("Fac" %in% colnames(DF))) stop("FacName not found in DF")
  if(!("V1" %in% colnames(DF))) stop("V1 not found in DF")
  if(!("V2" %in% colnames(DF))) stop("V2 not found in DF")
  
  
  if(!is.factor(DF$Fac)) DF$Fac <- factor(DF$Fac)
  
  p1 <- cowplot::plot_grid(plotlist = lapply(levels(DF$Fac), function(xN){
    
    ggplot(subset(DF, Fac == xN), aes(V1, V2)) + geom_bin2d(bins = bins) + theme_bw() +
      theme(legend.position = "bottom") +
      ggtitle(paste0(xN)) +
      scale_fill_viridis() + guides(colour = guide_legend(override.aes = list(size=2, alpha=1)))
    
  }),
  ncol = ncol)
  
  print(p1)
  
}
                      



