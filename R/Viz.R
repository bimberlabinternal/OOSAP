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
                      

#' @title PlotExprThrTab
#'
#' @description This function from a seurat object, takes the gene expression of a single gene and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SeuratObj, A Seurat object.
#' @param GeneName, the Gene name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotExprThrTab <- function(SeuratObj, GeneName = NULL, CutThresh = NULL, PlotBar = T,
                       MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                       col_vector = col_vector){
  
  if(is.null(GeneName)) stop("GeneName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(GeneName) != 1) stop ("GeneName Length != 1")
  
  
  # tempTab <- table(SeuratObj@assays$RNA@data[GeneName, ]>CutThresh, SeuratObj@meta.data[,MetaDataName])
  tempTab <- table(GetAssayData(object = SeuratObj, slot = 'data')[GeneName, ]>CutThresh, 
                   SeuratObj@meta.data[,MetaDataName])
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(melt(tempTab)) +
      geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
      theme_bw()  + scale_fill_manual(values=col_vector) +
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90)) +
      ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    
    tempDF  <- as.data.frame(cbind(Expr = SeuratObj@assays$RNA@data[GeneName, ], CellN = colnames(SeuratObj@assays$RNA@data)), stringsAsFactors = F)
    tempDF$Expr <- as.numeric(tempDF$Expr)
    
    vp <- VlnPlot(SeuratObj, features = GeneName, group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=Expr)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .3) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
    
    
    
  }
  
  
  
}


#' @title PlotHistDensOptima
#'
#' @description This function plots histogram+density and finds all optima from a df usually from FetchData() of Seurat. 
#' @param dftemp, a dataframe with column var1 to be plotted. 
#' @param col_vector, color vector.
#' @param nn, a color from color_vector.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotHistDensOptima <- function(dftemp, Print = F, col_vector = col_vector, nn = NULL, Title = ""){
  dens <- density(dftemp$var1, n = 1000)
  peaks <- dens$x[find_peaks(dens$y, m=20)]
  dips <- dens$x[find_peaks(-1*dens$y, m=20)]
  
  if(is.null(nn)) nn=1
  
  
  p1 <- ggplot(dftemp, aes(x=var1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.05) +
    geom_density(alpha=.2, fill=col_vector[nn]) +
    geom_vline(aes(xintercept=mean(var1)), color="blue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=mean(var1) - 2*sd(var1)), color="dodgerblue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=mean(var1) + 2*sd(var1)), color="dodgerblue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=peaks[1] ), color="plum", size=.8) +
    geom_vline(aes(xintercept=dips[1] ), color="red", size=1)  +
    theme_bw() + ggtitle(Title)
  
  
  if((length(peaks) >= 2)) p1 <- p1 + geom_vline(aes(xintercept=peaks[2] ), color="plum", size=.8)
  if((length(peaks) >= 3)) p1 <- p1 + geom_vline(aes(xintercept=peaks[3] ), color="plum", size=.8)
  if((length(peaks) >= 4)) p1 <- p1 + geom_vline(aes(xintercept=peaks[4] ), color="plum", size=.8)
  
  if((length(dips) >= 2)) p1 <- p1 + geom_vline(aes(xintercept=dips[2] ), color="brown", size=.8)
  if((length(dips) >= 3)) p1 <- p1 + geom_vline(aes(xintercept=dips[3] ), color="brown", size=.8)
  if((length(dips) >= 4)) p1 <- p1 + geom_vline(aes(xintercept=dips[4] ), color="brown", size=.8)
  
  if(Print) print(p1) else return(p1)
}


#' @title PlotFeatThrTab
#'
#' @description This function from a seurat object, takes the feature selected and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SeuratObj, A Seurat object.
#' @param FeatName, the Feature name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotFeatThrTab <- function(SeuratObj, FeatName = NULL, CutThresh = NULL, PlotBar = T,
                       MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                       col_vector = col_vector){

  if(is.null(FeatName)) stop("FeatName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(FeatName) != 1) stop ("FeatName Length != 1")
  
  tempDF <- FetchData(SeuratObj, vars = c(FeatName, MetaDataName))
  colnames(tempDF) <- c("var1", "var2")
  tempTab <- table(tempDF$var1>CutThresh, tempDF$var2)
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    

    
    vp <- VlnPlot(SeuratObj, features = FeatName, 
                  group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=var1)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .03) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
    
    
    
  }
  
  
  
}


