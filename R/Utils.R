#' @importFrom grDevices colorRampPalette colors hcl recordPlot
#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist filter kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
range01 <- function(x, MaxN = NULL, MinN = NULL){
  if(is.null(MaxN)) MaxN = max(x)
  if(is.null(MinN)) MinN = min(x)
  
  (x - MinN)/(MaxN - MinN)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
range01b <- function(x, MaxN = 10, MinN = -10){
  if(is.null(MaxN)) MaxN = max(x)
  if(is.null(MinN)) MinN = min(x)
  
  (x - MinN)/(MaxN - MinN)
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
MDSmyDF <- function(dfx, labelsDF, factorV, title = "MDS Plot", col_vector){
  
  
  if(length(factorV) == nrow(dfx)){
    mds <- cmdscale(as.matrix(dist(dfx)))
    colnames(mds) <- c("MDS1", "MDS2")
    
    mds <- cbind(mds, labelsDF) #samples as rows add a cols of labels
    
    
    p1 <- ggplot(mds, aes(x = MDS1, y = MDS2)) +
      theme_bw() +
      geom_hline(yintercept = 0, color = "gray70") +
      geom_vline(xintercept = 0, color = "gray70") +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      coord_cartesian(xlim = c(min(mds[,1])-5, max(mds[,1])+5)) +
      scale_color_manual(values = col_vector, name="Samples")
    
    # the graphic with ggrepel
    p1 + geom_text_repel(aes(y = MDS2 + 0.25), label = factorV) +
      ggtitle(paste("MDS of:",title ,sep=" "))+
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    print("rotate t() your dataX")
  }
  
  
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
PCAmyDF <- function (dfx, labels, factorV, title = "PCA Plot", scale, center, col_vector, namePointBL = F) {
  if(class(labels) == "function") {
    print("no labels, using factor as names")
    labels = as.character(factorV)
  }
  if(length(as.character(factorV)) != length(labels)) {
    print("labels length != factor length, using factor as names")
    labels = as.character(factorV)
  }
  
  dfx.pca <- prcomp(t(dfx), scale.=scale, center = center)
  
  MinXY <- min(c(round(min(dfx.pca$rotation[,1:2]) - abs(min(dfx.pca$rotation[,1:2])*0.5),1), -1) )
  MaxXY <- max(c(round(max(dfx.pca$rotation[,1:2]) + abs(max(dfx.pca$rotation[,1:2])*0.5),1),  1) )
  
  
  if(namePointBL){
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      scale_color_manual(values = col_vector, name="Samples") +
      geom_text_repel(aes(y = PC2, label = labels))  +
      ggtitle(paste("PCA of:",title ,sep=" ")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  } else {
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      scale_color_manual(values = col_vector, name="Samples")  +
      ggtitle(paste("PCA of:",title ,sep=" "))+
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  }
  
  
  
  
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @import data.table
transposedt <- function(dt, varlabel="myVar") {
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
LogAdd <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
UniformSampleDF_FacPor <- function(x, ClassF, p){
  nr <- NROW(x)
  size <- (nr * p) %/% length(unique(ClassF))
  idx <- lapply(split(seq_len(nr), ClassF), function(.x) sample(.x, size))
  unlist(idx)
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
WichIn1not2 <- function(Clus1N = c(1), DataT = "", Clus2N = c(2)){
  Gs1 <- subset(DataT, cluster %in% Clus1N)$gene 
  Gs2 <- subset(DataT, cluster %in% Clus2N)$gene
  Gs1[which(!Gs1 %in% Gs2)]
  
}

#' @title quickTabulate
#'
#' @description A description
#' @param sparce.matrix, A Seurat object.
#' @return histo_numers
#' @keywords 
#' @export
quickTabulate <- function(spMat){
  histo_numers <- matrix(c(0:max(spMat), rep(0, max(spMat)+1)), ncol = 2)
  histo_numers[1:max(spMat)+1, 2] <- tabulate(as.matrix(spMat))
  histo_numers[1, 2] <- sum(spMat == 0)
  return(histo_numers)
}

#' @title is.even
#'
#' @description A description
#' @param x, numbers
#' @keywords 
#' @export
is.even <- function(x) x %% 2 == 0

#' @title is.odd
#'
#' @description A description
#' @param x, numbers.
#' @return histo_numers
#' @keywords 
#' @export
is.odd <- function(x) x %% 2 != 0


#' @title PlotAvgExpr
#'
#' @description A description
#' @param x, numbers.
#' @return histo_numers
#' @keywords 
#' @export
PlotAvgExpr <- function(GenesNames2Show, X_avg, Y_avg, features=NULL, Xlab="Xlab", Ylab="Ylab"){
  
  X_avg$gene <- rownames(X_avg)
  Y_avg$gene <- rownames(Y_avg)
  
  if(is.null(features)) features = rownames(Y_avg)
  
  avg.combo.cells <- merge(X_avg[features,], Y_avg[features,], by = "gene")
  
  colnames(avg.combo.cells) <- c("gene", "X", "Y")
  
  
  avg.combo.cells$gene3 <- avg.combo.cells$gene
  
  avg.combo.cells$gene2 <- ifelse(avg.combo.cells$gene %in% GenesNames2Show, "show", "hide")
  
  avg.combo.cells[which(avg.combo.cells$gene2=="show"),]$gene3 <- avg.combo.cells[which(avg.combo.cells$gene2=="show"),]$gene
  
  avg.combo.cells[which(avg.combo.cells$gene2=="hide"),]$gene3 <- NA
  
  
  ggplot(avg.combo.cells, aes(X, Y)) + geom_point() + 
    geom_text(aes(label=gene3), size=3, colour="dodgerblue",
              vjust=0, hjust=-0.1) +
    ggtitle("Day 0 Vs Day 8 LymphAxLN : Tcell_Other") + xlab(Xlab) + ylab(Ylab) + 
    theme_bw() +
    geom_point(data=subset(avg.combo.cells, gene2 == "show"), aes(x=X, y=Y), colour="dodgerblue", size=2)
  
  
}