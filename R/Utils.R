#' @importFrom grDevices colorRampPalette colors hcl recordPlot
#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table



#' @title Range01
#'
#' @description Scales the range of input to be between 0 and 1
#' @return A numeric vector
#' @keywords range
range01 <- function(x, MaxN = NULL, MinN = NULL){
  if(is.null(MaxN)) MaxN = max(x)
  if(is.null(MinN)) MinN = min(x)
  
  (x - MinN)/(MaxN - MinN)
}

#' @title Range01b
#'
#' @description Scales the range of input to be between 0 and 1, w/ a prior Min Max 
#' @return A numeric vector
#' @keywords range
range01b <- function(x, MaxN = 10, MinN = -10){
  range01(x, MaxN = MaxN, MinN = MinN)
  
}


utils::globalVariables(
  names = c('MDS1', 'MDS2', 'geom_text_repel'),
  package = 'OOSAP',
  add = TRUE
)

#' @title MDSmyDF
#'
#' @description Compute Multidimensional scaling on a dataframe
#' @return MDS plot of data or data
#' @keywords MDS, ggplot2
MDSmyDF <- function(dfx, labelsDF, factorV, title = "MDS Plot", col_vector, returnMDS = F){
  
  
  if(length(factorV) != nrow(dfx)) stop("rotate t() your dataX")
  
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
    p1 <- p1 + geom_text_repel(aes(y = MDS2 + 0.25), label = factorV) +
      ggtitle(paste("MDS of:",title ,sep=" "))+
      theme(plot.title = element_text(hjust = 0.5))
    
  
  if(returnMDS) return(mds) else p1
  
  
}

#' @title PCAmyDF
#'
#' @description Compute PCA on a dataframe
#' @return PCA plot of data or data
#' @keywords PCA, ggplot2
#' @import ggplot2
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

#' @title transposedt
#'
#' @description Transpose a data.table
#' @return A transposed data.table
#' @keywords transpose, t
#' @param dt The datatable
#' @import data.table
TransposeDT <- function(dt, varlabel="myVar") {
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}



#' @title LogAdd
#'
#' @description Sum a vector of log valued
#' @return A vector.
#' @keywords log, sum
LogAdd <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}

#' @title SetIfNull
#'
#' @description checks is.null, used in various Seurat based functions.
#' @return A vector x
#' @keywords is.null, Seurat
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}


#' @title UniformSampleDF_FacPor
#'
#' @description uniformly sample a dataframe based on factor and porportion
#' @return index of the uniformly sampled samples
#' @keywords sample, uniform
UniformSampleDF_FacPor <- function(x, ClassF, p){
  nr <- NROW(x)
  size <- (nr * p) %/% length(unique(ClassF))
  idx <- lapply(split(seq_len(nr), ClassF), function(.x) sample(.x, size))
  unlist(idx)
}

#' @title WichIn1not2
#'
#' @description compares unique association between two labels from A dataftrame usually from GO annotation with a column named cluster
#' @return unique genes to the comparison
#' @keywords cluster, gene, 
WichIn1not2 <- function(Clus1N = c(1), DataT = "", Clus2N = c(2)){
  Gs1 <- subset(DataT, cluster %in% Clus1N)$gene 
  Gs2 <- subset(DataT, cluster %in% Clus2N)$gene
  Gs1[which(!Gs1 %in% Gs2)]
  
}

#' @title quickTabulate
#' @description tablulates a matrix
#' @param spMat, A sparse matrix
#' @return histo_numers
quickTabulate <- function(spMat){
  histo_numers <- matrix(c(0:max(spMat), rep(0, max(spMat)+1)), ncol = 2)
  histo_numers[1:max(spMat)+1, 2] <- tabulate(as.matrix(spMat))
  histo_numers[1, 2] <- sum(spMat == 0)
  return(histo_numers)
}

#' @title is.even
#'
#' @description logical returns T if even
#' @param x, numbers
#' @export
is.even <- function(x) x %% 2 == 0

#' @title is.odd
#'
#' @description logical returns T if odd
#' @param x, numbers.
#' @return histo_numers
#' @export
is.odd <- function(x) x %% 2 != 0


#' @title PlotAvgExpr
#'
#' @description Plots average expression of each cluster/group for all associated cells
#' @param x, numbers.
#' @return histo_numers
#' @keywords 
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
  
  
  ggplot(avg.combo.cells, aes(X, Y)) + geom_point() + 
    geom_text(aes(label=gene3), size=3, colour=HighColor,
              vjust=0, hjust=-0.1) +
    ggtitle(Title) + xlab(Xlab) + ylab(Ylab) + 
    theme_bw() +
    geom_point(data=subset(avg.combo.cells, gene2 == "show"), aes(x=X, y=Y), colour="dodgerblue", size=2)
  
  
}



#' @title RunUMAP.Matrix
#' @description Return UMAP 2D with similar parameters as Seurat
#' @param DGEmat, A matrix rows are cells
#' @return 2D UMAP rows are cellss
#' @keywords 
#' @export
#' @import reticulate
RunUMAP.Matrix <- function(
  #originally from Seurat pacakge, 
  DGEmat,
  assay = NULL,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  reduction.key = 'UMAP_',
  verbose = TRUE,
  ...
) {
  
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  umap_output <- umap$fit_transform(as.matrix(x = DGEmat))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- rownames(DGEmat)
  
  return(umap_output)
}