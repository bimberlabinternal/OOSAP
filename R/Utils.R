#' @importFrom grDevices colorRampPalette colors hcl recordPlot
#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table

#' @title frac_to_numeric
#' @description Converts fraction to numeric
#' @return A vector
#' @keywords numeric
frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))



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


utils::globalVariables(
  names = c('X', 'Y', 'gene2', 'gene3'),
  package = 'OOSAP',
  add = TRUE
)





#' @title RunUMAP.Matrix
#' @description Return UMAP 2D with similar parameters as Seurat
#' @param DGEmat, A matrix rows are cells
#' @return 2D UMAP rows are cellss
#' @export
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

.GetCCGenes <- function(){
  # Cell cycle genes were obtained from the Seurat example (See regev_lab_cell_cycle_genes.txt)
  # and stored using use_data(internal = T) (https://github.com/r-lib/usethis and use_data)
  # cc.genes
  # g2m.genes.orig

  return(cc.genes)
}

.GetSPhaseGenes <- function(){
  return (.GetCCGenes()[1:43])
}

.GetG2MGenes <- function() {
  return(unique(c(g2m.genes.orig, .GetCCGenes()[44:97])))
}


.WriteLogMsg <- function(msg, prefixTime = TRUE, file = 'OOSAP.log.txt') {
  if (prefixTime) {
    msg <- paste0(Sys.time(), ' ', msg)
  }

  write(msg, file = file, append = T)
}


#' @title find_peaks
#' @description Returns the maxima of points. To get minima, -1*x
#' @param x, A vector of numbers, if small or not very smooth, use a smoothing function. Like density(x).
#' @param m, An integer that acts as a loose hand for resolution.
#' @return vector of peaks positions.
find_peaks <- function (x, m = 4){
  #https://github.com/stas-g/findPeaks
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


.InferPerplexityFromSeuratObj <- function(seuratObj, perplexity = 30) {
  return(.InferPerplexity(ncol(seuratObj), perplexity))
}


.InferPerplexity <- function(sampleNumber, perplexity = 30) {
  if (sampleNumber - 1 < 3 * perplexity) {
    print(paste0('Perplexity is too large for the number of samples: ', sampleNumber))
    perplexityNew <- floor((sampleNumber - 1) / 3)
    print(paste0('lowering from ', perplexity, ' to: ', perplexityNew))
    perplexity <- perplexityNew
  }

  return(perplexity)
}

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @export
#'
#' @import Rlabkey
DownloadOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL) {
	return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = F))
}

#' @title DownloadOutputFile
#' @description Downloads an output file tracked in LabKey to the local filesystem.
#' @param outputFileId The rowid of the sequence outputfiles on the webserver
#' @param outFile The local path to write this file
#' @param overwrite If true, any pre-existing local copy will be replaced.
#' @param pathTranslator There are instances where the goal is to download some file located relative to the actual outputfile. For example, 10x data has an outputfile for the Loupe file, but it might be useful to download the raw counts. To accomplish this, a function can be provided that accepts the primary file remote path, and returns a modified path to download.
#' @export
#'
#' @import Rlabkey
DownloadOutputDirectoryFromOutputFile <- function(outputFileId, outFile, overwrite = T, pathTranslator = NULL) {
  if (is.null(pathTranslator)) {
    stop('When attempting to download a folder relative to an outputfile, you must provide a pathTranslator function')
  }
  return(.DownloadOutputFileOrDirectory(outputFileId = outputFileId, outFile = outFile, overwrite = overwrite, pathTranslator = pathTranslator, asDirectory = T))
}

.DownloadOutputFileOrDirectory <- function(outputFileId, outFile, asDirectory, overwrite = T, pathTranslator = NULL) {
  if (is.na(outputFileId)) {
    stop('Output file ID cannot be NA')
  }

  if (file.exists(outFile) & !overwrite) {
    print(paste0("File exists, will not overwrite: ", outFile))
    return(outFile)
  }

  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
    colFilter=makeFilter(c("rowid", "EQUAL", outputFileId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  if (nrow(rows) != 1) {
    stop(paste0('More than one matching file found, this should not occur.  RowId: ', outputFileId))
  }

  wb <- rows[['workbook_workbookid']]
  if (is.na(wb) || is.null(wb)){
    wb <- ''
  }

  remotePath <- rows[['dataid_webdavurlrelative']]
  if (!is.null(pathTranslator)) {
    remotePath <- pathTranslator(remotePath)
  }

	if (asDirectory) {
		success <- labkey.webdav.downloadFolder(
			baseUrl=lkBaseUrl,
			folderPath=paste0(lkDefaultFolder,wb),
			remoteFilePath = remotePath,
			overwrite = overwrite,
			localBaseDir = outFile
		)

		if (!success) {
			stop(paste0('labkey.webdav.downloadFolder failed for file: ', remotePath))
		}
	} else {
		success <- labkey.webdav.get(
			baseUrl=lkBaseUrl,
			folderPath=paste0(lkDefaultFolder,wb),
			remoteFilePath = remotePath,
			overwrite = overwrite,
			localFilePath = outFile
		)

    if (!success || !file.exists(outFile)) {
      stop(paste0('labkey.webdav.get failed for file: ', remotePath))
    }
	}

  return(outFile)
}