#' @import Seurat

# Note: this file contains methods that are essentially copy/paste from Seurat3 to add fixes.  We should minimize these and
# remove them from here once Seurat itself is fixed.

# Eisa: do we still need this, or is Seurat3 fixed?
#' @title Cell Cycle Scoring function.
#'
#' @description Old name: CellCycleScoring.
#' @param s.features, A vector of genes associated to S-phase.
#' @param g2m.features, A vector of genes associated to G2M-phases.
#' @param set.ident, Boolean, if T sets cell cycle as identity.
#' @return A modified Seurat object.
#' @keywords SerIII, cellcycle, AddModuleScore
#' @export
CellCycleScoring <- function (object,
s.features,
g2m.features,
set.ident = FALSE) {
    enrich.name <- 'Cell Cycle'
    genes.list <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
    object <- AddModuleScoreAvg(
    SeurObj = object,
    genes.list = genes.list,
    enrich.name = enrich.name,
    ctrl.size = min(vapply(X = genes.list, FUN = length, FUN.VALUE = numeric(1)))
    )

    cc.columns <- grep(pattern = enrich.name, x = colnames(x = object@meta.data))
    cc.scores <- object@meta.data[, cc.columns]

    gc(verbose = FALSE)
    assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
        if (all(scores < 0)) {
            return(null)
        } else {
            return(c(first, second)[which(x = scores == max(scores))])
        }
    }
    )
    cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
    colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
    rownames(x = cc.scores) <- cc.scores$rownames
    cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]

    object$S.Score <- cc.scores$S.Score
    object$G2M.Score <- cc.scores$G2M.Score
    object$Phase <- cc.scores$Phase

    if (set.ident) {
        object$old.or.idents <- Idents(object = object)
        Idents(object) <- cc.scores$Phase
    }
    return(object)
}


# Eisa: Do we need this, or can we call Seurat:::LengthCheck?
#' @title A Title
#'
#' @description A description
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
LengthCheck <- function(values, cutoff = 0) {
    #Seurat v2 fx
    # Check the length of components of a list
    #
    # @param values A list whose components should be checked
    # @param cutoff A minimum value to check for
    #
    # @return a vector of logicals
    #
    return(vapply(
    X = values,
    FUN = function(x) {
        return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
    ))
}


#' @title AddModuleScoreAvg
#'
#' @description Instead of the status quo score of Seurat which is 1 score 1 gene, this function takes a list of genes and computes per set, the average of the individual scores.
#' @param object, A Seurat object.
#' @param genes.list, Gene list to obtain a score for
#' @param genes.pool, Gene list to base as the pool; NULL = all.
#' @param n.bin, number of bins to evaluate score across; default 25.
#' @param seed.use, seed use.
#' @param ctrl.size, control gene set size.
#' @param enrich.name, A name for the assesment
#' @param random.seed, random seed
#' @return A modified Seurat object.
#' @keywords SerIII_AddModuleScore
#' @export
#' @importFrom Hmisc cut2
#' @importFrom Matrix colMeans rowMeans
AddModuleScoreAvg <- function(
#May-2019 version

#this is a modified version of the AddModuleScore
#returnScore = F/T controls the output.
#if T, just the scores are returned,
#if F, the scores are put in the Seurat obj and the Seurat SeurObj is returned.
#Also this FX is modified to work for Seurat V3
SeurObj,
genes.list = NULL,
genes.pool = NULL,
n.bin = 25,
seed.use = 1,
ctrl.size = 100,
enrich.name = "Cluster",
random.seed = 1, returnScore = F) {


    set.seed(seed = random.seed)
    genes.old <- genes.list


    if (is.null(x = genes.list)) {
        stop("Missing input gene list")
    }

    genes.list <- lapply(
    X = genes.list,
    FUN = function(x) {
        return(intersect(x = x, y = rownames(SeurObj)))
    }
    )

    cluster.length <- length(x = genes.list)

    if (!all(Seurat:::LengthCheck(values = genes.list))) {
        warning(paste(
        'Could not find enough genes in the SeurObj from the following gene lists:',
        paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
        'Attempting to match case...'
        ))

        genes.list <- lapply(
        X = genes.old,
        FUN = CaseMatch, match = rownames(SeurObj)
        )
    }

    if (!all(Seurat:::LengthCheck(values = genes.list))) {
        stop(paste(
        'The following gene lists do not have enough genes present in the SeurObj:',
        paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
        'exiting...'
        ))
    }
    if (is.null(x = genes.pool)) {
        genes.pool = rownames(SeurObj)
    }
    data.avg <- Matrix::rowMeans(SeurObj@assays$RNA@data[genes.pool, ])
    data.avg <- data.avg[order(data.avg)]
    data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
    ))
    names(x = data.cut) <- names(x = data.avg)
    ctrl.use <- vector(mode = "list", length = cluster.length)
    for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        for (j in 1:length(x = genes.use)) {
            ctrl.use[[i]] <- c(
            ctrl.use[[i]],
            names(x = sample(
            x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
            size = ctrl.size,
            replace = FALSE
            ))
            )
        }
    }

    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
    ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = SeurObj@assays$RNA@data)
    )

    for (i in 1:length(ctrl.use)) {
        genes.use <- ctrl.use[[i]]
        ctrl.scores[i, ] <- Matrix::colMeans(x = SeurObj@assays$RNA@data[genes.use, ])
    }
    genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = SeurObj@assays$RNA@data)
    )
    for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        data.use <- SeurObj@assays$RNA@data[genes.use, , drop = FALSE]
        genes.scores[i, ] <- Matrix::colMeans(x = data.use)
    }
    genes.scores.use <- genes.scores - ctrl.scores
    rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
    genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
    rownames(x = genes.scores.use) <- colnames(x = SeurObj@assays$RNA@data)

    for (colName in colnames(genes.scores.use)) {
        SeurObj[[colName]] <- genes.scores.use[colnames(SeurObj), colName]
    }

    gc(verbose = FALSE)

    if(!returnScore){

        return(SeurObj)

    } else {
        SeurObj@meta.data$cID <- rownames(SeurObj@meta.data)

        return(SeurObj@meta.data[, c("cID", colnames(genes.scores.use))] )
    }
}


#' @title A Title
#'
#' @description A description
#' @param object, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @importFrom fitdistrplus fitdist
#' @importFrom cluster clara
HTODemux2 <- function(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  verbose = TRUE
) {
  # This is a hack around Seurat's method.  If this is improved, shift to use that:

  #initial clustering
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = GetAssayData(object = object, assay = assay)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
        x = t(x = GetAssayData(object = object, assay = assay)),
        k = ncenters,
        samples = nsamples
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    },
    stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
  )
  #average hto signals per cluster
  #work around so we don't average all the RNA levels which takes time
  average.expression <- AverageExpression(
    object = object,
    assays = c(assay),
    verbose = FALSE
  )[[assay]]

  #TODO: checking for any HTO negative in all clusters:

  #if (sum(average.expression == 0) > 0) {
  #  stop("Cells with zero counts exist as a cluster.")
  #}

  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]

    minNonZero <- which.min(x = average.expression[iter,average.expression[iter, ] > 0])
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[[minNonZero]]
    )]
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    }
  }
  # now assign cells to HTO based on discretized values
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
    }
  )])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
  )
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )
  colnames(x = classification.metadata) <- paste(
    assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
    sep = '_'
  )
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  # object@meta.data$hash.ID <- Idents(object)
  object$hash.ID <- Idents(object = object)
  return(object)
}
