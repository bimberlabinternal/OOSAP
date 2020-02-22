#' @import Seurat

# Note: this file contains methods that are essentially copy/paste from Seurat3 to add fixes.  We should minimize these and
# remove them from here once Seurat itself is fixed.




# This method looks different than the Seurat::HTODemux so keep for now

#' @title HTODemux2
#'
#' @description A description
#' @param object, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @importFrom fitdistrplus fitdist
#' @importFrom cluster clara
#' @importFrom Matrix t
HTODemux2 <- function(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  verbose = TRUE,
  plotDist = FALSE
) {
  # This is a hack around Seurat's method to make it more robust and fault tolerant.  If this is improved, shift to use that:
  slot <- "counts"

  #initial clustering
  data <- GetAssayData(object = object, assay = assay)
  assayData <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot
  )[, colnames(x = object)]

  ncenters <- (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = data),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
      	x = t(x = data),
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
    slot = slot,
    verbose = FALSE
  )[[assay]]

  if (sum(average.expression == 0) > 0) {
    warning("Cells with zero counts exist as a cluster.")
  }

  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay, slot = slot)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (hto in naturalsort::naturalsort(rownames(x = data))) {
    values <- assayData[hto, colnames(object)]

    minNonZero <- which.min(x = average.expression[hto,average.expression[hto, ] > 0])
    maxNonZero <- which.max(x = average.expression[hto,average.expression[hto, ] > 0])

    # This indicates the only non-zero cluster is the primary one.
    # It's not especially likely, but could occur when there are only 2 possible HTOs
    if (minNonZero == maxNonZero){
      print('Min. non-zero value is the same as max value, using cutoff of 1 read')
      cutoff <- 1
    } else {
      values.use <- values[WhichCells(
        object = object,
        idents = levels(x = Idents(object = object))[[minNonZero]]
      )]

      doSkip <- FALSE
      tryCatch({
        fit <- suppressWarnings(fitdist(data = values.use, distr = "nbinom"))
        if (plotDist) {
          print(plot(fit))
        }

        cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
      }, error = function(e) {
        print(paste0('Error fitting nbinom, skipping: ', hto))
        print(e)
        saveRDS(values.use, file = paste0('./', hto, '.nbinom.data.rds'))
        doSkip <- TRUE
      })

      if (doSkip) {
        next
      }
    }

    discrete[hto, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      print(paste0("Cutoff for ", hto, " : ", cutoff, " reads"))
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

  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'

  object$hash.ID <- Idents(object = object)
  return(object)
}
