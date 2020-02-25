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
  positive.quantile = 0.98,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  verbose = TRUE,
  plotDist = FALSE
) {
  # This is a hack around Seurat's method to make it more robust and fault tolerant.  If this is improved, shift to use that:

  #initial clustering
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
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
    slot = 'counts',
    verbose = FALSE
  )[[assay]]

  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (hto in naturalsort::naturalsort(rownames(x = data))) {
    values <- counts[hto, colnames(object)]

    # Take the bottom 2 clusters (top 2 assumed to be HTO and doublet) as background.
    maxPossibleBackgroundCols <- max(nrow(data) - 2, 1)
    numBackgroundCols <- min(2, maxPossibleBackgroundCols)
    backgroundIndices <- order(average.expression[hto, ])[1:numBackgroundCols]

    if (sum(average.expression[hto, backgroundIndices]) == 0) {

      allPossibleBackgroundIndices <- order(average.expression[hto, ])[1:maxPossibleBackgroundCols]
      for (i in 1:maxPossibleBackgroundCols) {
        print('Expanding clusters until non-zero background obtained')
        backgroundIndices <- allPossibleBackgroundIndices[1:i]
        if (sum(average.expression[hto, backgroundIndices]) > 0) {
          break
        }
      }
    }

    if (verbose) {
      print(paste0('Will select bottom ', numBackgroundCols, ' columns as background'))
      print(paste0('Background clusters for ', hto, ': ', paste0(backgroundIndices, collapse = ',')))
    }

    if (sum(average.expression[hto, backgroundIndices]) == 0) {
      #TODO: unclear what to do with this?
      print('The background clusters have zero reads, defaulting to threshold of 1')
      cutoff <- 1
    } else {
      values.use <- values[WhichCells(
        object = object,
        idents = levels(x = Idents(object = object))[backgroundIndices]
      )]

      if (verbose) {
        print(paste0('total cells for background: ', length(values.use)))
      }

      cutoff <- NULL
      tryCatch(expr = {
        fit <- suppressWarnings(fitdist(data = values.use, distr = "nbinom"))
        if (plotDist) {
          print(plot(fit))
        }

        cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
      }, error = function(e) {
        saveRDS(values.use, file = paste0('./', hto, '.fail.nbinom.rds'))
      })

      if (is.null(cutoff)) {
        print(paste0('Skipping HTO due to failure to fit distribution: ', hto))
        next
      }
    }

    if (verbose) {
      print(paste0("Cutoff for ", hto, " : ", cutoff, " reads"))
    }
    discrete[hto, names(x = which(x = values > cutoff))] <- 1
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

#' @rawNamespace import(Matrix, except = c('tail', 'head'))
MULTIseqDemux2 <- function(
  object,
  assay = "HTO",
  quantile = 0.7,
  autoThresh = FALSE,
  maxiter = 5,
  qrange = seq(from = 0.1, to = 0.9, by = 0.05),
  verbose = TRUE
) {
  if (is.na(assay) || is.null(assay)) {
    assay <- DefaultAssay(object = object)
  }
  multi_data_norm <- t(x = GetAssayData(
    object = object,
    slot = "data",
    assay = assay
  ))
  if (autoThresh) {
    iter <- 1
    negatives <- c()
    neg.vector <- c()
    while (iter <= maxiter) {
      # Iterate over q values to find ideal barcode thresholding results by maximizing singlet classifications
      bar.table_sweep.list <- list()
      n <- 0
      for (q in qrange) {
        n <- n + 1
        # Generate list of singlet/doublet/negative classifications across q sweep
        bar.table_sweep.list[[n]] <- ClassifyCells(data = multi_data_norm, q = q)
        names(x = bar.table_sweep.list)[n] <- paste0("q=" , q)
      }

      # Determine which q values results in the highest pSinglet
      threshold.results1 <- Seurat:::FindThresh(call.list = bar.table_sweep.list)
      res_round <- threshold.results1$res
      res.use <- res_round[res_round$Subset == "pSinglet", ]
      q.use <- res.use[which.max(res.use$Proportion),"q"]

      round.calls <- ClassifyCells(data = multi_data_norm, q = q.use)
      #remove negative cells
      neg.cells <- names(x = round.calls)[which(x = round.calls == "Negative")]
      called.cells <- sum(round.calls != "Negative")
      neg.vector <- c(neg.vector, rep(x = "Negative", length(x = neg.cells)))
      negatives <- c(negatives, neg.cells)

      print(ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) +
        geom_line() +
        theme(
          legend.position = "right"
        ) +
        geom_vline(xintercept=q.use, lty=2) +
        scale_color_manual(values=c("red","black","blue")) +
        ggtitle(paste0('Multi-Seq Iteration: ', iter, ', quantile: ', q.use, ', input: ', length(round.calls), ', called: ', called.cells))
      )

      if (length(x = neg.cells) == 0) {
        print(paste0('no negatives, breaking after iteration: ', iter))
        break
      }

      multi_data_norm <- multi_data_norm[-which(x = rownames(x = multi_data_norm) %in% neg.cells), ]
      iter <- iter + 1
    }
    names(x = neg.vector) <- negatives
    demux_result <- c(round.calls,neg.vector)
    demux_result <- demux_result[rownames(x = object[[]])]
  } else{
    demux_result <- ClassifyCells(data = multi_data_norm, q = quantile)
  }
  demux_result <- demux_result[rownames(x = object[[]])]
  object[['MULTI_ID']] <- factor(x = demux_result)
  Idents(object = object) <- "MULTI_ID"
  bcs <- colnames(x = multi_data_norm)
  bc.max <- bcs[apply(X = multi_data_norm, MARGIN = 1, FUN = which.max)]
  bc.second <- bcs[unlist(x = apply(
  X = multi_data_norm,
  MARGIN = 1,
  FUN = function(x) {
    return(which(x == Seurat:::MaxN(x)))
  }
  ))]
  doublet.names <- unlist(x = lapply(
  X = 1:length(x = bc.max),
  FUN = function(x) {
    return(paste(sort(x = c(bc.max[x], bc.second[x])), collapse =  "_"))
  }
  ))
  doublet.id <- which(x = demux_result == "Doublet")
  MULTI_classification <- as.character(object$MULTI_ID)
  MULTI_classification[doublet.id] <- doublet.names[doublet.id]
  object$MULTI_classification <- factor(x = MULTI_classification)
  return(object)
}


#' @importFrom KernSmooth bkde
#' @importFrom stats approxfun quantile
ClassifyCells <- function(data, q) {
  ## Generate Thresholds: Gaussian KDE with bad barcode detection, outlier trimming
  ## local maxima estimation with bad barcode detection, threshold definition and adjustment
  n_cells <- nrow(x = data)
  bc_calls <- vector(mode = "list", length = n_cells)
  n_bc_calls <- numeric(length = n_cells)
  for (i in 1:ncol(x = data)) {
    model <- NULL
    tryCatch(expr = {
      model <- approxfun(x = bkde(x = data[, i], kernel = "normal"))
    }, error = function(e) {
      print(paste0("Unable to fit model for ", rownames(x = data)[i], ", for ", q, "..."))
      saveRDS(data[, i], file = paste0('./', hto, '.fail.approxfun.rds'))
    })

    # This is changed relative to seurat
    if (is.null(x = model)) {
      print(paste0("Unable to fit model for ", rownames(x = data)[i], ", for ", q, ", skipping"))
      next
    }

    x <- seq.int(
      from = quantile(x = data[, i], probs = 0.001),
      to = quantile(x = data[, i], probs = 0.999),
      length.out = 100
    )
    extrema <- Seurat:::LocalMaxima(x = model(x))
    if (length(x = extrema) <= 1) {
      print(paste0("No extrema/threshold found for ", colnames(x = data)[i]))
      next
    }
    low.extremum <- min(extrema)
    high.extremum <- max(extrema)
    thresh <- (x[high.extremum] + x[low.extremum])/2
    ## Account for GKDE noise by adjusting low threshold to most prominent peak
    low.extremae <- extrema[which(x = x[extrema] <= thresh)]
    new.low.extremum <- low.extremae[which.max(x = model(x)[low.extremae])]
    thresh <- quantile(x = c(x[high.extremum], x[new.low.extremum]), probs = q)

    ## Find which cells are above the ith threshold
    cell_i <- which(x = data[, i] >= thresh)
    n <- length(x = cell_i)
    if (n == 0) { ## Skips to next BC if no cells belong to the ith group
      next
    }
    bc <- colnames(x = data)[i]
    if (n == 1) {
      bc_calls[[cell_i]] <- c(bc_calls[[cell_i]], bc)
      n_bc_calls[cell_i] <- n_bc_calls[cell_i] + 1
      if (is.na(n_bc_calls[cell_i])) {
        stop('NA cell, pos 1')
      }
    } else {
      # have to iterate, lame
      for (cell in cell_i) {
        bc_calls[[cell]] <- c(bc_calls[[cell]], bc)
        n_bc_calls[cell] <- n_bc_calls[cell] + 1
        if (is.na(n_bc_calls[cell])) {
          stop('NA cell, pos 2')
        }
      }
    }
  }
  calls <- character(length = n_cells)

  for (i in 1:n_cells) {
    if (is.na(n_bc_calls[i]) || is.null(n_bc_calls[i]) || !is.finite(n_bc_calls[i])){
      print(head(n_bc_calls))
      stop(paste0('bad value: ', n_bc_calls[i], ' ', i))
    }

    if (n_bc_calls[i] == 0) { calls[i] <- "Negative"; next }
    if (n_bc_calls[i] > 1) { calls[i] <- "Doublet"; next }
    if (n_bc_calls[i] == 1) { calls[i] <- bc_calls[[i]] }
  }
  names(x = calls) <- rownames(x = data)
  return(calls)
}
