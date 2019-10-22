

#' @title ClassifySGS_Multiple
#'
#' @description This will classify cells in a seurat object using a list of gene sets
#' @param seuratObj The seurat object
#' @param moduleScoreGeneLists A named list of character vectors, where each vector is the set of genes to score
#' @param saveFilePath If provided, the datframe will be saved here
#' @return A dataframe with the results per cell, for all gene sets
#' @export
ClassifySGS_Multiple <- function(seuratObj,
	moduleScoreGeneLists = NULL,
	saveFilePath = NULL
) {
	if (is.null(moduleScoreGeneLists)) {
		stop('No gene lists provided')
	}

	if (is.null(names(moduleScoreGeneLists))) {
		stop('moduleScoreGeneLists must be a names list')
	}

	finalResults <- NA
	for (setName in names(moduleScoreGeneLists)) {
		print(paste0('Running set: ', setName))
		genes <- moduleScoreGeneLists[[setName]]
		results <- .DoModuleScoreClassify(seuratObj, columnPrefix = setName, geneList = genes)
		if (!all(is.na(results))) {
			if (all(is.na(finalResults))) {
				finalResults <- results
			} else {
				finalResults <- merge(finalResults, results, by = 'CellBarcode', all.x = T, all.y = T)
			}
		}
	}

	if (!all(is.na(finalResults)) && !is.null(saveFilePath)) {
		write.csv(finalResults, file = saveFilePath, row.names = FALSE)
	}

	return(finalResults)
}


#' @title ClassifySGSAndApply
#'
#' @description This will classify cells in a seurat object using a variant of Seurat's AddModuleScore and save the results to the seurat object
#' @param geneList A character vector of the genes to test
#' @param geneSetName A name to use for this gene list.  This will be used in the header of the dataframe
#' @param saveFilePath If provided, the dataframe will be saved here
#' @param doPlot If true, a QC plot will be created
#' @param positivityThreshold If provided, a boolean column will be added to the dataframe scoring cells as positive or negative, depending on if their score is greater than this threshold
#' @param reduction The reduction (i.e. tsne or umap) that will be used when plotting
#' @return The modified seurat object
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
ClassifySGSAndApply <- function(seuratObj,
geneList = NULL,
geneSetName = NULL,
saveFilePath = NULL,
doPlot = T,
positivityThreshold = NULL,
reduction = NULL
) {
	df <- ClassifySGS(seuratObj, geneList, geneSetName, saveFilePath, doPlot, positivityThreshold)
	if (!all(is.na(df))) {
		if (sum(colnames(seuratObj) != df$CellBarcode) > 0) {
			stop('Cell barcodes do not match between seurat object and table')
		}

		plots <- list()
		for (col in colnames(df)[colnames(df) != 'CellBarcode']) {
			seuratObj[[col]] <- df[col]

			if (doPlot && (all(grepl(pattern = '.Score', x = col)) || all(grepl(pattern = '.Call', x = col)))) {
				plots[[col]] <- Seurat::FeaturePlot(seuratObj, features = c(col), reduction = reduction)
			}
		}

		if (length(plots) > 1) {
			print(cowplot::plot_grid(plotlist = plots, ncol = 2))
		} else if (length(plots) == 1){
			print(plots[[1]])
		}
	}

	return(seuratObj)
}

#' @title ClassifySGS
#'
#' @description This will classify cells in a seurat object using a variant of Seurat's AddModuleScore
#' @param geneList A character vector of the genes to test
#' @param geneSetName A name to use for this gene list.  This will be used in the header of the dataframe
#' @param saveFilePath If provided, the dataframe will be saved here
#' @param doPlot If true, a QC plot will be created
#' @param positivityThreshold If provided, a boolean column will be added to the dataframe scoring cells as positive or negative, depending on if their score is greater than this threshold
#' @return A dataframe with the results per cell
#' @import ggplot2
#' @export
ClassifySGS <- function(seuratObj,
	geneList = NULL,
	geneSetName = NULL,
	saveFilePath = NULL,
	doPlot = T,
	positivityThreshold = NULL
) {

	if (is.null(geneList)) {
		stop('No geneList provided')
	}

	if (is.null(geneSetName)) {
		stop('No geneSetName provided')
	}

	results <- .DoModuleScoreClassify(seuratObj, columnPrefix = geneSetName, geneList = geneList)
	if (!all(is.na(results))) {
		fn <- paste0(geneSetName, '.Score')
		if (!is.null(positivityThreshold)) {
				print(paste0('Scoring cells using threshold: ', positivityThreshold))
				results[paste0(geneSetName, '.Call')] <- results[fn] > positivityThreshold

				print(table(Call = results[paste0(geneSetName, '.Call')]))
		}

		if (doPlot) {
			P1 <- ggplot(results, aes_string(x = fn, y = '..density..')) +
				#geom_histogram(colour = "black", fill = "white", binwidth = .1) +
				geom_density(alpha = .2, fill = "#E69F00") +
				theme(
					legend.position="bottom",
					legend.direction="horizontal",
					legend.title = element_blank(),
					axis.text.x = element_text(angle = 90)
				) +
				theme_bw() +
				ggtitle(paste0('Score: ', geneSetName)) + ylab("Number of Cells") + xlab("Score")

			if (!is.null(positivityThreshold)) {
				P1 <- P1 + geom_vline(aes(xintercept=positivityThreshold), color="blue", linetype="dashed", size=1)
			}

			print(P1)
		}

		if (!is.null(saveFilePath)) {
			write.csv(results, file = saveFilePath, row.names = FALSE)
		}
	}

	return(results)
}


.DoModuleScoreClassify <- function(seuratObj, columnPrefix, geneList) {
	if (length(geneList) == 0) {
		warning("length of gene list is 0")
		return(NA)
	}

	if (!is.character(geneList)) {
		stop(paste0("Must provide a character vector of gene IDs.  was: ", typeof(geneList)))
	}

	#TODO: is there a quantitative way to best estimate n.bin based on number of genes?
	logBinSize <- F
	for (binSize in c(25, 30 , 50)) {
		df <- try(
		CalculateModuleScoreAvg(seuratObj = seuratObj,
			genes.list = geneList,
			n.bin = binSize
		),
		silent = T
		)

		if (class(df) == "try-error") {
			print(paste0("Bin size ", binSize, " was problematic"))
			logBinSize <- T
		} else {
			if (logBinSize) {
				print(paste0("Bin size ", binSize, " worked!"))
			}

			names(df)[names(df) == 'Score'] <- paste0(columnPrefix, '.Score')

			return(df)
		}
	}

	print(paste0('CalculateModuleScoreAvg failed for gene set: ', geneList))

	return(NA)
}


#' @title CalculateModuleScoreAvg
#'
#' @description Instead of the status quo score of Seurat AddModuleScore which is 1 score 1 gene, this function takes a list of genes and computes per set, the average of the individual scores.
#' @param seuratObj, A Seurat object.
#' @param genes.list, A charater vector of genes to score
#' @param genes.pool, Gene list to base as the pool; NULL = all.
#' @param n.bin, number of bins to evaluate score across; default 25.
#' @param ctrl.size, control gene set size.
#' @param seed.use, random seed
#' @param assay, The seurat assay to use.  Defaults to the result of DefaultAssay()
#' @return A dataframe with CellBarcode and Score for this gene set
#' @export
#' @importFrom Hmisc cut2
#' @importFrom Matrix colMeans rowMeans
CalculateModuleScoreAvg <- function(
seuratObj,
genes.list = NULL,
genes.pool = NULL,
n.bin = 25,
ctrl.size = 100,
assay = NULL,
seed.use = 1
) {
	set.seed(seed = seed.use)
	genes.old <- genes.list

	if (is.null(x = genes.list)) {
		stop("Missing input gene list")
	}

	if (is.null(assay)) {
		assay <- Seurat::DefaultAssay(seuratObj)
	}

	genes.list <- intersect(x = genes.list, y = rownames(seuratObj))
	if (length(genes.list) == 0) {
		print(paste0('No matching genes found in the seurat object from the gene list, attempting to match case.'))
		genes.list <- Seurat::CaseMatch(genes.old, match = rownames(seuratObj))
	}

	if (length(genes.list) == 0) {
		print(paste0('No matching genes found in the seurat object from the gene list, aborting'))
		return(NA)
	}

	if (is.null(x = genes.pool)) {
		genes.pool = rownames(seuratObj)
	}

	data.avg <- Matrix::rowMeans(Seurat::GetAssayData(seuratObj, assay = assay)[genes.pool, ])
	data.avg <- data.avg[order(data.avg)]
	data.cut <- as.numeric(x = Hmisc::cut2(
	x = data.avg,
	m = round(x = length(x = data.avg) / n.bin)
	))
	names(x = data.cut) <- names(x = data.avg)

	genes.use <- genes.list
	ctrl.use <- character()
	for (j in 1:length(x = genes.use)) {
		ctrl.use <- c(ctrl.use, names(x = sample(
		x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
		size = ctrl.size,
		replace = FALSE
		)
		))
	}

	ctrl.use <- unique(ctrl.use)
	ctrl.scores <- matrix(
	data = numeric(length = 1L),
	nrow = length(x = ctrl.use),
	ncol = ncol(x = Seurat::GetAssayData(seuratObj, assay = assay))
	)

	data <- Seurat::GetAssayData(seuratObj, assay = assay)
	ctrl.scores <- Matrix::colMeans(x = data[ctrl.use, , drop = F])
	genes.scores <- Matrix::colMeans(x = data[genes.list, , drop = FALSE])

	return(data.frame(CellBarcode = colnames(x = Seurat::GetAssayData(seuratObj, assay = assay)), Score = (genes.scores - ctrl.scores)))
}