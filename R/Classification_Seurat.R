
#' @title A Title
#'
#' @description A description
#' @export
.ClassifyCellsCustom <- function(classifierFile,
	ClassifierNames = NULL,
	testing.data,
	log10T = T,
	assay = NULL
) {
	if (!file.exists(classifierFile)) {
		stop(paste0("File does not exist: ", classifierPath))
	}

	MultiClassifierResults <- readRDS(classifierFile)
	print(names(MultiClassifierResults))

	Available.Classifiers <- names(MultiClassifierResults)[grepl("classifier", names(MultiClassifierResults))]

	print(Available.Classifiers)

	if (is.null(ClassifierNames)) {
		ClassifierNames <- Available.Classifiers
	}

	if (length(which(ClassifierNames %in% Available.Classifiers)) < length(ClassifierNames)) {
		print("Names in ClassifierNames did not match whats available")
		ClassifierNames <- Available.Classifiers
	}

	MultiClassifier <- MultiClassifierResults[Available.Classifiers]

	if (class(testing.data)[1] == "seurat") {
		print("Converting Seurat Sparse Matrix to one with 0s .... ")

		if (is.null(assay)) {
			assay = DefaultAssay(testing.data)
		}

		testing.data <- Matrix::as.matrix(Matrix::t(Seurat::GetAssayData(testing.data, assay = assay, slot = 'data')))
	} else {
		if (! (class(testing.data)[1] %in% c("dgCMatrix", "matrix", "Matrix"))) {
			print("DGE mat not Seurat or dgCMatrix")
			print("Make sure columns are cells and rows are genes and convert to expected format for speed")
			testing.data <- Matrix::as.matrix((testing.data))
		}
	}

	if (log10T) {
		print("log-transforming")
		testing.data <- log10(testing.data + 1)
	}

	print("Starting Classification")

	testing.data.yhat <- data.frame(lapply(Available.Classifiers, function(ClassifX){
		print(ClassifX)

		tempClassifier <- MultiClassifier[[ClassifX]]

		tempClassifier$coefnames[! (tempClassifier$coefnames %in% colnames(testing.data))]

		y.hat <- stats::predict(tempClassifier, newdata = testing.data[,])

		return(y.hat)
	}))

	colnames(testing.data.yhat) <- Available.Classifiers

	NotLevelName <- levels(testing.data.yhat[, 1])[grep("Not", levels(testing.data.yhat[, 1]))]

	testing.data.yhat$CountNot <- apply(testing.data.yhat, 1, function(x) sum(grepl(NotLevelName, x)))
	testing.data.yhat$CountTot <- rep(length(Available.Classifiers), nrow(testing.data.yhat))
	testing.data.yhat$NotProb <- testing.data.yhat$CountNot / testing.data.yhat$CountTot

	rownames(testing.data.yhat) <- rownames(testing.data)

	return(list(
		yhat.DF = testing.data.yhat,
		classification.levels = levels(testing.data.yhat[, 1]),
		MultiClassifierResults = MultiClassifierResults
	))
}


#' @title ClassifySeurat
#'
#' @description This will classify cells in a seurat object using the specified classifiers
#' @return A modified Seurat object.
#' @export
ClassifySeurat <- function(seuratObj,
	trainedClassifiers = NULL,
	moduleScoreGeneLists = NULL,
	savePrefix = NULL,
	doGarnettClassify = F,
	garnett.path = "./data/Garnett/pbmc_classification.txt",
	assay = NULL
) {
	ret <- list()

	if (doGarnettClassify) {
		ret$Garnett <- list()
	}

	if (is.null(assay)) {
		assay = DefaultAssay(seuratObj)
	}

	if (!is.null(moduleScoreGeneLists)) {
		ret$SeuratGeneScore <- .DoModuleScoreGeneListClassify(seuratObj, moduleScoreGeneLists = moduleScoreGeneLists, csvFile = paste0(savePrefix, ".SGS.csv"))
	}

	if (!is.null(trainedClassifiers)) {
		ret$MCE.yhat <- list()
		ret$MCE <- list()
		for (name in names(trainedClassifiers)) {
			ret$MCE.yhat[name] <- .DoMCEClassify(
				seuratObj,
				assay = assay,
				classifierFile = trainedClassifiers[name],
				savePath = paste(savePrefix, ".mce.", name, ".rds")
			)

			ret$MCE[name] <- 1 - ret$MCE.yhat[name]$NotProb
		}
	}

	return(ret)
}


.DoMCEClassify <- function(seuratObj, assay, classifierFile, savePath = NULL) {
	testData <- log10(Matrix::as.matrix(Matrix::t(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'data')))+1)

	if (!is.null(savePath) && file.exists(savePath)) {
		print(paste0("Classification already done ... loading: ", savePath))
		return(readRDS(savePath))
	}

	ret <- .ClassifyCellsCustom(classifierFile = classifierFile, testing.data = testData, log10T = F)

	if (!is.null(savePath)) {
		saveRDS(ret, savePath)
	}

	return(ret)
}


ScoreScheme <- list(
 	CD8T = c(Lymph = '> 0.5', CD8T = '> 0.5', CD4T = '< 0.5', NK = '< 0.5', B = '< 0.5'),
	NK = c(Lymph = '> 0.5', CD8T = '< 0.5', CD4T = '< 0.5', NK = '> 0.5', B = '< 0.5'),
	LymphNotTBNK = c(Lymph = '> 0.5', CD8T = '< 0.5', CD4T = '< 0.5', NK = '< 0.5', B = '< 0.5'),
	NotLymphTBNK_MCE = c(Lymph = '< 0.5', CD8T = '< 0.5', CD4T = '< 0.5', NK = '< 0.5', B = '< 0.5')
)

.ScoreCellTypeFromMCE <- function(data, scoreScheme = ScoreScheme, savePath = NULL){
	ret <- data.frame()

	for (type in names(scoreScheme)) {
		rules <- scoreScheme[[names]]
		colName <- paste0(names,'.Call')
		ret[colName] <- logical()

		for (field in names(rules)) {
			for (i in 1:nrow(data)) {
				val <- data[i,field]
				if (!eval(parse_expr(paste0(field, rules[field])))) {
					data[i, colName] <- FALSE
					break
				}
			}
		}

		ret[i, colName] <- FALSE
	}

	if (!is.null(savePath)) {
		write.table(ret, file = savePath, sep = '\t', row.names = F, quote = F)
	}

	return(ret)
}


.DoModuleScoreGeneListClassify <- function(seuratObj, moduleScoreGeneLists, csvFile = NULL) {
	print("Starting Seurat's AddModule Scoring for GeneSets")
	SeuratGeneScore <- list()

	for (geneList in names(moduleScoreGeneLists)) {
		print(geneList)
		genes <- moduleScoreGeneLists[[geneList]]
		if (!is.list(genes)) {
			genes <- list(genes)
		}

		if (length(genes) == 0) {
			warning("length of gene list is 0")
			next
		}

		#TODO: is there a quantitative way to best estimate n.bin based on number of genes?
		logBinSize <- F
		for (binSize in c(25, 30 , 50)) {
			MS.temp <- try(
				AddModuleScoreAvg(SeurObj = seuratObj,
					genes.list = genes,
					genes.pool = NULL,
					n.bin = binSize,
					seed.use = 123,
					ctrl.size = 100,
					enrich.name = geneList,
					random.seed = 1, returnScore = T
				),
				silent = T
			)

			if (class(MS.temp) == "try-error") {
				print(paste0("Bin size ", binSize, " was problematic"))
				logBinSize <- T
			} else {
				if (logBinSize)
					print(paste0("Bin size ", binSize, " worked!"))

				SeuratGeneScore[[geneList]] <- MS.temp
				break
			}
		}

		if (!(geneList %in% names(SeuratGeneScore))) {
			print(paste0('AddModuleScoreAvg failed for gene set: ', geneList))
		}
	}

	#Seurate gene score (SGS)
	SGS.DF <- as.data.frame(SeuratGeneScore)
	print(colnames(SeuratGeneScore))
	SGS.DF <- SGS.DF[, ! grepl(".cID", colnames(SGS.DF)), drop = F]
	SGS.DF$CellBarcode <- rownames(SGS.DF)

	if (!is.null(csvFile)) {
		write.csv(SGS.DF, file = csvFile, row.names = FALSE)
	}

	return(SGS.DF)
}


