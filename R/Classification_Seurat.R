#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
SaveDimRedux_SERIII <- function(seuratObj, reductions=c("pca", "tsne", "umap"),
                                save.path=NA, maxPCAcomps=10, nameID=""){

  if(is.na(save.path)){
    stop("save.path is NA")
  }

  require(data.table)

  tempDT <- data.table(cbind(1:ncol(serObj), colnames(serObj)))
  rownames(tempDT) <- colnames(serObj)
  colnames(tempDT) <- c("nID", "cID")

  #tempDT <- merge(tempDT, )

  if("tsne" %in% reductions) {
    if(is.null(serObj@reductions$tsne)) print("tsne slot NULL") else {

    }
    tsneDT <- data.table(serObj@reductions$tsne@cell.embeddings)
    tsneDT$cID <- rownames(serObj@reductions$tsne@cell.embeddings)
    tempDT <- merge(tempDT, tsneDT, by="cID")
  }
  if("pca" %in% reductions) {
    if(is.null(serObj@reductions$pca)) print("pca slot NULL") else {
      pcaDT <- data.table(serObj@reductions$pca@cell.embeddings[,1:maxPCAcomps])
      pcaDT$cID <- rownames(serObj@reductions$pca@cell.embeddings)
      tempDT <- merge(tempDT, pcaDT, by="cID")
    }

  }
  if("umap" %in% reductions) {
    if(is.null(serObj@reductions$umap)) print("umap slot NULL") else {
      umapDT <- data.table(serObj@reductions$umap@cell.embeddings)
      umapDT$cID <- rownames(serObj@reductions$umap@cell.embeddings)
      tempDT <- merge(tempDT, umapDT, by="cID")
    }

  }

  print("saving DimRedux")
  save.path <- paste(save.path, "/", nameID, "_DimReduxComps.csv", sep="")
  write.csv(tempDT, file = save.path, row.names=TRUE, col.names = TRUE)


}






#' @title AddModuleScore_SERIII
#'
#' @description Using Seurat v2, AddModuleScore and converting to v3.
#'     this worked in v2 but their new version 3 breaks so this replaces it
#' @param object, A Seurat object.
#' @param genes.list, Gene list to obtain a score for
#' @param genes.pool, Gene list to base as the pool; NULL = all.
#' @param n.bin, number of bins to evaluate score across; default 25.
#' @param seed.use, seed use.
#' @param ctrl.size, control gene set size.
#' @param use.k, use k ... decap
#' @param enrich.name, A name for the assesment
#' @param random.seed, random seed
#' @return A modified Seurat object.
#' @keywords SerIII_AddModuleScore
#' @export
#' @examples
AddModuleScore_SERIII <- function(
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
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
ClassifyCellsCustom_SERIII <- function(Classifier.rds.path = "",
                                       ClassifierNames="",
                                       testing.data, log10T=T, returnTraining=F){
  #March-2019 version


  if(!file.exists(Classifier.rds.path)){
    stop("Check classifier path ....")
  }

  MultiClassifierResults <- readRDS(Classifier.rds.path)
  Available.Classifiers  <- names(MultiClassifierResults)[grepl("classifier", names(MultiClassifierResults))]

  print(Available.Classifiers)

  if(ClassifierNames[1]=="") ClassifierNames <- Available.Classifiers

  if(length(which(ClassifierNames %in% Available.Classifiers)) < length(ClassifierNames)) {
    print("Names in ClassifierNames did not match whats available")
    ClassifierNames <- Available.Classifiers
  }

  MultiClassifier        <- MultiClassifierResults[Available.Classifiers]


  if(class(testing.data)[1]=="seurat") {
    print("Converting Seurat Sparse Matrix to one with 0s .... ")

    testing.data <- Matrix::as.matrix(Matrix::t(testing.data@assays$RNA@data))

  } else {
    if(!(class(testing.data)[1] %in% c("dgCMatrix", "matrix", "Matrix")) ) {
      print("DGE mat not Seurat or dgCMatrix")
      print("Make sure columns are cells and rows are genes and convert to expected format for speed")
      testing.data <- Matrix::as.matrix((testing.data))

    }
  }


  if(log10T) {
    print("log-transforming")
    testing.data <- log10(testing.data+1)
  }


  print("Starting Classification")

  testing.data.yhat <- data.frame(lapply(Available.Classifiers, function(ClassifX){

    print(ClassifX)

    tempClassifier <- MultiClassifier[[ClassifX]]

    tempClassifier$coefnames[!(tempClassifier$coefnames %in% colnames(testing.data))]


    y.hat <- predict(tempClassifier, newdata = testing.data[,])




    return(y.hat)
  }))

  colnames(testing.data.yhat) <- Available.Classifiers

  NotLevelName <- levels(testing.data.yhat[,1])[grep("Not", levels(testing.data.yhat[,1]))]

  testing.data.yhat$CountNot <- apply(testing.data.yhat, 1, function(x) sum(grepl(NotLevelName, x)))
  testing.data.yhat$CountTot <- rep(length(Available.Classifiers), nrow(testing.data.yhat))

  testing.data.yhat$NotProb <- testing.data.yhat$CountNot/testing.data.yhat$CountTot



  rownames(testing.data.yhat) <- rownames(testing.data)


  return(list(yhat.DF = testing.data.yhat,
              classification.levels = levels(testing.data.yhat[,1]),
              Available.Classifiers = Available.Classifiers,
              log10T = log10T,
              DGEcolNames = colnames(testing.data),
              DGErowNames = rownames(testing.data),
              MultiClassifierResults = ifelse(returnTraining, MultiClassifierResults, "NULL")))


}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
predict_SERIII <- function(ProcSERobj.path = NULL, PatternOfProcSERobj="_proc.rds",
                           classification.path = NULL, file.select = NULL,
                           TrainedClassifiers.path = "../PBMC3k/data",
                           save.fig.path = NULL, col_vector=NULL, returnLS = F, GarnettClassify=F,
                           Garnett.path = "./data/Garnett/pbmc_classification.txt", MCEClassify=T,
                           ModuleScoreGeneListClassify=F, ModuleScoreGeneLists=NULL,
                           RhesusConvDavid.path = "./data/Rhesus/David6.8_ConvertedRhesus_ENSMMUG.txt",
                           RhesusConvDavid = F, ENSMB.tag="ENSMM", yhat.save.path = NA, cleanName=T){


  if(is.null(ProcSERobj.path)){
    stop("path does not exists")
  }

  ClassifiersLS <- list()
  ClassifiersLS$MCEyhat <- list()
  ClassifiersLS$MCEyhat$CD8T <- list()
  ClassifiersLS$MCEyhat$CD4T <- list()
  ClassifiersLS$MCEyhat$NK <- list()
  ClassifiersLS$MCEyhat$B <- list()
  ClassifiersLS$MCEyhat$Lymph <- list()

  if(GarnettClassify) ClassifiersLS$Garnett <- list()

  if(ModuleScoreGeneListClassify) ClassifiersLS$SeuratGeneScore <- list()

  SERObjects_processed.paths <- list.files(ProcSERobj.path, full.names = T, pattern = PatternOfProcSERobj)

  if(length(SERObjects_processed.paths)==0){
    print("no files found...")
  }else {

    if(is.null(classification.path)) classification.path <- ProcSERobj.path
    if(is.null(save.fig.path)) save.fig.path <- ProcSERobj.path
    if(is.null(col_vector)) col_vector <- colors(distinct = T)

    if(!is.null(file.select)) SERObjects_processed.paths <- SERObjects_processed.paths[grepl(file.select, SERObjects_processed.paths)]

    for(SERObj.path in SERObjects_processed.paths){

      print(basename(SERObj.path))
      tempSER <- readRDS(SERObj.path)
      #tempSER <- SeuratObjs



      tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))


      # ModuleScoreGeneLists <- CTL_Immune_GeneList(QuickGO.path="/Volumes/Maggie/Work/OHSU/Eisa/R/scRNASeq/data/QuickGO")





      if(ModuleScoreGeneListClassify){

        print("Starting Seurat's AddModule Scoring for GeneSets")
        #tempLSScores <- list()
        ClassifiersLS$SeuratGeneScore[[tempName]] <- list()


        for(GeneList in names(ModuleScoreGeneLists)){
          # GeneList = names(ModuleScoreGeneLists)[1]
          print(GeneList)



          if(length(ModuleScoreGeneLists[[GeneList]])>0){

            MS.temp <-  try(AddModuleScore_SERIII(SeurObj=tempSER,
                                                  genes.list = list(ModuleScoreGeneLists[[GeneList]]),
                                                  genes.pool = NULL,
                                                  n.bin = 25,
                                                  seed.use = 123,
                                                  ctrl.size = 100,
                                                  enrich.name = GeneList,
                                                  random.seed = 1, returnScore = T), silent = T)

            if(class(MS.temp)=="try-error") {

              #TODO: is there a quantitative way to best estimate n.bin based on number of genes?

              print("bin size 25 was problematic... trying 30")
              MS.temp <-  try(AddModuleScore_SERIII(SeurObj=tempSER,
                                                    genes.list = list(ModuleScoreGeneLists[[GeneList]]),
                                                    genes.pool = NULL,
                                                    n.bin = 30,
                                                    seed.use = 123,
                                                    ctrl.size = 100,
                                                    enrich.name = GeneList,
                                                    random.seed = 1, returnScore = T), silent = T)

              if(class(MS.temp)=="try-error") {

                #TODO: is there a quantitative way to best estimate n.bin based on number of genes?

                print("bin size 30 was problematic... trying 50")
                MS.temp <-  try(AddModuleScore_SERIII(SeurObj =tempSER,
                                                      genes.list = list(ModuleScoreGeneLists[[GeneList]]),
                                                      genes.pool = NULL,
                                                      n.bin = 50,
                                                      seed.use = 123,
                                                      ctrl.size = 100,
                                                      enrich.name = GeneList,
                                                      random.seed = 1, returnScore = T), silent = T)

                if(class(MS.temp)=="try-error") {

                  warning("bin size 25, 30 and 50 failed")

                } else {

                  print("bin size 50 worked!")
                  ClassifiersLS$SeuratGeneScore[[tempName]][[GeneList]] <- MS.temp

                }

              } else {
                print("bin size 30 worked!")
                ClassifiersLS$SeuratGeneScore[[tempName]][[GeneList]] <- MS.temp
              }

            } else {
              #print("bin size 25 worked!")
              ClassifiersLS$SeuratGeneScore[[tempName]][[GeneList]] <- MS.temp
            }

          } else {
            warning("length of gene list is <= 0 ")

          }
        }



        #Seurate gene score (SGS)
        SGS.DF <- as.data.frame(ClassifiersLS$SeuratGeneScore[[tempName]])
        SGS.DF <- SGS.DF[,!grepl(".cID", colnames(SGS.DF))]
        SGS.DF$cID <- rownames(SGS.DF)

        ClassifiersLS$SeuratGeneScore[[tempName]]$SGS.DF <- SGS.DF#ClassifiersLS$SeuratGeneScore

        if(!is.na(yhat.save.path)){

          yhat.s.path <- paste(yhat.save.path, "/", basename(SERObj.path), "_SGS.csv", sep="")
          write.csv(SGS.DF, file = yhat.s.path, row.names=TRUE, col.names = T)



        }
      }



      #CD8 T cells

      if(MCEClassify){

        if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))){

          #one can directly give the Seurat object to the ClassifyCellsCustom_SERIII()
          #since looping, its faster to compute the non-sparse log once

          X.SerObj.temp <- log10(Matrix::as.matrix(Matrix::t(tempSER@assays$RNA@data))+1)

          if(RhesusConvDavid){

            if(!file.exists(RhesusConvDavid.path)) print("David file not found") else {
              colnames(X.SerObj.temp)[grepl(ENSMB.tag, colnames(X.SerObj.temp))] <- RhesusGeneDavidConv(ColNames2Conv=colnames(X.SerObj.temp),
                                                                                                        RhesusConvDavid.path=RhesusConvDavid.path)
            }


          }

          dim(X.SerObj.temp)


          #for now because I did that with the training set, when new training is done,
          #I will not do this, as its just fine to keep the dash and avoid dupes
          if(cleanName) colnames(X.SerObj.temp) <- gsub("-", "", colnames(X.SerObj.temp))


          ClassifiersLS$MCEyhat$CD8T[[tempName]]  <- list(MCE=ClassifyCellsCustom_SERIII(
            Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_CD8T.rds", sep=""),
            testing.data = X.SerObj.temp, log10T=F))

          ClassifiersLS$MCEyhat$CD8T[[tempName]]$MCE$log10T = T

          #either save results as tsv or csv
          #ClassifiersLS$MCEyhat$CD8T[[tempName]]$yhat.DF or entire object
          #but its best so compute all inferences, save as rds, but as a csv, save a complete sheet of all classifiers
          saveRDS(ClassifiersLS$MCEyhat$CD8T[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))


        } else {
          print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))
          print("Already done ... loading for LS")
          ClassifiersLS$MCEyhat$CD8T[[tempName]] <- readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))
        }




        #CD4T cells

        if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))){

          if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

          ClassifiersLS$MCEyhat$CD4T[[tempName]]  <- list(MCE=ClassifyCellsCustom_SERIII(
            Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_CD4T.rds", sep=""),
            testing.data = X.SerObj.temp, log10T=F))
          ClassifiersLS$MCEyhat$CD4T[[tempName]]$MCE$log10T = T


          saveRDS(ClassifiersLS$MCEyhat$CD4T[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))


        } else {
          print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))
          print("Already done ... loading for LS")
          ClassifiersLS$MCEyhat$CD4T[[tempName]] <- readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))
        }


        #NK cells

        if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))){

          if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

          ClassifiersLS$MCEyhat$NK[[tempName]]  <- list(MCE=ClassifyCellsCustom_SERIII(
            Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_NK.rds", sep=""),
            testing.data = X.SerObj.temp, log10T=F))
          ClassifiersLS$MCEyhat$NK[[tempName]]$MCE$log10T = T

          #either save results as tsv or csv

          saveRDS(ClassifiersLS$MCEyhat$NK[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))


        } else {
          print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))
          print("Already done ... loading for LS")
          ClassifiersLS$MCEyhat$NK[[tempName]] <- readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))
        }



        #B cells

        if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))){

          if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

          ClassifiersLS$MCEyhat$B[[tempName]]  <- list(MCE=ClassifyCellsCustom_SERIII(
            Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_B.rds", sep=""),
            testing.data = X.SerObj.temp, log10T=F))
          ClassifiersLS$MCEyhat$B[[tempName]]$MCE$log10T = T


          saveRDS(ClassifiersLS$MCEyhat$B[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))


        } else {
          print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))
          print("Already done ... loading for LS")
          ClassifiersLS$MCEyhat$B[[tempName]] <- readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))
        }



        #Lymph cells

        if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))){

          if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

          ClassifiersLS$MCEyhat$Lymph[[tempName]]  <- list(MCE=ClassifyCellsCustom_SERIII(
            Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_Lymph.rds", sep=""),
            testing.data = X.SerObj.temp, log10T=F))

          ClassifiersLS$MCEyhat$Lymph[[tempName]]$MCE$log10T = T


          saveRDS(ClassifiersLS$MCEyhat$Lymph[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))
        } else {
          print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))
          print("Already done ... loading for LS")
          ClassifiersLS$MCEyhat$Lymph[[tempName]] <- readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))
        }

        print("All classifiers are done with this file ... ")



        yhat.Combo <- as.data.frame(cbind(1-ClassifiersLS$MCEyhat$CD8T[[tempName]]$MCE$yhat.DF$NotProb,
                                          1-ClassifiersLS$MCEyhat$CD4T[[tempName]]$MCE$yhat.DF$NotProb,
                                          1-ClassifiersLS$MCEyhat$NK[[tempName]]$MCE$yhat.DF$NotProb,
                                          1-ClassifiersLS$MCEyhat$B[[tempName]]$MCE$yhat.DF$NotProb,
                                          1-ClassifiersLS$MCEyhat$Lymph[[tempName]]$MCE$yhat.DF$NotProb))

        print(dim(yhat.Combo))
        colnames(yhat.Combo) <- c("CD8T", "CD4T", "NK", "B", "Lymph")
        rownames(yhat.Combo) <- colnames(tempSER@assays$RNA@data)

        Classificatio.meta.data <- list()

        #TODO?: The prob association is set to 50% changing it will affect classification stats e.g. precision/recall etc
        Classificatio.meta.data$CD8T_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] > 0.5 &
                   xR["CD8T"] > 0.5 &
                   xR["CD4T"] < 0.5 &
                   xR["NK"] < 0.5 &
                   xR["B"] < 0.5, 1, 0)
        })

        Classificatio.meta.data$CD4T_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] > 0.5 &
                   xR["CD8T"] < 0.5 &
                   xR["CD4T"] > 0.5 &
                   xR["NK"] < 0.5 &
                   xR["B"] < 0.5, 1, 0)
        })


        Classificatio.meta.data$B_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] > 0.5 &
                   xR["CD8T"] < 0.5 &
                   xR["CD4T"] < 0.5 &
                   xR["NK"] < 0.5 &
                   xR["B"] > 0.5, 1, 0)
        })


        Classificatio.meta.data$NK_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] > 0.5 &
                   xR["CD8T"] < 0.5 &
                   xR["CD4T"] < 0.5 &
                   xR["NK"] > 0.5 &
                   xR["B"] < 0.5, 1, 0)
        })


        Classificatio.meta.data$LymphNotTBNK_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] > 0.5 &
                   xR["CD8T"] < 0.5 &
                   xR["CD4T"] < 0.5 &
                   xR["NK"] < 0.5 &
                   xR["B"] < 0.5, 1, 0)
        })

        Classificatio.meta.data$NotLymphTBNK_MCE <- apply(yhat.Combo, 1, function(xR){
          ifelse(xR["Lymph"] < 0.5 &
                   xR["CD8T"] < 0.5 &
                   xR["CD4T"] < 0.5 &
                   xR["NK"] < 0.5 &
                   xR["B"] < 0.5, 1, 0)
        })

        Classificatio.meta.data <- as.data.frame(Classificatio.meta.data)



        if(!is.na(yhat.save.path)){

          yhat.s.path <- paste(yhat.save.path, "/", basename(SERObj.path), "_MCEmeta.csv", sep="")
          if(file.exists(yhat.s.path)){
            print(yhat.s.path); print("already exists")
          } else {
            write.csv(Classificatio.meta.data, file = yhat.s.path, row.names=TRUE, col.names = T)

          }

          yhat.s.path <- paste(yhat.save.path, "/", basename(SERObj.path), "_MCEyhat.csv", sep="")

          if(file.exists(yhat.s.path)){
            print(yhat.s.path); print("already exists")
          } else {
            write.csv(yhat.Combo, file = yhat.s.path, row.names=TRUE, col.names = T)

          }

        }


        ClassifiersLS$classificationLS[[tempName]] <- list(Classificatio.meta.data = Classificatio.meta.data,
                                                           yhat.Combo=yhat.Combo)


      } # END MCE_classification



    }



  }
  if(returnLS) return(ClassifiersLS)



}




