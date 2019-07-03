#' @import Seurat
#' @import glmnet
#' @import rpart
#' @import rpart.plot
#' @import ranger
#' @importFrom stats predict

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
GetTraingStats <- function(ClassifierSet=NA, AvailableTrainConfMatrix="", AvailableTestConfMatrix=""){

  tempDF1 <- data.frame(lapply(AvailableTrainConfMatrix, function(conftabN){
    c(ClassifierSet[[conftabN]]$byClass[c(1:4, 6:10)],
      ClassifierSet[[conftabN]]$overall[c(1,3,4)])
  }))
  colnames(tempDF1) <- AvailableTrainConfMatrix

  tempDF2 <- data.frame(lapply(AvailableTestConfMatrix, function(conftabN){
    c(ClassifierSet[[conftabN]]$byClass[c(1:4, 6:10)],
      ClassifierSet[[conftabN]]$overall[c(1,3,4)])
  }))
  colnames(tempDF2) <- AvailableTestConfMatrix

  data.frame(cbind(tempDF1, tempDF2))[,sort(c(AvailableTestConfMatrix, AvailableTrainConfMatrix))]
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
GetEnsmbleVarImportance <- function(ClassiferSet, AvailableClassifiers="", ScaleVarImp=T){

  if(AvailableClassifiers[1]=="") AvailableClassifiers =  names(ClassifierSet)

  FoundVarImps <- lapply(AvailableClassifiers, function(clN){
    print(clN)

    tempOut <-   try(varImp(ClassiferSet[[clN]],scale = ScaleVarImp), silent = T)
    if(class(tempOut)[1]=="try-error") {
      return(NULL)
      print("No VarImp...")
    } else{
      tempOut2 <- tempOut$importance
      if(ncol(tempOut2)<2){
        tempOut2 <- data.frame(cbind(tempOut2, rep(NA, length(tempOut2))))
      }
      return(data.frame(t(tempOut2)))
    }

  })
  names(FoundVarImps) <- AvailableClassifiers
  FoundVarImps[sapply(FoundVarImps, is.null)] <- NULL


  FoundVarImpsDF <- as.data.frame(data.table::rbindlist(FoundVarImps))
  FoundVarImpsDF <- FoundVarImpsDF[is.odd(1:nrow(FoundVarImpsDF)),]
  rownames(FoundVarImpsDF) <- names(FoundVarImps)
  FoundVarImpsDF
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
ClassifyCellsCustom <- function(Classifier.rds.path = "", ClassifierNames="", testing.data, log10T=T, returnTraining=F){

  if(file.exists(Classifier.rds.path)){
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

      testing.data <- Matrix::as.matrix(t(testing.data@data))
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

      y.hat <- stats::predict(tempClassifier, newdata = testing.data)


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

  } else {
    print("Check classifier path ....")
  }
}



#' @title MultiClassifier_Cells
#'
#' @description This function takes in Seurat objects and performs machine learning; generalized for common use.
#'     However as with any approach, especially in machine learning and statistical inference world, robust assessments are needed. 
#'     This is a first step.
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords Machine Learning, Seurat, classification
#' @importFrom fastAdaboost adaboost
#' @importFrom stats predict
#' @export
MultiClassifier_Cells <- function (object_train,
                                  object_test,
                                  data2use = "data",
                                  assay = "RNA",
                                  training.genes = NULL,
                                  Y_Train_True = NULL,
                                  Y_Test_True = NULL,
                                  log10p1=F,
                                  asinhp1 = F,
                                  do_Adaboost = F,
                                  do_Caret_LinSVM = F,
                                  do_ranger_RF = F,
                                  do_Caret_NaiveBayes = F,
                                  do_Caret_NNet = F,
                                  do_Caret_radialSVM = F,
                                  do_ElasticNet = F,
                                  do_tree = F,
                                  do_tree_gini = F,
                                  do_Caret_pcaNNet = F,
                                  do_Caret_rangerRF = F,
                                  do_glmStepAIC=F,
                                  do_stepLDA=F,
                                  crossValReps = 3,
                                  NcrossVal = 10, 
                                  preProc = c( "scale", "center", "nzv", "zv"))
{


  results_ls <- list()
  
  training.classes <- as.vector(x = Y_Train_True)
  #training.genes   <- #SetIfNull(x = training.genes, default = rownames(x = object_train@assays$RNA@data))
  training.data    <- as.data.frame(Matrix::t(Seurat::GetAssayData(object = object_train, assay = assay, slot = data2use)))[,training.genes]
  
  
  testing.classes <- as.vector(x = Y_Test_True)
  testing.genes   <- training.genes #SetIfNull(x = training.genes, default = rownames(x = object_test@assays$RNA@data))
  testing.data    <- as.data.frame(Matrix::t(Seurat::GetAssayData(object = object_test, assay = assay, slot = data2use)))[,testing.genes ]
  



  if(log10p1) training.data <- log10(training.data + 1)
  if(log10p1) testing.data  <- log10(testing.data + 1)
  
  if(asinhp1) training.data <- asinh(training.data)
  if(asinhp1) testing.data  <- asinh(testing.data)
  




  training.data$Class <- factor(x = training.classes)



  #colnames(training.data) <- gsub("-", "", colnames(training.data))
  #colnames(testing.data) <- gsub("-", "", colnames(testing.data))

  results_ls$test_data           <- testing.data
  results_ls$test_Y              <- testing.classes
  results_ls$train_data          <- training.data
  results_ls$train_Y             <- training.classes



  #############template
  # if(do_newModel){
  #
  #
  #   print("running newModel")
  #   TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
  #
  #
  #  newModel_classifier <- train(Class ~ ., data = training.data,
  #                                method = "newModel",
  #                                trControl= TrainingParameters,
  #                                preProcess = c("scale","center"),
  #                                na.action = na.omit)
  #
  #
  #   newModel_yhat <- stats::predict(newModel_classifier, testing.data)
  #
  #
  #   results_ls$Caret_newModel_train_yhat      <- stats::predict(newModel_classifier, newdata=training.data[,which(colnames(training.data) !="Class")])
  #
  #
  #   results_ls$Caret_newModel_classifier      <- newModel_classifier
  #   results_ls$Caret_newModel_test_yhat       <- newModel_yhat
  #   results_ls$Caret_newModel_test_ConfMat    <- confusionMatrix(newModel_yhat, (Y_Test_True))
  #   results_ls$Caret_newModel_test_ConfTab    <- table(pred=newModel_yhat, truth=Y_Test_True)
  #
  #   results_ls$Caret_newModel_train_ConfMat   <- confusionMatrix((results_ls$Caret_newModel_train_yhat), (training.data$Class))
  #   results_ls$Caret_newModel_train_ConfTab   <- table(pred=results_ls$Caret_newModel_train_yhat, truth=training.data$Class)
  #
  #
  #   colnames(training.data) <- colnames(results_ls$train_data)
  #   colnames(testing.data) <- colnames(results_ls$test_data)
  #
  #
  # }

  #### Caret based classifiers
  
  TrainingParameters <- trainControl(method = "repeatedcv", 
                                     number = NcrossVal, 
                                     repeats=crossValReps, 
                                     ## Estimate class probabilities
                                     classProbs = TRUE,
                                     ## Evaluate performance using 
                                     ## the following function
                                     summaryFunction = twoClassSummary)
  
  
  print("getting pre-processing values")
  preProcValues <- preProcess(training.data, method = preProc);preProcValues
  
  print("Transforming training and testing data with preprocessing values")
  
  trainTransformed <- stats::predict(preProcValues, training.data)
  testTransformed <- stats::predict(preProcValues, testing.data)
  
  #heatmap(as.matrix(testTransformed))
  results_ls$TrainingParameters        <- TrainingParameters
  results_ls$preProcValues             <- preProcValues
  
  
  
  
  if(do_Caret_LinSVM) {
    # training model with SVM
    print("running Caret LinSVM")
    
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    SVM_linear_classifier <- train(Class ~ ., data = trainTransformed,
                                   method = "svmLinear",
                                   trControl= TrainingParameters,
                                   #preProcess = preProc,
                                   na.action = na.omit, 
                                   metric = "ROC")
    
    SVM_Lin_yhat <- stats::predict(SVM_linear_classifier, testTransformed)
    results_ls$Caret_LinSVM_train_yhat      <- stats::predict(SVM_linear_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    #stats::predict(SVM_linear_classifier, newdata=trainTransformed[,training.genes])
    
    
    results_ls$Caret_LinSVM_classifier      <- SVM_linear_classifier
    results_ls$Caret_LinSVM_test_yhat       <- SVM_Lin_yhat
    results_ls$Caret_LinSVM_test_ConfMat    <- confusionMatrix(SVM_Lin_yhat, (Y_Test_True))
    results_ls$Caret_LinSVM_test_ConfTab    <- table(pred=SVM_Lin_yhat, truth=Y_Test_True)
    
    results_ls$Caret_LinSVM_train_ConfMat   <- confusionMatrix((results_ls$Caret_LinSVM_train_yhat), (training.data$Class))
    results_ls$Caret_LinSVM_train_ConfTab   <- table(pred=results_ls$Caret_LinSVM_train_yhat, truth=training.data$Class)
    
  }
  
  if(do_Caret_radialSVM) {
    # training model with SVM
    print("running Caret radialSVM")
    
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    SVM_radial_classifier <- train(Class ~ ., data = trainTransformed,
                                   method = "svmRadial",
                                   trControl= TrainingParameters,
                                   #preProcess = preProc,
                                   na.action = na.omit, 
                                   metric = "ROC")
    
    SVM_Radial_yhat <- stats::predict(SVM_radial_classifier, testTransformed)
    results_ls$Caret_radialSVM_train_yhat      <- stats::predict(SVM_radial_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_radialSVM_classifier      <- SVM_radial_classifier
    results_ls$Caret_radialSVM_test_yhat       <- SVM_Radial_yhat
    results_ls$Caret_radialSVM_test_ConfMat    <- confusionMatrix(SVM_Radial_yhat, (Y_Test_True))
    results_ls$Caret_radialSVM_test_ConfTab    <- table(pred=SVM_Radial_yhat, truth=Y_Test_True)
    
    results_ls$Caret_radialSVM_train_ConfMat   <- confusionMatrix((results_ls$Caret_radialSVM_train_yhat), (trainTransformed$Class))
    results_ls$Caret_radialSVM_train_ConfTab   <- table(pred=results_ls$Caret_radialSVM_train_yhat, truth=trainTransformed$Class)
    
  }
  
  if(do_stepLDA){
    
    
    print("running stepLDA")
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps, 
    #                                    ## Estimate class probabilities
    #                                    classProbs = TRUE,
    #                                    ## Evaluate performance using 
    #                                    ## the following function
    #                                    summaryFunction = twoClassSummary)
    
    
    stepLDA_classifier <- train(Class ~ ., data = trainTransformed,
                                method = "stepLDA",
                                trControl= TrainingParameters,
                                #preProcess = preProc,
                                na.action = na.omit)
    
    
    stepLDA_yhat <- stats::predict(stepLDA_classifier, testTransformed)
    
    
    results_ls$Caret_stepLDA_train_yhat      <- stats::predict(stepLDA_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_stepLDA_classifier      <- stepLDA_classifier
    results_ls$Caret_stepLDA_test_yhat       <- stepLDA_yhat
    results_ls$Caret_stepLDA_test_ConfMat    <- confusionMatrix(stepLDA_yhat, (Y_Test_True))
    results_ls$Caret_stepLDA_test_ConfTab    <- table(pred=stepLDA_yhat, truth=Y_Test_True)
    
    results_ls$Caret_stepLDA_train_ConfMat   <- confusionMatrix((results_ls$Caret_stepLDA_train_yhat), (trainTransformed$Class))
    results_ls$Caret_stepLDA_train_ConfTab   <- table(pred=results_ls$Caret_stepLDA_train_yhat, truth=trainTransformed$Class)
    
    
    # colnames(training.data) <- colnames(results_ls$train_data)
    # colnames(testing.data) <- colnames(results_ls$test_data)
    
    
  }
  
  if(do_glmStepAIC){
    
    
    print("running glmStepAIC")
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    
    glmStepAIC_classifier <- train(Class ~ ., data = trainTransformed,
                                   method = "glmStepAIC",
                                   trControl= TrainingParameters,
                                   #preProcess = preProc,
                                   na.action = na.omit)
    
    
    glmStepAIC_yhat <- stats::predict(glmStepAIC_classifier, testTransformed)
    
    
    results_ls$Caret_glmStepAIC_train_yhat      <- stats::predict(glmStepAIC_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_glmStepAIC_classifier      <- glmStepAIC_classifier
    results_ls$Caret_glmStepAIC_test_yhat       <- glmStepAIC_yhat
    results_ls$Caret_glmStepAIC_test_ConfMat    <- confusionMatrix(glmStepAIC_yhat, (Y_Test_True))
    results_ls$Caret_glmStepAIC_test_ConfTab    <- table(pred=glmStepAIC_yhat, truth=Y_Test_True)
    
    results_ls$Caret_glmStepAIC_train_ConfMat   <- confusionMatrix((results_ls$Caret_glmStepAIC_train_yhat), (trainTransformed$Class))
    results_ls$Caret_glmStepAIC_train_ConfTab   <- table(pred=results_ls$Caret_glmStepAIC_train_yhat, truth=trainTransformed$Class)
    
    
    # colnames(training.data) <- colnames(results_ls$train_data)
    # colnames(testTransformed) <- colnames(results_ls$test_data)
    
    
  }
  
  
  
  if(do_tree){
    
    
    print("running tree based class with rpart")
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    
    
    rpart_classifier <- train(x=trainTransformed[,setdiff(colnames(trainTransformed), c("Class"))],
                              y=trainTransformed$Class,
                              method = "rpart",
                              trControl= TrainingParameters,
                              #preProcess =preProc,
                              na.action = na.omit,
                              parms = list(split = "information"),
                              tuneLength = 10, metric = "ROC")
    
    prp(rpart_classifier$finalModel, box.palette = "Blues", tweak = 1.2)
    results_ls$Caret_rpart_treeViz <- recordPlot()
    
    
    rpart_yhat <- stats::predict(rpart_classifier, testTransformed)
    
    
    results_ls$Caret_rpart_train_yhat      <- stats::predict(rpart_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_rpart_classifier      <- rpart_classifier
    results_ls$Caret_rpart_test_yhat       <- rpart_yhat
    results_ls$Caret_rpart_test_ConfMat    <- confusionMatrix(rpart_yhat, (Y_Test_True))
    results_ls$Caret_rpart_test_ConfTab    <- table(pred=rpart_yhat, truth=Y_Test_True)
    
    results_ls$Caret_rpart_train_ConfMat   <- confusionMatrix((results_ls$Caret_rpart_train_yhat), (trainTransformed$Class))
    results_ls$Caret_rpart_train_ConfTab   <- table(pred=results_ls$Caret_rpart_train_yhat, truth=trainTransformed$Class)
    
    
    #colnames(training.data) <- colnames(results_ls$train_data)
    #colnames(testing.data) <- colnames(results_ls$test_data)
    
    
  }
  
  if(do_tree_gini){
    
    
    print("running tree based class with rpart w/ gini index")
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    
    
    
    rpartGini_classifier <- train(x=trainTransformed[,setdiff(colnames(trainTransformed), c("Class"))],
                                  y=trainTransformed$Class,
                                  method = "rpart",
                                  trControl= TrainingParameters,
                                  #preProcess = preProc,
                                  na.action = na.omit,
                                  parms = list(split = "gini"),
                                  tuneLength = 10,
                                  metric = "ROC")
    
    prp(rpartGini_classifier$finalModel, box.palette = "Blues", tweak = 1.2)
    results_ls$Caret_rpartGini_treeViz <- recordPlot()
    
    
    rpartGini_yhat <- stats::predict(rpartGini_classifier, testTransformed)
    
    
    results_ls$Caret_rpartGini_train_yhat      <- stats::predict(rpartGini_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_rpartGini_classifier      <- rpartGini_classifier
    results_ls$Caret_rpartGini_test_yhat       <- rpartGini_yhat
    results_ls$Caret_rpartGini_test_ConfMat    <- confusionMatrix(rpartGini_yhat, (Y_Test_True))
    results_ls$Caret_rpartGini_test_ConfTab    <- table(pred=rpartGini_yhat, truth=Y_Test_True)
    
    results_ls$Caret_rpartGini_train_ConfMat   <- confusionMatrix((results_ls$Caret_rpartGini_train_yhat), (trainTransformed$Class))
    results_ls$Caret_rpartGini_train_ConfTab   <- table(pred=results_ls$Caret_rpartGini_train_yhat, truth=trainTransformed$Class)
    
    #results_ls
    
    # colnames(training.data) <- colnames(results_ls$train_data)
    # 
    # colnames(testing.data) <- colnames(results_ls$test_data)
    # 
    
  }
  
  if(do_Caret_NaiveBayes) {
    print("running Naive Bayes")
    
    # training model with SVM
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    #Naive Bayes algorithm
    NaiveBayes_classifier <- train(Class ~ ., data = trainTransformed,
                                   method = "nb",
                                   #preProcess=preProc,
                                   trControl= TrainingParameters,
                                   na.action = na.omit, 
                                   silent = T)
    
    NaiveBayes_yhat <- stats::predict(NaiveBayes_classifier, testTransformed)
    
    results_ls$NaiveBayes_train_yhat      <- stats::predict(NaiveBayes_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    results_ls$NaiveBayes_classifier      <- NaiveBayes_classifier
    results_ls$NaiveBayes_test_yhat       <- NaiveBayes_yhat
    results_ls$NaiveBayes_test_ConfMat    <- confusionMatrix(NaiveBayes_yhat, (Y_Test_True))
    results_ls$NaiveBayes_test_ConfTab    <- table(pred=NaiveBayes_yhat, truth=Y_Test_True)
    
    results_ls$NaiveBayes_train_ConfMat   <- confusionMatrix((results_ls$NaiveBayes_train_yhat), (trainTransformed$Class))
    results_ls$NaiveBayes_train_ConfTab   <- table(pred=results_ls$NaiveBayes_train_yhat, truth=trainTransformed$Class)
    
  }
  
  if(do_Caret_NNet) {
    print("running NNet")
    
    # training model with SVM
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    #Naive Bayes algorithm
    NNet_classifier <- train(Class ~ ., data = trainTransformed,
                             method = "nnet",
                             # preProcess=preProc,
                             trControl= TrainingParameters,
                             na.action = na.omit, 
                             silent = T, metric = "ROC")
    
    NNet_yhat <- stats::predict(NNet_classifier, testTransformed)
    results_ls$NNet_train_yhat      <- stats::predict(NNet_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$NNet_classifier      <- NNet_classifier
    results_ls$NNet_test_yhat       <- NNet_yhat
    results_ls$NNet_test_ConfMat    <- confusionMatrix(NNet_yhat, (Y_Test_True))
    results_ls$NNet_test_ConfTab    <- table(pred=NNet_yhat, truth=Y_Test_True)
    
    results_ls$NNet_train_ConfMat   <- confusionMatrix((results_ls$NNet_train_yhat), (trainTransformed$Class))
    results_ls$NNet_train_ConfTab   <- table(pred=results_ls$NNet_train_yhat, truth=trainTransformed$Class)
    
  }
  
  if(do_Caret_pcaNNet) {
    # training model with SVM
    print("running Caret pcaNNet ")
    
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    pcaNNetCaret_classifier <- train(Class ~ ., data = trainTransformed,
                                     method = "pcaNNet",
                                     trControl= TrainingParameters,
                                     #preProcess = preProc,
                                     na.action = na.omit, metric = "ROC", silent = T)
    
    pcaNNetCaret_yhat <- stats::predict(pcaNNetCaret_classifier, testTransformed)
    results_ls$Caret_pcaNNet_train_yhat      <- stats::predict(pcaNNetCaret_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_pcaNNet_classifier      <- pcaNNetCaret_classifier
    results_ls$Caret_pcaNNet_test_yhat       <- pcaNNetCaret_yhat
    results_ls$Caret_pcaNNet_test_ConfMat    <- confusionMatrix(pcaNNetCaret_yhat, (Y_Test_True))
    results_ls$Caret_pcaNNet_test_ConfTab    <- table(pred=pcaNNetCaret_yhat, truth=Y_Test_True)
    
    results_ls$Caret_RFranger_train_ConfMat   <- confusionMatrix((results_ls$Caret_pcaNNet_train_yhat), (trainTransformed$Class))
    results_ls$Caret_RFranger_train_ConfTab   <- table(pred=results_ls$Caret_pcaNNet_train_yhat, truth=trainTransformed$Class)
    
  }
  
  
  if(do_Caret_rangerRF) {
    # training model with SVM
    print("running Caret RF::ranger ")
    
    # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps)
    
    rangerRFCaret_classifier <- train(Class ~ ., data = trainTransformed,
                                      method = "ranger",
                                      trControl= TrainingParameters,
                                      # preProcess = preProc,
                                      na.action = na.omit)
    
    rangerRFCaret_yhat <- stats::predict(rangerRFCaret_classifier, testTransformed)
    results_ls$Caret_RFranger_train_yhat      <- stats::predict(rangerRFCaret_classifier, trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$Caret_RFranger_classifier      <- rangerRFCaret_classifier
    results_ls$Caret_RFranger_test_yhat       <- rangerRFCaret_yhat
    results_ls$Caret_RFranger_test_ConfMat    <- confusionMatrix(rangerRFCaret_yhat, (Y_Test_True))
    results_ls$Caret_RFranger_test_ConfTab    <- table(pred=rangerRFCaret_yhat, truth=Y_Test_True)
    
    results_ls$Caret_RFranger_train_ConfMat   <- confusionMatrix((results_ls$Caret_RFranger_train_yhat), (trainTransformed$Class))
    results_ls$Caret_RFranger_train_ConfTab   <- table(pred=results_ls$Caret_RFranger_train_yhat, truth=trainTransformed$Class)
    
  }
  
  
  if(do_ElasticNet){
    
    print("running ElasticNet Reg.")
    print("Elasting net is part of Garnett and its optimized .. try that ")
    # # TrainingParameters <- trainControl(method = "repeatedcv", number = NcrossVal, repeats=crossValReps,
    # #                                    returnResamp="all",
    # #                                    classProbs=TRUE)
    # 
    # # training.data.Y <- factor(as.character(training.data$Class))
    # # training.data.X <- training.data[,setdiff(colnames(training.data), c("Class"))]
    # 
    # 
    # ElasticNet_classifier <- try(train(x=trainTransformed[,setdiff(colnames(trainTransformed), c("Class"))],
    #                                    y=trainTransformed$Class,
    #                                           method = "glmnet",
    #                                           trControl= TrainingParameters,
    #                                           #preProcess = preProc,
    #                                           na.action = na.omit), silent = T)
    # #print(class(ElasticNet_classifier))
    # if(class(ElasticNet_classifier)[1]=="try-error") {
    #   print("first model failed, tring with alternate alpha and lambda...")
    #   ElasticNet_classifier <- try(caret::train(training.data.X, training.data.Y,
    #                                             method = "glmnet",
    #                                             metric = "ROC",
    #                                             trControl= TrainingParameters,
    #                                             preProcess = preProc,
    #                                             na.action = na.omit,
    #                                             tuneGrid = expand.grid(alpha = c(0.1, 1),
    #                                                                    lambda = c(0.001, 0.0001))),
    #                                silent = T)
    # }
    # if(class(ElasticNet_classifier)[1]=="try-error") {
    #   warning("ElasticNet classifier failed ")
    # } else {
    # 
    #   ElasticNet_yhat <- stats::predict(ElasticNet_classifier, testing.data)
    #   results_ls$Caret_ElasticNet_train_yhat      <- stats::predict(ElasticNet_classifier, newdata=training.data[,training.genes])
    # 
    # 
    #   results_ls$Caret_ElasticNet_classifier      <- ElasticNet_classifier
    #   results_ls$Caret_ElasticNet_test_yhat       <- ElasticNet_yhat
    #   results_ls$Caret_ElasticNet_test_ConfMat    <- confusionMatrix(ElasticNet_yhat, (Y_Test_True))
    #   results_ls$Caret_ElasticNet_test_ConfTab    <- table(pred=ElasticNet_yhat, truth=Y_Test_True)
    # 
    #   results_ls$Caret_ElasticNet_train_ConfMat   <- confusionMatrix((results_ls$Caret_ElasticNet_train_yhat), (training.data$Class))
    #   results_ls$Caret_ElasticNet_train_ConfTab   <- table(pred=results_ls$Caret_ElasticNet_train_yhat, truth=training.data$Class)
    # 
    # }
    
  }
  
  
  #### Non Caret functions
  
  if(do_Adaboost){
    
    print("running Adaboost")
    
    Adaboost_classifier <- fastAdaboost::adaboost(Class ~ ., trainTransformed, NcrossVal)
    
    y_hat_Adaboost                      <- stats::predict(Adaboost_classifier, newdata=testTransformed)
    results_ls$Adaboost_train_yhat      <- stats::predict(Adaboost_classifier, newdata=trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    results_ls$Adaboost_classifier      <- Adaboost_classifier
    results_ls$Adaboost_test_yhat       <- y_hat_Adaboost
    results_ls$Adaboost_test_ConfMat    <- confusionMatrix((y_hat_Adaboost$class), (Y_Test_True))
    results_ls$Adaboost_test_ConfTab    <- table(pred=y_hat_Adaboost$class, truth=Y_Test_True)
    
    results_ls$Adaboost_train_ConfMat   <- confusionMatrix((results_ls$Adaboost_train_yhat$class), (trainTransformed$Class))
    results_ls$Adaboost_train_ConfTab   <- table(pred=results_ls$Adaboost_train_yhat$class, truth=trainTransformed$Class)
    
    
  }
  
  
  
  
  
  
  
  if(do_ranger_RF) {
    # training
    print("running ranger-RF")
    
    RRF_classifier <- ranger::ranger(data = trainTransformed,
                                     dependent.variable.name = "Class",
                                     classification = TRUE,
                                     write.forest = TRUE)
    
    
    
    RRF_yhat <- stats::predict(RRF_classifier, testTransformed)
    results_ls$RRF_train_yhat      <- stats::predict(RRF_classifier, data=trainTransformed[,which(colnames(trainTransformed) !="Class")])
    
    
    results_ls$RRF_classifier      <- RRF_classifier
    results_ls$RRF_test_yhat       <- RRF_yhat
    results_ls$RRF_test_ConfMat    <- confusionMatrix(RRF_yhat$predictions, (Y_Test_True))
    results_ls$RRF_test_ConfTab    <- table(pred=RRF_yhat$predictions, truth=Y_Test_True)
    
    results_ls$RRF_train_ConfMat   <- confusionMatrix((results_ls$RRF_train_yhat$predictions), (trainTransformed$Class))
    results_ls$RRF_train_ConfTab   <- table(pred=results_ls$RRF_train_yhat$predictions, truth=trainTransformed$Class)
    
  }
  
  
  
  
  
  
  classifierNames <- names(results_ls)[grep("classifier", names(results_ls))]
  
  
  
  results_ls$ReSamp <- resamples(lapply(classifierNames, function(clsf){
    results_ls[[clsf]]
  }))
  
  results_ls$ReSampSum <- summary(results_ls$ReSamp)
  
  bwplot(results_ls$ReSamp)
  results_ls$ReSampbwplot <- recordPlot()
  
  dotplot(results_ls$ReSamp)
  results_ls$ReSampdotplot <- recordPlot()
  
  results_ls$modelCor <- modelCor(results_ls$ReSamp)
  
  splom(results_ls$ReSamp)
  results_ls$ReSampsplom <- recordPlot()
  
  
  
  
  return(results_ls)

  }


