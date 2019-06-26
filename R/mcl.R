#' @import Seurat


# Functions related to Seurat but not directly in pipeline or in development



## Make Combo Seur Obj From Several Pre-processed SerObj.rds

#' @title QuickSerCombObjs
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
QuickSerCombObjs <- function(save.fig.path="./Figs",
                             working.serobjs.path="./data/10X/ser/proc", returnComboObj=F){

  CleaningLS <- list()
  CleaningLS$SeurObjs <- list()
  getwd()



  if(!dir.exists(save.fig.path)) dir.create(save.fig.path, recursive = T)

  CombObj.path <- paste(working.serobjs.path, "/SeurComboObj.rds", sep="")

  if(!file.exists( CombObj.path)){

    CleaningLS$all_Seurat_files  <- list.files(paste(getwd(), working.serobjs.path, sep=""),
                                               full.names = T,
                                               pattern = "SeuratObj.rds")

    CleaningLS$SeurObj.rds_files <-  CleaningLS$all_Seurat_files[grep("SeuratObj.rds", CleaningLS$all_Seurat_files)]



    SerObj.files <- CleaningLS$SeurObj.rds_files#[1:4]

    for(fN in 1:length(SerObj.files)){
      # fN=1
      print(fN)
      print("reading in Ser obj")
      exptNum <- gsub("-", "_", gsub("_SeuratObj.rds", "", basename(SerObj.files[fN])))

      seuratObjs <- readRDS(SerObj.files[fN])

      seuratObjs[['OrigFileName']] <- c(exptNum)

      CleaningLS$SeurObjs[[fN]] <- seuratObjs
    }; remove(seuratObjs)

    print("all loaded in now preping for merging..")
    length(CleaningLS$SeurObjs)

    names(CleaningLS$SeurObjs) <- paste("ID", as.character(unlist(lapply(SerObj.files, function(xN){

      gsub("-", "_", gsub("_SeuratObj.rds", "", basename(xN)))
    }))), sep="_")


    TempA <- CleaningLS$SeurObjs[[names(CleaningLS$SeurObjs)[1]]]
    TempB <- c(CleaningLS$SeurObjs[[names(CleaningLS$SeurObjs)[1]]])

    for(ij in 3:length(CleaningLS$SeurObjs)){
      #ij = 2
      TempB <- append(TempB, CleaningLS$SeurObjs[[ij]])
    }


    print("merging ser objs")
    SeurComboObj <- merge(x = TempA,
                          y = TempB,
                          add.cell.ids = gsub("_", "", gsub("-", "_", gsub("\\.", "-", names(CleaningLS$SeurObjs)))),
                          do.normalize = F, project = "214")

    print("completed merging ser objs, now basic pre-processing")
    SeurComboObj <- NormalizeData(object = SeurComboObj)
    SeurComboObj <- FindVariableFeatures(object = SeurComboObj)
    SeurComboObj <- ScaleData(object = SeurComboObj)
    SeurComboObj <- RunPCA(object = SeurComboObj)
    SeurComboObj <- FindNeighbors(object = SeurComboObj)
    SeurComboObj <- FindClusters(object = SeurComboObj)
    #SeurComboObj <- RunTSNE(object = SeurComboObj, check_duplicates = F, verbose=T)
    #SeurComboObj <- RunUMAP(object = SeurComboObj)

    print("saving combo obj")
    saveRDS(SeurComboObj, CombObj.path)
  } else {
    print("Already Made.... delete or rename to remake")
    if(returnComboObj) SeurComboObj <- readRDS(CombObj.path)
  }
  if(returnComboObj) return(SeurComboObj)
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
ColorTheme <- function(){
  scaleyellowred <- grDevices::colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)

  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[-4]

  col_vector2 <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52",
                   "#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                   '#e6194b', '#3cb44b', '#ffe119', '#4363d8',
                   '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                   '#bcf60c', '#fabebe', '#008080', '#e6beff',
                   '#9a6324', '#fffac8', '#800000', '#aaffc3',
                   '#808000', '#ffd8b1', '#000075', '#808080',
                   '#ffffff', '#000000')


  col_vector2 <- as.character(sapply(1:20, function(xN){
    c(col_vector[xN], col_vector2[xN])
  }))
  col_vector  <- c(col_vector2, col_vector[21:length(col_vector)])

  col_vector <- col_vector[c(1:4, 6:length(col_vector))]

  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}

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
quickTabulate <- function(spMat){
  histo_numers <- matrix(c(0:max(spMat), rep(0, max(spMat)+1)), ncol = 2)
  histo_numers[1:max(spMat)+1, 2] <- tabulate(as.matrix(spMat))
  histo_numers[1, 2] <- sum(spMat == 0)
  return(histo_numers)
}
