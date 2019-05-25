


# Functions related to Seurat but not directly in pipeline or in development



## Make Combo Seur Obj From Several Pre-processed SerObj.rds

#' @title QuickSerCombObjs
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
QuickSerCombObjs <- function(save.fig.path="./Figs",
                             working.serobjs.path="./data/10X/ser/proc", returnComboObj=F){

  require(Matrix)
  require(Seurat)

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
#' @examples
findElbow <- function(y, plot = FALSE, returnIndex = TRUE, ignore.concavity=F, min.y=NA, min.x=NA) {

  # minor modification to debug specic scenarios when fail to find elbow
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.

  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  #' @examples
  #' tmp <- findElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
  #' 	plot = TRUE) # wandering
  #' tmp <- findElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
  #' 	plot = TRUE) # late rise
  #' tmp <- findElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
  #' 	plot = TRUE) # gradual, no obvious break
  #'
  #' # Not the usual way to choose the number of PCs:
  #' library("chemometrics")
  #' data(glass)
  #' pca <- prcomp(glass)
  #' eigensum <- sum(pca$sdev * pca$sdev)
  #' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
  #' cs <- cumsum(vv)
  #' tmp <- findElbow(vv, plot = TRUE)
  #' tmp <- findElbow(cs, plot = TRUE)
  #'

  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }

  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if(any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if(any((u < 0.00001) || (u > 1))) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      if(ix > iy)  ans <- iy
      else ans <- ix
    } else {
      ## Intersecting point is on the line, use the formula
      ix <- x1 + u * (x2 - x1)
      iy <- y1 + u * (y2 - y1)
      ans <- lineMagnitude(px, py, ix, iy)
    }
    ans
  }

  # End of helper functions by PB

  ### Now for the actual findElbow function!

  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).


  y <- sort(y, decreasing = T)

  # Add an index to argument values for easy plotting
  DF <- data.frame(x = 1:length(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]

  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.

  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if(ignore.concavity) concave <- TRUE
  if (!concave) stop("Your curve doesn't appear to be concave")

  # Calculate the orthogonal distances

  if(is.na(min.x)){
    if(!is.na(min.y)){
      if(!length(which(DF$y<=min.y))<1){
        min.x = min(DF[which(DF$y<=min.y), ]$x)
      } else {
        print("min.y greater than smallest y")
        min.x = 2
      }
    } else {
      print("min.x and min.y are NA")
      min.x = 2
    }

  }

  use     <- min.x:(nrow(DF)-1)
  elbowd  <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- rep(NA, nrow(DF))
  DF$dist[use]  <- elbowd # c(NA, elbowd, NA) # first & last points don't have a distance

  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
         main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
             DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")
    points(DF$x[edm], DF$y[edm], pch = 20)
  }

  if (returnIndex) return(which.max(DF$dist)) else return(DF)

} # end of findElbow


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
ColorTheme <- function(){
  require(RColorBrewer)

  scaleyellowred <- colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
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



  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
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
#' @examples
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
#' @examples
transposedt <- function(dt, varlabel="myVar") {
  require(data.table)
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
#' @examples
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
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
#' @examples
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
#' @examples
UniformSampleDF_FacPor <- function(x, ClassF, p){
  nr <- NROW(x)
  size <- (nr * p) %/% length(unique(ClassF))
  idx <- lapply(split(seq_len(nr), ClassF), function(.x) sample(.x, size))
  unlist(idx)
}


