

#' @title bimodal_posDist1D
#' @description 1D bimodal Dist given parameters
#' @export
bimodal_posDist1D <- function (n=1000, cpct = 0.4, mu1 = log(1), mu2 = log(100), sig1 = log(3), sig2=log(2), p=.5) {
  # BMD <- bimodal_posDist1D(n=1000, cpct = 0.4, mu1 = log(1), mu2 = log(100), sig1 = log(3), sig2=log(2))
  # hist(log(BMD), breaks=100)
  
  y0 <- rlnorm(n*p,mean=mu1, sd = sig1)
  y1 <- rlnorm(n*(1-p),mean=mu2, sd = sig2)
  
  # print(min(c(y0, y1)))
  
  if (min(c(y0, y1))<1) {
   
    y0 <- y0 + abs(min(c(y0, y1))) + 3
    y1 <- y1 + abs(min(c(y0, y1))) + 3
  }
  
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag 
}

#' @title bimodal_posDisHD
#' @description ndim-HD bimodal Dist given parameters
#' @export
bimodal_posDisHD <- function(ndim = 5, n=1000, cpct = 0.4, mu1 = log(1), mu2 = log(100), sig1 = log(3), 
                             sig2=log(2), p=.5, JitterPeaks = T, 
                             removeDups = T, tSNEplot = T){
  
  if(JitterPeaks) {
    tempX <- as.data.frame( lapply(1:ndim, function(xN){
      bimodal_posDist1D(n=n, cpct = cpct, mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2=sig2, p=p)
    }), col.names = paste0("col", 1:ndim))
    
  } else {
    tempX <- as.data.frame( lapply(1:ndim, function(xN){
      bimodal_posDist1D(n=n, cpct = cpct, mu1 = round(mean(rnorm(1000, mu1, 3))), 
                        mu2 = round(mean(rnorm(1000, mu2, 3))), 
                        sig1 = sig1, 
                        sig2=sig2, 
                        p=p)
    }), col.names = paste0("col", 1:ndim))
  }
  
  if(removeDups) tempX <- tempX[!duplicated(tempX),]
  if(tSNEplot) plot(Rtsne::Rtsne(tempX,)$Y)
  
  return(tempX)

 
}


# tempDFX   <- bimodal_posDisHD(tSNEplot = F)
# tempMeans <- apply(tempDFX, 2, mean)
# tempSds   <- apply(tempDFX, 2, sd)
# 
# tempDFX$class <- "neg"
# 
# # tempPosID <- rownames(subset(tempDFX, 
# #                                col1 < tempMeans[1] + tempSds[1] & 
# #                                col2 < tempMeans[2] + tempSds[2] & 
# #                                col3 > tempMeans[3] + tempSds[3] & 
# #                                col4 > tempMeans[4] + tempSds[4] & 
# #                                col5 > tempMeans[5] + tempSds[5]))
# 
# tempPosID <- rownames(subset(tempDFX, 
#                               col1 < tempMeans[1] &#+ tempSds[1] & 
#                                col2 > tempMeans[2] #+ tempSds[2] #& 
#                                #col3 > tempMeans[3] + tempSds[3]# & 
#                                #col4 > tempMeans[4] + tempSds[4] & 
#                                #col5 > tempMeans[5] + tempSds[5]
#                                ))
# 
# 
# tempDFX[tempPosID,]$class <- "pos"
# 
# plot(Rtsne::Rtsne(tempDFX[,1:5],)$Y, col = factor(tempDFX$class), pch=20)

