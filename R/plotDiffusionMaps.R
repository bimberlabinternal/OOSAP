
#' @title Plots destiny's Diffusion Map
#' @param DiffMap cell names
#' @param ndimsToPlot number of dimensions to plot
#' @import patchwork
#' @export
#' 

PlotDiffComp <- function(diffMap, ndimsToPlot){
  dcPlots <- list()
  for(i in 1:ndimsToPlot){
    if((i %% 2) != 0){
      plotData <- data.frame(DCx = destiny::eigenvectors(diffMap)[, i], 
                             DCy = destiny::eigenvectors(diffMap)[, i+1],
                             var1 = rownames(diffMap@eigenvectors))
      
      #colnames(var1)["var1"] <- variable_col
      
      plot1 <- ggplot(plotData)+
        aes(x = DCx, y = DCy, colour = var1)+
        geom_point(alpha=.2)  + 
        xlab(paste0("Diffusion component ", i)) + 
        ylab(paste0("Diffusion component ", i+1)) +
        theme_bw()+ 
        #scale_color_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))
      
      
      print(plot1)
      dcPlots[[paste0("p", i, ".", i+1)]] <- plot1
      
    }
  }
  print(patchwork::wrap_plots(dcPlots), ncol = 2)
  
}


#' @title Plots destiny's Diffusion Map. This function includes subplots that shows the trajectory of the cells by cell Ids on each componenets. 
#' @param DiffMap cell names
#' @param ndimsToPlot number of dimensions to plot
#' @import patchwork
#' @export
#' 

PlotDiffComp2 <- function(diffMap, ndimsToPlot){
  #dcPlots <- list()
  for(i in 1:ndimsToPlot){
    if((i %% 2) != 0){
      plotData <- data.frame(DCx = destiny::eigenvectors(diffMap)[,i],
                             DCy = destiny::eigenvectors(diffMap)[,i+1],
                             var1 = factor(rownames(diffMap@eigenvectors), 
                                           levels = gtools::mixedsort(levels(factor(rownames(diffMap@eigenvectors))))
                             ))
      
      #colnames(var1)["var1"] <- variable_col
      
      plot1 <- ggplot(plotData)+
        aes(x = DCx, y = DCy, colour = var1)+
        geom_point(alpha=.2)  + 
        xlab(paste0("Diffusion component ", i)) + 
        ylab(paste0("Diffusion component ", i+1)) +
        theme_bw()+ 
        #scale_color_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))
      
      plot1x <- ggplot(plotData)+
        aes(x = rank(DCx), y = var1, colour = var1)+
        ggbeeswarm::geom_quasirandom(groupOnX = FALSE)+ 
        xlab(paste0("Diffusion component ", i))+ 
        ylab('Timepoint')+
        theme_bw()+ 
        #scale_color_manual(values=col_vector) +
        theme(legend.position="none",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())
      
      plot1y <- ggplot(plotData)+
        aes(x = rank(DCy), y = var1, colour = var1)+
        ggbeeswarm::geom_quasirandom(groupOnX = FALSE)+ 
        xlab(paste0("Diffusion component ", i+1))+ 
        ylab('Timepoint')+
        theme_bw()+ 
        #scale_color_manual(values=col_vector) +
        theme(legend.position="none",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank())
      
      pGrid <- plot1 / (plot1x | plot1y) 
      
      print(pGrid)
    }
  }
}
