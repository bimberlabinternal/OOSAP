
#' @import RColorBrewer


#' @title A Title
#'
#' @description A description
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





