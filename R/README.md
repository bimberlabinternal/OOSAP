# OOSAP Developmental guidelines

Here, we house individual R files that contain the functions in this package. They are split based on utility. 

## Naming Schema

### Functions and objects

A function needs to have documentation, and we use two packages Roxygen 2 and Devtools. For version control and maintaining, between local and Github, whenever a function is changed with a brief summary. 

```{r }
library(devtools)
library(roxygen2)
```

#### Seurat functions

Seurat functions are functions that input a Seurat object, and returns an appended/changed Seurat object. All functions in this package are for Seurat 3, unless specified. 

  
  
##### Template function and documentation

```{r }

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
#' @examples
#' MyFxName(SeurObj=SO)
MyFxName <- function(SeurObj, ...){
  #Laste changes: 
  #03/01/2019: changed SeurObj$ExampleInput from NULL to "NA" :: EM
  #01/01/2019: added retrun() :: BB
  
  SeurObj$ExampleInput <- "NA"
  
  return(SeurObj)
}

```

