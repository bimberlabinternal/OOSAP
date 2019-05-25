# OOSAP Developmental guidelines

Here, we house individual R files that contain the functions in this package. They are split based on utility. 

## Naming Schema

### Functions and objects

A function needs to have documentation, and we use two packages Roxygen 2 and Devtools,

```{r }
library(devtools)
library(roxygen2)
```

#### Seurat functions

Seurat functions are functions that input a Seurat object, and returns an appended/changed Seurat object.



  -suffix _SERII: refers to functions that only work on Seurat V 2.x These are legacy and are (need to be) changed to work for V 3.x

  -suffix _SERIII: refers to functions that only work on Seurat V 3.x; since future versions may arrive (usually sooner than expected)
  
  
##### Template function and documentation

```{r }


#' Function Name.
#' 
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @examples
#' MyFxName_SER(SeurObj=SO)

MyFxName_SERIII <- function(SeurObj, ...){
  #Laste change date: 01/01/2019

  return(SeurObj)
}

```

