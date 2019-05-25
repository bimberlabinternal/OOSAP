# OOSAP

OHSU-ONPRC Single-cell Analysis Package (alpha 1)

## Developmental guidelines

### Naming Schema

#### Functions and objects

##### Seurat functions

Seurat functions are functions that input a Seurat object, and returns an appended/changed Seurat object.

-suffix _SERII: refers to functions that only work on Seurat V 2.x These are legacy and are (need to be) changed to work for V 3.x

-suffix _SERIII: refers to functions that only work on Seurat V 3.x; since future versions may arrive (usually sooner than expected)

```{r }
# template function

MyFxName_SER## <- function(SeurObj, ...){

  return(SeurObj)
}

```

