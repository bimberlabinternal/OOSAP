# OOSAP :: alpha - development

OHSU-ONPRC Single-cell Analysis Package (alpha 1)

The goal of this package is bringing together tools and piplines for the analysis of single-cell RNASeq data.

## Installation

```{r }
# install.packages("devtools")

library(devtools)
install_github("bimberlabinternal/OOSAP")

```

### Developmental progress

TODO:
    
    Make function and object names uniform across pacakge.
    limit exported functions to pipeline starting points.
    Clean up and speed up functions.
    Move towards S4 S3 structure.
    Create test cases to demo main functions.
    finalize documentation.
    


Latest Updates:

    May/28/2019 Merged with B.B. Seurat Fxs. _SERIII prefix removed. Minor documentation added.
    May/25/2019 All functions in current pipeline moved to this package.
