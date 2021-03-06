% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleR_Seurat.R
\name{RunSingleR}
\alias{RunSingleR}
\title{generate SingleR object}
\usage{
RunSingleR(
  seuratObj = NULL,
  datasets = c("hpca", "blueprint", "dice", "monaco"),
  assay = NULL,
  resultTableFile = NULL,
  singlerSavePrefix = NULL,
  minFraction = 0.01,
  showHeatmap = TRUE,
  maxCellsForHeatmap = 20000
)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{datasets}{One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, and immgen. See cellDex package for available datasets.}

\item{assay}{The assay in the seurat object to use}

\item{resultTableFile}{If provided, a table of results will be saved here}

\item{singlerSavePrefix}{If provided, the SingleR results will be saved to RDS here}

\item{minFraction}{If provided, any labels present with fraction of this or fewer across cells will be converted to Unknown}

\item{showHeatmap}{If true, heatmaps will be generated showing the SingleR calls}

\item{maxCellsForHeatmap}{The heatmap will only be plotted if the total cells is below this number}
}
\value{
The modified seurat object
}
\description{
Compute SingleR classification on a Seurat object
}
