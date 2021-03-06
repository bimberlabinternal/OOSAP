% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_III.R
\name{HasStepRun}
\alias{HasStepRun}
\title{HasStepRun}
\usage{
HasStepRun(seuratObj, name, forceReCalc = F, printResult = T)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{name}{The name of the step to mark complete}
}
\value{
A modified Seurat object.
}
\description{
An internal method to mark steps as complete, to avoid recalculating
}
