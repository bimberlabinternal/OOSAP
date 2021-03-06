% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_III.R
\name{MarkStepRun}
\alias{MarkStepRun}
\title{MarkStepRun}
\usage{
MarkStepRun(seuratObj, name, saveFile = NULL)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{name}{The name of the step to test}

\item{saveFile}{If provided, the seurat object will be saved as RDS to this location}
}
\value{
A modified Seurat object.
}
\description{
An internal method to determine is a step has been marked as complete
}
