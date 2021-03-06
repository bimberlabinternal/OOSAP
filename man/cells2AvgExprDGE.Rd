% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseudoBulk.R
\name{cells2AvgExprDGE}
\alias{cells2AvgExprDGE}
\title{cells2AvgExprDGE}
\usage{
cells2AvgExprDGE(
  so = NULL,
  group.by = NULL,
  downSample = F,
  sampleNmax = NULL,
  slot = "data",
  genes = NULL
)
}
\arguments{
\item{so}{A Seurat Object}

\item{group.by}{a parameter in the Seurat Objects's meta to sample againts factor-levels}

\item{downSample}{If True random sampling is performed to return equal number of cells across group.by. if sampleNmax = NULL this is set by min cells found across all group.by factor levels}

\item{sampleNmax}{maximum number of cells to sample from, if greater than any given group.by factor level since relace = F the maximum number avaiable for that level is taken. leave as NULL if equal number of samples across all levels is desired determined by the lowest number of cells across group.by factor-levels.}

\item{slot}{Passed directly to Seurat::AverageExpression slot parameter.}

\item{genes}{A vector of genes desired, leave as NULL for all genes.}
}
\description{
cells2AvgExprDGE
}
