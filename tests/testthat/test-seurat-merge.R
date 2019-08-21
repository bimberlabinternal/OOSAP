context("scRNAseq")

test_that("Seurat-merge works as expected", {
  data <- list(
  'Set1' = '../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/10xCounts/CellRanger3/raw_feature_bc_matrix'
  )

  seuratObjs <- list()
  for (datasetName in names(data)) {
    seuratObjs[[datasetName]] <- ReadAndFilter10xData(testthat::test_path(data[[datasetName]]), datasetName, emptyDropNIters=1000)
  }

  expect_equal(length(seuratObjs), 2)

  seuratObj <- MergeSeuratObjs(seuratObjs, data)
  rm(seuratObjs)

  # NOTE: this might not be deterministic.  expect about 7987
  expect_equal(7987, ncol(seuratObj), tolerance = 10)
})