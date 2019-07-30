context("Seurat")

data <- list(
  'Set1' = '../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/10xCounts/CellRanger3/raw_feature_bc_matrix'
)

seuratObjs <- list()
for (datasetName in names(data)) {
  seuratObjs[[datasetName]] <- ReadAndFilter10xData(testthat::test_path(data[[datasetName]]), datasetName, emptyDropNIters=1000)
}

test_that("object count correct", {
  expect_equal(length(seuratObjs), 2)
})

seuratObj <- MergeSeuratObjs(seuratObjs, data)
rm(seuratObjs)

# NOTE: this might not be deterministic.  expect about 7987
cellDiff <- abs(ncol(seuratObj) - 7987)
test_that("cell count within tolerances", {
  expect_true(cellDiff < 10)
})