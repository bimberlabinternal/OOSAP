skip_on_cran()

context("Seurat")

outDir <- './'
outPrefix <- paste0(outDir, 'testData')
resolutionToUse <- 0.6

data <- list(
  'Set1' = '../testdata/10XCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/10xCounts/CellRanger3/raw_feature_bc_matrix'
)

seuratObjs <- list()
for (datasetName in names(data)) {
  seuratObjs[[datasetName]] <- ReadAndFilter10xData(data[[datasetName]], datasetName)
}

test_that("object count correct", {
  expect_equal(length(seuratObjs), 2)
})


seuratObj <- MergeSeuratObjs(seuratObjs, data)
rm(seuratObjs)

test_that("cell count correct", {
  expect_equal(ncol(seuratObj), 7985)
})