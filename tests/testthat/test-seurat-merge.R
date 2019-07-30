skip_on_cran()

context("Seurat")

data <- list(
  'Set1' = '../testdata/10XCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/10xCounts/CellRanger3/raw_feature_bc_matrix'
)

print('working dir:')
print(getwd())
print(list.files(getwd()))
print(list.files('../'))

seuratObjs <- list()
for (datasetName in names(data)) {
  print(paste0('using path: ', testthat::test_path(data[[datasetName]])))
  seuratObjs[[datasetName]] <- ReadAndFilter10xData(testthat::test_path(data[[datasetName]]), datasetName)
}

test_that("object count correct", {
  expect_equal(length(seuratObjs), 2)
})


seuratObj <- MergeSeuratObjs(seuratObjs, data)
rm(seuratObjs)

# NOTE: this might not be deterministic
test_that("cell count correct", {
  expect_equal(ncol(seuratObj), 7987)
})