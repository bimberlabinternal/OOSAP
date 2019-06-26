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

vgFile <- 'variableGenes.txt'
seuratObj <- ProcessSeurat1(seuratObj, doCellCycle = T, variableGeneTable = vgFile)

test_that("variable genes not saved", {
  expect_equal(file.exists(vgFile), T)
})

test_that("variable gene list not expected length", {
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)
})

seuratObj <- FindClustersAndDimRedux(seuratObj)

test_that("cell count correct", {
  expect_equal(ncol(seuratObj), 7985)
})

unlink(vgFile)

mf <- paste0(outPrefix, '.markers.txt')
md <- paste0(outPrefix, '.markers.rds')
FindMarkers(seuratObj, resolutionToUse, outFile = mf, saveFileMarkers = md, testsToUse = c('wilcox', 't'))

test_that("marker list not expected length", {
  expect_equal(nrow(utils::read.table(mf, sep = '\t', header = T)), 1065)
})

unlink(md)
unlink(mf)

sf <- paste0(outPrefix, '.summary.txt')
WriteSummaryMetrics(seuratObj, file = sf)

test_that("summary file not expected length", {
  expect_equal(nrow(utils::read.table(sf, sep = '\t', header = T)), 2)
})

unlink(sf)

dr <- paste(outPrefix, ".DimReduxComps.csv", sep="")
SaveDimRedux(seuratObj, file = dr)

