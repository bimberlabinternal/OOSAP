context("Seurat")

outDir <- './'
outPrefix <- paste0(outDir, 'testData')
resolutionToUse <- 0.6

seuratObj <- ReadAndFilter10xData('../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', 'Set1', emptyDropNIters=5000)

test_that("cell count correct", {
  expect_equal(ncol(seuratObj), 3353)
})

vgFile <- 'variableGenes.txt'
seuratObj <- ProcessSeurat1(seuratObj, doCellCycle = T, variableGeneTable = vgFile, doCellFilter = T)

test_that("variable genes not saved", {
  expect_equal(file.exists(vgFile), T)
})

test_that("variable gene list not expected length", {
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)
})

seuratObj <- FindClustersAndDimRedux(seuratObj)

test_that("cell count correct after filter", {
  expect_equal(ncol(seuratObj), 1557)
})

unlink(vgFile)

mf <- paste0(outPrefix, '.markers.txt')
md <- paste0(outPrefix, '.markers.rds')
Find_Markers(seuratObj, resolutionToUse = resolutionToUse, outFile = mf, saveFileMarkers = md, testsToUse = c('wilcox', 't'))

test_that("marker list not expected length", {
  expect_equal(nrow(utils::read.table(mf, sep = '\t', header = T)), 201)
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

test_that("dim redux file not expected length", {
  expect_equal(nrow(utils::read.table(dr, sep = '\t', header = T)), ncol(seuratObj))
})

unlink(dr)

# Test SingleR:
singleRFile <- paste0(outPrefix, '.singleR.rds')
seuratObj <- SingleRmySerObj(SeurObj = seuratObj, SavePath = singleRFile, ReturnSeurObj = T)

test_that("cell type counts not as expected", {
  expect_equal(94, sum(seuratObj$SingleR_Labels1 == 'NK_cell'))
  expect_equal(13, sum(seuratObj$SingleR_Labels1 == 'T_cell:CD4+_central_memory'))

  expect_equal(465, sum(seuratObj$SingleR_Labels2 == 'CD8+ Tcm'))
  expect_equal(911, sum(seuratObj$SingleR_Labels2 == 'CD8+ Tem'))
})

unlink(singleRFile)
