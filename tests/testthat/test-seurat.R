context("scRNAseq")

test_that("Serat processing works as expected", {
  outDir <- './'
  outPrefix <- paste0(outDir, 'testData')
  resolutionToUse <- 0.6

  seuratObj <- ReadAndFilter10xData('../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', 'Set1', emptyDropNIters=5000)

  expect_equal(ncol(seuratObj), 3353)

  vgFile <- 'variableGenes.txt'
  seuratObj <- ProcessSeurat1(seuratObj, doCellCycle = T, variableGeneTable = vgFile, doCellFilter = T)

  expect_equal(file.exists(vgFile), T)
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)

  seuratObj <- FindClustersAndDimRedux(seuratObj)
  expect_equal(ncol(seuratObj), 1557)

  unlink(vgFile)

  mf <- paste0(outPrefix, '.markers.txt')
  md <- paste0(outPrefix, '.markers.rds')
  Find_Markers(seuratObj, resolutionToUse = resolutionToUse, outFile = mf, saveFileMarkers = md, testsToUse = c('wilcox', 't'))

  expect_equal(nrow(utils::read.table(mf, sep = '\t', header = T)), 201, tolerance = 5)

  unlink(md)
  unlink(mf)

  sf <- paste0(outPrefix, '.summary.txt')
  WriteSummaryMetrics(seuratObj, file = sf)

  expect_equal(nrow(utils::read.table(sf, sep = '\t', header = T)), 2)

  unlink(sf)

  dr <- paste(outPrefix, ".DimReduxComps.csv", sep="")
  SaveDimRedux(seuratObj, file = dr)
  expect_equal(nrow(utils::read.table(dr, sep = '\t', header = T)), ncol(seuratObj))

  unlink(dr)
})