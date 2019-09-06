context("scRNAseq")

test_that("Serat processing works as expected", {
  outDir <- './'
  outPrefix <- paste0(outDir, 'testData')
  resolutionToUse <- 0.6

  seuratObj <- ReadAndFilter10xData('../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', 'Set1', emptyDropNIters=5000)
  #expectedSeuratObj <- readRDS('../testdata/seuratOutputSS.rds')

  expect_equal(ncol(seuratObj), 3353, tolerance = 5)

  vgFile <- 'variableGenes.txt'
  seuratObj <- ProcessSeurat1(seuratObj, doCellCycle = T, variableGeneTable = vgFile, doCellFilter = T)

  expect_equal(file.exists(vgFile), T)
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)

  seuratObj <- FindClustersAndDimRedux(seuratObj)
  expect_equal(ncol(seuratObj), 1557)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 7)

  expect_equal(length(rownames(seuratObj@assays$RNA@scale.data)), length(rownames(seuratObj@assays$RNA@counts)))

  #Note: Seurat::PercentageFeatureSet returns 0-100.  our code is currently a fraction (0-1.0)
  expect_true(max(seuratObj$p.mito) < 1.0)
  expect_true(max(seuratObj$p.mito) > 0)

  seuratObj0 <- FindClustersAndDimRedux(seuratObj, minDimsToUse = 12, forceReCalc = T)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 7)
  rm(seuratObj0)

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

  #Note: if the expectations change, save this output as a reference:
  #seuratObjSS <- seuratObj[1:100]
  #saveRDS(seuratObjSS, file = '../testdata/seuratOutputSS.rds')
})

test_that("Serat SCTransform works as expected", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
  seuratObjSCT <- OOSAP::CreateSeuratObj(seuratData = seuratObj@assays$RNA@counts, project = 'Set1')

  seuratObjSCT <- ProcessSeurat1(seuratObjSCT, doCellCycle = F, useSCTransform = T)

  expect_equal(length(rownames(seuratObjSCT@assays$SCT@scale.data)), length(rownames(seuratObjSCT@assays$SCT@counts)))
  expect_equal(ncol(seuratObjSCT), ncol(seuratObj))

  #seuratObjSS <- seuratObjSCT[1:100]
  #saveRDS(seuratObjSS, file = '../testdata/seuratObjSCT.rds')
})