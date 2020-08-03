context("scRNAseq")

test_that("Highly Activated Cells Called", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
  
  f <- './testSeuratHA.SGS.csv'
  ret <- ClassifySGS(
    seuratObj,
    geneSetName = 'HighlyActivated',
    geneList = OOSAP::Phenotyping_GeneList()[['HighlyActivated']],
    positivityThreshold = 0.5,
    saveFilePath = f, doPlot = T
  )
  
  expect_equal(names(ret), c("CellBarcode", "HighlyActivated.Score", "HighlyActivated.Call"))
  
  df <- read.table(f, header = T, sep = ',')
  expect_equal(names(ret), c("CellBarcode", "HighlyActivated.Score", "HighlyActivated.Call"))
  expect_equal(nrow(df), ncol(seuratObj))
  expect_equal(ncol(df), 3)
  expect_equal(max(df['HighlyActivated.Score']), 1.26, tolerance = 0.001)
  expect_equal(min(df['HighlyActivated.Score']), -0.279, tolerance = 0.001)
  expect_equal(sum(df[['HighlyActivated.Call']]), 58)
  
  unlink(f)

  #Print single FeaturePlot
  seuratObj2 <- ClassifySGSAndApply(
    seuratObj,
    geneSetName = 'HighlyActivated',
    geneList = OOSAP::Phenotyping_GeneList()[['HighlyActivated']],
    doPlot = T
  )
  
  #This will print two FeaturePlots using plot_grid
  seuratObj2 <- ClassifySGSAndApply(
    seuratObj,
    geneSetName = 'HighlyActivated',
    geneList = OOSAP::Phenotyping_GeneList()[['HighlyActivated']],
    positivityThreshold = 0.5,
    doPlot = T
  )
  
  ts <- seuratObj2$HighlyActivated.Score
  names(ts) <- NULL
  expect_equal(df[['HighlyActivated.Score']], ts)

  # This should write metrics to a file:
  sf <- 'summary.txt'
  WriteSummaryMetrics(seuratObj2, file = sf)
  df <- utils::read.table(sf, sep = '\t', header = T)
  head(df)
  expect_equal(nrow(df), 3)

  totalActivated <- sum(seuratObj2$HighlyActivated.Call)
  print(paste0('Total activated: ', totalActivated))
  expect_equal(totalActivated, 58)

  fractionActivated <- totalActivated / ncol(seuratObj2)
  print(paste0('Fraction activated: ', fractionActivated))
  expect_equal(fractionActivated, 0.0373, tolerance = 0.001)
  expect_equal(df$Value[df$MetricName == 'FractionActivated'], fractionActivated, tolerance = 0.001)

  unlink(sf)

})


test_that("Multi-gene set scoring works", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')

  genes <- OOSAP::Phenotyping_GeneList()
  genes <- genes[names(genes) != 'Testis']
  
  f <- './testSeurat.SGS.csv'
  ret <- ClassifySGS_Multiple(
      seuratObj,
      moduleScoreGeneLists = genes,
      saveFilePath = f
  )
  
  df <- read.table(f, header = T, sep = ',')
  expect_equal(ncol(ret), length(genes) + 1)
  expect_equal(nrow(df), ncol(seuratObj))
  expect_equal(max(df['Dendritic.Score']), 5.480072, tolerance = 0.0001)
  expect_equal(min(df['Dendritic.Score']), -1.229455, tolerance = 0.0001)

  unlink(f)
})


# test_that("MCE Runs Correctly", {
#     library(OOSAP)
#     setwd('/projects/oosap/tests/testthat')
#     seuratObj <- readRDS('../testdata/seuratOutput.rds')
#
#     ret <- OOSAP:::ClassifySeurat(
#         seuratObj,
#         trainedClassifiers = list(
#             B = '/projects/oosap/classifiers/MCR_LS_B.rds',
#             CD4T = '/projects/oosap/classifiers/MCR_LS_CD4T.rds',
#             CD8T = '/projects/oosap/classifiers/MCR_LS_CD8T.rds',
#             Lymph = '/projects/oosap/classifiers/MCR_LS_Lymph.rds',
#             NK = '/projects/oosap/classifiers/MCR_LS_NK.rds'
#         ),
#         doGarnettClassify = F,
#         savePrefix = './testSeurat2'
#     )
# })