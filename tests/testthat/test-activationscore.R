context("scRNAseq")


test_that("Highly Activated Cells Called", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')

  ret <- OOSAP:::ClassifySeurat(
      seuratObj,
      trainedClassifiers = NULL,
      doGarnettClassify = F,
      moduleScoreGeneLists = OOSAP::Phenotyping_GeneList(),
      savePrefix = './testSeurat'
  )
  
  f <- 'testSeurat.SGS.csv'
  df <- read.table(f, header = T, sep = ',')
  expect_equal(names(ret), c('SeuratGeneScore'))
  expect_equal(nrow(df), ncol(seuratObj))
  expect_equal(ncol(df), 74)
  expect_equal(max(df['HighlyActivated.HighlyActivated1']), 1.620322)
  expect_equal(min(df['HighlyActivated.HighlyActivated1']), -0.4411742)

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