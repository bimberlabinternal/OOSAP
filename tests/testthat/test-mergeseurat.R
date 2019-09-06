context("scRNAseq")

test_that("merge seurat works as expected on panc8 data", {
  # #Testing with panc8 data
  # options(future.globals.maxSize = 10000 * 1024^2)
  # library(Seurat)
  # library(SeuratData)
  # library(ggplot2)
  # data("panc8")
  # seuratObjs <- SplitObject(panc8, split.by = "tech")
  # method = "simple" # simple or cca
  # useAllFeatures = F
  # nVariableFeatures = 2000
  # includeCellCycleGenes = T
  # assay = NULL
  # cc.genes <- OOSAP:::cc.genes
  # g2m.genes.orig <- OOSAP:::g2m.genes.orig
  # spike.genes = unique(c(OOSAP::Phenotyping_GeneList()$TCellCanonical,
  #                      OOSAP::Phenotyping_GeneList()$TCellSecondary,
  #                      OOSAP::Phenotyping_GeneList()$TCellTranscription,
  #                      OOSAP::Phenotyping_GeneList()$CD8Canonical,
  #                      OOSAP::Phenotyping_GeneList()$MAIT,
  #                      OOSAP::Phenotyping_GeneList()$CD8Subphenos1,
  #                      OOSAP::Phenotyping_GeneList()$CD4Canonical,
  #                      OOSAP::Phenotyping_GeneList()$CD4Subphenos1,
  #                      OOSAP::Phenotyping_GeneList()$NKCanonical,
  #                      OOSAP::Phenotyping_GeneList()$HighlyActivated))
  # metadata = NULL
  # maxCCAspaceDim = 20
  # normalization.method = "LogNormalize"
  # maxPCs2Weight = 20
    #saveRDS(seuratObj, file = '../testdata/singleR.rds')
})

