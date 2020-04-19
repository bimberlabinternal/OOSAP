context("scRNAseq")

test_that("Seurat-merge works as expected", {
  data <- list(
  'Set1' = '../testdata/10xCounts/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/10xCounts/CellRanger3/raw_feature_bc_matrix'
  )

  seuratObjs <- list()
  for (datasetName in names(data)) {
    seuratObjs[[datasetName]] <- ReadAndFilter10xData(testthat::test_path(data[[datasetName]]), datasetName, emptyDropNIters=1000)
  }

  expect_equal(length(seuratObjs), 2)

  #Simple method
  seuratObj <- MergeSeuratObjs(seuratObjs, data, method = NULL)
  expect_equal("RNA", Seurat::DefaultAssay(seuratObj))
  expect_equal("simple", seuratObj@misc$MergeMethod)

  #Gene IDs preserved:
  expect_equal(nrow(seuratObj), length(seuratObj@assays$RNA@meta.features$GeneId))
  geneIds <- GetGeneIds(seuratObj, c('HES4', 'CALML6'))
  names(geneIds) <- NULL
  expect_equal(geneIds, c('ENSMMUG00000001817', 'ENSMMUG00000012392'))

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)

  # NOTE: this might not be deterministic.  expect about 7987
  print(paste0('cells: ', ncol(seuratObj)))
  expect_equal(7987, ncol(seuratObj), tolerance = 10)
  
  #Invalid method
  expect_error(MergeSeuratObjs(seuratObjs, data, method = 'bad'))
               
  #CCA method
  seuratObj <- MergeSeuratObjs(seuratObjs, NULL, method = 'cca')
  expect_equal("Integrated", Seurat::DefaultAssay(seuratObj))
  expect_equal("cca", seuratObj@misc$MergeMethod)
  expect_equal(7987, ncol(seuratObj), tolerance = 10)
  expect_equal(nrow(seuratObj), 2354)

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)

  #CCA method + spike genes
  spikeGenes <- unique(c(OOSAP::Phenotyping_GeneList()$TCellCanonical,
                         OOSAP::Phenotyping_GeneList()$TCellSecondary,
                         OOSAP::Phenotyping_GeneList()$TCellTranscription,
                         OOSAP::Phenotyping_GeneList()$CD8Canonical,
                         OOSAP::Phenotyping_GeneList()$MAIT,
                         OOSAP::Phenotyping_GeneList()$CD8Subphenos1,
                         OOSAP::Phenotyping_GeneList()$CD4Canonical,
                         OOSAP::Phenotyping_GeneList()$CD4Subphenos1,
                         OOSAP::Phenotyping_GeneList()$NKCanonical,
                         OOSAP::Phenotyping_GeneList()$HighlyActivated))
  
  expectedSpike = intersect(rownames(seuratObjs[[1]]), spikeGenes)

  #Not all genes will be present:
  expect_equal(length(intersect(rownames(seuratObj), expectedSpike)), 44)


  #Repeat with spike genes
  seuratObj <- MergeSeuratObjs(seuratObjs, NULL, method = 'cca', spike.genes = spikeGenes)

  expect_equal(nrow(seuratObj), 2358)
  
  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)
  
  #All spike genes should be present
  expect_equal(sort(intersect(rownames(seuratObj), expectedSpike)), sort(expectedSpike))
  
})

