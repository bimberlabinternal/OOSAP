context("scRNAseq")

test_that("SingleR works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    results <- 'singleR.txt'
    singleRPrefix <- 'singleR.results'
    seuratObj <- RunSingleR(seuratObj = seuratObj, resultTableFile = results, singlerSavePrefix = singleRPrefix)

    print(table(seuratObj$SingleR_Labels))
    print(table(seuratObj$SingleR_Labels_Fine))

    sr1 <- paste0(singleRPrefix, '.singleR.rds')
    expect_true(file.exists(sr1))
    unlink(sr1)

    sr2 <- paste0(singleRPrefix, '.singleR.fine.rds')
    expect_true(file.exists(sr2))
    unlink(sr2)

    expect_equal(101, sum(seuratObj$SingleR_Labels == 'NK_cell'))
    expect_equal(0, sum(seuratObj$SingleR_Labels == 'B_cell'))
    expect_equal(0, sum(seuratObj$SingleR_Labels == 'Neutrophils'))
    expect_equal(22, sum(seuratObj$SingleR_Labels == 'Unknown'))

    expect_equal(170, sum(seuratObj$SingleR_Labels_Fine == 'T_cell:CD8+_Central_memory'))
    expect_equal(68, sum(seuratObj$SingleR_Labels_Fine == 'T_cell:CD8+_naive'))

    expect_equal(ncol(seuratObj), nrow(read.table(results, sep = '\t', header = T)))

    unlink(results)

    DimPlot_SingleRClassLabs(seuratObj)

    Tabulate_SingleRClassLabs(seuratObj)

    #saveRDS(seuratObj, file = '../testdata/singleR.rds')
})

