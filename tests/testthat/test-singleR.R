context("Seurat")

seuratObj <- readRDS('../testdata/testSeurat.rds')

results <- 'singleR.txt'
seuratObj <- RunSingleR(seuratObj = seuratObj, resultTableFile = results)

print(table(seuratObj$SingleR_Labels))
print(table(seuratObj$SingleR_Labels_Fine))

test_that("cell type counts not as expected", {
    expect_equal(120, sum(seuratObj$SingleR_Labels == 'NK_cell'))
    expect_equal(2, sum(seuratObj$SingleR_Labels == 'B_cell'))
    expect_equal(7, sum(seuratObj$SingleR_Labels == 'Neutrophils'))

    expect_equal(188, sum(seuratObj$SingleR_Labels_Fine == 'T_cell:CD8+_Central_memory'))
    expect_equal(3, sum(seuratObj$SingleR_Labels_Fine == 'T_cell:CD4+_central_memory'))
})

test_that("cell type counts not as expected", {
    expect_equal(3886, nrow(read.table(results, sep = '\t', header = T)))
})

unlink(results)

DimPlot_SingleRClassLabs(seuratObj)

Tabulate_SingleRClassLabs(seuratObj)