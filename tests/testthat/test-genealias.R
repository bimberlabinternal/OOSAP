context("scRNAseq")

test_that("All aliases preserved", {
    homologs <- QueryEnsemblHumanHomologs(c('ENSMMUG00000040244', 'CD8'))
    expect_equal(homologs$Label, c('TRAV-1', 'CD8'))

    cdGenes <- RenameGenesUsingCD(c('PTPR', '12345', 'DPP4'))
    expect_equal(cdGenes, c('PTPR', '12345', 'DPP4 (CD26,ADCP2)'))

})