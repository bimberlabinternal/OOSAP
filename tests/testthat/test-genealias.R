context("scRNAseq")

test_that("All aliases preserved", {
    # ENSMMUG00000040244: resolved by external_gene_name
    # ENSMMUG00000008350: resolved by homolog
    # CD8: Already  symbol, not resolved, return original
    homologs <- QueryEnsemblSymbolAndHumanHomologs(c('ENSMMUG00000040244', 'ENSMMUG00000008350', 'CD8'), dataset = 'mmulatta_gene_ensembl', version = '97')
    expect_equal(homologs$Label, c('TRAV1-1', 'MDK(hs)', 'CD8'))

    # This is a private method, but test directly anyway
    cdg <- OOSAP:::RenameGenesUsingCD(c('PTPR', '12345', 'DPP4', 'ENSMMUG00000040244'))
    expect_equal(cdg, c('PTPR', '12345', 'DPP4 (CD26,ADCP2)', 'ENSMMUG00000040244'))

    aliased <- AliasGeneNames(c('ENSMMUG00000040244', 'ENSMMUG00000008350', 'CD8', 'DPP4'), ensemblVersion = '97')
    expect_equal(aliased, c('TRAV1-1', 'MDK(hs)', 'CD8', 'DPP4 (CD26,ADCP2)'))

    #verify concat works when it returns two hits:
    aliased <- AliasGeneNames(c('ENSMMUG00000029821'), ensemblVersion = '97')
    expect_equal(aliased, c('HSPA1A,HSPA1B(hs)'))
})