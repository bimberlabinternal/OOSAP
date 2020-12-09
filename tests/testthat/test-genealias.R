context("scRNAseq")

test_that("All aliases preserved", {
    biomaRt::biomartCacheClear()

    # ENSMMUG00000040244: resolved by external_gene_name
    # ENSMMUG00000008350: resolved by homolog
    # CD8: Already  symbol, not resolved, return original
    homologs <- QueryEnsemblSymbolAndHumanHomologs(c('ENSMMUG00000042504', 'ENSMMUG00000008350', 'CD8'), dataset = 'mmulatta_gene_ensembl')
    expect_equal(homologs$Label, c('APOBEC3A', 'ENSMMUG00000008350', 'CD8'))

    #This will fail the regexp and should not get queried
    homologs <- QueryEnsemblSymbolAndHumanHomologs(c('foo'), dataset = 'mmulatta_gene_ensembl')
    expect_equal(homologs$Label, c('foo'))

    # This is a private method, but test directly anyway
    cdg <- OOSAP:::RenameGenesUsingCD(c('PTPR', '12345', 'DPP4', 'ENSMMUG00000040244'))
    expect_equal(cdg, c('PTPR', '12345', 'DPP4 (CD26,ADCP2)', 'ENSMMUG00000040244'))

    aliased <- AliasGeneNames(c('ENSMMUG00000042504', 'CD8', 'DPP4'))
    expect_equal(aliased, c('APOBEC3A', 'CD8', 'DPP4 (CD26,ADCP2)'))

    #TODO: not working on newest gene build. Need another example.
    #verify concat works when it returns two hits:
    #aliased <- AliasGeneNames(c('ENSMMUG00000041076'))
    #expect_equal(aliased, c('HSPA1A,HSPA1B(hs)'))

    biomaRt::biomartCacheClear()
})