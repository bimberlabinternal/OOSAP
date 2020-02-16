context("scRNAseq")

test1a =
  c("ENSMMUG00000000016",
     "ENSMMUG00000015664",
     "ENSMMUG00000010387",
     "ENSMMUG00000060218",
     "ENSMMUG00000000065",
     "ENSMMUG00000059499",
     "ENSMMUG00000050963",
     "ENSMMUG00000059669",
     "ENSMMUG00000022935")

test1b = 
 c("ADSS",
   "ZBTB18",
   "BRP44",
   "ENSMMUG00000060218",
   "CCDC181",
   "TNLG2A",
   "SNORD24",
   "COP1",
   "BRINP2")

test2a = 
  c("ENSMMUG00000013380",
    "ENSMMUG00000050963",
    "ENSMMUG00000059669")

test2b = 
 c("MROH9",
   "SNORD24",
   "COP1")

test3a = 
  c("ENSMMUG00000031200",
    "ENSMMUG00000017368",
    "ENSMMUG00000018891",
    "ENSMMUG00000021266",
    "ENSMMUG00000023311",
    "ENSMMUG00000031191",
    "ENSMMUG00000000021",
    "ENSMMUG00000009226",
    "ENSMMUG00000031107",
    "ENSMMUG00000013779",
    "ENSMMUG00000007266",
    "ENSMMUG00000026913",
    "ENSMMUG00000027103",
    "ENSMMUG00000026930")

test3b = 
  c("ENSMMUG00000031200",
    "ENSMMUG00000017368",
    "ENSMMUG00000018891",
    "ENSMMUG00000021266",
    "ENSMMUG00000023311",
    "ENSMMUG00000031191",
    "CATSPERE",
    "MAP3K21",
    "ENSMMUG00000031107",
    "ENSMMUG00000013779",
    "EEF1AKNMT",
    "mml-mir-214",
    "mml-mir-199a-1",
    "mml-mir-488")

expect_test1a_ensembl = 
  c("ADSS2",
    "ZBTB18",
    "BRP44",
    "XCL1",
    "CCDC181",
    "TNFSF18",
    "SNORD24",
    "COP1",
    "BRINP2")

expect_test2a_string = 
  c("C1orf129",
    "RPL7A",
    "RFWD2")

expect_test3a_david = 
  c("LOC100427314",
    "LOC694380",
    "LOC100427785",
    "LOC701769",
    "ZNF669",
    "ZNF695",
    "C1H1orf101",
    "KIAA1804",
    "LOC695277",
    "LOC100423131",
    "METTL13",
    "MIR214",
    "MIR199A-1",
    "MIR488")





test_that("Gene aliasing code works as expected", {
  testEnsembl <- aliasENSEMBL(ensemblIds = test1a,
                              geneSymbols = NULL,
                              attributes = c('hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'),
                              replaceUnmatched = F)

  expect_equal(as.character(testEnsembl$Coalesced.ENSEMBL), expect_test1a_ensembl)

  testString <- aliasSTRINGdb(ensemblIds = NULL,
                              geneSymbols = test2a,
                              speciesId = 9606)
  expect_equal(as.character(testString$Coalesced.STRING), expect_test2a_string)

  davidEmail <- Sys.getenv('DAVID_EMAIL')
  if (is.na(davidEmail) || davidEmail == '') {
    stop('DAVID_EMAIL environment variable must be set!')
  }

  testDavid <- aliasDAVID(ensemblIds = test3a,
                          geneSymbols = NULL,
                          email = davidEmail)
  expect_equal(as.character(testDavid$Coalesced.DAVID), expect_test3a_david)
  expect_length(as.character(testDavid$Coalesced.DAVID), length(expect_test3a_david))

  testAliasTable <- aliasTable(geneSymbols = test2b, ensemblIds = test2a, davidEmail = davidEmail)
  expect_length(as.character(testAliasTable$Consensus), length(expect_test2b))

})