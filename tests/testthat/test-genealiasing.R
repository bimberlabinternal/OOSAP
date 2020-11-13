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
   "MPC2",
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
    "MPC2",
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
	biomaRt::biomartCacheClear()

  #Ensembl
  testEnsemblWithId <- TranslateToEnsembl(ensemblIds = test1a, geneSymbols = NULL)
  expect_equal(as.character(testEnsemblWithId$EnsemblId), test1a)  #Ensure order preserved
  expect_equal(as.character(testEnsemblWithId$external_gene_name), c('ADSS2','ZBTB18','MPC2','XCL1','CCDC181','TNFSF18','SNORD24','COP1','BRINP2'))
  
  testEnsemblWithSymbol <- TranslateToEnsembl(ensemblIds = NA, geneSymbols = expect_test1a_ensembl)
  expect_equal(as.character(testEnsemblWithSymbol$GeneSymbol), expect_test1a_ensembl)
  expect_true(all(is.na(testEnsemblWithSymbol$EnsemblId)))
  expect_equal(as.character(testEnsemblWithSymbol$external_gene_name), c('ADSS2','ZBTB18','MPC2','XCL1','CCDC181','TNFSF18', NA,'COP1','BRINP2'))
  expect_equal(as.character(testEnsemblWithSymbol$ensembl_gene_id), c('ENSMMUG00000000016','ENSMMUG00000015664','ENSMMUG00000010387','ENSMMUG00000060218','ENSMMUG00000000065','ENSMMUG00000059499',NA,'ENSMMUG00000059669','ENSMMUG00000022935'))

  testEnsemblWithBoth <- TranslateToEnsembl(ensemblIds = test1a, geneSymbols = expect_test1a_ensembl)
  expect_equal(as.character(testEnsemblWithBoth$GeneSymbol), expect_test1a_ensembl)
  expect_equal(as.character(testEnsemblWithBoth$EnsemblId), test1a)
  expect_equal(as.character(testEnsemblWithBoth$external_gene_name), c('ADSS2','ZBTB18','MPC2','XCL1','CCDC181','TNFSF18', 'SNORD24','COP1','BRINP2'))
  expect_equal(as.character(testEnsemblWithBoth$ensembl_gene_id), c('ENSMMUG00000000016','ENSMMUG00000015664','ENSMMUG00000010387','ENSMMUG00000060218','ENSMMUG00000000065','ENSMMUG00000059499','ENSMMUG00000050963','ENSMMUG00000059669','ENSMMUG00000022935'))
  
  #STRINGdb
  testString <- TranslateToStringDb(ensemblIds = NULL, geneSymbols = expect_test1a_ensembl)
  expect_equal(as.character(testString$GeneSymbol), expect_test1a_ensembl)
  expect_true(all(is.na(testString$EnsemblId)))
  expect_equal(testString$STRING_id, c('9606.ENSP00000355493','9606.ENSP00000351539','9606.ENSP00000356820','9606.ENSP00000356792','9606.ENSP00000356780','9606.ENSP00000385470', NA,'9606.ENSP00000356641,9606.ENSP00000364858','9606.ENSP00000354481'))
  
  #DAVID
  davidEmail <- Sys.getenv('DAVID_EMAIL')
  if (is.na(davidEmail) || davidEmail == '') {
    stop('DAVID_EMAIL environment variable must be set!')
  }

  testDavidId <- TranslateToDAVID(ensemblIds = test3a, geneSymbols = NULL, email = davidEmail)
  expect_equal(as.character(testDavidId$EnsemblId), test3a)
  expect_equal(as.character(testDavidId$DAVID.Symbol), c('LOC100427314','LOC694380','LOC100427785','LOC701769','ZNF669','ZNF695','C1H1orf101','KIAA1804','LOC695277','LOC100423131','METTL13','MIR214','MIR199A-1','MIR488'))
  
  #Not currently possible, so returns nothing
  testDavidSymbol <- TranslateToDAVID(ensemblIds = NA, geneSymbols = test3b, email = davidEmail)
  expect_equal(as.character(testDavidSymbol$GeneSymbol), test3b)
  expect_true(all(is.na(testDavidSymbol$DAVID.Id)))
  
  testDavidBoth <- TranslateToDAVID(ensemblIds = test3a, geneSymbols = test3b, email = davidEmail)
  expect_equal(as.character(testDavidBoth$GeneSymbol), test3b)
  expect_equal(as.character(testDavidBoth$DAVID.Symbol), c('LOC100427314','LOC694380','LOC100427785','LOC701769','ZNF669','ZNF695','C1H1orf101','KIAA1804','LOC695277','LOC100423131','METTL13','MIR214','MIR199A-1','MIR488'))
  
  #Now all combined:
  testAliasTable <- TranslateGeneNames(ensemblIds = test3a, geneSymbols = test3b, davidEmail = davidEmail)
  expect_equal(as.character(testAliasTable$EnsemblId), test3a)
  expect_equal(as.character(testAliasTable$GeneSymbol), test3b)
  expect_equal(as.character(testAliasTable$DAVID.Symbol), c('LOC100427314','LOC694380','LOC100427785','LOC701769','ZNF669','ZNF695','C1H1orf101','KIAA1804','LOC695277','LOC100423131','METTL13','MIR214','MIR199A-1','MIR488'))
  #expect_equal(as.character(testAliasTable$ensembl_gene_id), c('ENSMMUG00000000016','ENSMMUG00000015664','ENSMMUG00000010387','ENSMMUG00000060218','ENSMMUG00000000065','ENSMMUG00000059499','ENSMMUG00000050963','ENSMMUG00000059669','ENSMMUG00000022935'))

	biomaRt::biomartCacheClear()
})