skip_on_cran()

context("CellHashing")

DoTest <- function(barcodeFile, summaryFile, whitelistFile) {

    barcodeData <- ProcessCiteSeqCount(bFile=barcodeFile)

    if (nrow(barcodeData) == 0) {
        stop('No passing HTOs')
    }

    if (ncol(barcodeData) == 0) {
        stop('No passing cells')
    }

    generateQcPlots(barcodeData)

    sc <- generateCellHashCallsSeurat(barcodeData)

    mc <- generateCellHashCallsMultiSeq(barcodeData)

    dt <- processEnsemblHtoCalls(mc, sc, barcodeData, outFile = summaryFile)

    if (!is.null(whitelistFile) && !is.null(summaryFile)){
        generateSummaryForExpectedBarcodes(dt, whitelistFile = whitelistFile, outputFile = summaryFile, barcodeData = barcodeData)
    }

    return(list(barcodeData = barcodeData, dt = dt))
}


tests <- list(
    '282-1' = list(input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt', htos = c(1:3, 8, 10, 12), gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt'),
    '283' = list(input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt', htos = c(2:6), gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv'),
    'NewFormat' = list(input = '../testdata/cellHashing/umi_count', htos = c(1:4), gexBarcodeFile =NULL)
)

for (testName in names(tests)) {
    print(paste0('Running test: ', testName))
    test <- tests[[testName]]

    summaryFile <- NULL
    if (!is.null(test$gexBarcodeFile)) {
        summaryFile <- 'summary.txt'
    }

    l <- DoTest(test$input, summaryFile, test$gexBarcodeFile)
    barcodeData <- l$barcodeData

    expectedHtos <- sort(paste0('HTO-', test$htos))
    actualHtosMatrix <- sort(unname(simplifyHtoNames(rownames(barcodeData))))

    test_that("HTOs not equal", {
        expect_equal(expectedHtos, actualHtosMatrix)
    })

    dt <- l$dt
    print(str(dt))
    print(nrow(dt))

    if (!is.null(test$gexBarcodeFile)) {

    }
}

