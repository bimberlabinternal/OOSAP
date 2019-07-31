context("CellHashing")

DoTest <- function(barcodeFile, callsFile, summaryFile, whitelistFile, doRowFilter = F) {

    barcodeData <- ProcessCiteSeqCount(bFile=barcodeFile, doRowFilter = doRowFilter)

    if (nrow(barcodeData) == 0) {
        stop('No passing HTOs')
    }

    if (ncol(barcodeData) == 0) {
        stop('No passing cells')
    }

    GenerateQcPlots(barcodeData)

    sc <- GenerateCellHashCallsSeurat(barcodeData)

    mc <- GenerateCellHashCallsMultiSeq(barcodeData)

    dt <- ProcessEnsemblHtoCalls(mc, sc, barcodeData, outFile = callsFile)

    if (!is.null(whitelistFile) && !is.null(summaryFile)){
      GenerateSummaryForExpectedBarcodes(dt, whitelistFile = whitelistFile, outputFile = summaryFile, barcodeData = barcodeData)
    }

    return(list(barcodeData = barcodeData, dt = dt))
}


tests <- list(
    '282-1' = list(
        input = '../testdata/cellHashing/282-1-HTO_cellHashingRawCounts.txt',
        htos = c(2:3, 8, 10, 12),
        gexBarcodeFile = '../testdata/cellHashing/282-1-whitelist.txt',
        CalledCells = 21263,
        Singlet = 12176,
        MultiSeq = 1206,
        Seurat = 13550,
        CallRows = 21263,
        DoRowFilter = T
    ),
    '283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt', htos = c(2:6),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 5346,
        Singlet = 2697,
        MultiSeq = 853,
        Seurat = 2900,
        CallRows = 5346,
        DoRowFilter = T
    )
    # 'NewFormat' = list(
    #     input = '../testdata/cellHashing/umi_count',
    #     htos = c(1:4),
    #     gexBarcodeFile = NULL,
    #     CalledCells = 5346,
    #     Singlet = 2697,
    #     MultiSeq = 853,
    #     Seurat = 2900,
    #     CallRows = 100,
    #     DoRowFilter = F
    # )
)

for (testName in names(tests)) {
    print(paste0('Running test: ', testName))
    test <- tests[[testName]]

    callsFile <- paste0(testName, '-calls.txt')
    summaryFile <- NULL
    if (!is.null(test$gexBarcodeFile)) {
        summaryFile <- paste0(testName, '-summary.txt')
    }

    l <- DoTest(test$input, callsFile=callsFile, summaryFile=summaryFile, whitelistFile=test$gexBarcodeFile, doRowFilter = test$DoRowFilter)
    barcodeData <- l$barcodeData

    expectedHtos <- sort(paste0('HTO-', test$htos))
    actualHtosMatrix <- sort(unname(simplifyHtoNames(rownames(barcodeData))))

    test_that("HTOs not equal", {
        expect_equal(expectedHtos, actualHtosMatrix)
    })

    dt <- l$dt

    test_that("Called cells not equal", {
        expect_equal(test[['CalledCells']], nrow(dt))
    })

    test_that("Called Singlets not equal", {
        expect_equal(test[['Singlet']], nrow(dt[dt$HTO_Classification == 'Singlet',]))
    })

    test_that("Seurat called not equal", {
        expect_equal(test[['Seurat']], sum(dt$Seurat))
    })

    test_that("MultiSeq called not equal", {
        expect_equal(test[['MultiSeq']], sum(dt$MultiSeq))
    })

    d <- read.table(callsFile, header = T, sep = '\t')
    test_that("Call file rowcount not correct", {
        expect_equal(test[['CallRows']], nrow(d))
    })
    unlink(callsFile)

    if (!is.null(summaryFile)) {
        d <- read.table(summaryFile, header = T, sep = '\t')
        test_that("Metric rowcount not correct", {
            expect_equal(19, nrow(d))
        })
       unlink(summaryFile)
    }
}

