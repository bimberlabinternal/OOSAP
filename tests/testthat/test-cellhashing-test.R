skip_on_cran()

context("CellHashing")

DoTest <- function(barcodeFile, summaryFile, whitelistFile, doRowFilter = F) {

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

    dt <- ProcessEnsemblHtoCalls(mc, sc, barcodeData, outFile = summaryFile)

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
        Singlet = 8155,
        MultiSeq = 1206,
        Seurat = 9529,
        DoRowFilter = T
    ),
    '283' = list(
        input = '../testdata/cellHashing/283-cellbarcodeToHTO.calls.citeSeqCounts.txt', htos = c(2:6),
        gexBarcodeFile = '../testdata/cellHashing/283-validCellIndexes.csv',
        CalledCells = 5346,
        Singlet = 2697,
        MultiSeq = 853,
        Seurat = 2900,
        DoRowFilter = T
    )
    # TODO: revisit the example data
    # 'NewFormat' = list(
    #     input = '../testdata/cellHashing/umi_count',
    #     htos = c(1:4),
    #     gexBarcodeFile = NULL,
    #     CalledCells = 5346,
    #     Singlet = 2697,
    #     MultiSeq = 853,
    #     Seurat = 2900,
    #     DoRowFilter = F
    # )
)

for (testName in names(tests)) {
    print(paste0('Running test: ', testName))
    test <- tests[[testName]]

    summaryFile <- NULL
    if (!is.null(test$gexBarcodeFile)) {
        summaryFile <- 'summary.txt'
    }

    l <- DoTest(test$input, summaryFile, test$gexBarcodeFile, doRowFilter = test$DoRowFilter)
    barcodeData <- l$barcodeData

    expectedHtos <- sort(paste0('HTO-', test$htos))
    actualHtosMatrix <- sort(unname(simplifyHtoNames(rownames(barcodeData))))

    test_that("HTOs not equal", {
        expect_equal(expectedHtos, actualHtosMatrix)
    })

    dt <- l$dt
    print(nrow(dt))
    print(nrow(dt[dt$HTO_Classification == 'Singlet',]))
    print(sum(dt$Seurat))
    print(sum(dt$MultiSeq))

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

    if (!is.null(summaryFile)) {
        unlink(summaryFile)
    }
}

