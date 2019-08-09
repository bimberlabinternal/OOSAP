#' @include LabKeySettings.R
#' @import Rlabkey

#' @title GetCdnaRecords
#'
#' @export
GetCdnaRecords <- function(readsetId, type = 'GEX') {
  fieldName <- switch(type,
                      'GEX' = 'readsetId',
                      'TCR' = 'enrichedReadsetId',
                      'HTO' = 'hashingReadsetId'
  )

  df <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="tcrdb",
    queryName="clones",
    showHidden=TRUE,
    #colSelect=c(''),
    colFilter=makeFilter(c(fieldName, "EQUAL", readsetId)),
    containerFilter=NULL,
    colNameOpt='rname'
  )

  return(df)
}

utils::globalVariables(
  names = c('SortOrder'),
  package = 'OOSAP',
  add = TRUE
)
#' @title QueryAndApplyCdnaMetadata
#'
#' @param seuratObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr %>% group_by_at summarise_at arrange
QueryAndApplyCdnaMetadata <- function(seuratObj,
                                      fieldSelect = c('rowid', 'sortid/population', 'sortid/stimid/animalId', 'sortid/stimid/date', 'sortid/stimid/stim'),
                                      fieldNames = c('cDNA ID', 'Population', 'SubjectId', 'SampleDate', 'Stim'), overwriteExisting = F) {
  if (length(fieldSelect) != length(fieldNames)) {
    stop('The length of fields must equal the length of fieldNames')
  }

  if (sum(duplicated(fieldSelect)) > 0) {
    stop('All values for fieldSelect must be unique')
  }

  #spike in readset, since we need this to join dataframes
  fieldSelect <- tolower(fieldSelect)
  readsetAdded <- F
  if (!('readsetid' %in% fieldSelect)) {
    readsetAdded <- T
    fieldSelect <- c('readsetid', fieldSelect)
    fieldNames <- c('Readset', fieldNames)
  }
  readsetIdx <- which('readsetid' %in% fieldSelect)
  readsetLabel <- fieldNames[readsetIdx]

  #also HTO:
  htoAdded <- F
  if (!('sortid/hto' %in% fieldSelect)) {
    htoAdded <- T
    fieldSelect <- c('sortid/hto', fieldSelect)
    fieldNames <- c('HTO', fieldNames)
  }
  htoIdx <- which('sortid/hto' %in% fieldSelect)
  htoLabel <- fieldNames[htoIdx]

  #Download info, based on BarcodePrefix:
  outputFiles <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colFilter=makeFilter(c('rowId', "IN", paste0(unique(seuratObj$BarcodePrefix), collapse = ';'))),
    colSelect="rowid,readset",
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(outputFiles) <- c('BarcodePrefix', readsetLabel)
  outputFiles$BarcodePrefix <- as.character(outputFiles$BarcodePrefix)
  print(paste0('total outputfile rows: ', nrow(outputFiles)))
  if (nrow(outputFiles) != length(unique(seuratObj$BarcodePrefix))) {
    missing <- sort(c(setdiff(unique(seuratObj$BarcodePrefix), unique(outputFiles$BarcodePrefix)), setdiff(unique(outputFiles$BarcodePrefix), unique(seuratObj$BarcodePrefix))))
    warning(paste0('Did not find output file record for all prefixes!  Missing: ', paste0(missing, collapse = ',')))
  }

  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="tcrdb",
    queryName="cdnas",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("readsetId", "IN", paste0(unique(outputFiles$Readset), collapse = ";"))),
    colSelect=paste0(fieldSelect, collapse = ','),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(rows) <- fieldNames
  print(paste0('total cDNA rows: ', nrow(rows)))

  rows <- merge(rows, outputFiles, by = c(readsetLabel), all.x = T)
  rows <- unique(rows)

  #Force unique values
  groupVars <- c(htoLabel, 'BarcodePrefix')
  colSummary <- fieldNames[!(fieldNames %in% groupVars)]
  rows <- rows %>% group_by_at(groupVars) %>% summarise_at(colSummary, function(x){
    paste(sort(unique(x)), collapse=',')
  })

  origBarcodes <- colnames(seuratObj)
  df <- data.frame(HTO = as.character(seuratObj$HTO), BarcodePrefix = as.character(seuratObj$BarcodePrefix), Barcode = origBarcodes, SortOrder = 1:length(origBarcodes))
  names(df) <- c(htoLabel, 'BarcodePrefix', 'Barcode', 'SortOrder')
  df <- merge(df, rows, by = c(htoLabel, 'BarcodePrefix'), all.x = T)
  df <- dplyr::arrange(df, SortOrder)
  df <- df[colnames(df) != 'SortOrder']
  rownames(df) <- df$Barcode
  df <- df[colnames(df) != 'Barcode']

  if (nrow(df) != ncol(seuratObj)) {
    stop('Length of original seurat object and metadata not equal. Something went wrong merging')
  }

  if (sum(origBarcodes != rownames(df)) > 0) {
    stop('BarcodePrefix does not match for all cells.  Something went wrong merging')
  }

  #strings to factors:
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor)

  #now apply the result:
  for (idx in 1:length(fieldSelect)) {
    fieldName <- fieldNames[idx]
    fieldKey <- fieldSelect[idx]
    if (readsetAdded && fieldName == readsetLabel) {
      next
    }

    if (htoAdded && fieldName == htoLabel) {
      next
    }

    if (!overwriteExisting && (fieldName %in% names(seuratObj@meta.data))) {
      print(paste0('column already exists, skipping: ', fieldName))
      next
    }

    print(paste0('Adding column: ', fieldName, ' (', fieldKey, ')'))
    v <- df[[fieldName]]
    names(v) <- names(df)
    seuratObj[[fieldName]] <- v
  }

  return(seuratObj)
}

GenerateDataToCompareBarcodeSets <- function(workbooks, savePath = './') {
  metrics <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="quality_metrics",
    viewName="",
    colSort="-rowid",
    colFilter=makeFilter(
      c("category", "IN", "Cell Hashing;Cell Hashing Concordance"),
      c("workbook/workbookId", "IN", paste0(workbooks, collapse = ';'))
    ),
    containerFilter=NULL,
    colNameOpt="rname"
  )

  summary <- data.frame(
    Name = character(),
    GEX_ReadsetId = integer(),
    HTO_Reads = integer(),
    GEX_Reads = integer(),
    GEX_FractionOfInputCalled = numeric(),
    GEX_InputBarcodes = numeric(),
    GEX_TotalCalledNotInInput = numeric(),
    TCR_FractionOfInputCalled = numeric(),
    TCR_InputBarcodes = numeric(),
    TCR_TotalCalledNotInInput = numeric(),
    GEX_WhitelistFile = character(),
    TCR_WhitelistFile = character(),
    HTO_Top_BarcodesFile = character(),
    GEX_CallsFile = character(),
    TCR_CallsFile = character(),
    ExpectedHTOs = character()
  )

  for (wb in workbooks){
    print(paste0('processing: ', wb))
    localPath <- file.path(savePath, wb)
    if (!dir.exists(localPath)){
      dir.create(localPath)
    }

    cDNAs <- labkey.selectRows(
      baseUrl=lkBaseUrl,
      folderPath=paste0(lkDefaultFolder, wb),
      schemaName="tcrdb",
      queryName="cdnas",
      viewName="",
      colSort="-rowid",
      colSelect = 'rowid,readsetid,readsetid/name,hashingreadsetid,enrichedreadsetid,hashingreadsetid/totalforwardReads,readsetid/totalforwardReads,sortId/hto',
      containerFilter=NULL,
      colNameOpt="rname"
    )
    print(paste0('total cDNA records: ', nrow(cDNAs)))

    callFiles <- labkey.selectRows(
      baseUrl=lkBaseUrl,
      folderPath=paste0(lkDefaultFolder, '/', wb),
      schemaName="sequenceanalysis",
      queryName="outputfiles",
      viewName="",
      colSort="-rowid",
      colSelect="rowid,name,description,readset,readset/name,category,dataid/RowId,workbook,dataid/WebDavUrlRelative,dataid/WebDavUrlRelative,created",
      colFilter=makeFilter(c("category", "IN", "Cell Hashing Calls (VDJ);Cell Hashing Calls;10x GEX Cell Hashing Calls")),
      containerFilter=NULL,
      colNameOpt="rname"
    )

    htoSummary <- cDNAs %>% group_by(readsetid) %>% summarise(ExpectedHTOs = paste0(sort(unique(sortid_hto)), collapse = ","))
    uniqueRs <- c()
    for (i in 1:nrow(cDNAs)) {
      row <- cDNAs[i,]
      if (row$readsetid %in% uniqueRs) {
        next
      }

      uniqueRs <- c(uniqueRs, row$readsetid)

      n <- gsub(x = row[['readsetid_name']], pattern = '-GEX', replacement = '')
      toAdd <- data.frame(Name = n, GEX_ReadsetId = row[['readsetid']], HTO_Reads = row[['hashingreadsetid_totalforwardreads']], GEX_Reads = row[['readsetid_totalforwardreads']])

      #metrics:
      htoMetrics <- metrics[metrics$readset == row[['hashingreadsetid']],]
      if (nrow(htoMetrics) > 0) {
        m <- htoMetrics[htoMetrics$metricname == 'Singlet',]
        if (nrow(m) > 0) {
          m <- m[m$dataid == max(m$dataid),]

          toAdd$HTO_Singlet <- m$metricvalue
        }
      } else {
        toAdd$HTO_Singlet <- NA
      }
      toAdd <- AppendMetrics(row, toAdd, metrics, 'GEX', 'readsetid')
      toAdd <- AppendMetrics(row, toAdd, metrics, 'TCR', 'enrichedreadsetid')

      #call files:
      toAdd <- DownloadCallFile(wb, callFiles, row[['readsetid']], toAdd, 'GEX_WhitelistFile', savePath, T)
      toAdd <- DownloadCallFile(wb, callFiles, row[['enrichedreadsetid']], toAdd, 'TCR_WhitelistFile', savePath, T)
      toAdd <- DownloadCallFile(wb, callFiles, row[['hashingreadsetid']], toAdd, 'HTO_Top_BarcodesFile', savePath, F)

      toAdd <- DownloadCallFile(wb, callFiles, row[['readsetid']], toAdd, 'GEX_CallsFile', savePath, F)
      toAdd <- DownloadCallFile(wb, callFiles, row[['enrichedreadsetid']], toAdd, 'TCR_CallsFile', savePath, F)

      #now merge HTOs:
      toAdd <- merge(toAdd, htoSummary, by.x = c('GEX_ReadsetId'), by.y = c('readsetid'), all.x = T)

      summary <- rbind(summary, toAdd)
    }
  }

  return(summary)
}

AppendMetrics <- function(row, toAdd, metrics, type, field) {
  for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
    toAdd[[paste0(type, '_', name)]] <- NA
  }

  xMetrics <- metrics[metrics$readset == row[[field]],]
  if (nrow(xMetrics) > 0) {
    latest <- xMetrics[xMetrics$dataid == max(xMetrics$dataid),]

    for (name in c('FractionOfInputCalled', 'InputBarcodes', 'TotalCalledNotInInput')) {
      toAdd[[paste0(type, '_', name)]] <- latest[latest$metricname == name,]$metricvalue
    }
  }

  return(toAdd)
}

DownloadCallFile <- function(wb, callFiles, readsetId, toAdd, fieldName, localPath, downloadWhitelist) {
  toAdd[fieldName] <- NA

  cf <- callFiles[callFiles$readset == readsetId,]
  if (nrow(cf) > 1) {
    print(paste0('Multiple call files, using latest: ', toAdd$Name, ' ', fieldName))
  }

  #use the most recent (highest rowId)
  if (nrow(cf) > 0) {
    row <- cf[cf$rowid == max(cf$rowid),]

    suffix <- '.calls'
    if (downloadWhitelist) {
      suffix <- '.validCellIndexes'
    }

    fn <- paste0(row['readset_name'], '.', row['readset'], '.', row['rowid'], '.', row$category, suffix, '.txt')
    fn <- gsub('\\(', '_', fn)
    fn <- gsub('\\)', '_', fn)
    fn <- gsub(' ', '_', fn)
    fn <- gsub('__', '_', fn)

    remotePath <- row[['dataid_webdavurlrelative']]
    if (downloadWhitelist) {
      remotePath <- paste0(dirname(remotePath), '/validCellIndexes.csv')
    }

    lp <- paste0(localPath, '/', wb, '/', fn)
    if (file.exists(lp)) {
      print(paste0('file exists, reusing: ', lp))
    } else {
      success <- labkey.webdav.get(
        baseUrl=lkBaseUrl,
        folderPath=paste0(lkDefaultFolder, wb),
        remoteFilePath = remotePath,
        overwrite = T,
        localFilePath = lp
      )

      if (!success) {
        warning(paste0('Unable to download whitelist for readset: ', row['readset_name'], ' from: ', remotePath))
      }
    }

    toAdd[fieldName] <- lp
  }

  return(toAdd)
}

ProcessSet <- function(df, set1, set2, type1, type2) {
  for (name1 in names(set1)) {
    h <- set1[[name1]]

    for (name2 in names(set2)) {
      g <- set2[[name2]]
      i <- length(intersect(h, g))

      df <- rbind(df, data.frame(dataset1 = name1, type1 = type1, dataset2 = name2, type2 = type2, intersect = i, length1 = length(h), length2 = length(g)))
    }
  }

  return(df)
}

utils::globalVariables(
  names = c('dataset1', 'type1', 'type2', 'max_intersect', 'max_intersect_by_type'),
  package = 'OOSAP',
  add = TRUE
)
#' @title CompareCellBarcodeSets
#'
#' @description This iterates the provided workbooks, identifying every 10x cDNA library and the GEX/TCR/HTO libraries.  It generates summaries of the cell barcode intersect between them, which can help debug sample swaps.
#' @param workbooks A vector of workbook IDs
#' @export
#' @importFrom dplyr %>% mutate group_by
CompareCellBarcodeSets <- function(workbooks, savePath = './') {
  summary <- GenerateDataToCompareBarcodeSets(workbooks, savePath)

  htoBC <- list()
  gexBC <- list()
  tcrBC <- list()

  for (i in 1:nrow(summary)) {
    row <- summary[i,]
    name <- as.character(row$Name)

    if (!is.na(row$HTO_Top_BarcodesFile) & !is.na(row$GEX_WhitelistFile) & !is.na(row$TCR_WhitelistFile)){
      bc <- read.table(row$HTO_Top_BarcodesFile, header = T, sep = '\t')
      htoBC[[name]] <- bc$CellBarcode

      bc <- read.table(row$GEX_WhitelistFile, header = F, sep = '\t')$V1
      gexBC[[name]] <- bc

      bc <- read.table(row$TCR_WhitelistFile, header = F, sep = '\t')$V1
      tcrBC[[name]] <- bc

    } else {
      print(paste0('missing one ore more files: ', name, ':'))
      if (is.na(row$HTO_Top_BarcodesFile)){
        print('HTO missing')
      }

      if (is.na(row$GEX_WhitelistFile)){
        print('GEX missing')
      }

      if (is.na(row$TCR_WhitelistFile)){
        print('TCR missing')
      }
    }
  }

  df <- data.frame(dataset1 = character(), type1 = character(), dataset2 = character(), type2 = character(), intersect = integer(), length1 = integer(), length2 = integer())

  df <- ProcessSet(df, htoBC, gexBC, 'HTO', 'GEX')
  df <- ProcessSet(df, htoBC, tcrBC, 'HTO', 'TCR')

  df <- ProcessSet(df, gexBC, htoBC, 'GEX', 'HTO')
  df <- ProcessSet(df, gexBC, tcrBC, 'GEX', 'TCR')

  df <- ProcessSet(df, tcrBC, htoBC, 'TCR', 'HTO')
  df <- ProcessSet(df, tcrBC, gexBC, 'TCR', 'GEX')

  df$fraction <- df$intersect / df$length1
  df <- df[order(df$dataset1, df$type1, df$type2, -df$intersect),]
  df <- df %>% group_by(dataset1, type1, type2) %>% mutate(max_intersect_by_type = max(intersect))
  df <- df %>% group_by(dataset1, type1) %>% mutate(max_intersect = max(intersect))

  write.table(df, file = file.path(savePath, 'cell_barcode_comparisons.txt'), sep = '\t', row.names = F, quote = F)

  self <- df[df$dataset1 == df$dataset2, c('dataset1', 'type1', 'type2', 'intersect', 'fraction')]
  names(self) <- c('dataset1', 'type1', 'type2', 'self_intersect', 'self_intersect_fraction')

  df2 <- df[df$intersect == df$max_intersect_by_type,]
  write.table(df2, file = 'top_intersect_by_type.txt', sep = '\t', row.names = F, quote = F)

  df3 <- df2[df2$dataset1 != df2$dataset2,]
  df3 <- merge(df3, self, by = c('dataset1', 'type1', 'type2'), all.x = T, all.y = F)
  #filter ties:
  df3 <- df3[df3$intersect != df3$self_intersect,]
  write.table(df3, file = file.path(savePath, 'conflicting_intersect.txt'), sep = '\t', row.names = F, quote = F)

  #Now look for instances where the raw HTO calls find more HTOs than either TCR or GEX
  dfSummary <- NA
  for (i in 1:nrow(summary)) {
    row <- summary[i,]
    name <- as.character(row$Name)

    if (is.na(row$HTO_Top_BarcodesFile)){
      print(paste0('No HTO call file for: ', name))
    }

    df <- data.frame(HTO = character(), Type = character(), Count.HTO = integer(), Fraction.HTO = numeric(), Count.Compare = integer(), Fraction.Compare = numeric(), Difference = numeric())
    #GEX:
    dfG <- CompareHtosByCall(name, row$HTO_Top_BarcodesFile, row$GEX_CallsFile, 'GEX')
    if (!all(is.na(dfG))){
      df <- rbind(df, dfG)
    }

    #TCR
    dfT <- CompareHtosByCall(name, row$HTO_Top_BarcodesFile, row$TCR_CallsFile, 'TCR')
    if (!all(is.na(dfT))){
      df <- rbind(df, dfT)
    }

    if (all(is.na(dfSummary))) {
      dfSummary <- df
    } else {
      dfSummary <- rbind(dfSummary, df)
    }
  }

  dfSummary <- merge(dfSummary, summary[c('Name', 'ExpectedHTOs')], all.x = T, by.x = c('Dataset'), by.y = c('Name'))
  dfSummary$Unexpected <- apply(dfSummary, 1, function(r){
    htos <- unlist(strsplit(r['ExpectedHTOs'], ','))

    return(!(r['HTO'] %in% htos))
  })


  write.table(dfSummary, file = file.path(savePath, paste0('htoCompare.txt')), sep = '\t', row.names = F, quote = F)
  write.table(dfSummary[dfSummary$Unexpected & dfSummary$Fraction.HTO > 0.025,], file = file.path(savePath, paste0('UnexpectedHTOs.txt')), sep = '\t', row.names = F, quote = F)

  return(dfSummary)
}

CompareHtosByCall <- function(name, htoCallFile, compareFile, type) {
  if (is.na(compareFile)){
    print(paste0('Missing ', type, ' call file for: ', name))
    return(NA)
  }

  t1 <- read.table(htoCallFile, header = T, sep = '\t')
  t1 <- t1 %>% group_by(HTO) %>% summarise(Count = n())
  t1 <- t1[!(t1$HTO %in% c('Doublet', 'Negative')),]
  t1$Fraction <- t1$Count / sum(t1$Count)

  t2 <- read.table(compareFile, header = T, sep = '\t')
  t2 <- t2 %>% group_by(HTO) %>% summarise(Count = n())
  t2 <- t2[!(t2$HTO %in% c('Doublet', 'Negative')),]
  t2$Fraction <- t2$Count / sum(t2$Count)

  df <- merge(x = t1, y = t2, by = 'HTO', all = T, suffixes = c('.HTO', '.Compare'))
  df$Count.HTO[is.na(df$Count.HTO)] <- 0
  df$Fraction.HTO[is.na(df$Fraction.HTO)] <- 0

  df$Count.Compare[is.na(df$Count.Compare)] <- 0
  df$Fraction.Compare[is.na(df$Fraction.Compare)] <- 0

  df$Type <- type
  df$Difference <- abs(df$Fraction.HTO - df$Fraction.Compare)
  df$Dataset <- name

  df <- df[c('Dataset', 'HTO', 'Type', 'Count.HTO', 'Fraction.HTO', 'Count.Compare', 'Fraction.Compare', 'Difference')]

  return(df)
}

