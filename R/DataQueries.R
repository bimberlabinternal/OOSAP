#' @import Rlabkey

Rlabkey::labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")

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
    folderPath="/Labs/Bimber/",
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

#' @title QueryAndApplyCdnaMetadata
#'
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @export
#' @importFrom dplyr %>% group_by_at summarise_at
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
    folderPath="/Labs/Bimber/",
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

  rows <- labkey.selectRows(
    folderPath="/Labs/Bimber/",
    schemaName="tcrdb",
    queryName="cdnas",
    viewName="",
    colSort="-rowid",
    colFilter = makeFilter(c("readsetId", "IN", paste0(unique(outputFiles$Readset, collapse = ";")))),
    colSelect=paste0(fieldSelect, collapse = ','),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  names(rows) <- fieldNames

  rows <- merge(rows, outputFiles, by = c(readsetLabel), all.x = T)
  rows <- unique(rows)

  #Force unique values
  groupVars <- c(htoLabel, 'BarcodePrefix')
  colSummary <- fieldNames[!(fieldNames %in% groupVars)]
  rows <- rows %>% group_by_at(groupVars) %>% summarise_at(colSummary, function(x){
    paste(sort(unique(x)), collapse=',')
  })

  df <- data.frame(HTO = as.character(seuratObj$HTO), BarcodePrefix = as.character(seuratObj$BarcodePrefix))
  names(df) <- c(htoLabel, 'BarcodePrefix')
  df <- merge(df, rows, by = c(htoLabel, 'BarcodePrefix'), all.x = T)
  if (nrow(df) != ncol(seuratObj)) {
    stop('Length of original seurat object and metadata not equal. Something went wrong merging')
  }

  if (sum(seuratObj$BarcodePrefix != df$BarcodePrefix) > 0) {
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
    seuratObj[[fieldName]] <- df[[fieldName]]
  }

  return(seuratObj)
}
