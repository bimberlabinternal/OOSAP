context("scRNAseq")

test_that("Cite-Seq works", {
  #Reduce size, touch up data:
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
	seuratObj <- seuratObj[1:1000,1:100]

	#this is technically cell hashing, but shouldnt matter:
	countDir <- '../testdata/cellHashing/umi_count'
	citeseqData <- Seurat::Read10X(data.dir = countDir, gene.column=1)
	citeseqData1 <- citeseqData[,1:40]
	colnames(citeseqData1) <- colnames(seuratObj)[1:40]
	inputPath1 <- './input1'
	if (dir.exists(inputPath1)) {
	  unlink(inputPath1, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath1, citeseqData1)
	
	citeseqData2 <- citeseqData[,60:100]
	colnames(citeseqData2) <- colnames(seuratObj)[60:100]
	inputPath2 <- './input2'
	if (dir.exists(inputPath2)) {
	  unlink(inputPath2, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath2, citeseqData2)
	
	#First w/o prefix:
	seuratObjCite <- OOSAP:::AppendCiteSeq(seuratObj = seuratObj, countMatrixDir = inputPath1, barcodePrefix = NULL)
	expect_equal(colnames(seuratObjCite@assays$ADT), colnames(seuratObj@assays$RNA))
	expect_equal(rownames(seuratObjCite@assays$ADT), rownames(citeseqData1)[rownames(citeseqData1) != 'unmapped'])
  data <- Seurat::GetAssayData(seuratObjCite, assay = 'ADT', slot = 'counts')
  expect_equal(max(data[,41:100]), 0) #These have no data
  expect_equal(max(data[,1:40]), 6) 
	
	#Rename cells with a barcodeprefix:
	seuratObj@meta.data$BarcodePrefix <- c(rep('12345', 50), rep('67890', 50))
	seuratObj <- Seurat::RenameCells(object = seuratObj, new.names = paste0(seuratObj@meta.data$BarcodePrefix, '_', colnames(seuratObj)))
	
	#With prefix, no assay present
	seuratObjCite <- OOSAP:::AppendCiteSeq(seuratObj = seuratObj, countMatrixDir = inputPath1, barcodePrefix = '12345')
	expect_equal(colnames(seuratObjCite@assays$ADT), colnames(seuratObj@assays$RNA))
	expect_equal(rownames(seuratObjCite@assays$ADT), rownames(citeseqData1)[rownames(citeseqData1) != 'unmapped'])
	data <- Seurat::GetAssayData(seuratObjCite, assay = 'ADT', slot = 'counts')
	expect_equal(max(data[,41:100]), 0) #These have no data
	expect_equal(max(data[,1:40]), 6) 
	
	#Now add again, existing assay present:
	seuratObjCite2 <- OOSAP:::AppendCiteSeq(seuratObj = seuratObjCite, countMatrixDir = inputPath2, barcodePrefix = '67890', minRowSum = 0)
	d1 <- Seurat::GetAssayData(seuratObjCite, 'ADT', slot = 'counts')
	d2 <- Seurat::GetAssayData(seuratObjCite2, 'ADT', slot = 'counts')
	expect_equal(d2[,1:50], d1[,1:50]) #The original data, should be identical
	expect_equal(max(d2[,41:59]), 0) #These have no data
	expect_equal(max(d2[,60:100]), 1) 
	
	unlink(inputPath1, recursive = T)
	unlink(inputPath2, recursive = T)
})