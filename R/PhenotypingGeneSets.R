

#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
Phenotyping_GeneList <- function(QuickGO.path="./data/QuickGO"){

  #This Function is used to store lists of interesting gene sets. They can be useful for phenotyping various cell type and their states.

  #This is not a comprehensive list (yet); please contribute, if you can, with ref.

  #Some are literature derived, some are data-driven. However, ideally, if ref exists it should be here.

  ## Refs:
  #1) https://www.rndsystems.com/research-area/lymphoid-lineage-markers
  #2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1479811/
  #3) https://www.rndsystems.com/research-area/red-blood-cells--rbcs
  #4) https://www.bio-rad-antibodies.com/flow-myeloid-cell-immunophenotyping.html#1
  #5) https://www.bio-rad-antibodies.com/macrophage-m1-m2-tam-tcr-cd169-cd-markers-antibodies.html
  #6) 10.1126/SCIENCE.AAH4573
  #7) Yan 2017
  #8) Hermann 2018


  getwd()



  SGS.LS <- list()
  #### IxNG to IFNG
  SGS.LS$WBC                      <- c("PTPRC")

    #Lymphocytes generally are though of as T cells, B cells, and NK cells.
    SGS.LS$Lymphoid               <- c("PTPRC", "CD2", "CD5", "RPL21", "RPL23", "RPL27", "RPL31", "RPL35", "RPS14", "RPS21", "RPS24", "RPS3A", "RPS6", "TGT", "TPI1")

      SGS.LS$TCellCanonical       <- c("CD3G", "CD3D", "CD3E") #canonical TCR complex
      SGS.LS$TCellSecondary       <- c("CD2", "CD5", "CD6")
      SGS.LS$TCellTranscription   <- c("MAL", "LAT", "TCRIM", "CD28", "ZAP70", "FYN", "GATA3", "LEF1", "TCF7", "RUNX2", "STAT4", "SATB1")
        SGS.LS$CD8Canonical       <- c("CD8A", "CD8B")
        SGS.LS$MAIT               <- c("CD8A", "CD8B", "CD4", "KLRB1", "DPP4")
        SGS.LS$CD8Subphenos1      <- c("IL2RB", "KLRC1", "KLRG1", "PRF1", "GNLY", "GZMC", "GZMH", "TBX21")
        SGS.LS$CD4Canonical       <- c("CD4", "IL7R")
        SGS.LS$CD4Subphenos1      <- c("NOSIP", "RGS10", "FOXP3", "IL2RA", "LTB")
        SGS.LS$CD4Subphenos2      <- c("ANK3", "MXI1",  "CTSB")

      SGS.LS$BCellCanonical     <- c("CD19", "MS4A1", "CD79A")
      SGS.LS$BCellSecondary     <- c("CD20", "CD21", "CD22", "FCGR2B") #[2]
      SGS.LS$BCellCanonicalV3   <- c("CD24", "CD38", "CD72", "CD74", "CD79B", "CD83", "CD86") #[2]

      SGS.LS$NKCanonical        <- c("NCAM1", "CD8B", "GZMK", "LYAR", "NKG7", "GZMA", "GNLY", "FGFBP2", "FCGR3A", "CCL4", "GZMH")


    SGS.LS$Myeloid            <- c("ITGAM", "FUT4", "ANPEP", "CD14", "ITGAX","FCGR3A") #[4]

      SGS.LS$Macrophage         <- c("CD14", "ITGAX","FCGR1A", "CD68", "TFRC", "CD86", "CD163", "TLR2", "TLR4") #[5]
      SGS.LS$Monocytes          <- c("CSF1R", "CD14", "FCGR1A", "CD68", "S100A12", "MS4A7", "CKB", "LILRA3")
      SGS.LS$MonocytesFCGR3A    <- c("HES4", "CDKN1C", "FCGR3A", "MS4A7", "CKB", "LILRA3", "IFITM3", "MS4A4A", "LRRC25")
      SGS.LS$MonocytesCD34p     <- c("CD34", "THY1", "ENG", "KIT", "PROM1" )


    SGS.LS$HighlyActivated    <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "IFNG", "CD83", "CD82", "PLEK", "RGCC") #"ENSMMUG00000013779"
    #SGS.LS$HighlyActivated2   <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "IFNG") #the original from B.B.
    #SGS.LS$HighlyActivated3   <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "CD83", "CD82", "PLEK", "RGCC") #same as first without IFNG due to typo in original set
    # SGS.LS$HighlyActivated2   <- c("NR4A2", "NR4A1", "LMNA", "SDC4", "SLA", "JARID2",
    #                                "IRF4", "RGS16", "ENSMMUG00000045612", "TXNIP", "ENSMMUG00000031417", "LTB",
    #                                "GZMK", "IRF8", "c-fos", "RGS1", "KLF10", "ARID5B",
    #                                "CD83", "AHI1", "MYC", "SPAG9", "ID2", "TNFAIP3",
    #                                "HIVEP3", "BACH2", "ZFP36", "GZMA", "ENSMMUG00000041662", "ISG15",
    #                                "KLRB1", "SPINK2", "ENSMMUG00000042057", "SAMHD1", "APOL2", "TCF7",
    #                                "CD74", "S100A6", "S100A4", "IFI6", "CCL5", "ARL5B",
    #                                "CDKN1A", "CD160", "CENPF")
    SGS.LS$HighlyActivated2   <- c("NR4A1", "KLF10", "IRF8", "NR4A2", "RGS1", "AHI1", "ZFP36",
                                   "ID2", "BACH2", "MYC", "TNFAIP3", "KLRB1", "TCF7", "TXNIP", "LTB")
    SGS.LS$HighlyActivated3 <- unique(c(SGS.LS$HighlyActivated, SGS.LS$HighlyActivated2))


    SGS.LS$LessActivated      <- c("LTB", "IL7R", "PTPRCAP", "GZMK")
    SGS.LS$Pheno1             <- c("MT1M", "MGEA5")





    SGS.LS$Erythrocyte        <- c("HBB", "GYPA", "BCAM", "CD36", "EPO", "HPX", "SLC14A1", "LEPR") #[3]
    SGS.LS$Megakaryocytes     <- c("PF4", "GP9", "ITGA2B", "TMEM40", "LY6G6F","SEPT5", "PTCRA", "TREML1", "CLDN5", "HGD")
    SGS.LS$Eosinophils        <- c("CCR3")
    SGS.LS$PlasmaCells        <- c("SDC1")
    SGS.LS$Basophil           <- c("LMO4", "ENPP3", "IL3RA")
    SGS.LS$Neutrophil         <- c("CEBPE", "S100A8", "S100A9", "FUT4")
    SGS.LS$Dendritic          <- c("CD74")
    SGS.LS$Granulocyte        <- c("CSF2RB", "CSF3R", "IL1R2", "IL1RN", "IL8RB", "IL13RA1", "FPR1", "MME")

    # Blood Monocytes, DCs, and etc #[6]
    SGS.LS$Vilani_DC1_CD141CLEC9A       <- c("CLEC9A", "C1ORF54", "HLA-DPA1", "CADM1", "CAMK2D")
    SGS.LS$Vilani_DC2_CD1C_A            <- c("CD1C", "FCER1A", "CLEC10A", "ADAM8", "CD1D")
    SGS.LS$Vilani_DC3_CD1C_B            <- c("S100A9", "S100A8", "VCAN", "LYZ", "ANXA1")
    SGS.LS$Vilani_DC4_CD1CnegCD141neg   <- c("FCGR3A", "FTL", "SERPINA1", "LST1", "AIF1")
    SGS.LS$Vilani_DC5_AXLposSIGLEC6pos  <- c("AXL", "PPP1R14A", "SIGLEC6", "CD22", "DAB2")
    SGS.LS$Vilani_DC6_pDC               <- c("GZMB", "IGJ", "AK128525", "SERPINF1", "ITM2C")

    SGS.LS$Vilani_Mono1_classical_CD14high_CD16neg <- c("CD14", "VCAN", "S100A8", "S100A9", "FCN1", "ITGB2", "LRP1", "CSF3R", "TKT", "LYZ", "APLP2", "FPR1", "CD36", "S100A12", "CLEC4E", "ITGAM", "SLC2A3", "CTSD", "NEAT1", "PTAFR", "TREM1", "NAIP", "NCF1", "FCGR2A", "SCPEP1", "CTSA", "NLRP3", "ACSL1", "SDCBP", "SLC11A1", "IRS2", "VNN2", "DPYD", "CLEC7A", "BST1", "PLBD1", "PYGL", "QPCT", "BC013828", "CD163", "AQP9", "PELI1", "FAM198B", "GAS7", "STAB1", "CDA", "DOK3", "IRAK3", "PLAUR", "AL137655", "LILRA6", "TLR4", "AX747598", "TLR2", "AGTRAP", "CRISPLD2", "CCR1", "NFAM1", "ETS2", "RAB27A", "BNIP3L", "HPSE", "PER1", "MEGF9", "CD300E", "CYP1B1", "FCAR", "SOD2", "UPP1", "IER3", "C5AR1", "NLRP12", "SMA", "DMXL2", "NCF1B", "CREB5", "CR1", "ALDH1A1", "ASGR1", "FNDC3B", "DUSP6", "TOM1", "CDC42EP3", "ZBTB16", "DYSF", "KCNE3", "CD93", "CEBPD", "FCGR1A", "PLEKHM1", "CPM", "MPP7", "AK302511", "IL1B", "PFKFB3", "PLD3", "SMA3", "F13A1", "G0S2", "LOC100133161", "PHF21A", "TLR8", "CLMN", "TNFAIP3")
    SGS.LS$Vilani_Mono2_nonclassical_CD14posCD16high <- c("LAIR2", "ASAH1", "APOBEC3A", "TSPAN14", "LIPA", "CYTIP", "SIGLEC10", "LILRB1", "EMR1", "TTYH3", "CAMKK2", "CX3CR1", "C3AR1", "BC013828", "RASGEF1B", "BIRC3", "PLIN2", "CD300C", "CD83", "XYLT1", "KLF2", "FBP1")
    SGS.LS$Vilani_Mono3_undef1 <- c("G0S2", "NAMPT", "NEAT1", "AL137655", "CSF3R", "FCGR3B", "SRGN", "TREM1", "TNFRSF10C", "MXD1", "SOD2", "CXCR2", "SLC25A37", "S100A8", "FPR1", "ITM2B", "MNDA", "VNN2", "SDCBP", "CXCR1", "S100A9", "AQP9", "SORL1", "ACSL1", "AX747598", "R3HDM4", "NCF1", "IFITM2", "FCGR2A", "XPO6", "GCA", "C5AR1", "TKT", "PELI1", "SLC2A3", "CLEC4E", "MMP25", "GLUL", "CD14", "LOC388312", "NCF1C", "VMP1", "RTN3", "ACTN1", "PTAFR", "S100A12", "SEC14L1", "DQ574721", "LITAF", "TLR2", "SHKBP1", "LIMK2", "LOC100505702", "PYGL", "RNF24", "DNAJC25-GNG10", "IL8", "FPR2", "LOC731275", "SLC12A6", "IL1R2", "VNN3", "CFD", "VCAN", "BC013828", "NAIP", "ZBTB16", "BCL2A1", "FAM129A", "PLAUR", "FNDC3B", "FP15737", "SEPX1", "LOC100133161", "PER1", "FBXL5", "IL17RA", "TLR4", "IGF2R", "ITGAM", "HIST1H2AC", "LRP1", "KREMEN1", "C12ORF35", "PRRG4", "CR1", "RAB27A", "LOC100505815", "BST1", "NUMB", "USP15", "CDA", "IER3", "ACADSB", "DYSF", "PXN", "PDP2", "TNFRSF1A", "LRG1", "LOC91948", "FLJ45445", "SMAP2", "LOC643802", "NINJ1", "ABTB1", "CCNY", "TMEM154", "CCR1", "CARD8", "TACC3", "TMEM71", "PTGS2", "HPSE", "C3ORF72", "FAM157A", "AK130076", "CD163", "NBEAL2", "IL1RAP", "GK", "AZGP1P1", "DOK3", "PROK2", "FAM115C", "QPCT", "ALPL", "BEST1", "CES3", "CREB5", "SPAG9", "GPR97", "TBL1X", "FAM198B", "FCAR", "PHF21A", "IRS2", "CYP1B1", "NCF1B", "BC048113", "BACH1", "AX747405", "RCBTB2", "CEBPD", "ALPK1", "LAT2", "OSBPL8", "PCNX", "LPPR2", "CCPG1", "DOCK5", "TUBA4A", "F2RL3", "NCF4", "FAM157B",
                                    "TECPR2", "SLA", "TM6SF1", "CRISPLD2", "FAS", "PADI4", "RUFY1", "AK302511", "PDE4B", "AK091866", "DQ580909", "FAM126B", "LRP10", "PADI2", "TRIB1", "ZDHHC18", "F5", "PDLIM7", "RBM47", "SIRPA", "ARHGAP26", "DSTYK", "TLR6", "FBXL13", "LOC649305", "P2RY8", "HBP1", "SGSM1", "ABCA1", "SEMA4D", "ABHD5", "MRS2P2")
    SGS.LS$Vilani_Mono4_undef2 <- c("PRF1", "GNLY", "KLRC4-KLRK1", "TCRBV3S1", "CTSW", "CCL5", "KLRD1", "FGFBP2", "NKG7", "IL2RB", "SPON2", "HOPX", "GZMA", "CST7", "ZAP70", "GPR56", "SYNE2", "KLRF1", "GZMH", "IL32", "TXK", "IFITM1", "IKZF3", "LCK", "TC2N", "S1PR5", "S100A8", "MCTP2", "S100A12", "CD96", "SAMD3", "TRGC2", "TTC38", "PXN", "S100A9", "SH2D1B", "LAIR2", "SYNE1", "PRKCH", "RARRES3", "PIK3R1", "CCL4", "PARP8", "TGFBR3", "GSTM1", "CD2", "CD247", "PDE4D", "PRDM1", "CBLB", "GIMAP1", "BC013828", "DENND2D", "GZMM", "SKAP1", "TMEM41A", "KLRB1", "PLEKHG3", "FCRL6", "PYHIN1", "AAK1", "CCR1", "IRS2", "STAT4", "IL18RAP", "INADL", "DIP2A", "LOC388692", "FAIM3", "CD160", "PAPD5", "PAM", "PIK3IP1", "PRSS23", "PVRIG", "VNN2", "CREB5", "CCND2", "RORA", "ATXN7", "PTPN4", "LIMK2", "SEPX1", "KLF12", "TRDC", "AK094156", "NCR3", "KIF21B", "PTGDR", "IER3", "ITK", "BTN3A2", "CPD", "NCAM1", "ZBTB16", "RAB27A", "RUNX3", "SLC25A37", "SLFN13", "GCA", "RASA3", "IPCEF1", "SCML4", "NID1", "PADI4", "S1PR1", "ZBTB38", "FCGR1A", "PARP15", "ETS1", "LAT", "TRPM2", "FNDC3B", "CCL3", "CLEC4D", "OPTN", "RASSF3", "LOC100216546", "IL1B", "GBP5", "ENC1", "KLRG1", "SYTL3", "BC051736", "TRAPPC10", "LIN54", "LOC374443", "ZNF44", "F2R", "TFDP2", "CEP78", "CXCR2", "G0S2", "GABARAPL1", "TUBD1", "PDPR", "DQ573668", "FXYD6-FXYD2", "BRF2", "SLAMF6", "CREM", "TGIF1", "SLFN5", "ARHGAP24", "ZMYM5", "ZNF276", "SUPV3L1", "FAM190B", "LPIN1")





  #[7]
  SGS.LS$Testis$Yan2017_GonadSomatic <- c("KIT", "POU5F1", "NANOG", "DAZL",
                                          "DDX4", "STRA8", "SYCP1", "ZP3",
                                          "WT1", "AMH", "ARX", "CYP17A1")

  #[8]

  SGS.LS$Testis$Hermann2018_Human_Pertitubular_Sertoli_Leydig <- c("DCN", "SPARCL1", "C11orf96", "APOD", "BST2", "PLAC9", "SMOC2", "KCNMB1", "BGN", "MFAP4", "INMT", "SFRP4", "TNS1", "COL1A1", "LMOD1", "ANGPTL1", "OGN", "DPT", "CYBRD1", "CARMN", "EMILIN1", "PODN", "TMEM100", "THBD", "CYR61", "CTSK", "PTRF", "ZFP36", "KRT13", "RGN", "CILP", "ENG", "CPED1", "AEBP1", "C1S", "TCF21", "CBLN4", "TMEM47", "MYC", "COL1A2", "CAV1", "NR2F2", "SELENBP1", "TGFB1I1", "TSHZ2", "OSR2", "APOC1", "TCEAL3", "LRP1", "SH3BGRL", "TIMP3", "RASD1", "MIR503HG", "OSR1", "CTGF", "COL15A1", "C1R", "TM4SF1", "PLPP3", "MGP", "MMP2", "S100A6", "NUPR1", "HNMT", "CLMP", "LCNL1", "FXYD6", "FOS", "ISLR", "HHIP-AS1", "NEAT1", "ABCA8", "FLNA", "ARX", "IFITM3", "JUNB", "EGR1", "RAB34", "F10", "TSHZ1", "COL6A1", "TIMP1", "EPDR1", "LAMC3", "COL6A2", "NID1", "FBLN5", "ALDH1A3", "LGALS3BP", "CFH", "TCEAL1", "PRPS2", "ADIRF", "MAOB", "SSPN", "GYPC", "PDGFRB", "PRELP", "SERPINF1", "WBP5", "ACTA2", "PRDX2", "LTBP2", "GADD45G", "PDGFRA", "FSTL1", "MYL9", "KRT17", "FKBP9", "FRZB", "SMIM3", "PKDCC", "TNFSF12", "IL34", "ALDH1A1", "FLRT2", "IFI6", "ITGA9", "FAM127A", "FBLIM1", "COL18A1", "PLAGL1", "LTBP4", "BTG2", "MAGED2", "RARRES2", "FOSB", "TAGLN", "ARMCX1", "GPNMB", "SEPP1", "VEGFB", "IGF1", "LGALS1", "IGFBP7", "COL6A3", "ARMCX3", "TCEAL8", "EHD2", "SPARC", "CALD1", "FHL2", "BHMT2", "TBX2", "F3", "CBR3", "MARVELD1", "PTGIS", "PDLIM3", "PPP1R14A", "NGFRAP1", "CRISPLD2", "CYP26B1", "CAMK2N1", "GALM", "RHOB", "PDLIM2",
                                                                   "IL32", "GSTM5", "GSN", "TPM2", "PTMS", "MYLK", "MAF", "NBL1", "APOBEC3C", "CHGA", "TCEAL7", "APOE", "HTRA1", "LINC01420", "COMT", "LAMB1", "VWA1", "KANK2", "CBX6", "PLA2G5", "IDS", "MRGPRF", "ALDH2", "ANGPTL4", "IFI44L", "B2M", "TXNIP", "WNT2B", "SAT1", "TSPO", "SVIL", "TCEAL4", "HSPB6", "ANKDD1A", "PEMT", "ELN", "CPE", "PROK1", "CDH11", "AHNAK", "MFGE8", "HEXB", "AR", "CD63", "MMP23B", "MYH11", "CD151", "GSTP1", "IGFBP4", "AK3", "COX7A1", "IGF2", "ST6GALNAC4", "SNED1", "LAG3", "JUN", "DPEP1", "AKR1C3", "IFI27", "FAM210B", "AOC3", "HIC1", "CCDC8", "TRIP6", "SYNPO2", "FAM114A1", "LTBR", "MAGEH1", "ZFP36L2", "PTCH2", "IGFBP6", "THY1", "ANKRD35", "RPL39", "CSRP1", "SERPING1", "IRF2", "FXYD1", "MXRA8", "COLEC11", "NDRG2", "KRT19", "VCAN", "CCDC85B", "MYOCD", "TNS2", "RPL36A", "TPM1", "VIM", "FBN1", "EMP3", "NR2F2-AS1", "EGFR", "MAP1B", "TMSB4X", "CNN2", "PLEKHH2", "RPL10", "CD99", "CCND2", "RAB31", "TMSB10", "NEXN", "NDN", "RPS4X", "ASS1", "ITM2B", "GNB5", "S100A4", "RHOBTB3", "MT1X", "SCN7A", "HGSNAT", "PHYHD1", "MGST3", "PMP22", "LINC01197", "LY6E", "MEG3", "DHRS7", "ABHD14B", "CD81", "SLC22A17", "EDNRA", "CFD", "PPP1R15A", "PSME1", "CPQ", "DUSP1", "CAV2", "IER2", "SLC25A6", "IFITM2", "LGALS3", "AKR1C2", "KLHL42", "EEF1A1", "COL3A1", "ACTB", "INSR", "PLTP", "ITM2A", "PDGFRL", "EID1", "NCKAP5L", "MYL6", "CST3", "RPS8", "ITGB1", "ADCY5", "RPL41", "SLC40A1", "RPL27A", "LAMP2", "PLOD1", "MALAT1", "PLEC", "TRIM56",
                                                                   "SPTBN1", "SGCE", "RPL13A", "RPL10A", "FTL", "HSPG2", "RPL23", "PSAP", "GAMT", "CRTAP", "SFRP1", "RPL34", "VKORC1", "FCGRT", "PTGDS", "RPL31", "EMP2", "RPL7", "ST6GALNAC6", "FBLN2", "RPS27A", "RPS9", "RPL6", "SLC44A1", "ANXA5", "ABLIM1", "THSD4", "DDIT4", "RPL5", "MGLL", "SMIM4", "HYI", "HLA-E", "BDH2", "MARCKSL1", "FBXO17", "KCTD12", "RPL30", "RPS3", "FKBP10", "ARL6IP5", "RHOC", "AKR1C1", "RPL9", "CSTB", "RPL12", "HLA-C", "PARVA", "RPLP1", "RP11-14N7.2", "RPS3A", "SLC51A", "RPL21", "FN1", "MATN2", "RPL3", "MSN", "NFKBIA", "RPS17", "CYP11A1", "RPS23", "TPT1", "RPL22", "NDUFA4", "RPL26", "ITM2C", "SNHG18", "PGRMC1", "PTCH1", "NFIB", "TIMP2", "PRRT2", "PLXDC2", "BLVRB", "ZFP36L1", "RPS15A", "PID1", "RPS14", "PPIB", "NACA", "COL4A2", "GADD45B", "ATP6AP2", "LTBP1", "NFIC", "RPS4Y1", "SERINC5", "PNRC1", "FTH1", "MXRA7", "RPL18", "ITIH5", "RPS13", "LAMA2", "LAYN", "KLF2", "SERPINH1", "SOD3", "CITED2", "NDUFB11", "ANXA2", "CCDC3", "RPS15", "SLC38A2", "TAX1BP3", "NPDC1", "RPL24", "LDHB", "RPL36", "ZNF358", "RPL14", "CEBPD", "CPNE3", "CIRBP", "RPL37", "EEF1D", "IFI27L2", "RPS24", "COX8A", "RPL11", "RPL18A", "DSTN", "ARHGDIB", "GADD45A", "RPL23A", "RPL35A", "CDO1", "CYB5R3", "ATRX", "LRP10", "EPHX1", "RPL27", "EIF4A2", "AXL", "TUBB", "MORF4L2", "CEBPB", "RPL29", "LEPROT", "RPS18", "PLD3", "NUCKS1", "MYH9", "RPS16", "RPS11", "WFDC1", "RPS28", "EIF2S3", "HNRNPA1", "DPYSL3", "FAM198B", "RPL19", "DBP", "CNN3", "TMEM43",
                                                                   "TPM4", "RPL8", "PPIA", "MSRB2", "PTN", "RPS27", "CAMK1", "TXN", "NFIA", "BICC1", "HSPB1", "CTSC", "ECE1", "RPL32", "HSPA1B", "CRYAB", "ATP6V0E1", "NME4", "GAS6", "C3orf70", "TRBC2", "DDAH2", "FAU", "PEBP1", "VAMP5", "ST3GAL4", "SMTN", "KDSR", "LAMB2", "GNB2L1", "CCND1", "HLA-B", "PCBD1", "EIF2AK2", "RPS27L", "ANXA7", "PHACTR2", "RPS25", "SMOC1", "ANXA6", "PDLIM7", "LTBP3", "TPI1", "SCPEP1", "TMEM98", "AK1", "MYADM", "SERTAD1", "LMNA", "RPS5", "NT5DC2", "TMEM14C", "AHDC1", "APP", "TSC22D3", "FAM46A", "RP11-553L6.5", "TMEM109", "MT-ND3", "RPL28", "SNHG7", "DNASE1L3", "EBF1", "MLXIP", "SAT2", "SHISA5", "TOMM7", "SRSF5", "RPLP2", "RPS7", "RASSF2", "PPDPF", "SRSF8", "PDE5A", "H1F0", "RPS2", "C12orf57", "GNG11", "MSRB3", "SELM", "RPS19", "CDC42EP1", "CPM", "NREP", "S100A13", "RNF213", "ALDOA", "MSI2", "HNRNPA0", "SLC29A1", "GSTA4", "NENF", "WIPF3", "SERPINE2", "GAPDH", "SDCBP", "SLC25A3", "FOXS1", "RPL13", "VAMP2", "SCP2", "ZFAS1", "GK5", "BTF3", "GRHPR", "BNC2", "RPS12", "PTMA", "MRFAP1", "PRDX1", "TSC22D4", "PLEKHA4", "TBC1D2B", "EEF2", "HSD17B11", "MT-CYB", "NPTX2", "CCNL1", "NDUFA1", "PRNP", "STAT3", "TUBA1A", "SERPINA5", "CBR1", "PGAM1", "VASP", "FGFR1", "SNX21", "JTB", "PCOLCE", "SLC25A5", "ATP5G2", "BEX4", "MIF", "RPS6", "MT-ATP6", "FAM50A", "AES", "N4BP2L2", "RBM3", "CCNI", "HCFC1R1", "COX7B", "RPS21", "SNHG8", "NDUFS5", "NME3", "JUND", "LINC01116", "CLEC11A", "HSD17B10", "SMARCA1", "KDELR1", "RUNX1T1", "IGFBP2",
                                                                   "ATP5I", "GLTSCR2", "NR1D2", "FUS", "PSRC1", "RPL4", "PEG3", "RPL17", "SEMA3D", "RGMA", "MIR99AHG", "RBMX", "CCNG1", "HLA-A", "OAZ1", "SLC25A4", "TBCA", "SNTA1", "KLF9", "RPLP0", "MT-ND5", "MTCH1", "MZT2A", "CREB3L2", "ATP1B1", "SON", "NISCH", "MYL12A", "RPS10", "LUC7L3", "EEF1B2", "ADH5", "PHLDA3", "CTSD", "ST5", "DMD", "DNPH1", "LAMTOR1", "GSTO1", "SH3KBP1", "ESD", "ZNF580", "PCYOX1", "CHCHD2", "GCSH", "EFEMP2", "SOD1", "DDR2", "SMARCD3", "ACTG1", "NSMF", "NEDD8", "FBXO21", "MT-ND4", "CMIP", "PFN1", "RPL7A", "CD59", "SRP14", "EIF3E", "TACC1", "RPL15", "MAGEF1", "VCL", "GALNT10", "PLS3", "ASAH1", "MT-ND1", "PLP2", "TSPAN8", "PFKL", "CA11", "RBBP7", "NAP1L1", "DDX17", "EDF1", "JAM3", "ZMAT3", "SYNE2", "RILPL2", "MGST1", "TSPAN4", "WLS", "TMEM256", "PNISR", "ZFYVE21", "ERGIC1", "FUCA2", "SDHD", "PPIC", "CFL1", "SQSTM1", "CUTA", "CKB", "SEPW1", "H1FX", "GLO1", "PKM", "RPL35", "VGLL4", "CLIC1", "SUCLG2", "HSP90AB1", "RPS29", "PDK4", "NPM1", "ADI1", "DST", "MT-ND2", "EIF3L", "CYB5A", "PRKAR2B", "ACTN1", "ZBTB7A", "HMGN3", "SH3BGRL3", "PYGB", "KMT2A", "HMGB1", "KIDINS220", "NFIX", "S100A11", "AC013461.1", "QTRT1", "NOTCH2", "TLN1", "HSP90B1", "MT-CO3", "CTSL", "LMO4", "RAC1", "COL16A1", "PRDX3", "CALM1", "ATP5G1", "ATP2B4", "SEC22B", "ENAH", "ATP5G3", "MYLIP", "CNP", "STX10", "ATP5D", "PRRC2B", "DYNLRB1", "ADAR", "LAPTM4A", "ILK", "TMEM14A", "TMEM165", "LGMN", "CLTB", "CBX7", "DIO3OS", "FERMT2", "IER5", "PHB2", "SERINC1", "PLGRKT", "MDK", "RASL11B")

  SGS.LS$Testis$Hermann2018_Human_Pertitubular_Sertoli_Leydig_top10 <- c("DCN", "TIMP1", "ACTA2", "MYL9", "TAGLN", "IGFBP7", "TPM2", "APOE", "TMSB4X", "PTGDS")
  SGS.LS$Testis$Hermann2018_Human_Perivascular_Macrophage_top10     <- c("CD74", "SPARCL1", "HLA-DRA", "HLA-DPA1", "MGP", "B2M", "TMSB10", "VIM", "RPL10", "TMSB4X")
  SGS.LS$Testis$Hermann2018_Human_UndiffSg1_top10                   <- c("SMS", "ASB9", "CST3", "PFN1", "PABPC4", "RPL18A", "RPS19", "CCNI", "MT-CO1", "MT-CO3")
  SGS.LS$Testis$Hermann2018_Human_UndiffSg2_top10                   <- c("ZNF428", "HMGA1", "RPS4X", "GNB2L1", "RPS12", "EEF1B2", "RPL18", "CCNI", "RPSA", "RPS2")
  SGS.LS$Testis$Hermann2018_Human_DifferentiatingSg1_top10          <- c("DMRT1", "PRAME", "IGFBP2", "TKTL1", "CIRBP", "ANP32B", "PTMA", "HIST1H4C", "HMGB1", "RPS2")
  SGS.LS$Testis$Hermann2018_Human_PrelepSct_top10                   <- c("TEX101", "ZCWPW1", "SCML1", "SMC1B", "MEIOB", "SYCP3", "RAD51AP2", "HIST1H4C", "SMC3", "SYCP1")
  SGS.LS$Testis$Hermann2018_Human_LeptoteneZygotene_top10           <- c("TDRG1", "LINC01206", "KB-1592A4.15", "TEX101", "TPTE", "RCN2", "MT-ND4", "MT-ND3", "MT-ATP6", "RP11-620J15.3")
  SGS.LS$Testis$Hermann2018_Human_PachSct_top10                     <- c("CETN3", "CCDC112", "RP11-620J15.3", "OVOS2", "CALM2", "RNFT1", "HSPA2", "SPATA8", "TMEM99", "C15orf48")
  SGS.LS$Testis$Hermann2018_Human_Diplotene_2ndSct_top10            <- c("AURKA", "CCNA1", "CLGN", "PRKCDBP", "CIAPIN1", "PPP3R2", "RP5-1023B21.1", "ASRGL1", "ZMYND10", "SPATA8")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid1_top10             <- c("RP5-1023B21.1", "PRKCDBP", "ISOC2", "TBPL1", "CCNA1", "LYAR", "AURKA", "PPP3R2", "CCDC42", "GLIPR1L1")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid2_top10             <- c("CDRT15", "LINC00643", "CFAP126", "C11orf97", "MNS1", "C9orf116", "MCHR2-AS1", "GOLGA6L2", "LYRM5", "CCDC110")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid3_top10             <- c("FAM24A", "EQTN", "LYZL1", "C4orf17", "LRRC39", "TJP3", "NLRP1", "RP11-639B1.1", "RP13-39P12.3", "GOLGA6L2")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid4_top10             <- c("DPP10-AS3", "LRRC39", "FAM24A", "RP11-639B1.1", "RP13-39P12.3", "PEX5L-AS2", "C7orf61", "LYZL1", "LACE1", "GOLGA6L2")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid5_top10             <- c("TMEM190", "CCDC168", "TEX29", "OLAH", "TEX43", "ACRV1", "ERICH2", "SPACA1", "FAM209B", "FAM209A")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid6_top10             <- c("LINC00919", "FAM71A", "FAM71B", "PRSS58", "RNF151", "HEMGN", "CCDC179", "C1orf100", "C17orf105", "TNP1")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid7_top10             <- c("PRM2", "PRM1", "LELP1", "AC007557.1", "TSSK6", "C10orf62", "GLUL", "SPATA3", "TNP1", "OAZ3")

  SGS.LS$Testis$Hermann2018_Human_UndiffSg_Marker                     <- c("NANOS3", "GFRA1", "ETV5", "ID4", "NANOS2", "SOHLH1", "RHOXF1", "MAGEA4", "ZBTB16", "BCL6B", "UCHL1", "FOXO1", "DMRT1", "TAF4B", "PIWIL4", "PIWIL2")
  SGS.LS$Testis$Hermann2018_Human_DifferentiatingSgSct_Marker         <- c("KIT", "SOHLH2", "STRA8", "SCML2", "SYCP3", "DMC1", "PEG10", "GAL3ST1", "TOP2B", "SPO11", "MEIOB", "DAZL", "ATR", "NEUROG3", "SYCP1", "MYBL1", "HORMAD1", "PIWIL1")
  SGS.LS$Testis$Hermann2018_Human_RoundSpermatid_Marker               <- c("ACR", "PGK2", "SAXO1", "ACRV1", "CATSPER1", "CATSPER3", "HSFY2", "RBM5", "CHD5", "TEKT1", "TXNDC8", "GAPDHS", "PRM1", "CA2", "PRM2", "PRM3", "TNP1", "TNP2")
  SGS.LS$Testis$Hermann2018_Human_Perivascular_Macrophage_Marker      <- c("CLU", "PECAM1", "CDH5", "AIF1", "CD68", "CX3CR1", "WT1")
  SGS.LS$Testis$Hermann2018_Human_Pertitubular_Sertoli_Leydig_Marler  <- c("GATA4", "ACTA2", "MYH11", "CYP11A1", "NR5A1", "STAR")
  SGS.LS$Testis$Hermann2018_Human_Spermatogenesis_Markers             <- c("NANOS3", "GFRA1", "ETV5", "ID4", "NANOS2", "SOHLH1", "RHOXF1", "MAGEA4", "ZBTB16", "BCL6B", "UCHL1", "FOXO1", "DMRT1", "TAF4B", "PIWIL4", "PIWIL2", "KIT", "SOHLH2", "STRA8", "SCML2", "SYCP3", "DMC1", "PEG10", "GAL3ST1", "TOP2B", "SPO11", "MEIOB", "DAZL", "ATR", "NEUROG3", "SYCP1", "MYBL1", "HORMAD1", "PIWIL1", "ACR", "PGK2", "SAXO1", "ACRV1", "CATSPER1", "CATSPER3", "HSFY2", "RBM5", "CHD5", "TEKT1", "TXNDC8", "GAPDHS", "PRM1", "CA2", "PRM2", "PRM3", "TNP1", "TNP2")

  SGS.LS$Testis$Guo18_Sertoli             <- c("SOX9", "AMH", "WFDC2", "BEX2", "PRND")
  SGS.LS$Testis$Guo18_Endothelial             <- c("VWF", "NOTCH4", "HES1", "JAG1", "PECAM1", "MAML1")
  SGS.LS$Testis$Guo18_Myoid          <- c("MYH11", "ACTA2", "PTCH1", "PTCH2", "GLI1", "IGFBP6")
  SGS.LS$Testis$Guo18_Leydig          <- c("DLK1", "IGF1", "INHBA", "VIT", "IGFBP3", "IGFBP5")



  AnnotationFiles <- list.files(QuickGO.path, pattern = "QuickGO", full.names = T)
  GeneLists <- list()

  for(AnnFile in AnnotationFiles){
    # AnnFile = AnnotationFiles[1]
    tempName = gsub(".txt", "", gsub("_","",gsub("annotations-", "", gsub(".tsv","",gsub("-","",basename(AnnFile))))))
    GeneLists$Extra[[tempName]] <-  read.table(
      AnnFile,
      sep="\t", header=TRUE, row.names = NULL, fill = TRUE )
  }

  SGS.LS$QuickGOgenes <- as.character(data.table::rbindlist(lapply(GeneLists$Extra, function(setX){
    subset(setX, TAXON.ID == 9606)[,c("GO.NAME", "SYMBOL")] #10090 = mouse ; 9606 = human
  }))$SYMBOL)

  return(SGS.LS)
}



#' @title A Title
#'
#' @description A description
#' @param SeurObj, A Seurat object.
#' @return A modified Seurat object.
#' @keywords SerIII_template
#' @export
RhesusGeneDavidConv <- function(ColNames2Conv, RhesusConvDavid.path, ENSMB.tag = "ENSMM", returnFull=F){

  #Seurat 3 cant update the gene names !!
  # see https://github.com/satijalab/seurat/issues/1207

  print("Reading in David Data...")
  noGeneSYM <- ColNames2Conv[grepl(ENSMB.tag, ColNames2Conv)]


  David6.8ConvTable <- data.frame(read.csv(RhesusConvDavid.path, sep = "\t", header = T))
  rownames(David6.8ConvTable) <- David6.8ConvTable$From
  David6.8ConvTable <- David6.8ConvTable[noGeneSYM, ]
  length(unique(noGeneSYM)); length((noGeneSYM))
  rownames(David6.8ConvTable) <- noGeneSYM

  David6.8ConvTable$Final <- as.character(David6.8ConvTable$To)

  David6.8ConvTable$Final[which(is.na(David6.8ConvTable$To))] <- rownames(David6.8ConvTable)[which(is.na(David6.8ConvTable$To))]

  duplicatedGeneNames <- names(which(table(David6.8ConvTable$Final)>1))

  #change the second duplace name to add a .2
  #perhaps can avg the expr?
  for(geneN in duplicatedGeneNames){

    David6.8ConvTable$Final[grep(geneN , David6.8ConvTable$Final)][2] <- paste(geneN, ".2", sep="")

  }


  if(returnFull) return(David6.8ConvTable) else return(David6.8ConvTable$Final)


}
