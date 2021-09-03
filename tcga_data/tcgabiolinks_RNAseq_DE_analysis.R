library(TCGAbiolinks)
library(SummarizedExperiment)

outFolder <- file.path("TCGABIOLINKS_RNASEQ_DE_ANALYSIS")
dir.create(outFolder, recursive = TRUE)

# some refs
# https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Transcriptomic_analysis
# https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/
# https://support.bioconductor.org/p/63981/#64004  # -> explain limma/voom on tcga values

# http://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html

#https://stat.ethz.ch/pipermail/bioconductor/2013-November/056309.html
#If FPKM is really all you have, then convert the values to a log2 scale 
#and do an ordinary limma analysis as you would for microarray data, using 
#eBayes() with trend=TRUE. 

## RNASEQ ANALYSIS FOR MICROARRAY
## https://support.bioconductor.org/p/67590/

##############################################
##############################################  LEGACY = TRUE workflow (hg19)
##############################################

if(buildTable){

  # You can define a list of samples to query and download providing relative TCGA barcodes.
  listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                   "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                   "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                   "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                   "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")
  
  # Query platform Illumina HiSeq with a list of barcode 
  query <- GDCquery(project = "TCGA-BRCA", 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    experimental.strategy = "RNA-Seq",
                    platform = "Illumina HiSeq",
                    file.type = "results",
                    barcode = listSamples, 
                    legacy = TRUE)
  
  # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
  GDCdownload(query)
  
  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
  # rsem.genes.results as values
  BRCARnaseqSE <- GDCprepare(query)
  
  outFile <- file.path(otuFolder, "BRCARnaseqSE.Rdata")
  save(BRCARnaseqSE, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(otuFolder, "BRCARnaseqSE.Rdata")
  BRCARnaseqSE <- get(load(outFile))
}


BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
dim(BRCAMatrix)
# [1] 19947    10
stopifnot(ncol(BRCAMatrix) == length(listSamples))


# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)

# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
# save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = BRCARnaseq_CorOutliers,  
                                      geneInfo =  geneInfo)  ### ! from legacy=FALSE, done geneLength

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT")) # Solid Tissue Normal
 
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP")) # PRIMARY SOLID TUMOR

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])
# pipeline =  ("limma" or "edgeR") default: "edgeR"
# method = 'glmLRT' (1) or 'exactTest' (2) used for edgeR
# (1) Fit a negative binomial generalized log-linear model to the read counts for each gene
# (2) Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.


##############################################
##############################################  LEGACY = TRUE home-made limma-voom DE analysis
##############################################
# https://f1000research.com/articles/5-1408
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

library(edgeR)
library(limma)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
dim(BRCAMatrix)
# [1] 19947    10
stopifnot(ncol(BRCAMatrix) == length(listSamples))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(BRCAMatrix),
                                   typesample = c("NT")) # Solid Tissue Normal

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(BRCAMatrix), 
                                   typesample = c("TP")) # PRIMARY SOLID TUMOR

BRCAMatrix <- BRCAMatrix[,c(samplesNT, samplesTP)]

minCpm_filt <- 0

cpm_mat <- cpm(BRCAMatrix)
rows_to_keep <- which(rowSums(cpm_mat) > minCpm_filt)
filt_mat <- BRCAMatrix[rows_to_keep,]

# 1 for symbol, 2 for entrezID
i_symb <- 1
geneSymbols <- sapply(strsplit(x=rownames(filt_mat), split="\\|"), function(x)x[[i_symb]])
dup_genes <- geneSymbols[duplicated(geneSymbols)]
rows_to_keep2 <- which(!geneSymbols %in% dup_genes)
filt_mat2 <- filt_mat[rows_to_keep2,]
geneSymbols2 <- sapply(strsplit(x=rownames(filt_mat2), split="\\|"), function(x)x[[i_symb]])
stopifnot(!duplicated(geneSymbols2))
rownames(filt_mat2) <- geneSymbols2

samples_groups <- sapply(colnames(filt_mat2), function(x) {
  if(x %in% samplesNT) return("normal")
  if(x %in% samplesTP) return("tumor")
  sys.exit()
})

seqData <- DGEList(filt_mat2, group=samples_groups, genes = geneSymbols2)
seqData <- calcNormFactors(seqData)

boxplot(cpm(seqData, log=TRUE), las=2, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")

my_group_design <- factor(samples_groups, levels = c("normal", "tumor"))
my_design <- model.matrix( ~ my_group_design)

par(mfrow=c(1,2))
v <- voom(seqData, my_design, plot=TRUE)
vfit <- lmFit(v, my_design)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

DE_topTable <- topTable(efit, coef=ncol(v$design), number=Inf, sort.by="p") ## if not 0+ in design -> coef=2

dataDEGs$genes <- rownames(dataDEGs)
DE_topTable$genes <- rownames(DE_topTable)
cmp_dt <- merge(dataDEGs, DE_topTable, by="genes")

plot(cmp_dt$logFC.x, cmp_dt$logFC.y)
plot(cmp_dt$PValue, cmp_dt$P.Value)# !!!!!!???????

##### with the other way to design the design matrix
# https://bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
# my_design2 <- model.matrix( ~ 0+my_group_design)
# contr.matrix <- makeContrasts(
#   NormalVsTumor = my_group_designnormal-my_group_designtumor,
#   levels = colnames(my_design2))
# 
# v2 <- voom(seqData, my_design2, plot=TRUE)
# vfit2 <- lmFit(v2, my_design2)
# vfit2 <- contrasts.fit(vfit2, contrasts=contr.matrix)
# efit2 <- eBayes(vfit2)
# DE_topTable2 <- topTable(efit2,  number=Inf, sort.by="p", coef=NULL) 
# DE_topTable2$genes <- rownames(DE_topTable2)
# 
# 
# cmp_dt <- merge(DE_topTable, DE_topTable2, by="genes")
# 
# plot(cmp_dt$adj.P.Val.x, cmp_dt$adj.P.Val.y)
# stopifnot(sum(abs(cmp_dt$adj.P.Val.x - cmp_dt$adj.P.Val.y))< 10^-6)



dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            pipeline = "limma",
                            voom = TRUE)
dataDEGs$genes <- rownames(dataDEGs)
DE_topTable$genes <- rownames(DE_topTable)
cmp_dt <- merge(dataDEGs, DE_topTable, by="genes")

plot(cmp_dt$logFC.x, cmp_dt$logFC.y)
plot(cmp_dt$adj.P.Val.x, cmp_dt$adj.P.Val.y)# !!!!!!???????

##############################################
##############################################  LEGACY = FALSE workflow (hg38)
##############################################

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
dir.create(DataDirectory,recursive = TRUE)
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")


if(buildTable) {
  query <- GDCquery(project = CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  samplesDown <- getResults(query,cols=c("cases"))
  
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                    typesample = "TP")
  
  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                    typesample = "NT")
  dataSmTP_short <- dataSmTP[1:10]
  dataSmNT_short <- dataSmNT[1:10]
  
  queryDown <- GDCquery(project = CancerProject, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        barcode = c(dataSmTP_short, dataSmNT_short))
  
  GDCdownload(query = queryDown,
              directory = DataDirectory)
  
  dataPrep <- GDCprepare(query = queryDown, 
                         save = TRUE, 
                         directory =  DataDirectory,
                         save.filename = FileNameData)
  
  outFile <- file.path(outFolder, "dataPrep.Rdata")
  save(dataPrep, file=outFile)
  cat(paste0("... written: ", dataPrep, "\n"))
  
} else {
  outFile <- file.path(outFolder, "dataPrep.Rdata")
  dataPrep <- get(load(outFile))
  
}



dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") ### ! from legacy=FALSE, done gcContent

boxplot(dataPrep, outline = FALSE)

boxplot(dataNorm, outline = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
                            mat2 = dataFilt[,dataSmNT_short],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,dataSmTP_short],dataFilt[,dataSmNT_short])











########################################################
######################################################## ENRICHMENT ANALYSIS
########################################################

library(TCGAbiolinks)
# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)

genelist <- rownames(dataDEGsFiltLevel)
library(EnsDb.Hsapiens.v86)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= genelist, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",geneIDs1$SYMBOL)

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist, 
                        nBar = 10)

sys.exit(0)


########################################################
######################################################## SURVIVAL ANALYSIS
########################################################


clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
TCGAanalyze_survival(clin.gbm,
                     "gender",
                     main = "TCGA Set\n GBM",height = 10, width=10)


library(TCGAbiolinks)
# Survival Analysis SA

clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
dataBRCAcomplete <- log2(BRCA_rnaseqv2)

tokenStop<- 1

tabSurvKMcomplete <- NULL

for( i in 1: round(nrow(dataBRCAcomplete)/100)){
  message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataBRCAcomplete,
                                    Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]



