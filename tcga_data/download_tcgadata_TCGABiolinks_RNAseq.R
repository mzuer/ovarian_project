require(TCGAbiolinks)
require(SummarizedExperiment)
# Rscript download_tcgadata_TCGABiolinks_RNAseq.R

# Ovarian Serous Cystadenocarcinoma (TCGA-OV)	TCGA-OV

startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGADATA_TCGABIOLINKS_RNASEQ"
dir.create(outFolder)

buildData <- TRUE

if(buildData) {
  # https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#legacy_archive_examples
query_rna_hg38_ <- GDCquery(project = "TCGA-OV", 
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts",
                      experimental.strategy = "RNA-Seq")


# https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#legacy_archive_examples
query_rna_hg19_ <- GDCquery(project = "TCGA-OV", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results",
                      experimental.strategy = "RNA-Seq")



cat("... download for hg38\n")
GDCdownload(query_rna_hg38_)
ov_RNAseqSE_hg38 <- GDCprepare(query_rna_hg38_)
outFile <- file.path(outFolder, "TCGAbiolinks_ov_RNAseqSE_hg38.Rdata")
save(ov_RNAseqSE_hg38, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


cat("... download for hg19\n")
GDCdownload(query_rna_hg19_)
ov_RNAseqSE_hg19 <- GDCprepare(query_rna_hg19_)
outFile <- file.path(outFolder, "TCGAbiolinks_ov_RNAseqSE_hg19.Rdata")
save(ov_RNAseqSE_hg19, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFolder, "TCGAbiolinks_ov_RNAseqSE_hg19.Rdata")
  ov_RNAseqSE_hg19 <- get(load(outFile))
  outFile <- file.path(outFolder, "TCGAbiolinks_ov_RNAseqSE_hg38.Rdata")
  ov_RNAseqSE_hg38 <- get(load(outFile))
}


ov_rnaseq_mat_h19 <- assay(ov_RNAseqSE_hg19,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
dim(ov_RNAseqSE_hg19)

ov_rnaseq_mat_h38 <- assay(ov_RNAseqSE_hg38,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
dim(ov_RNAseqSE_hg38)

# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006701
# TCGAquery_recount2 function to download tumor and normal ovary samples 
# from the Recount2 platform as Ranged Summarized Experiment (RSE) objects.
# Recount2 resource contains reads, some of them soft-clipped, aligned to Gencode version
# 25 hg38 using the splice-aware Rail-RNA aligner. Moreover, the RSE shows coverage counts 
# instead of standard read count matrices. Since most methods are adapted to read count matrices,
# there are some highly recommended transformations to perform before commencing with DEA
# https://support.bioconductor.org/p/9136299/#9136308
# recount has a function scale_counts to create a count matrix appropriate for running DESeq2,
# edgeR, limma-voom or other count based packages. From their vignette, in the code just before running DESeq2:
# https://bioconductor.org/packages/release/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html
ovary_rec2_gtex <- TCGAquery_recount2(project="gtex", tissue = "ovary")
ovary_rec2_gtex_scaled <- scale_counts(ovary_rec2_gtex$gtex_ovary)

ovary_rec2_tcga <- TCGAquery_recount2(project="tcga", tissue = "ovary")
ovary_rec2_tcga_scaled <- scale_counts(ovary_rec2_tcga$tcga_ovary)


## Extract counts and filter out lowly expressed geens
gtex_counts <- assays(ovary_rec2_gtex_scaled)$counts
gtex_filter <- rowMeans(gtex_counts) > 0.5
gtex_counts <- gtex_counts[gtex_filter,]
dim(gtex_counts)
gtex_counts[1:5,1:5]
# 37644   108
tcga_counts <- assays(ovary_rec2_tcga_scaled)$counts
tcga_filter <- rowMeans(tcga_counts) > 0.5
tcga_counts <- tcga_counts[tcga_filter,]
dim(tcga_counts)
tcga_counts[1:5,1:5]
# 38670   430

common_genes <- intersect(rownames(gtex_counts), rownames(tcga_counts))
gtex_counts <- gtex_counts[common_genes,]
tcga_counts <- tcga_counts[common_genes,]
stopifnot(dim(gtex_counts)[1] == dim(tcga_counts)[1])

all_counts <- cbind(gtex_counts, tcga_counts)
stopifnot(rownames(all_counts) == common_genes)
# [1] 35120   538

library("limma")
library("edgeR")

## Build DGEList object
dge <- DGEList(counts = all_counts)

## Calculate normalization factors
dge <- calcNormFactors(dge)

## Specify our design matrix
samples_groups <- c(rep("gtex", ncol(gtex_counts)), rep("tcga", ncol(tcga_counts))) 
my_group_design <- factor(samples_groups, levels = c("gtex", "tcga"))
my_design <- model.matrix( ~ my_group_design)


## Explore the data
# plotMDS(dge), labels = substr(colData(rse_gene_scaled)$prenatal, 1, 2))
# plotMDS(dge, labels = substr(colData(rse_gene_scaled)$sex, 1, 1))
# tapply(
#   colData(rse_gene_scaled)$rin, colData(rse_gene_scaled)$prenatal,
#   summary
# )

## Run voom
v <- voom(dge, my_design, plot = TRUE)

## Run remaining parts of the DE analysis
fit <- lmFit(v, my_design)
fit <- eBayes(fit)

## Run remaining parts of the DE analysis
fit <- lmFit(v, my_design)
fit <- eBayes(fit)
## Visually explore DE results
limma::plotMA(fit)
limma::volcanoplot(fit)
## Extract data from limma-voom results
top <- topTable(fit,
                number = Inf, sort.by = "none",
                coef = ncol(v$design))




 # is possible to use the function scale_counts from the recount package. 
# After that, we merged the two prepared gene count matrices, normalized for GC-content and applied the quantile filtering with a 25% cut-off. The data were then loaded into the TCGAanalyze_DEA function for comparison of normal samples versus cancer samples using the limma-voom pipeline. 

sys.exit(0)
# https://bioconductor.org/packages/release/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html
# library("clusterProfiler")
# library("org.Hs.eg.db")
# 
# ## Remember that limma_res had ENSEMBL IDs for the genes
# head(rownames(limma_res))
# 
# ## [1] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938"
# 
# ## Perform enrichment analysis for Biological Process (BP)
# ## Note that the argument is keytype instead of keyType in Bioconductor 3.5
# enrich_go <- enrichGO(
#   gene = rownames(limma_res)[limma_res$padj < 0.001],
#   OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
#   pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
#   universe = rownames(limma_res)
# )
# 
# ## Visualize enrichment results
# dotplot(enrich_go, font.size = 7)


BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")




ov_rnaseq_hg38 <- getResults(query_rna_hg38_)
colnames(ov_rnaseq_hg38)
head(ov_rnaseq_hg38$sample_type)
table(ov_rnaseq_hg38$sample_type)
# Primary Tumor Recurrent Tumor 
# 374               5 

ov_rnaseq_hg19 <- getResults(query_rna_hg19_)
colnames(ov_rnaseq_hg19)
head(ov_rnaseq_hg19$sample_type)
table(ov_rnaseq_hg19$sample_type)
# normalized_results OR results -> same result
# Primary Tumor Recurrent Tumor 
# 304       5 

outFile <- file.path(outFolder, "TCGAbiolinks_ov_rnaseq_hg38.Rdata")
save(ov_rnaseq_hg38, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, "TCGAbiolinks_ov_rnaseq_hg19.Rdata")
save(ov_rnaseq_hg19, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

} else {
outFile <- file.path(outFolder, "TCGAbiolinks_ov_rnaseq_hg38.Rdata")
ov_rnaseq_hg38 <- get(load(outFile))

outFile <- file.path(outFolder, "TCGAbiolinks_ov_rnaseq_hg19.Rdata")
ov_rnaseq_hg19 <- get(load(outFile))
}
names(ov_rnaseq_hg19)
# [1] "id"                    "data_format"           "access"               
# [4] "cases"                 "file_name"             "submitter_id"         
# [7] "data_category"         "type"                  "file_size"            
# [10] "platform"              "state_comment"         "tags"                 
# [13] "updated_datetime"      "md5sum"                "file_id"              
# [16] "data_type"             "state"                 "experimental_strategy"
# [19] "file_state"            "version"               "data_release"         
# [22] "project"               "center_id"             "center_center_type"   
# [25] "center_code"           "center_name"           "center_namespace"     
# [28] "center_short_name"     "sample_type"           "is_ffpe"              
# [31] "cases.submitter_id"    "sample.submitter_id"  
GDCdownload(query_rna_hg38_)

BRCARnaseqSE <- GDCprepare(query)






out_dt <- assay(ov_rnaseq_hg19, "raw_counts")









ov.exp_27_clinic <- ov_rnaseq_hg19@colData

sample_type <- sapply(ov.exp_27_clinic$barcode, substr, 14,15)
table(sample_type)
# 01  02  11 
# 583  18  12
# 12 normal tissue
# TCGA-02-0001-01C-01D-0192-01
### 0001 is the participant then 01 indicates that it is tumor 
# (tumor from 01-09; 10-19 normal and control 20-29) (C = the third vial)

table(ov.exp_27_clinic$tumor_stage)
# not reported 
# 594 
table(ov.exp_27_clinic$tumor_grade)
# not reported 
# 594 
table(ov.exp_27_clinic$classification_of_tumor)
# not reported 
# 594 
table(ov.exp_27_clinic$tissue_type)
# Not Reported 
# 613 
# View(ov.exp_27_clinic)
table(ov.exp_27_clinic$primary_diagnosis)
# Serous cystadenocarcinoma, NOS Serous surface papillary carcinoma 
# 590                                  4 
table(ov.exp_27_clinic$figo_stage)
# Stage IA   Stage IB   Stage IC  Stage IIA  Stage IIB  Stage IIC Stage IIIA Stage IIIB 
# 4          3         10          3          4         23          8         24 
# Stage IIIC   Stage IV 
# 424         86 
#FIGO: the International Federation of Gynecology and Obstetrics.

### the variables I keep for output dt
figo_stages <- setNames(ov.exp_27_clinic$figo_stage, rownames(ov.exp_27_clinic))
sample_types <- setNames(ov.exp_27_clinic$sample_type, rownames(ov.exp_27_clinic))
diagnosis <- setNames(ov.exp_27_clinic$primary_diagnosis, rownames(ov.exp_27_clinic))


#***********************************
# output the mapping of the probes 
ov.exp_27_probes_dt <- data.frame(rowData(ov.exp_27))
stopifnot(rownames(ov.exp_27_probes_dt) == ov.exp_27_probes_dt$Composite.Element.REF)
rownames(ov.exp_27_probes_dt) <- NULL

outFile <- file.path(outFolder, "TCGAbiolinks_OV_DNAmet27_hg38_probesDT.txt")
write.table(ov.exp_27_probes_dt, file =outFile, col.names = TRUE, row.names=FALSE, sep="\n", quote=F)
cat(paste0("... written: ", outFile, "\n"))



#***********************************
# output the methylation values 
out_dt <- t(assay(ov.exp_27))
dim(out_dt)
out_dt[1:5,1:5]
# 613x27578
stopifnot(rownames(out_dt) %in% rownames(ov.exp_27_clinic))
stopifnot(rownames(out_dt) %in% names(figo_stages))
stopifnot(rownames(out_dt) %in% names(diagnosis))
stopifnot(rownames(out_dt) %in% names(sample_types))

# remove the full NA
sum(apply(out_dt, 2, function(x) all(is.na(x))))
# [1] 2597
to_keep <- which(apply(out_dt, 2, function(x) !all(is.na(x))))
stopifnot(length(to_keep) + sum(apply(out_dt, 2, function(x) all(is.na(x)))) == ncol(out_dt))
out_dt_filtered <- out_dt[,to_keep]

probe_mads <- apply(out_dt_filtered, 2, mad, na.rm=TRUE)
head(sort(probe_mads, decreasing = TRUE))
# cg17965019 cg14159672 cg22881914 cg04797323 cg23495733 cg10516359 
# 0.5941993  0.5613968  0.5511603  0.5432971  0.5390746  0.5218498 
# compare to what obtained in python:
# mad_probes.sort_values(ascending=False)
# Out[199]: 
#   cg17965019    0.398750
# cg14159672    0.381125
# cg22881914    0.372497
# cg04797323    0.370243
# cg23495733    0.364115
# -> looks good

col_order <- names(sort(probe_mads, decreasing=TRUE))
stopifnot(setequal(col_order, colnames(out_dt_filtered)))
out_dt_filtAndSort <- out_dt_filtered[,col_order]
stopifnot(setequal(colnames(out_dt_filtAndSort), colnames(out_dt_filtered)))
stopifnot(colnames(out_dt_filtAndSort) == col_order)
out_dt_filtAndSort[1:5,1:5]

patient_dt <- data.frame(
  sampleID = rownames(out_dt_filtAndSort),
  sample_types = as.character(sample_types[rownames(out_dt_filtAndSort)]),
  diagnosis = as.character(diagnosis[rownames(out_dt_filtAndSort)]),
  figo_stages = as.character(figo_stages[rownames(out_dt_filtAndSort)]),
  stringsAsFactors = FALSE
)

stopifnot(patient_dt$sampleID == rownames(out_dt_filtAndSort))

rownames(out_dt_filtAndSort) <- NULL

ov_met27_data_dt <- cbind(patient_dt, out_dt_filtAndSort)
ov_met27_data_dt[1:5,1:5]

outFile <- file.path(outFolder, "TCGAbiolinks_OV_DNAmet27_hg38_metDT.txt")
write.table(ov_met27_data_dt, file =outFile, col.names = TRUE, row.names=FALSE, sep="\n", quote=F)
cat(paste0("... written: ", outFile, "\n"))


##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))

