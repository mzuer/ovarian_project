require(TCGAbiolinks)
require(SummarizedExperiment)
require(recount)
# Rscript download_tcga_gtex_recount2.R

# https://github.com/LieberInstitute/recount-brain/blob/master/SupplementaryTable2.csv


startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGA_GTEX_RECOUNT2"
dir.create(outFolder)

dataset_rowMeans_filter <- 0

tcga_purity_thresh <- 0.6
# taken from Lucchetta1 et al. 2019 https://bmccancer.biomedcentral.com/track/pdf/10.1186/s12885-019-5965-x
# https://github.com/ELELAB/LUAD_LUSC_TCGA_comparison/blob/master/6-recount/unifiedLUAD_18062018.R
# purityinfo.R.luad<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)

ovary_rec2_gtex <- TCGAquery_recount2(project="gtex", tissue = "ovary")
ovary_rec2_gtex_scaled <- scale_counts(ovary_rec2_gtex$gtex_ovary)

ovary_rec2_tcga <- TCGAquery_recount2(project="tcga", tissue = "ovary")
ovary_rec2_tcga_scaled <- scale_counts(ovary_rec2_tcga$tcga_ovary)

head(ovary_rec2_tcga_scaled@colData[,"cgc_sample_id"])
table(ovary_rec2_tcga_scaled@colData[,"cgc_case_primary_site"])
# Ovary: 430
table(ovary_rec2_tcga_scaled@colData[,"cgc_case_clinical_stage"])
# Stage IC  Stage IIA  Stage IIB  Stage IIC Stage IIIA Stage IIIB Stage IIIC   Stage IV 
# 1          3          5         18          7         18        311         64 
table(ovary_rec2_tcga_scaled@colData[,"cgc_sample_sample_type"])
# Primary Tumor Recurrent Tumor 
# 422               8 

table(ovary_rec2_tcga_scaled@colData[,"cgc_case_tumor_status"])
# TUMOR FREE WITH TUMOR 
# 96        283 

tcgacols <- names(ovary_rec2_tcga_scaled@colData)[grepl("cgc_", names(ovary_rec2_tcga_scaled@colData))]
tcga_sampleAnnot <- data.frame(ovary_rec2_tcga_scaled@colData[, tcgacols])

# https://github.com/LieberInstitute/recount-brain/blob/master/SupplementaryTable2.csv
annot_dt <- read.delim("GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt", header=T)
# ovary_rec2_gtex_scaled <- add_predictions(ovary_rec2_gtex_scaled)
stopifnot(as.character(ovary_rec2_gtex_scaled@colData[, "smtsd"]) =="Ovary") # Anatomic Stie of tissue 2
gtex_sampleAnnot_tmp <- data.frame(ovary_rec2_gtex_scaled@colData[, c("sampid", "smtsd")])
gtex_sampleAnnot_tmp$SUBJID <- gsub("(.+-.+?)-.+", "\\1",gtex_sampleAnnot_tmp$sampid)
gtex_sampleAnnot_tmp2 <- merge(gtex_sampleAnnot_tmp, annot_dt, by="SUBJID", all.x=TRUE, all.y=FALSE)
stopifnot(!any(is.na(gtex_sampleAnnot_tmp2$AGE)))
stopifnot(nrow(gtex_sampleAnnot_tmp2) == nrow(gtex_sampleAnnot_tmp))
gtex_sampleAnnot <- gtex_sampleAnnot_tmp2[match(gtex_sampleAnnot_tmp$sampid, gtex_sampleAnnot_tmp2$sampid),]
stopifnot(gtex_sampleAnnot$sampid == gtex_sampleAnnot_tmp$sampid)
rownames(gtex_sampleAnnot) <- rownames(gtex_sampleAnnot_tmp)

## Extract counts and filter out lowly expressed genes
gtex_counts_all <- assays(ovary_rec2_gtex_scaled)$counts
stopifnot(gtex_counts_all >= 0)
gtex_filter <- rowMeans(gtex_counts_all) >= dataset_rowMeans_filter
sum(gtex_filter)
# 37644

tcga_counts_all <- assays(ovary_rec2_tcga_scaled)$counts
stopifnot(tcga_counts_all >= 0)
tcga_filter <- rowMeans(tcga_counts_all) >= dataset_rowMeans_filter
sum(tcga_filter)
# 38670

tcga_symbs <- sapply(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
tcga_ids <- as.character(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(tcga_counts_all) == tcga_ids)
stopifnot(length(tcga_ids) == length(tcga_symbs))
stopifnot(length(tcga_ids) == length(tcga_filter))   
stopifnot(!duplicated(tcga_ids))
tcga_g2s <- setNames(tcga_symbs[tcga_filter], tcga_ids[tcga_filter]) ##### keep only the IDs that have passed the count filter
any(duplicated(as.character(tcga_g2s)))
tcga_g2s_filt <- tcga_g2s[!duplicated(tcga_g2s)]
tcga_g2s_filt <- tcga_g2s_filt[!is.na(tcga_g2s_filt)]


gtex_symbs <- sapply(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
gtex_ids <- as.character(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(gtex_counts_all) == gtex_ids)
stopifnot(length(gtex_ids) == length(gtex_symbs))
stopifnot(length(gtex_ids) == length(gtex_filter))
stopifnot(!duplicated(gtex_ids))
gtex_g2s <- setNames(gtex_symbs[gtex_filter], gtex_ids[gtex_filter]) ##### keep only the IDs that have passed the count filter
any(duplicated(as.character(gtex_g2s)))
gtex_g2s_filt <- gtex_g2s[!duplicated(gtex_g2s)]
gtex_g2s_filt <- gtex_g2s_filt[!is.na(gtex_g2s_filt)]

# I need to do like this because I can have a given gene symbol 
# mapped to  ENSG in tcga and to a different ENSG in gtex
# so with this table I ensure exact gene ID gene Symbol matching
common_dt = data.frame(
  geneSymb=c(as.character(tcga_g2s_filt), as.character(gtex_g2s_filt)),
  geneID=c(names(tcga_g2s_filt),names(gtex_g2s_filt)),stringsAsFactors = FALSE
)
common_dt$filtid <- paste0(common_dt$geneID, "_", common_dt$geneSymb)
id_to_keep <- names(table(common_dt$filtid))[which(as.numeric(table(common_dt$filtid)) == 2)]
intersect_dt <- common_dt[common_dt$filtid %in% id_to_keep,]
stopifnot(nrow(intersect_dt) == nrow(unique(intersect_dt)) * 2)
intersect_dt <- unique(intersect_dt)

common_symbs <- as.character(intersect_dt$geneSymb)
common_ids <- as.character(intersect_dt$geneID)
stopifnot(length(common_symbs) == length(common_ids))
stopifnot(!duplicated(common_ids))
stopifnot(!duplicated(common_symbs))
length(common_ids)
# 20698
g2s <- setNames(as.character(intersect_dt$geneSymb), as.character(intersect_dt$geneID))

new_rownames1 <- as.character(gtex_g2s_filt[common_ids])
new_rownames2 <- as.character(tcga_g2s_filt[common_ids])
stopifnot(length(new_rownames2) == length(new_rownames1))
stopifnot(new_rownames2==new_rownames1)
stopifnot(new_rownames2==g2s[common_ids])


sub_gtex_counts <- gtex_counts_all[common_ids,]
sub_tcga_counts <- tcga_counts_all[common_ids,]
stopifnot(rownames(sub_gtex_counts) == rownames(sub_tcga_counts))

stopifnot(colnames(sub_gtex_counts) == rownames(gtex_sampleAnnot))
stopifnot(!duplicated(gtex_sampleAnnot$sampid))
colnames(sub_gtex_counts) <- gtex_sampleAnnot$sampid

stopifnot(colnames(sub_tcga_counts) == rownames(tcga_sampleAnnot))
stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
colnames(sub_tcga_counts) <- tcga_sampleAnnot$cgc_sample_id

all_counts <- cbind(sub_gtex_counts, sub_tcga_counts)
stopifnot(rownames(all_counts)==names(g2s[common_ids]))
rownames(all_counts) <- g2s[common_ids]
dim(all_counts)
stopifnot(length(new_rownames1) == nrow(all_counts))
stopifnot(ncol(sub_gtex_counts) + ncol(sub_tcga_counts) == ncol(all_counts))

outFile <- file.path(outFolder, paste0("all_counts_dsCountsGe_",dataset_rowMeans_filter,".Rdata"))
save(all_counts, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "tcga_sampleAnnot.Rdata")
save(tcga_sampleAnnot, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "gtex_sampleAnnot.RData")
save(gtex_sampleAnnot, file=outFile)
cat(paste0("... written ", outFile, "\n"))



##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))



require(TCGAbiolinks)

#the following lines are only to query the TCGA data 
query_ov<- GDCquery(project = "TCGA-OV",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts")


samplesOv <- getResults(query_ov,cols=c("cases"))


dataSmTP_ov <- TCGAquery_SampleTypes(barcode = samplesOv,
                                       typesample = "TP")

dataSmNT_ov <- TCGAquery_SampleTypes(barcode = samplesOv,
                                       typesample = "NT")



query.luad2 <- GDCquery(project = "TCGA-LUAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts",
                        barcode = c(dataSmTP.luad, dataSmNT.luad))


GDCdownload(query=query.luad2)

dataPrep1.luad <- GDCprepare(query = query.luad2, 
                             save = TRUE )

purityinfo_ov <-TCGAtumor_purity(xx, 0, 0, 0, 0, 0.6)


#####getting samples with more than 60% tumor purity and removing discordant samples

tcga.barcodes<-c(dataSmTP.luad, dataSmNT.luad)
dataSmTP.luad.pure.R<-purityinfo.R.luad$pure
tcga.pure_barcodes<-c(dataSmTP.luad.pure.R, dataSmNT.luad)


uuid.tcga<-barcodeToUUID(tcga.pure_barcodes)$case_id
uuid.tcga.normal<-barcodeToUUID(dataSmNT.luad)$case_id
uuid.tcga.cancer<-barcodeToUUID(dataSmTP.luad.pure.R)$case_id

