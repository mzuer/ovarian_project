require(TCGAbiolinks)
require(SummarizedExperiment)
require(recount)
# Rscript download_tcga_brca_recount2.R

startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGA_BRCA_RECOUNT2"
dir.create(outFolder)


tcga_purity_thresh <- 0.6 # cf Lucchetta et al. 2019


############################################## 
############################################## prepare TCGA data
############################################## 


breast_rec2_tcga <- TCGAquery_recount2(project="tcga", tissue = "breast")
breast_rec2_tcga_scaled <- scale_counts(breast_rec2_tcga$tcga_breast)
tcga_counts_raw_all <- assays(breast_rec2_tcga_scaled)$counts


stopifnot(breast_rec2_tcga_scaled@colData[,"cgc_case_primary_site"] == "Breast")
annot_cols <- c(grep("^cgc_", names(breast_rec2_tcga_scaled@colData)), grep("^gdc_", names(breast_rec2_tcga_scaled@colData)))
tcgacols <- names(breast_rec2_tcga_scaled@colData)[annot_cols]
tcga_sampleAnnot <- data.frame(breast_rec2_tcga_scaled@colData[, tcgacols])

table(tcga_sampleAnnot[,"cgc_case_histological_diagnosis"])
# Infiltrating Carcinoma NOS    Infiltrating Ductal Carcinoma   Infiltrating Lobular Carcinoma 
# 1                              896                              209 
# Medullary Carcinoma            Metaplastic Carcinoma Mixed Histology (please specify) 
# 8                               12                               37 
# Mucinous Carcinoma                   Other  specify 
# 18                               54 
table(tcga_sampleAnnot[,"cgc_case_pathologic_stage"])# 
# Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III Stage IIIA Stage IIIB 
# 105         94          6          6        409        296          2        173         33 
# Stage IIIC   Stage IV  Stage Tis    Stage X 
# 69         22          1         14 
table(tcga_sampleAnnot[,"cgc_file_disease_type"])# 
# Breast Invasive Carcinoma 
# 1246 
table(tcga_sampleAnnot[,"cgc_sample_sample_type"])# 
# Metastatic       Primary Tumor Solid Tissue Normal 
# 7                1127                 112 
table(tcga_sampleAnnot[,"cgc_file_experimental_strategy"] == "RNA-Seq")# 

#### FILTER HERE -> TAKE ONLY PRIMARY TUMOR !!!
tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$cgc_sample_sample_type == "Primary Tumor",]
stopifnot(nrow(tcga_sampleAnnot) > 0)

# retrieve ER+ ER- status
# as indicated here: https://www.biostars.org/p/279048/#279069
# download from https://portal.gdc.cancer.gov/legacy-archive/search/f, 9/9/21
full_annot_dt <- read.delim("nationwidechildrens.org_clinical_patient_brca.txt", header=TRUE, stringsAsFactors = FALSE)
head(full_annot_dt$bcr_patient_barcode)
head(full_annot_dt$er_status_by_ihc)
head(full_annot_dt$nte_er_status)
head(full_annot_dt$bcr_patient_uuid)
# [1] "bcr_patient_uuid"                     "CDE_ID:"                             
# [3] "6E7D5EC6-A469-467C-B748-237353C23416" "55262FCB-1B01-4480-B322-36570430C917"
# [5] "427D0648-3F77-4FFC-B52C-89855426D647" "C31900A4-5DCD-4022-97AC-638E86E889E4"
head(full_annot_dt$bcr_patient_barcode)
# [1] "bcr_patient_barcode" "CDE_ID:2673794"      "TCGA-3C-AAAU"        "TCGA-3C-AALI"       
# [5] "TCGA-3C-AALJ"        "TCGA-3C-AALK"       
head(full_annot_dt$patient_id)
# [1] "patient_id" "CDE_ID:"    "AAAU"       "AALI"       "AALJ"       "AALK"      

tcga_sampleAnnot$patient_barcode <- substr(start=1, stop=12, tcga_sampleAnnot$cgc_sample_id)
stopifnot(tcga_sampleAnnot$patient_barcode %in% full_annot_dt$bcr_patient_barcode)
stopifnot(!duplicated(full_annot_dt$bcr_patient_barcode))
sampleAnnot <- setNames(full_annot_dt$er_status_by_ihc, full_annot_dt$bcr_patient_barcode)
tcga_sampleAnnot$ER_status <- sampleAnnot[as.character(tcga_sampleAnnot$patient_barcode)]
stopifnot(!is.na(tcga_sampleAnnot$ER_status))
table(tcga_sampleAnnot$ER_status)
# [Not Evaluated]   Indeterminate        Negative        Positive 
# 48               2             251             826 

stopifnot(rownames(tcga_sampleAnnot) %in% colnames(tcga_counts_raw_all))

ERpos_samples <- rownames(tcga_sampleAnnot)[tcga_sampleAnnot$ER_status == "Positive"] 
nERpos <- length(ERpos_samples)
cat(nERpos, "\n")
# 826
ERneg_samples <- rownames(tcga_sampleAnnot)[tcga_sampleAnnot$ER_status == "Negative"] 
nERneg <- length(ERneg_samples)
cat(nERneg, "\n")
# 251

tcga_counts_raw_all <- tcga_counts_raw_all[,c(ERpos_samples, ERneg_samples)]
stopifnot(!is.na(tcga_counts_raw_all))  
  
tcga_sampleAnnot <- tcga_sampleAnnot[c(ERpos_samples, ERneg_samples),]


### PURITY FILTER FOR THE TUMOR SAMPLES
# taken from Lucchetta1 et al. 2019 https://bmccancer.biomedcentral.com/track/pdf/10.1186/s12885-019-5965-x
# https://github.com/ELELAB/LUAD_LUSC_TCGA_comparison/blob/master/6-recount/unifiedLUAD_18062018.R
# purityinfo.R.luad<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)

tcga_sampleAnnot$labs_for_purity <- tcga_sampleAnnot$cgc_sample_id
purityinfo_brca <-TCGAtumor_purity(tcga_sampleAnnot$labs_for_purity, 0, 0, 0, 0, tcga_purity_thresh)
#  $pure_barcodes attribute as a vector of pure samples and
# $filtered attribute as filtered samples with no purity info
stopifnot(is.null(purityinfo_brca$filtered)) # ensure all samples had purity info
stopifnot(purityinfo_brca$pure_barcodes %in% tcga_sampleAnnot$labs_for_purity)
nPureSamp <- length(purityinfo_brca$pure_barcodes)
cat(paste0("... samples to keep after purity filter (", tcga_purity_thresh, ") = ",
           round(sum(tcga_sampleAnnot$labs_for_purity%in% purityinfo_brca$pure_barcodes )/nrow(tcga_sampleAnnot) * 100, 2)), " %\n")
# ... samples to keep after purity filter (0.6) = 84.83  %
tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$labs_for_purity%in% purityinfo_brca$pure_barcodes,]
stopifnot(rownames(tcga_sampleAnnot) %in% colnames(tcga_counts_raw_all))
# keep only the ones that pass purity filter
tcga_counts_raw_pF <- tcga_counts_raw_all[, colnames(tcga_counts_raw_all) %in% rownames(tcga_sampleAnnot)]
stopifnot(ncol(tcga_counts_raw_pF) == nPureSamp)
stopifnot(nrow(tcga_sampleAnnot) == nPureSamp)
stopifnot(colnames(tcga_counts_raw_pF) %in% c(ERpos_samples, ERneg_samples))


# for tcga -> takes the one filtered
tcga_counts_all <- tcga_counts_raw_pF
stopifnot(tcga_counts_all >= 0)

# take the first symbol of the list... don't know how to do better...
tcga_symbs <- sapply(data.frame(breast_rec2_tcga_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
tcga_ids <- as.character(data.frame(breast_rec2_tcga_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(tcga_counts_all) == tcga_ids)
stopifnot(length(tcga_ids) == length(tcga_symbs))
stopifnot(!duplicated(tcga_ids))

gene_dt = data.frame(
  geneSymb=tcga_symbs,
  geneID=tcga_ids,
  stringsAsFactors = FALSE
)

stopifnot(rownames(tcga_counts_all) == gene_dt$geneID)





outFile <- file.path(outFolder, paste0("all_counts_onlyPF_", tcga_purity_thresh, ".Rdata"))
save(tcga_counts_all, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "gene_dt.Rdata")
save(gene_dt, file=outFile)
cat(paste0("... written ", outFile, "\n"))




##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))

stop("---ok\n")
