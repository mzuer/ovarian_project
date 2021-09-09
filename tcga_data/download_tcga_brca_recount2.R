require(TCGAbiolinks)
require(SummarizedExperiment)
require(recount)
require(matrixStats)
# Rscript download_tcga_brca_recount2.R

startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGA_BRCA_RECOUNT2"
dir.create(outFolder)


tcga_purity_thresh <- 0.6 # cf Lucchetta et al. 2019

# ER status
# https://www.biostars.org/p/279048/
# duplicated samples
# https://www.biostars.org/p/311017/

compute_meds <- TRUE

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
stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_counts_raw_all))
tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$cgc_sample_sample_type == "Primary Tumor",]
stopifnot(nrow(tcga_sampleAnnot) > 0)
tcga_counts_raw_all <- tcga_counts_raw_all[,rownames(tcga_sampleAnnot)]
stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_counts_raw_all))

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


### DO THIS BEFORE COMPUTING THE MEDIAN !
# take the first symbol of the list... don't know how to do better...
tcga_symbs <- sapply(data.frame(breast_rec2_tcga_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
tcga_ids <- as.character(data.frame(breast_rec2_tcga_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(tcga_counts_raw_all) == tcga_ids)
stopifnot(length(tcga_ids) == length(tcga_symbs))
stopifnot(!duplicated(tcga_ids))

gene_dt = data.frame(
  geneSymb=tcga_symbs,
  geneID=tcga_ids,
  stringsAsFactors = FALSE
)
# I think lot of non coding no gene symb
stopifnot(rownames(tcga_counts_raw_all) == gene_dt$geneID)
cat(paste0("... filter NA gene symbols: ", sum(is.na(gene_dt$geneSymb)), "\n"))
# 32511
gene_dt <- gene_dt[!is.na(gene_dt$geneSymb),]
stopifnot(gene_dt$geneID %in% rownames(tcga_counts_raw_all))
tcga_counts_raw_all <- tcga_counts_raw_all[gene_dt$geneID,]
stopifnot(!is.na(tcga_counts_raw_all))
dim(tcga_counts_raw_all)
#[1] 25526  1127

### subset ER+/- before the median
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

cat("any duplicated(tcga_sampleAnnot$cgc_sample_id)", "\n")
cat(any(duplicated(tcga_sampleAnnot$cgc_sample_id)), "\n")
# stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
# this is not TRUE ! I have duplicated samples
# see discussion here: https://www.biostars.org/p/311017/
# when duplicated, I take the one with highest median value
nGenesCheck <- nrow(tcga_counts_raw_all)
withdup_count_dt <- as.data.frame(t(tcga_counts_raw_all))
stopifnot(rownames(withdup_count_dt) == rownames(tcga_sampleAnnot))
withdup_count_dt$samp_id <- tcga_sampleAnnot$cgc_sample_id
sum(duplicated(withdup_count_dt$samp_id))
# [1] 23
nSampsCheck <- length(unique(withdup_count_dt$samp_id))
dup_samples <- unique(withdup_count_dt$samp_id[duplicated(withdup_count_dt$samp_id)])

if(compute_meds) {
  no_dup_data <- withdup_count_dt[!withdup_count_dt$samp_id %in% dup_samples,]
  
  dup_data <- withdup_count_dt[withdup_count_dt$samp_id %in% dup_samples,]
  stopifnot(length(unique(dup_data$samp_id)) == length(dup_samples))
  ### !!! TODO
  cat(paste0("compute median dup samples\n"))
  dup_data_duprm <- do.call(rbind, by(dup_data, dup_data$samp_id, function(sub_dt) {
    # each row = a sample, keep the sample that has highest median value
    tmp_dt <- data.frame(sub_dt)
    tmp_dt$samp_id <- NULL
    tmp_dt <- log10(tmp_dt+1)
    stopifnot(!is.na(tmp_dt))
    stopifnot(ncol(tmp_dt) == ncol(dup_data)-1)
    stopifnot(nrow(tmp_dt) == nrow(sub_dt))
    save(tmp_dt, file="tmp_dt.Rdata")
    cat(paste0("written\n"))
    all_meds <- apply(tmp_dt, 1, median)
    all_meds <- rowMedians(as.matrix(tmp_dt))
    stopifnot(length(all_meds) == nrow(sub_dt))
    stopifnot(!is.na(all_meds))
    to_keep <- which.max(all_meds)
    stopifnot(length(to_keep) == 1)
    # need this hack otherwise rownames will be samp_id
    out_dt <- sub_dt[to_keep,, drop=FALSE]
    out_dt$full_id <- rownames(sub_dt)[to_keep]
    out_dt
  }))
  # need this hack otherwise rownames will be samp_id (because of by())
  stopifnot(!duplicated(dup_data_duprm$full_id))
  rownames(dup_data_duprm) <- dup_data_duprm$full_id
  dup_data_duprm$full_id <- NULL
  stopifnot(nrow(dup_data_duprm) == length(dup_samples))
  stopifnot(!duplicated(dup_data_duprm$samp_id))
  cat(paste0("... ok meds\n"))
  
  stopifnot(colnames(dup_data_duprm) == colnames(no_dup_data))
  
  new_tcga_counts_raw_all <- rbind(dup_data_duprm, no_dup_data)
  new_tcga_counts_raw_all$samp_id <- NULL
  new_tcga_counts_raw_all <- t(new_tcga_counts_raw_all)
  stopifnot(rownames(new_tcga_counts_raw_all) == rownames(tcga_counts_raw_all))
  stopifnot(colnames(new_tcga_counts_raw_all) %in% colnames(tcga_counts_raw_all))
  outFile <- file.path(outFolder,"new_tcga_counts_raw_all.Rdata" )
  save(new_tcga_counts_raw_all, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder,"new_tcga_counts_raw_all.Rdata" )
  new_tcga_counts_raw_all <- get(load(outFile))
}


## update the count data
tcga_counts_raw_all <- new_tcga_counts_raw_all

stopifnot(nrow(tcga_counts_raw_all) == nGenesCheck) 
stopifnot(ncol(tcga_counts_raw_all) == nSampsCheck) 

## update annot data
stopifnot(colnames(tcga_counts_raw_all)%in% rownames(tcga_sampleAnnot))
tcga_sampleAnnot <- tcga_sampleAnnot[colnames(tcga_counts_raw_all),] 
stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
stopifnot(colnames(tcga_counts_raw_all) == rownames(tcga_sampleAnnot))

ERpos_samples <- rownames(tcga_sampleAnnot)[tcga_sampleAnnot$ER_status == "Positive"] 
nERpos <- length(ERpos_samples)
cat(nERpos, "\n")
ERneg_samples <- rownames(tcga_sampleAnnot)[tcga_sampleAnnot$ER_status == "Negative"] 
nERneg <- length(ERneg_samples)
cat(nERneg, "\n")
stopifnot(c(ERpos_samples, ERneg_samples) %in% colnames(tcga_counts_raw_all))
stopifnot(colnames(tcga_counts_raw_all) %in% c(ERpos_samples, ERneg_samples))
tcga_counts_raw_all <- tcga_counts_raw_all[,c(ERpos_samples, ERneg_samples)]

stopifnot(c(ERpos_samples, ERneg_samples) %in% rownames(tcga_sampleAnnot))
stopifnot(rownames(tcga_sampleAnnot) %in% c(ERpos_samples, ERneg_samples))
tcga_sampleAnnot <- tcga_sampleAnnot[c(ERpos_samples, ERneg_samples),]

stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_counts_raw_all))

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

ERpos_samples <- ERpos_samples[ERpos_samples %in% rownames(tcga_sampleAnnot)]
nERpos <- length(ERpos_samples)
ERneg_samples <- ERneg_samples[ERneg_samples %in% rownames(tcga_sampleAnnot)]
nERneg <- length(ERneg_samples)
stopifnot(colnames(tcga_counts_raw_pF) == c(ERpos_samples, ERneg_samples))

# for tcga -> takes the one filtered
tcga_counts_all <- tcga_counts_raw_pF
stopifnot(tcga_counts_all >= 0)

stopifnot(rownames(tcga_counts_all) == gene_dt$geneID)
stopifnot(nrow(tcga_counts_all) == nrow(gene_dt))
stopifnot(colnames(tcga_counts_all) == c(ERpos_samples, ERneg_samples))
stopifnot(ncol(tcga_counts_all) == nERneg+nERpos)

stopifnot(colnames(tcga_counts_all) == rownames(tcga_sampleAnnot))
stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
colnames(tcga_counts_all) <- tcga_sampleAnnot$cgc_sample_id


outFile <- file.path(outFolder, "tcga_sampleAnnot.Rdata")
save(tcga_sampleAnnot, file=outFile)
cat(paste0("... written ", outFile, "\n"))


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
