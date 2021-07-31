require(TCGAbiolinks)
require(SummarizedExperiment)
# Rscript download_tcgadata_TCGABiolinks.R

startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGADATA_TCGABIOLINKS"
dir.create(outFolder)

######## 27k
query.met2 <- GDCquery(project = "TCGA-OV", 
                      legacy = FALSE,
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 27")

ov_res_27 <- getResults(query.met2)
colnames(ov_res_27)
head(ov_res_27$sample_type)
table(ov_res_27$sample_type)
# Primary Tumor     Recurrent Tumor Solid Tissue Normal 
# 582                  19                  12 
GDCdownload(query.met2)
ov.exp_27 <- GDCprepare(query = query.met2, save = F)
dim(ov.exp_27)
# 27578   613

ov_dnaMet_27 <- ov.exp_27
outFile <- file.path(outFolder, "TCGAbiolinks_OV_DNAmet27_hg38.Rdata")
save(ov_dnaMet_27, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

ov.exp_27_clinic <- ov.exp_27@colData

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
write.table(ov.exp_27_probes_dt, file =outFile, col.names = TRUE, row.names=FALSE, sep="\t", quote=F)
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
write.table(ov_met27_data_dt, file =outFile, col.names = TRUE, row.names=FALSE, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))


############################ 450 k => not enough data
# legacy = aligned to hg19

# TCGAbiolinks:::getProjectSummary("TCGA-OV")
# 
# query.met <- GDCquery(project = "TCGA-OV", 
#                       legacy = FALSE,
#                       data.category = "DNA Methylation",
#                       platform = "Illumina Human Methylation 450")
# 
# ov_res <- getResults(query.met)
# 
# 
# colnames(ov_res)
# head(ov_res$sample_type)
# table(ov_res$sample_type)
# # Primary Tumor 
# # 10 
# 
# GDCdownload(query.met)
# ov.exp <- GDCprepare(query = query.met, save = F)
# dim(ov.exp)
# # 485577x10
# ov.exp_clinic <- ov.exp@colData
# table(ov.exp_clinic$sample_type)
# # Primary Tumor 
# # 10 
# sample_type <- sapply(ov.exp_clinic$barcode, substr, 14,15)
# stopifnot(sample_type == "01")
# 
# 
# 
# table(ov.exp_clinic$tumor_stage)
# # not reported
# table(ov.exp_clinic$tumor_grade)
# table(ov.exp_clinic$classification_of_tumor)
# table(ov.exp_clinic$tissue_type)
# # all not reported 
# 
# ov_dnaMet_450 <- ov.exp
# save(ov_dnaMet_450, file=file.path(outFolder, "TCGAbiolinks_OV_DNAmet450_hg38.Rdata"))
