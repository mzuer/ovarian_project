require(TCGAbiolinks)
require(SummarizedExperiment)
require(recount)
# Rscript download_tcga_gtex_recount2.R

# https://github.com/LieberInstitute/recount-brain/blob/master/SupplementaryTable2.csv


startTime <- Sys.time()

cat(paste0("> START ", startTime, "\n"))

outFolder <- "DOWNLOAD_TCGA_GTEX_RECOUNT2"
dir.create(outFolder)


tcga_purity_thresh <- 0.6 # cf Lucchetta et al. 2019

compute_meds <- TRUE


# update
# no filter based on duplicated ensembl ID but remove the NA gene symbol (I guess might be the non-protein coding)


# are tere all high grades ??????

############################################## 
############################################## prepare TCGA data
############################################## 
ovary_rec2_gtex <- TCGAquery_recount2(project="gtex", tissue = "ovary")
ovary_rec2_gtex_scaled <- scale_counts(ovary_rec2_gtex$gtex_ovary)

ovary_rec2_tcga <- TCGAquery_recount2(project="tcga", tissue = "ovary")
ovary_rec2_tcga_scaled <- scale_counts(ovary_rec2_tcga$tcga_ovary)

stopifnot(ovary_rec2_tcga_scaled@colData[,"cgc_case_primary_site"] == "Ovary")


annot_cols <- c(grep("^cgc_", names(ovary_rec2_tcga_scaled@colData)), grep("^gdc_", names(ovary_rec2_tcga_scaled@colData)))
tcgacols <- names(ovary_rec2_tcga_scaled@colData)[annot_cols]
tcga_sampleAnnot <- data.frame(ovary_rec2_tcga_scaled@colData[, tcgacols])

stopifnot(tcga_sampleAnnot[,"gdc_cases.project.name"] == "Ovarian Serous Cystadenocarcinoma")
stopifnot(tcga_sampleAnnot[,"cgc_file_disease_type"] == "Ovarian Serous Cystadenocarcinoma")
stopifnot(tcga_sampleAnnot[,"cgc_case_histological_diagnosis"] == "Serous Cystadenocarcinoma")
stopifnot(tcga_sampleAnnot[,"cgc_file_experimental_strategy"] == "RNA-Seq")# 


####~~~ 1st FILTER HERE -> TAKE ONLY PRIMARY TUMOR !!!
stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_counts_raw_all))
#*** update here
tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$cgc_sample_sample_type == "Primary Tumor",]
####~~~2nd ADDED FILTER -> discard FFPE samples
# (seems to be a good practice, e.g.: https://www.biostars.org/p/308192/)
stopifnot(tcga_sampleAnnot$cgc_sample_is_ffpe == ifelse(tcga_sampleAnnot$gdc_cases.samples.is_ffpe, "YES", "NO"))
cat(paste0("... remove FFPE: ", sum(tcga_sampleAnnot$gdc_cases.samples.is_ffpe), "\n"))
#*** update here
tcga_sampleAnnot <- tcga_sampleAnnot[! tcga_sampleAnnot$gdc_cases.samples.is_ffpe,]
stopifnot(nrow(tcga_sampleAnnot) > 0)
#*** update here
tcga_counts_raw_all <- tcga_counts_raw_all[,rownames(tcga_sampleAnnot)]
stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_counts_raw_all))


###~~~3d FILTER - protein coding only? DO THIS BEFORE COMPUTING THE MEDIAN !  -> to remove the NA gene Symbols (non coding ???)
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
#*** update here:
tcga_counts_raw_all <- tcga_counts_raw_all[gene_dt$geneID,]
stopifnot(!is.na(tcga_counts_raw_all))
dim(tcga_counts_raw_all)


####~~~ 4th FILTER:  PURITY FILTER FOR THE TUMOR SAMPLES - do this before removing duplicated ones !!!
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
#*** update here
# ~~~~ KEEP ONLY THE ONES THAT PASSED THE PURITY FILTER
tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$labs_for_purity%in% purityinfo_brca$pure_barcodes,]
stopifnot(tcga_sampleAnnot$cgc_sample_id %in% colnames(tcga_counts_raw_all))
tcga_counts_raw_all <- tcga_counts_raw_all[,colnames(tcga_counts_raw_all) %in% tcga_sampleAnnot$cgc_sample_id]  
stopifnot(!is.na(tcga_counts_raw_all))  
# 5th filter -> e.g. ER status
# not needed here


###~~~ 6th FILTER: DUPLICATED SAMPLES
###~~~ FILTER 6: KEEP MAX MEDS FOR DUPLICATED SAMPLES
cat("any duplicated(tcga_sampleAnnot$cgc_sample_id)", "\n")
cat(any(duplicated(tcga_sampleAnnot$cgc_sample_id)), "\n")
# stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
# this is not TRUE ! I have duplicated samples
# see discussion here: https://www.biostars.org/p/311017/
# when duplicated, I take the one with highest median value
nGenesCheck <- nrow(tcga_counts_raw_all)
withdup_count_dt <- as.data.frame(t(tcga_counts_raw_all))
stopifnot(rownames(withdup_count_dt) == rownames(tcga_sampleAnnot))

##### UPDATE -> they might be not only duplicated analytes, but also duplicated vials
# I noticed it when doing survival analysis...
###withdup_count_dt$samp_id <- tcga_sampleAnnot$cgc_sample_id
###sum(duplicated(withdup_count_dt$samp_id))
# [1] 23
# check all should Solid Tumor; see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
stopifnot(substr(start=14,stop=15,x=tcga_sampleAnnot$cgc_sample_id) %in% c("01"))# because filter primary tumor
# update: could differ by vial
tcga_sampleAnnot$id_for_survival <- gsub("(^.+?-.+?-.+?)-.+", "\\1", tcga_sampleAnnot$cgc_sample_id)
cat("any duplicated(tcga_sampleAnnot$id_for_survival)", "\n")
cat(any(duplicated(tcga_sampleAnnot$id_for_survival)), "\n")
head(colnames(tcga_counts_raw_all))
head( tcga_sampleAnnot$id_for_survival)

withdup_count_dt$samp_id <- tcga_sampleAnnot$id_for_survival
sum(duplicated(withdup_count_dt$samp_id))

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


##*** update the count data
tcga_counts_raw_all <- new_tcga_counts_raw_all

stopifnot(nrow(tcga_counts_raw_all) == nGenesCheck) 
stopifnot(ncol(tcga_counts_raw_all) == nSampsCheck) 

##*** update annot data
stopifnot(colnames(tcga_counts_raw_all)%in% rownames(tcga_sampleAnnot))
tcga_sampleAnnot <- tcga_sampleAnnot[colnames(tcga_counts_raw_all),] 
stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
stopifnot(!duplicated(tcga_sampleAnnot$id_for_survival))
stopifnot(colnames(tcga_counts_raw_all) == rownames(tcga_sampleAnnot))

tcga_counts_all <- tcga_counts_raw_all
stopifnot(tcga_counts_all >= 0)

stopifnot(rownames(tcga_counts_all) == gene_dt$geneID)
stopifnot(nrow(tcga_counts_all) == nrow(gene_dt))
stopifnot(colnames(tcga_counts_all) == c(ERpos_samples, ERneg_samples))
stopifnot(ncol(tcga_counts_all) == nERneg+nERpos)

stopifnot(colnames(tcga_counts_all) == rownames(tcga_sampleAnnot))
stopifnot(!duplicated(tcga_sampleAnnot$cgc_sample_id))
stopifnot(!duplicated(tcga_sampleAnnot$id_for_survival))
colnames(tcga_counts_all) <- tcga_sampleAnnot$cgc_sample_id

stopifnot(substr(start=1,stop=12,x=colnames(tcga_counts_all)) == tcga_sampleAnnot$id_for_survival)


############################################## 
############################################## prepare GTEX data
############################################## 


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


# take the first symbol of the list... don't know how to do better...
tcga_symbs <- sapply(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
tcga_ids <- as.character(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$gene_id)
# need add it here because now I do the protein coding filter earlier (filter 3)
stopifnot(gene_dt$geneSymb %in% tcga_symbs)
stopifnot(gene_dt$geneID %in% tcga_ids)
tcga_symbs <- tcga_symbs[tcga_symbs %in% gene_dt$geneSymb]
tcga_ids <- tcga_symbs[tcga_symbs %in% gene_dt$geneID]

stopifnot(rownames(tcga_counts_all) == tcga_ids)
stopifnot(length(tcga_ids) == length(tcga_symbs))
stopifnot(!duplicated(tcga_ids))
## here filter -> want the non NA
to_keep <- !is.na(tcga_symbs)
tcga_symbs <- tcga_symbs[to_keep]
tcga_ids <- tcga_ids[to_keep]
tcga_g2s <- setNames(tcga_symbs, tcga_ids)


gtex_symbs <- sapply(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
gtex_ids <- as.character(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(gtex_counts_all) == gtex_ids)
stopifnot(length(gtex_ids) == length(gtex_symbs))
stopifnot(!duplicated(gtex_ids))
## here filter -> want the non NA
to_keep <- !is.na(gtex_symbs)
gtex_symbs <- gtex_symbs[to_keep]
gtex_ids <- gtex_ids[to_keep]
gtex_g2s <- setNames(gtex_symbs, gtex_ids)

# now I do a bit different
# cf post here: https://www.biostars.org/p/352492/#352535
# gene names are more ambiguous, so I will go forward
# to the downstream analyses with the ENSG symbols
# (so I dont discard the gene symbol duplicated genes)
gtex_g2s_filt <- gtex_g2s
tcga_g2s_filt <- tcga_g2s


# I need to do like this because I can have a given gene symbol 
# mapped to  ENSG in tcga and to a different ENSG in gtex
# so with this table I ensure exact gene ID gene Symbol matching
# in this version, I allow an ENSG to map multiple symbols
# but should be the case in both datasets
common_dt = data.frame(
  geneSymb=c(as.character(tcga_g2s_filt), as.character(gtex_g2s_filt)),
  geneID=c(names(tcga_g2s_filt),names(gtex_g2s_filt)),
  mapping_source = c(rep("TCGA", length(gtex_g2s_filt)), rep("GTEX", length(gtex_g2s_filt))),
  stringsAsFactors = FALSE
)
tmp_dt <- common_dt
tmp_dt$mapping_source <- NULL

common_dt$filtid <- paste0(common_dt$geneID, "_", common_dt$geneSymb)
id_to_keep <- names(table(common_dt$filtid))[which(as.numeric(table(common_dt$filtid)) == 2)]
intersect_dt <- common_dt[common_dt$filtid %in% id_to_keep,]
out_intersect_dt <- intersect_dt
intersect_dt$mapping_source  <- NULL
stopifnot(nrow(intersect_dt) == nrow(unique(intersect_dt)) * 2)
intersect_dt <- unique(intersect_dt)
stopifnot(!any(is.na(intersect_dt)))

common_symbs <- as.character(intersect_dt$geneSymb)
common_ids <- as.character(intersect_dt$geneID)
stopifnot(!is.na(common_symbs))
stopifnot(!is.na(common_ids))
stopifnot(length(common_symbs) == length(common_ids))
stopifnot(!duplicated(common_ids))
sum(duplicated(common_symbs)) # 302
# stopifnot(!duplicated(common_symbs)) ### this is not necessary true in this version !
length(common_ids)
# 20698
g2s <- setNames(as.character(intersect_dt$geneSymb), as.character(intersect_dt$geneID))


## => do not do the renaming of row names, keep ensembl ID
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
stopifnot(nrow(all_counts)==length(g2s))
dim(all_counts)
stopifnot(ncol(sub_gtex_counts) + ncol(sub_tcga_counts) == ncol(all_counts))

outFile <- file.path(outFolder, paste0("all_counts_onlyPF_", tcga_purity_thresh, ".Rdata"))
save(all_counts, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "tcga_sampleAnnot.Rdata")
save(tcga_sampleAnnot, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "gtex_sampleAnnot.RData")
save(gtex_sampleAnnot, file=outFile)
cat(paste0("... written ", outFile, "\n"))


outFile <- file.path(outFolder, "out_intersect_dt.RData")
save(out_intersect_dt, file=outFile)
cat(paste0("... written ", outFile, "\n"))


##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))

stop("---ok\n")

####################################################################################################
### THRASH
####################################################################################################
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


rownames(tcga_counts_raw_pF)<-gsub("\\..*", "",rownames(tcga_counts_raw_pF))
library(biomaRt)
# get list of all protein-coding genes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                     filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
head(listOfGenes)
dim(listOfGenes)

tcga_counts_raw_pF <- subset(tcga_counts_raw_pF,rownames(tcga_counts_raw_pF) %in% listOfGenes$ensembl_gene_id)
dim(tcga_counts_raw_pF)

xx=TCGAanalyze_Normalization(tabDF = tcga_counts_raw_pF,
                             geneInfo = geneInfoHT,
                             method = "gcContent")


### was a first idea to do this filtering by dataset
# but as mentioned in limma 1,2,3
# Whilst it is of interest to examine genes that are expressed 
# in one condition but not in another, some genes are unexpressed throughout *all samples*
dataset_rowMeans_filter <- sqrt(0.5) #
tcga_filter <- rowMeans(tcga_counts_all) >= dataset_rowMeans_filter
sum(tcga_filter)
# 38670
gtex_filter <- rowMeans(gtex_counts_all) >= dataset_rowMeans_filter
sum(gtex_filter)
# 37644

stopifnot(length(tcga_ids) == length(tcga_filter))   
stopifnot(length(gtex_ids) == length(gtex_filter))
tcga_g2s <- setNames(tcga_symbs[tcga_filter], tcga_ids[tcga_filter]) ##### keep only the IDs that have passed the count filter
gtex_g2s <- setNames(gtex_symbs[gtex_filter], gtex_ids[gtex_filter]) ##### keep only the IDs that have passed the count filter


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



# take the first symbol of the list... don't know how to do better...
tcga_symbs <- sapply(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
tcga_ids <- as.character(data.frame(ovary_rec2_tcga_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(tcga_counts_all) == tcga_ids)
stopifnot(length(tcga_ids) == length(tcga_symbs))
stopifnot(!duplicated(tcga_ids))
tcga_g2s <- setNames(tcga_symbs, tcga_ids) ##### count filter now done on merged data
any(duplicated(as.character(tcga_g2s)))
tcga_g2s_filt <- tcga_g2s[!duplicated(tcga_g2s)]
tcga_g2s_filt <- tcga_g2s_filt[!is.na(tcga_g2s_filt)]


# the issue is: I cannot make the merging on gene symbols, because duplicated
# I can do the merging on

gtex_symbs <- sapply(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$symbol, function(x)x[[1]])
gtex_ids <- as.character(data.frame(ovary_rec2_gtex_scaled@rowRanges@elementMetadata)$gene_id)
stopifnot(rownames(gtex_counts_all) == gtex_ids)
stopifnot(length(gtex_ids) == length(gtex_symbs))
stopifnot(!duplicated(gtex_ids))
gtex_g2s <- setNames(gtex_symbs, gtex_ids) ##### count filter now done on merged data
any(duplicated(as.character(gtex_g2s)))
gtex_g2s_filt <- gtex_g2s[!duplicated(gtex_g2s)]
gtex_g2s_filt <- gtex_g2s_filt[!is.na(gtex_g2s_filt)]

# do not do the renaming based on ensembl 2 gene symb matching

new_rownames1 <- as.character(gtex_g2s_filt[common_ids])
new_rownames2 <- as.character(tcga_g2s_filt[common_ids])
stopifnot(length(new_rownames2) == length(new_rownames1))
stopifnot(new_rownames2==new_rownames1)
stopifnot(new_rownames2==g2s[common_ids])

rownames(all_counts) <- g2s[common_ids]

stopifnot(length(new_rownames1) == nrow(all_counts))



# update: I think it makes more sense to filter based on the 2 datasets
# (e.g. in Campbell and Yau, they filter the 2 conditions together, see also limma 1,2,3 https://f1000research.com/articles/5-1408)
full_countMeans_filter <- 0
# now done in later script