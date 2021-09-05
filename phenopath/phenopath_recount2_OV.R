
# Rscript phenopath_TCGA_OV.R

require(recount)
require(TCGAbiolinks)
require(biomaRt)
require(phenopath)

# I load the data retrieved from recount2 scripts
# so far:
# purity filter
# scale() from 

# TODO: 
# - gc_content normalization (TCGAbiolinks)
# - filter lowly expressed genes

purityFilter <- 0.6

inFolder <- file.path("..","tcga_data","DOWNLOAD_TCGA_GTEX_RECOUNT2")

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

outFolder <- "PHENOPATH_RECOUNT2_TCGA_OV"
dir.create(outFolder, recursive=TRUE)

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

# load input data
ov_data_raw <- get(load(file.path(inFolder, paste0("all_counts_onlyPF_", purityFilter, ".Rdata"))))
dim(ov_data_raw)
# 25226x533
stopifnot(ov_data_raw >= 0)


gtex_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_GTEX_RECOUNT2/gtex_sampleAnnot.RData"))
tcga_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_GTEX_RECOUNT2/tcga_sampleAnnot.Rdata"))

# 1- gc content normalization
stopifnot(!duplicated(rownames(gsub("\\..*", "",rownames(ov_data_raw)))))
rownames(ov_data_raw)<-gsub("\\..*", "",rownames(ov_data_raw))
# get list of all protein-coding genes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                     filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
head(listOfGenes)
dim(listOfGenes)

ov_data_forNorm <- subset(ov_data_raw,rownames(ov_data_raw) %in% listOfGenes$ensembl_gene_id)
dim(ov_data_forNorm)
# [1] 19177   533


### normalize as in Lucchetta et al. 2019
ov_dat_gcNorm <- TCGAanalyze_Normalization(tabDF = ov_data_forNorm,
                             geneInfo = geneInfoHT,
                             method = "gcContent")


# filter to keep only the highly variable genes
# Campbell 2018: 
# for BRAC, whose variance in log(TPM+1) expression was greater than 1 and whose median absolute deviation was greater than 0
# for COAD: s whose median absolute deviation in log(TPM+1) expression was greater than sqrt(0.5)

mad_thresh <- sqrt(0.5)

ov_data_gcNorm_log <- log2(ov_dat_gcNorm + 1)

all_mad <- apply(ov_data_gcNorm_log, 1, mad)
to_keep <- all_mad > mad_thresh
cat(paste0(sum(to_keep), "/", length(to_keep), "\n"))

stopifnot(length(to_keep) == nrow(ov_data_gcNorm_log))

ov_data_filtered <- ov_data_gcNorm_log[to_keep,]
stopifnot(nrow(ov_data_filtered) == length(to_keep))


### look how different are the TCGA and the GTEX data

gtex_data_log <- ov_data_filtered[,grep("^GTEX", colnames(ov_data_filtered))]
tcga_data_log <- ov_data_filtered[,grep("^TCGA", colnames(ov_data_filtered))]

gtex_density <- density(gtex_data_log)
tcga_density <- density(tcga_data_log)

dev.off()
plot(gtex_density, main="data density (log2(.+1))", 
     xlim=range(c(gtex_density$x, tcga_density$x)),
     ylim=range(c(gtex_density$y, tcga_density$y)))
lines(tcga_density, col="red")
legend("topright", legend=c("GTEX", "TCGA"), 
       bty="n",
       col=c("black", "red"), lty=c(1,1))

# phenopath tuto: sim$y is the N×G matrix of gene expression (N=100 cells and G=40 genes)
# so I have to take the transpose of my ov_data_log_matrix to have samples in line and genes in columns
ov_data_filteredT <- t(ov_data_filtered)

nGTEX <- sum(grepl("^GTEX", rownames(ov_data_filteredT)))
nTCGA <- sum(grepl("^TCGA", rownames(ov_data_filteredT)))
stopifnot(nGTEX+nTCGA == nrow(ov_data_filteredT))

pca_ov <- prcomp(ov_data_filteredT)
pca_ov_lowrepr <- pca_ov$x
stopifnot(nrow(pca_ov_lowrepr) == nGTEX+nTCGA)

dev.off()
plot(x = pca_ov_lowrepr[,1],
     y = pca_ov_lowrepr[,2],
     xlab="PC1", ylab="PC2",
     main="TCGA+GTEX OV (log2(.+1))",
     pch=16,
     col=1+as.numeric(grepl("TCGA", rownames(pca_ov_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")

# do I have batch effect within TCGA (plate)
tcga_pca_ov <- prcomp(t(tcga_data_log))
tcga_pca_ov_lowrepr <- tcga_pca_ov$x
stopifnot(nrow(tcga_pca_ov_lowrepr) == nTCGA)

stopifnot(tcga_annot_dt$cgc_sample_id == colnames(tcga_data_log))
tcga_annot_dt$cgc_sample_tissue_source_site_fact <- as.factor(tcga_annot_dt$cgc_sample_tissue_source_site)

dev.off()
plot(x = tcga_pca_ov_lowrepr[,1],
     y = tcga_pca_ov_lowrepr[,2],
     xlab="PC1", ylab="PC2",
     main="TCGA OV (log2(.+1))",
     pch=16,
     col = tcga_annot_dt$cgc_sample_tissue_source_site_fact)
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect


# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(ov_data_filteredT) == nGTEX+nTCGA)
mycovar <- as.numeric(grepl("^GTEX", rownames(ov_data_filteredT)))
stopifnot(sum(mycovar) == nGTEX)
mycovar[mycovar == 1] <- "GTEX"
mycovar[mycovar == 0] <- "TCGA"
mycovar <- factor(mycovar, levels=c("GTEX", "TCGA"))
stopifnot(!is.na(mycovar))
fit <- phenopath(ov_data_filteredT, mycovar, elbo_tol = 1e-6, thin = 40)


##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))

# TCGA annot of interest

# [1] "cgc_file_platform"
# Illumina HiSeq 
# 430 
# [1] "cgc_file_investigation"
# TCGA-OV 
# 430 
# [1] "cgc_file_disease_type"
# Ovarian Serous Cystadenocarcinoma 
#                               430 
# [1] "cgc_sample_sample_type"
#   Primary Tumor Recurrent Tumor 
#             422               8 
# [1] "cgc_case_vital_status"
# Alive  Dead 
#   192   237 
# [1] "cgc_case_primary_site"
# Ovary 
#   430 
# [1] "cgc_case_clinical_stage"
# 
#   Stage IC  Stage IIA  Stage IIB  Stage IIC Stage IIIA Stage IIIB Stage IIIC   Stage IV 
#          1          3          5         18          7         18        311         64 
# [1] "cgc_case_gender"
# FEMALE 
#    430 
# [1] "cgc_drug_therapy_pharmaceutical_therapy_type"
# 
#               Chemotherapy            Hormone Therapy              Immunotherapy 
#                        376                          6                          2 
# Targeted Molecular therapy 
# 7 





# 
# for(i in colnames(tcga_annot_dt)) {
#   x <- table(tcga_annot_dt[,i])
#   if(length(x) <= 10) {
#     print(i)
#     print(x)
#   }
# }



# sc$plate <- droplevels(sc$plate)
# batch <- sc$plate
# design <- model.matrix(~ 1, data = pData(sc))
# exprs_combat <- ComBat(exprs(sc), batch, design)
# exprs(sc) <- exprs_combat
# Subset to those with mutations:

# make each sample to have the same range of expression
dim(tcga_data_log)
# [1] 25224   430
tcga_data_log_norm <- apply(tcga_data_log, 2, function(x) (x- min(x)) / (max(x)- min(x)))
gtex_data_log_norm <- apply(gtex_data_log, 2, function(x) (x- min(x)) / (max(x)-min(x)))

ov_data_log_norm <- cbind(tcga_data_log_norm, gtex_data_log_norm)
pca_ov_norm <- prcomp(t(ov_data_log_norm))
pca_ov_norm_lowrepr <- pca_ov_norm$x
stopifnot(nrow(pca_ov_norm_lowrepr) == nGTEX+nTCGA)

dev.off()
plot(x = pca_ov_norm_lowrepr[,1],
     y = pca_ov_norm_lowrepr[,2],
     xlab="PC1", ylab="PC2",
     main="TCGA+GTEX OV (log2(.+1))",
     pch=16,
     col=1+as.numeric(grepl("TCGA", rownames(pca_ov_norm_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")


