
# Rscript phenopath_TCGA_OV.R

# Collado-Torres6 et al., developed recount2, where >4.4 trillion RNA-seq raw reads from >70,000 human 
# samples deposited in TCGA, GTEx and Short Read Archive (SRA) were processed uniformly using the Rail-RNA 
# pipeline7. Raw sequencing reads from all studies were processed identically to maximize cross-sample compa-
#   rability.
# Vivian9 et al., developed Toil, a robust workflow which processed >20,000 RNA-seq samples from four 
# large scale studies (including TCGA and GTEx) in under four days using 32,000 commercial cloud-computing 
# preemptable cores. 
# We show that although the 
# majority of protein coding genes show arguably acceptable concordance between pipelines, for >12% of protein 
# coding genes in the genome, different pipelines produce gene expression estimates that are different by more than 
# four-fold.

# The differences among these pipelines arise from diverse imple-
#   mentation choices including statistical and algorithmic methods, software versions, and run-time parameters. 
# The discordant abundance and fold-change estimates revealed here do not imply any technical errors. Rather 
# they highlight inherent uncertainty in processing noisy and complex data. Nonetheless, the end result is that for 
# the DQ genes reported here (Supplementary Table S2b), we can have little confidence in the abundance estimates 
# produced by any RNA-seq processing pipeline. For critical applications such as biomedical research and clinical 
# practice, a concerted, community-wide effort will be needed to develop gold-standards for estimating the mRNA 
# abundance of these genes

# Unifying cancer and normal RNA sequencing data from different sources
# -> provide GTEX + TCGA harmonize datasets 
# but not available for ovarian cancer (no normal samples from TCGA)
# https://github.com/mskcc/RNAseqDB
# https://www.nature.com/articles/sdata201861

# might still need to correct batch effect:
# From https://www.nature.com/articles/nbt.3838
# Although all recount2 samples have been processed and summarized with a single pipeline, 
# so-called 'batch' effects could occur and should be considered in downstream analyses, particularly when comparing among studies. 

# ## scale_counts: Scale counts by taking into account the total coverage per sample
# [I think implicitly correct for gene length]
# cf.  supp note of the article 
# https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3838/MediaObjects/41587_2017_BFnbt3838_MOESM35_ESM.pdf
# they do limma voom directly after scale_count without further correction

# from phenopath tuto
# The expression data should represent something comparable to logged counts, such as log2(TPM+1)

# from Campbell and Yau 2018:
# for COAD dataset: 4,801 genes whose MAD in log(TPM+1) expression > sqrt(0.5)
# for BRCA dataset: 4,579 genes whose variance in log(TPM+1) expression was greater than 1 and whose MAD > 0.
# I think this is log2: cf. phenopath tuto + Campbell thesis (https://kieranrcampbell.github.io/blog/files/krc-thesis.pdf)

# NB - maybe worth to check
# https://github.com/kieranrcampbell/tcga_phenotime/blob/master/analysis/exploratory_OV_tcga.Rmd
# >>>> Exploratory analysis of TCGA OV data
# If we select tau_beta from gamma(10, 0.01) then we get a prior on beta like
# ```{r prior-on-beta}
# qplot(rnorm(10000, 0, 1 / sqrt(rgamma(10000, 1000, 1))), geom = "density")
# ```
# which seems reasonable.
# https://github.com/kieranrcampbell/phenopath_analyses/blob/master/analysis/OV/clvm_analysis.Rmd
# >>>> Immune-mutational burden interactions in ovarian cancer using Variational Bayes Phenotime"
# https://github.com/kieranrcampbell/tcga_phenotime/blob/master/analysis/OV_metabolomics.Rmd
# >>>> Exploratory analysis of TCGA OV data - metabolomics
# Load genes associated with "lipid metabolism":
#   ```{r load-lipid-met}
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# go_lipid <- "GO:0006629" # lipid metabolic process
# xx <- as.list(GOBPCHILDREN)
# all_lipid_children <- xx[[go_lipid]]
# lipid_annotation <- getBM("ensembl_transcript_id", filters = "go_id",
#                           values = all_lipid_children, mart = ensembl)
# Let's remove the batch effect:
# 
# ```{r remove-batch}
# sc$plate <- droplevels(sc$plate)
# batch <- sc$plate
# design <- model.matrix(~ 1, data = pData(sc))
# exprs_combat <- ComBat(exprs(sc), batch, design)
# exprs(sc) <- exprs_combat
# Subset to those with mutations:
#   
#   ```{r sce-lipid2}
# sc_lipid <- sce_lipid[, !is.na(sce_lipid$n_mutations)]
# sc_lipid <- sc_lipid[rowVars(exprs(sc_lipid)) > 0, ]
# ```
# 
# Quick heatmap of lipids:
#   
#   ```{r heatmap}
# prop_cells_exprs <- rowMeans(exprs(sc_lipid) > 0)
# sc_lipid <- sc_lipid[prop_cells_exprs > 0.2,]# & prop_cells_exprs < 1]
# exprs_for_heatmap <- exprs(sc_lipid) 
# # exprs_for_heatmap[exprs_for_heatmap > 5] <- 5
# col_side_col <- plyr::mapvalues(sc_lipid$censored, from = c(TRUE, FALSE),
#                                 to = brewer.pal(3, "Set1")[1:2])
# heatmap.2(exprs_for_heatmap, trace = "none", col = "viridis",
#           ColSideColors = col_side_col)



outFolder <- "PHENOPATH_TCGA_OV"
dir.create(outFolder, recursive=TRUE)

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

# load input data
ov_data <- get(load("../tcga_data/DOWNLOAD_TCGA_GTEX_RECOUNT2/all_counts_dsCountsGe_0.Rdata"))
dim(ov_data)
# 25224x538
stopifnot(ov_data >= 0)
ov_data_log <- log2(ov_data + 1)

gtex_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_GTEX_RECOUNT2/gtex_sampleAnnot.RData"))
tcga_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_GTEX_RECOUNT2/tcga_sampleAnnot.Rdata"))

### look how different are the TCGA and the GTEX data

gtex_data_log <- ov_data_log[,grep("^GTEX", colnames(ov_data_log))]
tcga_data_log <- ov_data_log[,grep("^TCGA", colnames(ov_data_log))]

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
ov_data_logT <- t(ov_data_log)

nGTEX <- sum(grepl("^GTEX", rownames(ov_data_logT)))
nTCGA <- sum(grepl("^TCGA", rownames(ov_data_logT)))
stopifnot(nGTEX+nTCGA == nrow(ov_data_logT))


pca_ov <- prcomp(ov_data_logT)
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

tcga_annot_dt$cgc_sample_tissue_source_site_fact <- as.factor(tcga_annot_dt$cgc_sample_tissue_source_site)

dev.off()
plot(x = tcga_pca_ov_lowrepr[,1],
     y = tcga_pca_ov_lowrepr[,2],
     xlab="PC1", ylab="PC2",
     main="TCGA OV (log2(.+1))",
     pch=16,
     col = tcga_annot_dt$cgc_sample_tissue_source_site_fact)
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
# apparently there is no batch effect


for(i in colnames(tcga_annot_dt)) {
  x <- table(tcga_annot_dt[,i])
  if(length(x) <= 10) {
    print(i)
    print(x)
  }
}



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
