
# Rscript phenopath_recount2_OV_pcaScaleFalse.R

require(recount)
require(TCGAbiolinks)
require(biomaRt)
require(phenopath)
require(ggplot2)

pcaplot <- function(pca_dt, pctoplot, summ_dt,...) {
  stopifnot(length(pctoplot) == 2)
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  plot(x = pca_dt[, pctoplot[1]],
       y = pca_dt[, pctoplot[2]],
       xlab=paste0("PC", pctoplot[1], " (", var1, " % variance explained)"),
       ylab=paste0("PC", pctoplot[2], " (", var2, " % variance explained)"),
       pch = 16,
       cex=0.7,
       ...
  )
}

runPheno <- TRUE
runNorm <- TRUE
#recount
# https://f1000research.com/articles/6-1558/v1
# why not using gene symbols: https://www.biostars.org/p/352492/#352535

# see workflow from
# https://github.com/ELELAB/LUAD_LUSC_TCGA_comparison/blob/master/6-recount/unifiedLUAD_Rail_18062018.R
# then
# https://github.com/ELELAB/LUAD_LUSC_TCGA_comparison/blob/master/6-recount/LUAD/DE_unified_LUAD.R
# from this article https://bmccancer.biomedcentral.com/track/pdf/10.1186/s12885-019-5965-x.pdf


# I load the data retrieved from recount2 scripts
# so far:
# purity filter
# scale() from 

# TODO: 
# - gc_content normalization (TCGAbiolinks)
# - filter lowly expressed genes

plotType <- "png"
myHeight <- myWidth <- 400

# filter to keep only the highly variable genes
# Campbell 2018: 
# for BRAC, whose variance in log(TPM+1) expression was greater than 1 and whose median absolute deviation was greater than 0
# for COAD: s whose median absolute deviation in log(TPM+1) expression was greater than sqrt(0.5)
mad_thresh <- sqrt(0.5)


purityFilter <- 0.6

inFolder <- file.path("..","tcga_data","DOWNLOAD_TCGA_GTEX_RECOUNT2")

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

outFolder <- "PHENOPATH_RECOUNT2_TCGA_OV_PCASCALEFALSE"
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
if(runNorm) {
  ov_dat_gcNorm <- TCGAanalyze_Normalization(tabDF = ov_data_forNorm,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")
  outFile <- file.path(outFolder, "ov_data_gcNorm.Rdata")
  save(ov_dat_gcNorm, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "ov_data_gcNorm.Rdata")
  ov_dat_gcNorm <- get(load(outFile))
}

ov_data_gcNorm_log <- log2(ov_dat_gcNorm + 1)

nGTEX <- sum(grepl("^GTEX", colnames(ov_data_gcNorm_log)))
nTCGA <- sum(grepl("^TCGA", colnames(ov_data_gcNorm_log)))

################################### 
################################### PCA ON not norm DATA -> for comparison
################################### 
ov_data_notNorm_log <- log2(ov_data_forNorm + 1)

# ov_data_notNorm_log_no0 <- ov_data_notNorm_log[rowSums(ov_data_notNorm_log) > 0,]
ov_data_notNorm_log_no0 <- ov_data_notNorm_log

pca_ov <- prcomp(t(ov_data_notNorm_log_no0), scale=FALSE)
pca_ov_lowrepr <- pca_ov$x
stopifnot(nrow(pca_ov_lowrepr) == nGTEX+nTCGA)

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_tcga_gtex.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()

pcaplot(pca_dt=pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_ov),
        main="TCGA+GTEX OV notNorm (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_ov_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_tcga_gtex.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_ov),
        main="TCGA+GTEX OV notNorm (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_ov_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))


#### check if I have some outliers within TCGA or GTEX cohorts
gtex_raw_data_log <- ov_data_notNorm_log[,grep("^GTEX", colnames(ov_data_notNorm_log))]
tcga_raw_data_log <- ov_data_notNorm_log[,grep("^TCGA", colnames(ov_data_notNorm_log))]
stopifnot(ncol(gtex_raw_data_log) + ncol(tcga_raw_data_log) == ncol(ov_data_notNorm_log))

### FOR PCA -> I need to remove genes that sum to 0 otherwise cannot scale data
# gtex_raw_data_log <- gtex_raw_data_log[rowSums(gtex_raw_data_log) > 0,]
# tcga_raw_data_log <- tcga_raw_data_log[rowSums(tcga_raw_data_log) > 0,]


# do I have batch effect within TCGA (plate)
# The default is FALSE for consistency with S, but in general scaling is advisable
# also scaled in https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
tcga_pca_ov <- prcomp(t(tcga_raw_data_log), scale=FALSE)
tcga_pca_ov_lowrepr <- tcga_pca_ov$x
stopifnot(nrow(tcga_pca_ov_lowrepr) == nTCGA)

stopifnot(tcga_annot_dt$cgc_sample_id == colnames(tcga_raw_data_log))
tcga_annot_dt$cgc_sample_tissue_source_site_fact <- as.factor(tcga_annot_dt$cgc_sample_tissue_source_site)


#dev.off()
outFile <- file.path(outFolder, paste0("in_raw_data_tcga_pca_12.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=tcga_pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(tcga_pca_ov),
        col = tcga_annot_dt$cgc_sample_tissue_source_site_fact,
        main="TCGA OV notNorm (log2(.+1))")
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))

#dev.off()
outFile <- file.path(outFolder, paste0("in_raw_data_tcga_pca_23.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=tcga_pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(tcga_pca_ov),
        col = tcga_annot_dt$cgc_sample_tissue_source_site_fact,
        main="TCGA OV notNorm (log2(.+1))")
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))


# do I have batch effect within gtex
gtex_pca_ov <- prcomp(t(gtex_raw_data_log), scale=FALSE)
gtex_pca_ov_lowrepr <- gtex_pca_ov$x
stopifnot(nrow(gtex_pca_ov_lowrepr) == nGTEX)


#dev.off()
outFile <- file.path(outFolder, paste0("in_raw_data_gtex_pca_12.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(gtex_pca_ov),
        main="GTEX OV notNorm (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# here there is an issue with one guy !!!
dev.off()
cat(paste0("... written: ", outFile, "\n"))

#dev.off()
outFile <- file.path(outFolder, paste0("in_raw_data_gtex_pca_23.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(gtex_pca_ov),
        main="GTEX OV notNorm (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# here there is an issue with one guy !!!
dev.off()
cat(paste0("... written: ", outFile, "\n"))

rm(gtex_raw_data_log)
rm(tcga_raw_data_log)

################################### 
################################### PCA ON NORM DATA 
################################### 

#### check if I have some outliers within TCGA or GTEX cohorts
gtex_data_log <- ov_data_gcNorm_log[,grep("^GTEX", colnames(ov_data_gcNorm_log))]
tcga_data_log <- ov_data_gcNorm_log[,grep("^TCGA", colnames(ov_data_gcNorm_log))]

stopifnot(ncol(gtex_data_log) + ncol(tcga_data_log) == ncol(ov_data_gcNorm_log))

### FOR PCA -> I need to remove genes that sum to 0 otherwise cannot scale data
# tcga_data_log <- tcga_data_log[rowSums(tcga_data_log) > 0,]
# gtex_data_log <- gtex_data_log[rowSums(gtex_data_log) > 0,]

# do I have batch effect within TCGA (plate)
# The default is FALSE for consistency with S, but in general scaling is advisable
# also scaled in https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
tcga_pca_ov <- prcomp(t(tcga_data_log), scale=FALSE)
tcga_pca_ov_lowrepr <- tcga_pca_ov$x
stopifnot(nrow(tcga_pca_ov_lowrepr) == nTCGA)

stopifnot(tcga_annot_dt$cgc_sample_id == colnames(tcga_data_log))
tcga_annot_dt$cgc_sample_tissue_source_site_fact <- as.factor(tcga_annot_dt$cgc_sample_tissue_source_site)

#dev.off()
outFile <- file.path(outFolder, paste0("in_data_tcga_pca_12.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=tcga_pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(tcga_pca_ov),
        col = tcga_annot_dt$cgc_sample_tissue_source_site_fact,
        main="TCGA OV (log2(.+1))")
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))

#dev.off()
outFile <- file.path(outFolder, paste0("in_data_tcga_pca_23.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=tcga_pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(tcga_pca_ov),
        col = tcga_annot_dt$cgc_sample_tissue_source_site_fact,
        main="TCGA OV (log2(.+1))")
mtext(text=paste0("nTCGA=",nTCGA, " (color-coded by tissue source site)"), side=3)
legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))


# do I have batch effect within gtex
gtex_pca_ov <- prcomp(t(gtex_data_log), scale=FALSE)
gtex_pca_ov_lowrepr <- gtex_pca_ov$x
stopifnot(nrow(gtex_pca_ov_lowrepr) == nGTEX)


#dev.off()
outFile <- file.path(outFolder, paste0("in_data_gtex_pca_12.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(gtex_pca_ov),
     main="GTEX OV (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# here there is an issue with one guy !!!
dev.off()
cat(paste0("... written: ", outFile, "\n"))

# outlier on PC2 !!!
################################## remove outlier in GTEX

#dev.off()
outFile <- file.path(outFolder, paste0("in_data_gtex_pca_23.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(gtex_pca_ov),
        main="GTEX OV (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))


pc2_outlier <- names(which.max(abs(gtex_pca_ov_lowrepr[,2])))
stopifnot(length(pc2_outlier) == 1)
gtex_annot_dt[gtex_annot_dt$sampid == pc2_outlier,]

# remove it from everywhere
ov_data_gcNorm_log <- ov_data_gcNorm_log[,colnames(ov_data_gcNorm_log) != pc2_outlier]
stopifnot(ncol(ov_data_gcNorm_log) == nGTEX + nTCGA -1)
nGTEX <- nGTEX -1 
gtex_data_log <- ov_data_gcNorm_log[,grep("^GTEX", colnames(ov_data_gcNorm_log))]
stopifnot(ncol(gtex_data_log) == nGTEX)
gtex_annot_dt <- gtex_annot_dt[gtex_annot_dt$sampid != pc2_outlier,]
stopifnot(gtex_annot_dt$sampid == colnames(gtex_data_log))

########################### REDO THE GTEX PCA WITHOUT THE OUTLIER
# do I have batch effect within gtex
gtex_pca_ov <- prcomp(t(gtex_data_log), scale=FALSE)
gtex_pca_ov_lowrepr <- gtex_pca_ov$x
stopifnot(nrow(gtex_pca_ov_lowrepr) == nGTEX)

#dev.off()
outFile <- file.path(outFolder, paste0("in_data_gtex_pca_12_no_outlier.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(gtex_pca_ov),
        main="GTEX OV (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# here there is an issue with one guy !!!
dev.off()
cat(paste0("... written: ", outFile, "\n"))

#dev.off()
outFile <- file.path(outFolder, paste0("in_data_gtex_pca_23_no_outlier.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
pcaplot(pca_dt=gtex_pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(gtex_pca_ov),
        main="GTEX OV (log2(.+1))")
mtext(text=paste0("nGTEX=",nGTEX, ""), side=3)
#legend("topleft", legend=c("color-code source site"), bty="n")
# apparently there is no batch effect
dev.off()
cat(paste0("... written: ", outFile, "\n"))

rm(gtex_data_log)
rm(tcga_data_log)
########################### KEEP ONLY VARIABLE GENES
# median absolute deviation in log(TPM+1)
# filter to keep only the highly variable genes
### this should be done after having checked for outliers...
nGenes <- nrow(ov_data_gcNorm_log)
all_mad <- apply(ov_data_gcNorm_log, 1, mad)
to_keep <- all_mad > mad_thresh
stopifnot(length(to_keep) == nGenes)
cat(paste0(sum(to_keep), "/", length(to_keep), "\n"))

plot(density(all_mad), main="all gene MADs (GTEX+TCGA)")
abline(v = mad_thresh, col="red")
legend("topright",legend=c(
  paste0("threshold = ", round(mad_thresh, 2)),
  paste0("to keep: ", sum(to_keep), "/", length(to_keep), "\n")), bty="n")

stopifnot(length(to_keep) == nrow(ov_data_gcNorm_log))
stopifnot(nGenes == nrow(ov_data_gcNorm_log))

ov_data_filtered <- ov_data_gcNorm_log[to_keep,]
stopifnot(nrow(ov_data_filtered) == sum(to_keep))
nGenes <- sum(to_keep)

outFile <- file.path(outFolder, paste0("ov_data_filtered.Rdata"))
save(ov_data_filtered, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


### look how different are the TCGA and the GTEX data
##### distribution of the counts - before MAD filter
gtex_density_all <- density(ov_data_gcNorm_log[,grep("^GTEX", colnames(ov_data_gcNorm_log))])
tcga_density_all <- density(ov_data_gcNorm_log[,grep("^TCGA", colnames(ov_data_gcNorm_log))])

# dev.off()
outFile <- file.path(outFolder, paste0("all_counts_log_tcga_gtex_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(gtex_density_all, main="data density (log2(.+1))", 
     xlim=range(c(gtex_density_all$x, tcga_density_all$x)),
     ylim=range(c(gtex_density_all$y, tcga_density_all$y)))
mtext(side=3, text="(all)")
lines(tcga_density_all, col="red")
legend("topright", legend=c("GTEX", "TCGA"), 
       bty="n",
       col=c("black", "red"), lty=c(1,1))
dev.off()
cat(paste0("... written: ", outFile, "\n"))

##### distribution of the counts - after MAD filter
gtex_density_madF <- density(ov_data_filtered[,grep("^GTEX", colnames(ov_data_filtered))])
tcga_density_madF <- density(ov_data_filtered[,grep("^TCGA", colnames(ov_data_filtered))])

# dev.off()
outFile <- file.path(outFolder, paste0("MADfiltered_counts_log_tcga_gtex_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(gtex_density_madF, main="data density (log2(.+1))", 
     xlim=range(c(gtex_density_madF$x, tcga_density_madF$x)),
     ylim=range(c(gtex_density_madF$y, tcga_density_madF$y)))
mtext(side=3, text="(MAD filtered)")
lines(tcga_density_madF, col="red")
legend("topright", legend=c("GTEX", "TCGA"), 
       bty="n",
       col=c("black", "red"), lty=c(1,1))
dev.off()
cat(paste0("... written: ", outFile, "\n"))


##### distribution of the mads - before filtering
gtex_density_all <- density(apply(ov_data_gcNorm_log[,grep("^GTEX", colnames(ov_data_gcNorm_log))], 1, mad))
tcga_density_all <- density(apply(ov_data_gcNorm_log[,grep("^TCGA", colnames(ov_data_gcNorm_log))], 1, mad))

# dev.off()
outFile <- file.path(outFolder, paste0("all_mads_tcga_gtex_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(gtex_density_all, main="data density (log2(.+1))", 
     xlim=range(c(gtex_density_all$x, tcga_density_all$x)),
     ylim=range(c(gtex_density_all$y, tcga_density_all$y)))
lines(tcga_density_all, col="red")
mtext(side=3, text="(MAD filtered)")
legend("topright", legend=c("GTEX", "TCGA"), 
       bty="n",
       col=c("black", "red"), lty=c(1,1))
dev.off()
cat(paste0("... written: ", outFile, "\n"))


##### distribution of the mads - after filtering
gtex_density_madF <- density(apply(ov_data_filtered[,grep("^GTEX", colnames(ov_data_filtered))], 1, mad))
tcga_density_madF <- density(apply(ov_data_filtered[,grep("^TCGA", colnames(ov_data_filtered))], 1, mad))

# dev.off()
outFile <- file.path(outFolder, paste0("MADfiltered_mads_tcga_gtex_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(gtex_density_madF, main="data density (log2(.+1))", 
     xlim=range(c(gtex_density_madF$x, tcga_density_madF$x)),
     ylim=range(c(gtex_density_madF$y, tcga_density_madF$y)))
lines(tcga_density_madF, col="red")
mtext(side=3, text="(MAD filtered)")
legend("topright", legend=c("GTEX", "TCGA"), 
       bty="n",
       col=c("black", "red"), lty=c(1,1))
dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################################################
############################################################## PCA on merged data
##############################################################


# phenopath tuto: sim$y is the N×G matrix of gene expression (N=100 cells and G=40 genes)
# so I have to take the transpose of my ov_data_log_matrix to have samples in line and genes in columns
ov_data_filteredT <- t(ov_data_filtered)
stopifnot(nGenes == ncol(ov_data_filteredT))
stopifnot(nGTEX+nTCGA == nrow(ov_data_filteredT))

pca_ov <- prcomp(ov_data_filteredT, scale=FALSE)
pca_ov_lowrepr <- pca_ov$x
stopifnot(nrow(pca_ov_lowrepr) == nGTEX+nTCGA)

pc2 <- pca_ov_lowrepr[,2]
pc1 <- pca_ov_lowrepr[,1]

outFile <- file.path(outFolder, paste0("in_data_pca_12_tcga_gtex.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()

pcaplot(pca_dt=pca_ov_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_ov),
        main="TCGA+GTEX OV (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_ov_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_data_pca_23_tcga_gtex.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_ov_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_ov),
        main="TCGA+GTEX OV (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_ov_lowrepr))))
mtext(text=paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA), side=3)
legend("topleft", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################################################
############################################################## LET'S DO PHENO
##############################################################

# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(ov_data_filteredT) == nGTEX+nTCGA)
stopifnot(ncol(ov_data_filteredT) == nGenes)
mycovar <- as.numeric(grepl("^GTEX", rownames(ov_data_filteredT)))
stopifnot(sum(mycovar) == nGTEX)
mycovar[mycovar == 1] <- "GTEX"
mycovar[mycovar == 0] <- "TCGA"
mycovar <- factor(mycovar, levels=c("GTEX", "TCGA"))
stopifnot(!is.na(mycovar))

if(runPheno) {
  ov_phenopath_fit <- phenopath(ov_data_filteredT, mycovar, elbo_tol = 1e-6, thin = 40)
  
  outFile <- file.path(outFolder, paste0("ov_phenopath_fit.Rdata"))
  save(ov_phenopath_fit, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("ov_phenopath_fit.Rdata"))
  ov_phenopath_fit <- get(load(outFile))
}

system.exit(0)


# it is important to check convergence with a call to plot_elbo(ov_phenopath_fit) to ensure the ELBO is approximately flat:
  
  plot_elbo(ov_phenopath_fit)
  
  dev.off()
  qplot(pc1, trajectory(ov_phenopath_fit)) +
    xlab("pc1") + ylab("Phenopath z")
  
  dev.off()
  qplot(pc2, trajectory(ov_phenopath_fit)) +
    xlab("pc2") + ylab("Phenopath z")
  

  gene_names <- colnames(ov_data_filteredT)
  
  df_beta <- data.frame(beta = interaction_effects(ov_phenopath_fit),
                        beta_sd = interaction_sds(ov_phenopath_fit),
                        is_sig = significant_interactions(ov_phenopath_fit),
                        gene = gene_names)
  
  # from 
  df_beta$gene <- fct_relevel(df_beta$gene, gene_names)
  
  ggplot(df_beta, aes(x = gene, y = beta, color = is_sig)) + 
    geom_point() +
    geom_errorbar(aes(ymin = beta - 2 * beta_sd, ymax = beta + 2 * beta_sd)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylab(expression(beta)) +
    scale_color_brewer(palette = "Set2", name = "Significant")
  
  
  which_largest <- which.max(df_beta$beta)
  
  df_large <- data.frame(
    y = ov_data_filteredT[, which_largest],
    x = mycovar,
    z = trajectory(ov_phenopath_fit)
  )
  dev.off()
  
  ggplot(df_large, aes(x = z, y = y, color = x)) +
    geom_point() +
    scale_color_brewer(palette = "Set1") +
    stat_smooth()
  
  
gene_lab_dt <- get(load(file.path(inFolder, "out_intersect_dt.RData")))
save(out_intersect_dt, file=outFile)
cat(paste0("... written ", outFile, "\n"))


gene_lab_dt$geneID_short <-gsub("\\..*", "",gene_lab_dt$geneID)

  
stopifnot(colnames(ov_data_filteredT))
stopifnot(colnames(ov_data_filteredT) %in% gene_lab_dt$geneID_short)  

unique(sapply(colnames(ov_data_filteredT), function(x) sum(x %in% gene_lab_dt$geneID_short)))
# [1] 1   # there is no ambiguity in gene name
sub_gene_lab_dt <- gene_lab_dt[gene_lab_dt$geneID_short %in% colnames(ov_data_filteredT), ]
sub_gene_lab_dt$mapping_source <- NULL
sub_gene_lab_dt <- unique(sub_gene_lab_dt)
any(duplicated(sub_gene_lab_dt$geneID_short))
any(duplicated(sub_gene_lab_dt$geneSymb))

sub_gene_lab_dt[sub_gene_lab_dt$geneID_short == df_beta$gene[which_largest],]
# PLXNC1 -> tumor suppressor gene !


stopifnot(rownames(ov_data_filteredT) == c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id))


all_traj <- setNames(trajectory(ov_phenopath_fit), rownames(ov_data_filteredT))
stopifnot(names(all_traj) == c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id))

stg_ords <-c(  "Stage IC", "Stage IIA" ,"Stage IIB" , "Stage IIC","Stage IIIA", "Stage IIIB","Stage IIIC",   "Stage IV")

tcga_traj_dt <- data.frame(
  tcga_samp = tcga_annot_dt$cgc_sample_id,
  tcga_stage = tcga_annot_dt$cgc_case_clinical_stage,
  pseudotime = all_traj[tcga_annot_dt$cgc_sample_id],
  stringsAsFactors = FALSE
)
sum(is.na(tcga_traj_dt$tcga_stage))
# 3
tcga_traj_dt <- tcga_traj_dt[!is.na(tcga_traj_dt$tcga_stage),]
stopifnot(!is.na(tcga_traj_dt))
tcga_traj_dt$tcga_stage <- factor(tcga_traj_dt$tcga_stage, levels=stg_ords)
stopifnot(!is.na(tcga_traj_dt$tcga_stage))

ggplot(tcga_traj_dt, aes(x= tcga_stage, y= pseudotime) )+
  geom_boxplot()



age_ords <- sort(unique(gtex_annot_dt$AGE))
  
gtex_traj_dt <- data.frame(
    gtex_samp = gtex_annot_dt$sampid,
  gtex_age = gtex_annot_dt$AGE,
  pseudotime = all_traj[gtex_annot_dt$sampid],
  stringsAsFactors = FALSE
)
sum(is.na(gtex_traj_dt$gtex_age))
# 3
gtex_traj_dt <- gtex_traj_dt[!is.na(gtex_traj_dt$gtex_age),]
stopifnot(!is.na(gtex_traj_dt))
gtex_traj_dt$gtex_age <- factor(gtex_traj_dt$gtex_age, levels=age_ords)
stopifnot(!is.na(gtex_traj_dt$gtex_age))

ggplot(gtex_traj_dt, aes(x= gtex_age, y= pseudotime) )+
  geom_boxplot()




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
pca_ov_norm <- prcomp(t(ov_data_log_norm), scale=FALSE)
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


