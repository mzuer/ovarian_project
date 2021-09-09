
# Rscript phenopath_recount2_TCGA_BRCA.R

require(recount)
require(TCGAbiolinks)
require(biomaRt)
require(phenopath)
require(ggplot2)
require(ggrepel)
require(ggsci)
require(viridis)
require(dplyr)
require(limma)
require(edgeR)
require(matrixStats)
require(goseq)
genome <- "hg19"
id <- "ensGene"
go_cat <- "GO:BP"
plotNgos <- 10
topCorrThresh <- 0.5  # 0.5 in Campbell and Yau 2018; to select most corr. for GO enrichment
topBetaThresh <- 0.5  # 0.5 in Campbell and Yau 2018; to select highest beta for GO enrichment

meanExprThresh_limma <- 0.5

startTime <-  Sys.time()

source("phenopath_my_utils.R")


betaU <-"\u03B2" 
lambdaU <- "\u03BB"
chiU <- "\u03C7"

runPheno <- T
runNorm <- T


plotType <- "png"
myHeight <- myWidth <- 400
myHeightGG <- 6
myWidthGG <- 6

# filter to keep only the highly variable genes
# Campbell 2018: 
# for BRAC, whose variance in log(TPM+1) expression was greater than 1 and whose median absolute deviation was greater than 0
var_thresh <- 1
mad_thresh <- 0


purityFilter <- 0.6

inFolder <- file.path("..","tcga_data","DOWNLOAD_TCGA_BRCA_RECOUNT2")

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

outFolder <- "PHENOPATH_RECOUNT2_TCGA_BRCA"
dir.create(outFolder, recursive=TRUE)


# load input data
brca_data_raw <- get(load(file.path(inFolder, paste0("all_counts_onlyPF_", purityFilter, ".Rdata"))))
dim(brca_data_raw)
# 25526   893

stopifnot(brca_data_raw >= 0)


tcga_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_BRCA_RECOUNT2/tcga_sampleAnnot.Rdata"))

# 1- gc content normalization
# xx=rownames(brca_data_raw)
# xx=xx[!grepl("_PAR_Y$", xx)]
# stopifnot(!duplicated(gsub("\\..*", "",xx)))
## I have some weird ENSGxx_PAR_Y genes -> remove them

## check the weird PAR_Y genes
tocheckgenes <- rownames(brca_data_raw)[grepl("_PAR_Y$", rownames(brca_data_raw))]
mygene=tocheckgenes[1]
for(mygene in tocheckgenes) {
  othergene <- gsub("_PAR_Y", "", mygene)
  if(!othergene %in% rownames(brca_data_raw)) next
  tmp_dt <- data.frame(brca_data_raw[mygene,],
                       brca_data_raw[othergene,])
  colnames(tmp_dt) <- c(mygene, othergene)
  tmp_dt <- log2(tmp_dt+1)
  outFile <- file.path(outFolder, paste0("check_dupgenes_", mygene, "_", othergene, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(tmp_dt)
  dev.off()
  cat(paste0("written ... ", outFile, "\n"))  
}

httr::set_config(httr::config(ssl_verifypeer = FALSE))  ### added to access ensembl biomart connection
tokeep <- ! grepl("_PAR_Y$", rownames(brca_data_raw))
cat(paste0("... remove weird genes: ", sum(!tokeep), "/", length(tokeep), "\n"))
brca_data_raw <- brca_data_raw[tokeep,]
stopifnot(!duplicated(gsub("\\..*", "",rownames(brca_data_raw))))
rownames(brca_data_raw)<-gsub("\\..*", "",rownames(brca_data_raw))
              # get list of all protein-coding genes
              mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
              listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                                   filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
              head(listOfGenes)
              dim(listOfGenes)

              brca_data_forNorm <- subset(brca_data_raw,rownames(brca_data_raw) %in% listOfGenes$ensembl_gene_id)
              dim(brca_data_forNorm)
# [1] 19177   533

### normalize as in Lucchetta et al. 2019
if(runNorm) {
  brca_dat_gcNorm <- TCGAanalyze_Normalization(tabDF = brca_data_forNorm,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")
  outFile <- file.path(outFolder, "brca_data_gcNorm.Rdata")
  save(brca_dat_gcNorm, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "brca_data_gcNorm.Rdata")
  brca_dat_gcNorm <- get(load(outFile))
}

brca_data_gcNorm_log <- log2(brca_dat_gcNorm + 1)

nTCGA <- sum(grepl("^TCGA", colnames(brca_data_gcNorm_log)))
stopifnot(nTCGA == ncol(brca_data_gcNorm_log))
ERpos_samples <- colnames(brca_data_gcNorm_log)[colnames(brca_data_gcNorm_log) %in%
                                                  tcga_annot_dt$cgc_sample_id[tcga_annot_dt$ER_status == "Positive"]]
nERpos <- length(ERpos_samples)

ERneg_samples <- colnames(brca_data_gcNorm_log)[colnames(brca_data_gcNorm_log) %in%
                                                  tcga_annot_dt$cgc_sample_id[tcga_annot_dt$ER_status == "Negative"]]
nERneg <- length(ERneg_samples)

stopifnot(colnames(brca_data_gcNorm_log) == c(ERpos_samples, ERneg_samples))


################################### 
################################### PCA ON not norm DATA -> for comparison
################################### 
brca_data_notNorm_log <- log2(brca_data_forNorm + 1)

brca_data_notNorm_log_no0 <- brca_data_notNorm_log[rowSums(brca_data_notNorm_log) > 0,]

pca_brca <- prcomp(t(brca_data_notNorm_log_no0), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nTCGA)

stopifnot(rownames(pca_brca_lowrepr) == tcga_annot_dt$cgc_sample_id)

mycolvect <- 1+as.numeric(tcga_annot_dt$ER_status == "Positive")

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_brca_lowrepr))))
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

### here I see a cluster
# but not on the norm data.. ?! 




################################### 
################################### PCA ON  norm but not filt DATA -> for comparison
################################### 

brca_data_gcNorm_log_no0 <- brca_data_gcNorm_log[rowSums(brca_data_gcNorm_log) > 0,]

pca_brca <- prcomp(t(brca_data_gcNorm_log_no0), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nTCGA)

stopifnot(rownames(pca_brca_lowrepr) == tcga_annot_dt$cgc_sample_id)

mycolvect <- 1+as.numeric(tcga_annot_dt$ER_status == "Positive")

outFile <- file.path(outFolder, paste0("in_norm_data_pca_12_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_^norm_pca_23_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_brca_lowrepr))))
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

#### ADD SOME FILTERING HERE

# For input to PhenoPath we used the 4,579 genes whose variance in log(TPM+1)
# expression was greater than 1 and whose median absolute deviation was greater than 0

var_exprs <- rowVars(brca_data_gcNorm_log)
mad_exprs <- rowMads(brca_data_gcNorm_log)
to_keep <-  var_exprs > var_thresh & mad_exprs > mad_thresh

brca_data_gcNorm_log <- brca_data_gcNorm_log[to_keep,]

cat(paste0("... keep highly variable genes ", nrow(brca_data_gcNorm_log), "/", length(to_keep), "\n"))

################################### 
################################### PCA ON NORM + filtered DATA 
################################### 

pca_brca <- prcomp(t(brca_data_gcNorm_log), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nTCGA)

stopifnot(rownames(pca_brca_lowrepr) == tcga_annot_dt$cgc_sample_id)

mycolvect <- 1+as.numeric(tcga_annot_dt$ER_status == "Positive")

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_12_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_23_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_brca_lowrepr))))
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))



stop("--ok\n")
# ################################## remove outlier with mclust [not needed with the norm  + filt ???]





We need some QC to remove outlier cells:
  
  ```{r qc}
sce <- plotPCA(sce, ntop = nrow(sce), return_SCESet = TRUE, ncomponents = 3)
set.seed(123L)
mc <- Mclust(redDim(sce)[,c(1,3)], G = 2)
sce$cluster <- mc$classification
sce$Cluster <- factor(mc$classification)
plotReducedDim(sce, colour_by = 'Cluster', ncomponents = 3)
pca_plot <- last_plot()
saveRDS(pca_plot, file = "../../data/BRCA/brca_pca_plot.rds")
to_keep_index <- which.max(table(sce$cluster))
samples_to_keep <- which(sce$cluster == to_keep_index)
# samples_to_keep <- sce$pct_exprs_top_100_features < 2.4
sce <- sce[, samples_to_keep]
sce_gene <- sce_gene[, samples_to_keep]
```




rm(gtex_data_log)
rm(tcga_data_log)

# stop("--ok\n")

########################### KEEP ONLY VARIABLE GENES
# median absolute deviation in log(TPM+1)
# filter to keep only the highly variable genes
### this should be done after having checked for outliers...
nGenes <- nrow(brca_data_gcNorm_log)
all_mad <- apply(brca_data_gcNorm_log, 1, mad)
to_keep <- all_mad > mad_thresh
stopifnot(length(to_keep) == nGenes)
cat(paste0(sum(to_keep), "/", length(to_keep), "\n"))

plot(density(all_mad), main="all gene MADs (GTEX+TCGA)")
abline(v = mad_thresh, col="red")
legend("topright",legend=c(
  paste0("threshold = ", round(mad_thresh, 2)),
  paste0("to keep: ", sum(to_keep), "/", length(to_keep), "\n")), bty="n")

stopifnot(length(to_keep) == nrow(brca_data_gcNorm_log))
stopifnot(nGenes == nrow(brca_data_gcNorm_log))

brca_data_filtered <- brca_data_gcNorm_log[to_keep,]
stopifnot(nrow(brca_data_filtered) == sum(to_keep))
nGenes <- sum(to_keep)

outFile <- file.path(outFolder, paste0("brca_data_filtered.Rdata"))
save(brca_data_filtered, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


### look how different are the TCGA and the GTEX data
##### distribution of the counts - before MAD filter
gtex_density_all <- density(brca_data_gcNorm_log[,grep("^GTEX", colnames(brca_data_gcNorm_log))])
tcga_density_all <- density(brca_data_gcNorm_log[,grep("^TCGA", colnames(brca_data_gcNorm_log))])

# dev.off()
outFile <- file.path(outFolder, paste0("all_counts_log_ERpos_ERneg_density.", plotType))
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
gtex_density_madF <- density(brca_data_filtered[,grep("^GTEX", colnames(brca_data_filtered))])
tcga_density_madF <- density(brca_data_filtered[,grep("^TCGA", colnames(brca_data_filtered))])

# dev.off()
outFile <- file.path(outFolder, paste0("MADfiltered_counts_log_ERpos_ERneg_density.", plotType))
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
gtex_density_all <- density(apply(brca_data_gcNorm_log[,grep("^GTEX", colnames(brca_data_gcNorm_log))], 1, mad))
tcga_density_all <- density(apply(brca_data_gcNorm_log[,grep("^TCGA", colnames(brca_data_gcNorm_log))], 1, mad))

# dev.off()
outFile <- file.path(outFolder, paste0("all_mads_ERpos_ERneg_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(gtex_density_all, main="data density (log2(.+1))",
     xlim=range(c(gtex_density_all$x, tcga_density_all$x)),
     ylim=range(c(gtex_density_all$y, tcga_density_all$y)))
lines(tcga_density_all, col="red")
mtext(side=3, text="(all)")
legend("topright", legend=c("GTEX", "TCGA"),
       bty="n",
       col=c("black", "red"), lty=c(1,1))
dev.off()
cat(paste0("... written: ", outFile, "\n"))

# stop("OK\n")

##### distribution of the mads - after filtering
gtex_density_madF <- density(apply(brca_data_filtered[,grep("^GTEX", colnames(brca_data_filtered))], 1, mad))
tcga_density_madF <- density(apply(brca_data_filtered[,grep("^TCGA", colnames(brca_data_filtered))], 1, mad))

# dev.off()
outFile <- file.path(outFolder, paste0("MADfiltered_mads_ERpos_ERneg_density.", plotType))
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
# brca_data_filtered <- get(load("PHENOPATH_RECOUNT2_TCGA_brca/brca_data_filtered.Rdata"))

# phenopath tuto: sim$y is the N×G matrix of gene expression (N=100 cells and G=40 genes)
# so I have to take the transpose of my brca_data_log_matrix to have samples in line and genes in columns
brca_data_filteredT <- t(brca_data_filtered)
stopifnot(nGenes == ncol(brca_data_filteredT))
stopifnot(nGTEX+nTCGA == nrow(brca_data_filteredT))

pca_brca <- prcomp(brca_data_filteredT, scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nGTEX+nTCGA)

pc1 <- pca_brca_lowrepr[,1]
pc2 <- pca_brca_lowrepr[,2]
pc3 <- pca_brca_lowrepr[,3]

outFile <- file.path(outFolder, paste0("in_data_pca_12_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()

pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA+GTEX BRCA (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_brca_lowrepr))))
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_data_pca_23_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA+GTEX BRCA (log2(.+1))",
        col=1+as.numeric(grepl("TCGA", rownames(pca_brca_lowrepr))))
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################################################
############################################################## LET'S DO PHENO
##############################################################

### how to set the covar ?? explained here: https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
# # 4.2 Inference with PhenoPath
# First we must decide how to pass in the covariate information (ie the stimulant applied) to the software as the x
# values. Here we will give cells exposed to LPS a value of 1 and cells exposed to PAM a value of -1.
# This means the overall pathway loading λ is the average change for LPS and PAM cells, 
# while if the β parameter is positive it means the gene is more upregulated over pseudotime under LPS and 
# if β is negative it means the gene is more upregulated under PAM11 
# Instead we could encode LPS to 1 and PAM to 0, in which case the pathway loading λ would be the change under PAM and λ+β
# the change under LPS stimulation..

# so here give a value of 1 in tumor and a value of -1 in normal
# This means the overall pathway loading λ is the average change for normal and tumor
# if β > 0 =>  the gene is more upregulated over pseudotime under tumor  
# if β  < 0 =>   the gene is more upregulated under normal 


# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(brca_data_filteredT) == nGTEX+nTCGA)
stopifnot(ncol(brca_data_filteredT) == nGenes)
mycovar <- 2 * grepl("^TCGA", rownames(brca_data_filteredT)) - 1
stopifnot(sum(mycovar== -1) == nGTEX)
stopifnot(sum(mycovar== 1) == nTCGA)
stopifnot(!is.na(mycovar))

if(runPheno) {
  brca_phenopath_fit <- phenopath(brca_data_filteredT, mycovar, elbo_tol = 1e-6, thin = 40)
  
  outFile <- file.path(outFolder, paste0("brca_phenopath_fit.Rdata"))
  save(brca_phenopath_fit, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("brca_phenopath_fit.Rdata"))
  brca_phenopath_fit <- get(load(outFile))
}
#   brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_brca/brca_phenopath_fit.Rdata"))


######################################################## 
##############  PHENOPATH RESULT ANLAYSIS
######################################################## 

####### check the model
# it is important to check convergence with a call to plot_elbo(brca_phenopath_fit) to ensure the ELBO is approximately flat:

outFile <- file.path(outFolder, paste0("phenopath_elbo_fit.", plotType))
p <- plot_elbo(brca_phenopath_fit)
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

  
####### retrieve the pseudotimes
# maximum a-posteriori (MAP) estimates of the pseudotimes using the trajectory()

brca_pseudotimes <-  trajectory(brca_phenopath_fit)
 
############## correlation with PCs ##############
### look if the pseudotimes correlate with the PCs
# tumor=1, normal=-1

# dev.off()
outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc1, 
     y = brca_pseudotimes,
     xlab = "PC1", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc1, brca_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc2.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc2, 
     y = brca_pseudotimes,
     xlab = "PC2", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc2, brca_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc3.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc3, 
     y = brca_pseudotimes,
     xlab = "PC3", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc3, brca_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

############## LOOK A BIT AT THE BETA VALUES ##############
# density of the beta with color code signifcance

gene_names <- colnames(brca_data_filteredT)

df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = gene_names)

outFile <- file.path(outFolder, paste0("beta_values_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(df_beta$beta))
lines(density(df_beta$beta[df_beta$is_sig]), col=2)
# lines(density(df_beta$beta[!df_beta$is_sig]), col=3)
legend("topright", legend=c("all dist.", "only signif. dist."), lty=1, col=c(1,2), bty="n")
mtext(side=3, text=paste0("# signif: ", sum(df_beta$is_sig), "/", nrow(df_beta)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


############## top significant interactions pseudotime x covar ##############

# look at the top bottom and up interaction effects
nTop <- 10
  
# show the top genes with signif interactions
df_beta <- df_beta[order(df_beta$beta, decreasing = TRUE),]
topPosGenes <- df_beta[1:nTop,]
topPosGenes$dir <- "pos"
df_beta <- df_beta[order(df_beta$beta, decreasing = FALSE),]
topNegGenes <- df_beta[1:nTop,]
topNegGenes$dir <- "neg"

dfBeta_topGenes <- rbind(topPosGenes, topNegGenes)
stopifnot(!duplicated(dfBeta_topGenes$gene))
dfBeta_topGenes$gene <- as.character(dfBeta_topGenes$gene)

gene_lab_dt <- get(load(file.path(inFolder, "out_intersect_dt.RData")))
gene_lab_dt$geneID_short <-gsub("\\..*", "",gene_lab_dt$geneID)
gene_lab_dt <- unique(gene_lab_dt[,c("geneSymb", "geneID_short")])
stopifnot(!duplicated(gene_lab_dt$geneID_short))
ens2genes <- setNames(gene_lab_dt$geneSymb, gene_lab_dt$geneID_short)
stopifnot(dfBeta_topGenes$gene %in% gene_lab_dt$geneID_short)
stopifnot(dfBeta_topGenes$gene %in% names(ens2genes))
dfBeta_topGenes$geneSymb <- ens2genes[dfBeta_topGenes$gene]
stopifnot(!duplicated(dfBeta_topGenes$geneSymb))
dfBeta_topGenes$geneSymb <- factor(dfBeta_topGenes$geneSymb, levels=dfBeta_topGenes$geneSymb)
dfBeta_topGenes$dir <- factor(dfBeta_topGenes$dir, levels=c("pos", "neg"))
  
p <-  ggplot(dfBeta_topGenes, aes(x = geneSymb, y = beta, color = is_sig)) +
    geom_point() +
    facet_wrap(. ~ dir ,scales="free") +
    ggtitle(paste0("Genes with top ", betaU), subtitle=paste0("(", nTop, " lowest and hightest)"))+
    geom_errorbar(aes(ymin = beta - 2 * beta_sd, ymax = beta + 2 * beta_sd)) +
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylab(expression(beta)) +
    scale_color_brewer(palette = "Set2", name = "Significant")

outFile <- file.path(outFolder, paste0("genes_with_top_and_bottom_n", nTop, "_beta.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


############## top significant interactions pseudotime x covar ##############
# which gene has highest interaction effect ?

stopifnot(df_beta$gene %in% names(ens2genes))
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- ens2genes[df_beta$gene]

# the same for the lowest interaction effect ?

p <- plot_iGeneExpr(igene= which.max(df_beta$beta),
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
               valuedt=df_beta,
               valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.max(df_beta$beta),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



p <- plot_iGeneExpr(igene= which.min(df_beta$beta),
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
               valuedt=df_beta,
                           valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.min(df_beta$beta),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))





############## look trajectory and phenotypes (tumor stages, age, etc.) ##############

table(tcga_annot_dt$cgc_case_days_to_death)
table(tcga_annot_dt$cgc_slide_percent_tumor_cells)


stopifnot(rownames(brca_data_filteredT) == c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id))


all_traj <- setNames(brca_pseudotimes, rownames(brca_data_filteredT))
stopifnot(names(all_traj) == c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id))

############### TCGA data
stopifnot(tcga_annot_dt$gdc_cases.samples.sample_type == tcga_annot_dt$cgc_sample_sample_type)

# > cgc_case_age_at_diagnosis
# Error: object 'cgc_case_age_at_diagnosis' not found
# > cgc_case_days_to_death
# Error: object 'cgc_case_days_to_death' not found
# > cgc_case_clinical_stage
# > cgc_durg_therapy_drug_name
# # Error: object 'cgc_durg_therapy_drug_name' not found -> too many
# > cgc_follow_up_days_to_deathg
# Error: object 'cgc_follow_up_days_to_deathg' not found
# > histologic name and gdc_cases.project .name grep serous
# Error: unexpected symbol in "histologic name"
# > sampletypeID


myplotlab <- "tumor stage"
stg_ords <-c(  "Stage IC", "Stage IIA" ,"Stage IIB" , "Stage IIC","Stage IIIA", "Stage IIIB","Stage IIIC",   "Stage IV")
p <- plot_pheno_catego(tcga_annot_dt, plotvar= "cgc_case_clinical_stage", plotlab=myplotlab, varords=stg_ords)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# plot_pheno_catego(tcga_annot_dt, plotvar= "cgc_drug_therapy_drug_name", plotlab="drug therapy", varords=NULL)
## too many categories


myplotlab <- "year of birth"
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.demographic.year_of_birth", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "days to death"
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_case_days_to_death", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "diag. days to birth"
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.diagnoses.days_to_birth", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "age at diag."
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.diagnoses.age_at_diagnosis", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct normal cells"
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_normal_cells", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct neutrophil infilt."
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_neutrophil_infiltration", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct monocyte infilt."
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_monocyte_infiltration", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "slide pct lymphocyt infilt."
p <- plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_lymphocyte_infiltration", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))







############### GTEX data

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


p <- ggplot(gtex_traj_dt, aes(x= gtex_age, y= pseudotime) )+
  geom_boxplot(notch = F, outlier.shape=NA)+
  geom_jitter(aes(col=gtex_age),alpha=0.7,position=position_jitterdodge())+
  ggtitle(paste0("Pseudotime by age class"), subtitle = paste0("(GTEX data)"))+
  scale_color_nejm()+
  ylab("PP Pseudotimes")+
  xlab("GTEX donnor age class")+
  myG_theme +
  labs(color="")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank() )

outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_age_GTEX.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


############## further look at the signif. interactions ##############

# extract the interaction parameters 
# -  feature =>  The feature (usually gene)
# - covariate => The covariate, specified from the order originally supplied to the call to phenopath
# - interaction_effect_size => The effect size of the interaction (β
# - significant => Boolean for whether the interaction effect is significantly different from 0
# - chi => The precision of the ARD prior on β [Automatic Relevance Determination]
# - pathway_loading => The pathway loading λ

int_dt <- interactions(brca_phenopath_fit)
stopifnot(as.character(int_dt$feature) %in% names(ens2genes))
int_dt$featureSymb <- ens2genes[as.character(int_dt$feature) ]

# plot the posterior ARD variances (1/χ) against the posterior interaction effect sizes (β)
# colouring them by which are found to be significant and annotating the top few genes:

chi_cutoff <- sort(int_dt$chi)[10]


p <- ggplot(int_dt, aes(x = interaction_effect_size, y = 1 / chi,
                 color = significant_interaction)) +
  ggtitle("", subtitle=paste0("labs: ", chiU, " top10"))+
  xlab(paste0("posterior interaction effect sizes (", lambdaU, ")"))+
  ylab(paste0("posterior ARD variances (1/", chiU, ")"))+
  geom_point() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_text_repel(data = dplyr::filter(int_dt, chi < chi_cutoff),
                  aes(label = featureSymb)) +
  labs(color="significant")+
  scale_colour_brewer(palette = "Set1")+
  theme(plot.subtitle = element_text(hjust=0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))


outFile <- file.path(outFolder, paste0("posteriorARDchi_vs_posteriorEffectSizeBeta.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

# plot the “landscape” of interactions, where we plot the interaction effect size against the pathway score.
# The further up the y-axis a gene is, the more it is upregulated under tumor rather than normal (and vice-versa),
# while the further along the x-axis a gene is, the more it is upregulated over pseudotime regardless of the covariate condition

# y = covariate-pseudotime interaction (beta)
# x = gene regulation over pseudotime (pathway loading)
# x<0 & y>0 [left,top]: gene downreg. along pseudotime, cond=-1 increases downregulation
# x<0 & y<0 [left, bottom]: gene downreg. along pseudotime, cond=+1 increases downregulation
# x>0 & y>0 [right, top]: gene upreg. along pseudotime, cond=+1 increases upregulation
# x>0 & y<0 [right, bottom]: gene upreg. along pseudotime, cond=-1 increases upregulation

p <- ggplot(int_dt, aes(x = pathway_loading, y = interaction_effect_size,
                 color = significant_interaction)) +
  ggtitle("", subtitle=paste0("labs: ", chiU, " top10"))+
  geom_point() +
  ylab(paste0("posterior interaction effect sizes (", betaU, ")"))+
  xlab(paste0("pathway loading (", lambdaU, ")"))+
  geom_text_repel(data = dplyr::filter(int_dt, chi < chi_cutoff),
                  aes(label = featureSymb), size = 5) +
  scale_colour_brewer(palette = "Set1")  +
  labs(color="significant")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(plot.subtitle = element_text(hjust=0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

outFile <- file.path(outFolder, paste0("posteriorEffectSizeBeta_vs_pathwayloadingLambda.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

############## look at the gene with top and bottom pathway loading ##############

int_dt <- as.data.frame(int_dt)
tmp1 <- as.data.frame(int_dt)
tmp2 <- df_beta
tmp1 <- tmp1[order(tmp1$interaction_effect_size),]
tmp2 <- tmp2[order(tmp2$beta),]
stopifnot(tmp1$featureSymb == tmp2$geneSymb)


p <- plot_iGeneExpr(igene= which.min(int_dt$pathway_loading),
                                 exprdt=brca_data_filteredT,
                                 pseudot=brca_pseudotimes,
                                 covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
                                 valuedt=int_dt,
                                             valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.min(int_dt$pathway_loading),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr(igene= which.max(int_dt$pathway_loading),
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
               valuedt=int_dt,
               valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr_gg2(igene= which.max(int_dt$pathway_loading),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

############## top significant pathway loading ##############

# look at the top bottom and up interaction effects
nTop <- 10

# show the top genes with signif interactions
int_dt <- as.data.frame(int_dt)
int_dt <- int_dt[order(int_dt$pathway_loading, decreasing = TRUE),]
topPosGenes <- int_dt[1:nTop,]
topPosGenes$dir <- "pos"
int_dt <- int_dt[order(int_dt$pathway_loading, decreasing = FALSE),]
topNegGenes <- int_dt[1:nTop,]
topNegGenes$dir <- "neg"

dfLambda_topGenes <- rbind(topPosGenes, topNegGenes)
stopifnot(!duplicated(dfLambda_topGenes$feature))
dfLambda_topGenes$feature <- as.character(dfLambda_topGenes$feature)


stopifnot(dfLambda_topGenes$feature %in% gene_lab_dt$geneID_short)
stopifnot(dfLambda_topGenes$feature %in% names(ens2genes))
dfLambda_topGenes$featureSymb <- ens2genes[dfLambda_topGenes$feature]
stopifnot(!duplicated(dfLambda_topGenes$featureSymb))
dfLambda_topGenes$featureSymb <- factor(dfLambda_topGenes$featureSymb, levels=dfLambda_topGenes$featureSymb)
dfLambda_topGenes$dir <- factor(dfLambda_topGenes$dir, levels=c("pos", "neg"))

p <-  ggplot(dfLambda_topGenes, aes(x = featureSymb, y = pathway_loading, color = significant_interaction)) + 
  geom_point() +
  facet_wrap(. ~ dir ,scales="free") +
  ggtitle(paste0("Genes with top pathway loading (", lambdaU, ")"), subtitle=paste0("(", nTop, " lowest and highest)"))+
  # geom_errorbar(aes(ymin = beta - 2 * beta_sd, ymax = beta + 2 * beta_sd)) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab(expression(beta)) +
  scale_color_brewer(palette = "Set2", name = "Significant")

outFile <- file.path(outFolder, paste0("genes_with_top_and_bottom_n", nTop, "_lambda.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

############## PC dots with color-coded by pseudotime gradient ##############

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_brca),
            condvect = gsub("(^.+?)-.+", "\\1", rownames(pca_brca_lowrepr)),
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA+GTEX BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_pseudotimeGrad_ERpos_ERneg.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,2), summ_dt=summary(pca_brca),
            condvect = gsub("(^.+?)-.+", "\\1", rownames(pca_brca_lowrepr)),
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA+GTEX BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_pseudotimeGrad_ERpos_ERneg.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

############## gene set enrichment on top and bottom beta value genes ##############


# A GO enrichment analysis of the genes upregulated along pseudotime whose upregulation was increased by LPS stimulation
# showed enrichment for immune system processes.
# do GO enrichment analysis of the top upregulated along pseudotime for cond=+1 et cond=-1 separately
# GO category vs. p-value with colorcond by condvar

# compare to common c 
# A GO enrichment analysis of upregulated genes confirms the latent trajectory encodes immune pathway activation in each tumour

# in the case of COAD -> was common trajectory
# in the case of LS/PAM -> was different

# I think the following would be of interest
# the highest/lowest betas => this would give the ones with highest interaction effects up/down regulated normal/tumor
# the highest/lowest lambda => the one most up/down regulated across pseudotimes


# What is pseudotime?  - GO association
# from: https://github.com/kieranrcampbell/phenopath_analyses/blob/master/analysis/BRCA/clvm_analysis_er_pos.Rmd


pp_input_data <- brca_data_filteredT

all_genes <- colnames(pp_input_data) # ensemblID
corr_expr_pt <- apply(pp_input_data, 2, cor, brca_pseudotimes)
cdf <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt)
stopifnot(!is.na(cdf))

outFile <- file.path(outFolder, paste0("corr_gene_expr_pseudot_distribution_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cdf$correlation), main="corr. gene expr. and pseudotimes")
mtext(side=3, text="(all samples; # = ", nrow(cdf), ")")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

stopifnot(length(brca_pseudotimes) == nrow(pp_input_data))

tcga_samp <- grep("^TCGA", rownames(pp_input_data))
stopifnot(length(tcga_samp) == nTCGA)
corr_expr_pt_tcga <- apply(pp_input_data[tcga_samp,], 2, cor, brca_pseudotimes[tcga_samp])
cdf_tcga <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt_tcga)
stopifnot(!is.na(cdf_tcga))

gtex_samp <- grep("^GTEX", rownames(pp_input_data))
stopifnot(length(gtex_samp) == nGTEX)
corr_expr_pt_gtex <- apply(pp_input_data[gtex_samp,], 2, cor, brca_pseudotimes[gtex_samp])
to_keep <- which(!is.na(corr_expr_pt_gtex))  ### there is 2 genes with all 0 values !!! -> cannot compute corr
cdf_gtex <- data.frame(feature = all_genes[to_keep], 
                       correlation = corr_expr_pt_gtex[to_keep])
stopifnot(!is.na(cdf_gtex))


outFile <- file.path(outFolder, paste0("corr_gene_expr_pseudot_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cdf_tcga$correlation), main="corr. gene expr. and pseudotimes", col=2)
lines(density(cdf_gtex$correlation), col=1)
legend("topright", legend=c("GTEX", "TCGA"), pch=16, col=c(1,2), bty="n")
mtext(side=3, text=paste0("(# GTEX=", nGTEX, "; # TCGA=", nTCGA, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


stopifnot(topCorrThresh >= 0)
upreg_genes <- cdf$feature[cdf$correlation > topCorrThresh]
downreg_genes <- cdf$feature[cdf$correlation < -topCorrThresh]
upg <- 1 * (all_genes %in% upreg_genes)
downg <- 1 * (all_genes %in% downreg_genes)
names(upg) <- names(downg) <- all_genes
pwfup <- nullp(upg, genome, id) # nullp = Probability Weighting Function from goseq
goup <- goseq(pwfup, genome, id, test.cats = go_cat)
pwfdown <- nullp(downg, genome, id)
godown <- goseq(pwfdown, genome, id, test.cats = go_cat)


corrExpr_gos <- rbind(
  parse_go(goup, "Up-regulated", length(upreg_genes), plotntop=plotNgos),
  parse_go(godown, "Down-regulated", length(downreg_genes), plotntop=plotNgos)
)

corrExpr_gos$term <- factor(corrExpr_gos$term, levels = corrExpr_gos$term[order(corrExpr_gos$log10qval)])

p <- ggplot(corrExpr_gos, aes(x = term, y = log10qval)) +
  geom_point() +
  ggtitle(paste0("Top ", plotNgos, " GOs"), subtitle = paste0("(abs. corr. with Pseudotime > ", topCorrThresh, ")"))+
  facet_wrap(~ type, scales = "free", nrow = 2) +
  coord_flip() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  ylab(expression(paste("-", log[10], " q-value"))) +
  theme(
    plot.title = element_text(hjust=0.5),
    plot.subtitle = element_text(hjust=0.5),
    legend.position = c(0, 0), legend.direction = "horizontal",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)) 

outFile <- file.path(outFolder, paste0("topCorrExpr_topGOs.", plotType))
ggsave(p, filename = outFile, height=myHeightGG*1.4, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


  
### Do the same for top beta:

stopifnot(topBetaThresh >= 0)

df_beta$ensembl_gene_id <- df_beta$gene
upreg_genes <- df_beta$ensembl_gene_id[df_beta$beta > topBetaThresh & df_beta$is_sig]
downreg_genes <- df_beta$ensembl_gene_id[df_beta$beta < -topBetaThresh & df_beta$is_sig]

upg <- 1 * (all_genes %in% upreg_genes)
downg <- 1 * (all_genes %in% downreg_genes)
names(upg) <- names(downg) <- all_genes
pwfup <- nullp(upg, genome, id)
goup <- goseq(pwfup, genome, id, test.cats = go_cat)
pwfdown <- nullp(downg, genome, id)
godown <- goseq(pwfdown, genome, id, test.cats = go_cat)

betacorr_gos <- rbind(
  parse_go(goup, "Up-regulated", length(upreg_genes), plotNgos),
  parse_go(godown, "Down-regulated", length(downreg_genes), plotNgos)
)


betacorr_gos$term <- factor(betacorr_gos$term, levels = betacorr_gos$term[order(betacorr_gos$log10qval)])

p <- ggplot(betacorr_gos, aes(x = term, y = log10qval, color = type)) +
  ggtitle(paste0("Top ", plotNgos," GOs"), 
          subtitle = paste0("abs(", betaU, ")>", topBetaThresh, 
                            " & signif. (#up=",length(upreg_genes) ,"; #down=", length(downreg_genes) , ")"))+
  geom_point() +
  facet_wrap(~ type, scales = "free", nrow = 2) +
  coord_flip() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  ylab(expression(paste(log[10], " q-value"))) +
  theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 11))

outFile <- file.path(outFolder, paste0("topBeta_topGOs.", plotType))
ggsave(p, filename = outFile, height=myHeightGG*1.4, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


# for reactome pathway analysis
# https://github.com/kieranrcampbell/phenopath_revisions/blob/master/analysis/brca_reanalysis/clvm_analysis.Rmd

############## complementary info to std gene DE ##############

# look at genes that are differentially expressed (large difference of gene expression normal vs. tumor)
# and that have different trajectories (high interaction effects)
# Whilst many of these 92 genes are differentially expressed between MSI groups, including MLH2 and TGFBR2 (Fig. 5c), 
# PhenoPath is able to resolve the dynamic contribution to these expression differences


# see e.g. http://research.libd.org/recountWorkflow/articles/recount-workflow.html
# brca_data_raw comes from here:
# ovary_gtex <- scale_counts(TCGAquery_recount2(project="gtex", tissue = "ovary")$gtex_brcaary)
# ovary_tcga <- scale_counts(TCGAquery_recount2(project="tcga", tissue = "ovary")$tcga_brcaary)
# # I skip some steps of processing, but I have in the end 
# # gene x sample matrix with GTEX samples starting with GTEX- and TCGA samples with TCGA-
# ovary_all_data <- cbind(assays(ovary_gtex)$counts, assays(ovary_tcga)$counts)

brca_data_raw <- get(load(file.path(inFolder, paste0("all_counts_onlyPF_", purityFilter, ".Rdata"))))

# filter because I have removed the outlier !
stopifnot(tcga_annot_dt$cgc_sample_id %in% colnames(brca_data_raw))
stopifnot(gtex_annot_dt$sampid %in% colnames(brca_data_raw))
stopifnot(!duplicated(tcga_annot_dt$cgc_sample_id))
stopifnot(!duplicated(gtex_annot_dt$sampid))
brca_data_raw <- brca_data_raw[, c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id) ]

# remove the weird gene
tokeep <- ! grepl("_PAR_Y$", rownames(brca_data_raw))
cat(paste0("... remove weird genes: ", sum(!tokeep), "/", length(tokeep), "\n"))
brca_data_raw <- brca_data_raw[tokeep,]
stopifnot(!duplicated(gsub("\\..*", "",rownames(brca_data_raw))))


## Extract counts and filter out lowly expressed geens

to_keep <- rowMeans(brca_data_raw) > meanExprThresh_limma
cat(paste0("... keep sufficient expr: ", sum(to_keep), "/", length(tokeep), "\n"))


## Build DGEList object
dge <- DGEList(counts = brca_data_raw[to_keep, ])
## Calculate normalization factors
dge <- calcNormFactors(dge)
# plotMDS(dge)
samples_groups <- c(rep("normal", nrow(gtex_annot_dt)), rep("tumor", nrow(tcga_annot_dt)))
my_group_design <- factor(samples_groups, levels = c("normal", "tumor"))
my_design <- model.matrix( ~ my_group_design)
## Run voom
v <- voom(dge, my_design, plot = TRUE)
## Run remaining parts of the DE analysis
fit <- lmFit(v, my_design)
efit <- eBayes(fit)

## Visually explore DE results
# limma::plotMA(efit)
# limma::volcanoplot(efit)

DE_topTable <- topTable(efit, coef=ncol(v$design), number=Inf, sort.by="p") ## if not 0+ in design -> coef=2


stopifnot(!duplicated(rownames(gsub("\\..*", "",rownames(brca_data_raw)))))

DE_topTable$ensemblID <- gsub("\\..*", "",rownames(DE_topTable))
stopifnot(!duplicated(DE_topTable$ensemblID))

## look for genes that has high logFC + has high beta value

nTopDEvalueThresh <- 800
topDEvalueThresh <- sort(abs(DE_topTable$logFC), decreasing = TRUE)[nTopDEvalueThresh]
topDEvalueThresh

topLogFC_genes <- DE_topTable$ensemblID[abs(DE_topTable$logFC) >= topDEvalueThresh]

nTopBetaValueThresh <- 800
topBetaValueThresh <- sort(abs(df_beta$beta), decreasing = TRUE)[nTopBetaValueThresh]
topBetaValueThresh
topBeta_genes <- df_beta$gene[abs(df_beta$beta) >= topBetaValueThresh]

intersect(topLogFC_genes, topBeta_genes)

# not convincing -> up to 500, no intersect

############## gene limma voom pvalue vs. phenopath beta value  ##############

df_beta$ensemblID <- df_beta$gene

merged_dt <- merge(df_beta, DE_topTable, by="ensemblID", all.x = TRUE, all.y=FALSE)
stopifnot(!is.na(merged_dt))

outFile <- file.path(outFolder, "merged_dt.Rdata")
save(merged_dt, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "DE_topTable.Rdata")
save(DE_topTable, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("limma_adjPval_vs_phenopath_beta.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  y=-log10(merged_dt$adj.P.Val),
  x=df_beta$beta,
  pch=16,
  col=2+2*merged_dt$is_sig,
  cex.lab=1.2,
  cex.axis=1.2,
  ylab="limma adj. Pval [-log10]",
  xlab=paste0("Phenopath ", betaU)
)
legend("topright", legend=c("signif.", "not signif."),
       pch=16,bty="n",
       col = c(2+2*c(TRUE, FALSE)))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



outFile <- file.path(outFolder, paste0("limma_logFC_vs_phenopath_beta.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  y=merged_dt$logFC,
  x=df_beta$beta,
  pch=16,
  col=2+2*merged_dt$is_sig,
  cex.lab=1.2,
  cex.axis=1.2,
  ylab="limma logFC",
  xlab=paste0("Phenopath ", betaU)
)
legend("topright", legend=c("signif.", "not signif."),
       pch=16,bty="n",
       col = c(2+2*c(TRUE, FALSE)))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))
stop("--ok\n")

####################################################################################################
### THRASH
####################################################################################################

