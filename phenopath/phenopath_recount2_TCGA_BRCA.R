
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
require(mclust)
require(ggExtra)

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

runPheno <- F
runNorm <- FALSE


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


tcga_annot_dt <- get(load(file.path(inFolder, "tcga_sampleAnnot.Rdata")))

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
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_13_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=mycolvect)
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

outFile <- file.path(outFolder, paste0("in_norm_pca_23_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_norm_pca_13_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
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

stopifnot(colnames(brca_data_gcNorm_log) == c(ERpos_samples, ERneg_samples))

var_exprs <- setNames(var_exprs,colnames(brca_data_gcNorm_log)) 
mad_exprs <- setNames(mad_exprs,colnames(brca_data_gcNorm_log)) 

# dev.off()
outFile <- file.path(outFolder, paste0("all_vars_ERpos_ERneg_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(var_exprs), main="all vars (log2(.+1))")
mtext(side=3, text="(all)")
abline(v=var_thresh, col="red")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("all_mads_ERpos_ERneg_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(mad_exprs), main="all mads (log2(.+1))")
mtext(side=3, text="(all)")
abline(v=var_thresh, col="red")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

nGenes <- sum(to_keep)

brca_data_gcNorm_log <- brca_data_gcNorm_log[to_keep,]
cat(paste0("... keep highly variable genes ", nrow(brca_data_gcNorm_log), "/", length(to_keep), "\n"))

brca_data_filtered <- brca_data_gcNorm_log
outFile <- file.path(outFolder, paste0("brca_data_filtered.Rdata"))
save(brca_data_filtered, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


################################### 
################################### PCA ON NORM + filtered DATA 
################################### 

pca_brca <- prcomp(t(brca_data_filtered), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nTCGA)

pc1 <- pca_brca_lowrepr[,1]
pc2 <- pca_brca_lowrepr[,2]
pc3 <- pca_brca_lowrepr[,3]

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
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_13_ERpos_ERneg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))


# ################################## remove outlier with mclust [not needed with the norm  + filt ???]



mc <- Mclust(pca_brca_lowrepr[,c(1,3)], G = 2)

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_13_ERpos_ERneg_withClust.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mc$classification)
mtext(text=paste0("nERneg=",nERneg, "; nERpos=",nERpos), side=3)
legend("topleft", legend=c("clust1", "clust2"), pch=16, col=c(1,2), bty="n")
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

# so here give a value of 1 in ER+ and a value of -1 in ER-
# This means the overall pathway loading λ is the average change for normal and tumor
# if β > 0 =>  the gene is more upregulated over pseudotime under tumor  
# if β  < 0 =>   the gene is more upregulated under normal 

brca_data_filteredT <- t(brca_data_filtered)
stopifnot(nGenes == ncol(brca_data_filteredT))
stopifnot(nERpos+nERneg == nrow(brca_data_filteredT))

# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(brca_data_filteredT) == nERneg+nERpos)
stopifnot(ncol(brca_data_filteredT) == nGenes)
mycovar <- 2 * ( rownames(brca_data_filteredT) %in% ERpos_samples) - 1
stopifnot(sum(mycovar== -1) == nERneg)
stopifnot(sum(mycovar== 1) == nERpos)
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
add_corr(pc1, brca_pseudotimes, cond_plus1="ER+", cond_minus1="ER-")
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
add_corr(pc2, brca_pseudotimes, cond_plus1="ER+", cond_minus1="ER-")
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
add_corr(pc3, brca_pseudotimes, cond_plus1="ER+", cond_minus1="ER-")
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
plot(density(df_beta$beta[df_beta$is_sig]), col=2)
lines(density(df_beta$beta), col=1)
# lines(density(df_beta$beta[!df_beta$is_sig]), col=3)
legend("topright", legend=c("all dist.", "only signif. dist."), lty=1, col=c(1,2), bty="n")
mtext(side=3, text=paste0("# signif: ", sum(df_beta$is_sig), "/", nrow(df_beta)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


############## top significant interactions pseudotime x covar ##############

# look at the top bottom and up interaction effects
nTop <- 10
  
# show the top genes with signif interactions
tmp_beta <- df_beta[order(df_beta$beta, decreasing = TRUE),]
topPosGenes <- tmp_beta[1:nTop,]
topPosGenes$dir <- "pos"
tmp_beta <- df_beta[order(df_beta$beta, decreasing = FALSE),]
topNegGenes <- tmp_beta[1:nTop,]
topNegGenes$dir <- "neg"

dfBeta_topGenes <- rbind(topPosGenes, topNegGenes)
stopifnot(!duplicated(dfBeta_topGenes$gene))
dfBeta_topGenes$gene <- as.character(dfBeta_topGenes$gene)

gene_lab_dt <- get(load(file.path(inFolder, "gene_dt.Rdata")))
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
               covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
               valuedt=df_beta,
               valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.max(df_beta$beta),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



p <- plot_iGeneExpr(igene= which.min(df_beta$beta),
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
               valuedt=df_beta,
                           valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.min(df_beta$beta),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


# plot selected genes from Campbell and Yau 2018
selected_symbs <- c("ESR1", "FOXC1", "FBP1")
gs="FOXC1"
stopifnot(df_beta$gene == colnames(brca_data_filteredT))
stopifnot(df_beta$gene %in% names(ens2genes))
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- ens2genes[df_beta$gene]

for(gs in selected_symbs) {
  if(!gs %in% ens2genes) {
    cat(paste0("warning: ", gs, " not available\n"))
    next
  }
  i_gs <- which(df_beta$geneSymb == gs)
  ens_gs <- df_beta$gene[i_gs]
  stopifnot(length(i_gs) == 1)
  tmp_dt <- df_beta
  tmp_dt <- tmp_dt[order(abs(tmp_dt$beta), decreasing = TRUE),]
  
  gs_beta_rank <- which(tmp_dt$gene == ens_gs)
  stopifnot(length(gs_beta_rank) == 1)
  gs_signif <- as.character(df_beta$is_sig[i_gs])
  gs_beta <- round(df_beta$beta[i_gs], 3)
  
  p <- plot_iGeneExpr_gg2(igene= i_gs,
                          exprdt=brca_data_filteredT,
                          pseudot=brca_pseudotimes,
                          covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
                          valuedt=df_beta,
                          valuecol="beta", symbcol="geneSymb", 
                          subtit=paste0(betaU, " rank=", gs_beta_rank, "; signif=",gs_signif ))
  
  outFile <- file.path(outFolder, paste0(gs, "_geneExpr_along_pseudotime_gg2.", plotType))
  ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
  cat(paste0("... written: ", outFile, "\n"))
  
}




############## look trajectory and phenotypes (tumor stages, age, etc.) ##############

table(tcga_annot_dt$cgc_case_days_to_death)
table(tcga_annot_dt$cgc_slide_percent_tumor_cells)


stopifnot(rownames(brca_data_filteredT) == c(tcga_annot_dt$cgc_sample_id))


all_traj <- setNames(brca_pseudotimes, rownames(brca_data_filteredT))
stopifnot(names(all_traj) == c(tcga_annot_dt$cgc_sample_id))

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


table(tcga_annot_dt[,"cgc_case_histological_diagnosis"])
# Infiltrating Carcinoma NOS    Infiltrating Ductal Carcinoma   Infiltrating Lobular Carcinoma 
# 1                              896                              209 
# Medullary Carcinoma            Metaplastic Carcinoma Mixed Histology (please specify) 
# 8                               12                               37 
# Mucinous Carcinoma                   Other  specify 
# 18                               54 
table(tcga_annot_dt[,"cgc_case_pathologic_stage"])# 
# Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III Stage IIIA Stage IIIB 
# 66         71          6          6        281        218          1        129         25 
# Stage IIIC   Stage IV  Stage Tis    Stage X 
# 50         15          1         13 
# cgc_drug_therapy_drug_name
# table(tcga_annot_dt[,"cgc_drug_therapy_drug_name"])# too many

myplotlab <- "pathologic stage"
stg_ords <-c( "Stage I",  "Stage IA","Stage IB", "Stage II", "Stage IIA" ,"Stage IIB" , "Stage III","Stage IIIA", "Stage IIIB",
              "Stage IIIC",   "Stage IV", "Stage Tis", "Stage X")
p <- plot_pheno_catego(tcga_annot_dt,  pt_traj = all_traj,plotvar= "cgc_case_pathologic_stage", plotxlab=paste0("TCGA BRCA ", myplotlab), varords=stg_ords)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "histological diagnosis"
curr_plotvar <- "cgc_case_histological_diagnosis"
tmp_dt <- tcga_annot_dt[!is.na(tcga_annot_dt[,curr_plotvar]),]
tmp_traj <- all_traj[names(all_traj) %in% tmp_dt$cgc_sample_id]
p <- plot_pheno_catego(tmp_dt,  pt_traj = tmp_traj, plotvar= curr_plotvar, plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



# plot_pheno_catego(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_drug_therapy_drug_name", plotlab="drug therapy", varords=NULL)
## too many categories


myplotlab <- "year of birth"
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "gdc_cases.demographic.year_of_birth", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "days to death"
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_case_days_to_death", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "diag. days to birth"
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "gdc_cases.diagnoses.days_to_birth", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "age at diag."
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "gdc_cases.diagnoses.age_at_diagnosis", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct normal cells"
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_slide_percent_normal_cells", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct neutrophil infilt."
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_slide_percent_neutrophil_infiltration", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct monocyte infilt."
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_slide_percent_monocyte_infiltration", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "slide pct lymphocyt infilt."
p <- plot_pheno_continuous(tcga_annot_dt,  pt_traj = all_traj, plotvar= "cgc_slide_percent_lymphocyte_infiltration", plotxlab=paste0("TCGA BRCA ", myplotlab))
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
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
                                 covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
                                 valuedt=int_dt,
                                             valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= which.min(int_dt$pathway_loading),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr(igene= which.max(int_dt$pathway_loading),
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
               valuedt=int_dt,
               valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr_gg2(igene= which.max(int_dt$pathway_loading),
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "ER+", "ER-"),
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

myconds <- ifelse(rownames(pca_brca_lowrepr) %in% ERpos_samples, "ER+", 
                  ifelse(rownames(pca_brca_lowrepr) %in% ERneg_samples, "ER-", NA))
stopifnot(!is.na(myconds))

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_brca),
            condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nER+=",nERpos, "; nER-=",nERneg))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_pseudotimeGrad_ERpos_ERneg.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,2), summ_dt=summary(pca_brca),
                 condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nER+=",nERpos, "; nER-=",nERneg))

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

ERpos_samp <- which(rownames(pp_input_data) %in% ERpos_samples)
stopifnot(length(ERpos_samp) == nERpos)
stopifnot(rownames(pp_input_data[ERpos_samp]) == rownames(pp_input_data[ERpos_samples]))
corr_expr_pt_ERpos <- apply(pp_input_data[ERpos_samp,], 2, cor, brca_pseudotimes[ERpos_samp])
cdf_ERpos <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt_ERpos)
stopifnot(!is.na(cdf_ERpos))

ERneg_samp <- which(rownames(pp_input_data) %in% ERneg_samples)
stopifnot(length(ERneg_samp) == nERneg)
stopifnot(rownames(pp_input_data[ERneg_samp]) == rownames(pp_input_data[ERneg_samples]))
corr_expr_pt_ERneg <- apply(pp_input_data[ERneg_samp,], 2, cor, brca_pseudotimes[ERneg_samp])
to_keep <- which(!is.na(corr_expr_pt_ERneg))  ### there is 2 genes with all 0 values !!! -> cannot compute corr
cdf_ERneg <- data.frame(feature = all_genes[to_keep], 
                       correlation = corr_expr_pt_ERneg[to_keep])
stopifnot(!is.na(cdf_ERneg))


outFile <- file.path(outFolder, paste0("corr_gene_expr_pseudot_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cdf_ERpos$correlation), main="corr. gene expr. and pseudotimes", col=2)
lines(density(cdf_ERneg$correlation), col=1)
legend("topright", legend=c("ER-", "ER+"), pch=16, col=c(1,2), bty="n")
mtext(side=3, text=paste0("(# ER-=", nERneg, "; # ER+=", nERpos, ")"))
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
# https://github.com/kieranrcampbell/phenopath_revisions/blob/master/analysis/brca_reanalysis/clvm_analysis.Rmd
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
stopifnot(colnames(brca_data_raw)==c(ERpos_samples, ERneg_samples))
stopifnot(colnames(brca_data_raw) == tcga_annot_dt$cgc_sample_id)
stopifnot(!duplicated(tcga_annot_dt$cgc_sample_id))

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
samples_groups <- c(rep("ERpos", length(ERpos_samples)), rep("ERneg", length(ERneg_samples)))
my_group_design <- factor(samples_groups, levels = c("ERpos", "ERneg"))
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

## campbell extracts pvals as follows:
# https://github.com/kieranrcampbell/phenopath_revisions/blob/master/analysis/brca_reanalysis/clvm_analysis.Rmd
qvals <- p.adjust(efit$p.value[,2], method = 'BH')
stopifnot(qvals[rownames(DE_topTable)] == DE_topTable$adj.P.Val )
# ouf... this is equivalent

stopifnot(rownames(efit) %in% rownames(DE_topTable))
stopifnot(rownames(DE_topTable) %in% rownames(efit))
plot(qvals[rownames(DE_topTable)],DE_topTable$adj.P.Val )

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

stopifnot(!duplicated(merged_dt$gene))
stopifnot(!duplicated(df_beta$gene))

outFile <- file.path(outFolder, paste0("limma_adjPval_vs_phenopath_beta.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  y=-log10(merged_dt$adj.P.Val),
  x=merged_dt$beta,
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

##### NICER WAY with colors as in 
# https://github.com/kieranrcampbell/phenopath_revisions/blob/master/analysis/brca_reanalysis/clvm_analysis.Rmd


cols <- RColorBrewer::brewer.pal(3, "Set2")
cols2 <- c("#c5e2d9", cols[2])
p <- ggplot(merged_dt, aes(x = beta, y =-log10(adj.P.Val), color = is_sig)) + 
  geom_point() +
  ylab(expression(paste("limma voom -", log[10], "(adj. P-value)"))) + 
  xlab(expression(paste("Phenopath ", beta))) +
  # theme(legend.position = 'top') +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)) +
  #theme(legend.title = element_text(size = 10),
  #      legend.text = element_text(size = 9)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cols2, name = "Interaction") 
p <- ggExtra::ggMarginal(p, margins = "y", type = "histogram", size = 10)

outFile <- file.path(outFolder, paste0("nicer_limma_adjPval_vs_phenopath_beta.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("limma_logFC_vs_phenopath_beta.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  y=merged_dt$logFC,
  x=merged_dt$beta,
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


############## selected expression comparison  ##############


# SNAI1 vs FBP1 expression supp fig18
genesymb1 <- "FBP1"
stopifnot(genesymb1 %in% ens2genes)
gene_id1 <- names(ens2genes)[ens2genes == genesymb1]
stopifnot(length(gene_id1) == 1)
stopifnot(gene_id1 %in% rownames(brca_data_filtered))

genesymb2 <- "SNAI1"
stopifnot(genesymb %in% ens2genes)
gene_id2 <- names(ens2genes)[ens2genes == genesymb2]
stopifnot(length(gene_id2) == 1)
stopifnot(gene_id2 %in% rownames(brca_data_filtered))

plot_dt <- data.frame(t(brca_data_filtered[c(gene_id1,gene_id2),]))
colnames(plot_dt) <- c("gene1", "gene2")
stopifnot(rownames(plot_dt) == c(ERpos_samples, ERneg_samples))
plot_dt$ERstat <- c(rep("ER+", nERpos), rep("ER-", nERneg))

p <- ggplot(plot_dt, aes(x = gene1, y = gene2, color = ERstat)) +
  geom_point(alpha = 0.6) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = "ER status") +
  xlab("FBP1 expression "~log[2]~"(.+1)") + 
  ylab(expression("SNAIL expression "~log[2]~"(.1)")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  theme_bw()+
  theme(legend.position = "top")# +
# theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 

outFile <- file.path(outFolder, paste0(genesymb2, "_vs_", genesymb1, "_expr.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))
stop("--ok\n")

####################################################################################################
### THRASH
####################################################################################################



```

Time for some reactome fun: NOT EVAL'D

```{r reactome-time, eval = FALSE}
pathways <- c("1643713", "1226099",
              "2219528", "2644603",
              "3304351", "4791275")
id_to_name <- as.list(reactomePATHID2NAME)
pathway_names <- id_to_name[pathways]
pathways_to_genes <- as.list(reactomePATHID2EXTID)
gene_list <- pathways_to_genes[pathways]
pathway_names <- sapply(pathway_names, function(pn) {
  gsub("Homo sapiens: ", "", pn)
})
names(gene_list) <- pathway_names
mart <- useMart("ensembl", "hsapiens_gene_ensembl")
to_ensembl <- function(gl) {
  bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = c("entrezgene"),
      values = as.numeric(gl),
      mart = mart)
  return(bm$ensembl_gene_id)
}
gene_list_ensembl <- lapply(gene_list, to_ensembl)
```

And graph the results for various metrics:

```{r graph-for-metrics, eval = FALSE, eval = FALSE}
sce <- readRDS("../../data/BRCA/sce_brca_gene_level.rds")
all_genes <- unique(unlist(gene_list_ensembl))
all_genes <- all_genes[all_genes %in% rowData(sce)$ensembl_gene_id]
mm <- match(all_genes, rowData(sce)$ensembl_gene_id)
is_basal <- sce$PAM50_mRNA == "Basal-like"
Y <- t(exprs(sce)[mm, ])
colnames(Y) <- all_genes
df_gex <- Y %>% 
  as_data_frame() %>% 
  dplyr::mutate(is_basal, z = colData(sce)[['tmap']]) %>% 
  gather(gene, expression, -is_basal, -z)
df_pheno <- bind_rows(
  lapply(names(gene_list_ensembl), function(n) {
    data_frame(pathway = n, gene = gene_list_ensembl[[n]])
  })
)
df <- inner_join(df_gex, df_pheno, by = 'gene')
df2 <- filter(df, !is.na(is_basal)) %>% 
  group_by(gene, pathway, is_basal) %>% 
  summarise(gex = mean(expression))
df3 <- df %>% # filter(df, !is.na(is_basal)) %>% 
  group_by(pathway, z) %>% 
  summarise(gex = median(expression))
ggplot(df3, aes(x = z, y = gex)) + 
  geom_point(alpha = 0.7) + 
  facet_wrap(~ pathway, scales = 'free_y') +
  ylab("Median pathway expression") +
  stat_smooth(color = 'red', se = F)
df4 <- df %>% # filter(df, !is.na(is_basal)) %>% 
  group_by(pathway, z) %>% 
  summarise(gex = var(expression))
# ggplot(df4, aes(x = z, y = gex)) + 
#   geom_point(alpha = 0.7) + 
#   facet_wrap(~ pathway, scales = 'free_y') +
#   ylab("Variance in pathway expression") +
#   stat_smooth(color = 'red', se = F)
```





# Crossover bit for chris

```{r df-beta}
df_beta <- mutate(df_beta,
                  crossover = - alpha / beta)
filter(df_beta, is_sig) %>% 
  ggplot(aes(x = crossover)) + 
  geom_histogram(fill = "#74a9cf", color = "grey90", bins = 30) +
  xlab("Crossover point") + ylab("Number of genes") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10))
cross_plot <- last_plot()
ggsave("../../figs/supplementary_crossover.png", width = 6, height = 5)
```



Crossover GO analysis:

```{r crossover-go}
library(goseq)
genome <- "hg19"
id <- "ensGene"
all_genes <- rowData(sce)$ensembl_gene_id
crossover_genes <- filter(df_beta, 
                          is_sig & crossover > 0.4) %>% 
  extract2("ensembl_gene_id") 
crossg <- 1 * (all_genes %in% crossover_genes)
names(crossg) <- all_genes
pwfcross <- nullp(crossg, genome, id)
gocross <- goseq(pwfcross, genome, id, test.cats = "GO:BP")
pgo <- parse_go(gocross, "crossover", length(all_genes))
```

Crossover plots:

```{r crossover-gene-plots}
tmap <- pcavi$m_z
cross_df <- dplyr::select(df_beta, gene, crossover, hgnc_symbol)
top_genes <- dplyr::filter(df_beta, is_sig) %>% 
  arrange(desc(abs(beta))) %>% 
  extract2("gene") %>% head(n=12)
df_gex <- t(exprs(sce))[, top_genes, drop=FALSE] %>% 
  as_data_frame() %>% 
  dplyr::mutate(phenotime = tmap, x = as.character(x)) %>% 
  gather(gene, expression, -phenotime, -x)
df_gex$x[is.na(df_gex$x)] <- "NA"
df_gex$x <- plyr::mapvalues(df_gex$x, from = sort(unique(df_gex$x)), to = c("ER negative", "ER positive"))
df_gex$x <- factor(df_gex$x, levels = c("ER negative", "ER positive"))
df_gex <- inner_join(df_gex, cross_df, by = "gene")
ggplot(df_gex, aes(x = phenotime, y = expression, color = x)) + 
  geom_point(alpha = 0.1) +
  facet_wrap(~ hgnc_symbol, scales = "free_y", ncol = 3) + 
  theme(legend.position = "top", strip.text.x = element_text(size = 9)) +
  stat_smooth(se = FALSE, method = "lm", size = 2) + 
  scale_color_brewer(name = "", palette = 'Set1') +
  geom_vline(aes(xintercept = crossover), linetype = 2) +
  xlab("Pathway score") +
  ylab(expression(Expression ~ log[2] ~ "(TPM + 1)")) 
# saveRDS(last_plot(), "../../figs/brca/crossover_thesis.rds")
# 
# ggsave("../../figs/supplementary_crossover_2.png", width = 10, height = 9)
```



# Plots for paper

```{r more-vis}
textinfo <- frame_data(
  ~x, ~y, ~label,
  0.6, 0.15, "Gene upregulated\nER+ increases upregulation",
  0.6, -0.15, "Gene upregulated\nER+ decreases upregulation",
  -0.7, 0.15, "Gene downregulated\nER+ decreases downregulation",
  -0.7, -0.15, "Gene downregulated\nER+ increases downregulation"
)
cols <- RColorBrewer::brewer.pal(3, "Set2")
cols2 <- c("#c5e2d9", cols[2])
outline_cols = c("#c5e2d9", 'black')
ggplot(df_beta, aes(x = c, y = beta)) + 
  geom_point(shape = 21, aes(fill = is_sig_graph, color = is_sig_graph), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_fill_manual(values = cols2, name = "Interaction") +
  scale_color_manual(values = outline_cols, name = "Interaction") +
  geom_text_repel(data = dplyr::filter(df_beta, is_sig, abs(beta) > 0.7),
                 aes(label = hgnc_symbol), color = 'black',
                 size = 3) +
  ylab("Covariate-pseudotime interaction") +
  xlab("Gene regulation over pseudotime") +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) +
  geom_text(data = textinfo, aes(x = x, y = y, label = label), 
            color = 'black', size = 3, fontface = "bold")
lplot <- last_plot()
lplot




New heatmap plots:

```{r, cache = FALSE}
sce_gene <- readRDS("../../data/BRCA/sce_brca_gene_level.rds")
sce_gene <- updateSCESet(sce_gene)
sce_gene <- sce_gene[, colnames(sce)]
angiogenesis_genes_all <- read_csv("../../data/BRCA/ag_hallmark_geneset.txt", skip = 2, col_names = FALSE)[[1]]
df_beta <- mutate(df_beta, is_angiogenesis = hgnc_symbol %in% angiogenesis_genes_all)
growth_factor_grep <- c("^FGF", "^VEGF", "^EGF", "^TGF")
growth_factors <- unlist(lapply(growth_factor_grep, grep, rowData(sce)$hgnc_symbol, value = TRUE))
growth_factors <- growth_factors[growth_factors != "EFGLAM"]
df_beta <- mutate(df_beta, is_gf = hgnc_symbol %in% growth_factors)
ag_to_plot <- filter(df_beta, is_gf) %>%
  arrange(desc(c)) %>%
  head(n = 10) %>%
  .$hgnc_symbol
ag_to_plot <- c(ag_to_plot, "VEGFC")
ag_inds <- match(ag_to_plot, rowData(sce_gene)$hgnc_symbol)
marker_mat <- t(exprs(sce_gene))[, ag_inds, drop=FALSE]
colnames(marker_mat) <- ag_to_plot
marker_df <-  marker_mat %>% 
  as_data_frame() %>% 
  dplyr::mutate(pseudotime = tmap)
mx <- as.matrix(dplyr::select(marker_df, -pseudotime))
d <- dist(t(mx))
hc <- hclust(d)
gene_order <- rev(colnames(mx)[hc$order])
winsorize <- function(x, lower = -2.5, upper = 2.5) {
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}
marker_df_2 <- mutate(marker_df, pseudotime_order = rank(pseudotime)) %>% 
  gather(gene, expression, -pseudotime_order, -pseudotime)
marker_df_2$gene <- factor(marker_df_2$gene, levels = gene_order)
marker_df_2 <- group_by(marker_df_2, gene) %>% 
  mutate(norm_expression = winsorize((expression - mean(expression)) / sd(expression)))
ggplot(marker_df_2, aes(x = pseudotime_order, y = gene, fill = norm_expression)) +
         geom_raster(interpolate = FALSE) +
  scale_fill_viridis(name = "Expression\nz-score") +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7, margin = margin(r = -0.5, l = 0, unit = "cm")),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10, margin = margin(t = 0, unit = "cm")),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.direction = "vertical",
        legend.margin = margin(l = -.5, unit = "cm")) +
  labs(y = "Gene", x = "Pseudotime order")
heatmap_plot <- last_plot()
```






# Demanding supervisor plots

## Vascular growth factors

```{r vasc-growth-fact}
x_str <- plyr::mapvalues(x, from = sort(unique(x)),
                         to = c("ER-", "ER+"))
tmap <- sce$tmap
vgf <- c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "FGF2", "CXCL8") 
mm <- match(vgf, rowData(sce_gene)$hgnc_symbol)
vgf_exprs <- t(exprs(sce_gene)[mm, ]) %>% as_data_frame()
names(vgf_exprs) <- vgf
vgf_exprs <- mutate(vgf_exprs, z = tmap, ER_status = x_str)
vgf_exprs_tidy <- gather(vgf_exprs, gene, expression, -z, -ER_status)
ggplot(vgf_exprs_tidy, aes(x = z, y = expression, color = ER_status)) +
  geom_point(alpha = 0.6) + facet_wrap(~ gene, scales = "free_y") +
  xlab("Pathway score") + ylab("Expression") +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = "ER status") 
ggsave("../../figs/vgf_genes.png", width = 8, height = 6)
```

"Can we plot these brca genes?"

```{r brca-genes}
brca_genes <- c("CEP55", "ESR1", "FOXA1", "FOXC1", "KRT17", "MAPT", "MELK", "MMP11", "NAT1", "SFRP1", "UBE2C", "UBE2T")
mm <- match(brca_genes, rowData(sce_gene)$hgnc_symbol)
brca_exprs <- t(exprs(sce_gene)[mm, ]) %>% as_data_frame()
names(brca_exprs) <- brca_genes
brca_exprs <- mutate(brca_exprs, z = tmap, ER_status = x_str)
brca_exprs_tidy <- gather(brca_exprs, gene, expression, -z, -ER_status)
ggplot(brca_exprs_tidy, aes(x = z, y = expression, color = ER_status)) +
  geom_point(alpha = 0.6) + facet_wrap(~ gene, scales = "free_y") +
  xlab("Pathway score") + ylab("Expression") +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = "ER status") 
ggsave("../../figs/brca_genes.png", width = 9, height = 6)
```

"Can we look at FBP1 vs SNAIL expression?"

```{r fbp1-vs-snail}
genes <- c("FBP1", "SNAI1")
mm <- match(genes, rowData(sce_gene)$hgnc_symbol)
fs_exprs <- t(exprs(sce_gene)[mm, ]) %>% as_data_frame()
names(fs_exprs) <- genes
fs_exprs <- mutate(fs_exprs, ER_status = x_str)
ggplot(fs_exprs, aes(x = FBP1, y = SNAI1, color = ER_status)) +
  geom_point(alpha = 0.6) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = "ER status") +
  xlab("FBP1 expression "~log[2]~"(TPM+1)") + 
  ylab(expression("SNAIL expression "~log[2]~"(TPM+1)"))
ggsave("../../figs/fbp1-vs-snail.png", width = 4, height = 4)
```

