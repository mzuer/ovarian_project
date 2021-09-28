
# Rscript phenopath_recount2_TCGA_BRCA_lumAlumB.R

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
require(plyr)
require(edgeR)
require(matrixStats)
require(goseq)
require(reshape2)
require(mclust)
require(ggExtra)
require(ReactomePA)
require(reactome.db)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(tidyr)
library(gridExtra)
library(ggplot2)
require(survival)
require(survminer)

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
runNorm <- F


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

inFolder <- file.path("..","tcga_data","DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2")

# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

outFolder <- "PHENOPATH_RECOUNT2_TCGA_BRCA_LUMALUMB"
dir.create(outFolder, recursive=TRUE)


# load input data
brca_data_raw <- get(load(file.path(inFolder, paste0("all_counts_onlyPF_", purityFilter, ".Rdata"))))
dim(brca_data_raw)
# [1] 25526   882

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

cat(paste0("... grep brca_data_raw"))              
grep("ENSG00000124216", rownames(brca_data_raw))


httr::set_config(httr::config(ssl_verifypeer = FALSE))  ### added to access ensembl biomart connection
tokeep <- ! grepl("_PAR_Y$", rownames(brca_data_raw))
cat(paste0("... remove weird genes: ", sum(!tokeep), "/", length(tokeep), "\n"))
brca_data_raw <- brca_data_raw[tokeep,]
dim(brca_data_raw)
# [1] 25503   882
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
# [1] 19159   533

cat(paste0("... grep brca_data_forNorm"))              
grep("ENSG00000124216", rownames(brca_data_forNorm))

geneIDs_beforeNorm <- rownames(brca_data_forNorm)

outFile <- file.path(outFolder, "brca_data_forNorm.Rdata")
save(brca_data_forNorm, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


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
cat(paste0("... grep brca_dat_gcNorm"))              
grep("ENSG00000124216", rownames(brca_dat_gcNorm))

# the function reorder the gene...
# reput in same order (just to be sure)
geneIDs_afterNorm <- rownames(brca_dat_gcNorm)
stopifnot(geneIDs_afterNorm %in% geneIDs_beforeNorm)
kept_genes <- geneIDs_beforeNorm[geneIDs_beforeNorm %in% geneIDs_afterNorm]
stopifnot(setequal(kept_genes, rownames(brca_dat_gcNorm)))
stopifnot(length(kept_genes) == nrow(brca_dat_gcNorm))
brca_dat_gcNorm <- brca_dat_gcNorm[kept_genes,]
cat(paste0("... grep brca_dat_gcNorm reordererd"))              
grep("ENSG00000124216", rownames(brca_dat_gcNorm))

brca_data_gcNorm_log <- log2(brca_dat_gcNorm + 1)

nTCGA <- sum(grepl("^TCGA", colnames(brca_data_gcNorm_log)))
stopifnot(nTCGA == ncol(brca_data_gcNorm_log))

stopifnot(tcga_annot_dt$PAM50 %in% c("lumA", "lumB"))

lumA_samples <- colnames(brca_data_gcNorm_log)[colnames(brca_data_gcNorm_log) %in%
                                                  tcga_annot_dt$cgc_sample_id[tcga_annot_dt$PAM50 == "lumA"]]
nLumA <- length(lumA_samples)

lumB_samples <- colnames(brca_data_gcNorm_log)[colnames(brca_data_gcNorm_log) %in%
                                                  tcga_annot_dt$cgc_sample_id[tcga_annot_dt$PAM50 == "lumB"]]
nLumB <- length(lumB_samples)

stopifnot(colnames(brca_data_gcNorm_log) == c(lumA_samples, lumB_samples))


################################### 
################################### PCA ON not norm DATA -> for comparison
################################### 
brca_data_notNorm_log <- log2(brca_data_forNorm + 1)

brca_data_notNorm_log_no0 <- brca_data_notNorm_log[rowSums(brca_data_notNorm_log) > 0,]

pca_brca <- prcomp(t(brca_data_notNorm_log_no0), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nTCGA)

stopifnot(rownames(pca_brca_lowrepr) == tcga_annot_dt$cgc_sample_id)

mycolvect <- 1+as.numeric(tcga_annot_dt$PAM50 == "lumA")

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_13_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA notNorm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
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

mycolvect <- 1+as.numeric(tcga_annot_dt$PAM50 == "lumA")

outFile <- file.path(outFolder, paste0("in_norm_data_pca_12_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_norm_pca_23_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_norm_pca_13_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA norm notFilt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))


#### ADD SOME FILTERING HERE

# For input to PhenoPath we used the 4,579 genes whose variance in log(TPM+1)
# expression was greater than 1 and whose median absolute deviation was greater than 0


cat(paste0("... grep brca_data_gcNorm_log"))              
grep("ENSG00000124216", rownames(brca_data_gcNorm_log))

var_exprs <- rowVars(brca_data_gcNorm_log)
mad_exprs <- rowMads(brca_data_gcNorm_log)
to_keep <-  var_exprs > var_thresh & mad_exprs > mad_thresh

stopifnot(colnames(brca_data_gcNorm_log) == c(lumA_samples, lumB_samples))

var_exprs <- setNames(var_exprs,colnames(brca_data_gcNorm_log)) 
mad_exprs <- setNames(mad_exprs,colnames(brca_data_gcNorm_log)) 

# dev.off()
outFile <- file.path(outFolder, paste0("all_vars_lumA_lumB_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(var_exprs), main="all vars (log2(.+1))")
mtext(side=3, text="(all)")
abline(v=var_thresh, col="red")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("all_mads_lumA_lumB_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(mad_exprs), main="all mads (log2(.+1))")
mtext(side=3, text="(all)")
abline(v=var_thresh, col="red")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

nGenes <- sum(to_keep)

brca_data_notFiltered <- brca_data_gcNorm_log
brca_data_gcNorm_log <- brca_data_gcNorm_log[to_keep,]
cat(paste0("... keep highly variable genes ", nrow(brca_data_gcNorm_log), "/", length(to_keep), "\n"))

cat(paste0("... grep brca_data_gcNorm_log"))              
grep("ENSG00000124216", rownames(brca_data_gcNorm_log))


brca_data_filtered <- brca_data_gcNorm_log
outFile <- file.path(outFolder, paste0("brca_data_filtered.Rdata"))
save(brca_data_filtered, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

# stop("--ok\n")
# SNAI did not pass the filter


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

mycolvect <- 1+as.numeric(tcga_annot_dt$PAM50 == "lumA")

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_12_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_23_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_13_lumA_lumB.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
legend("topleft", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))


# ################################## remove outlier with mclust [not needed with the norm  + filt ???]



mc <- Mclust(pca_brca_lowrepr[,c(1,3)], G = 2)

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_13_lumA_lumB_withClust.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="TCGA BRCA Norm+Filt (log2(.+1))",
        col=mc$classification)
mtext(text=paste0("nLumB=",nLumB, "; nLumA=",nLumA), side=3)
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

# so here give a value of 1 in lumA and a value of -1 in lumB
# This means the overall pathway loading λ is the average change for normal and tumor
# if β > 0 =>  the gene is more upregulated over pseudotime under tumor  
# if β  < 0 =>   the gene is more upregulated under normal 

brca_data_filteredT <- t(brca_data_filtered)
stopifnot(nGenes == ncol(brca_data_filteredT))
stopifnot(nLumA+nLumB == nrow(brca_data_filteredT))

# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(brca_data_filteredT) == nLumB+nLumA)
stopifnot(ncol(brca_data_filteredT) == nGenes)
mycovar <- 2 * ( rownames(brca_data_filteredT) %in% lumA_samples) - 1
stopifnot(sum(mycovar== -1) == nLumB)
stopifnot(sum(mycovar== 1) == nLumA)
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
# brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_brca/brca_phenopath_fit.Rdata"))
brca_pseudotimes <-  trajectory(brca_phenopath_fit)


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
stopifnot(colnames(brca_data_raw)==c(lumA_samples, lumB_samples))
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
# 21636/25526

## Build DGEList object
dge <- DGEList(counts = brca_data_raw[to_keep, ])
## Calculate normalization factors
dge <- calcNormFactors(dge)

x <- c(rep(1, length(lumA_samples)), rep(-1, length(lumB_samples)))

### was done like:
# design <- model.matrix(~ x, colData(sce))
# were x=sce$x
# from Tuto -> sce$x == colData(sce) except that colData(sce) returns in data frame format
# usable by model.matrix -> not totally sure !!!!
design <- model.matrix(~ x, as.data.frame(x))
stopifnot(design[,2] == x)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- decideTests(fit)

outFile <- file.path(outFolder, paste0("limma_v2_venn_diag.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
vennDiagram(results)
mtext(side=3, text="limma v2")
dev.off()
cat(paste0("... written: ", outFile, "\n"))



limma_genes <- gsub("(.+)\\..+", "\\1", rownames(fit$coefficients))
stopifnot(!duplicated(limma_genes))
stopifnot(nrow(as.data.frame(fit)) == length(limma_genes))
qvals <- p.adjust(fit$p.value[,2], method = 'BH')

df_limma <- data_frame(coef = fit$coefficients[,2], 
                       pval = fit$p.value[,2],
                       qval = qvals)
stopifnot(nrow(df_limma) == length(limma_genes))
df_limma$feature <- limma_genes


tmp_int_dt <- as.data.frame(interactions(brca_phenopath_fit))
stopifnot(tmp_int_dt$feature == brca_phenopath_fit$feature_names)
tmp_int_dt$beta <- brca_phenopath_fit$m_beta[1,]
tmp_int_dt$alpha <- brca_phenopath_fit$m_alpha[1,]
tmp_int_dt$mu <- brca_phenopath_fit$m_mu
stopifnot(tmp_int_dt$beta == tmp_int_dt$interaction_effect_size)

limma_merged_dt <- merge(tmp_int_dt, df_limma, by="feature", all.x=T,all.y=F)
stopifnot(nrow(limma_merged_dt) == nrow(tmp_int_dt))
stopifnot(!is.na(limma_merged_dt))
stopifnot(setequal(limma_merged_dt$feature, tmp_int_dt$feature))
limma_merged_dt$log10qval <- -log10(limma_merged_dt$qval)
limma_merged_dt$limma_sig <- limma_merged_dt$qval < 0.05

p <- ggplot(limma_merged_dt, aes(x = coef, y = alpha, color = mu)) +
  geom_point(alpha = 0.5) +
  xlab("Limma voom coefficient") +
  ylab(expression(paste("Phenotime ", alpha))) +
  scale_color_viridis(name = expression(mu))

outFile <- file.path(outFolder, paste0("limma_v2_coef_vs_phenopath_alpha.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# There is some "twisting" there which may be due to using a mixed effects model rather than straight forward DE.

limma_merged_dt$is_sig_graph <- ifelse(limma_merged_dt$significant_interaction,
                                       "Significant", "Non-significant")


cols <- RColorBrewer::brewer.pal(3, "Set2")
cols2 <- c("#c5e2d9", cols[2])
p <- ggplot(limma_merged_dt, x =-log10(qval), aes(x = beta, y =-log10(qval), color = is_sig_graph)) + 
  geom_point() +
  ylab(expression(paste("Limma voom -", log[10], "(q-value)"))) + 
  xlab(expression(paste("PhenoPath ", beta))) +
  # theme(legend.position = 'top') +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)) +
  #theme(legend.title = element_text(size = 10),
  #      legend.text = element_text(size = 9)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cols2, name = "Interaction") 
p <- ggExtra::ggMarginal(p, margins = "y", type = "histogram", size = 10)


outFile <- file.path(outFolder, paste0("nicer_limma_v2_adjPval_vs_phenopath_beta.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, "limma_merged_dt.Rdata")
save(limma_merged_dt, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "fit.Rdata")
save(fit, file=outFile)
cat(paste0("... written: ", outFile, "\n"))




# stop("--ok \n") 





















 
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
add_corr(pc1, brca_pseudotimes, cond_plus1="lumA", cond_minus1="lumB")
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
add_corr(pc2, brca_pseudotimes, cond_plus1="lumA", cond_minus1="lumB")
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
add_corr(pc3, brca_pseudotimes, cond_plus1="lumA", cond_minus1="lumB")
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
    ggtitle(paste0("Genes with top ", betaU), subtitle=paste0("(", nTop, " lowest and highest)"))+
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
i_max <- which.max(df_beta$beta)
p <- plot_iGeneExpr(igene= i_max,
              exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
               valuedt=df_beta,
               valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_", df_beta$geneSymb[i_max], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_max,
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_", df_beta$geneSymb[i_max], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

i_min <-which.min(df_beta$beta)

p <- plot_iGeneExpr(igene= i_min,
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
               valuedt=df_beta,
                           valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_", df_beta$geneSymb[i_min], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_min,
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_", df_beta$geneSymb[i_min], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


# plot selected genes from Campbell and Yau 2018
selected_symbs <- c("ESR1", "FOXC1", "FBP1","FGD5", "METTL6", "CPT1A",
                    "DTX3", "MRPS23", "EIF2S2", "EIF6", "SLC2A10",
                    "PI3KCA", "ERBB2", "ERBB3", "ESR1",
                    "ERS1", "GATA3", "FOXA1", "XBP1") 
# 1st row from the Campbell
# 2nd row Gatza et al. "uniquely amplified in patients with highly proliferative luminal breast"
# 3d row "marker genes" from Harbeck et al. 2019 - lumB
# 4th row "marker genes" from Harbeck et al. 2019 - lumA [NB: ERS1 likely a typo, should be ESR1]

gs="FOXC1"
stopifnot(df_beta$gene == colnames(brca_data_filteredT))
stopifnot(df_beta$gene %in% names(ens2genes))
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- ens2genes[df_beta$gene]

## many of the 20 top beta value genes exhibit convergence
stopifnot(brca_phenopath_fit$m_beta[1,] == df_beta$beta)
stopifnot(brca_phenopath_fit$feature_names == df_beta$gene)

df_beta$alpha <- brca_phenopath_fit$m_alpha[1,]
df_beta$crossover <- -df_beta$alpha/df_beta$beta
stopifnot(!is.na(df_beta))

                                                for(gs in selected_symbs) {
                                                  if(!gs %in% df_beta$geneSymb) {
                                                    cat(paste0("warning: ", gs, " not available\n"))
                                                    next
                                                  }
                                                  i_gs <- which(df_beta$geneSymb == gs)
                                                  ens_gs <- df_beta$gene[i_gs]
                                                  stopifnot(ens_gs == names(ens2genes)[ens2genes==gs])
                                                  cp_gs <- df_beta$crossover[i_gs]
                                                  cat(paste0(gs, "- length i_gs = ",length(i_gs) ,"\n" ))
                                                  # stopifnot(length(i_gs) == 1)
                                                  tmp_dt <- df_beta
                                                  tmp_dt <- tmp_dt[order(abs(tmp_dt$beta), decreasing = TRUE),]

                                                  gs_beta_rank <- which(tmp_dt$gene == ens_gs)
                                                  stopifnot(length(gs_beta_rank) == 1)
                                                  gs_signif <- as.character(df_beta$is_sig[i_gs])
                                                  gs_beta <- round(df_beta$beta[i_gs], 3)

                                                  p <- plot_iGeneExpr_gg2(igene= i_gs,
                                                                          exprdt=brca_data_filteredT,
                                                                          pseudot=brca_pseudotimes,
                                                                          covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                                                                          valuedt=df_beta,
                                                                          valuecol="beta", symbcol="geneSymb",
                                                                          subtit=paste0(betaU, " rank=", gs_beta_rank, "; signif=",gs_signif ))

                                                  # add crossoverline
                                                  p <- p + geom_vline(xintercept = cp_gs, linetype=2)

                                                  outFile <- file.path(outFolder, paste0(gs, "_geneExpr_along_pseudotime_withCP_gg2.", plotType))
                                                  ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
                                                  cat(paste0("... written: ", outFile, "\n"))

                                                }


#### do the same for the top 9 (3x3 grid) pos beta 
top_gs <- df_beta$gene[order(df_beta$beta, decreasing = TRUE)][1:9]

all_p <- list()

for(i in 1:length(top_gs)) {
  ens_gs <- top_gs[i]
  if(!ens_gs %in% names(ens2genes)) {
    cat(paste0("warning: ", ens_gs, " not available\n"))
    next
  }
  i_gs <- which(df_beta$gene == ens_gs)
  stopifnot(ens_gs == df_beta$gene[i_gs])
  cat(paste0(i), "\n")
  # cat(paste0(gs), "\n")
  cat(paste0(ens_gs), "\n")
  # cat(paste0(which(ens2genes==gs)), "\n")
  # cat(paste0(names(ens2genes)[ens2genes==gs]), "\n")
  cp_gs <- df_beta$crossover[i_gs]
  
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
                          covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                          valuedt=df_beta,
                          valuecol="beta", symbcol="geneSymb", 
                          subtit=paste0(betaU, " rank=", gs_beta_rank, "; signif=",gs_signif ))
  
  # add crossoverline
  p <- p + geom_vline(xintercept = cp_gs, linetype=2)
  all_p[[i]] <- p
}
ag_allp <- do.call(gridExtra:::grid.arrange,c(all_p, ncol=3))


outFile <- file.path(outFolder, paste0("all_topPosBeta_geneExpr_along_pseudotime_withCP_gg2.", plotType))
ggsave(plot = ag_allp, filename = outFile, height=myHeightGG*1.9, width = myWidthGG*3.2)
cat(paste0("... written: ", outFile, "\n"))


# stop("--ok\n")

#### do the same for the top 12 (3x4 grid) neg beta 

#### do the same for the top 9 (3x3 grid) pos beta 
top_gs <- df_beta$gene[order(df_beta$beta, decreasing = FALSE)][1:9]

all_p <- list()

for(i in 1:length(top_gs)) {
  ens_gs <- top_gs[i]
  if(!ens_gs %in% names(ens2genes)) {
    cat(paste0("warning: ", ens_gs, " not available\n"))
    next
  }
  i_gs <- which(df_beta$gene == ens_gs)
  stopifnot(ens_gs == df_beta$gene[i_gs])
  cat(paste0(i), "\n")
  # cat(paste0(gs), "\n")
  cat(paste0(ens_gs), "\n")
  # cat(paste0(which(ens2genes==gs)), "\n")
  # cat(paste0(names(ens2genes)[ens2genes==gs]), "\n")
  cp_gs <- df_beta$crossover[i_gs]
  
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
                          covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                          valuedt=df_beta,
                          valuecol="beta", symbcol="geneSymb", 
                          subtit=paste0(betaU, " rank=", gs_beta_rank, "; signif=",gs_signif ))
  
  # add crossoverline
  p <- p + geom_vline(xintercept = cp_gs, linetype=2)
  all_p[[i]] <- p
}
ag_allp <- do.call(gridExtra:::grid.arrange,c(all_p, ncol=3))


outFile <- file.path(outFolder, paste0("all_topNegBeta_geneExpr_along_pseudotime_withCP_gg2.", plotType))
ggsave(plot = ag_allp, filename = outFile, height=myHeightGG*1.9, width = myWidthGG*3.2)
cat(paste0("... written: ", outFile, "\n"))




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
p <- plot_pheno_catego(tcga_annot_dt,  pt_traj = all_traj,
                       plotvar= "cgc_case_pathologic_stage", plotxlab=paste0("TCGA BRCA ", myplotlab),
                       varords=stg_ords)
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
# - significant_interaction => Boolean for whether the interaction effect is significantly different from 0
# - chi => The precision of the ARD prior on β [Automatic Relevance Determination]
# - pathway_loading => The pathway loading λ

int_dt <- interactions(brca_phenopath_fit)
stopifnot(as.character(int_dt$feature) %in% names(ens2genes))
int_dt$featureSymb <- ens2genes[as.character(int_dt$feature) ]
stopifnot(!is.na(int_dt$featureSymb))
# plot the posterior ARD variances (1/χ) against the posterior interaction effect sizes (β)
# colouring them by which are found to be significant and annotating the top few genes:

chi_cutoff <- sort(int_dt$chi)[11]


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
ggsave(plot = p, filename = outFile, height=myHeightGG, width=myWidthGG*1.4)
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
ggsave(plot = p, filename = outFile, height=myHeightGG, width=myWidthGG*1.4)
cat(paste0("... written: ", outFile, "\n"))



# stop("--ok\n")



int_dt$is_sig_graph <- 
  plyr::mapvalues(int_dt$significant_interaction, from = c(FALSE, TRUE),
                  to = c("Non-significant", "Significant"))

textinfo <- frame_data(
  ~x, ~y, ~label,
  0.6, 0.15, "Gene upregulated\nLumA increases upregulation",
  0.6, -0.15, "Gene upregulated\nLumA decreases upregulation",
  -0.7, 0.15, "Gene downregulated\nLumA decreases downregulation",
  -0.7, -0.15, "Gene downregulated\nLumA increases downregulation"
)
cols <- RColorBrewer::brewer.pal(3, "Set2")
cols2 <- c("#c5e2d9", cols[2])
outline_cols = c("#c5e2d9", 'black')
p <- ggplot(int_dt, aes(x = pathway_loading, y = interaction_effect_size)) + 
  geom_point(shape = 21, aes(fill = is_sig_graph, color = is_sig_graph), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_fill_manual(values = cols2, name = "Interaction") +
  scale_color_manual(values = outline_cols, name = "Interaction") +
  geom_text_repel(data = dplyr::filter(int_dt, significant_interaction, abs(interaction_effect_size) > 0.7),
                  aes(label = featureSymb), color = 'black',
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


outFile <- file.path(outFolder, paste0("posteriorEffectSizeBeta_vs_pathwayloadingLambda_nicer.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

# stop("ok\n")

############## look at the gene with top and bottom pathway loading ##############

int_dt <- as.data.frame(int_dt)
tmp1 <- as.data.frame(int_dt)
tmp2 <- df_beta
tmp1 <- tmp1[order(tmp1$interaction_effect_size),]
tmp2 <- tmp2[order(tmp2$beta),]
stopifnot(tmp1$featureSymb == tmp2$geneSymb)

i_min <- which.min(int_dt$pathway_loading)
p <- plot_iGeneExpr(igene= i_min,
                                 exprdt=brca_data_filteredT,
                                 pseudot=brca_pseudotimes,
                                 covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                                 valuedt=int_dt,
                                             valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_", int_dt$featureSymb[i_min], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_min,
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_", int_dt$featureSymb[i_min], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

i_max <- which.max(int_dt$pathway_loading)
p <- plot_iGeneExpr(igene= i_max,
               exprdt=brca_data_filteredT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
               valuedt=int_dt,
               valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_", int_dt$featureSymb[i_max], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr_gg2(igene= i_max,
                    exprdt=brca_data_filteredT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, "lumA", "lumB"),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_", int_dt$featureSymb[i_max], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
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
  ylab(expression(lambda)) +
  scale_color_brewer(palette = "Set2", name = "Significant")

outFile <- file.path(outFolder, paste0("genes_with_top_and_bottom_n", nTop, "_lambda.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

############## PC dots with color-coded by pseudotime gradient ##############  <<<<<<<<<<<<<<< FIG 6a

myconds <- ifelse(rownames(pca_brca_lowrepr) %in% lumA_samples, "lumA", 
                  ifelse(rownames(pca_brca_lowrepr) %in% lumB_samples, "lumB", NA))
stopifnot(!is.na(myconds))

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_brca),
            condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nLumA=",nLumA, "; nLumB=",nLumB))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_pseudotimeGrad_lumA_lumB.", plotType))
ggsave(p, filename = outFile, height=myHeightGG*0.9, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,2), summ_dt=summary(pca_brca),
                 condvect = myconds,
                 colvect=brca_pseudotimes,
                 mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
                 mysubtit = paste0("nLumA=",nLumA, "; nLumB=",nLumB))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_pseudotimeGrad_lumA_lumB.", plotType))
ggsave(p, filename = outFile, height=myHeightGG*0.9, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,3), summ_dt=summary(pca_brca),
                 condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("nLumA=",nLumA, "; nLumB=",nLumB))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_13_pseudotimeGrad_lumA_lumB.", plotType))
ggsave(p, filename = outFile, height=myHeightGG*0.9, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

############## gene set enrichment on top and bottom beta value genes ##############  <<<<<<<<<<<<<<< FIG 6c
# A GO enrichment analysis indicated that the genes driving the inferred pseudotemporal trajectory 
# were indeed enriched for vascular growth pathways (Fig. 6c).
# -> how is "driving" defined ???  is it with all genes ???

# from the code https://github.com/kieranrcampbell/phenopath_analyses/blob/058fbd77dea3be6de03bf57a319d272a0eb15c64/analysis/BRCA/clvm_analysis.Rmd
# https://github.com/kieranrcampbell/phenopath_revisions/blob/2c52ed25cedb0abb49465e6bf71b1b25640fd068/analysis/brca_reanalysis/clvm_analysis.Rmd
# https://github.com/kieranrcampbell/phenopath_analyses/blob/058fbd77dea3be6de03bf57a319d272a0eb15c64/analysis/BRCA/clvm_analysis_er_pos.Rmd
# upreg_genes <- filter(cdf, correlation > 0.5) %>% extract2("feature") 
# downreg_genes <- filter(cdf, correlation < -0.5) %>% extract2("feature")
# gos <- bind_rows(
#   parse_go(goup, "Up-regulated", length(upreg_genes)),
#   parse_go(godown, "Down-regulated", length(downreg_genes))
# )

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

lumA_samp <- which(rownames(pp_input_data) %in% lumA_samples)
stopifnot(length(lumA_samp) == nLumA)
stopifnot(rownames(pp_input_data[lumA_samp]) == rownames(pp_input_data[lumA_samples]))
corr_expr_pt_lumA <- apply(pp_input_data[lumA_samp,], 2, cor, brca_pseudotimes[lumA_samp])
cdf_lumA <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt_lumA)
stopifnot(!is.na(cdf_lumA))

lumB_samp <- which(rownames(pp_input_data) %in% lumB_samples)
stopifnot(length(lumB_samp) == nLumB)
stopifnot(rownames(pp_input_data[lumB_samp]) == rownames(pp_input_data[lumB_samples]))
corr_expr_pt_lumB <- apply(pp_input_data[lumB_samp,], 2, cor, brca_pseudotimes[lumB_samp])
to_keep <- which(!is.na(corr_expr_pt_lumB))  ### there is 2 genes with all 0 values !!! -> cannot compute corr
cdf_lumB <- data.frame(feature = all_genes[to_keep], 
                       correlation = corr_expr_pt_lumB[to_keep])
stopifnot(!is.na(cdf_lumB))


outFile <- file.path(outFolder, paste0("corr_gene_expr_pseudot_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cdf_lumA$correlation), main="corr. gene expr. and pseudotimes", col=2)
lines(density(cdf_lumB$correlation), col=1)
legend("topright", legend=c("lumB", "lumA"), pch=16, col=c(1,2), bty="n")
mtext(side=3, text=paste0("(# lumB=", nLumB, "; # lumA=", nLumA, ")"))
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
ggsave(p, filename = outFile, height=myHeightGG*1.4, width=myWidthGG*1.5)
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
  geom_segment(aes(y = min(log10qval - 1), yend = log10qval, x = term, xend = term),
               color = 'grey30', linetype = 3)+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(),
      axis.text = element_text(size = 10),
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
stopifnot(colnames(brca_data_raw)==c(lumA_samples, lumB_samples))
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
# 21636/25526

## Build DGEList object
dge <- DGEList(counts = brca_data_raw[to_keep, ])
## Calculate normalization factors
dge <- calcNormFactors(dge)
# plotMDS(dge)
samples_groups <- c(rep("lumA", length(lumA_samples)), rep("lumB", length(lumB_samples)))
my_group_design <- factor(samples_groups, levels = c("lumA", "lumB"))
my_design <- model.matrix( ~ my_group_design)
## Run voom
v <- voom(dge, my_design, plot = TRUE)
## Run remaining parts of the DE analysis
fit <- lmFit(v, my_design)
efit <- eBayes(fit)


results <- decideTests(efit)
outFile <- file.path(outFolder, paste0("limma_venn_diag.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
vennDiagram(results)
mtext(side=3, text="limma")
dev.off()
cat(paste0("... written: ", outFile, "\n"))


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
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
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
genesymb2 <- "SNAI1"

if(genesymb2 %in% rownames(brca_data_filtered) & genesymb1 %in% rownames(brca_data_filtered)) {
  stopifnot(genesymb1 %in% ens2genes)
  gene_id1 <- names(ens2genes)[ens2genes == genesymb1]
  stopifnot(length(gene_id1) == 1)
  stopifnot(gene_id1 %in% rownames(brca_data_filtered))
  
  stopifnot(genesymb2 %in% ens2genes)
  gene_id2 <- names(ens2genes)[ens2genes == genesymb2]
  stopifnot(length(gene_id2) == 1)
  stopifnot(gene_id2 %in% rownames(brca_data_filtered))
  
  plot_dt <- data.frame(t(brca_data_filtered[c(gene_id1,gene_id2),]))
  colnames(plot_dt) <- c("gene1", "gene2")
  stopifnot(rownames(plot_dt) == c(lumA_samples, lumB_samples))
  plot_dt$PAM50 <- c(rep("lumA", nLumA), rep("lumB", nLumB))
  
  p <- ggplot(plot_dt, aes(x = gene1, y = gene2, color = PAM50)) +
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
  
  outFile <- file.path(outFolder, paste0(genesymb2, "_vs_", genesymb1, "_expr_filtered_data.", plotType))
  ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}  else {
  stopifnot(genesymb1 %in% ens2genes)
  gene_id1 <- names(ens2genes)[ens2genes == genesymb1]
  stopifnot(length(gene_id1) == 1)
  stopifnot(gene_id1 %in% rownames(brca_data_notFiltered))
  
  stopifnot(genesymb2 %in% ens2genes)
  gene_id2 <- names(ens2genes)[ens2genes == genesymb2]
  stopifnot(length(gene_id2) == 1)
  stopifnot(gene_id2 %in% rownames(brca_data_notFiltered))
  
  plot_dt <- data.frame(t(brca_data_notFiltered[c(gene_id1,gene_id2),]))
  colnames(plot_dt) <- c("gene1", "gene2")
  stopifnot(rownames(plot_dt) == c(lumA_samples, lumB_samples))
  plot_dt$PAM50 <- c(rep("lumA", nLumA), rep("lumB", nLumB))
  
  p <- ggplot(plot_dt, aes(x = gene1, y = gene2, color = PAM50)) +
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
  
  outFile <- file.path(outFolder, paste0(genesymb2, "_vs_", genesymb1, "_expr_notFiltered_data.", plotType))
  ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
}



##################  <<<<<<<<<<<<<<< FIG 6b

# download from here https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_ANGIOGENESIS.html 11/9/21
angiogenesis_genes_all <- read.delim("HALLMARK_ANGIOGENESIS.txt", skip=2, header=F)[,1]
# gene_dt=get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))
# brca_phenopath_fit=get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))
# brca_data_filteredT=t(get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_data_filtered.Rdata")))

# df_beta <- data_frame(
#   c = pcavi$m_lambda,
#   beta = pcavi$m_beta[1,],
#   chi = pcavi$chi_exp[1,],
#   alpha = pcavi$m_alpha[1,],
#   pos_sd = sqrt(pcavi$s_beta[1,]),
#   gene = rownames(sce),
#   is_sig = as.vector(significant_interactions(pcavi, n = 3)),
#   hgnc_symbol = rowData(sce)$hgnc_symbol,
#   ensembl_gene_id = rowData(sce)$ensembl_gene_id
# )
int_dt <- as.data.frame(interactions(brca_phenopath_fit))
df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = colnames(brca_data_filteredT),
                      stringsAsFactors = FALSE)
stopifnot(df_beta$beta == int_dt$interaction_effect_size)
stopifnot(df_beta$gene == int_dt$feature)
stopifnot(df_beta$gene %in% names(ens2genes))
df_beta$lambda <- int_dt$pathway_loading
df_beta$gene_symb <- ens2genes[paste0(df_beta$gene)]
df_beta$is_angiogenesis <- df_beta$gene_symb %in% angiogenesis_genes_all
cat(paste0("angio genes av.: ", sum(df_beta$is_angiogenesis), "/", length(angiogenesis_genes_all), "\n"))
# 19/36
growth_factor_grep <- c("^FGF", "^VEGF", "^EGF", "^TGF")
growth_factors <- unlist(lapply(growth_factor_grep, grep, df_beta$gene_symb, value = TRUE))
growth_factors <- growth_factors[growth_factors != "EFGLAM"]

df_beta$is_gf <- df_beta$gene_symb  %in% growth_factors
cat(paste0("gowth factors av.: ", sum(df_beta$is_gf), "/", length(growth_factors), "\n"))
# 32/32 by definition...

# takes the top 10 of the growth factors sorted by lambda
ag_to_plot <- filter(df_beta, is_gf) %>%
  arrange(desc(lambda)) %>%
  head(n = 10) %>%
  .$gene_symb
if("VEGF" %in% ens2genes)
  ag_to_plot <- c(ag_to_plot, "VEGFC")

# retrieve the ones to plot
stopifnot(rownames(brca_data_filtered) %in% names(ens2genes))
pp_genes <- ens2genes[paste0(rownames(brca_data_filtered))]
stopifnot(!is.na(pp_genes))
stopifnot(rownames(brca_data_filtered) == names(pp_genes))
stopifnot(ag_to_plot %in% pp_genes)

ag_inds <- match(ag_to_plot, pp_genes)
stopifnot(pp_genes[ag_inds] == ag_to_plot)

marker_mat <- t(brca_data_filtered)[, ag_inds, drop=FALSE]
stopifnot(pp_genes[colnames(marker_mat)] == ag_to_plot)
colnames(marker_mat) <- ag_to_plot
marker_df <-  marker_mat %>% 
  as_data_frame() %>% 
  dplyr::mutate(pseudotime = trajectory(brca_phenopath_fit))

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
p <- ggplot(marker_df_2, aes(x = pseudotime_order, y = gene, fill = norm_expression)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_viridis(name = "Expression\nz-score") +
  # scale_x_continuous(expand = expansion(mult = c(1,1)))+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.background=element_rect(color="white"),
        legend.position = "right",
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

outFile <- file.path(outFolder, paste0("angio_genes_ptime_order.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



                                              # 
                                              # 
                                              # ###>>>>>>> fig 6a barplot
                                              # genes <- c("FGF2", "FBP1", "FOXC1")
                                              # stopifnot(genes %in% ens2genes)
                                              # genes_ens <- names(ens2genes)[ens2genes %in% genes]
                                              # stopifnot(genes_ens %in% rownames(brca_data_filtered))
                                              # 
                                              # lumB_samples <- tcga_annot_dt$cgc_sample_id[tcga_annot_dt$PAM50 == "lumB"]
                                              # lumA_samples <- tcga_annot_dt$cgc_sample_id[tcga_annot_dt$PAM50 == "lumA"]
                                              # 
                                              # stopifnot(colnames(brca_data_filtered) == c(lumA_samples, lumB_samples))
                                              # stopifnot(rownames(brca_data_filteredT) == c(lumA_samples, lumB_samples))
                                              # brca_pts <-  trajectory(brca_phenopath_fit)
                                              # names(brca_pts) <- c(lumA_samples, lumB_samples)
                                              # 
                                              # selected_expr_dt <- brca_data_filtered[genes_ens,]
                                              # m_sexpr_dt <- melt(t(selected_expr_dt))
                                              # colnames(m_sexpr_dt) <- c("samp", "gene","logexpr")
                                              # m_sexpr_dt$geneSymb <- ens2genes[as.character(m_sexpr_dt$gene)]
                                              # m_sexpr_dt$PAM50 <- ifelse(m_sexpr_dt$samp %in% lumB_samples, "lumB",
                                              #                                ifelse(m_sexpr_dt$samp %in% lumA_samples, "lumA",NA))
                                              # m_sexpr_dt$pseudotime <- brca_pts[m_sexpr_dt$samp]
                                              # stopifnot(!is.na(m_sexpr_dt))
                                              # 
                                              # # thresholds from Campbell - modified to have different "end" for comparison
                                              # # and modified for using apply (faster)
                                              # classify_cell <- function(row) {
                                              #   if(row['pseudotime'] > 0.5) return("end")
                                              #   if(row['pseudotime'] < -0.3 && row['PAM50'] == "lumB") return('beginning-lumB')
                                              #   if(row['pseudotime'] < -0.3 && row['PAM50'] == "lumA") return("beginning-lumA")
                                              #   return(NA)
                                              # }
                                              # classify_cell_mz <- function(row) {
                                              #   if(row['pseudotime'] > 0.5 && row['PAM50'] == "lumA") return('end-lumA')
                                              #   if(row['pseudotime'] > 0.5 && row['PAM50'] == "lumB") return('end-lumB')
                                              #   if(row['pseudotime'] < -0.3 && row['PAM50'] == "lumB") return('beginning-lumB')
                                              #   if(row['pseudotime'] < -0.3 && row['PAM50'] == "lumA") return("beginning-lumA")
                                              #   return(NA)
                                              # }
                                              # 
                                              # m_sexpr_dt$cell_Class <- sapply(seq_len(nrow(m_sexpr_dt)), function(i) classify_cell(m_sexpr_dt[i, ]))
                                              # m_sexpr_dt$cell_ClassMZ <- sapply(seq_len(nrow(m_sexpr_dt)), function(i) classify_cell_mz(m_sexpr_dt[i, ]))
                                              # 
                                              # stopifnot(grepl("beginning", m_sexpr_dt$cell_Class) == grepl("beginning",m_sexpr_dt$cell_ClassMZ ))
                                              # stopifnot(grepl("end", m_sexpr_dt$cell_Class) == grepl("end",m_sexpr_dt$cell_ClassMZ ))
                                              # 
                                              # ## agg mean expression
                                              # mean_byClass_dt <- aggregate(logexpr ~ geneSymb + cell_Class , data=m_sexpr_dt, FUN=mean)
                                              # mean_byClassMZ_dt <- aggregate(logexpr ~ geneSymb + cell_ClassMZ , data=m_sexpr_dt, FUN=mean)
# > zsc_byClass_dt <- aggregate(logexpr ~ geneSymb+cell_class, data = m_sexpr_dt, 
#                               +                              function(x) mean(winsorize((x - mean(x)) / sd(x))))
# > 
#   > zsc_byClassMZ_dt <- aggregate(logexpr ~ geneSymb+cell_class_mz, data = m_sexpr_dt, 
#                                   +                              function(x) mean(winsorize((x - mean(x)) / sd(x))))
# > range(zsc_byClass_dt$logexpr)
# [1] -0.04195008  0.02286164
# z-score not done, not sure which across which dimension z-score should have been taken

                                                # aggMet="mean"
                                                # myclass="Class"
                                                # for(aggMet in c("mean")) {
                                                #   for(myclass in c("Class", "ClassMZ")) {
                                                #     
                                                #     mylab <- ifelse(aggMet == "mean", "Mean", ifelse(aggMet=="zsc", "Mean z-score", NA))
                                                #     stopifnot(!is.na(mylab))
                                                #     
                                                #     plot_dt <- eval(parse(text = paste0(aggMet, "_by", myclass, "_dt")))
                                                #     
                                                #     plot_dt$geneSymb <- factor(plot_dt$geneSymb, levels=c("FGF2", "FBP1", "FOXC1"))
                                                #     stopifnot(!is.na(plot_dt$geneSymb))
                                                #     
                                                #     p <- ggplot(plot_dt, aes(x = geneSymb, y = logexpr, fill = geneSymb)) +
                                                #       geom_bar(stat = "identity", color = 'grey20') +
                                                #       scale_fill_brewer(palette = "Dark2") +
                                                #       # facet_wrap(~cell_class)+ 
                                                #       facet_wrap(c(paste0("cell_", myclass)), scales="free")+
                                                #       scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
                                                #       # scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
                                                #       theme(panel.background = element_blank(),
                                                #             legend.position = "None",
                                                #             axis.line = element_line(),
                                                #             axis.title.x = element_blank(),
                                                #             axis.title.y = element_text(size = 12, face = "bold"),
                                                #             axis.text.y = element_text(size = 10),
                                                #             axis.text.x = element_text(size = 8, angle = -90,hjust=0.5, face = "bold")) +
                                                #       labs(y = paste0(mylab, "\n expression"))
                                                #     
                                                #     if(myclass=="Class"){
                                                #       myWidthGG2 <- myWidthGG*2
                                                #       myHeightGG2 <- myHeightGG*0.7
                                                #     } else {
                                                #         myWidthGG2 <- myWidthGG*1.5
                                                #         myHeightGG2 <- myHeightGG*1.3
                                                #       }
                                                #     
                                                #     outFile <- file.path(outFolder, paste0("selectedGenes_", aggMet, "Expr_",myclass , "_barplot.", plotType))
                                                #     ggsave(p, filename = outFile, height=myHeightGG2, width=myWidthGG2)
                                                #     cat(paste0("... written: ", outFile, "\n"))
                                                #     
                                                #     
                                                #   }}


stopifnot(brca_phenopath_fit$m_beta[1,] == df_beta$beta)
stopifnot(brca_phenopath_fit$feature_names == df_beta$gene)

df_beta$alpha <- brca_phenopath_fit$m_alpha[1,]
df_beta$crossover <- -df_beta$alpha/df_beta$beta
stopifnot(!is.na(df_beta))
df_beta_sig <- df_beta[df_beta$is_sig,]

mytit <- paste0("Crossover points distribution")
subtit <- paste0("only signif. genes (n=", nrow(df_beta_sig), ")")
med_val <- median(df_beta_sig$crossover)

p <- ggplot(df_beta_sig, aes(x = crossover)) + 
  ggtitle(paste0(mytit), subtitle=paste0(subtit))+
geom_histogram(fill = "#74a9cf", color = "grey90", bins = 30) +
xlab("Crossover point") + ylab("Number of genes") + 
  geom_vline(xintercept = med_val, linetype=2)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
  # plot.title = element_text(hjust=0.5),
  # plot.subtitle = element_text(hjust=0.5),
  panel.background = element_blank(),
      axis.line=element_line(),
  axis.text = element_text(size = 9),
      axis.title = element_text(size = 10))

max_val <- max(ggplot_build(p)$data[[1]]$count)
p <- p + annotate("text", label=paste0("median = ", round(med_val, 3)),
                  x=median(df_beta_sig$crossover) -0.1,
                  y = max_val, hjust=1)

outFile <- file.path(outFolder, paste0("signifGenes_crossoverPointsDistribution.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.8)
cat(paste0("... written: ", outFile, "\n"))









#~~~~~ survival analysis

brca_data_filtered <- get(load(file.path(outFolder, "brca_data_filtered.Rdata")))
brca_data_filteredT <- t(brca_data_filtered)
brca_phenopath_fit <- get(load(file.path(outFolder, "brca_phenopath_fit.Rdata")))
brca_pseudotimes <-  trajectory(brca_phenopath_fit)

library(RTCGA)
library(RTCGA.clinical)

stopifnot(!duplicated(tcga_annot_dt$patient_barcode))

surv_dt <- survivalTCGA(BRCA.clinical, 
                        extract.cols="admin.disease_code")
# Show the first few lines
head(surv_dt)

tcga_annot_dt$labs_for_survival <- gsub("(^.+?-.+?-.+?)-.+", "\\1", tcga_annot_dt$cgc_sample_id)
stopifnot(tcga_annot_dt$labs_for_survival %in% surv_dt$bcr_patient_barcode)
stopifnot(tcga_annot_dt$patient_barcode == tcga_annot_dt$labs_for_survival)
stopifnot(!duplicated(tcga_annot_dt$labs_for_survival))

all_lumStatus <- ifelse(tcga_annot_dt$PAM50 == "lumA", "lumA",
                       ifelse(tcga_annot_dt$PAM50 == "lumB", "lumB", NA))
stopifnot(!is.na(all_lumStatus))
all_lumStatus <- setNames(all_lumStatus, tcga_annot_dt$labs_for_survival)
surv_dt <- surv_dt[surv_dt$bcr_patient_barcode %in% tcga_annot_dt$labs_for_survival,]
stopifnot(nrow(surv_dt) == nrow(tcga_annot_dt))

# add pseudotime info
all_ptimes <- brca_pseudotimes
stopifnot(rownames(brca_data_filteredT) == tcga_annot_dt$cgc_sample_id)
names(all_ptimes) <- tcga_annot_dt$labs_for_survival

surv_dt$pseudotime <- all_ptimes[paste0(surv_dt$bcr_patient_barcode)]
stopifnot(!is.na(surv_dt$pseudotime))

surv_dt$PAM50 <- all_lumStatus[paste0(surv_dt$bcr_patient_barcode)]
stopifnot(!is.na(surv_dt$PAM50))

### TAKE ONLY THE pseudotime > 0 ???
surv_dt <- surv_dt[surv_dt$pseudotime > 0,]

# let’s run a Cox PH model
# By default it’s going to treat lumB cancer as the baseline, because alphabetically it’s first.
coxph(Surv(times, patient.vital_status)~PAM50 + pseudotime, data=surv_dt)
# This tells us that compared to the baseline lumB group, lumA have ~0.7x increase in hazards, 
# and pseuodtime 1.12x worse survival. 
# Let’s create a survival curve, visualize it with a Kaplan-Meier plot, and show a table for the first 5 years survival rates.
sfit <- survfit(Surv(times, patient.vital_status)~PAM50 + pseudotime, data=surv_dt)
# time_range <- seq(0,800,2000)
time_range <- seq(from=min(surv_dt$times),to=max(surv_dt$times), by=200)
summary(sfit, times=time_range)

surv_dt$pt_bin <- cut(surv_dt$pseudotime, breaks=seq(0,1,0.25))
stopifnot(!is.na(surv_dt$pt_bin))

fit1 <- survfit(Surv(times,patient.vital_status) ~ pt_bin + strata(PAM50), data = surv_dt)
ggsurvplot(fit1, data = surv_dt, risk.table = FALSE)

fit2 <- survfit(Surv(times,patient.vital_status) ~ pt_bin + PAM50, data = surv_dt)
ggsurvplot(fit2, data = surv_dt, risk.table = FALSE)

# I don't understand the difference...

p <- ggsurvplot_facet(fit2, data = surv_dt, 
                      conf.int = TRUE, 
                      pval=TRUE,
                      legend="bottom",
                      legend.title="Pseudotime",
                      legend.labs =as.character(levels(surv_dt$pt_bin)),
                      facet.by="PAM50",
                      risk.table = FALSE)



outFile <- file.path(outFolder, paste0("survival_PAM50_by_pseudotime.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



############## REACTOME PATHWAY ENRICHMENT analysis  ##############
# *a pathway enrichment analysis using Reactome to discover whether any of the top 20 interacting genes (by β value) *
# not sure if take absolute beta in the article ??
# brca_data_filtered <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_data_filtered.Rdata"))
# brca_data_filteredT <- t(brca_data_filtered)
# brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))

gene_names <- colnames(brca_data_filteredT)
df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = gene_names)

###### FIRST WAY TO GET MATCHING IDS

httr::set_config(httr::config(ssl_verifypeer = FALSE))  ### added to access ensembl biomart connection

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes_entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype"),
                            filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
dim(listOfGenes_entrez)
# [1] 23069     3
stopifnot(df_beta$gene %in% listOfGenes_entrez$ensembl_gene_id)
listOfGenes_entrez <- listOfGenes_entrez[listOfGenes_entrez$ensembl_gene_id %in% df_beta$gene,]
dim(listOfGenes_entrez)
# [1] 6005    3
listOfGenes_entrez <- listOfGenes_entrez[!is.na(listOfGenes_entrez$entrezgene_id),]
dim(listOfGenes_entrez)
# 5982 3
# I need to have unique ensemblID
listOfGenes_entrez <- listOfGenes_entrez[!duplicated(listOfGenes_entrez$ensembl_gene_id),]
dim(listOfGenes_entrez)
# [1] 5982 3
any(duplicated(listOfGenes_entrez$entrezgene_id))
# TRUE
stopifnot(!duplicated(listOfGenes_entrez$ensembl_gene_id))
stopifnot(!is.na(listOfGenes_entrez$ensembl_gene_id))
listOfGenes_entrez1 <- setNames(listOfGenes_entrez$entrezgene_id, listOfGenes_entrez$ensembl_gene_id)
stopifnot(!is.na(listOfGenes_entrez1))

###### 2ND WAY TO GET MATCHING IDS
#columns(org.Hs.eg.db) # returns list of available keytypes
listOfGenes_entrez2 <- mapIds(org.Hs.eg.db,
                              keys=df_beta$gene, #Column containing Ensembl gene ids
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
listOfGenes_entrez2 <- listOfGenes_entrez2[!is.na(listOfGenes_entrez2)]
stopifnot(!is.na(listOfGenes_entrez2))
length(listOfGenes_entrez2)
# 5996
any(duplicated(listOfGenes_entrez2))
# TRUE # also duplicated entrez

stopifnot(names(listOfGenes_entrez1) %in% df_beta$gene)
stopifnot(names(listOfGenes_entrez2) %in% df_beta$gene)
stopifnot(!is.na(listOfGenes_entrez2))
stopifnot(!is.na(listOfGenes_entrez1))

# if one is longer, take the longest, otherwise take the one with least duplicated entrezIDs (less ambiguous)
if(length(listOfGenes_entrez1) != length(listOfGenes_entrez2)) {
  if(length(listOfGenes_entrez1) > length(listOfGenes_entrez2)) {
    ens2entrez <- listOfGenes_entrez1
    cat(paste0("... take list 1\n"))
  } else {
    ens2entrez <- listOfGenes_entrez2
    cat(paste0("... take list 2\n"))
  }
} else {
  if(length(unique(listOfGenes_entrez1)) > length(unique(listOfGenes_entrez2))) {
    ens2entrez <- listOfGenes_entrez1
    cat(paste0("... take list 1\n"))
  } else { # if == take list 2 (do not know better way)
    ens2entrez <- listOfGenes_entrez2
    cat(paste0("... take list 2\n"))
  }
}
cat(paste0("av. genes with matched entrez ID:", length(ens2entrez), "/", nrow(df_beta), "\n"))
cat(paste0("# duplicated entrezIDs: ", sum(duplicated(ens2entrez)), "\n"))

df_beta_entrez <- df_beta[df_beta$gene %in% names(ens2entrez),]
stopifnot(paste0(df_beta_entrez$gene) %in% names(ens2entrez))
df_beta_entrez$gene_entrez <- ens2entrez[paste0(df_beta_entrez$gene)]
stopifnot(!is.na(df_beta_entrez))

top20_posbeta <- df_beta_entrez[order(df_beta_entrez$beta, decreasing = TRUE),][1:20,"gene_entrez"]

top20_negbeta <- df_beta_entrez[order(df_beta_entrez$beta, decreasing = FALSE),][1:20,"gene_entrez"]

top20_absbeta <- df_beta_entrez[order(abs(df_beta_entrez$beta), decreasing = TRUE),][1:20,"gene_entrez"]

# will not use the cutoff because no one signif, but still want to look at the top ones
# reactome_pvaluecutoff <- 0.05
# from here: https://www.biostars.org/p/434669/
# seems ok to leave universe param out
cat(paste0("... run Reactome PA\n"))
top_pos_React <- enrichPathway(gene=top20_posbeta, pvalueCutoff = 1 , readable=TRUE)
top_neg_React <- enrichPathway(gene=top20_negbeta, pvalueCutoff = 1 , readable=TRUE)
top_abs_React <- enrichPathway(gene=top20_absbeta, pvalueCutoff = 1 , readable=TRUE)
cat(paste0("...... finished\n"))

keep_cols <- c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID")
out_topPos <- top_pos_React@result[,keep_cols]
out_topPos$pvalue <- round(out_topPos$pvalue, 4)
out_topPos$p.adjust <- round(out_topPos$p.adjust, 4)
out_topPos$qvalue <- round(out_topPos$qvalue, 4)

outFile <- file.path(outFolder, "reactomeEnrichment_topPosBeta.txt")
write.table(out_topPos, file = outFile, sep="\t", col.names=T, row.names=F, quote=F)
cat(paste0("... written: ", outFile, "\n"))

out_topNeg <- top_neg_React@result[,keep_cols]
out_topNeg$pvalue <- round(out_topNeg$pvalue, 4)
out_topNeg$p.adjust <- round(out_topNeg$p.adjust, 4)
out_topNeg$qvalue <- round(out_topNeg$qvalue, 4)

outFile <- file.path(outFolder, "reactomeEnrichment_topNegBeta.txt")
write.table(out_topNeg, file = outFile, sep="\t", col.names=T, row.names=F, quote=F)
cat(paste0("... written: ", outFile, "\n"))


out_topAbs <- top_abs_React@result[,keep_cols]
out_topAbs$pvalue <- round(out_topAbs$pvalue, 4)
out_topAbs$p.adjust <- round(out_topAbs$p.adjust, 4)
out_topAbs$qvalue <- round(out_topAbs$qvalue, 4)
outFile <- file.path(outFolder, "reactomeEnrichment_topAbsBeta.txt")
write.table(out_topAbs, file = outFile, sep="\t", col.names=T, row.names=F, quote=F)
cat(paste0("... written: ", outFile, "\n"))






##################################################
cat(paste0("****** DONE"))
cat(paste0(startTime, " - ", Sys.time(),  "\n"))
stop("--ok\n")

####################################################################################################
### THRASH
####################################################################################################
dt1=get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/limma_merged_dt.Rdata"))
dt2=get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/merged_dt.Rdata"))
stopifnot(dt1$feature==dt2$ensemblID)
stopifnot(round(dt1$pval,4)==round(dt2$P.Value,4))
stopifnot(round(dt1$qval,4)==round(dt2$adj.P.Val,4))

############## survival analysis  ##############
# brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))
# brca_pseudotimes <-  trajectory(brca_phenopath_fit)
# brca_data_filtered <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_data_filtered.Rdata"))
# brca_data_filteredT <- t(brca_data_filtered)
#tcga_annot_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_BRCA_RECOUNT2/tcga_sampleAnnot.Rdata"))
library(RTCGA)
library(RTCGA.clinical)

stopifnot(!duplicated(tcga_annot_dt$patient_barcode))

surv_dt <- survivalTCGA(BRCA.clinical, 
                     extract.cols="admin.disease_code")
# Show the first few lines
head(surv_dt)

tcga_annot_dt$labs_for_survival <- gsub("(^.+?-.+?-.+?)-.+", "\\1", tcga_annot_dt$cgc_sample_id)
stopifnot(tcga_annot_dt$labs_for_survival %in% surv_dt$bcr_patient_barcode)
stopifnot(tcga_annot_dt$patient_barcode == tcga_annot_dt$labs_for_survival)
stopifnot(!duplicated(tcga_annot_dt$labs_for_survival))

all_lumStatus <- ifelse(tcga_annot_dt$PAM50 == "lumA", "lumA",
                       ifelse(tcga_annot_dt$PAM50 == "lumB", "lumB", NA))
stopifnot(!is.na(all_lumStatus))
all_lumStatus <- setNames(all_lumStatus, tcga_annot_dt$labs_for_survival)
surv_dt <- surv_dt[surv_dt$bcr_patient_barcode %in% tcga_annot_dt$labs_for_survival,]
stopifnot(nrow(surv_dt) == nrow(tcga_annot_dt))

# add pseudotime info
all_ptimes <- brca_pseudotimes
stopifnot(rownames(brca_data_filteredT) == tcga_annot_dt$cgc_sample_id)
names(all_ptimes) <- tcga_annot_dt$labs_for_survival

surv_dt$pseudotime <- all_ptimes[paste0(surv_dt$bcr_patient_barcode)]
stopifnot(!is.na(surv_dt$pseudotime))

surv_dt$PAM50 <- all_lumStatus[paste0(surv_dt$bcr_patient_barcode)]
stopifnot(!is.na(surv_dt$PAM50))

### TAKE ONLY THE pseudotime > 0 ???
surv_dt <- surv_dt[surv_dt$pseudotime > 0,]

# let’s run a Cox PH model
# By default it’s going to treat lumB cancer as the baseline, because alphabetically it’s first.
coxph(Surv(times, patient.vital_status)~PAM50 + pseudotime, data=surv_dt)
# This tells us that compared to the baseline lumB group, lumA have ~0.7x increase in hazards, 
# and pseuodtime 1.12x worse survival. 
# Let’s create a survival curve, visualize it with a Kaplan-Meier plot, and show a table for the first 5 years survival rates.
sfit <- survfit(Surv(times, patient.vital_status)~PAM50 + pseudotime, data=surv_dt)
time_range <- seq(0,800,2000)
summary(sfit, times=time_range)

surv_dt$pt_bin <- cut(surv_dt$pseudotime, breaks=seq(0,1,0.25))
stopifnot(!is.na(surv_dt$pt_bin))

fit1 <- survfit(Surv(times,patient.vital_status) ~ pt_bin + strata(PAM50), data = surv_dt)
ggsurvplot(fit1, data = surv_dt, risk.table = FALSE)

fit2 <- survfit(Surv(times,patient.vital_status) ~ pt_bin + PAM50, data = surv_dt)
ggsurvplot(fit2, data = surv_dt, risk.table = FALSE)

# I don't understand the difference...

p <- ggsurvplot_facet(fit2, data = surv_dt, 
                 conf.int = TRUE, 
                 pval=TRUE,
                 legend="bottom",
                 legend.title="Pseudotime",
                 legend.labs =as.character(levels(surv_dt$pt_bin)),
                 facet.by="PAM50",
                 risk.table = FALSE)



outFile <- file.path(outFolder, paste0("survival_PAM50_by_pseudotime.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



############## REACTOME PATHWAY ENRICHMENT analysis  ##############
# *a pathway enrichment analysis using Reactome to discover whether any of the top 20 interacting genes (by β value) *
# not sure if take absolute beta in the article ??
# brca_data_filtered <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_data_filtered.Rdata"))
# brca_data_filteredT <- t(brca_data_filtered)
# brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))

gene_names <- colnames(brca_data_filteredT)
df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = gene_names)

###### FIRST WAY TO GET MATCHING IDS
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes_entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype"),
                     filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
dim(listOfGenes_entrez)
# [1] 23069     3
stopifnot(df_beta$gene %in% listOfGenes_entrez$ensembl_gene_id)
listOfGenes_entrez <- listOfGenes_entrez[listOfGenes_entrez$ensembl_gene_id %in% df_beta$gene,]
dim(listOfGenes_entrez)
# [1] 6005    3
listOfGenes_entrez <- listOfGenes_entrez[!is.na(listOfGenes_entrez$entrezgene_id),]
dim(listOfGenes_entrez)
# 5982 3
# I need to have unique ensemblID
listOfGenes_entrez <- listOfGenes_entrez[!duplicated(listOfGenes_entrez$ensembl_gene_id),]
dim(listOfGenes_entrez)
# [1] 5982 3
any(duplicated(listOfGenes_entrez$entrezgene_id))
# TRUE
stopifnot(!duplicated(listOfGenes_entrez$ensembl_gene_id))
stopifnot(!is.na(listOfGenes_entrez$ensembl_gene_id))
listOfGenes_entrez1 <- setNames(listOfGenes_entrez$entrezgene_id, listOfGenes_entrez$ensembl_gene_id)
stopifnot(!is.na(listOfGenes_entrez1))

###### 2ND WAY TO GET MATCHING IDS
#columns(org.Hs.eg.db) # returns list of available keytypes
listOfGenes_entrez2 <- mapIds(org.Hs.eg.db,
                     keys=df_beta$gene, #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
listOfGenes_entrez2 <- listOfGenes_entrez2[!is.na(listOfGenes_entrez2)]
stopifnot(!is.na(listOfGenes_entrez2))
length(listOfGenes_entrez2)
# 5996
any(duplicated(listOfGenes_entrez2))
# TRUE # also duplicated entrez

stopifnot(names(listOfGenes_entrez1) %in% df_beta$gene)
stopifnot(names(listOfGenes_entrez2) %in% df_beta$gene)
stopifnot(!is.na(listOfGenes_entrez2))
stopifnot(!is.na(listOfGenes_entrez1))

# if one is longer, take the longest, otherwise take the one with least duplicated entrezIDs (less ambiguous)
if(length(listOfGenes_entrez1) != length(listOfGenes_entrez2)) {
  if(length(listOfGenes_entrez1) > length(listOfGenes_entrez2)) {
    ens2entrez <- listOfGenes_entrez1
    cat(paste0("... take list 1\n"))
  } else {
    ens2entrez <- listOfGenes_entrez2
    cat(paste0("... take list 2\n"))
  }
} else {
  if(length(unique(listOfGenes_entrez1)) > length(unique(listOfGenes_entrez2))) {
    ens2entrez <- listOfGenes_entrez1
    cat(paste0("... take list 1\n"))
  } else { # if == take list 2 (do not know better way)
    ens2entrez <- listOfGenes_entrez2
    cat(paste0("... take list 2\n"))
  }
}
cat(paste0("av. genes with matched entrez ID:", length(ens2entrez), "/", nrow(df_beta), "\n"))
cat(paste0("# duplicated entrezIDs: ", sum(duplicated(ens2entrez)), "\n"))

df_beta_entrez <- df_beta[df_beta$gene %in% names(ens2entrez),]
stopifnot(paste0(df_beta_entrez$gene) %in% names(ens2entrez))
df_beta_entrez$gene_entrez <- ens2entrez[paste0(df_beta_entrez$gene)]
stopifnot(!is.na(df_beta_entrez))

top20_posbeta <- df_beta_entrez[order(df_beta_entrez$beta, decreasing = TRUE),][1:20,"gene_entrez"]

top20_negbeta <- df_beta_entrez[order(df_beta_entrez$beta, decreasing = FALSE),][1:20,"gene_entrez"]

top20_absbeta <- df_beta_entrez[order(abs(df_beta_entrez$beta), decreasing = TRUE),][1:20,"gene_entrez"]

# will not use the cutoff because no one signif, but still want to look at the top ones
# reactome_pvaluecutoff <- 0.05
# from here: https://www.biostars.org/p/434669/
# seems ok to leave universe param out
cat(paste0("... run Reactome PA\n"))
top_pos_React <- enrichPathway(gene=top20_posbeta, pvalueCutoff = 1 , readable=TRUE)
top_neg_React <- enrichPathway(gene=top20_negbeta, pvalueCutoff = 1 , readable=TRUE)
top_abs_React <- enrichPathway(gene=top20_absbeta, pvalueCutoff = 1 , readable=TRUE)
cat(paste0("...... finished\n"))

keep_cols <- c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID")
out_topPos <- top_pos_React@result[,keep_cols]
out_topPos$pvalue <- round(out_topPos$pvalue, 4)
out_topPos$p.adjust <- round(out_topPos$p.adjust, 4)
out_topPos$qvalue <- round(out_topPos$qvalue, 4)

outFile <- file.path(outFolder, "reactomeEnrichment_topPosBeta.txt")
write.table(out_topPos, file = outFile, sep="\t", col.names=T, row.names=T, quote=F)
cat(paste0("... written: ", outFile, "\n"))

out_topNeg <- top_neg_React@result[,keep_cols]
out_topNeg$pvalue <- round(out_topNeg$pvalue, 4)
out_topNeg$p.adjust <- round(out_topNeg$p.adjust, 4)
out_topNeg$qvalue <- round(out_topNeg$qvalue, 4)

outFile <- file.path(outFolder, "reactomeEnrichment_topNegBeta.txt")
write.table(out_topNeg, file = outFile, sep="\t", col.names=T, row.names=T, quote=F)
cat(paste0("... written: ", outFile, "\n"))


out_topAbs <- top_abs_React@result[,keep_cols]
out_topAbs$pvalue <- round(out_topAbs$pvalue, 4)
out_topAbs$p.adjust <- round(out_topAbs$p.adjust, 4)
out_topAbs$qvalue <- round(out_topAbs$qvalue, 4)
outFile <- file.path(outFolder, "reactomeEnrichment_topAbsBeta.txt")
write.table(out_topAbs, file = outFile, sep="\t", col.names=T, row.names=T, quote=F)
cat(paste0("... written: ", outFile, "\n"))


# geneRatio = ratio of input genes that are annotated in a term
# BgRatio = ratio of all genes that are annotated in this term
# https://www.sciencedirect.com/science/article/pii/S2666675821000667
# fold enrichment = defined as the ratio of the frequency of input genes annotated in a term to the frequency
# of all genes annotated to that term, and it is easy to calculate by dividing geneRatio by BgRatio
# 
# Let is suppose I have a collection of genesets called : HALLMARK
# Now let is suppose there is a specific geneset there called: E2F_targets
# 
# BgRatio = M/N.
# 
# M = size of the geneset (eg size of the E2F_targets); 
# (# of genes within that distribution that are annotated (either directly or indirectly) to the node of interest
#   # of genes from the universe annotated to gene set
# 
# N = size of all of the unique genes in the collection of genesets (example the HALLMARK collection);
# # (is the total number of genes in the background distribution (universe)
# 
#  GeneRatio is k/n.
# k = size of the overlap of 'a vector of gene id' you input with the specific geneset (eg E2F_targets),
# only unique genes; (the number of genes within that list n, which are annotated to the node.
#                                                                                                                                                                                                                                  
# n = size of the overlap of 'a vector of gene id' you input with all the members of 
# the collection of genesets (eg the HALLMARK collection),only unique genes;
# is the size of the list of genes of interest 

# BgRatio = M/N => size given set set/size universe
# GeneRatio = k/n => # input genes in a given set/# input genes in universe


### look at the pathways they looked at
pathways <- c("1643713", "1226099",
              "2219528", "2644603",
              "3304351", "4791275")
pathways <- paste0("R-HSA-", pathways)
id_to_name <- as.list(reactomePATHID2NAME)
pathway_names <- id_to_name[pathways]
pathways_to_genes <- as.list(reactomePATHID2EXTID)
gene_list <- pathways_to_genes[pathways]
pathway_names <- sapply(pathway_names, function(pn) {
  gsub("Homo sapiens: ", "", pn)
})
# R-HSA-1643713                                      R-HSA-1226099 
# "Signaling by EGFR in Cancer"                     "Signaling by FGFR in disease" 
# R-HSA-2219528                                      R-HSA-2644603 
# "PI3K/AKT Signaling in Cancer"                    "Signaling by NOTCH1 in Cancer" 
# R-HSA-3304351                                      R-HSA-4791275 
# "Signaling by TGF-beta Receptor Complex in Cancer"                       "Signaling by WNT in cancer" 

# I don't have these pathways in my analyses... not conclusive
na.omit(top_pos_React[pathways,])
# ID                   Description GeneRatio  BgRatio     pvalue  p.adjust     qvalue
# R-HSA-2644603 R-HSA-2644603 Signaling by NOTCH1 in Cancer      1/10 58/10856 0.05218144 0.1108856 0.07209278
# geneID Count
# R-HSA-2644603  MAML2     1
na.omit(top_neg_React[pathways,])
# 0
na.omit(top_neg_React[pathways,])
# 0


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

```{r crossovlumBgo}
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

```{r crossovlumBgene-plots}
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
df_gex$x <- plyr::mapvalues(df_gex$x, from = sort(unique(df_gex$x)), to = c("lumB", "lumA"))
df_gex$x <- factor(df_gex$x, levels = c("lumB", "lumA"))
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
                         to = c("lumB", "lumA"))
tmap <- sce$tmap
vgf <- c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "FGF2", "CXCL8") 
mm <- match(vgf, rowData(sce_gene)$hgnc_symbol)
vgf_exprs <- t(exprs(sce_gene)[mm, ]) %>% as_data_frame()
names(vgf_exprs) <- vgf
vgf_exprs <- mutate(vgf_exprs, z = tmap, PAM50 = x_str)
vgf_exprs_tidy <- gather(vgf_exprs, gene, expression, -z, -PAM50)
ggplot(vgf_exprs_tidy, aes(x = z, y = expression, color = PAM50)) +
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
brca_exprs <- mutate(brca_exprs, z = tmap, PAM50 = x_str)
brca_exprs_tidy <- gather(brca_exprs, gene, expression, -z, -PAM50)
ggplot(brca_exprs_tidy, aes(x = z, y = expression, color = PAM50)) +
  geom_point(alpha = 0.6) + facet_wrap(~ gene, scales = "free_y") +
  xlab("Pathway score") + ylab("Expression") +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = "ER status") 
ggsave("../../figs/brca_genes.png", width = 9, height = 6)
```




gather() does the reverse of spread(). gather() collects a set of column names and
places them into a single “key” column. It also collects the cells of those columns 
and places them into a single value column. You can use gather() to tidy table4.
table4  # cases

## Source: local data frame [3 x 3]
## 
##       country   1999   2000
## 1 Afghanistan    745   2666
## 2      Brazil  37737  80488
## 3       China 212258 213766

gather(table4, "year", "cases", 2:3)

## Source: local data frame [6 x 3]
## 
##       country year  cases
## 1 Afghanistan 1999    745
## 2      Brazil 1999  37737
## 3       China 1999 212258
## 4 Afghanistan 2000   2666
## 5      Brazil 2000  80488
## 6       China 2000 213766


