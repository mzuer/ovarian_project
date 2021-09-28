
# Rscript phenopath_recount2_johansson_BRCA_lumAlumB.R

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

runPheno <- T

cond1 <- "LumA"
cond2 <- "LumB"


plotType <- "png"
myHeight <- myWidth <- 400
myHeightGG <- 6
myWidthGG <- 6

# filter to keep only the most variables -> try with hard-set threshold 5000
# cannot use same threshold as tcga counts (lower values)
nMostVarToKeep <- 5000


purityFilter <- 0.6


# setwd("~/Documents/FREITAS_LAB/ovarian_project/phenopath")

outFolder <- paste0("PHENOPATH_RECOUNT2_JOHANSSON_BRCA_", toupper(cond1), toupper(cond2))
dir.create(outFolder, recursive=TRUE)


# load input data

brca_data_raw <- read.delim(file.path("../../breast_project/data/johansson_data_relative_ratios_to_pool.csv"), sep=",")
dim(brca_data_raw)
# [1] 9995 47

# I will keep the 5'000 most var genes
prot_values <- brca_data_raw
stopifnot(!duplicated(brca_data_raw$gene_symbol))
stopifnot(!duplicated(brca_data_raw$ensembl_id))
rownames(prot_values) <- prot_values$gene_symbol
prot_values$gene_symbol <- prot_values$ensembl_id <- NULL
stopifnot(prot_values >=0)
prot_values_log2 <- log2(prot_values+1) # phenopath expects log2 tpm+1


prot_vars <- apply(prot_values_log2, 1, var)
stopifnot(!is.na(prot_vars))
prot_vars <- sort(prot_vars, decreasing = TRUE)
tokeep <- names(prot_vars)[1:nMostVarToKeep]
stopifnot(tokeep %in% rownames(prot_values_log2))

prot_values_log2_filt <- prot_values_log2[paste0(tokeep),]
stopifnot(nrow(prot_values_log2_filt) == nMostVarToKeep)


### load annotation
annot_dt <- read.delim(file.path("../../breast_project/data/johansson_tumor_annot.csv"), sep=",")
stopifnot(setequal(annot_dt$Tumor.ID, colnames(prot_values_log2_filt)))

stopifnot(cond1 %in% annot_dt$PAM50.subtype)
stopifnot(cond2 %in% annot_dt$PAM50.subtype)

cond1_samples <- annot_dt$Tumor.ID[annot_dt$PAM50.subtype == cond1]
cond2_samples <- annot_dt$Tumor.ID[annot_dt$PAM50.subtype == cond2]
  
nCond1 <- length(cond1_samples)
nCond2 <- length(cond2_samples)

stopifnot(cond1_samples %in% colnames(prot_values_log2_filt))
stopifnot(cond2_samples %in% colnames(prot_values_log2_filt))

prot_values_log2_filt <- prot_values_log2_filt[,c(cond1_samples, cond2_samples)]
stopifnot(colnames(prot_values_log2_filt) == c(cond1_samples, cond2_samples))

nGenes <- nrow(prot_values_log2_filt)


################################### 
################################### PCA ON NORM + filtered DATA 
################################### 

pca_brca <- prcomp(t(prot_values_log2_filt), scale=TRUE)
pca_brca_lowrepr <- pca_brca$x
stopifnot(nrow(pca_brca_lowrepr) == nCond1+nCond2)

pc1 <- pca_brca_lowrepr[,1]
pc2 <- pca_brca_lowrepr[,2]
pc3 <- pca_brca_lowrepr[,3]

mycolvect <- 1+as.numeric(annot_dt$PAM50.subtype == cond1)

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_12_", cond1, "_", cond2, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,2), summ_dt=summary(pca_brca),
        main="Proteo Johan. BRCA Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("#", cond1,  "=", nCond1, "; #", cond2, "=", nCond2), side=3)
legend("topleft", legend=c(cond2, cond1), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_23_", cond1, "_", cond2, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(2,3), summ_dt=summary(pca_brca),
        main="Proteo Johan. BRCA Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("#", cond1,  "=", nCond1, "; #", cond2, "=", nCond2), side=3)
legend("topleft", legend=c(cond2, cond1), pch=16, col=c(1,2), bty="n")
dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("in_normFilt_data_pca_13_", cond1, "_", cond2, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#dev.off()
pcaplot(pca_dt=pca_brca_lowrepr, pctoplot=c(1,3), summ_dt=summary(pca_brca),
        main="Proteo Johan. BRCA Filt (log2(.+1))",
        col=mycolvect)
mtext(text=paste0("#", cond1,  "=", nCond1, "; #", cond2, "=", nCond2), side=3)
legend("topleft", legend=c(cond2, cond1), pch=16, col=c(1,2), bty="n")
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

# so here give a value of 1 in cond1 and a value of -1 in cond2
# This means the overall pathway loading λ is the average change for normal and tumor
# if β > 0 =>  the gene is more upregulated over pseudotime under tumor  
# if β  < 0 =>   the gene is more upregulated under normal 

prot_values_log2_filtT <- t(prot_values_log2_filt)
stopifnot(nGenes == ncol(prot_values_log2_filtT))
stopifnot(nCond2+nCond1 == nrow(prot_values_log2_filtT))

# phenopath needs a cell-by-gene matrix = N×G matrix of gene expression (N=samples, G = genes)
stopifnot(nrow(prot_values_log2_filtT) == nCond1+nCond2)
stopifnot(ncol(prot_values_log2_filtT) == nGenes)
mycovar <- 2 * ( rownames(prot_values_log2_filtT) %in% cond1_samples) - 1
stopifnot(sum(mycovar== -1) == nCond2)
stopifnot(sum(mycovar== 1) == nCond1)
stopifnot(!is.na(mycovar))

if(runPheno) {
  brca_phenopath_fit <- phenopath(prot_values_log2_filtT, mycovar, elbo_tol = 1e-6, thin = 40)
  
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
add_corr(pc1, brca_pseudotimes, cond_plus1=cond1, cond_minus1=cond2)
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
add_corr(pc2, brca_pseudotimes, cond_plus1=cond1, cond_minus1=cond2)
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
add_corr(pc3, brca_pseudotimes, cond_plus1=cond1, cond_minus1=cond2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

############## LOOK A BIT AT THE BETA VALUES ##############
# density of the beta with color code signifcance

gene_names <- colnames(prot_values_log2_filtT)

df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = gene_names)

if(sum(df_beta$is_sig) > 2) {
  outFile <- file.path(outFolder, paste0("beta_values_distribution.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(density(df_beta$beta[df_beta$is_sig]), col=2)
  lines(density(df_beta$beta), col=1)
  # lines(density(df_beta$beta[!df_beta$is_sig]), col=3)
  legend("topright", legend=c("all dist.", "only signif. dist."), lty=1, col=c(1,2), bty="n")
  mtext(side=3, text=paste0("# signif: ", sum(df_beta$is_sig), "/", nrow(df_beta)))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


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
dfBeta_topGenes$geneSymb <- dfBeta_topGenes$gene
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
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- df_beta$gene

# the same for the lowest interaction effect ?
i_max <- which.max(df_beta$beta)
p <- plot_iGeneExpr(igene= i_max,
              exprdt=prot_values_log2_filtT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, cond1, cond2),
               valuedt=df_beta,
               valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_", df_beta$geneSymb[i_max], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_max,
                    exprdt=prot_values_log2_filtT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, cond1, cond2),
                    valuedt=df_beta,
                    valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_", df_beta$geneSymb[i_max], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

i_min <-which.min(df_beta$beta)

p <- plot_iGeneExpr(igene= i_min,
               exprdt=prot_values_log2_filtT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, cond1, cond2),
               valuedt=df_beta,
                           valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_", df_beta$geneSymb[i_min], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_min,
                    exprdt=prot_values_log2_filtT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, cond1, cond2),
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
stopifnot(df_beta$gene == colnames(prot_values_log2_filtT))
# stopifnot(df_beta$gene %in% names(ens2genes))
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- df_beta$gene

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
                                                  # stopifnot(ens_gs == names(ens2genes)[ens2genes==gs])
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
                                                                          exprdt=prot_values_log2_filtT,
                                                                          pseudot=brca_pseudotimes,
                                                                          covarlab=ifelse(mycovar == 1, cond1, cond2),
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
  if(!ens_gs %in% df_beta$gene) {
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
                          exprdt=prot_values_log2_filtT,
                          pseudot=brca_pseudotimes,
                          covarlab=ifelse(mycovar == 1, cond1, cond2),
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
  if(!ens_gs %in% df_beta$gene) {
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
                          exprdt=prot_values_log2_filtT,
                          pseudot=brca_pseudotimes,
                          covarlab=ifelse(mycovar == 1, cond1, cond2),
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



all_traj <- setNames(brca_pseudotimes, rownames(prot_values_log2_filtT))




############## further look at the signif. interactions ##############

# extract the interaction parameters 
# -  feature =>  The feature (usually gene)
# - covariate => The covariate, specified from the order originally supplied to the call to phenopath
# - interaction_effect_size => The effect size of the interaction (β
# - significant_interaction => Boolean for whether the interaction effect is significantly different from 0
# - chi => The precision of the ARD prior on β [Automatic Relevance Determination]
# - pathway_loading => The pathway loading λ

int_dt <- interactions(brca_phenopath_fit)
# stopifnot(as.character(int_dt$feature) %in% names(ens2genes))
int_dt$featureSymb <- as.character(int_dt$feature)
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
  0.6, 0.15, "Gene upregulated\nCond2 increases upregulation",
  0.6, -0.15, "Gene upregulated\nCond2 decreases upregulation",
  -0.7, 0.15, "Gene downregulated\nCond2 decreases downregulation",
  -0.7, -0.15, "Gene downregulated\nCond2 increases downregulation"
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
                                 exprdt=prot_values_log2_filtT,
                                 pseudot=brca_pseudotimes,
                                 covarlab=ifelse(mycovar == 1, cond1, cond2),
                                 valuedt=int_dt,
                                             valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_", int_dt$featureSymb[i_min], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p <- plot_iGeneExpr_gg2(igene= i_min,
                    exprdt=prot_values_log2_filtT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, cond1, cond2),
                    valuedt=int_dt,
                    valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU))

outFile <- file.path(outFolder, paste0("lowestLambdaGene_", int_dt$featureSymb[i_min], "_expr_along_pseudotime_gg2.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG*0.9, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

i_max <- which.max(int_dt$pathway_loading)
p <- plot_iGeneExpr(igene= i_max,
               exprdt=prot_values_log2_filtT,
               pseudot=brca_pseudotimes,
               covarlab=ifelse(mycovar == 1, cond1, cond2),
               valuedt=int_dt,
               valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestLambdaGene_", int_dt$featureSymb[i_max], "_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr_gg2(igene= i_max,
                    exprdt=prot_values_log2_filtT,
                    pseudot=brca_pseudotimes,
                    covarlab=ifelse(mycovar == 1, cond1, cond2),
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


# stopifnot(dfLambda_topGenes$feature %in% gene_lab_dt$geneID_short)
# stopifnot(dfLambda_topGenes$feature %in% names(ens2genes))
dfLambda_topGenes$featureSymb <- dfLambda_topGenes$feature
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

############## PC dots with color-coded by pseudotime gradient ##############  <<<<<<<<<<<<<<< FIG 6a

myconds <- ifelse(rownames(pca_brca_lowrepr) %in% cond1_samples, cond1, 
                  ifelse(rownames(pca_brca_lowrepr) %in% cond2_samples, cond2, NA))
stopifnot(!is.na(myconds))

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_brca),
            condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("#", cond2, "=",nCond2, "; #", cond1, "=",nCond1))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_pseudotimeGrad_", cond1, "_", cond2, ".", plotType))
ggsave(p, filename = outFile, height=myHeightGG*0.9, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,2), summ_dt=summary(pca_brca),
                 condvect = myconds,
                 colvect=brca_pseudotimes,
                 mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
                 mysubtit = paste0("#", cond2, "=",nCond2, "; #", cond1, "=",nCond1))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_pseudotimeGrad_", cond1, "_", cond2, ".", plotType))
ggsave(p, filename = outFile, height=myHeightGG*0.9, width=myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

p <- pcaplot_gg2(pca_dt=data.frame(pca_brca_lowrepr), pctoplot=c(1,3), summ_dt=summary(pca_brca),
                 condvect = myconds,
            colvect=brca_pseudotimes,
            mytit = paste0("TCGA BRCA notNorm (log2(.+1))"),
            mysubtit = paste0("#", cond2, "=",nCond2, "; #", cond1, "=",nCond1))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_13_pseudotimeGrad_", cond1, "_", cond2, ".", plotType))
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


pp_input_data <- prot_values_log2_filtT

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


cond1_samp <- which(rownames(pp_input_data) %in% cond1_samples)
stopifnot(length(cond1_samp) == nCond1)
stopifnot(rownames(pp_input_data[cond1_samp]) == rownames(pp_input_data[cond1_samples]))
corr_expr_pt_cond1 <- apply(pp_input_data[cond1_samp,], 2, cor, brca_pseudotimes[cond1_samp])
to_keep <- which(!is.na(corr_expr_pt_cond1))  ### there is 2 genes with all 0 values !!! -> cannot compute corr
cdf_cond1 <- data.frame(feature = all_genes[to_keep], 
                        correlation = corr_expr_pt_cond1[to_keep])
stopifnot(!is.na(cdf_cond1))


cond2_samp <- which(rownames(pp_input_data) %in% cond2_samples)
stopifnot(length(cond2_samp) == nCond2)
stopifnot(rownames(pp_input_data[cond2_samp]) == rownames(pp_input_data[cond2_samples]))
corr_expr_pt_cond2 <- apply(pp_input_data[cond2_samp,], 2, cor, brca_pseudotimes[cond2_samp])
cdf_cond2 <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt_cond2)
stopifnot(!is.na(cdf_cond2))

outFile <- file.path(outFolder, paste0("corr_gene_expr_pseudot_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cdf_cond2$correlation), main="corr. gene expr. and pseudotimes", col=2)
lines(density(cdf_cond1$correlation), col=1)
legend("topright", legend=c(cond2, cond1), pch=16, col=c(1,2), bty="n")
mtext(side=3, text=paste0("(# ", cond1, "=", nCond1, "; # ", cond2, "=", nCond2, ")"))
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





############## REACTOME PATHWAY ENRICHMENT analysis  ##############
# *a pathway enrichment analysis using Reactome to discover whether any of the top 20 interacting genes (by β value) *
# not sure if take absolute beta in the article ??
# prot_values_log2_filt <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/prot_values_log2_filt.Rdata"))
# prot_values_log2_filtT <- t(prot_values_log2_filt)
# brca_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_BRCA/brca_phenopath_fit.Rdata"))

gene_names <- colnames(prot_values_log2_filtT)
df_beta <- data.frame(beta = interaction_effects(brca_phenopath_fit),
                      beta_sd = interaction_sds(brca_phenopath_fit),
                      is_sig = significant_interactions(brca_phenopath_fit),
                      gene = gene_names)

###### FIRST WAY TO GET MATCHING IDS

httr::set_config(httr::config(ssl_verifypeer = FALSE))  ### added to access ensembl biomart connection

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes_entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                            filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
dim(listOfGenes_entrez)
# [1] 23069     3
 stopifnot(any(df_beta$gene %in% listOfGenes_entrez$hgnc_symbol))
listOfGenes_entrez <- listOfGenes_entrez[listOfGenes_entrez$hgnc_symbol %in% df_beta$gene,]
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
                              keytype="ENTREZID", # for johansson
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

