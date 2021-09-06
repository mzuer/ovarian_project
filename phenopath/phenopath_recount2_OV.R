
# Rscript phenopath_recount2_OV.R

require(recount)
require(TCGAbiolinks)
require(biomaRt)
require(phenopath)
require(ggplot2)
require(ggrepel)
require(ggsci)
require(viridis)
require(dplyr)
require(goseq)
genome <- "hg19"
id <- "ensGene"
go_cat <- "GO:BP"
plotNgos <- 10
topCorrThresh <- 0.5  # 0.5 in Campbell and Yau 2018; to select most corr. for GO enrichment
topBetaThresh <- 0.5  # 0.5 in Campbell and Yau 2018; to select highest beta for GO enrichment

library("limma")
library("edgeR")
meanExprThresh_limma <- 0.5


parse_go <- function(go, type, n_tested, plotntop=12) {
  go <- go %>% 
    mutate(qval = p.adjust(over_represented_pvalue)) %>% 
    mutate(log10qval = -log10(qval),
           prop_in_cat = numInCat / n_tested) %>% 
    head(n = plotntop) %>% 
    mutate(type = type) %>% 
    tbl_df() %>% 
    arrange(desc(log10qval))
  go
}

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

pcaplot_gg <- function(pca_dt, pctoplot, summ_dt,colvect, collab="pseudotime") {
  stopifnot(length(pctoplot) == 2)
  # colnames(pca_dt) <- paste0("col", 1:ncol(pca_dt))
  pca_dt$colvect <- colvect
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  p <- ggplot(pca_dt, aes_string(x = paste0("PC", pctoplot[1]),
                            y = paste0("PC", pctoplot[2]),
                            col = paste0("colvect")))+
    xlab(paste0("PC", pctoplot[1], " (", var1, " % variance explained)"))+
  ylab(paste0("PC", pctoplot[2], " (", var2, " % variance explained)"))+
    geom_point() +
    ggtitle(paste0(""), subtitle=paste0(""))+
    labs(color=paste0(collab))+
    scale_color_viridis() +
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  
       return(p)
}
# pcaplot_gg(pca_dt=data.frame(pca_ov_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_ov),colvect=ov_pseudotimes)


pcaplot_gg2 <- function(pca_dt, pctoplot, summ_dt,colvect, condvect, mytit="", mysubtit="", collab="pseudotime") {
  stopifnot(length(pctoplot) == 2)
  # colnames(pca_dt) <- paste0("col", 1:ncol(pca_dt))
  pca_dt$colvect <- colvect
  pca_dt_nocond <- pca_dt
  pca_dt$condvect <- condvect
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  p <- ggplot(pca_dt, aes_string(x = paste0("PC", pctoplot[1]),
                                 y = paste0("PC", pctoplot[2])))+
    geom_point(data = pca_dt_nocond, fill = "grey80", color = "grey80", size = 3) +
    geom_point(aes(fill = colvect), shape = 21, color = 'grey20', size = 3) +
    scale_fill_viridis(name = "Pseudotime", option = "C") + 
    xlab(paste0("PC", pctoplot[1], " (", var1, " % variance explained)"))+
    ylab(paste0("PC", pctoplot[2], " (", var2, " % variance explained)"))+
    ggtitle(paste0(mytit), subtitle=paste0(mysubtit))+
    labs(color=paste0(collab))+
    facet_wrap(~ condvect) +
    theme(legend.position = c(.4,.08),
          legend.direction = "horizontal") +
    theme(strip.background = element_rect(fill = "grey95", color = "grey30", size = 1),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10))
  return(p)
}



myG_theme <-   theme( 
  plot.title = element_text(hjust = 0.5, face = "bold", size=16),
  plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
  panel.grid = element_blank(),
  panel.grid.major.y = element_line(colour = "grey"),
  panel.grid.minor.y = element_line(colour = "grey"),
  axis.line.x= element_line(size = .2, color = "black"),
  axis.line.y = element_line(size = .2, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
  axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
  # axis.ticks.x = element_blank(),
  axis.title.y = element_text(color="black", size=14),
  axis.title.x = element_text(color="black", size=14),
  panel.border = element_blank(),
  panel.background = element_rect(fill = "transparent"),
  legend.background =  element_rect(),
  legend.text = element_text(size=12),
  legend.key = element_blank(),
  legend.title = element_text(face="bold", size=12)
)
plot_iGeneExpr <- function(igene, exprdt, pseudot, covarlab, valuedt, 
                           valuecol="beta", symbcol="geneSymb", subtit="") {
  stopifnot(is.numeric(igene))
  stopifnot(nrow(exprdt) == length(pseudot))
  stopifnot(length(covarlab) == length(pseudot))
  gene_lab <- valuedt[,paste0(symbcol)][igene]
  
  plot_dt <- data.frame(
    y = exprdt[, igene],
    x = covarlab,
    z = pseudot
  )
  p <- ggplot(plot_dt, aes(x = z, y = y, color = x))+
    geom_point() +
    ylab("Expression (GCnorm + log2(.+1))") +
    xlab("PP Pseudotimes")+
    ggtitle(paste0("",gene_lab ), subtitle=paste0(subtit," (", round(valuedt[,paste0(valuecol)][igene], 3), ")"))+
    labs(color="")+
    scale_color_brewer(palette = "Set1") +
    stat_smooth()+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  return(p)
}

betaU <-"\u03B2" 
lambdaU <- "\u03BB"
chiU <- "\u03C7"

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
myHeightGG <- 6
myWidthGG <- 6

# filter to keep only the highly variable genes
# Campbell 2018: 
# for BRAC, whose variance in log(TPM+1) expression was greater than 1 and whose median absolute deviation was greater than 0
# for COAD: s whose median absolute deviation in log(TPM+1) expression was greater than sqrt(0.5)
mad_thresh <- sqrt(0.5)


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
# xx=rownames(ov_data_raw)
# xx=xx[!grepl("_PAR_Y$", xx)]
# stopifnot(!duplicated(gsub("\\..*", "",xx)))
## I have some weird ENSGxx_PAR_Y genes -> remove them
tokeep <- ! grepl("_PAR_Y$", rownames(ov_data_raw))
cat(paste0("... remove weird genes: ", sum(!tokeep), "/", length(tokeep), "\n"))
ov_data_raw <- ov_data_raw[tokeep,]
stopifnot(!duplicated(gsub("\\..*", "",rownames(ov_data_raw))))
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

ov_data_notNorm_log_no0 <- ov_data_notNorm_log[rowSums(ov_data_notNorm_log) > 0,]


pca_ov <- prcomp(t(ov_data_notNorm_log_no0), scale=TRUE)
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
gtex_raw_data_log <- gtex_raw_data_log[rowSums(gtex_raw_data_log) > 0,]
tcga_raw_data_log <- tcga_raw_data_log[rowSums(tcga_raw_data_log) > 0,]


# do I have batch effect within TCGA (plate)
# The default is FALSE for consistency with S, but in general scaling is advisable
# also scaled in https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
tcga_pca_ov <- prcomp(t(tcga_raw_data_log), scale=TRUE)
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
gtex_pca_ov <- prcomp(t(gtex_raw_data_log), scale=TRUE)
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
tcga_data_log <- tcga_data_log[rowSums(tcga_data_log) > 0,]
gtex_data_log <- gtex_data_log[rowSums(gtex_data_log) > 0,]

# do I have batch effect within TCGA (plate)
# The default is FALSE for consistency with S, but in general scaling is advisable
# also scaled in https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
tcga_pca_ov <- prcomp(t(tcga_data_log), scale=TRUE)
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
gtex_pca_ov <- prcomp(t(gtex_data_log), scale=TRUE)
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

gtex_data_log <- gtex_data_log[rowSums(gtex_data_log) > 0,]

########################### REDO THE GTEX PCA WITHOUT THE OUTLIER
# do I have batch effect within gtex
gtex_pca_ov <- prcomp(t(gtex_data_log), scale=TRUE)
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
# ov_data_filtered <- get(load("PHENOPATH_RECOUNT2_TCGA_OV/ov_data_filtered.Rdata"))

# phenopath tuto: sim$y is the N×G matrix of gene expression (N=100 cells and G=40 genes)
# so I have to take the transpose of my ov_data_log_matrix to have samples in line and genes in columns
ov_data_filteredT <- t(ov_data_filtered)
stopifnot(nGenes == ncol(ov_data_filteredT))
stopifnot(nGTEX+nTCGA == nrow(ov_data_filteredT))

pca_ov <- prcomp(ov_data_filteredT, scale=TRUE)
pca_ov_lowrepr <- pca_ov$x
stopifnot(nrow(pca_ov_lowrepr) == nGTEX+nTCGA)

pc1 <- pca_ov_lowrepr[,1]
pc2 <- pca_ov_lowrepr[,2]
pc3 <- pca_ov_lowrepr[,3]

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
stopifnot(nrow(ov_data_filteredT) == nGTEX+nTCGA)
stopifnot(ncol(ov_data_filteredT) == nGenes)
mycovar <- 2 * grepl("^TCGA", rownames(ov_data_filteredT)) - 1
stopifnot(sum(mycovar== -1) == nGTEX)
stopifnot(sum(mycovar== 1) == nTCGA)
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
#   ov_phenopath_fit <- get(load("PHENOPATH_RECOUNT2_TCGA_OV/ov_phenopath_fit.Rdata"))


######################################################## 
##############  PHENOPATH RESULT ANLAYSIS
######################################################## 

####### check the model
# it is important to check convergence with a call to plot_elbo(ov_phenopath_fit) to ensure the ELBO is approximately flat:

outFile <- file.path(outFolder, paste0("phenopath_elbo_fit.", plotType))
p <- plot_elbo(ov_phenopath_fit)
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

  
####### retrieve the pseudotimes
# maximum a-posteriori (MAP) estimates of the pseudotimes using the trajectory()

ov_pseudotimes <-  trajectory(ov_phenopath_fit)
 
############## correlation with PCs ##############
### look if the pseudotimes correlate with the PCs
# tumor=1, normal=-1

# dev.off()
add_corr <- function(x,y, legPos="topright") {
  corval <- as.numeric(round(cor.test(x,y, method="spearman")$estimate, 3))
  legend(legPos, legend=c("TCGA", "GTEX",
    paste0("SCC = ", corval)),
    pch=c(16, 16, -1),
    col =c(1+3, -1+3, "black"), bty="n")
}
outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc1, 
     y = ov_pseudotimes,
     xlab = "PC1", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc1, ov_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc2.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc2, 
     y = ov_pseudotimes,
     xlab = "PC2", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc2, ov_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("corr_pseudotimes_pc3.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = pc3, 
     y = ov_pseudotimes,
     xlab = "PC3", ylab="PP Pseudotimes z",
     col = mycovar+3,
     pch=16, cex=0.7,
     cex.lab = 1.2, cex.axis=1.2)
add_corr(pc3, ov_pseudotimes)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

############## LOOK A BIT AT THE BETA VALUES ##############
# density of the beta with color code signifcance

gene_names <- colnames(ov_data_filteredT)

df_beta <- data.frame(beta = interaction_effects(ov_phenopath_fit),
                      beta_sd = interaction_sds(ov_phenopath_fit),
                      is_sig = significant_interactions(ov_phenopath_fit),
                      gene = gene_names)

outFile <- file.path(outFolder, paste0("beta_values_distribution.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(df_beta$beta))
lines(density(df_beta$beta[df_beta$is_sig]), col=2)
# lines(density(df_beta$beta[!df_beta$is_sig]), col=3)
legend("topright", legend=c("all dist.", "only signif. dist."), col=c(1,2), bty="n")
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
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))
  
  
############## top significant interactions pseudotime x covar ##############
# which gene has highest interaction effect ?

stopifnot(df_beta$gene %in% names(ens2genes))
stopifnot(!duplicated(df_beta$gene))
# if duplicated I would need to do some manual curation...
df_beta$geneSymb <- ens2genes[df_beta$gene]

# the same for the lowest interaction effect ?
  
p <- plot_iGeneExpr(igene= which.max(df_beta$beta), 
               exprdt=ov_data_filteredT,
               pseudot=ov_pseudotimes, 
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"), 
               valuedt=df_beta, 
               valuecol="beta", symbcol="geneSymb", subtit=paste0("highest ", betaU)) 

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr(igene= which.min(df_beta$beta), 
               exprdt=ov_data_filteredT,
               pseudot=ov_pseudotimes, 
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"), 
               valuedt=df_beta, 
                           valuecol="beta", symbcol="geneSymb", subtit=paste0("lowest ", betaU))

outFile <- file.path(outFolder, paste0("lowestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





############## look trajectory and phenotypes (tumor stages, age, etc.) ##############

table(tcga_annot_dt$cgc_case_days_to_death)
table(tcga_annot_dt$cgc_slide_percent_tumor_cells)


stopifnot(rownames(ov_data_filteredT) == c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id))


all_traj <- setNames(ov_pseudotimes, rownames(ov_data_filteredT))
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

plot_pheno_catego <- function(annot_dt, plotvar, plotlab, varords=NULL) {
  tcga_traj_dt <- data.frame(
    tcga_samp = annot_dt$cgc_sample_id,
    tcga_plotvar = annot_dt[,plotvar],
    pseudotime = all_traj[annot_dt$cgc_sample_id],
    stringsAsFactors = FALSE
  )
  na_txt <- paste0(sum(!is.na(tcga_traj_dt$tcga_plotvar)),
  "/", length(tcga_traj_dt$tcga_plotvar))
  tcga_traj_dt <- tcga_traj_dt[!is.na(tcga_traj_dt$tcga_plotvar),]
  stopifnot(!is.na(tcga_traj_dt))
  if(!is.null(varords)){
    tcga_traj_dt$tcga_plotvar <- factor(tcga_traj_dt$tcga_plotvar, levels=varords) 
  }
  stopifnot(!is.na(tcga_traj_dt$tcga_plotvar))
  
  p <- ggplot(tcga_traj_dt, aes(x= tcga_plotvar, y= pseudotime) )+
    geom_boxplot(notch = F, outlier.shape=NA)+
    geom_jitter(aes(col=tcga_plotvar),alpha=0.7,position=position_jitterdodge())+
    ggtitle(paste0("Pseudotime by ", plotlab), 
            subtitle = paste0("(TCGA data; av.: ", na_txt, ")"))+
    scale_color_nejm()+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    ylab("PP Pseudotimes")+
    xlab(paste0("TCGA HGSC  ", plotlab))+
    myG_theme +
    labs(color="")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank() )
  return(p)
}
myplotlab <- "tumor stage"
stg_ords <-c(  "Stage IC", "Stage IIA" ,"Stage IIB" , "Stage IIC","Stage IIIA", "Stage IIIB","Stage IIIC",   "Stage IV")
plot_pheno_catego(tcga_annot_dt, plotvar= "cgc_case_clinical_stage", plotlab=myplotlab, varords=stg_ords)

outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# plot_pheno_catego(tcga_annot_dt, plotvar= "cgc_drug_therapy_drug_name", plotlab="drug therapy", varords=NULL)
## too many categories 

plot_pheno_continuous <- function(annot_dt, plotvar, plotlab) {
  tcga_traj_dt <- data.frame(
    tcga_samp = annot_dt$cgc_sample_id,
    tcga_plotvar = annot_dt[,plotvar],
    pseudotime = all_traj[annot_dt$cgc_sample_id],
    stringsAsFactors = FALSE
  )
  na_txt <- paste0(sum(!is.na(tcga_traj_dt$tcga_plotvar)),
                   "/", length(tcga_traj_dt$tcga_plotvar))
  tcga_traj_dt <- tcga_traj_dt[!is.na(tcga_traj_dt$tcga_plotvar),]
  stopifnot(!is.na(tcga_traj_dt))

  p <- ggplot(tcga_traj_dt, aes(x= tcga_plotvar, y= pseudotime) )+
    geom_point() +
    ggtitle(paste0("Pseudotime by ", plotlab), 
            subtitle = paste0("(TCGA data; av.: ", na_txt, ")"))+
    ylab("PP Pseudotimes")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    xlab(paste0("TCGA HGSC  ", plotlab))+
    stat_smooth()+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  return(p)
}
myplotlab <- "year of birth"
plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.demographic.year_of_birth", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "days to death"
plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_case_days_to_death", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "diag. days to birth"
plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.diagnoses.days_to_birth", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "age at diag."
plot_pheno_continuous(tcga_annot_dt, plotvar= "gdc_cases.diagnoses.age_at_diagnosis", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct normal cells"
plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_normal_cells", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct neutrophil infilt."
plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_neutrophil_infiltration", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

myplotlab <- "slide pct monocyte infilt."
plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_monocyte_infiltration", plotlab=myplotlab)
outFile <- file.path(outFolder, paste0("boxplot_pseudotimes_by_", paste0(gsub(" ", "_",myplotlab)),"_TCGA.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


myplotlab <- "slide pct lymphocyt infilt."
plot_pheno_continuous(tcga_annot_dt, plotvar= "cgc_slide_percent_lymphocyt_infiltration", plotlab=myplotlab)
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

int_dt <- interactions(ov_phenopath_fit)
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
                                 exprdt=ov_data_filteredT,
                                 pseudot=ov_pseudotimes, 
                                 covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"), 
                                 valuedt=int_dt, 
                                             valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("lowest ", lambdaU)) 

outFile <- file.path(outFolder, paste0("lowestLambdaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


p <- plot_iGeneExpr(igene= which.max(int_dt$pathway_loading), 
               exprdt=ov_data_filteredT,
               pseudot=ov_pseudotimes, 
               covarlab=ifelse(mycovar == 1, "TCGA", "GTEX"), 
               valuedt=int_dt, 
               valuecol="pathway_loading", symbcol="featureSymb", subtit=paste0("highest ", lambdaU))

outFile <- file.path(outFolder, paste0("highestBetaGene_expr_along_pseudotime.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


############## top significant pathway loading ##############

# look at the top bottom and up interaction effects
nTop <- 10

# show the top genes with signif interactions
int_dt <- int_dt[order(int_dt$pathway_loading, decreasing = TRUE),]
topPosGenes <- int_dt[1:nTop,]
topPosGenes$dir <- "pos"
int_dt <- int_dt[order(int_dt$pathway_loading, decreasing = FALSE),]
topNegGenes <- int_dt[1:nTop,]
topNegGenes$dir <- "neg"

dfLambda_topGenes <- rbind(topPosGenes, topNegGenes)
stopifnot(!duplicated(dfLambda_topGenes$feature))
dfLambda_topGenes$feature <- as.character(dfLambda_topGenes$feature)


stopifnot(dfLambda_topGenes$feature %in% gene_lab_dt$featureID_short)
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
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

############## PC dots with color-coded by pseudotime gradient ##############

p <- pcaplot_gg2(pca_dt=data.frame(pca_ov_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_ov),
            condvect = gsub("(^.+?)-.+", "\\1", rownames(pca_ov_lowrepr)),
            colvect=ov_pseudotimes,
            mytit = paste0("TCGA+GTEX OV notNorm (log2(.+1))"),
            mysubtit = paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_23_pseudotimeGrad_tcga_gtex.", plotType))
ggsave(p, filename = outFile, height=myHeight, width=myWidth*1.2)
cat(paste0("... written: ", outFile, "\n"))


p <- pcaplot_gg2(pca_dt=data.frame(pca_ov_lowrepr), pctoplot=c(1,2), summ_dt=summary(pca_ov),
            condvect = gsub("(^.+?)-.+", "\\1", rownames(pca_ov_lowrepr)),
            colvect=ov_pseudotimes,
            mytit = paste0("TCGA+GTEX OV notNorm (log2(.+1))"),
            mysubtit = paste0("nGTEX=",nGTEX, "; nTCGA=",nTCGA))

outFile <- file.path(outFolder, paste0("in_raw_data_pca_12_pseudotimeGrad_tcga_gtex.", plotType))
ggsave(p, filename = outFile, height=myHeight, width=myWidth*1.2)
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


pp_input_data <- ov_data_filteredT

all_genes <- colnames(pp_input_data) # ensemblID
corr_expr_pt <- apply(pp_input_data, 2, cor, ov_pseudotimes)
cdf <- data.frame(feature = all_genes, 
                  correlation = corr_expr_pt)
stopifnot(!is.na(cdf))

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
  ggtitle(paste0("GO enrichment - top ", plotNgos), subtitle = paste0("abs. corr. with Pseudotime > ", topCorrThresh, ")"))+
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
  parse_go(goup, "Up-regulated", length(upreg_genes), 10),
  parse_go(godown, "Down-regulated", length(downreg_genes), 10)
)


betacorr_gos$term <- factor(betacorr_gos$term, levels = betacorr_gos$term[order(betacorr_gos$log10qval)])
ggplot(betacorr_gos, aes(x = term, y = log10qval, color = type)) +
  geom_point() +
  facet_wrap(~ type, scales = "free", nrow = 2) +
  coord_flip() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  ylab(expression(paste(log[10], " q-value"))) +
  theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 11))

# for reactome pathway analysis
# https://github.com/kieranrcampbell/phenopath_revisions/blob/master/analysis/brca_reanalysis/clvm_analysis.Rmd

############## complementary info to std gene DE ##############

# look at genes that are differentially expressed (large difference of gene expression normal vs. tumor)
# and that have different trajectories (high interaction effects)
# Whilst many of these 92 genes are differentially expressed between MSI groups, including MLH2 and TGFBR2 (Fig. 5c), 
# PhenoPath is able to resolve the dynamic contribution to these expression differences


# see e.g. http://research.libd.org/recountWorkflow/articles/recount-workflow.html
# ov_data_raw comes from here:
# ovary_gtex <- scale_counts(TCGAquery_recount2(project="gtex", tissue = "ovary")$gtex_ovary)
# ovary_tcga <- scale_counts(TCGAquery_recount2(project="tcga", tissue = "ovary")$tcga_ovary)
# # I skip some steps of processing, but I have in the end 
# # gene x sample matrix with GTEX samples starting with GTEX- and TCGA samples with TCGA-
# ovary_all_data <- cbind(assays(ovary_gtex)$counts, assays(ovary_tcga)$counts)

ov_data_raw <- get(load(file.path(inFolder, paste0("all_counts_onlyPF_", purityFilter, ".Rdata"))))

# filter because I have removed the outlier !
stopifnot(tcga_annot_dt$cgc_sample_id %in% colnames(ov_data_raw))
stopifnot(gtex_annot_dt$sampid %in% colnames(ov_data_raw))
stopifnot(!duplicated(tcga_annot_dt$cgc_sample_id))
stopifnot(!duplicated(gtex_annot_dt$sampid))
ov_data_raw <- ov_data_raw[, c(gtex_annot_dt$sampid, tcga_annot_dt$cgc_sample_id) ]

# remove the weird gene
tokeep <- ! grepl("_PAR_Y$", rownames(ov_data_raw))
cat(paste0("... remove weird genes: ", sum(!tokeep), "/", length(tokeep), "\n"))
ov_data_raw <- ov_data_raw[tokeep,]
stopifnot(!duplicated(gsub("\\..*", "",rownames(ov_data_raw))))


## Extract counts and filter out lowly expressed geens

to_keep <- rowMeans(ov_data_raw) > meanExprThresh_limma
cat(paste0("... keep sufficient expr: ", sum(to_keep), "/", length(tokeep), "\n"))


## Build DGEList object
dge <- DGEList(counts = ov_data_raw[to_keep, ])
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


stopifnot(!duplicated(rownames(gsub("\\..*", "",rownames(ov_data_raw)))))

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

outFile <- file.path(outFolder, paste0("limma_adjPval_vs_phenopath_beta.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  y=-log10(merged_dt$adj.P.Val),
  x=df_beta$beta,
  pch=16,
  col=2+2*merged_dt$is_sig,
  cex.lab=1.2,
  cex.axis=1.2,
  ylab="limm adj Pval [-log10]",
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
pca_ov_norm <- prcomp(t(ov_data_log_norm), scale=TRUE)
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





df_large <- data.frame(
  y = ov_data_filteredT[, which_lowest],
  x = mycovar,
  z = ov_pseudotimes
)

#dev.off()
df_large$covar_lab <- ifelse(mycovar == 1, "TCGA", "GTEX")

p <- ggplot(df_large, aes(x = z, y = y, color = covar_lab))+
  geom_point() +
  ylab("Expression (GCnorm + log2(.+1))") +
  xlab("PP Pseudotimes")+
  ggtitle(paste0("",beta_topGene ), subtitle=paste0("highest beta (", round(df_beta$beta[which_lowest], 3), ")"))+
  labs(color="")+
  scale_color_brewer(palette = "Set1") +
  stat_smooth()+
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5))


which_largest <- which.max(df_beta$beta)
beta_topGene <- df_beta$geneSymb[which_largest]

df_large <- data.frame(
  y = ov_data_filteredT[, which_largest],
  x = mycovar,
  z = ov_pseudotimes
)

#dev.off()
df_large$covar_lab <- ifelse(mycovar == 1, "TCGA", "GTEX")

p <- ggplot(df_large, aes(x = z, y = y, color = covar_lab))+
  geom_point() +
  ylab("Expression (GCnorm + log2(.+1))") +
  xlab("PP Pseudotimes")+
  ggtitle(paste0("",beta_topGene ), subtitle=paste0("highest beta (", round(df_beta$beta[which_largest], 3), ")"))+
  labs(color="")+
  scale_color_brewer(palette = "Set1") +
  stat_smooth()+
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5))


df_curve <- frame_data(
  ~ ER_status, ~ x, ~ xend, ~ y, ~ yend, ~ curvature,
  "ER-positive", -50, 20, 27, 40, 0.1,
  "ER-negative", 55, 60, -35, 40, -0.1
)
geom_curve(aes(x = x, y = y, xend = xend, yend = yend),
           data = filter(df_curve, ER_status == "ER-negative"), color = 'black',
           curvature = df_curve$curvature[2], arrow = arrow(length = unit(0.3, "cm"), type = "open"), 
           size = 1.5) +
  
  
df_curve$ER_status <- factor(df_curve$ER_status, levels = df_curve$ER_status)
plt <- filter(df_pca, ER_status != "indeterminate") %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(data = df_no_er, fill = "grey80", color = "grey80", size = 3) +
  geom_point(aes(fill = pseudotime), shape = 21, color = 'grey20', size = 3) +
  scale_fill_viridis(name = "Pseudotime", option = "C") + 
  facet_wrap(~ ER_status) +
  theme(legend.position = c(.4,.08),
        legend.direction = "horizontal") +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend),
             data = filter(df_curve, ER_status == "ER-negative"), color = 'black',
             curvature = df_curve$curvature[2], arrow = arrow(length = unit(0.3, "cm"), type = "open"), 
             size = 1.5) +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend),
             data = filter(df_curve, ER_status == "ER-positive"), color = 'black',
             curvature = df_curve$curvature[1], arrow = arrow(length = unit(0.3, "cm"), type = "open"), 
             size = 1.5) +
  theme(strip.background = element_rect(fill = "grey95", color = "grey30", size = 1),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  xlim(-72, 83) + ylim(-50, 50)
```



Time for some reactome fun:
  
  ```{r reactome-time}
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
  
  ```{r graph-for-metrics, eval = FALSE}
sce <- readRDS("../../data/BRCA/sce_brca_gene_level.rds")
all_genes <- unique(unlist(gene_list_ensembl))
all_genes <- all_genes[all_genes %in% fData(sce)$ensembl_gene_id]
mm <- match(all_genes, fData(sce)$ensembl_gene_id)
is_basal <- sce$PAM50_mRNA == "Basal-like"
Y <- t(exprs(sce)[mm, ])
colnames(Y) <- all_genes
df_gex <- Y %>% 
  as_data_frame() %>% 
  dplyr::mutate(is_basal, z = pData(sce)[['tmap']]) %>% 
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



