### I have tried to reproduce the vignette:
# but I use here glioblastoma data scRNA data from the conquer project
# (Soneson and Robinson (2017); re-quantification of the original dataset) 
# http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE84465.rds

require(phenopath)
require(scater)
require(MultiAssayExperiment)
library(matrixStats)
library(dplyr)

library(ggplot2)


###############################################################
#########################################  PREPARE THE DATA
###############################################################

# We will now parse the data into a form suitable for the scater package, 
# an excellent package for handling single-cell gene expression data. 

# First, read in the MultiAssayExperiment:
mae <- readRDS(file.path("data", "GSE84465.rds"))

# retrieve counts, transcript-per-million (TPM) values and the phenotypic (cell-specific) data 
# and convert it into an SCESet used by scater. We’ll set the default “expression” values to log2(TPM+1)

cts <- assays(experiments(mae)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(mae)[["gene"]])[["TPM"]]
phn <- colData(mae)

################### prepare SCESet and filter the data of interest

sce <- newSCESet(countData = cts, 
                 phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
tpm(sce) <- tpms
exprs(sce) <- log2(tpm(sce) + 1)

# We’re only interested in cells exposed to LPS or PAM, so we parse these from the description column of the SCESet and subset the data accordingly:
  
  is_lps_pam <- grepl("LPS|PAM", sce$description)
sce <- sce[, is_lps_pam]

# Finally, we need to parse the capture time and stimulant from the description column of the SCESet and add them as new columns:
  
  split <- strsplit(as.character(sce$description), "_", fixed = TRUE)
stimulant <- sapply(split, `[`, 1)
time <- sapply(split, `[`, 2)
sce$stimulant <- stimulant
sce$time <- time

# Finally, let’s get MGI symbols for the genes so we actually know what they are:
  
  suppressPackageStartupMessages(library(biomaRt))
ensembl_gene_ids <- sapply(strsplit(featureNames(sce), ".", fixed = TRUE), `[`, 1)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
            filters = "ensembl_gene_id",
            values = ensembl_gene_ids,
            mart = mart)

fData(sce)$mgi_symbol <- rep(NA, nrow(sce))

mm2 <- match(bm$ensembl_gene_id, ensembl_gene_ids)
fData(sce)$mgi_symbol[mm2] <- bm$mgi_symbol


################### Quality control and removal of low-quality cells
  
  sce <- calculateQCMetrics(sce)

#  total number of genes expressed (total_features) against the total number of counts to each cell:
  
  plotPhenoData(sce, aes(x = total_counts, y = total_features))

# We see there are quite a few cells with low counts and features. We’ll remove these via threholds:
  
  sce$to_keep <- sce$total_counts > 5e5 & sce$total_features > 5e3
plotPhenoData(sce, aes(x = total_counts, y = total_features, color = to_keep)) +
  labs(subtitle = "QC pass: total_features > 5000 and total_counts > 50000")

# and subset to the post-qc’d cells:
  
  sce_qc <- sce[, sce$to_keep]

  ################### some additional check
  
  
# In the original publication (Shalek et al. (2014)) 
# the author identified a subset of “cluster-disrupted” cells that were removed.
# These were identified as having low Lyz1 expression and high Serpinb6b expression. 
# Let’s have a look at the co-expression of these two:
  
  Lyz1_index <- grep("Lyz1", fData(sce_qc)$mgi_symbol)
SerpinB6b_index <- grep("SerpinB6b", fData(sce_qc)$mgi_symbol, ignore.case = TRUE)

Lyz1 <- exprs(sce_qc)[Lyz1_index,]
Serpinb6b <- exprs(sce_qc)[SerpinB6b_index,]

qplot(Lyz1, Serpinb6b)

# Accepting cells with Lyz1 expression greater than 0 and Serpbinb6b expression less than 2.5 seems reasonable. 
# Let’s see how this would look:
  
  Serpinb6b_threshold <- 2.5
Lyz1_threshold <- 0

to_keep <- Lyz1 > Lyz1_threshold & Serpinb6b < Serpinb6b_threshold

qplot(Lyz1, Serpinb6b, color = to_keep) +
  geom_vline(xintercept = Lyz1_threshold, linetype = 2) +
  geom_hline(yintercept = Serpinb6b_threshold, linetype = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(subtitle = "Non-cluster-disrupted: Serpinb6b > 2.5 and Lyz1 > 0")

# Let’s now subset the data appropriately:
  
  sce_qc2 <- sce_qc[, to_keep]

  ################### check for batch effects
  
  
# Finally, technical variation can have a large effect on single-cell RNA-seq data. 
# Unfortunately we don’t know the experimental design, but one of the key signs of batch effects 
# is large variation in the number of genes expressed across cells (Hicks et al. (2017)). 
# Let’s see how this affects the principal components of the data:
#   
  plotQC(sce_qc2, type = 'find', var = 'total_features', ntop = 2e3)

# We see this has a huge effect on the overall variation, contributing to the first principal component.
# We can remove this effect using the handy normaliseExprs function in scater:
  
  m <- model.matrix(~ sce_qc2$total_features)
sce_qc2 <- normaliseExprs(sce_qc2, design = m)
exprs(sce_qc2) <- norm_exprs(sce_qc2)

# Let’s tidy up all the SCESets we have lying around before we’re ready for the PhenoPath analysis:
  
  sce <- sce_qc2
rm(sce_qc, sce_qc2)
print(sce)



###############################################################
######################################### Covariate-adjusted pseudotime analysis with PhenoPath
###############################################################

################### Preparing the SCESet for input to PhenoPath
# 
# It’s an open question in the field precisely what genes to use in any pseudotime fit.
# In this work we opt for the most variable genes (in log)
# 
# Here we’ll create a new SCESet that consists of the 7500 most variable genes:

gene_vars <- rowVars(exprs(sce))
var_cutoff <- sort(gene_vars, decreasing = TRUE)[7500]
sce_hvg <- sce[gene_vars >= var_cutoff, ]

# We just have a couple of more things to tidy up before we can fit the model with PhenoPath:
  
#   If any MGI symbol is NA, set it to the corresponding ensembl gene ID
# For some reason the gene Rasgefb1 that’s important to our analysis isn’t annotated, so let’s fix that:
  
  is_na_mgi_symbol <- is.na(fData(sce_hvg)$mgi_symbol)
fData(sce_hvg)$mgi_symbol[is_na_mgi_symbol] <- featureNames(sce_hvg)[is_na_mgi_symbol]

is_Rasgefb1 <- match("ENSMUSG00000029333.14", featureNames(sce_hvg))
fData(sce_hvg)$mgi_symbol[is_Rasgefb1] <- "Rasgefb1"

###################  Inference with PhenoPath
# 
# First we must decide how to pass in the covariate information (ie the stimulant applied) to the software as the x
# values. Here we will give cells exposed to LPS a value of 1 and cells exposed to PAM a value of -1.
# This means the overall pathway loading λ is the average change for LPS and PAM cells, 
# while if the β parameter is positive it means the gene is more upregulated over pseudotime under LPS and 
# if β is negative it means the gene is more upregulated under PAM11 
# Instead we could encode LPS to 1 and PAM to 0, in which case the pathway loading λ would be the change under PAM and λ+β
# the change under LPS stimulation..
# In R we construct this via

x <- 2 * (sce_hvg$stimulant == "LPS") - 1

# By default PhenoPath initialises to the first principal component of the data. 
# However, variational inference is non-convex and we can easily end up in a local maximum in which pseudotime essentially 
# runs backwards in time.
# Simply for convenience sake, we’ll initialise the latent space with the first principal component “flipped” so that the pseudotimes
# will run forwards in time22 Since all pseudotime trajectories are essentially equivalent up to a parity transformation, 
# this won’t affect any of the benchmarking with existing software.:
# ????
  
  scale_vec <- function(x) (x - mean(x)) / sd(x)
pc1 <- prcomp(t(exprs(sce_hvg)), scale = TRUE)$x[,1]
pc1 <- scale_vec(pc1)
time_numeric <- as.numeric(gsub("h", "", sce$time))
pc1 <- pc1 * sign(cor(pc1, time_numeric))

# let's fit phenopath
  
fit <- phenopath(sce_hvg, x, z_init = pc1)


###################  Monitoring convergence
# By default the model is considered converged when the change in the ELBO falls below 10−5
# %. The user should plot the elbo using the plot_elbo command to ensure the lower bound as sufficiently converged:
  
  plot_elbo(fit)

###############################################################
######################################### Interpreting the results
  ###############################################################
  ################### Pseudotime recapitulates capture time

# First, let’s check that the pseudotimes roughly correspond to the true capture times. 
# Note that we can get maximum a-posteriori (MAP) estimates of the pseudotimes using the trajectory function.


## Warning: package 'dplyr' was built under R version 3.4.1
zdf <- data_frame(z = trajectory(fit), time = sce$time)

ggplot(zdf, aes(x = time, y = z, fill = time)) + 
  geom_violin(alpha = 0.8) +
  theme(legend.position = "none") + 
  scale_fill_brewer(palette = "Set1") +
  xlab("Capture time") +
  ylab("Pathway score\n(pseudotime)") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) 

5.2 Identifying significant interactions

We can extract the interaction parameters using the interactions function. This returns a data_frame with the following entries:
  
  feature The feature (usually gene)
covariate The covariate, specified from the order originally supplied to the call to phenopath
interaction_effect_size The effect size of the interaction (β
                                                            
                                                            from the statistical model)
significant Boolean for whether the interaction effect is significantly different from 0
chi The precision of the ARD prior on β
pathway_loading The pathway loading λ

, showing the overall effect for each gene marginalised over the covariate

Let’s have a look at this for our dataset. We’ll swap out the ensembl gene IDs with the MGI symbols to make things a little easier to read:
  
  ints <- interactions(fit)
ints$feature <- fData(sce_hvg)$mgi_symbol
ints[,1:3]

## # A tibble: 7,500 x 3
##    feature   covariate interaction_effect_size
##      <chr>       <chr>                   <dbl>
##  1   Gnai3 covariate_1             -0.21224978
##  2    Narf covariate_1              0.18136389
##  3    Cav2 covariate_1             -0.40507337
##  4   Scmh1 covariate_1              0.07543741
##  5    Xpo6 covariate_1             -0.28933790
##  6    Tfe3 covariate_1             -0.32925474
##  7   Gna12 covariate_1             -0.11623157
##  8    Dlat covariate_1             -0.21062996
##  9    Sdhd covariate_1             -0.23311038
## 10   Ccnd2 covariate_1              0.38030757
## # ... with 7,490 more rows

ints[,4:6]

## # A tibble: 7,500 x 3
##    significant_interaction       chi pathway_loading
##                      <lgl>     <dbl>           <dbl>
##  1                   FALSE 10.121751      -0.1130153
##  2                   FALSE 10.184667      -0.8729432
##  3                   FALSE  9.570431       0.1637015
##  4                   FALSE 10.319062       0.6004500
##  5                   FALSE  9.940135       0.9765032
##  6                   FALSE  9.823192      -0.1777573
##  7                   FALSE 10.286568      -1.1952394
##  8                   FALSE 10.132019      -1.1192479
##  9                   FALSE 10.085082      -1.0163325
## 10                   FALSE  9.702828       1.7546201
## # ... with 7,490 more rows

A nice way to visualise this is to plot the posterior ARD variances (1/χ
) against the posterior interaction effect sizes (β
                                                  
), colouring them by which are found to be significant and annotating the top few genes:
  
  library(ggrepel)
chi_cutoff <- sort(ints$chi)[10]

ggplot(ints, aes(x = interaction_effect_size, y = 1 / chi, 
                 color = significant_interaction)) +
  geom_point() +
  geom_text_repel(data = dplyr::filter(ints, chi < chi_cutoff), 
                  aes(label = feature)) +
  scale_colour_brewer(palette = "Set1")

We can also plot the “landscape” of interactions, where we plot the interaction effect size against the pathway score. The further up the y
-axis a gene is, the more it is upregulated under LPS rather than PAM (and vice-versa), while the further along the x

-axis a gene is, the more it is upregulated over pseudotime regardless of stimulant applied.

ggplot(ints, aes(x = pathway_loading, y = interaction_effect_size, 
                 color = significant_interaction)) +
  geom_point() +
  geom_text_repel(data = dplyr::filter(ints, chi < chi_cutoff), 
                  aes(label = feature), size = 5) +
  scale_colour_brewer(palette = "Set1")  +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)) 

The Tnf gene has the largest interaction effect size - let’s plot it over pseudotime coloured by the stimulant applied:
  
  tnf_index <- grep("^Tnf$", fData(sce_hvg)$mgi_symbol)
sce_hvg$phenopath_pseudotime <- trajectory(fit)

plotExpression(sce_hvg, 
               features = tnf_index, 
               x = "phenopath_pseudotime",
               colour_by = "stimulant",
               show_violin = FALSE)

We see that under PAM it’s upregulated, while under LPS it’s downregulated.
6 Technical

sessionInfo()

## R version 3.4.0 (2017-04-21)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] bindrcpp_0.2                ggrepel_0.6.5              
##  [3] dplyr_0.7.4                 phenopath_0.99.4           
##  [5] matrixStats_0.52.2          biomaRt_2.33.3             
##  [7] scater_1.5.0                ggplot2_2.2.1              
##  [9] Biobase_2.37.2              BiocGenerics_0.23.1        
## [11] MultiAssayExperiment_1.3.35 BiocStyle_2.5.39           
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.4.0              edgeR_3.19.1              
##  [3] tidyr_0.7.1                viridisLite_0.2.0         
##  [5] shiny_1.0.3                assertthat_0.2.0          
##  [7] stats4_3.4.0               vipor_0.4.5               
##  [9] GenomeInfoDbData_0.99.1    yaml_2.1.14               
## [11] progress_1.1.2             RSQLite_1.1-2             
## [13] backports_1.1.0            lattice_0.20-35           
## [15] glue_1.1.1                 limma_3.33.9              
## [17] digest_0.6.12              GenomicRanges_1.29.12     
## [19] XVector_0.17.0             colorspace_1.3-2          
## [21] htmltools_0.3.6            httpuv_1.3.3              
## [23] Matrix_1.2-10              plyr_1.8.4                
## [25] XML_3.98-1.9               pkgconfig_2.0.1           
## [27] bookdown_0.4               zlibbioc_1.23.0           
## [29] purrr_0.2.3                xtable_1.8-2              
## [31] scales_0.4.1               tibble_1.3.4              
## [33] IRanges_2.11.12            SummarizedExperiment_1.7.5
## [35] lazyeval_0.2.0             magrittr_1.5              
## [37] mime_0.5                   memoise_1.1.0             
## [39] evaluate_0.10              beeswarm_0.2.3            
## [41] shinydashboard_0.6.1       tools_3.4.0               
## [43] data.table_1.10.4          prettyunits_1.0.2         
## [45] stringr_1.2.0              S4Vectors_0.15.5          
## [47] munsell_0.4.3              locfit_1.5-9.1            
## [49] DelayedArray_0.3.19        AnnotationDbi_1.39.1      
## [51] compiler_3.4.0             GenomeInfoDb_1.13.4       
## [53] rlang_0.1.2                rhdf5_2.21.2              
## [55] grid_3.4.0                 RCurl_1.95-4.8            
## [57] tximport_1.5.0             rjson_0.2.15              
## [59] bitops_1.0-6               rmarkdown_1.6             
## [61] gtable_0.2.0               DBI_0.7                   
## [63] reshape2_1.4.2             R6_2.2.2                  
## [65] gridExtra_2.2.1            knitr_1.16                
## [67] bindr_0.1                  rprojroot_1.2             
## [69] stringi_1.1.5              ggbeeswarm_0.5.3          
## [71] Rcpp_0.12.13

saveRDS(sce_hvg, "~/Desktop/delete-me.rds")

References

Hicks, Stephanie C, F William Townes, Mingxiang Teng, and Rafael A Irizarry. 2017. “Missing Data and Technical Variability in Single-Cell Rna-Sequencing Experiments.” BioRxiv. Cold Spring Harbor Labs Journals, 025528.

Patro, Rob, Geet Duggal, Michael I Love, Rafael A Irizarry, and Carl Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression.” Nature Methods 14 (4). Nature Research: 417–19.

Shalek, Alex K, Rahul Satija, Joe Shuga, John J Trombetta, Dave Gennert, Diana Lu, Peilin Chen, et al. 2014. “Single-Cell RNA-seq Reveals Dynamic Paracrine Control of Cellular Variation.” Nature 510 (7505): 363–69.

Soneson, Charlotte, and Mark D Robinson. 2017. “Bias, Robustness and Scalability in Differential Expression Analysis of Single-Cell Rna-Seq Data.” BioRxiv. Cold Spring Harbor Labs Journals, 143289.
