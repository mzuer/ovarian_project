

# for rnaseq qqnorm
# madnorm for microarray

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

quantNorm <- function(x) {qqnorm(x, plot.it=F)$x}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

madNorm <- function(x) {
  rna_rnaseqDT_tmp <- x
  # median center
  all_meds <- apply(rna_rnaseqDT_tmp, 1, median)
  rna_rnaseqDT_tmp <- sweep(rna_rnaseqDT_tmp, 1, all_meds, "-")  # substract from each row their corresponding median
  # mad norm
  # In order to use the MAD as a consistent estimator for the estimation of the standard deviation σ, one takes
  # σ ^ = k ⋅ MAD where k is a constant scale factor, which depends on the distribution.[1]
  # For normally distributed data k is taken to be: ~1.4826
  all_mads <-1.4826 * apply(abs(rna_rnaseqDT_tmp), 1, median)
  rna_madnorm_rnaseqDT <- sweep(rna_rnaseqDT_tmp, 1, all_mads, "/") 
  stopifnot( all ( dim(x) == dim(rna_madnorm_rnaseqDT)))
  rownames(rna_madnorm_rnaseqDT) <- rownames(x)
  colnames(rna_madnorm_rnaseqDT) <- colnames(x)
  return(rna_madnorm_rnaseqDT)
}
