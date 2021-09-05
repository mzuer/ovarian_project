
# Run CAVI for CLVM and save result

# Format:
# Rscript clvm [input sceset] [output file] [pc to initialise] 

library(clvm)
library(scater)

scale_vec <- function(x) (x - mean(x)) / sd(x)

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[1]
output_file <- args[2]

pc_initialise <- as.numeric(args[3])

sce <- readRDS(input_rds)

y <- scale(t(exprs(sce)))



ER <- scale_vec(1 * (sce$ER_status == "positive"))
PR <- scale_vec(1 * (sce$PR_status == "positive"))

HER2 <- sce$IHC_HER2
is_equiv <- HER2 == "Equivocal" | HER2 == "Indeterminate" | is.na(HER2)

h2 <- HER2
h2[is.na(h2)] <- "NA"
h2 <- scale_vec(1 * (h2 == "Positive"))

x <- cbind(ER, PR, h2)

pcavi <- clvm(y, x)

saveRDS(pcavi, file = output_file)


