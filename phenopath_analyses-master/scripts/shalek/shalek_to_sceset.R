
library(tidyverse)
library(scater)


raw <- read_tsv("data/shalek/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz")

genes <- raw[[1]]
raw <- raw[,-1]

cell_names <- colnames(raw)

pam_cells <- grep("^PAM", cell_names)
pic_cells <- grep("^PIC", cell_names)
lps_cells <- grep("^LPS", cell_names)

cells_oi <- c(pam_cells, pic_cells, lps_cells)
raw <- raw[, cells_oi]
cell_names_pd <- names(raw)

names_split <- strsplit(names(raw), "_")
stimulant <- sapply(names_split, `[`, 1)
time <- sapply(names_split, `[`, 2)

pdata <- data.frame(stimulant, time)
rownames(pdata) <- cell_names_pd

tpmdata <- as.matrix(raw)
rownames(tpmdata) <- genes

sce <- newSCESet(tpmData = tpmdata, 
                 phenoData = new("AnnotatedDataFrame", pdata))

is_exprs(sce) <- exprs(sce) > 0

sce <- calculateQCMetrics(sce)

saveRDS(sce, file = "data/shalek/sce_shalek.rds")

