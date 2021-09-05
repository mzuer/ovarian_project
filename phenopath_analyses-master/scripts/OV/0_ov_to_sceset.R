library(scater)
library(tidyverse)
library(magrittr)


library(RTCGA.clinical)
data("OV.clinical")

library(RTCGA.mutations)
data("OV.mutations")

## Read in Kallisto quantified data
ov <- read_tsv("data/OV/TCGA_OV_tpm.tsv.gz")
sample_names <- names(ov)[-1]
ov <- data.frame(ov)
names(ov)[1] <- "feature_id"
rownames(ov) <- ov$feature_id; ov$feature_id <- NULL
names(ov) <- sample_names

ov_counts <- read_tsv("data/OV/TCGA_OV_counts.tsv.gz")
sample_names <- names(ov_counts)[-1]
ov_counts <- data.frame(ov_counts)
names(ov_counts)[1] <- "feature_id"
rownames(ov_counts) <- ov_counts$feature_id; ov_counts$feature_id <- NULL
names(ov_counts) <- sample_names

stopifnot(all.equal(dim(ov), dim(ov_counts)))
stopifnot(all.equal(rownames(ov), rownames(ov_counts)))
stopifnot(all.equal(colnames(ov), colnames(ov_counts)))

## Find feature names (ensembl transcript id + hgnc symbol) and gene types
id_split <- strsplit(rownames(ov), "|", fixed = TRUE)
feature_names <- sapply(id_split, function(x) paste0(x[1], "_", x[6]))
gene_type <- sapply(id_split, `[`, 8)
ensembl_gene_id <- sapply(id_split, `[`, 2)
rownames(ov) <- feature_names
rownames(ov_counts) <- feature_names

## Match CGHubAnalysisID with comparible ID in OV.clinical

id_map <- read_csv("data/TCGA_ID_MAP.csv") %>% filter(Disease == "OV")

to_tcga_barcode <- function(x) tolower(paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "-"))
id_map %<>% mutate(patient_barcode = sapply(AliquotBarcode, to_tcga_barcode))

## Get clinical data associated with IDs
ov_clinical <- OV.clinical %>% 
  filter(patient.bcr_patient_barcode %in% id_map$patient_barcode) %>% 
  select(patient_barcode = patient.bcr_patient_barcode,
         patient.days_to_death,
         patient.age_at_initial_pathologic_diagnosis, 
         patient.days_to_last_followup, patient.days_to_tumor_progression,
         patient.days_to_tumor_recurrence,
         patient.stage_event.clinical_stage,
         patient.tumor_stage) %>% 
  mutate(censored = is.na(patient.days_to_death))

ov_clinical %<>% inner_join(select(id_map, patient_barcode, CGHubAnalysisID, AliquotBarcode),
                            by = c("patient_barcode"))

## Tidy up classes
for(covariate_index in grep("days", names(ov_clinical)))
  ov_clinical[[covariate_index]] <- as.numeric(ov_clinical[[covariate_index]])

plate <- factor(sapply(strsplit(ov_clinical$AliquotBarcode, "-"), `[`, 6))
centre <- factor(sapply(strsplit(ov_clinical$AliquotBarcode, "-"), `[`, 7))
ov_clinical$plate <- factor(plate)
ov_clinical$centre <- factor(centre)

common_cghub_ids <- intersect(colnames(ov), ov_clinical$CGHubAnalysisID)

## Put clinical data in order corresponding to common_cghub_ids
ov_clinical <- ov_clinical[match(common_cghub_ids, ov_clinical$CGHubAnalysisID), ]
ov <- ov[, match(common_cghub_ids, names(ov))]
ov_counts <- ov_counts[, match(common_cghub_ids, names(ov_counts))]
stopifnot(all.equal(colnames(ov), ov_clinical$CGHubAnalysisID))

ov_clinical_pd <- data.frame(ov_clinical)
rownames(ov_clinical_pd) <- ov_clinical_pd$CGHubAnalysisID

sce <- newSCESet(countData = as.matrix(ov_counts),
                 phenoData = AnnotatedDataFrame(ov_clinical_pd))
tpm(sce) <- as.matrix(ov)
exprs(sce) <- log2(as.matrix(ov) + 1)

# sce <- sce[, sce$centre != 31] # get rid of the crappy centre

fData(sce)$gene_type <- gene_type
fData(sce)$ensembl_gene_id <- ensembl_gene_id



# Somatic Mutations -------------------------------------------------------
to_patient_barcode <- function(bcr_barcode) {
  tolower(paste0(strsplit(bcr_barcode, "-")[[1]][1:3], collapse = "-"))
}

patient_mutations <- OV.mutations %>% 
  filter(bcr_patient_barcode != "bcr_patient_barcode") %>% # I mean come on
  group_by(bcr_patient_barcode) %>% 
  summarise(n_mutations = n(), n_somatic = sum(Mutation_Status == "Somatic"),
            n_unknown = sum(Mutation_Status == "Unknown")) %>% 
  mutate(patient_barcode = sapply(bcr_patient_barcode, to_patient_barcode))

pdata_df <- select(pData(sce), patient_barcode)
pdata_df <- left_join(pdata_df, patient_mutations, by = "patient_barcode")

stopifnot(all.equal(sce$patient_barcode, pdata_df$patient_barcode))

pData(sce) <- cbind(pData(sce), select(pdata_df, n_mutations, n_somatic, n_unknown))

saveRDS(sce, file = "data/OV/sce_ov.rds")
