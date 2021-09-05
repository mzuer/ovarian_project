library(scater)
library(tidyverse)
library(magrittr)
library(readxl)


library(RTCGA.clinical)
data("COAD.clinical")

library(RTCGA.mutations)
data("COAD.mutations")

## Read in Kallisto quantified data
## Downloaded from https://osf.io/95wnv/
coad <- read_tsv("data/COAD/TCGA_COAD_tpm.tsv.gz")
sample_names <- names(coad)[-1]
coad <- data.frame(coad)
names(coad)[1] <- "feature_id"
rownames(coad) <- coad$feature_id; coad$feature_id <- NULL
names(coad) <- sample_names

coad_counts <- read_tsv("data/COAD/TCGA_COAD_counts.tsv.gz")
sample_names <- names(coad_counts)[-1]
coad_counts <- data.frame(coad_counts)
names(coad_counts)[1] <- "feature_id"
rownames(coad_counts) <- coad_counts$feature_id; coad_counts$feature_id <- NULL
names(coad_counts) <- sample_names

stopifnot(all.equal(dim(coad), dim(coad_counts)))
stopifnot(all.equal(rownames(coad), rownames(coad_counts)))
stopifnot(all.equal(colnames(coad), colnames(coad_counts)))

## Find feature names (ensembl transcript id + hgnc symbol) and gene types
id_split <- strsplit(rownames(coad), "|", fixed = TRUE)
feature_names <- sapply(id_split, function(x) paste0(x[1], "_", x[6]))
gene_type <- sapply(id_split, `[`, 8)
ensembl_gene_id <- sapply(id_split, `[`, 2)
rownames(coad) <- feature_names
rownames(coad_counts) <- feature_names

## Match CGHubAnalysisID with comparible ID in OV.clinical

id_map <- read_csv("data/TCGA_ID_MAP.csv") %>% filter(Disease == "COAD")

## Parse and add data to id_map
to_tcga_barcode <- function(x) tolower(paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "-"))
id_map %<>% mutate(patient_barcode = sapply(AliquotBarcode, to_tcga_barcode))
sample_split <- sapply(strsplit(id_map$AliquotBarcode, "-", fixed = TRUE), `[`, 4)
sample <- as.numeric(sapply(sample_split, substr, 1, 2))
vial <- sapply(sample_split, substr, 3, 3)

#' If the vial is in 1:10 then it's a tumour and otherwise
#' it's normal, so classify as such:

sample_type <- sapply(sample %in% 1:10, ifelse, "tumour", "normal")

id_map %<>% mutate(sample_type, vial)

## Get clinical data associated with IDs
coad_clinical <- COAD.clinical %>% 
  filter(patient.bcr_patient_barcode %in% id_map$patient_barcode) %>% 
  select(patient_barcode = patient.bcr_patient_barcode,
         patient.days_to_death,
         patient.age_at_initial_pathologic_diagnosis, 
         patient.days_to_last_followup,
         days_to_sample_procurement = patient.biospecimen_cqcf.days_to_sample_procurement,
         patient.samples.sample.days_to_collection,
         stage = patient.stage_event.pathologic_stage,
         nuclei_percent = as.numeric(patient.biospecimen_cqcf.tumor_samples.tumor_sample.tumor_nuclei_percent),
         weight = as.numeric(patient.biospecimen_cqcf.tumor_samples.tumor_sample.tumor_weight),
         neoplasm_subdivision = patient.anatomic_neoplasm_subdivision,
         histological_type = patient.clinical_cqcf.histological_type,
         polyps_present = patient.colon_polyps_present,
         t_stage = patient.stage_event.tnm_categories.pathologic_categories.pathologic_t,
         n_stage = patient.stage_event.tnm_categories.pathologic_categories.pathologic_n,
         m_stage = patient.stage_event.tnm_categories.pathologic_categories.pathologic_m,
         msi_status = patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status) %>% 
  mutate(censored = is.na(patient.days_to_death))

coad_clinical %<>% inner_join(select(id_map, patient_barcode, CGHubAnalysisID, AliquotBarcode, sample_type, vial),
                            by = c("patient_barcode"))

## Tidy up classes
for(covariate_index in grep("days", names(coad_clinical)))
  coad_clinical[[covariate_index]] <- as.numeric(coad_clinical[[covariate_index]])

plate <- factor(sapply(strsplit(coad_clinical$AliquotBarcode, "-"), `[`, 6))
coad_clinical$plate <- factor(plate)

common_cghub_ids <- intersect(colnames(coad), coad_clinical$CGHubAnalysisID)

## Put clinical data in order corresponding to expression data
coad_clinical <- coad_clinical[match(common_cghub_ids, coad_clinical$CGHubAnalysisID), ]
coad <- coad[, match(common_cghub_ids, names(coad))]
coad_counts <- coad_counts[,match(common_cghub_ids, names(coad_counts)) ]
stopifnot(all.equal(colnames(coad), coad_clinical$CGHubAnalysisID))

coad_clinical_pd <- data.frame(coad_clinical)
rownames(coad_clinical_pd) <- coad_clinical_pd$CGHubAnalysisID

sce <- newSCESet(tpmData = as.matrix(coad), countData = coad_counts,
                 phenoData = AnnotatedDataFrame(coad_clinical_pd))

tpm(sce) <- as.matrix(coad)
exprs(sce) <- log2(tpm(sce) + 1)

fData(sce)$gene_type <- gene_type
fData(sce)$ensembl_gene_id <- ensembl_gene_id


# Collapse some of the plates ---------------------------------------------

sce$short_plate <- as.character(sce$plate)
to_collapse <- names(which(table(sce$plate) < 30))
sce$short_plate[which(sce$short_plate %in% to_collapse)] <- "other"

# Somatic Mutations -------------------------------------------------------
to_patient_barcode <- function(bcr_barcode) {
  tolower(paste0(strsplit(bcr_barcode, "-")[[1]][1:3], collapse = "-"))
}

patient_mutations <- COAD.mutations %>% 
  filter(bcr_patient_barcode != "bcr_patient_barcode") %>% # I mean come on
  group_by(bcr_patient_barcode) %>% 
  summarise(n_mutations = n(), n_somatic = sum(Mutation_Status == "Somatic"),
            n_unknown = sum(Mutation_Status == "Unknown")) %>% 
  mutate(patient_barcode = sapply(bcr_patient_barcode, to_patient_barcode))

pdata_df <- select(pData(sce), patient_barcode)
pdata_df <- left_join(pdata_df, patient_mutations, by = "patient_barcode")

stopifnot(all.equal(sce$patient_barcode, pdata_df$patient_barcode))

pData(sce) <- cbind(pData(sce), select(pdata_df, n_mutations, n_somatic, n_unknown))


# Add in MSI data ----------------------------------

msi_data <- read_excel("data/COAD/nm.4191-S3.xlsx")
names(msi_data) <- paste0("nm_", sapply(names(msi_data), function(x) gsub(" ", "_", x)))
msi_data <- filter(msi_data, nm_Tumor_Type == "COAD") %>% 
  mutate(nm_Sample_Name = tolower(nm_Sample_Name)) %>% 
  rename(patient_barcode = nm_Sample_Name)

pd <- pData(sce)

pd2 <- left_join(pd, msi_data, by = "patient_barcode")

stopifnot(all.equal(pd2$patient_barcode, pd$patient_barcode))
pData(sce) <- pd2


# Add new survival data ---------------------------------------------------

coad_survival <- read_tsv("data/COAD/coadread_tcga_clinical_data.tsv")[c("Patient ID", "Overall Survival (Months)", "Overall Survival Status")]
names(coad_survival) <- gsub(" ", "_", names(coad_survival))
names(coad_survival) <- gsub("(", "", names(coad_survival), fixed = TRUE)
names(coad_survival) <- gsub(")", "", names(coad_survival), fixed = TRUE)

coad_survival <- dplyr::rename(coad_survival, patient_barcode = Patient_ID) %>% 
  dplyr::mutate(patient_barcode = tolower(patient_barcode))

## remove duplicated rows
coad_survival <- coad_survival[-which(duplicated(coad_survival$patient_barcode)), ]

pd3 <- pData(sce)
pd3 <- left_join(pd3, coad_survival, by = "patient_barcode")

stopifnot(all.equal(pd3$patient_barcode, pData(sce)$patient_barcode))
pData(sce) <- pd3

saveRDS(sce, file = "data/COAD/sce_coad.rds")
