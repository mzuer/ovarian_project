library(scater)
library(tidyverse)
library(magrittr)

sanitize_names <- function(n) {
  n <- gsub("(", "", n, fixed = TRUE)
  n <- gsub(")", "", n, fixed =  TRUE)
  n <- gsub(" ", "_", n, fixed = TRUE)
  n <- gsub("-", "_", n, fixed  =TRUE)
  n
}


library(RTCGA.clinical)
data("BRCA.clinical")

library(RTCGA.mutations)
data("BRCA.mutations")

## Read in Kallisto quantified data
brca <- read_tsv("data/BRCA/TCGA_BRCA_tpm.tsv.gz")
sample_names <- names(brca)[-1]
brca <- data.frame(brca)
names(brca)[1] <- "feature_id"
rownames(brca) <- brca$feature_id; brca$feature_id <- NULL
names(brca) <- sample_names

brca_counts <- read_tsv("data/BRCA/TCGA_BRCA_counts.tsv.gz")
sample_names <- names(brca_counts)[-1]
brca_counts <- data.frame(brca_counts)
names(brca_counts)[1] <- "feature_id"
rownames(brca_counts) <- brca_counts$feature_id; brca_counts$feature_id <- NULL
names(brca_counts) <- sample_names

stopifnot(all.equal(dim(brca), dim(brca_counts)))
stopifnot(all.equal(rownames(brca), rownames(brca_counts)))
stopifnot(all.equal(colnames(brca), colnames(brca_counts)))

## Find feature names (ensembl transcript id + hgnc symbol) and gene types
id_split <- strsplit(rownames(brca), "|", fixed = TRUE)
feature_names <- sapply(id_split, function(x) paste0(x[1], "_", x[6]))
gene_type <- sapply(id_split, `[`, 8)
ensembl_gene_id <- sapply(id_split, `[`, 2)
rownames(brca) <- feature_names
rownames(brca_counts) <- feature_names

## Match CGHubAnalysisID with comparible ID in brca.clinical

id_map <- read_csv("data/TCGA_ID_MAP.csv") %>% filter(Disease == "BRCA")

to_tcga_barcode <- function(x) tolower(paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "-"))
id_map %<>% mutate(patient_barcode = sapply(AliquotBarcode, to_tcga_barcode))

## Get clinical data associated with IDs
brca_clinical <- BRCA.clinical %>% 
  filter(patient.bcr_patient_barcode %in% id_map$patient_barcode) %>% 
  select(patient_barcode = patient.bcr_patient_barcode,
         patient.days_to_death,
         patient.age_at_initial_pathologic_diagnosis, 
         patient.days_to_last_followup, 
         ER_status = patient.breast_carcinoma_estrogen_receptor_status,
         PR_status = patient.breast_carcinoma_progesterone_receptor_status) %>% 
  mutate(censored = is.na(patient.days_to_death))

## Add in HER2 and survival info
brca_clin_from_cbio <- read_tsv("data/BRCA/brca_tcga_clinical_data.tsv")
names(brca_clin_from_cbio) <- sapply(names(brca_clin_from_cbio), sanitize_names)
brca_clin_from_cbio <- select(brca_clin_from_cbio,
                              Disease_Free_Months, IHC_HER2, Overall_Survival_Months,
                              Patient_ID) %>% 
  mutate(patient_barcode = tolower(Patient_ID))

## Strip out duplicates
mm <- match(unique(brca_clin_from_cbio$patient_barcode), brca_clin_from_cbio$patient_barcode)
brca_clin_from_cbio <- brca_clin_from_cbio[mm, ]

brca_clinical %<>% inner_join(brca_clin_from_cbio, by = "patient_barcode")

brca_clinical %<>% inner_join(select(id_map, patient_barcode, CGHubAnalysisID, AliquotBarcode),
                            by = c("patient_barcode"))

## Tidy up classes
for(covariate_index in grep("days", names(brca_clinical)))
  brca_clinical[[covariate_index]] <- as.numeric(brca_clinical[[covariate_index]])

plate <- factor(sapply(strsplit(brca_clinical$AliquotBarcode, "-"), `[`, 6))
centre <- factor(sapply(strsplit(brca_clinical$AliquotBarcode, "-"), `[`, 7))
brca_clinical$plate <- factor(plate)
brca_clinical$centre <- factor(centre)

common_cghub_ids <- intersect(colnames(brca), brca_clinical$CGHubAnalysisID)

## Put clinical data in order corresponding to common_cghub_ids
brca_clinical <- brca_clinical[match(common_cghub_ids, brca_clinical$CGHubAnalysisID), ]
brca <- brca[, match(common_cghub_ids, names(brca))]
brca_counts <- brca_counts[, match(common_cghub_ids, names(brca_counts))]
stopifnot(all.equal(colnames(brca), brca_clinical$CGHubAnalysisID))

brca_clinical_pd <- data.frame(brca_clinical)
rownames(brca_clinical_pd) <- brca_clinical_pd$CGHubAnalysisID

sce <- newSCESet(countData = as.matrix(brca_counts),
                 phenoData = AnnotatedDataFrame(brca_clinical_pd))
exprs(sce) <- log2(as.matrix(brca) + 1)
tpm(sce) <- as.matrix(brca)

rm(brca_counts, brca, brca_clinical_pd, brca_clinical)
gc()

fData(sce)$gene_type <- gene_type
fData(sce)$ensembl_gene_id <- ensembl_gene_id



# Somatic Mutations -------------------------------------------------------
to_patient_barcode <- function(bcr_barcode) {
  tolower(paste0(strsplit(bcr_barcode, "-")[[1]][1:3], collapse = "-"))
}

patient_mutations <- BRCA.mutations %>%
  filter(bcr_patient_barcode != "bcr_patient_barcode") %>% # I mean come on
  group_by(bcr_patient_barcode) %>%
  summarise(n_mutations = n(), n_somatic = sum(Mutation_Status == "Somatic"),
            n_unknown = sum(Mutation_Status == "Unknown")) %>%
  mutate(patient_barcode = sapply(bcr_patient_barcode, to_patient_barcode))

patient_mutations <- filter(patient_mutations, 
                            patient_barcode %in% sce$patient_barcode)

## end

pdata_df <- dplyr::select(pData(sce), patient_barcode)
pdata_df <- left_join(pdata_df, patient_mutations, by = "patient_barcode")

mm <- match(sce$patient_barcode, pdata_df$patient_barcode)
pdata_df <- pdata_df[mm, ]

stopifnot(all.equal(sce$patient_barcode, pdata_df$patient_barcode))

pData(sce) <- cbind(pData(sce), dplyr::select(pdata_df, n_mutations, n_somatic, n_unknown))


# Add in phenotypic data --------------------------------------------------

# library(readxl)
# xl <- read_excel("data/BRCA/supplementary.xls", skip = 1)
# xl$patient_barcode <- tolower(xl$`Complete TCGA ID`)
# names(xl) <- sapply(names(xl), function(n) gsub(" ", "_", n))
# xl <- dplyr::select(xl, patient_barcode, HER2_Final_Status,
#                     Tumor, Node, Metastasis, AJCC_Stage, PAM50_mRNA, 
#                     RPPA_Clusters, Days_to_date_of_Death)

# pdata_df <- pData(sce)
# pd <- left_join(pdata_df, xl, by = "patient_barcode")

# stopifnot(all.equal(sce$patient_barcode, pd$patient_barcode))
# pData(sce) <- pd

saveRDS(sce, file = "data/BRCA/sce_brca.rds")
