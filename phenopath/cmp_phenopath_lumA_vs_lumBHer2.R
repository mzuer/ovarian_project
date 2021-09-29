# look at the intersect between top beta lumA-her2 and lumb-her2

vs_lumB_pp_ <-  get(load("PHENOPATH_RECOUNT2_TCGA_BRCA_LUMALUMB/brca_phenopath_fit.Rdata"))
vs_lumB_dt <-  get(load("PHENOPATH_RECOUNT2_TCGA_BRCA_LUMALUMB/brca_data_filtered.Rdata"))
vs_lumB_pp <- trajectory(vs_lumB_pp_)
gene_names_lumB <- rownames(vs_lumB_dt)
tmp_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2/gene_dt.Rdata"))
tmp_dt$geneID <- gsub("(.+)\\..+", "\\1", tmp_dt$geneID)
tmp_dt <- tmp_dt[tmp_dt$geneID %in% gene_names_lumB,]
tmp_dt <- unique(tmp_dt)
stopifnot(!duplicated(tmp_dt$geneID))
ens2gene_lumB <- setNames(tmp_dt$geneSymb, tmp_dt$geneID)

vs_her2_pp_ <-  get(load("PHENOPATH_RECOUNT2_TCGA_BRCA_LUMAHER2/brca_phenopath_fit.Rdata"))
vs_her2_dt <-  get(load("PHENOPATH_RECOUNT2_TCGA_BRCA_LUMAHER2/brca_data_filtered.Rdata"))
vs_her2_pp <- trajectory(vs_her2_pp_)
gene_names_her2 <- rownames(vs_her2_dt)
tmp_dt <- get(load("../tcga_data/DOWNLOAD_TCGA_BRCA_LUMAHER2_RECOUNT2/gene_dt.Rdata"))
tmp_dt$geneID <- gsub("(.+)\\..+", "\\1", tmp_dt$geneID)
tmp_dt <- tmp_dt[tmp_dt$geneID %in% gene_names_her2,]
tmp_dt <- unique(tmp_dt)
stopifnot(!duplicated(tmp_dt$geneID))
ens2gene_her2 <- setNames(tmp_dt$geneSymb, tmp_dt$geneID)

stopifnot(gene_names_lumB %in% names(ens2gene_lumB))
df_beta_lumB <- data.frame(beta = interaction_effects(vs_lumB_pp_),
                      beta_sd = interaction_sds(vs_lumB_pp_),
                      is_sig = significant_interactions(vs_lumB_pp_),
                      gene = gene_names_lumB,
                      geneSymb = ens2gene_lumB[paste0(gene_names_lumB)],
                      stringsAsFactors = FALSE)

stopifnot(gene_names_her2 %in% names(ens2gene_her2))
df_beta_her2 <- data.frame(beta = interaction_effects(vs_her2_pp_),
                           beta_sd = interaction_sds(vs_her2_pp_),
                           is_sig = significant_interactions(vs_her2_pp_),
                           gene = gene_names_her2,
                           geneSymb = ens2gene_her2[paste0(gene_names_her2)],
                           stringsAsFactors = FALSE)

nTop <- 50
topBeta_up_lumB <- df_beta_lumB[order(df_beta_lumB$beta, decreasing = TRUE), "geneSymb"][1:nTop]
topBeta_up_her2 <- df_beta_her2[order(df_beta_lumB$beta, decreasing = TRUE), "geneSymb"][1:nTop]

topBeta_down_lumB <- df_beta_lumB[order(df_beta_lumB$beta, decreasing = F), "geneSymb"][1:nTop]
topBeta_down_her2 <- df_beta_her2[order(df_beta_lumB$beta, decreasing = F), "geneSymb"][1:nTop]

intersect(topBeta_up_lumB, topBeta_up_her2)
intersect(topBeta_down_lumB, topBeta_down_her2)

