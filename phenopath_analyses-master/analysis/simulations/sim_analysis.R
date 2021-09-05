library(scater)
library(limma)
library(edgeR)
library(clvm)
library(ggbeeswarm)
library(cowplot)
library(viridis)

kallisto_dirs <- dir("data/simulations/quant", full.names = TRUE)

samples <- sapply(strsplit(kallisto_dirs, "_"), `[`, 2)

sce <- readKallistoResults(directories = kallisto_dirs, samples = samples)

fdata <- read_csv("data/simulations/gene_pars.csv")
pdata <- read_csv("data/simulations/pdata.csv")

sce$sample_index <- as.numeric(sampleNames(sce))
mm <- match(sce$sample_index, pdata$sample)
pdata <- pdata[mm, ]
pData(sce) <- cbind(pData(sce), pdata)

sce <- plotPCA(sce, colour_by = "x", ncomponents = 3, return_SCESet = TRUE)
plotPCA(sce, colour_by = "pst", ncomponents = 3)

fData(sce) <- cbind(fData(sce), fdata)

# limma voom

dge <- DGEList(counts(sce))
dge <- calcNormFactors(dge)

design <- model.matrix(~ (x == 1), pData(sce))
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- decideTests(fit)

qvals <- p.adjust(fit$p.value[,2], method = 'BH')

df_limma <- data_frame(limma_coef = fit$coefficients[,2], 
                       limma_qval = qvals,
                       feature = featureNames(sce))
df_limma <- mutate(df_limma, limma_signif = limma_qval < 0.05)


# PhenoPath

y <- scale(t(exprs(sce)))
x <- cbind(pData(sce)[[ 'x' ]])


pcavi <- clvm(y, x)

pc1 <- prcomp(t(exprs(sce)))$x[,1]

pp_df <- data_frame(pp_alpha = pcavi$m_alpha[1,],
                    pp_beta = pcavi$m_beta[1,],
                    pp_c = pcavi$m_c,
                    pp_signif = significant_interactions(pcavi, n = 2),
                    feature = featureNames(sce),
                    alpha = fData(sce)[['alpha']],
                    beta = fData(sce)[['beta']],
                    c = fData(sce)[['c']],
                    regime = fData(sce)[['regime']])

pp_df <- inner_join(pp_df, df_limma, by = 'feature')

regime_texts <- c("Differential\nexpression only",
                  "Pseudotime\nregulation only",
                  "Pseudotime and\ncovariate interactions",
                  "Differential\nexpression, pseudotime\nand covariate interactions")

pp_df$regime_txt <- plyr::mapvalues(pp_df$regime, 
                                    from = c("de", "pst", "pst_beta", "de_pst_beta"),
                                    to = regime_texts)

pp_df$regime_txt <- factor(pp_df$regime_txt, levels = regime_texts)

pp_pst <- data_frame(phenopath = pcavi[['m_t']],
                     pc1, pst = pData(sce)[['pst']])

limma_plot <- ggplot(pp_df, aes(x = regime_txt, y = abs(limma_coef), color = limma_signif)) + 
  geom_quasirandom(dodge.width = 0.5, alpha = 0.8) +
  xlab("Simulation regime") +
  ylab("Absolute value of limma coefficient") +
  scale_color_brewer(palette = "Set2",
                     name = "Limma voom significant at 5% FDR") +
  theme_cowplot(font_size = 11) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10))

phenopath_plot <- ggplot(pp_df, aes(x = regime_txt, y = abs(pp_beta), color = pp_signif)) + 
  geom_quasirandom(dodge.width = 0.5, alpha = 0.8) +
  xlab("Simulation regime") +
  ylab(expression(paste("|", beta, "|"))) +
  scale_color_brewer(palette = "Set2",
                     name = expression(paste("PhenoPath significant at 2", sigma))) +
  theme_cowplot(font_size = 11) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size = 8))

pc1_plot <- ggplot(pp_pst, aes(x = pst, y = pc1, color = factor(x))) + 
  geom_point(alpha = 0.8) +
  scale_colour_brewer(palette = "Set1", name = "x") +
  theme_cowplot(font_size = 11) +
  xlab("True pseudotime") + ylab("PC1")

z_plot <- ggplot(pp_pst, aes(x = pst, y = phenopath, color = factor(x))) + 
  geom_point(alpha = 0.8) +
  scale_colour_brewer(palette = "Set1", name = "x") +
  theme_cowplot(font_size = 11) +
  xlab("True pseudotime") + ylab("PhenoPath z") +
  theme(legend.position = "none")

comp_plot <- plot_grid(pc1_plot, z_plot, rel_widths = c(6,5))

pca_df <- as_data_frame(redDim(sce)[,1:2]) %>% 
  mutate(x = factor(pData(sce)[['x']]), pst = pp_pst[['pst']])

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = pst, shape = x)) +
  geom_point(size = 2.5) + scale_color_viridis(name = "True\npseudotime") +
  theme_cowplot(font_size = 11)


left_grid <- plot_grid(pca_plot, comp_plot, ncol = 1)
right_grid <- plot_grid(limma_plot, phenopath_plot, ncol = 1)

plot_grid(left_grid, right_grid, nrow = 1)
ggsave("figs/simulations_all.png", width = 11, height = 8)

plot_grid(limma_plot, phenopath_plot, nrow = 1)

ggsave("figs/simulations.png", width = 11.1, height = 3)
