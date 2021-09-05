library(ggplot2)
library(cowplot)

coad <- readRDS("data/COAD/coad_pca_plot.rds")
brca <- readRDS("data/BRCA/brca_pca_plot.rds")


plot_grid(coad, brca, labels = 'AUTO', ncol = 1)

ggsave("figs/supplementary_technical.png", width = 8, height = 10)
