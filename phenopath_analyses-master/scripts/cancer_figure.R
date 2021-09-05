library(ggplot2)
library(cowplot)
library(splines)


brca_limma_plot <- readRDS("figs/brca/limma_plot.rds")
brca_goplot <- readRDS("figs/brca/goplot.rds")
brca_gene_plot <- readRDS("figs/brca/gene_plot.rds")
brca_lplot <- readRDS("figs/brca/lplot.rds")
brca_crossplot <- readRDS("figs/brca/cross_plot.rds")

## Reduce font size of cross plot
brca_crossplot <- brca_crossplot + 
  theme(legend.position = "none", strip.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 11))
  

lsize <- 11
top_right <- plot_grid(NULL,
                       plot_grid(NULL, brca_limma_plot, NULL, nrow = 1, rel_widths = c(1,6, 1)), 
                       brca_goplot, NULL, rel_heights = c(0.3, 4, 4, 1),
                       ncol = 1, labels = c("", "", "C", ""),
                       label_size = lsize)

top_grid <- plot_grid(brca_lplot, top_right, ncol = 2, 
                      labels = c("A", "B"), label_size = lsize,
                      rel_widths = c(2,1))

bottom_grid <- plot_grid(brca_gene_plot, brca_crossplot, 
                       rel_widths = c(3,2), nrow = 1,
                       labels = c("D", "E"), label_size = lsize)


p <- plot_grid(top_grid, bottom_grid, 
               rel_heights = c(3,1.8), label_size = lsize, ncol = 1)

ggsave("figs/brca_figure.png", p, width = 11, height = 10)



# Extra figure for thesis -------------------------------------------------

crossover_thesis <- readRDS("figs/brca/crossover_thesis.rds")

crossover_thesis <- crossover_thesis + 
  theme(legend.position = "none", strip.text.x = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 11))

plot_grid(plot_grid(NULL, brca_crossplot, NULL, rel_widths = c(1,6,1),
                    labels = c("", "A", ""), nrow = 1, label_size = lsize),
          crossover_thesis, labels = c("", "B"),
          ncol = 1, label_size = lsize, rel_heights = c(1,3))

ggsave("figs/brca/crossover_thesis.png", width = 6, height = 10)


# Reformatted COAD plot ---------------------------------------------------


coad_limma_plot <- readRDS("figs/coad/limma_plot.rds")
coad_goplot <- readRDS("figs/coad/goplot.rds")
coad_gene_plot <- readRDS("figs/coad/gene_plot.rds")
coad_lplot <- readRDS("figs/coad/lplot.rds")
coad_tregs <- readRDS("figs/coad/tregs.rds")

upper_grid <- plot_grid(coad_lplot, coad_tregs, 
                        labels = c("A", "B"), nrow = 1, label_size = lsize,
                        rel_widths = c(3,1))

bottom_grid <- plot_grid(plot_grid(NULL, coad_goplot, NULL, rel_heights = c(1,6,1), ncol = 1), 
                         coad_gene_plot, coad_limma_plot, nrow = 1,
                         label_size = lsize, labels = c("C", "D", "E"), 
                         rel_widths = c(3,2,2))

pc <- plot_grid(upper_grid, bottom_grid, ncol = 1, rel_heights = c(5,2.5))

ggsave("figs/coad_figure.png", pc, width = 11, height = 9)
