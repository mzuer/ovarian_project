

betaU <-"\u03B2" 
lambdaU <- "\u03BB"
chiU <- "\u03C7"

add_corr <- function(x,y, legPos="topright", cond_plus1="TCGA", cond_minus1="GTEX") {
  corval <- as.numeric(round(cor.test(x,y, method="spearman")$estimate, 3))
  legend(legPos, legend=c(paste0(cond_plus1), paste0(cond_minus1),
                          paste0("SCC = ", corval)),
         pch=c(16, 16, -1),
         col =c(1+3, -1+3, "black"), bty="n")
}


parse_go <- function(go, type, n_tested, plotntop=12) {
  go <- go %>% 
    mutate(qval = p.adjust(over_represented_pvalue)) %>% 
    mutate(log10qval = -log10(qval),
           prop_in_cat = numInCat / n_tested) %>% 
    head(n = plotntop) %>% 
    mutate(type = type) %>% 
    tbl_df() %>% 
    arrange(desc(log10qval))
  go
}

pcaplot <- function(pca_dt, pctoplot, summ_dt,...) {
  stopifnot(length(pctoplot) == 2)
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  plot(x = pca_dt[, pctoplot[1]],
       y = pca_dt[, pctoplot[2]],
       xlab=paste0("PC", pctoplot[1], " (", var1, " % variance explained)"),
       ylab=paste0("PC", pctoplot[2], " (", var2, " % variance explained)"),
       pch = 16,
       cex=0.7,
       ...
  )
}

pcaplot_gg <- function(pca_dt, pctoplot, summ_dt,colvect, collab="pseudotime") {
  stopifnot(length(pctoplot) == 2)
  # colnames(pca_dt) <- paste0("col", 1:ncol(pca_dt))
  pca_dt$colvect <- colvect
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  p <- ggplot(pca_dt, aes_string(x = paste0("PC", pctoplot[1]),
                                 y = paste0("PC", pctoplot[2]),
                                 col = paste0("colvect")))+
    xlab(paste0("PC", pctoplot[1], " (", var1, " % variance explained)"))+
    ylab(paste0("PC", pctoplot[2], " (", var2, " % variance explained)"))+
    geom_point() +
    ggtitle(paste0(""), subtitle=paste0(""))+
    labs(color=paste0(collab))+
    scale_color_viridis() +
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  
  return(p)
}
# pcaplot_gg(pca_dt=data.frame(pca_ov_lowrepr), pctoplot=c(2,3), summ_dt=summary(pca_ov),colvect=ov_pseudotimes)


pcaplot_gg2 <- function(pca_dt, pctoplot, summ_dt,colvect, condvect, mytit="", mysubtit="", collab="pseudotime") {
  stopifnot(length(pctoplot) == 2)
  # colnames(pca_dt) <- paste0("col", 1:ncol(pca_dt))
  pca_dt$colvect <- colvect
  pca_dt_nocond <- pca_dt
  pca_dt$condvect <- condvect
  var1 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[1])], 2)
  var2 <- round(100*summ_dt$importance["Proportion of Variance", paste0("PC", pctoplot[2])], 2)
  p <- ggplot(pca_dt, aes_string(x = paste0("PC", pctoplot[1]),
                                 y = paste0("PC", pctoplot[2])))+
    geom_point(data = pca_dt_nocond, fill = "grey80", color = "grey80", size = 3) +
    geom_point(aes(fill = colvect), shape = 21, color = 'grey20', size = 3) +
    scale_fill_viridis(name = "Pseudotime", option = "C") + 
    xlab(paste0("PC", pctoplot[1], " (", var1, " % variance explained)"))+
    ylab(paste0("PC", pctoplot[2], " (", var2, " % variance explained)"))+
    ggtitle(paste0(mytit), subtitle=paste0(mysubtit))+
    labs(color=paste0(collab))+
    facet_wrap(~ condvect) +
    theme(legend.position = c(.4,.08),
          legend.direction = "horizontal") +
    theme(strip.background = element_rect(fill = "grey95", color = "grey30", size = 1),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10))
  return(p)
}



myG_theme <-   theme( 
  plot.title = element_text(hjust = 0.5, face = "bold", size=16),
  plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
  panel.grid = element_blank(),
  panel.grid.major.y = element_line(colour = "grey"),
  panel.grid.minor.y = element_line(colour = "grey"),
  axis.line.x= element_line(size = .2, color = "black"),
  axis.line.y = element_line(size = .2, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
  axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
  # axis.ticks.x = element_blank(),
  axis.title.y = element_text(color="black", size=14),
  axis.title.x = element_text(color="black", size=14),
  panel.border = element_blank(),
  panel.background = element_rect(fill = "transparent"),
  legend.background =  element_rect(),
  legend.text = element_text(size=12),
  legend.key = element_blank(),
  legend.title = element_text(face="bold", size=12)
)
plot_iGeneExpr <- function(igene, exprdt, pseudot, covarlab, valuedt, 
                           valuecol="beta", symbcol="geneSymb", subtit="") {
  stopifnot(is.numeric(igene))
  stopifnot(nrow(exprdt) == length(pseudot))
  stopifnot(length(covarlab) == length(pseudot))
  gene_lab <- valuedt[,paste0(symbcol)][igene]
  
  plot_dt <- data.frame(
    y = exprdt[, igene],
    x = covarlab,
    z = pseudot
  )
  p <- ggplot(plot_dt, aes(x = z, y = y, color = x))+
    geom_point() +
    ylab("Expression (GCnorm + log2(.+1))") +
    xlab("PP Pseudotimes")+
    ggtitle(paste0("",gene_lab ), subtitle=paste0(subtit," (", round(valuedt[,paste0(valuecol)][igene], 3), ")"))+
    labs(color="")+
    scale_color_brewer(palette = "Set1") +
    stat_smooth()+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  return(p)
}



plot_iGeneExpr_gg2 <- function(igene, exprdt, pseudot, covarlab, valuedt, 
                               valuecol="beta", symbcol="geneSymb", subtit="") {
  stopifnot(is.numeric(igene))
  stopifnot(nrow(exprdt) == length(pseudot))
  stopifnot(length(covarlab) == length(pseudot))
  gene_lab <- valuedt[,paste0(symbcol)][igene]
  
  plot_dt <- data.frame(
    y = exprdt[, igene],
    x = covarlab,
    z = pseudot
  )
  plot_dt_nocond <- plot_dt
  plot_dt_nocond$x <- NULL
  
  p <- ggplot(plot_dt, aes(x = z, y = y))+
    geom_point(data = plot_dt_nocond, fill = "grey80", color = "grey80", size = 3) +
    geom_point(aes(fill = x), shape = 21, color = 'grey20', size = 3) +
    facet_wrap(~ x) +
    ylab("Expression (GCnorm + log2(.+1))") +
    xlab("PP Pseudotimes")+
    ggtitle(paste0("",gene_lab ), subtitle=paste0(subtit," (", round(valuedt[,paste0(valuecol)][igene], 3), ")"))+
    labs(fill="")+
    scale_fill_brewer(palette = "Set1") +
    stat_smooth()+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  return(p)
}

plot_pheno_catego <- function(annot_dt, pt_traj,plotvar, plotxlab, plotylab="PP Pseudotimes", varords=NULL) {
  tcga_traj_dt <- data.frame(
    tcga_samp = annot_dt$cgc_sample_id,
    tcga_plotvar = annot_dt[,plotvar],
    pseudotime = pt_traj[annot_dt$cgc_sample_id],
    stringsAsFactors = FALSE
  )
  na_txt <- paste0(sum(!is.na(tcga_traj_dt$tcga_plotvar)),
                   "/", length(tcga_traj_dt$tcga_plotvar))
  tcga_traj_dt <- tcga_traj_dt[!is.na(tcga_traj_dt$tcga_plotvar),]
  stopifnot(!is.na(tcga_traj_dt))
  if(!is.null(varords)){
    tcga_traj_dt$tcga_plotvar <- factor(tcga_traj_dt$tcga_plotvar, levels=varords)
  }
  stopifnot(!is.na(tcga_traj_dt$tcga_plotvar))
  
  p <- ggplot(tcga_traj_dt, aes(x= tcga_plotvar, y= pseudotime) )+
    geom_boxplot(notch = F, outlier.shape=NA)+
    geom_jitter(aes(col=tcga_plotvar),alpha=0.7,position=position_jitterdodge())+
    ggtitle(paste0("Pseudotime by ", plotxlab),
            subtitle = paste0("(TCGA data; av.: ", na_txt, ")"))+
    scale_color_nejm()+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    xlab(paste0(plotxlab))+
    ylab(paste0(plotylab))+
    myG_theme +
    labs(color="")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank() )
  return(p)
}

plot_pheno_continuous <- function(annot_dt, pt_traj, plotvar, plotxlab, plotylab="PP Pseudotimes") {
  tcga_traj_dt <- data.frame(
    tcga_samp = annot_dt$cgc_sample_id,
    tcga_plotvar = annot_dt[,plotvar],
    pseudotime = pt_traj[annot_dt$cgc_sample_id],
    stringsAsFactors = FALSE
  )
  na_txt <- paste0(sum(!is.na(tcga_traj_dt$tcga_plotvar)),
                   "/", length(tcga_traj_dt$tcga_plotvar))
  tcga_traj_dt <- tcga_traj_dt[!is.na(tcga_traj_dt$tcga_plotvar),]
  stopifnot(!is.na(tcga_traj_dt))
  
  p <- ggplot(tcga_traj_dt, aes(x= tcga_plotvar, y= pseudotime) )+
    geom_point() +
    ggtitle(paste0("Pseudotime by ", plotxlab),
            subtitle = paste0("(TCGA data; av.: ", na_txt, ")"))+
    ylab(paste0(plotylab))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    xlab(paste0(plotxlab))+
    stat_smooth()+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  return(p)
}
