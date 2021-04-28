# ---------------------------------------------------------------
## Scripts to get and plot variant allele frequencies of germline mutations (identified in WES data) in merged clone scWGS data (Supplementary Figure 2)
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_load.R")

JK153_cnv_germ_cov <- list(b = load.germ.cov.vcf("./data/cellranger_dna/JK153_mergebams/clone_b.germline.vcf.gz"),
                           c = load.germ.cov.vcf("./data/cellranger_dna/JK153_mergebams/clone_c.germline.vcf.gz"),
                           d = load.germ.cov.vcf("./data/cellranger_dna/JK153_mergebams/clone_d.germline.vcf.gz"),
                           non_mal = load.germ.cov.vcf("./data/cellranger_dna/JK153_mergebams/clone_non-malignant.germline.vcf.gz"))

JK153_cnv_germ_cov_plots <- lapply(names(JK153_cnv_germ_cov), function(clone) {
  mat <- JK153_cnv_germ_cov[[clone]]
  
  #Fix chromosomes
  mat$CHROM <- factor(mat$CHROM, levels = paste0("chr", c(1:22, "X", "Y")))
  
  #Add bin info as a factor
  mat$bin_rank <- 1:nrow(mat)
  
  #Melt allele fractions per bin
  mat_melt <- melt(mat[,c("ref_af", "alt_af", "bin_rank")], id.vars = "bin_rank", value.name = "af")
  
  ggplot(mat_melt, aes(x = bin_rank, y = af)) +
    geom_point(colour = "grey60", size = 0.5, show.legend = FALSE) +
    geom_vline(data = data.frame(pos = cumsum(table(mat$CHROM))), 
               aes(xintercept = pos), colour = "black", size = 0.1) +
    labs(x = NULL, y = "variant alelle fraction", title = clone) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank())
})

cowplot::plot_grid(plotlist = JK153_cnv_germ_cov_plots, ncol = 1)