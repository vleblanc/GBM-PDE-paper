# ---------------------------------------------------------------
## Functions to make scWGS figures
## Author: Veronique LeBlanc
# ---------------------------------------------------------------


#------------------------------------------------
##Function to plot PCA/tSNE/UMAP with custom colours
#------------------------------------------------

plot.red.dim <- function(object, dim_use, exprs_use, colour_by = NULL, col_by_exp = FALSE, col_vals = NULL, shape_by = NULL, size_by = NULL, 
                         plot_title = NULL, leg_title = NULL) {
  require(hues)
  require(ggplot2)
  require(RColorBrewer)
  
  if(dim_use == "PCA") {
    #proj <- attr(object@reducedDims$PCA, "projection")
    proj <- object@reducedDims$PCA[,1:2]
    colnames(proj) <- c("X1", "X2")
  } else {
    proj <- object@reducedDims[[dim_use]]
  }
  
  if(col_by_exp) {
    df_to_plot <- data.frame(proj, colData(object), 
                             row.names = colnames(object))
    
    if(colour_by %in% rownames(assay(object, exprs_use))) {
      df_to_plot <- cbind(df_to_plot, assay(object, exprs_use)[colour_by,])
      colnames(df_to_plot)[ncol(df_to_plot)] <- colour_by
      leg_title <- paste(colour_by, "\nexpression")
    }
    
  } else {
    df_to_plot <- data.frame(proj, colData(object), 
                             row.names = colnames(object))
  }
  
  if(!col_by_exp & !is.factor(df_to_plot[, colour_by])) {
    if(is.numeric(df_to_plot[, colour_by])) {
      col_by_exp <- T
    } else {
      df_to_plot[, colour_by] <- as.factor(df_to_plot[, colour_by])
    }
  }
  
  x_lab <- ifelse(dim_use == "PCA", "PC1",
                  ifelse(dim_use == "tSNE", "tSNE dimension 1",
                         ifelse(dim_use == "UMAP", "UMAP dimension 1", "X1")))
  
  y_lab <- ifelse(dim_use == "PCA", "PC2",
                  ifelse(dim_use == "tSNE", "tSNE dimension 2",
                         ifelse(dim_use == "UMAP", "UMAP dimension 2", "X2")))
  
  p <- ggplot(df_to_plot, aes_string(x = "X1", y = "X2", colour = colour_by, shape = shape_by)) +
    geom_point(alpha = 0.6, stroke = 0) +
    labs(title = plot_title, x = x_lab, y = y_lab) +
    theme_bw()
  #theme(line = element_blank(),
  #      rect = element_blank(),
  #      axis.text = element_blank())
  
  if(col_by_exp) {
    p <- p + scale_colour_gradientn(colors = c("gray", "navy"), guide = guide_colourbar(title = leg_title))
  } else {
    if(is.null(col_vals)) col_vals <- setNames(iwanthue(length(levels(df_to_plot[, colour_by]))), levels(df_to_plot[, colour_by]))
    
    p <- p + scale_colour_manual(values = col_vals, guide = guide_legend(title = leg_title))
  }
  
  return(p)
}


#------------------------------------------------
##Function to plot median copy number across clusters
#------------------------------------------------

plot.clustmeds <- function(clustmeds, object, plot_type = c("dot", "heatmap")) {
  require(ggplot2)
  require(pheatmap)
  
  if(plot_type == "dot") {
    #Melt
    clustmeds_melt <- melt(as.matrix(clustmeds), varnames = c("bin", "dbclust"), value.name = "med_cn")
    
    #Add bin rank with correct order for plotting
    clustmeds_melt$bin <- factor(clustmeds_melt$bin, 
                                 levels = rownames(object)[rownames(object) %in% clustmeds_melt$bin])
    clustmeds_melt$bin_rank <- as.numeric(clustmeds_melt$bin)
    
    #Add chromosome info
    clustmeds_melt$chr <- factor(sapply(strsplit(as.character(clustmeds_melt$bin), "_"), "[", 1),
                                 levels = paste0("chr", c(1:22, "X", "Y")))
    
    #Add copy number status
    clustmeds_melt$status <- ifelse(clustmeds_melt$med_cn == 2, "neutral",
                                    ifelse(clustmeds_melt$med_cn > 2, "gain", "loss"))
    
    ggplot(clustmeds_melt, aes(x = bin_rank, y = med_cn)) +
      geom_point(aes(colour = status), size = 0.25, show.legend = FALSE) +
      scale_colour_manual(values = c("gain" = "red", "loss" = "blue", "neutral" = "gray45")) +
      geom_vline(data = data.frame(pos = cumsum(table(clustmeds_melt$chr[clustmeds_melt$dbclust == levels(clustmeds_melt$dbclust)[1]]))), 
                 aes(xintercept = pos), colour = "grey", size = 0.1) +
      labs(x = NULL, y = "copy number") +
      #coord_cartesian(ylim = c(0, 10)) +
      scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank()) +
      facet_wrap(~ dbclust, ncol = 1)
    
  } else {
    
    breaks = c(0, 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 7.99, 9.99, 11.99, 13.99, max(clustmeds, na.rm = T))
    color = c("#2F6FAB", "#92AFCF", "#F0F0F0", "#F8D6BD", "#F5BE92", "#F2A56C", "#F0884C", "#EB3B25", "#B2271B", "#761E11", "#3E1408")
    
    column_annot <- data.frame(chr = factor(sapply(strsplit(rownames(clustmeds), "_"), "[", 1), 
                                            levels = paste0("chr", c(1:22, "X", "Y"))),
                               row.names = rownames(clustmeds))
    
    pheatmap(t(clustmeds), show_rownames = TRUE, show_colnames = FALSE,
             cluster_cols = FALSE, cluster_rows = FALSE,
             #annotation_row = as.data.frame(colData(object)[,cell_annots, drop = F]),
             annotation_col = column_annot,
             annotation_colors = list(chr = setNames(rep(c("grey30", "grey80"), 12), paste0("chr", c(1:22, "X", "Y")))),
             gaps_col = cumsum(table(column_annot$chr)),
             breaks = breaks, color = color, 
             na_col = "grey70", border_color = NA)
  }
}



#------------------------------------------------
##Function to plot scDNA heatmap (with or without clustering)
#------------------------------------------------
plot.scDNA.cnv <- function(object, assay_use, cell_annots, annot_cols, heat_type = c("cn", "conf", "events"), 
                           cluster_cells = T, tree = NULL, dist_fun = "manhattan", out_file, height = 10, width = 10, ...) {
  require(pheatmap)
  
  cn <- t(assay(object, assay_use))
  
  #Cluster cells based on Manhattan distance and complete linkage (as per 10X)
  if(cluster_cells) {
    if(!is.null(tree)) {
      cells_clust <- tree
    } else {
      cells_dist <- dist(cn, method = dist_fun)
      cells_clust <- hclust(cells_dist)
    }
  } else {
    cells_clust = F
  }
  
  if(heat_type == "cn") {
    breaks = c(0, 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 7.99, 9.99, 11.99, 13.99, max(cn, na.rm = T))
    color = c("#2F6FAB", "#92AFCF", "#F0F0F0", "#F8D6BD", "#F5BE92", "#F2A56C", "#F0884C", "#EB3B25", "#B2271B", "#761E11", "#3E1408")
  }
  
  if(heat_type == "conf") {
    breaks = c(0, seq(0.99, 34.99, by = 1), max(cn[!is.na(cn)]))
    color = colorRampPalette(c("white", "royalblue4"))(36)
  } 
  
  if(heat_type == "events") {
    breaks = NULL
    color = c("white", "black")
  }
  
  pheatmap(cn, show_rownames = F, show_colnames = F,
           cluster_cols = F, cluster_rows = cells_clust,
           annotation_row = as.data.frame(colData(object)[,cell_annots, drop = F]),
           annotation_col = data.frame(chr = factor(sapply(strsplit(colnames(cn), "_"), "[", 1), 
                                                    levels = names(annot_cols$chr)),
                                       row.names = colnames(cn)),
           annotation_colors = annot_cols,
           breaks = breaks, color = color, 
           na_col = "grey70", border_color = NA,
           filename = out_file, height = height, width = width, ...)
}