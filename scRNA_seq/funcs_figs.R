#library(SC3)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


#------------------------------------------------
##Function to plot tSNE with custom colours
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
##Function to plot cluster tree
#------------------------------------------------

clust.tree <- function(object, method = c("pca", "exprs", "nn_seurat", "nn_shekhar"), group = "gclust",
                       dim_red_use = "PCA", weight_vals = NULL, pcs_use = NULL, genes_use = NULL) {
  require(Seurat)
  require(ape)
  
  if(method == "pca") {
    data_use <- data.frame(row.names = c(1:pcs_use))
    
    for (i in levels(colData(object)[, group])) {
      pca_tmp <- object@reducedDims[[dim_red_use]][colData(object)[, group] == i, 1:pcs_use]
      data_tmp <- apply(pca_tmp, 2, mean)
      data_use <- cbind(data_use, data_tmp)
    }
    colnames(data_use) <- paste0("c", levels(colData(object)[, group]))
    
    #eigenvals <- attr(object@reducedDims$PCA, "eigen")[1:pcs_use]
    weights <- weight_vals/sum(weight_vals)
    
    data_dist <- Seurat::CustomDistance(my.mat = data_use, my.function = weighted.euclidean.distance, w = weights)
    
    data_tree <- hclust(data_dist)
    data_tree_phylo <- ape::as.phylo(data_tree)
    return(list(hclust = data_tree, plot = plot.phylo(data_tree_phylo, direction = "downwards")))
  }
  
  
  if(method == "exprs") {
    if(is.null(genes_use)) {
      genes_use <- rownames(object)
    }
    
    avg_exp <- data.frame(row.names = genes_use)
    for (i in levels(object$gclust)) {
      exp_tmp <- exprs(object)[genes_use, object$gclust == i]
      exp_tmp <- apply(exp_tmp, 1, mean)
      avg_exp <- cbind(avg_exp, exp_tmp)
    }
    colnames(avg_exp) <- paste0("c", levels(object$gclust))
    
    data_dist <- dist(t(x = avg_exp))
    
    data_tree <- as.phylo(x = hclust(d = data_dist))
    plot.phylo(x = data_tree, direction = "downwards")
  }
  
  
  if(method == "nn_seurat") {
    data_use <- object@reducedDimension[, 1:pcs_use]
    adj <- get.edges(data_use, nn = 30, do_jaccard = TRUE)
    
    num_clusts <- length(levels(object$gclust))
    
    data_dist <- matrix(data = 0, nrow = num_clusts, ncol = num_clusts)
    rownames(data_dist) <- levels(object$gclust)
    colnames(data_dist) <- levels(object$gclust)
    
    for (i in 1:(num_clusts - 1)) {
      for (j in (i + 1):num_clusts) {
        adj_tmp <- adj[rownames(adj) %in% colnames(object)[object$gclust == i],
                       rownames(adj) %in% colnames(object)[object$gclust == j]]
        d <- mean(adj_tmp)
        
        if (is.na(d)) {
          data_dist[i, j] <- 0
        } else {
          data_dist[i, j] <- d
        }
      }
    }
    diag(data_dist) <- 1
    data.dist <- dist(data.dist)
    
    data_tree <- as.phylo(hclust(data_dist))
    plot.phylo(data.tree, direction = "downwards")
  }
  
  
  if(method == "nn_shekhar") {
    clust_dists <- compute.cluster.distances(object, clust_col = "gclust", reduction_use = "pca", dist_type = "nn", 
                                             pcs_use = 1:pcs_use)
    
    diag(clust_dists) <- 1
    data_dist <- dist(clust_dists)
    
    data_tree <- as.phylo(hclust(data_dist))
    plot.phylo(data_tree, direction = "downwards")
  }
}

weighted.euclidean.distance <- function(x, y, w) {
  v.dist <- sum(sqrt(x = w * (x - y) ^ 2))
  return(v.dist)
}





singleR.draw.heatmap.custom <- function (SingleR, cells.use = NULL, types.use = NULL, clusters = NULL, top.n = 40, 
                                         normalize = T, order.by.clusters = F, cells_order = NULL, silent = F, ...) {
  require(pheatmap)
  
  scores = SingleR$scores
  if (!is.null(cells.use)) {
    scores = scores[cells.use, ]
  }
  if (!is.null(types.use)) {
    scores = scores[, types.use]
  }
  m = apply(t(scale(t(scores))), 2, max)
  thres = sort(m, decreasing = TRUE)[min(top.n, length(m))]
  data = as.matrix(scores)
  if (normalize == T) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data - mmin)/(mmax - mmin)
    data = data^3
  }
  data = data[, m > (thres - 1e-06)]
  data = t(data)
  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = "Clusters"
    rownames(clusters) = colnames(data)
  }
  additional_params = list(...)
  if (is.null(additional_params$ors)) {
    annotation_colors = NA
  }
  else {
    annotation_colors = additional_params$annotation_colors
  }
  clustering_method = "ward.D2"
  if (order.by.clusters == T) {
    data = data[, order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters), , drop = F]
    
    clust_pos <- as.vector(table(clusters))
    for(i in 2:length(clust_pos)) {
      clust_pos[i] <- clust_pos[i-1] + clust_pos[i]
    }
    
    pheatmap(data, border_color = NA, show_colnames = FALSE, 
             clustering_method = clustering_method, fontsize_row = 9, 
             annotation_col = clusters, cluster_cols = F, silent = silent, 
             annotation_colors = annotation_colors, gaps_col = clust_pos, ...)
  }
  else if (!is.null(cells_order)) {
    data = data[, cells_order]
    clusters = clusters[cells_order, , drop = F]
    pheatmap(data, border_color = NA, show_colnames = FALSE, 
             clustering_method = clustering_method, fontsize_row = 9, 
             annotation_col = clusters, cluster_cols = F, silent = silent, 
             annotation_colors = annotation_colors, ...)
  }
  else {
    if (!is.null(clusters)) {
      pheatmap(data, border_color = NA, show_colnames = FALSE, 
               clustering_method = clustering_method, fontsize_row = 9, 
               annotation_col = clusters, silent = silent, annotation_colors = annotation_colors, ...)
    }
    else {
      pheatmap(data[, sample(ncol(data))], border_color = NA, 
               show_colnames = FALSE, clustering_method = clustering_method, 
               fontsize_row = 9, silent = silent, annotation_colors = annotation_colors, ...)
    }
  }
}
