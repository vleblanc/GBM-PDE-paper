#------------------------------------------------
## Utility functions for scRNA-seq data in SingleCellExperiment objects
#------------------------------------------------

#------------------------------------------------
## Function to subset SingleCellExperiment and scale expression data (if wanted)
#------------------------------------------------

subset.and.scale <- function(full, subset_by = NULL, name = NULL, gene_min_cells = 1, gene_min_umi = 1,
                             scale = FALSE, indiv_scale = FALSE, indiv_scale_by = NULL) {
  if(!is.null(name)) {
    #Subset by relevant colData column
    out <- full[, which(colData(full)[, subset_by] %in% name)]
    colData(out) <- droplevels(colData(out))
    
    #Remove genes that are not expressed
    out <- out[Matrix::rowSums(counts(out) > 0) >= gene_min_cells, ]
    out <- out[Matrix::rowSums(counts(out)) >= gene_min_umi, ]
    
    #Re-normalize within tumour set (centring size factors)
    out <- scater::normalize(out) 
  } else {
    out <- full
  }
  
  #Add scaled expression data
  if(scale) {
    assay(out, "scaled_exprs") <- t(scale(t(as.matrix(logcounts(out))), center = TRUE, scale = TRUE))
  }
  
  #Add scaled expression (individually by tumour)
  if(indiv_scale) {
    assay(out, "indiv_scale") <- as.matrix(do.call(cbind.data.frame, 
                                                   lapply(levels(colData(out)[, indiv_scale_by]), function(sub) {
                                                     exprs_mat <- as.matrix(logcounts(out)[, colData(out)[, indiv_scale_by] == sub])
                                                     scaled_mat <- t(scale(t(exprs_mat), center = TRUE, scale = TRUE))
                                                   })))
  }
  
  return(out)
}




#------------------------------------------------
##Function to run UMAP on SingleCellExperiment assay and return object
#------------------------------------------------

run.umap <- function(object, k = 30, use_red_dim = T, dim_red_use = NULL, sig_pcs = NULL, assay_use = NULL, rows_use = NULL, ...) {
  require(umap)
  
  if(use_red_dim) {
    if(is.null(sig_pcs)) sig_pcs <- attr(object@reducedDims$PCA, "npcs")
    
    mat <- object@reducedDims[[dim_red_use]][, 1:sig_pcs]
  } else {
    if(is.null(rows_use)) rows_use <- rownames(object)
    
    mat <- t(assay(object, assay_use)[rows_use,])
  }
  
  umap_out <- umap(mat, n_neighbors = k, ...)
  
  reducedDim(object, "UMAP") <- umap_out$layout
  metadata(object)$UMAP <- umap_out[names(umap_out) %in% c("knn", "config")]
  
  return(object)
}