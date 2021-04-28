# ---------------------------------------------------------------
## Functions for analyses of 10X CNV data
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

#------------------------------------------------
##Function to get coarse-grained CNV info
#------------------------------------------------

get.coarse.cnv <- function(object) {
  require(GenomicRanges)
  require(data.table)
  
  #Get median window size across cells
  window_bins <- median(object$window_size)
  window_length <- window_bins * 20000
  
  #Convert bins to GRanges object
  bins_gr <- makeGRangesFromDataFrame(data.frame(chr = sapply(strsplit(rownames(object), "_"), "[", 1),
                                                 start = sapply(strsplit(rownames(object), "_"), "[", 2),
                                                 end = sapply(strsplit(rownames(object), "_"), "[", 3)))
  
  #Get coarse bins (based on window length)
  cbins_gr <- unlist(slidingWindows(range(bins_gr), width = window_length, step = window_length))
  
  #Get indices for 20kb bins to collapse to form coarse bins
  bin_overlaps <- findOverlaps(query = cbins_gr, subject = bins_gr, ignore.strand = T)
  cbins_indices <- queryHits(bin_overlaps)
  
  #Aggregate copy number information across bins for each coarse bin (much faster if converted to data table)
  cbins_cn_agg <- data.table(cbind(cbins_indices, as.matrix(assay(object, "copy_number"))))
  
  cbins_cn_agg <- cbins_cn_agg[, lapply(.SD, function(cbin_cn) {
    #If more than half of the bins are unmappable, set the whole bin as unmappable - otherwise, take the rounded mean
    if(sum(is.na(cbin_cn)) > (length(cbin_cn) / 2)) {
      agg_cn <- as.integer(NA)
    } else {
      agg_cn <- as.integer(round(mean(cbin_cn, na.rm = T)))
    }
    return(agg_cn)
  }), by = cbins_indices]
  
  cbins_cn_agg <- data.frame(cbins_cn_agg[,2:ncol(cbins_cn_agg)])
  
  
  #Aggregate confidence information across bins for each coarse bin (much faster if converted to data table)
  cbins_conf_agg <- data.table(cbind(cbins_indices, as.matrix(assay(object, "copy_number"))))
  
  cbins_conf_agg <- cbins_conf_agg[, lapply(.SD, function(cbin_cn) {
    #If more than half of the bins are unmappable, set the whole bin as unmappable - otherwise, take the rounded mean
    if(sum(is.na(cbin_cn)) > (length(cbin_cn) / 2)) {
      agg_cn <- as.integer(NA)
    } else {
      agg_cn <- as.integer(round(mean(cbin_cn, na.rm = T)))
    }
    return(agg_cn)
  }), by = cbins_indices]
  
  cbins_conf_agg <- data.frame(cbins_conf_agg[,2:ncol(cbins_conf_agg)])
  
  #Fix column names and add bin ranges as row names
  colnames(cbins_conf_agg) <- colnames(cbins_cn_agg) <- gsub("[.]", "-", colnames(cbins_cn_agg))
  rownames(cbins_conf_agg) <- rownames(cbins_cn_agg) <- paste(as.character(seqnames(cbins_gr)), start(cbins_gr), end(cbins_gr), sep = "_")
  
  #Make new Single Cell Experiment object to return
  object_new <- SingleCellExperiment(list(copy_number = as.matrix(cbins_cn_agg), confidence = as.matrix(cbins_conf_agg)), 
                                     rowData = data.frame(chr = as.character(seqnames(cbins_gr)),
                                                          start_pos = start(cbins_gr),
                                                          end_pos = end(cbins_gr),
                                                          row.names = rownames(cbins_cn_agg)), 
                                     colData = colData(object))
  
  return(object_new)
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


