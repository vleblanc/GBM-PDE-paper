##-----------------------------------------------
# Functions to perform clustering and associated analyses
# Author: Veronique LeBlanc
##-----------------------------------------------


##-----------------------------------------------
# Master function to rull full clustering analysis on a SingleCellExperiment object
##-----------------------------------------------
run.full.cluster.analysis <- function(object, set, block = NULL, pca_iters = 1000, sig_pcs = NULL, batch_correct = FALSE, bio_thresh = 0.01, fdr_thresh = 0.05, 
                                      k = 30, de_method = c("t_test", "wilcox"), block_markers = TRUE,
                                      run_pca = TRUE, save_hvgs = TRUE, run_batch_correct = TRUE, run_gclust = TRUE, get_graph = TRUE, get_doublet_scores = TRUE, 
                                      calc_markers = TRUE, run_tSNE = TRUE, run_umap = TRUE, get_clust_scores = TRUE, 
                                      save_objects = TRUE, n_cores, rand_seed = 1234, ...) {
  source("/projects/vleblanc_prj/GBM_organoids/code/analysis/functions_analysis.R")
  require(scater)
  require(scran)
  
  set.seed(rand_seed)
  
  if(run_pca){
    metadata(object)$pca_iters <- pca_iters
    
    if(!is.null(block)) {
      if(!is.factor(block)){
        stop("\'block\' must be a factor")
      }
      
      metadata(object)$block <- block
      
      #-Get batch variances
      block_vars <- setNames(bplapply(levels(block), function(blc) {
        dat <- object[, block == blc]
        dat <- scater::normalize(dat)
        
        fit <- scran::trendVar(dat, method = "spline", min.mean = 0.001, use.spikes = FALSE)
        fit$trend <- makeTechTrend(x = dat)
        
        decomp <- decomposeVar(fit = fit) 
      }, BPPARAM = MulticoreParam(n_cores)), levels(block))
      
      #-Combine decomposed variances
      decomp <- combine.var.list(block_vars, method = "z", weighted = TRUE)
      
      #-Get HVGs
      hvgs <- rownames(decomp[decomp$bio > bio_thresh & decomp$FDR < fdr_thresh,])
      
      if(save_hvgs) {
        message(paste0("saving hvgs... (./results/clustering/", set, "/", set, "_hvgs.rds)"))
        saveRDS(hvgs, file = paste0("./results/clustering/", set, "/", set, "_hvgs.rds"))
      }
      
      #-Do PCA (including permutations to calculate number of PCs to keep)
      message("running pca and permutations...")
      object <- do.par.multi.batch.PCA(object, exprs_use = "logcounts", block = block, subset_row = hvgs,
                                       niters = pca_iters, approximate = TRUE, BPPARAM = MulticoreParam(n_cores, progressbar = TRUE), ...)
      if(is.null(sig_pcs)) sig_pcs <- attr(object@reducedDims$PCA, "npcs")
      
      if(save_objects) {
        message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
        saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
      }
      
    } else {
      #-Get Poisson trend
      fit <- trendVar(object, use.spikes = FALSE, method = "spline", min.mean = 0.001)
      fit$trend <- makeTechTrend(x = object)
      
      #-Decompose variances
      decomp <- decomposeVar(fit = fit)
      
      #-Get HVGs
      hvgs <- rownames(decomp[decomp$bio > bio_thresh & decomp$FDR < fdr_thresh,])
      
      if(save_hvgs){
        message(paste0("saving hvgs... (./results/clustering/", set, "/", set, "_hvgs.rds)"))
        saveRDS(hvgs, file = paste0("./results/clustering/", set, "/", set, "_hvgs.rds"))
      }
      
      #-Do PCA
      message("running pca and permutations...")
      object <- do.parallel.PCA(object, subset.row = hvgs, niters = pca_iters, approximate = TRUE, BPPARAM = MulticoreParam(n_cores), ...)
      rownames(object@reducedDims$PCA) <- colnames(object)
      sig_pcs <- attr(object@reducedDims$PCA, "npcs")
      
      if(save_objects){
        message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
        saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
      }
    } 
    
  } else {
    
    if(is.null(sig_pcs)) {
      sig_pcs <- attr(object@reducedDims$PCA, "npcs")
    }
    
  }
  
  
  if(batch_correct) {
    
    if(run_batch_correct) {
      #-Do mutual nearest neighbour batch correction
      message("doing batch correction...")
      pca_list <- lapply(levels(block), function(blc) {
        object@reducedDims$PCA[block == blc, 1:sig_pcs]
      })
      
      mnn_res <- fastMNN.list(pca_list, k = k, cos.norm = FALSE, auto.order = TRUE, pc.input = TRUE, compute.variances = TRUE)
      
      if(!identical(colnames(object), rownames(mnn_res$corrected))) mnn_res$corrected <- mnn_res$corrected[colnames(object),]
      
      object@reducedDims$MNN <- mnn_res$corrected
      
      if(save_objects){
        message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
        saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
      }
    }
    
    dim_red_use <- "MNN"
  } else {
    dim_red_use <- "PCA"
  }
  
  if(run_gclust){
    #-Do graph clustering
    message("doing graph clustering...")
    object <- do.graph.clustering(object, cells_use = NULL, pcs_use = 1:sig_pcs, num_nn = k,
                                  do_jaccard = TRUE, method = "Louvain", dim_red_use = dim_red_use)
    metadata(object)$k <- k
    
    if(save_objects){
      message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
      saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
    }
  }
  
  if(get_graph){
    graph <- do.graph.clustering(object, cells_use = NULL, pcs_use = 1:sig_pcs, num_nn = k, 
                                 do_jaccard = TRUE, method = "Louvain", dim_red_use = dim_red_use, return_graph = TRUE)
    
    if(save_objects){
      message(paste0("saving graph... (./results/clustering/", set, "/", set, "_graph.rds)"))
      saveRDS(graph, file = paste0("./results/clustering/", set, "/", set, "_graph.rds"))
    }
  }
  
  if(get_doublet_scores){
    message("calculating doublet scores...")
    dblt_scores <- mclapply(levels(object$sample), function(samp){
      message(samp)
      dat <- object[,object$sample == samp]
      scores <- doubletCells(dat, k = k)
      filter <- isOutlier(scores, nmads = 3, type = "higher")
      return(data.frame(score = scores, filter = filter, row.names = colnames(dat)))
    }, mc.cores = n_cores)
    
    dblt_scores <- do.call(rbind.data.frame, dblt_scores)
    dblt_scores <- dblt_scores[colnames(object),]
    object$dblt_score <- dblt_scores$score
    object$dblt_fil <- dblt_scores$filter
    
    if(save_objects){
      message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
      saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
    }
  }
  
  
  if(run_umap) {
    message("calculating UMAP...")
    umap <-  umap::umap(object@reducedDims[[dim_red_use]][, 1:sig_pcs], n_neighbours = k)
    
    reducedDim(object, "UMAP") <- umap$layout
    metadata(object)$UMAP <- umap[names(umap) %in% c("knn", "config")]
    
    if(save_objects){
      message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
      saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
    }
  }
  
  
  if(run_tSNE) {
    message("calculating tSNE...")
    tSNE <- Rtsne.multicore::Rtsne.multicore(object@reducedDims[[dim_red_use]][, 1:sig_pcs], pca = FALSE, 
                                             perplexity = 50, num_threads = n_cores*2, theta = 0.0,
                                             check_duplicates = FALSE, max_iter = 2000, verbose = FALSE)
    
    reducedDim(object, "tSNE") <- tSNE$Y
    metadata(object)$tSNE <- tSNE[names(tSNE) != "Y"]
    
    if(save_objects){
      message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
      saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
    }
  }
  
  
  if(calc_markers){
    if(block_markers) {
      block_mark <- block
    } else {
      block_mark <- NULL
    }
    
    #-Get marker genes
    message("calculating marker genes...")
    message("scran markers...")
    if(de_method == "t_test") scran_markers <- scran::findMarkers(object, clusters = object$gclust, block = block_mark, direction = "up")
    if(de_method == "wilcox") scran_markers <- scran::overlapExprs(object, groups = object$gclust, block = block_mark, direction = "up")
    
    scran_markers <- lapply(scran_markers, function(mat) {
      mat <- as.data.frame(mat[which(mat$FDR < 0.1),])
      #mat <- as.data.frame(mat[order(mat$FDR),], stringsAsFactors = FALSE)
    })
    
    if(save_objects){
      message(paste0("saving scran markers... (./results/clustering/", set, "/", set, "_scran_markers.xlsx)"))
      WriteXLS::WriteXLS(scran_markers, paste0("./results/clustering/", set, "/", set, "_scran_markers.xlsx"), row.names = TRUE)
    }
    
    message("hierarchical scran markers...")
    object <- calc.lev.scran.markers(object, clust_col = "gclust", exprs_use = "logcounts", de_method = de_method, 
                                     block = block_mark, direction = "up", n_cores = n_cores)
    
    if(save_objects){
      message(paste0("saving gene info... (./results/clustering/", set, "/", set, "_gene_info.xlsx)"))
      writexl::write_xlsx(as.data.frame(rowData(object)), paste0("./results/clustering/", set, "/", set, "_gene_info.xlsx"))
    }
    
    if(save_objects){
      message(paste0("saving object... (./data/Robjects/scRNA/clean/", set, ".rds)"))
      saveRDS(object, file = paste0("./data/Robjects/scRNA/clean/", set, ".rds"))
    }
  }
  
  
  if(get_clust_scores) {
    #-Run cluster permutation
    message("running cluster permutation analysis...")
    gclust_perm <- do.graph.clust.perm(object, iters = 500, pcs_use = 1:sig_pcs, do_jaccard = TRUE, 
                                       method = "Louvain", dim_red_use = dim_red_use, n_cores = n_cores, rand_seed = 12345)
    
    if(save_objects){
      message(paste0("saving cluster permutation analysis... (./results/clustering/", set, "/", set, "_gclust_perm.rds)"))
      saveRDS(gclust_perm, file = paste0("./results/clustering/", set, "/", set, "_gclust_perm.rds"))
    }
    
    #-Calculate stability and purity scores
    message("calculating stability and purity scores...")
    scores <- calc.stblty.prty(object, gclust_perm, "gclust")
    
    if(save_objects){
      message(paste0("saving scores... (./results/clustering/", set, "/", set, "_gclust_scores.rds)"))
      saveRDS(scores, file = paste0("./results/clustering/", set, "/", set, "_gclust_scores.rds"))
    }
  }
  
  if(!save_objects){
    return(object)
  }
}






#------------------------------------------------
##Modified scran functions accepting lists
#------------------------------------------------

combine.var.list <- function(all.results, method = "fisher", weighted = TRUE) {
  #all.results <- list(...)
  ref <- NULL
  for (x in all.results) {
    if (is.null(ref)) {
      ref <- rownames(x)
    }
    else if (!identical(ref, rownames(x))) {
      stop("gene identities should be the same for all arguments in '...'")
    }
  }
  to.average <- c("mean", "total", "bio", "tech")
  output <- vector("list", length(to.average))
  names(output) <- to.average
  all.means <- lapply(all.results, FUN = "[[", i = "mean")
  n.cells <- scran:::.extract_weightings(all.results, "num.cells")
  output$mean <- scran:::.weighted_average_vals(all.means, n.cells, 
                                                weighted)
  resid.df <- scran:::.extract_weightings(all.results, "resid.df")
  for (val in to.average[-1]) {
    cur.vals <- lapply(all.results, FUN = "[[", i = val)
    output[[val]] <- scran:::.weighted_average_vals(cur.vals, resid.df, 
                                                    weighted)
  }
  p.combine <- lapply(all.results, FUN = "[[", i = "p.value")
  comb.args <- list(method = method)
  if (weighted) {
    comb.args$weights <- resid.df
  }
  p.final <- do.call(combinePValues, c(p.combine, comb.args))
  output <- do.call(DataFrame, output)
  p.final <- unname(p.final)
  output$p.value <- p.final
  output$FDR <- p.adjust(p.final, method = "BH")
  rownames(output) <- ref
  metadata(output) <- list(num.cells = sum(unlist(n.cells)), 
                           resid.df = sum(unlist(resid.df)))
  return(output)
}



fastMNN.list <- function (batches, k = 20, cos.norm = TRUE, ndist = 3, d = 50, approximate = FALSE, 
                          irlba.args = list(), subset.row = NULL, auto.order = FALSE, 
                          pc.input = FALSE, compute.variances = FALSE, assay.type = "logcounts", 
                          get.spikes = FALSE, BNPARAM = NULL, BPPARAM = SerialParam()) {
  #batches <- list(...)
  nbatches <- length(batches)
  if (nbatches < 2L) {
    stop("at least two batches must be specified")
  }
  if (!pc.input) {
    out <- scran:::.SCEs_to_matrices(batches, assay.type = assay.type, 
                                     subset.row = subset.row, get.spikes = get.spikes)
    batches <- out$batches
    subset.row <- out$subset.row
    if (!is.null(subset.row) && !identical(subset.row, seq_len(nrow(batches[[1]])))) {
      batches <- lapply(batches, "[", i = subset.row, drop = FALSE)
    }
    if (cos.norm) {
      batches <- lapply(batches, FUN = cosineNorm, mode = "matrix")
    }
    pc.mat <- scran:::.multi_pca(batches, approximate = approximate, 
                                 irlba.args = irlba.args, d = d, use.crossprod = TRUE)
  }
  else {
    scran:::.check_batch_consistency(batches, byrow = FALSE)
    pc.mat <- batches
  }
  var.kept <- rep(1, nbatches)
  re.order <- NULL
  if (!is.logical(auto.order)) {
    re.order <- as.integer(auto.order)
    if (!identical(sort(re.order), seq_len(nbatches))) {
      stop("integer 'auto.order' must contain a permutation of 1:nbatches")
    }
    auto.order <- FALSE
  }
  use.order <- !is.null(re.order)
  mnn.pairings <- vector("list", nbatches - 1L)
  for (bdx in 2:nbatches) {
    if (auto.order) {
      if (bdx == 2L) {
        d.out <- scran:::.define_first_merge(pc.mat, k = k, BNPARAM = BNPARAM, 
                                             BPPARAM = BPPARAM)
        refdata <- pc.mat[[d.out$first]]
        curdata <- pc.mat[[d.out$second]]
        mnn.sets <- d.out$pairs
        precomp <- d.out$precomputed
        processed <- c(d.out$first, d.out$second)
      } else {
        d.out <- scran:::.define_next_merge(refdata, pc.mat, 
                                            processed, precomp, k = k, BNPARAM = BNPARAM, 
                                            BPPARAM = BPPARAM)
        curdata <- pc.mat[[d.out$other]]
        mnn.sets <- d.out$pairs
        processed <- c(processed, d.out$other)
      }
    } else {
      if (bdx == 2L) {
        ref.idx <- if (use.order) 
          re.order[1]
        else 1L
        refdata <- pc.mat[[ref.idx]]
        processed <- ref.idx
      }
      cur.idx <- if (use.order) 
        re.order[bdx]
      else bdx
      curdata <- pc.mat[[cur.idx]]
      mnn.sets <- scran:::find.mutual.nn(refdata, curdata, k1 = k, 
                                         k2 = k, BNPARAM = BNPARAM, BPPARAM = BPPARAM)
      processed <- c(processed, cur.idx)
    }
    ave.out <- scran:::.average_correction(refdata, mnn.sets$first, 
                                           curdata, mnn.sets$second)
    overall.batch <- colMeans(ave.out$averaged)
    if (compute.variances) {
      var.before <- scran:::.compute_intra_var(refdata, curdata, 
                                               pc.mat, processed)
    }
    refdata <- scran:::.center_along_batch_vector(refdata, overall.batch)
    curdata <- scran:::.center_along_batch_vector(curdata, overall.batch)
    if (compute.variances) {
      var.after <- scran:::.compute_intra_var(refdata, curdata, 
                                              pc.mat, processed)
      var.kept[seq_len(bdx)] <- var.kept[seq_len(bdx)] * 
        var.after/var.before
    }
    re.ave.out <- scran:::.average_correction(refdata, mnn.sets$first, 
                                              curdata, mnn.sets$second)
    curdata <- scran:::.tricube_weighted_correction(curdata, re.ave.out$averaged, 
                                                    re.ave.out$second, k = k, ndist = ndist, BNPARAM = BNPARAM, 
                                                    BPPARAM = BPPARAM)
    mnn.pairings[[bdx - 1L]] <- DataFrame(first = mnn.sets$first, 
                                          second = mnn.sets$second + nrow(refdata))
    refdata <- rbind(refdata, curdata)
  }
  if (auto.order || use.order) {
    ordering <- vector("list", nbatches)
    last <- 0L
    for (idx in processed) {
      ncells <- nrow(pc.mat[[idx]])
      ordering[[idx]] <- last + seq_len(ncells)
      last <- last + ncells
    }
    ordering <- unlist(ordering)
    refdata <- refdata[ordering, , drop = FALSE]
    relocate <- ordering
    relocate[ordering] <- seq_along(ordering)
    for (x in seq_along(mnn.pairings)) {
      current <- mnn.pairings[[x]]
      current$first <- relocate[current$first]
      current$second <- relocate[current$second]
      mnn.pairings[[x]] <- current
    }
    var.kept[processed] <- var.kept
  }
  if (!is.null(names(batches))) {
    batch.labels <- names(batches)
  }
  else {
    batch.labels <- seq_along(batches)
  }
  ncells <- vapply(pc.mat, FUN = nrow, FUN.VALUE = 0L)
  batch.ids <- Rle(batch.labels, ncells)
  output <- list(corrected = refdata, batch = batch.ids, pairs = mnn.pairings, 
                 order = batch.labels[processed])
  if (!pc.input) {
    output$rotation <- metadata(pc.mat)$rotation
  }
  if (compute.variances) {
    output$lost.var <- 1 - var.kept
  }
  return(output)
}






#------------------------------------------------
##Functions to run multi-batch PCA (from scran) and run permutations
#------------------------------------------------

do.par.multi.batch.PCA <- function(object, exprs_use, block = NULL, subset_row = NULL, min_rank = 5, max_rank = 100, niters = 50, threshold = 0.1, 
                                   approximate = FALSE, irlba_args = list(), run_permutations = TRUE, BPPARAM = SerialParam()) {
  if(is.null(block)) stop("batches must be supplied")
  
  if(!is.factor(block)) stop("\'block\' must be a factor")
  
  mat_list <- lapply(levels(block), function(blc) {
    assay(object, exprs_use)[, block == blc]
  })
  
  for (idx in seq_along(mat_list)) {
    current <- DelayedArray(mat_list[[idx]])
    if (!is.null(subset_row)) {
      current <- current[subset_row, , drop = FALSE]
    }
    mat_list[[idx]] <- current
  }
  
  message("object pca...")
  svd_out <- do.multi.batch.pca(mat_list, d = max_rank, approximate = approximate, subset_row = NULL, 
                                use_crossprod = TRUE, BPPARAM = BPPARAM)
  original_d2 <- svd_out@metadata$d^2
  
  if(run_permutations) {
    message("permutations...")
    permuted <- bplapply(rep(max_rank, niters), FUN = run.parallel.mb.pca, 
                         mat_list = mat_list, approximate = approximate, use_crossprod = TRUE, 
                         subset_row = NULL, irlba_args = list(), BPPARAM = BPPARAM)
    permutations <- do.call(cbind, permuted)
    
    prop <- rowMeans(permutations >= original_d2)
    above <- prop > threshold
    if (!any(above)) {
      npcs <- length(above)
    } else {
      npcs <- min(which(above)) - 1L
    }
    npcs <- scran:::.keep_rank_in_range(npcs, min_rank, length(original_d2))
  }
  
  x <- as.matrix(do.call(cbind.data.frame, mat_list))
  var_exp <- original_d2/(ncol(x) - 1)
  all_var <- sum(rowVars(x))
  
  mat_out <- do.call(rbind, svd_out)
  
  if(!identical(colnames(object), rownames(mat_out))) mat_out <- mat_out[colnames(object),]
  
  attr(mat_out, "rotation") <- metadata(svd_out)$rotation
  attr(mat_out, "d") <- metadata(svd_out)$d
  attr(mat_out, "percent_var") <- var_exp/all_var
  
  if(run_permutations) {
    attr(mat_out, "perm_percent_var") <- t(permutations)/(ncol(x) - 1L)/all_var
    attr(mat_out, "npcs") <- npcs
  }
  
  object@reducedDims$PCA <- mat_out
  
  return(object)
}



do.multi.batch.pca <- function(mat_list, subset_row = NULL, d = 50, approximate = FALSE, 
                               irlba_args = list(), use_crossprod = FALSE, BPPARAM = SerialParam()) {
  all_centers <- 0
  for (idx in seq_along(mat_list)) {
    current <- DelayedArray(mat_list[[idx]])
    if (!is.null(subset_row)) {
      current <- current[subset_row, , drop = FALSE]
    }
    centers <- DelayedMatrixStats::rowMeans2(current)
    all_centers <- all_centers + centers
    mat_list[[idx]] <- current
  }
  
  all_centers <- all_centers/length(mat_list)
  centered <- scaled <- mat_list
  
  for (idx in seq_along(mat_list)) {
    current <- mat_list[[idx]]
    current <- current - all_centers
    centered[[idx]] <- current
    current <- current/sqrt(ncol(current))
    scaled[[idx]] <- t(current)
  }
  
  if (d > min(ncol(scaled[[1]]), sum(vapply(scaled, FUN = nrow, FUN.VALUE = 0L)))) {
    stop("'d' is too large for the number of cells and genes")
  }
  
  if (use_crossprod) {
    #Set to true in multiBatchPCA
    svd_out <- scran:::.fast_svd(scaled, nv = d, irlba.args = irlba_args, 
                                 approximate = approximate, BPPARAM = BPPARAM)
  } else {
    combined <- as.matrix(do.call(rbind, scaled))
    if (!approximate) {
      svd_out <- svd(combined, nu = 0, nv = d)
    } else {
      svd_out <- do.call(irlba::irlba, c(list(A = combined, nu = 0, nv = d), irlba_args))
    }
  }
  
  final <- centered
  for (idx in seq_along(centered)) {
    final[[idx]] <- as.matrix(t(centered[[idx]]) %*% svd_out$v)
  }
  
  final <- as(final, "List")
  metadata(final) <- list(rotation = svd_out$v, d = svd_out$d)
  return(final)
}



run.parallel.mb.pca <- function(mat_list, ...) {
  mat_list_perm <- lapply(mat_list, function(mat) .Call(scran:::cxx_shuffle_matrix, mat))
  #mat_list_perm <- bplapply(mat_list, function(mat) t(apply(mat, 1, sample)), BPPARAM = MulticoreParam(2))
  out <- do.multi.batch.pca(mat_list = mat_list_perm, ..., BPPARAM = MulticoreParam(2))
  out@metadata$d^2
}



#------------------------------------------------
##Function that parallelizes PCA as implemented in scran
#------------------------------------------------

do.parallel.PCA <- function(object, assay.type = "logcounts", subset.row = NULL, get.spikes = FALSE, 
                            min.rank = 5, max.rank = 100, niters = 50, threshold = 0.1, 
                            approximate = FALSE, irlba.args = list(), run_permutations = TRUE, BPPARAM = SerialParam()) {
  
  subset.row <- scran:::.SCE_subset_genes(subset.row = subset.row, x = object, get.spikes = get.spikes)
  
  x <- assay(object, i = assay.type)
  
  x0 <- x
  if (!is.null(subset.row)) {
    subset.row <- scran:::.subset_to_index(subset.row, x, byrow = TRUE)
    x <- x[subset.row, , drop = FALSE]
  }
  y <- t(x)
  
  svd.out <- scran:::.centered_SVD(y, max.rank, approximate = approximate, 
                                   extra.args = irlba.args, keep.left = TRUE, keep.right = FALSE)
  original.d2 <- svd.out$d^2
  
  if(run_permutations){
    permuted <- bplapply(rep(max.rank, niters), FUN = scran:::.parallel_PA, 
                         y = y, approximate = approximate, extra.args = irlba.args, 
                         BPPARAM = BPPARAM)
    
    permutations <- do.call(cbind, permuted)
    prop <- rowMeans(permutations >= original.d2)
    above <- prop > threshold
    
    if (!any(above)) {
      npcs <- length(above)
    } else {
      npcs <- min(which(above)) - 1L
    }
    
    npcs <- scran:::.keep_rank_in_range(npcs, min.rank, length(original.d2))
  }
  
  out.val <- scran:::.svd_to_pca(svd.out, max.rank)
  
  var.exp <- original.d2/(ncol(x) - 1)
  all.var <- sum(DelayedMatrixStats::rowVars(DelayedArray(x)))
  
  attr(out.val, "percentVar") <- var.exp/all.var
  attr(out.val, "rotation") <- svd.out$v
  attr(out.val, "d") <- svd.out$d
  
  if(run_permutations){
    attr(out.val, "permuted.percentVar") <- t(permutations)/(ncol(x) - 1L)/all.var
    attr(out.val, "npcs") <- npcs
  }
  
  reducedDim(object, "PCA") <- out.val
  return(object)
}






#------------------------------------------------
##Function to perform Louvain-Jaccard clustering (from Shekhar et al)
#------------------------------------------------

do.graph.clustering <- function(object, cells_use = NULL, pcs_use = 1:10, num_nn = 30, do_jaccard = FALSE, method = "Louvain",
                                dim_red_use = "PCA", clust_name = "gclust", return_graph = FALSE) {
  if (do_jaccard){
    weights <- TRUE
    method_print <- paste0(method,"-","Jaccard")
  } else {
    weights <- NULL
    method_print <- method
  }
  
  message(paste0("Performing ", method_print, " clustering. Using ", num_nn, " nearest neighbors, and ", max(pcs_use), " PCs"))
  
  if (is.null(cells_use)){
    data_use <- object@reducedDims[[dim_red_use]][, pcs_use]
  } else {
    data_use <- object@reducedDims[[dim_red_use]][cells_use, pcs_use]
    
    object <- object[,cells_use]
  } 
  
  adj <- get.edges(data_use, nn = num_nn, do_jaccard = do_jaccard)
  
  g <- graph.adjacency(adj, mode = "undirected", weighted = weights)
  
  if(return_graph) {
    message("returning graph")
    return(g)
  }
  
  if (method == "Louvain") graph_out <- cluster_louvain(g)   
  if (method == "Infomap") graph_out <- cluster_infomap(g)
  
  clust_assign <- factor(graph_out$membership, levels = sort(unique(graph_out$membership)))
  names(clust_assign) <- graph_out$names
  
  k <- order(table(clust_assign), decreasing = TRUE)
  
  new_levels <- rep(1, length(unique(graph_out$membership)))
  new_levels[k] <- 1:length(unique(graph_out$membership))
  levels(clust_assign) <- new_levels
  clust_assign <- factor(clust_assign, levels = 1:length(unique(graph_out$membership)))
  
  message("Outputting clusters ..")
  colData(object)[, clust_name] <- clust_assign
  
  return(object) 
}


get.edges <- function(X, nn = 30 , do_jaccard = TRUE) {
  nearest <- nn2(X, X, k = nn+1, treetype = "bd", searchtype = "priority")
  message("Found nearest neighbors")
  
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] 
  #Convert to a similarity score
  nearest$nn.sim <- 1*(nearest$nn.dists >= 0 )
  
  edges <- melt(t(nearest$nn.idx)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
  edges$B <- edges$C
  edges$C <- 1
  
  #Remove repetitions
  edges <- unique(transform(edges, A = pmin(A,B), B = pmax(A,B)))
  
  if (do_jaccard){
    NN <- nearest$nn.idx
    jaccard_dist <- apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges$C <- jaccard_dist
    edges <- subset(edges, C != 0)
    edges$C <- edges$C/max(edges$C)
  }
  
  #adj <- matrix(0, nrow=nrow(X), ncol=nrow(X))
  adj <- Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse = TRUE)
  rownames(adj) <- rownames(X); colnames(adj) <- rownames(X)
  adj[cbind(edges$A, edges$B)] <- edges$C
  adj[cbind(edges$B, edges$A)] <- edges$C
  
  return(adj)
}



markers.binom <- function(object, clust_col, clust1, clust2 = NULL, genes_use = NULL, effect_size = log(2)) {
  if(is.null(genes_use)) genes_use <- rownames(object)
  
  cells1 <- rownames(colData(object)[colData(object)[, clust_col] == clust1, ])
  
  if (is.null(clust2)) {
    clust2 <- "rest"
    cells2 <- rownames(colData(object)[colData(object)[, clust_col] != clust1, ])
  } else {
    cells2 <- rownames(colData(object)[colData(object)[, clust_col] == clust2, ])
  }
  
  norm_mat <- logcounts(object)[genes_use, c(cells1, cells2), drop = FALSE]
  
  result <- binomcount.test(norm_mat, cells1, cells2, effect_size)
  
  if(length(cells1) == 1) {
    x <- exprs(object)[rownames(result), cells1]
    pos_frac1 <- ifelse(x > 0, 1, 0)
  } else {
    pos_frac1 <- apply(exprs(object)[rownames(result), cells1, drop = FALSE], 1, function(x) round(sum(x > 0)/length(x), 2))
  }
  
  if(!clust2 == "rest" & length(cells2) == 1){
    x <- exprs(object)[rownames(result), cells2]
    pos_frac2 <- ifelse(x > 0, 1, 0)
  } else {
    pos_frac2 <- apply(exprs(object)[rownames(result), cells2, drop = FALSE], 1, function(x) round(sum(x > 0)/length(x), 2))
  }
  
  if (clust2 == "rest"){
    genes_include <- pos_frac1 >= 0.1
  } else {
    genes_include <- (pos_frac1 >= 0.1) | (pos_frac2 >= 0.1)
  }
  
  result <- result[genes_include,]
  result <- result[order(abs(result$log_effect), decreasing=TRUE),]
  
  #Mean number of transcripts per cell
  if(length(cells1) == 1) {
    mean_trans1 <- counts(object)[rownames(result), cells1]
  } else {
    mean_trans1 <- apply(counts(object)[rownames(result), cells1, drop = FALSE], 1, function(x) round(mean(x),3))
  }
  
  if(!clust2 == "rest" & length(cells2) == 1) {
    mean_trans2 <- counts(object)[rownames(result), cells2]
  } else {
    mean_trans2 <- apply(counts(object)[rownames(result), cells2, drop = FALSE], 1, function(x) round(mean(x),3))
  }
  
  result[, paste0("mean_trans_", clust1)] = mean_trans1
  result[, paste0("mean_trans_", clust2)] = mean_trans2
  
  return(result)
}




binomcount.test <- function(norm_mat, cells1, cells2, effect_size) {
  #Test for enrichments in cluster #1
  if(length(cells2) == 1) {
    m <- sum(cells2 > 0)
  } else {
    m <- apply(norm_mat[, cells2, drop = FALSE], 1, function(x) sum(x > 0)) #Number of cells expressing marker in cluster #2
  }
  
  m1 <- m; m1[m == 0] <- 1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
  
  if(length(cells1) == 1) {
    n <- sum(cells1 > 0)
  } else {
    n <- apply(norm_mat[, cells1, drop = FALSE], 1, function(x) sum(x > 0)) #Number of cells expressing marker in cluster #1
  }
  
  #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
  pv1 <- pbinom(n, length(cells1), m1/length(cells2), lower.tail = FALSE) + dbinom(n, length(cells1), m1/length(cells2))
  
  log_fold_express <- log((n*length(cells2))/(m*length(cells1))) #log proportion of expressing cells
  d1 <- data.frame(log_effect = log_fold_express, pval = pv1)
  d1 <- subset(d1, log_effect >= effect_size)
  d1 <- d1[order(d1$pval, decreasing = FALSE),]
  
  #Enrichments in cells.2
  n1 <- n; n1[n == 0] <- 1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
  
  #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
  pv2 <- pbinom(m, length(cells2), n1/length(cells1), lower.tail=FALSE) + dbinom(m, length(cells2), n1/length(cells1))
  
  d2 <- data.frame(log_effect = log_fold_express, pval = pv2)
  d2 <- subset(d2, log_effect <= -effect_size)
  d2 <- d2[order(d2$pval, decreasing=FALSE),]
  
  d <- rbind(d1, d2)
  d <- d[order(d$pval, decreasing=FALSE),]
  return(d)
}





#------------------------------------------------
##Functions to do permutation analysis on Louvain-Jaccard clusters
#------------------------------------------------

graph.clust.perm <- function(object, cells_use = NULL, pcs_use = 1:10, do_jaccard = FALSE, method = "Louvain", dim_red_use = "PCA") {
  if (do_jaccard){
    weights <- TRUE
  } else {
    weights <- NULL
  }
  
  #Randomly select 85% of cells
  if (is.null(cells_use)){
    cells_sample <- sample(rownames(object@reducedDims$PCA), ceiling(0.85*nrow(object@reducedDims$PCA)))
  } else {
    cells_sample <- sample(cells_use, round(0.85*length(cells_use)))
  } 
  
  data_use <- object@reducedDims[[dim_red_use]][cells_sample[order(as.numeric(cells_sample))], pcs_use]
  
  #Build k-NN graph with a randomly chosen k
  num_nn <- round(runif(1, min = 15, max = 100))
  adj <- get.edges.perm(data_use, nn = num_nn, do_jaccard = do_jaccard)
  
  g <- graph.adjacency(adj, mode = "undirected", weighted = weights)
  if (method == "Louvain") graph_out <- cluster_louvain(g)   
  if (method == "Infomap") graph_out <- cluster_infomap(g)
  
  clust_assign <- factor(graph_out$membership, levels = sort(unique(graph_out$membership)))
  names(clust_assign) <- graph_out$names
  
  k <- order(table(clust_assign), decreasing = TRUE)
  
  new_levels <- rep(1, length(unique(graph_out$membership)))
  new_levels[k] <- 1:length(unique(graph_out$membership))
  levels(clust_assign) <- new_levels
  clust_assign <- factor(clust_assign, levels = 1:length(unique(graph_out$membership)))
  
  return(clust_assign) 
}



get.edges.perm <- function(X, nn = 30 , do_jaccard = TRUE) {
  nearest <- nn2(X, X, k = nn+1, treetype = "bd", searchtype = "priority")
  
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] 
  
  edges <- melt(t(nearest$nn.idx)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
  edges$B <- edges$C; edges$C <- 1
  
  #Remove repetitions
  edges <- unique(transform(edges, A = pmin(A,B), B = pmax(A,B)))
  
  #Remove 5% of edges
  edges_5pct <- round(0.05*nrow(edges))
  edges_prtrb <- edges[-sample.int(nrow(edges), edges_5pct),]
  
  #Add 5% spurious edges
  edges_add <- sapply(1:edges_5pct, function(pair) {
    cells <- sample(nrow(X), 2)
  })
  edges_add <- setNames(cbind(as.data.frame(t(edges_add)), as.integer(1)), c("A", "B", "C"))
  
  edges_prtrb <- rbind(edges_prtrb, edges_add)
  rownames(edges_prtrb) <- 1:nrow(edges_prtrb)
  
  if (do_jaccard){
    NN <- nearest$nn.idx
    jaccard_dist <- apply(edges_prtrb, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges_prtrb$C <- jaccard_dist
    edges_prtrb <- subset(edges_prtrb, C != 0)
    edges_prtrb$C <- edges_prtrb$C/max(edges_prtrb$C)
    
    #Add multiplicative noise to each edge
    mult <- runif(1, min = 0.6, max = 1.66)
    edges_prtrb$C <- edges_prtrb$C * mult
  }
  
  #adj <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  adj <- Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse = TRUE)
  rownames(adj) <- rownames(X); colnames(adj) <- rownames(X)
  adj[cbind(edges_prtrb$A, edges_prtrb$B)] <- edges_prtrb$C
  adj[cbind(edges_prtrb$B, edges_prtrb$A)] <- edges_prtrb$C
  return(adj)
}



do.graph.clust.perm <- function(object, iters, cells_use = NULL, pcs_use = 1:10, do_jaccard = FALSE, method = "Louvain",
                                dim_red_use = "PCA", n_cores = 4, rand_seed = 1) {
  require(doParallel)
  require(doRNG)
  
  set.seed(rand_seed)
  
  cl <- makeCluster(n_cores, outfile = "")
  registerDoParallel(cl)
  
  clust_assigns <- foreach(i = 1:iters, .packages = c('RANN', 'igraph', 'reshape2')) %dorng% {
    source("/projects/vleblanc_prj/GBM_organoids//code/analysis/functions_analysis.R")
    graph.clust.perm(object, cells_use = cells_use, pcs_use = pcs_use, do_jaccard = do_jaccard, method = method, dim_red_use = dim_red_use)
  }
  
  stopCluster(cl)
  
  return(clust_assigns)
}





#------------------------------------------------
##Function to calculate stability and purity of Louvain-Jaccard clusters
#------------------------------------------------


calc.stblty.prty <- function(object, clust_perm, clust_name) {
  #Get clusters
  clusts <- setNames(colData(object)[, clust_name], rownames(colData(object)))
  
  #Get number of trials
  trials <- length(clust_perm)
  
  #Calculate stability and purity for each cluster
  clust_vals <- sapply(levels(clusts), function(k) {
    #Calculate stability and purity values for each trial
    vals_i <- sapply(seq(1:trials), function(i) {
      
      #Get clusters for realization i
      trial_clusts <- clust_perm[[i]]
      
      ##STABILITY
      #Calculate Htot
      Htot <- sapply(levels(trial_clusts), function(j) {
        #Get p(j)
        p_j <- table(trial_clusts)[j] / length(trial_clusts)
        return(p_j * log(p_j))
      })
      
      Htot <- -sum(Htot, na.rm = TRUE)
      
      #Calculate Hk
      Hk <- sapply(levels(trial_clusts), function(j) {
        #Get cells in cluster k (original solution)
        cells_in_k <- names(clusts[clusts == k])
        
        #Get cells from cluster k for realization i
        cells_in_k_i <- trial_clusts[cells_in_k]
        
        #Get f(j)
        f_j <- table(cells_in_k_i)[j] / length(cells_in_k_i)
        return(f_j * log(f_j))
      })
      
      Hk <- -sum(Hk, na.rm = TRUE)
      stability_i <- Hk / Htot
      
      
      ##PURITY
      
      #Get cells in cluster k (original solution)
      cells_in_k <- names(clusts[clusts == k])
      
      #Get cells from cluster k for realization i
      cells_in_k_i <- trial_clusts[cells_in_k]
      
      #Get cluster with max cells from cluster k
      clust_max <- names(which(table(cells_in_k_i) == max(table(cells_in_k_i))))
      if(length(clust_max) > 1) clust_max <- clust_max[1]
      
      cells_in_clust_max <- trial_clusts[trial_clusts == clust_max]
      
      #Get purity for realization i
      purity_i <- sum(names(cells_in_k_i) %in% names(cells_in_clust_max)) / length(cells_in_clust_max)
      
      return(c(stability = stability_i, purity = purity_i))
    })
    
    vals <- data.frame(t(vals_i))
    
    stability <- 1 - ((1/trials) * sum(vals$stability))
    purity <- (1/trials) * sum(vals$purity)
    
    return(c(stability = stability, purity = purity))
  })
  
  clust_vals <- data.frame(t(clust_vals))
}



#------------------------------------------------
##Function to calculate cluster co-occurrence from permutation tests
#------------------------------------------------

gclust.perm.co.occur <- function(object, clust_perm, n_cores){
  co_mats <- mclapply(clust_perm, function(perm_iter) {
    dat <- data.frame(cell = colnames(object), clust = NA, row.names = colnames(object), stringsAsFactors = FALSE)
    dat[names(perm_iter), "clust"] <- perm_iter
    
    tbl <- Matrix::Matrix(table(dat), sparse = TRUE)
    co_occ <- Matrix::tcrossprod(tbl)
  }, mc.cores = n_cores)
  
  full_mat <- Reduce('+', co_mats)
  full_mat <- full_mat / length(clust_perm)
  full_mat <- full_mat[colnames(object), colnames(object)]
}





#------------------------------------------------
##Function to identify marker genes by levels using the scran approach
#------------------------------------------------

calc.lev.scran.markers <- function(object, clust_col = "gclust", exprs_use = "logcounts", de_method = c("t_test", "wilcox"), n_cores = 4, ...) {
  
  if(de_method == "t_test") markers_full <- scran::pairwiseTTests(assay(object, exprs_use), clusters = colData(object)[,clust_col],
                                                                  BPPARAM = MulticoreParam(n_cores), ...)
  if(de_method == "wilcox") markers_full <- scran::pairwiseWilcox(assay(object, exprs_use), clusters = colData(object)[,clust_col],
                                                                  BPPARAM = MulticoreParam(n_cores), ...)
  
  cl <- makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  num_genes <- nrow(object)
  inds <- split(seq_len(num_genes), sort(rep_len(seq_len(n_cores), num_genes)))
  
  markers <- foreach(l = seq_along(inds), .combine = rbind, .packages = "Matrix") %dopar% {
    genes <- rownames(object)[inds[[l]]]
    scores_inds <- sapply(genes, get.scran.markers.levs, res = markers_full, clusts = colData(object)[,clust_col])
    
    t(scores_inds)
  }
  
  markers <- as.data.frame(markers)
  
  if(identical(rownames(markers), rownames(object))) {
    f_data <- rowData(object)
    
    f_data[, paste0(clust_col, "_scran_levmarkers_level")] <- markers$level
    f_data[, paste0(clust_col, "_scran_levmarkers_score")] <- markers$score
    f_data[, paste0(clust_col, "_scran_levmarkers_clusts_high")] <- markers$clusts_high
    
    rowData(object) <- f_data
  }
  
  return(object)
}




get.scran.markers.levs <- function(gene, res, clusts, FDR_cutoff = 0.05) {
  if(!is.factor(clusts)) stop("\'clusts\' must be a factor")
  
  #For a gene, get results from all pairwise comparisons
  dat <- do.call(rbind.data.frame, lapply(res$statistics, function(pair) as.data.frame(pair[gene, , drop = FALSE])))
  dat <- cbind(res$pairs, dat)
  
  #Keep significant comparisons
  dat <- dat[which(dat$FDR < FDR_cutoff),]
  
  #Get level info
  lev <- NA
  clusts_high <- NA
  score <- NA
  
  num_clusts <- length(levels(clusts))
  clusts_named <- setNames(sort(1:(num_clusts-1), decreasing = TRUE), 1:(num_clusts-1))
  
  for(i in sort(seq_len(num_clusts-1), decreasing = TRUE)) {
    if(sum(table(dat$first) >= i) >= as.numeric(names(clusts_named[clusts_named == i]))) {
      clusts_high <- names(table(dat$first)[table(dat$first) >= i])
      
      if(length(clusts_high) == num_clusts) {
        clusts_high <- NA
        break
      }
      
      lev <- paste0("lev", names(clusts_named[clusts_named == i]))
      break
    }
  }
  
  if(nrow(dat) > 0) {
    score <- -log10(max(dat$FDR[dat$first %in% clusts_high], na.rm = TRUE))
  }
  
  c(level = lev, score = score, 
    clusts_high = paste(clusts_high, collapse = ","))
}