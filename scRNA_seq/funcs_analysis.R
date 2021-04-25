##-----------------------------------------------
# Functions to perform general analyses of 10X scRNA-seq data
# Author: Veronique LeBlanc
##-----------------------------------------------




#------------------------------------------------
##Function to run scUniFrac on already clustered data
#------------------------------------------------

sc.unifrac.multi <- function(object, group_by, clusts_name, tree, perm_iters, n_cores){
  require(doRNG)
  require(GUniFrac)
  
  groups <- setNames(colData(object)[,group_by], colnames(object))
  group_names = levels(groups)
  
  if(!is.factor(groups)) stop("\"group_by\" must be a factor")
  if(length(levels(groups)) == 2) stop("only two levels in the groups - use sc.unifrac.pair")
  
  clusts <- setNames(colData(object)[,clusts_name], colnames(object))
  num_clusts <- max(as.numeric(levels(clusts)))
  
  count_table <- table(groups, clusts)
  colnames(count_table) <- paste0("c", levels(clusts))
  
  tree_phylo <- as.phylo(tree$hclust)
  
  unifracs <- GUniFrac(count_table, tree_phylo, alpha = c(0, 0.5, 1))$unifracs
  
  dist_obs <- unifracs[, , "d_1"]
  rownames(dist_obs) <- colnames(dist_obs) <- rownames(count_table)
  
  cl <- parallel::makeCluster(n_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  perm_scores <- foreach(i = 1:perm_iters, .packages = c('scUnifrac', "GUniFrac")) %dorng% {
    groups_perm <- factor(setNames(sample(as.character(groups)), names(groups)),
                          levels = levels(groups))
    
    count_table_perm <- table(groups_perm, clusts)
    colnames(count_table_perm) <- paste0("c", levels(clusts))
    
    unifracs_perm <- GUniFrac(count_table_perm, tree_phylo, alpha = c(0, 0.5, 1))$unifracs
    
    return(unifracs_perm[, , "d_1"])
  }
  
  perm_scores_array <- array(unlist(perm_scores), dim = c(length(levels(groups)), length(levels(groups)), perm_iters))
  p_vals <- matrix(apply(apply(perm_scores_array, 3, ">=", dist_obs), 1, sum), nrow = length(levels(groups)), byrow = F)/perm_iters
  rownames(p_vals) <- colnames(p_vals) <- rownames(count_table)
  
  return(list(distance = dist_obs, count_table = count_table, pvalue = p_vals))
}





#------------------------------------------------
##Functions to calculate signature scores as per the Tirosh method (Neftel et al, Puram et al, etc)
#------------------------------------------------

get.sig.scores <- function(object, gene_sets, rank_method = c("aggregate", "mean"), n_cores = 2) {
  if(!is.list(gene_sets)) stop("\'gene_sets\' must be a list")
  
  if(n_cores == 1) {
    BPPARAM <- SerialParam()
  } else {
    if(n_cores > length(gene_sets)) n_cores <- length(gene_sets)
    BPPARAM <- MulticoreParam(n_cores)
  }
  
  #Keep only genes that are present in the dataset
  gene_sets <- lapply(gene_sets, function(gene_set) {
    genes_remove <- gene_set[!gene_set %in% rownames(object)]
    if(length(genes_remove) > 0) message(cat("removing the following genes: ", genes_remove))
    
    gene_set[gene_set %in% rownames(object)]
  })
  
  #If centred expression values are not present, calculate them and add them as an assay to the object
  if(!"centred_logcounts" %in% names(object@assays)) {
    assay(object, "centred_logcounts") <- t(scale(t(as.matrix(logcounts(object))), center = TRUE, scale = FALSE))
  }
  
  if(rank_method == "aggregate") {
    #Get aggregate expression for all genes and rank them
    agg_exp <- Matrix::rowSums(logcounts(object))
    agg_exp <- agg_exp[order(agg_exp, decreasing = TRUE)]
    
    #Split genes into 30 bins based on ranked aggregate expression values
    exp_bins <- data.frame(gene = names(agg_exp), bin = paste0("bin", cut(seq_along(agg_exp), 30, labels = FALSE)))
  }
  
  if(rank_method == "mean") {
    #Get mean expression for all genes and rank them
    mean_exp <- Matrix::rowMeans(logcounts(object))
    mean_exp <- mean_exp[order(mean_exp, decreasing = TRUE)]
    
    #Split genes into 30 bins based on ranked mean expression values
    exp_bins <- data.frame(gene = names(mean_exp), bin = paste0("bin", cut(seq_along(mean_exp), 30, labels = FALSE)))
  }
  
  #Get scores for each gene set
  scores <- setNames(bplapply(gene_sets, function(gene_set){
    #Get scores for the get set
    gene_set_scores <- colMeans(assay(object, "centred_logcounts")[gene_set,])
    
    #Get control set genes
    control_set <- c(sapply(gene_set, function(gene) {
      bin <- exp_bins[exp_bins$gene == gene, "bin"]
      bin_genes <- exp_bins[exp_bins$bin == bin, "gene"]
      genes_use <- sample(bin_genes, 100)
    }))
    
    #Get control set scores
    control_set_scores <- colMeans(assay(object, "centred_logcounts")[control_set,])
    
    #Get final scores
    scores <- gene_set_scores - control_set_scores
  }, BPPARAM = BPPARAM), names(gene_sets))
  
  #Combine scores into a dataframe
  scores <- as.data.frame(do.call(cbind, scores))
  
  #Get highest scoring gene set
  scores$max_set <- factor(apply(scores, 1, function(cell_scores) colnames(scores)[which.max(cell_scores)]),
                           levels = names(gene_sets))
  
  return(scores)
}


