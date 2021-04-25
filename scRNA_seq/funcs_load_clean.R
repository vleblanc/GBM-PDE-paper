##-----------------------------------------------
# Functions to load, clean, and normalize 10X scRNA-seq data
# Author: Veronique LeBlanc
##-----------------------------------------------


##-----------------------------------------------
# Function to load raw dataset if .rds file can be found, otherwise make it and save it
##-----------------------------------------------

load.raw.dat <- function(id, dirs, samp_levels, tum_levels, n_cores = 4, rand_seed = 1, ...) {
  require(parallel)
  require(doParallel)
  require(doRNG)
  require(DropletUtils)
  require(SingleCellExperiment)
  
  raw_dat_file <- paste0("./data/Robjects/scRNA/raw/", id, "_raw.rds")
  
  if(file.exists(raw_dat_file)){
    #Load .rds file if it exists
    raw_dat <- readRDS(raw_dat_file)
  } else {
    #If .rds file does not exist, make raw dataset and save as .rds file
    set.seed(rand_seed)
    
    #If there are fewer samples than cores used, reset the number of cores used
    if(length(dirs) < n_cores) {
      n_cores <- length(dirs)
    }
    
    #Get raw count matrices for each sample and identify cell-containing droplets
    cl <- parallel::makeCluster(n_cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    all_samples_indiv <- foreach(i = 1:length(dirs), .packages = c('DropletUtils', 'scater')) %dorng% {
      #Get sample info
      samp <- dirs[i]
      tum <- unlist(strsplit(samp, "/"))[4]
      tum <- unlist(strsplit(tum, "_"))[1]
      
      #If .rds file already exists, load it
      out_file = paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/raw_cells.rds")
      
      if(file.exists(out_file)) {
        raw_fil <- readRDS(out_file)
      } else {
        #Read raw count matrix
        raw <- DropletUtils::read10xCounts(unname(samp), col.names = TRUE)
        
        #Relabel rows with gene symbols (ensuring uniqueness for duplicated or missing symbols)
        rownames(raw) <- uniquifyFeatureNames(rowData(raw)$ID, rowData(raw)$Symbol)
        
        saveRDS(raw, file = out_file)
      }
      
      #Return counts matrix to combine into full dataset
      return(raw)
    }
    
    suppressMessages(stopCluster(cl))
    
    #Combine samples into a full dataset
    message("combining samples...")
    
    nsets <- length(all_samples_indiv)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)
    
    for(i in seq_len(nsets)) {
      samp <- all_samples_indiv[[i]]
      
      full_data[[i]] <- counts(samp)
      gene_info_list[[i]] <- rowData(samp)
      cell_info_list[[i]] <- data.frame(cell_name = gsub("-1", paste0("-", i), colnames(samp)),
                                        sample = names(dirs)[i])
    }
    
    if(nsets > 1 && length(unique(gene_info_list)) != 1L) {
      stop("gene information differs between runs")
    }
    
    gene_info <- gene_info_list[[1]]
    colnames(gene_info) <- c("id", "symbol")
    
    full_data <- do.call(cbind, full_data)
    
    cell_info <- do.call(rbind, cell_info_list)
    colnames(full_data) <- cell_info$cell_name
    
    raw_dat <- SingleCellExperiment(list(counts = full_data), rowData = gene_info, colData = cell_info)
    
    #Keep only genes that have at least one UMI
    keep_gene <- nexprs(raw_dat, detection_limit = 0, byrow = TRUE) >= 1
    raw_dat <- raw_dat[keep_gene, ]
    
    #Get gene (row) info
    message("adding gene info...")
    gtf <- read.table("/projects/vleblanc_prj/tools/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf", sep = "\t", skip = 5,
                      stringsAsFactors = FALSE)
    gtf <- gtf[gtf$V3 == "gene", c(1,4,5,9)]
    gtf$ensembl_gene_id <- sapply(strsplit(sapply(strsplit(as.character(gtf$V9), ";"), function(splt) splt[grepl("gene_id", splt)]), " "), "[", 2)
    gtf <- gtf[,c("V1", "V4", "V5", "ensembl_gene_id")]
    colnames(gtf)[1:3] <- c("chr", "start_pos", "end_pos")
    
    #Keep only genes in standard chromosomes
    std_chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                 "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")
    gtf <- gtf[gtf$chr %in% std_chr,]
    gtf$chr <- factor(gtf$chr, levels = std_chr)
    
    rowData(raw_dat) <- merge(rowData(raw_dat), gtf, by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
    rownames(raw_dat) <- uniquifyFeatureNames(rowData(raw_dat)$id, rowData(raw_dat)$symbol)
    
    #Add cell (column) info
    raw_dat <- mutate(raw_dat,
                      tumour = factor(sapply(strsplit(as.character(raw_dat$sample), "_"), "[", 1),
                                      levels = tum_levels),
                      region = droplevels(factor(sapply(strsplit(as.character(raw_dat$sample), "_"), "[", 2), 
                                                 levels = c("reg1", "reg2"))),
                      source = droplevels(factor(sapply(strsplit(as.character(raw_dat$sample), "_"), "[", 3), 
                                                 levels = c("tis", "org", "cell"))),
                      rep = factor(sapply(strsplit(as.character(raw_dat$sample), "_"), "[", 4)))
    raw_dat$sample <- factor(raw_dat$sample, levels = samp_levels)
    
    #Calculate QC metrics (treating mitochondrial genes as feature controls)
    message("calculating QC metrics...")
    is_mito <- grepl("MT-", rownames(raw_dat))
    is_ribo <- grepl("^RP[LS]", rowData(raw_dat)$symbol)
    feature_ctrls <- list(mito = rownames(raw_dat)[is_mito],
                          ribo = rownames(raw_dat)[is_ribo])
    raw_dat <- calculateQCMetrics(raw_dat, feature_controls = feature_ctrls, BPPARAM = MulticoreParam(n_cores))
    
    #Save .rds file
    message(paste0("saving raw object... (./data/Robjects/scRNA/clean/", id, ".rds)"))
    saveRDS(raw_dat, file = paste0("./data/Robjects/scRNA/raw/", id, "_raw.rds"))
  }
  return(raw_dat)  
}


##-----------------------------------------------
## Function to make QCd and normalized dataset
##-----------------------------------------------
make.qcd.dataset <- function(raw_dat, cell_min_genes = 250, cell_min_UMI = 500, cell_max_MT = NULL, 
                             gene_min_cells = 10, gene_min_UMI = 20, n_cores = 4, rand_seed = 1){
  require(scran)
  require(scater)
  require(irlba)
  
  set.seed(rand_seed)
  
  #Filter cells
  cells_keep <- colnames(raw_dat) %in% colnames(raw_dat)
  
  if(!is.null(cell_min_genes)) {
    cell_min_gene_fil <- raw_dat$total_features_by_counts >= cell_min_genes
    cells_keep <- cells_keep & cell_min_gene_fil
  }
  
  if(!is.null(cell_min_UMI)) {
    cell_min_UMI_fil <- raw_dat$total_counts >= cell_min_UMI
    cells_keep <- cells_keep & cell_min_UMI_fil
  }
  
  if(!is.null(cell_max_MT)) {
    cell_max_MT_fil <- raw_dat$pct_counts_MT <= cell_max_MT
    cells_keep <- cells_keep & cell_max_MT_fil
  }
  
  #Make object with cells that pass filters
  qc_dat <- raw_dat[, cells_keep]
  
  
  
  #Filter genes
  genes_keep <- rownames(qc_dat) %in% rownames(qc_dat)
  
  if(!is.null(gene_min_cells)) {
    gene_min_cells_fil <- rowData(qc_dat)$n_cells_by_counts >= gene_min_cells
    genes_keep <- genes_keep & gene_min_cells_fil
  }
  
  if(!is.null(gene_min_UMI)) {
    gene_min_UMI_fil <- rowData(qc_dat)$total_counts >= gene_min_UMI
    genes_keep <- genes_keep & gene_min_UMI_fil
  }
  
  #Make object with genes that pass filters
  qc_dat <- qc_dat[genes_keep, ]
  
  
  #print number of cells and genes removed by filter(s)
  print(setNames((dim(raw_dat) - dim(qc_dat)), c("Genes removed", "Cells removed")))
  
  #print number of cells and genes retained
  print(setNames(dim(qc_dat), c("Genes retained", "Cells retained")))
  
  #Normalize
  print("clustering cells...")
  qclust <- quickCluster(qc_dat, method = "igraph", min.mean = 0.1, BPPARAM = MulticoreParam(n_cores))
  
  print("computing sum factors...")
  qc_dat <- computeSumFactors(qc_dat, clusters = qclust, min.mean = 0.1, positive = FALSE, BPPARAM = MulticoreParam(n_cores))
  #If there are negative or 0 sum factors, return object and stop
  if(sum(sizeFactors(qc_dat) <= 0) > 0) {
    return(qc_dat)
    stop("0 or negative sum factors have been calculated; object has been saved")
  }
  
  print("normalizing...")
  qc_dat <- scater::normalize(qc_dat, log_exprs_offset = 1)
  
  #Calculate QC metrics (treating mitochondrial and ribosomal genes as feature controls)
  message("calculating QC metrics...")
  is_mito <- grepl("MT-", rownames(qc_dat))
  is_ribo <- grepl("^RP[LS]", rowData(qc_dat)$symbol)
  feature_ctrls <- list(mito = rownames(qc_dat)[is_mito],
                        ribo = rownames(qc_dat)[is_ribo])
  qc_dat <- calculateQCMetrics(qc_dat, feature_controls = feature_ctrls, BPPARAM = MulticoreParam(n_cores))
  
  return(qc_dat)
}