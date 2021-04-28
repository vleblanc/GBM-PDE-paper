# ---------------------------------------------------------------
## Functions to load 10X CNV data into a combined object
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

load.10x.cnv <- function(id, dirs, tum_levels, samp_levels, merged_calls = T, filter_cnvs = F, size_thresh = 1000000, conf_thresh = 10,
                         resume = F, return_unsplit_grmat = F, n_cores) {
  require(SingleCellExperiment)
  require(BiocParallel)
  
  if(n_cores == 1) {
    BPPARAM <- SerialParam()
  } else {
    if(n_cores > length(dirs)) {
      n_cores <- length(dirs)
    }
    BPPARAM <- MulticoreParam(n_cores)
  }
  
  #Get copy number and confidence matrices for individual samples
  all_samples_indiv <- bplapply(names(dirs), function(samp) {
    #Get file name for matrices
    mats_file <- paste0("./data/Robjects/scDNA/", samp, "_mats.rds")
    
    if(all(file.exists(mats_file) & resume)) {
      #If matrices exist, load and continue
      
      message(paste0("loading previously generated matrices for ", samp, "..."))
      mats <- readRDS(mats_file)
      return(mats)
      
    } else {
      mat_gr <- get.segments(dir = dirs[samp], merged_calls = merged_calls)
      
      mats <- get.mats(dir = dirs[samp], mat_gr = mat_gr)
    }
  }, BPPARAM = BPPARAM)
  
  #Combine samples into a full dataset
  message("combining samples...")
  
  if(return_unsplit_grmat) {
    all_grmats <- do.call(c, all_samples_indiv)
    return(all_grmats)
  }
  
  nsets <- length(all_samples_indiv)
  full_cn_dat <- vector("list", nsets)
  full_conf_dat <- vector("list", nsets)
  bin_info_list <- vector("list", nsets)
  cell_info_list <- vector("list", nsets)
  
  for(i in seq_len(nsets)) {
    samp <- all_samples_indiv[[i]]
    
    full_cn_dat[[i]] <- samp$cn
    full_conf_dat[[i]] <- samp$conf
    bin_info_list[[i]] <- data.frame(chr = factor(sapply(strsplit(rownames(samp$cn), "_"), "[", 1),
                                                  levels = paste0("chr", c(1:22, "X", "Y"))),
                                     start_pos = sapply(strsplit(rownames(samp$cn), "_"), "[", 2),
                                     end_pos = sapply(strsplit(rownames(samp$cn), "_"), "[", 3))
    cell_info_list[[i]] <- data.frame(barcode = gsub("-1", paste0("-", i), colnames(samp$cn)),
                                      sample = names(dirs)[i],
                                      cell_id = paste(samp$metrics, i, sep = "-"),
                                      samp$metrics[,3:ncol(samp$metrics)])
  }
  
  if(nsets > 1 && length(unique(bin_info_list)) != 1L) {
    stop("bin information differs between runs")
  }
  
  full_cn_dat <- do.call(cbind, full_cn_dat)
  full_conf_dat <- do.call(cbind, full_conf_dat)
  
  cell_info <- do.call(rbind, cell_info_list)
  colnames(full_cn_dat) <- colnames(full_conf_dat) <- cell_info$barcode
  
  sce <- SingleCellExperiment(list(copy_number = full_cn_dat, confidence = full_conf_dat), 
                              rowData = bin_info_list[[1]], colData = cell_info)
  
  #Add cell (column) info
  sce <- mutate(sce,
                tumour = factor(sapply(strsplit(as.character(sce$sample), "_"), "[", 1),
                                levels = tum_levels),
                region = droplevels(factor(sapply(strsplit(as.character(sce$sample), "_"), "[", 2), 
                                           levels = c("reg1", "reg2"))),
                source = droplevels(factor(sapply(strsplit(as.character(sce$sample), "_"), "[", 3), 
                                           levels = c("tis", "org"))))
  sce$sample <- factor(sce$sample, levels = samp_levels)
  
  #Save .rds file
  message(paste0("saving object... (./data/Robjects/scDNA/", id, ".rds)"))
  saveRDS(sce, file = paste0("./data/Robjects/scDNA/", id, ".rds"))
  
  return(sce)
}




get.segments <- function(dir, merged_calls = TRUE) {
  require(GenomicRanges)
  
  #Get per-cell metrics
  metrics <- read.csv(paste0(dir, "/per_cell_summary_metrics.csv"))
  ncells <- nrow(metrics)
  
  #Get 10X CNV output matrix
  if(merged_calls) {
    #Get merged CNV calls
    mat <- data.table::fread(paste0(dir, "/node_cnv_calls.bed"), data.table = F)
  } else {
    #Get unmerged CNV calls (unmappable regions have an event confidence of 0)
    mat <- data.table::fread(paste0(dir, "/node_unmerged_cnv_calls.bed"), data.table = F)
  }
  colnames(mat)[1] <- "chrom"
  #Add 1 to start position to make even bins
  mat$start <- mat$start + 1
  
  #Keep only CNV info for individual cells (info also included for clusters)
  mat <- mat[mat$id < ncells,]
  
  #Convert to GRanges object
  mat_gr <- makeGRangesFromDataFrame(mat, keep.extra.columns = T)
  
  #Add info
  mcols(mat_gr)$barcode <- metrics[mat_gr$id + 1, "barcode"]
  mcols(mat_gr)$size <- width(mat_gr)
  mcols(mat_gr)$ploidy <- round(metrics[mat_gr$id + 1, "mean_ploidy"])
  mcols(mat_gr)$cat <- factor(ifelse(mat_gr$copy_number == mat_gr$ploidy, "neutral", "cn_event"),
                              levels = c("neutral", "cn_event"))
  mcols(mat_gr)$event_dir <- ifelse(mat_gr$cat == "neutral", "neutral",
                                    ifelse(mat_gr$copy_number > mat_gr$ploidy, "gain", "loss"))
  mcols(mat_gr)$est_cnv_res_mb <- metrics[mat_gr$id + 1, "est_cnv_resolution_mb"]
  
  return(mat_gr)
}




get.mats <- function(dir, mat_gr) {
  require(GenomicRanges)
  
  #Get per-cell metrics
  metrics <- read.csv(paste0(dir, "/per_cell_summary_metrics.csv"))
  
  #Get mappable regions
  map_reg <- data.table::fread(paste0(dir, "/mappable_regions.bed"), data.table = F)
  colnames(map_reg)[1] <- "chrom"
  
  #Calculate total number of bases covered by mappable regions
  map_bases <- sum(map_reg$end - map_reg$start)
  
  #Add number of mappable bins to metrics table
  metrics$mappable_bins <- map_bases / 20000 #20kb bins
  
  #Add window size - number of bins needed to reach ~200 reads/bin
  metrics$window_size <- ceiling(200 / (metrics$num_mapped_dedup_reads / metrics$mappable_bins))
  
  #Get chromosome lengths
  chr_lengths <- read.table("/projects/vleblanc_prj/tools/refdata-GRCh38-1.0.0/fasta/chr_lengths.tsv")
  rownames(chr_lengths) <- chr_lengths[,1]
  
  seqlengths(mat_gr) <- setNames(chr_lengths[seqlevels(mat_gr), 2],
                                 seqlevels(mat_gr))
  
  #Get gaps at edges of chromosomes that were not covered and make new object with full ranges
  gaps <- gaps(mat_gr)
  gaps <- gaps[strand(gaps) == "*"]
  
  steps <- c(mat_gr, gaps)
  steps <- sort(steps)
  
  #Split all ranges into smallest step
  step <- min(width(mat_gr))
  steps <- unlist(slidingWindows(range(steps), width = step, step = step))
  
  #Split copy number object by cell
  mat_gr_bycell <- split(mat_gr, mat_gr$id)
  
  #Get copy number for each step for each cell
  cn_bycell <- GRangesList(lapply(mat_gr_bycell, function(cell) {
    #Get overlaps between query (steps) and subject (cell)
    overlaps <- findOverlaps(steps, cell, ignore.strand = T)
    
    #Make new GRanges object for cell with steps
    cell_steps <- steps
    
    #Get copy number values for each step
    cell_steps$copy_number <- rep(NA, length(steps))
    cell_steps$copy_number[queryHits(overlaps)] <- mcols(cell)[subjectHits(overlaps), "copy_number"]
    
    #Get event confidence values for each step
    cell_steps$event_confidence <- rep(NA, length(steps))
    cell_steps$event_confidence[queryHits(overlaps)] <- mcols(cell)[subjectHits(overlaps), "event_confidence"]
    
    return(cell_steps)
  }))
  
  #Collapse copy number and event confidence info into matrices
  cn_mat <- do.call(cbind.data.frame, lapply(cn_bycell, function(cell) mcols(cell)$copy_number))
  conf_mat <- do.call(cbind.data.frame, lapply(cn_bycell, function(cell) mcols(cell)$event_confidence))
  
  rownames(conf_mat) <- rownames(cn_mat) <- paste(as.character(seqnames(steps)), start(steps), end(steps), sep = "_")
  colnames(conf_mat) <- colnames(cn_mat) <- metrics$barcode
  
  mats <- list(cn = cn_mat, conf = conf_mat, metrics = metrics)
  
  return(mats)
}