# ---------------------------------------------------------------
## Scripts to load data into a combined object, QC filter cells and genes, and normalize expression
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

library(DropletUtils)
library(scater)

setwd("/projects/vleblanc_prj/GBM_organoids/")

# ---------------------------------------------------------------
## LOAD DATA
# ---------------------------------------------------------------

#Function to pull directories containing count matrices for each sample
mat_dir <- function(tumour, sample) paste0("./data/cellranger_v2_1/", tumour, "/", sample, "/outs/filtered_gene_bc_matrices/GRCh38/")

#Load raw counts into a SingleCellExperiment object
# Biological replicates are indicated by numbers in fourth position (e.g. tis_1 and tis_2 are two different pieces of tissue, org_1 and org_2 are two independent organoids)
# Technical replicates are indicated by periods (e.g. org_2.1 and org_2.2 were derived from the same cell suspension obtained from a single organoid)
# Cross-chip replicates are indicated by "r1" (for JK125) or "br"
all <- load.raw.dat(id = "all_200227", call_method = "none", 
                    n_cores = 15, rand_seed = 12345,
                    dirs = c(#JK124
                             JK124_reg1_tis_1 = mat_dir("JK124", "SI-GA-A6"), JK124_reg1_tis_2 =	mat_dir("JK124", "SI-GA-B6"),
                             JK124_reg1_org_1 =	mat_dir("JK124", "SI-GA-C6"), JK124_reg1_org_2.1 =	mat_dir("JK124", "SI-GA-D6"),
                             JK124_reg1_org_2.2 =	mat_dir("JK124", "SI-GA-E6"), JK124_reg2_tis_1 =	mat_dir("JK124", "SI-GA-F6"),
                             JK124_reg2_tis_2 =	mat_dir("JK124", "SI-GA-G6"), JK124_reg2_org_12mix = mat_dir("JK124", "SI-GA-H6"),
                             #JK125
                             JK125_reg1_tis_1.1 = mat_dir("JK125_run2", "SI-GA-A4"), JK125_reg1_tis_1.2 = mat_dir("JK125_run2", "SI-GA-A3"), 
                             JK125_reg1_org_1 = mat_dir("JK125_run2", "SI-GA-B3"), JK125_reg1_cell_1 = mat_dir("JK125_run2", "SI-GA-C3"), 
                             JK125_reg2_tis_1 = mat_dir("JK125_run2", "SI-GA-D3"), JK125_reg2_org_1 = mat_dir("JK125_run2", "SI-GA-E3"), 
                             JK125_reg2_org_2 = mat_dir("JK125_run2", "SI-GA-F3"), JK125_reg2_cell_1 = mat_dir("JK125_run2", "SI-GA-G3"), 
                             JK125_reg1_org_2_r1 = mat_dir("JK125_run1", "SI-GA-D1"), JK125_reg1_org_3_r1 = mat_dir("JK125_run1", "SI-GA-E1"), 
                             JK125_reg2_tis_2_r1 = mat_dir("JK125_run1", "SI-GA-B1"), JK125_reg2_org_3_r1 = mat_dir("JK125_run1", "SI-GA-F1"), 
                             JK125_reg2_org_4_r1 = mat_dir("JK125_run1", "SI-GA-G1"), JK125_reg2_org_5_r1 = mat_dir("JK125_run1", "SI-GA-H1"),
                             #JK126
                             JK126_reg1_tis_1.1 = mat_dir("JK126", "SI-GA-A5"), JK126_reg1_tis_1.2 = mat_dir("JK126", "SI-GA-B5"), 
                             JK126_reg1_org_2 = mat_dir("JK126", "SI-GA-C5"), JK126_reg1_org_13mix = mat_dir("JK126", "SI-GA-D5"),
                             JK126_reg2_tis_1 = mat_dir("JK126", "SI-GA-E5"), JK126_reg2_org_1 = mat_dir("JK126", "SI-GA-F5"), 
                             JK126_reg2_org_2 = mat_dir("JK126", "SI-GA-G5"), JK126_reg2_org_3 = mat_dir("JK126", "SI-GA-H5"), 
                             #JK134
                             JK134_reg1_tis_1 = mat_dir("JK134", "SI-GA-A2"), JK134_reg1_org_2 = mat_dir("JK134", "SI-GA-D2"), 
                             JK134_reg1_org_3 = mat_dir("JK134", "SI-GA-E2"), JK134_reg2_tis_1 = mat_dir("JK134", "SI-GA-B2"), 
                             JK134_reg2_org_1 = mat_dir("JK134", "SI-GA-F2"), JK134_reg2_org_2 = mat_dir("JK134", "SI-GA-G2"), 
                             JK134_reg2_org_3 = mat_dir("JK134", "SI-GA-H2"),
                             #JK136 (JK202 recurrence)
                             JK136_reg1_tis_1 = mat_dir("JK136", "SI-GA-A7"), JK136_reg1_org_1 = mat_dir("JK136", "SI-GA-B7"), 
                             JK136_reg1_org_2 = mat_dir("JK136", "SI-GA-C7"), JK136_reg1_org_3 = mat_dir("JK136", "SI-GA-D7"),
                             JK136_reg2_tis_1 = mat_dir("JK136", "SI-GA-E7"), JK136_reg2_tis_2_br = mat_dir("batch_run", "SI-GA-A2"),
                             JK136_reg2_org_1 = mat_dir("JK136", "SI-GA-F7"), JK136_reg2_org_2 = mat_dir("JK136", "SI-GA-G7"), 
                             JK136_reg2_org_3 = mat_dir("JK136", "SI-GA-H7"), JK136_reg2_org_4_br = mat_dir("batch_run", "SI-GA-H2"),
                             JK202_reg1_tis_1 = mat_dir("JK202", "SI-GA-E1"), JK202_reg1_org_1 = mat_dir("JK202", "SI-GA-F1"),
                             JK202_reg1_org_2 = mat_dir("JK202", "SI-GA-G1"), JK202_reg1_org_3 = mat_dir("JK202", "SI-GA-H1"),
                             #JK142 (JK196 recurrence)
                             JK142_reg1_tis_1 = mat_dir("JK142", "SI-GA-A8"), JK142_reg1_org_1 = mat_dir("JK142", "SI-GA-B8"), 
                             JK142_reg1_org_2 = mat_dir("JK142", "SI-GA-C8"), JK142_reg1_cell_1 = mat_dir("JK142", "SI-GA-D8"),
                             JK142_reg2_tis_1 = mat_dir("JK142", "SI-GA-E8"), JK142_reg2_tis_2.1_br = mat_dir("batch_run", "SI-GA-B2"),
                             JK142_reg2_tis_2.2_br = mat_dir("batch_run", "SI-GA-C2"), JK142_reg2_org_1 = mat_dir("JK142", "SI-GA-F8"), 
                             JK142_reg2_org_2 = mat_dir("JK142", "SI-GA-G8"), JK142_reg2_cell_1 = mat_dir("JK142", "SI-GA-H8"),
                             JK196_reg1_tis_1 = mat_dir("JK196", "SI-GA-A1"), JK196_reg1_tis_2.1_br = mat_dir("batch_run", "SI-GA-E2"),
                             JK196_reg1_tis_2.2_br = mat_dir("batch_run", "SI-GA-F2"), JK196_reg1_org_1 = mat_dir("JK196", "SI-GA-B1"),
                             JK196_reg1_org_2 = mat_dir("JK196", "SI-GA-C1"), JK196_reg1_org_3 = mat_dir("JK196", "SI-GA-D1"),
                             #JK152
                             JK152_reg1_tis_1 = mat_dir("JK152", "SI-GA-A10"), JK152_reg1_org_1 = mat_dir("JK152", "SI-GA-B10"), 
                             JK152_reg1_org_2 = mat_dir("JK152", "SI-GA-C10"), JK152_reg1_cell_1 = mat_dir("JK152", "SI-GA-D10"),
                             JK152_reg2_tis_1 = mat_dir("JK152", "SI-GA-E10"), JK152_reg2_org_1 = mat_dir("JK152", "SI-GA-F10"), 
                             JK152_reg2_org_2 = mat_dir("JK152", "SI-GA-G10"), JK152_reg2_cell_1 = mat_dir("JK152", "SI-GA-H10"),
                             #JK153
                             JK153_reg1_tis_1 = mat_dir("JK153", "SI-GA-A9"), JK153_reg1_org_1 = mat_dir("JK153", "SI-GA-B9"), 
                             JK153_reg1_org_2 = mat_dir("JK153", "SI-GA-C9"), JK153_reg1_cell_1 = mat_dir("JK153", "SI-GA-D9"),
                             JK153_reg2_tis_1 = mat_dir("JK153", "SI-GA-E9"), JK153_reg2_org_1 = mat_dir("JK153", "SI-GA-F9"), 
                             JK153_reg2_org_2 = mat_dir("JK153", "SI-GA-G9"), JK153_reg2_org_3 = mat_dir("JK153", "SI-GA-H9"),
                             #JK156
                             JK156_reg1_tis_1 = mat_dir("JK156", "SI-GA-A11"), JK156_reg1_org_1 = mat_dir("JK156", "SI-GA-B11"), 
                             JK156_reg1_org_2 = mat_dir("JK156", "SI-GA-C11"), JK156_reg1_org_3 = mat_dir("JK156", "SI-GA-D11"),
                             JK156_reg2_tis_1 = mat_dir("JK156", "SI-GA-E11"), JK156_reg2_tis_2_br = mat_dir("batch_run", "SI-GA-D2"),
                             JK156_reg2_org_1 = mat_dir("JK156", "SI-GA-F11"), JK156_reg2_org_3 = mat_dir("JK156", "SI-GA-G11"),
                             #JK163
                             JK163_reg1_tis_1 = mat_dir("JK163", "SI-GA-A12"), JK163_reg1_org_1 = mat_dir("JK163", "SI-GA-B12"),
                             JK163_reg1_org_2 = mat_dir("JK163", "SI-GA-C12"), JK163_reg1_cell_1 = mat_dir("JK163", "SI-GA-D12"),
                             JK163_reg2_tis_1 = mat_dir("JK163", "SI-GA-E12"), JK163_reg2_org_1 = mat_dir("JK163", "SI-GA-F12"),
                             JK163_reg2_org_2 = mat_dir("JK163", "SI-GA-G12"), JK163_reg2_cell_1 = mat_dir("JK163", "SI-GA-H12")),
                    samp_levels = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg1_org_1", "JK124_reg1_org_2.1", "JK124_reg1_org_2.2",
                                    "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK124_reg2_org_12mix",
                                    "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg1_org_1", 
                                    "JK125_reg1_org_2_r1", "JK125_reg1_org_3_r1", "JK125_reg1_cell_1", 
                                    "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", 
                                    "JK125_reg2_org_1", "JK125_reg2_org_2",
                                    "JK125_reg2_org_3_r1", "JK125_reg2_org_4_r1", "JK125_reg2_org_5_r1", 
                                    "JK125_reg2_cell_1", 
                                    "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg1_org_2", "JK126_reg1_org_13mix",
                                    "JK126_reg2_tis_1", "JK126_reg2_org_1", "JK126_reg2_org_2", "JK126_reg2_org_3", 
                                    "JK134_reg1_tis_1", "JK134_reg1_org_2", "JK134_reg1_org_3", 
                                    "JK134_reg2_tis_1", "JK134_reg2_org_1", "JK134_reg2_org_2", "JK134_reg2_org_3",
                                    "JK136_reg1_tis_1", "JK136_reg1_org_1", "JK136_reg1_org_2", "JK136_reg1_org_3",
                                    "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", 
                                    "JK136_reg2_org_1", "JK136_reg2_org_2", "JK136_reg2_org_3", "JK136_reg2_org_4_br",
                                    "JK202_reg1_tis_1", "JK202_reg1_org_1", "JK202_reg1_org_2", "JK202_reg1_org_3",
                                    "JK142_reg1_tis_1", "JK142_reg1_org_1", "JK142_reg1_org_2", "JK142_reg1_cell_1",
                                    "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", 
                                    "JK142_reg2_org_1", "JK142_reg2_org_2", "JK142_reg2_cell_1",
                                    "JK196_reg1_tis_1", "JK196_reg1_tis_2.1_br", "JK196_reg1_tis_2.2_br", 
                                    "JK196_reg1_org_1", "JK196_reg1_org_2", "JK196_reg1_org_3",
                                    "JK152_reg1_tis_1", "JK152_reg1_org_1", "JK152_reg1_org_2", "JK152_reg1_cell_1",
                                    "JK152_reg2_tis_1", "JK152_reg2_org_1", "JK152_reg2_org_2", "JK152_reg2_cell_1",
                                    "JK153_reg1_tis_1", "JK153_reg1_org_1", "JK153_reg1_org_2", "JK153_reg1_cell_1",
                                    "JK153_reg2_tis_1", "JK153_reg2_org_1", "JK153_reg2_org_2", "JK153_reg2_org_3",
                                    "JK156_reg1_tis_1", "JK156_reg1_org_1", "JK156_reg1_org_2", "JK156_reg1_org_3",
                                    "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK156_reg2_org_1", "JK156_reg2_org_3",
                                    "JK163_reg1_tis_1", "JK163_reg1_org_1", "JK163_reg1_org_2", "JK163_reg1_cell_1",
                                    "JK163_reg2_tis_1", "JK163_reg2_org_1", "JK163_reg2_org_2", "JK163_reg2_cell_1"),
                    tum_levels = c("JK124", "JK125", "JK126", "JK134", "JK136", "JK202",
                                   "JK142", "JK196", "JK152", "JK153", "JK156", "JK163"))


# ---------------------------------------------------------------
## QC filtering
# ---------------------------------------------------------------

#Identify outliers based on read counts
counts_drop <- isOutlier(all$total_counts, nmads = 3, log = TRUE, type = "lower", batch = all$sample) 

#Identify outliers based on gene counts
genes_drop <- isOutlier(all$total_features_by_counts, nmads = 3, log = TRUE, type = "lower", batch = all$sample) 

#Identify outliers based on proportion of reads in mitochondrial genes
mito_drop <- isOutlier(all$pct_counts_mito, nmads = 3, type = "higher", batch = all$sample) 

#Remove outliers
all_cell_fil <- all[, !(genes_drop | counts_drop | mito_drop)]
all_cell_fil <- calculateQCMetrics(all_cell_fil, BPPARAM = MulticoreParam(5))

#Normalize expression
all_cell_gene_fil <- make.qcd.dataset(all_cell_fil, cell_min_genes = NULL, cell_min_UMI = NULL, cell_max_MT = NULL, 
                                      gene_min_cells = 10, gene_min_UMI = 20, n_cores = 10, rand_seed = 12345)

#Calculate doublet scores and identify outliers
dblt_scores <- bplapply(levels(all_cell_gene_fill$sample), function(samp){
  dat <- all_cell_gene_fil[, all_cell_gene_fil$sample == samp]
  scores <- doubletCells(dat, k = 30)
  filter <- isOutlier(scores, nmads = 3, type = "higher")
  return(data.frame(score = scores, filter = filter, row.names = colnames(dat)))
}, BPPARAM = MulticoreParam(20))

dblt_scores <- do.call(rbind.data.frame, dblt_scores)
dblt_scores <- dblt_scores[colnames(all_cell_gene_fil),]
all_cell_gene_fil$dblt_score <- dblt_scores$score
all_cell_gene_fil$dblt_fil <- dblt_scores$filter

#Remove doublet outliers and re-normalize expression data
all_final <- all_cell_gene_fil[,!all_cell_gene_fil$dblt_fil]
all_final <- make.qcd.dataset(all_final, cell_min_genes = NULL, cell_min_UMI = NULL, cell_max_MT = NULL, 
                              gene_min_cells = 10, gene_min_UMI = 20, n_cores = 10, rand_seed = 12345)







##-----------------------------------------------
#Function to load raw dataset if .rds file can be found, otherwise make it and save it
##-----------------------------------------------

load.raw.dat <- function(id, dirs, samp_levels, tum_levels, call_method = c("none", "emptyDrops", "cellRanger"), 
                         p_val = 0.01, n_cores = 4, rand_seed = 1, ...) {
  require(doRNG)
  
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
    #clusterEvalQ(cl, .libPaths("/home/vleblanc/R/x86_64-pc-linux-gnu-library_gcc-7.2.0/3.5"))
    doParallel::registerDoParallel(cl, cores = n_cores)
    
    all_samples_indiv <- foreach(i = 1:length(dirs), .packages = c('DropletUtils', 'scater')) %dorng% {
      #Get sample info
      samp <- dirs[i]
      tum <- unlist(strsplit(samp, "/"))[4]
      tum <- unlist(strsplit(tum, "_"))[1]
      
      #If .rds file already exists, load it
      if(call_method == "none") {
        out_file = paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/raw_cells.rds")
      }
      
      if(call_method == "emptyDrops") {
        out_file = paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/raw_cells_emptyDrops.rds")
      }
      
      if(call_method == "cellRanger") {
        out_file = paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/raw_cells_cellRanger.rds")
      }
      
      if(file.exists(out_file)) {
        raw_fil <- readRDS(out_file)
      } else {
        #Read raw count matrix
        raw <- DropletUtils::read10xCounts(unname(samp), col.names = TRUE)
        
        #keep_barcode <- nexprs(raw, byrow = FALSE) >= 1
        #raw <- raw[, keep_barcode]
        
        #Relabel rows with gene symbols (ensuring uniqueness for duplicated or missing symbols)
        rownames(raw) <- uniquifyFeatureNames(rowData(raw)$ID, rowData(raw)$Symbol)
        
        if(call_method == "none") {
          raw_fil <- raw
        } else {
          #Make and print knee plot
          bcrank <- barcodeRanks(counts(raw))
          uniq <- !duplicated(bcrank$rank) #only showing unique points for plotting speed
          
          knee_file <- paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/knee.pdf")
          
          pdf(knee_file)
          suppressMessages(plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2))
          
          abline(h=bcrank$inflection, col="darkgreen", lty=2)
          abline(h=bcrank$knee, col="dodgerblue", lty=2)
          
          legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
          suppressMessages(dev.off())
          
          if(!call_method %in% c("emptyDrops", "cellRanger")) {
            stop("call_method needs to be one of \'emptyDrops\' or \'cellRanger\'")
          }
          
          if(call_method == "emptyDrops") {
            #Call cells and empty droplets
            cells_out <- emptyDrops(counts(raw), ...)
            
            #Check if any non-significant barcodes are true for limited (in which case more interations should be performed)
            if(table(sig = cells_out$FDR <= p_val, limited = cells_out$Limited)[1,2] > 0) {
              message("Warning: non-significant barcodes are TRUE for Limited")
            }
            
            #Print out info about cell numbers
            info_file <- paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/cell_info_emptyDrops.txt")
            cat(paste0(nrow(cells_out), " barcodes were considered in total. \n", 
                       sum(is.na(cells_out$FDR)), " barcodes were used to model ambient RNA. \n",
                       sum(cells_out$FDR <= p_val, na.rm = TRUE), " barcodes were identified as cells at a BH-adjusted p-value of ", p_val, ", ", 
                       sum(cells_out$FDR == 0, na.rm = TRUE), " of which passed the knee filter."),
                file = info_file)
            
            #Retain only barcodes identified as cells
            raw_fil <- raw[, which(cells_out$FDR <= p_val)]
          }
          
          if(call_method == "cellRanger") {
            cells_out <- defaultDrops(counts(raw), expected = 1000)
            
            #Print out info about cell numbers
            info_file <- paste0("./data/Robjects/scRNA/raw/", tum, "/", names(samp), "/cell_info_cellRanger.txt")
            cat(paste0(length(cells_out), " barcodes were considered in total. \n", 
                       sum(cells_out), " barcodes were identified as cells."),
                file = info_file)
            
            #Retain only barcodes identified as cells
            raw_fil <- raw[, which(cells_out)]
          }
        }
        
        saveRDS(raw_fil, file = out_file)
      }
      
      #Return counts matrix to combine into full dataset
      return(raw_fil)
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
    #raw_dat <- getBMFeatureAnnos(raw_dat, ids = rowData(raw_dat)$id, filters = "ensembl_gene_id",
    #                             attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene",
    #                                            "chromosome_name", "start_position", "end_position"),
    #                             dataset = "hsapiens_gene_ensembl", host = "mar2016.archive.ensembl.org")
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
##Make QCd dataset
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