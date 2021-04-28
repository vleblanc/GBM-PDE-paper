# ---------------------------------------------------------------
## Scripts to write barcode-containing csv files to run cellranger-dna bamslice to combine cells assigned to the same clone
## Author: Veronique LeBlanc
# ---------------------------------------------------------------


# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

JK136_cnv_crs <- readRDS("./data/Robjects/scDNA/JK136_cnv_crs.rds") # object from fig2_cn_clones.R

# Get barcode lists
JK136_clone_barcodes <- setNames(lapply(levels(JK136_cnv_crs$sample), function(samp) {
  setNames(lapply(levels(JK136_cnv_crs$clone)[!levels(JK136_cnv_crs$clone) %in% c("s-phase", "noisy")], function(clone) {
    #Get sample/clone-specific barcodes
    barcodes <- colnames(JK136_cnv_crs)[JK136_cnv_crs$clone == clone & JK136_cnv_crs$sample == samp]
    
    #Replace GEM groups to match with original BAM
    barcodes <- paste0(sapply(strsplit(barcodes, "-"), "[", 1), "-1")
  }), levels(JK136_cnv_crs$clone)[!levels(JK136_cnv_crs$clone) %in% c("s-phase", "noisy")])
}), levels(JK136_cnv_crs$sample))

# Write CSV files
lapply(names(JK136_clone_barcodes), function(samp) {
  clone_names <- names(JK136_clone_barcodes[[samp]])
  
  lapply(clone_names, function(clone) {
    write.table(JK136_clone_barcodes[[samp]][[clone]], file = paste0("./data/cellranger_dna/", samp, "/", clone, "_barcodes.csv"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  })
  
  csv_files <- data.frame(library_id = paste0("clone_", clone_names),
                          barcodes_csv = paste0("/projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/", samp, "/", clone_names, "_barcodes.csv"))
  write.csv(csv_files, file = paste0("./data/cellranger_dna/", samp, "/", "clones.csv"), row.names = FALSE, quote = FALSE)
})





# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

JK142_cnv_crs <- readRDS("./data/Robjects/scDNA/JK142_cnv_crs.rds") # object from fig2_cn_clones.R

# Get barcode lists
JK142_clone_barcodes <- setNames(lapply(levels(JK142_cnv_crs$sample), function(samp) {
  setNames(lapply(levels(JK142_cnv_crs$clone)[!levels(JK142_cnv_crs$clone) %in% c("s-phase", "noisy")], function(clone) {
    #Get sample/clone-specific barcodes
    barcodes <- colnames(JK142_cnv_crs)[JK142_cnv_crs$clone == clone & JK142_cnv_crs$sample == samp]
    
    #Replace GEM groups to match with original BAM
    barcodes <- paste0(sapply(strsplit(barcodes, "-"), "[", 1), "-1")
  }), levels(JK142_cnv_crs$clone)[!levels(JK142_cnv_crs$clone) %in% c("s-phase", "noisy")])
}), levels(JK142_cnv_crs$sample))

# Write CSV files
lapply(names(JK142_clone_barcodes), function(samp) {
  clone_names <- names(JK142_clone_barcodes[[samp]])
  
  lapply(clone_names, function(clone) {
    write.table(JK142_clone_barcodes[[samp]][[clone]], file = paste0("./data/cellranger_dna/", samp, "/", clone, "_barcodes.csv"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  })
  
  csv_files <- data.frame(library_id = paste0("clone_", clone_names),
                          barcodes_csv = paste0("/projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/", samp, "/", clone_names, "_barcodes.csv"))
  write.csv(csv_files, file = paste0("./data/cellranger_dna/", samp, "/", "clones.csv"), row.names = FALSE, quote = FALSE)
})







# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

JK153_cnv_crs <- readRDS("./data/Robjects/scDNA/JK153_cnv_crs.rds") # object from fig2_cn_clones.R

# Get barcode lists
JK153_clone_barcodes <- setNames(lapply(levels(JK153_cnv_crs$sample), function(samp) {
  setNames(lapply(levels(JK153_cnv_crs$clone)[!levels(JK153_cnv_crs$clone) %in% c("s-phase", "noisy")], function(clone) {
    #Get sample/clone-specific barcodes
    barcodes <- colnames(JK153_cnv_crs)[JK153_cnv_crs$clone == clone & JK153_cnv_crs$sample == samp]
    
    #Replace GEM groups to match with original BAM
    barcodes <- paste0(sapply(strsplit(barcodes, "-"), "[", 1), "-1")
  }), levels(JK153_cnv_crs$clone)[!levels(JK153_cnv_crs$clone) %in% c("s-phase", "noisy")])
}), levels(JK153_cnv_crs$sample))

# Write CSV files
lapply(names(JK153_clone_barcodes), function(samp) {
  clone_names <- names(JK153_clone_barcodes[[samp]])
  
  lapply(clone_names, function(clone) {
    write.table(JK153_clone_barcodes[[samp]][[clone]], file = paste0("./data/cellranger_dna/", samp, "/", clone, "_barcodes.csv"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  })
  
  csv_files <- data.frame(library_id = paste0("clone_", clone_names),
                          barcodes_csv = paste0("/projects/vleblanc_prj/GBM_organoids/data/cellranger_dna/", samp, "/", clone_names, "_barcodes.csv"))
  write.csv(csv_files, file = paste0("./data/cellranger_dna/", samp, "/", "clones.csv"), row.names = FALSE, quote = FALSE)
})