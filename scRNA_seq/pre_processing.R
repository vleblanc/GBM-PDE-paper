# ---------------------------------------------------------------
## Scripts to load data into a combined object, QC filter cells and genes, and normalize expression
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

library(scater)
library(scran)

source("./funcs_load_clean.R")

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