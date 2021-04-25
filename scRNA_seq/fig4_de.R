# ---------------------------------------------------------------
## Scripts to calculate scUniFrac distances, perform differential expression analyses, and plot figures in Figure 4 of the manuscript
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_figs")
source("./funcs_analysis")

library(ggplot2)
library(scater)
library(scran)


# ---------------------------------------------------------------
## Tissue cells
# ---------------------------------------------------------------

# Load tissue cells object
tis <- readRDS("./data/Robjects/scRNA/clean/tis_mito_dblt_fil_allsamps.rds") #generated in the fig3_cell_types.R script

# Get tree
tis_tree <- clust.tree(tis, method = "pca", weight_vals = attr(tis@reducedDims$PCA, "d")[1:tis_sig_pcs], 
                       pcs_use = tis_sig_pcs, group = "gclust")

# Get unique sample pairs
tis_reps <- expand.grid(samp1 = levels(tis$sample), samp2 = levels(tis$sample), stringsAsFactors = FALSE)
tis_reps <- tis_reps[tis_reps$samp1 != tis_reps$samp2,]

# Add replicate types
tis_reps$rep_type <- factor(apply(tis_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(all(grepl("JK136", samp1) & grepl("JK202", samp2)) | all(grepl("JK202", samp1) & grepl("JK136", samp2)) | 
     all(grepl("JK142", samp1) & grepl("JK196", samp2)) | all(grepl("JK196", samp1) & grepl("JK142", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "prim_rec_xchip" #primary/recurrent pairs (run in independent experiments [different 10X chips])
    } else {
      "prim_rec" #primary/recurrent pairs (run concurrently)
    }
  } else if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter" #samples from different patients
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "intra_xchip" #samples from different regions of the same tumour (run in independent experiments)
    } else {
      "intra" #samples from different regions of the same tumour (run concurrently)
    }
  } else if(all(grepl("[.]", samp1) & grepl("[.]", samp2))){
    "tech" #same cell suspension run on separate 10X chip wells
  } else if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
    "bio_xchip" #separate pieces of tissue from the same tumour region (run in independent experiments)
  } else {
    "bio" #separate pieces of tissue from the same tumour region (run concurrently)
  }
}), 
levels = c("same_samp", "tech", "bio", "bio_xchip", "intra", "intra_xchip", "prim_rec", "prim_rec_xchip", "inter"))

# Remove replicated pairs except for inter-patient pairs so that each patient can be plotted separately
rows_remove <- t(apply(tis_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & tis_reps$rep_type != "inter"

tis_reps <- tis_reps[!rows_remove,]


# Get scUniFrac distances between sample pairs
tis_samp_diss <- sc.unifrac.multi(tis, group_by = "sample", clusts_name = "gclust", tree = tis_tree, 
                                  perm_iters = 1000, n_cores = 10)

# Add distances to replicate data frame
tis_reps$distance <- apply(tis_reps, 1, function(pair) tis_samp_diss$distance[pair[1], pair[2]])

# Mark by unique tumour
tis_reps$samp1_tumour <- factor(sapply(strsplit(as.character(tis_reps$samp1), "_"), "[", 1), levels = levels(tis$tumour))





# ---------------------------------------------------------------
## PDO cells
# ---------------------------------------------------------------

# Load PDO cells object
org <- readRDS("./data/Robjects/scRNA/clean/org_mito_dblt_fil_allsamps.rds") #generated in the fig3_cell_types.R script

# Get tree
org_tree <- clust.tree(org, method = "pca", weight_vals = attr(org@reducedDims$PCA, "d")[1:org_sig_pcs], 
                       pcs_use = org_sig_pcs, group = "gclust")

# Get unique sample pairs
org_reps <- expand.grid(samp1 = levels(org$sample), samp2 = levels(org$sample), stringsAsFactors = FALSE)
org_reps <- org_reps[org_reps$samp1 != org_reps$samp2,]

# Add replicate types
org_reps$rep_type <- factor(apply(org_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(all(grepl("JK136", samp1) & grepl("JK202", samp2)) | all(grepl("JK202", samp1) & grepl("JK136", samp2)) | 
     all(grepl("JK142", samp1) & grepl("JK196", samp2)) | all(grepl("JK196", samp1) & grepl("JK142", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "prim_rec_xchip" #primary/recurrent pairs (run in independent experiments [different 10X chips])
    } else {
      "prim_rec" #primary/recurrent pairs (run concurrently)
    }
  } else if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter" #samples from different patients
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "intra_xchip" #samples from different regions of the same tumour (run in independent experiments)
    } else {
      "intra" #samples from different regions of the same tumour (run concurrently)
    }
  } else if(all(grepl("[.]", samp1) & grepl("[.]", samp2))){
    "tech" #same cell suspension run on separate 10X chip wells
  } else if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
    "bio_xchip" #separate pieces of orgsue from the same tumour region (run in independent experiments)
  } else {
    "bio" #separate pieces of orgsue from the same tumour region (run concurrently)
  }
}), 
levels = c("same_samp", "tech", "bio", "bio_xchip", "intra", "intra_xchip", "prim_rec", "prim_rec_xchip", "inter"))

# Remove replicated pairs except for inter-patient pairs so that each patient can be plotted separately
rows_remove <- t(apply(org_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & org_reps$rep_type != "inter"

org_reps <- org_reps[!rows_remove,]


# Get scUniFrac distances between sample pairs
org_samp_diss <- sc.unifrac.multi(org, group_by = "sample", clusts_name = "gclust", tree = org_tree, 
                                  perm_iters = 1000, n_cores = 10)

# Add distances to replicate data frame
org_reps$distance <- apply(org_reps, 1, function(pair) org_samp_diss$distance[pair[1], pair[2]])

# Mark by unique tumour
org_reps$samp1_tumour <- factor(sapply(strsplit(as.character(org_reps$samp1), "_"), "[", 1), levels = levels(org$tumour))





# ---------------------------------------------------------------
## BTICs
# ---------------------------------------------------------------

# Load BTICs object
btic <- readRDS("./data/Robjects/scRNA/clean/btic_mito_dblt_fil_allsamps.rds") #generated in the fig3_cell_types.R script

# Get tree
btic_tree <- clust.tree(btic, method = "pca", weight_vals = attr(btic@reducedDims$PCA, "d")[1:btic_sig_pcs], 
                       pcs_use = btic_sig_pcs, group = "gclust")

# Get unique sample pairs
btic_reps <- expand.grid(samp1 = levels(btic$sample), samp2 = levels(btic$sample), stringsAsFactors = FALSE)
btic_reps <- btic_reps[btic_reps$samp1 != btic_reps$samp2,]

# Add replicate types
btic_reps$rep_type <- factor(apply(btic_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter" #BTIC lines derived from samples from different patients
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    "intra" #BTIC lines derived from samples from different regions of the same tumour 
  }
}), levels = c("intra", "inter"))

# Remove replicated pairs except for inter-patient pairs so that each patient can be plotted separately
rows_remove <- t(apply(btic_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & btic_reps$rep_type != "inter"

btic_reps <- btic_reps[!rows_remove,]


# Get scUniFrac distances between sample pairs
btic_samp_diss <- sc.unifrac.multi(btic, group_by = "sample", clusts_name = "gclust", tree = btic_tree, 
                                  perm_iters = 1000, n_cores = 10)

# Add distances to replicate data frame
btic_reps$distance <- apply(btic_reps, 1, function(pair) btic_samp_diss$distance[pair[1], pair[2]])

# Mark by unique tumour
btic_reps$samp1_tumour <- factor(sapply(strsplit(as.character(btic_reps$samp1), "_"), "[", 1), levels = levels(btic$tumour))





# ---------------------------------------------------------------
## All cells
# ---------------------------------------------------------------

# Combined replicate tables from each sample type
all_reps <- rbind(data.frame(tis_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "tis"),
                  data.frame(org_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "org"),
                  data.frame(btic_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "btic"))

all_reps$samp1_patient <- droplevels(factor(ifelse(all_reps$samp1_tumour == "JK202", "JK136",
                                                   ifelse(all_reps$samp1_tumour == "JK196", "JK142", as.character(all_reps$samp1_tumour))),
                                            levels = levels(all_reps$samp1_tumour)))

# Plot (see Figure 4 and Supplementary Figure 7)
ggplot(all_reps, aes(x = rep_type, y = distance)) +
  geom_boxplot(aes(fill = source), outlier.shape = NA, position = position_dodge2(preserve = "single")) +
  geom_point(pch = 21, position = position_jitterdodge(), aes(group = source, fill = samp1_tumour)) +
  #stat_summary(fun.y = median, geom = "crossbar", width = 0.5, position = "dodge") +
  scale_fill_manual(values = c(tum_cols, source_cols)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))










# ---------------------------------------------------------------
## Differential expression analyses
# ---------------------------------------------------------------

# ---------------------------------------------------------------
## Tissue cells
# ---------------------------------------------------------------

# Get malignant cells
tis_mal <- tis[, which(tis$cell_type == "malignant")]
colData(tis_mal) <- droplevels(colData(tis_mal))
tis_mal <- tis_mal[Matrix::rowSums(logcounts(tis_mal) > 0) >= 10, ]
tis_mal <- scater::normalize(tis_mal)

# Calculate differential expression between regions of each tumour
tis_mal_reg_de_list <- setNames(lapply(levels(tis_mal$patient), function(pat){
  dat <- tis_mal[!rowData(tis_mal)$is_feature_control, tis_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  res <- as.data.frame(res[order(res$FDR),])
}), levels(tis_mal$patient))

# Remove JK142 (only 1 malignant cell in reg2)
tis_mal_reg_de_list <- tis_mal_reg_de_list[names(tis_mal_reg_de_list) != "JK142"]


# Calculate differential expression between primary/recurrent pairs
tis_mal_primrec_de_list <- setNames(lapply(c("JK136", "JK142"), function(pat){
  setNames(lapply(c("reg1", "reg2"), function(reg) {
    dat <- tis_mal[!rowData(tis_mal)$is_feature_control, tis_mal$patient == pat & tis_mal$reg_stg %in% c(reg, "rec")]
    colData(dat) <- droplevels(colData(dat))
    dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
    res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$reg_stg, BPPARAM = MulticoreParam(5))
    res <- res$statistics[[2]]
    res <- as.data.frame(res[order(res$FDR),])
  }), c("reg1", "reg2"))
}), c("JK136", "JK142"))


# ---------------------------------------------------------------
## PDO cells
# ---------------------------------------------------------------

# Get malignant cells
org_mal <- org[, which(org$cell_type == "malignant")]
colData(org_mal) <- droplevels(colData(org_mal))
org_mal <- org_mal[Matrix::rowSums(logcounts(org_mal) > 0) >= 10, ]
org_mal <- scater::normalize(org_mal)

# Calculate differential expression between regions of each tumour
org_mal_reg_de_list <- setNames(lapply(levels(org_mal$patient), function(pat){
  dat <- org_mal[!rowData(org_mal)$is_feature_control, org_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  res <- as.data.frame(res[order(res$FDR),])
}), levels(org_mal$patient))

#Remove JK124 (only 1 malignant cell in reg2)
org_mal_reg_de_list <- org_mal_reg_de_list[names(org_mal_reg_de_list) != "JK124"]


# Calculate differential expression between primary/recurrent pairs
org_mal_primrec_de_list <- setNames(lapply(c("JK136", "JK142"), function(pat){
  setNames(lapply(c("reg1", "reg2"), function(reg) {
    dat <- org_mal[!rowData(org_mal)$is_feature_control, org_mal$patient == pat & org_mal$reg_stg %in% c(reg, "rec")]
    colData(dat) <- droplevels(colData(dat))
    dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
    res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$reg_stg, BPPARAM = MulticoreParam(5))
    res <- res$statistics[[2]]
    res <- as.data.frame(res[order(res$FDR),])
  }), c("reg1", "reg2"))
}), c("JK136", "JK142"))



# ---------------------------------------------------------------
## BTICs
# ---------------------------------------------------------------

btic_mal <- btic

# Calculate differential expression between regions of each tumour
btic_mal_reg_de_list <- setNames(lapply(levels(btic_mal$patient), function(pat){
  dat <- btic_mal[!rowData(btic_mal)$is_feature_control, btic_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  res <- as.data.frame(res[order(res$FDR),])
}), levels(btic_mal$patient))

#Remove JK153 (only one region)
btic_mal_reg_de_list <- btic_mal_reg_de_list[names(btic_mal_reg_de_list) != "JK153"]



# ---------------------------------------------------------------
## All cells
# ---------------------------------------------------------------

# Combine regional DE results
all_mal_reg_de_res <- setNames(lapply(intersect(names(tis_mal_reg_de_list), names(org_mal_reg_de_list)),
                                      function(tum) {
                                        mrg <- merge(tis_mal_reg_de_list[[tum]], org_mal_reg_de_list[[tum]], by = 0)
                                        colnames(mrg) <- mgsub::mgsub(colnames(mrg), c(".x", ".y"), c("_tis", "_org"))
                                        
                                        if(tum %in% names(btic_mal_reg_de_list)) {
                                          mrg <- merge(mrg, btic_mal_reg_de_list[[tum]], by.x = "Row.names", by.y = 0)
                                          colnames(mrg)[8:10] <- paste0(colnames(mrg)[8:10], "_btic")
                                        } else {
                                          mrg <- cbind(mrg, data.frame(overlap_btic = NA, p.value_btic = NA, FDR_btic = NA))
                                        }
                                        
                                        mrg$tumour <- tum
                                        
                                        return(mrg)
                                      }),
                               intersect(names(tis_mal_reg_de_list), names(org_mal_reg_de_list)))

# Plot tissue/PDO cell overlaps
ggplot(do.call(rbind.data.frame, all_mal_reg_de_res), aes(x = overlap_tis, y = overlap_org)) +
  geom_point(colour = source_cols["org"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ tumour, ncol = 4)

# Plot tissue cell/BTIC overlaps
ggplot(do.call(rbind.data.frame, all_mal_reg_de_res), aes(x = overlap_tis, y = overlap_btic)) +
  geom_point(colour = source_cols["btic"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ tumour, ncol = 4)



# Combine primary/recurrent DE results
all_mal_primrec_de_res <- setNames(lapply(names(tis_mal_primrec_de_list), function(pat) {
  setNames(lapply(names(tis_mal_primrec_de_list[[pat]]), function(reg) {
    mrg <- merge(tis_mal_primrec_de_list[[pat]][[reg]], org_mal_primrec_de_list[[pat]][[reg]], by = 0)
    colnames(mrg) <- mgsub::mgsub(colnames(mrg), c(".x", ".y"), c("_tis", "_org"))
    
    mrg$patient <- pat
    mrg$region <- reg
    
    return(mrg)
  }), names(tis_mal_primrec_de_list[[pat]]))
}), names(tis_mal_primrec_de_list))

all_mal_primrec_de_res$JK142 <- all_mal_primrec_de_res$JK142[c("reg1")]

# Plot
ggplot(do.call(rbind.data.frame, list(do.call(rbind.data.frame, all_mal_primrec_de_res$JK136), all_mal_primrec_de_res$JK142$reg1)), 
       aes(x = overlap_tis, y = overlap_org)) +
  geom_point(colour = source_cols["org"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ patient + region)






# ---------------------------------------------------------------
## Enrichment plots
# ---------------------------------------------------------------

# ---------------------------------------------------------------
## Genes over-expressed in organoid cells compared to tissue cells
# ---------------------------------------------------------------

#Get multi-list Metascape enrichment results (performed on top 500 genes from each DE analysis)
all_mal_source_org_up_enrich <- list(group_res = read.csv("./results/gene_exp/all_mal_source_de_org_pertumreg_up_enrich/Enrichment_GO/_FINAL_GO.csv",
                                                          header = TRUE, stringsAsFactors = FALSE),
                                     indiv_res = read.csv("./results/gene_exp/all_mal_source_de_org_pertumreg_up_enrich/Enrichment_GO/GO_AllLists.csv",
                                                          header = TRUE, stringsAsFactors = FALSE))

#Add group ID to individual results and order by it
all_mal_source_org_up_enrich$indiv_res <- merge(all_mal_source_org_up_enrich$indiv_res, 
                                                all_mal_source_org_up_enrich$group_res[,c("GO", "GROUP_ID", "FirstInGroupByLogP")], 
                                                by = "GO")

all_mal_source_org_up_enrich$indiv_res <- all_mal_source_org_up_enrich$indiv_res[order(all_mal_source_org_up_enrich$indiv_res$GROUP_ID),]

#Add tumour/region info
all_mal_source_org_up_enrich$indiv_res$tumour <- factor(sapply(strsplit(all_mal_source_org_up_enrich$indiv_res$GeneList, "[.]"), "[", 1),
                                                        levels = levels(all$tumour))
all_mal_source_org_up_enrich$indiv_res$region <- factor(sapply(strsplit(all_mal_source_org_up_enrich$indiv_res$GeneList, "[.]"), "[", 2),
                                                        levels = levels(all$region))

all_mal_source_org_up_enrich$indiv_res$patient <- factor(ifelse(all_mal_source_org_up_enrich$indiv_res$tumour == "JK202", "JK136",
                                                                ifelse(all_mal_source_org_up_enrich$indiv_res$tumour == "JK196", "JK142",
                                                                       as.character(all_mal_source_org_up_enrich$indiv_res$tumour))),
                                                         levels = levels(tis$patient))
all_mal_source_org_up_enrich$indiv_res$reg_stg <- factor(ifelse(all_mal_source_org_up_enrich$indiv_res$tumour %in% c("JK202", "JK196"),
                                                                "rec", as.character(all_mal_source_org_up_enrich$indiv_res$region)),
                                                         levels = c("reg1", "reg2", "rec"))

#Add term ID to description and make into a factor for plotting
all_mal_source_org_up_enrich$indiv_res$term_anno <- paste0(all_mal_source_org_up_enrich$indiv_res$Description, " (",
                                                           all_mal_source_org_up_enrich$indiv_res$GO, ")")
all_mal_source_org_up_enrich$indiv_res$term_anno <- factor(all_mal_source_org_up_enrich$indiv_res$term_anno,
                                                           levels = rev(unique(all_mal_source_org_up_enrich$indiv_res$term_anno)))

# Plot 
ggplot(all_mal_source_org_up_enrich$indiv_res[all_mal_source_org_up_enrich$indiv_res$GROUP_ID %in% c(1:10) & 
                                                all_mal_source_org_up_enrich$indiv_res$FirstInGroupByLogP == 1,], 
       aes(x = patient, y = term_anno)) +
  geom_point(aes(size = -Log.q.value., colour = reg_stg), position = position_dodge2(width = 0.8, padding = 0.5)) +
  scale_colour_manual(values = reg_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# ---------------------------------------------------------------
## Genes under-expressed in organoid cells compared to tissue cells
# ---------------------------------------------------------------

#Get multi-list Metascape enrichment results (performed on top 500 genes from each DE analysis)
all_mal_source_org_down_enrich <- list(group_res = read.csv("./results/gene_exp/all_mal_source_de_org_pertumreg_up_enrich/Enrichment_GO/_FINAL_GO.csv",
                                                          header = TRUE, stringsAsFactors = FALSE),
                                     indiv_res = read.csv("./results/gene_exp/all_mal_source_de_org_pertumreg_up_enrich/Enrichment_GO/GO_AllLists.csv",
                                                          header = TRUE, stringsAsFactors = FALSE))

#Add group ID to individual results and order by it
all_mal_source_org_down_enrich$indiv_res <- merge(all_mal_source_org_down_enrich$indiv_res, 
                                                all_mal_source_org_down_enrich$group_res[,c("GO", "GROUP_ID", "FirstInGroupByLogP")], 
                                                by = "GO")

all_mal_source_org_down_enrich$indiv_res <- all_mal_source_org_down_enrich$indiv_res[order(all_mal_source_org_down_enrich$indiv_res$GROUP_ID),]

#Add tumour/region info
all_mal_source_org_down_enrich$indiv_res$tumour <- factor(sapply(strsplit(all_mal_source_org_down_enrich$indiv_res$GeneList, "[.]"), "[", 1),
                                                        levels = levels(all$tumour))
all_mal_source_org_down_enrich$indiv_res$region <- factor(sapply(strsplit(all_mal_source_org_down_enrich$indiv_res$GeneList, "[.]"), "[", 2),
                                                        levels = levels(all$region))

all_mal_source_org_down_enrich$indiv_res$patient <- factor(ifelse(all_mal_source_org_down_enrich$indiv_res$tumour == "JK202", "JK136",
                                                                ifelse(all_mal_source_org_down_enrich$indiv_res$tumour == "JK196", "JK142",
                                                                       as.character(all_mal_source_org_down_enrich$indiv_res$tumour))),
                                                         levels = levels(tis$patient))
all_mal_source_org_down_enrich$indiv_res$reg_stg <- factor(ifelse(all_mal_source_org_down_enrich$indiv_res$tumour %in% c("JK202", "JK196"),
                                                                "rec", as.character(all_mal_source_org_down_enrich$indiv_res$region)),
                                                         levels = c("reg1", "reg2", "rec"))

#Add term ID to description and make into a factor for plotting
all_mal_source_org_down_enrich$indiv_res$term_anno <- paste0(all_mal_source_org_down_enrich$indiv_res$Description, " (",
                                                           all_mal_source_org_down_enrich$indiv_res$GO, ")")
all_mal_source_org_down_enrich$indiv_res$term_anno <- factor(all_mal_source_org_down_enrich$indiv_res$term_anno,
                                                           levels = rev(unique(all_mal_source_org_down_enrich$indiv_res$term_anno)))

# Plot 
ggplot(all_mal_source_org_down_enrich$indiv_res[all_mal_source_org_down_enrich$indiv_res$GROUP_ID %in% c(1:10) & 
                                                all_mal_source_org_down_enrich$indiv_res$FirstInGroupByLogP == 1,], 
       aes(x = patient, y = term_anno)) +
  geom_point(aes(size = -Log.q.value., colour = reg_stg), position = position_dodge2(width = 0.8, padding = 0.5)) +
  scale_colour_manual(values = reg_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# ---------------------------------------------------------------
## Genes over-expressed in BTICs compared to tissue cells
# ---------------------------------------------------------------

#Get multi-list Metascape enrichment results (performed on top 500 genes from each DE analysis)
all_mal_source_btic_up_enrich <- list(group_res = read.csv("./results/gene_exp/all_mal_source_de_btic_pertumreg_up_enrich/Enrichment_GO/_FINAL_GO.csv",
                                                          header = TRUE, stringsAsFactors = FALSE),
                                     indiv_res = read.csv("./results/gene_exp/all_mal_source_de_btic_pertumreg_up_enrich/Enrichment_GO/GO_AllLists.csv",
                                                          header = TRUE, stringsAsFactors = FALSE))

#Add group ID to individual results and order by it
all_mal_source_btic_up_enrich$indiv_res <- merge(all_mal_source_btic_up_enrich$indiv_res, 
                                                all_mal_source_btic_up_enrich$group_res[,c("GO", "GROUP_ID", "FirstInGroupByLogP")], 
                                                by = "GO")

all_mal_source_btic_up_enrich$indiv_res <- all_mal_source_btic_up_enrich$indiv_res[order(all_mal_source_btic_up_enrich$indiv_res$GROUP_ID),]

#Add tumour/region info
all_mal_source_btic_up_enrich$indiv_res$tumour <- factor(sapply(strsplit(all_mal_source_btic_up_enrich$indiv_res$GeneList, "[.]"), "[", 1),
                                                        levels = levels(all$tumour))
all_mal_source_btic_up_enrich$indiv_res$region <- factor(sapply(strsplit(all_mal_source_btic_up_enrich$indiv_res$GeneList, "[.]"), "[", 2),
                                                        levels = levels(all$region))

all_mal_source_btic_up_enrich$indiv_res$patient <- factor(ifelse(all_mal_source_btic_up_enrich$indiv_res$tumour == "JK202", "JK136",
                                                                ifelse(all_mal_source_btic_up_enrich$indiv_res$tumour == "JK196", "JK142",
                                                                       as.character(all_mal_source_btic_up_enrich$indiv_res$tumour))),
                                                         levels = levels(tis$patient))
all_mal_source_btic_up_enrich$indiv_res$reg_stg <- factor(ifelse(all_mal_source_btic_up_enrich$indiv_res$tumour %in% c("JK202", "JK196"),
                                                                "rec", as.character(all_mal_source_btic_up_enrich$indiv_res$region)),
                                                         levels = c("reg1", "reg2", "rec"))

#Add term ID to description and make into a factor for plotting
all_mal_source_btic_up_enrich$indiv_res$term_anno <- paste0(all_mal_source_btic_up_enrich$indiv_res$Description, " (",
                                                           all_mal_source_btic_up_enrich$indiv_res$GO, ")")
all_mal_source_btic_up_enrich$indiv_res$term_anno <- factor(all_mal_source_btic_up_enrich$indiv_res$term_anno,
                                                           levels = rev(unique(all_mal_source_btic_up_enrich$indiv_res$term_anno)))

# Plot 
ggplot(all_mal_source_btic_up_enrich$indiv_res[all_mal_source_btic_up_enrich$indiv_res$GROUP_ID %in% c(1:10) & 
                                                all_mal_source_btic_up_enrich$indiv_res$FirstInGroupByLogP == 1,], 
       aes(x = patient, y = term_anno)) +
  geom_point(aes(size = -Log.q.value., colour = reg_stg), position = position_dodge2(width = 0.8, padding = 0.5)) +
  scale_colour_manual(values = reg_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# ---------------------------------------------------------------
## Genes under-expressed in BTICss compared to tissue cells
# ---------------------------------------------------------------

#Get multi-list Metascape enrichment results (performed on top 500 genes from each DE analysis)
all_mal_source_btic_down_enrich <- list(group_res = read.csv("./results/gene_exp/all_mal_source_de_btic_pertumreg_up_enrich/Enrichment_GO/_FINAL_GO.csv",
                                                            header = TRUE, stringsAsFactors = FALSE),
                                       indiv_res = read.csv("./results/gene_exp/all_mal_source_de_btic_pertumreg_up_enrich/Enrichment_GO/GO_AllLists.csv",
                                                            header = TRUE, stringsAsFactors = FALSE))

#Add group ID to individual results and order by it
all_mal_source_btic_down_enrich$indiv_res <- merge(all_mal_source_btic_down_enrich$indiv_res, 
                                                  all_mal_source_btic_down_enrich$group_res[,c("GO", "GROUP_ID", "FirstInGroupByLogP")], 
                                                  by = "GO")

all_mal_source_btic_down_enrich$indiv_res <- all_mal_source_btic_down_enrich$indiv_res[order(all_mal_source_btic_down_enrich$indiv_res$GROUP_ID),]

#Add tumour/region info
all_mal_source_btic_down_enrich$indiv_res$tumour <- factor(sapply(strsplit(all_mal_source_btic_down_enrich$indiv_res$GeneList, "[.]"), "[", 1),
                                                          levels = levels(all$tumour))
all_mal_source_btic_down_enrich$indiv_res$region <- factor(sapply(strsplit(all_mal_source_btic_down_enrich$indiv_res$GeneList, "[.]"), "[", 2),
                                                          levels = levels(all$region))

all_mal_source_btic_down_enrich$indiv_res$patient <- factor(ifelse(all_mal_source_btic_down_enrich$indiv_res$tumour == "JK202", "JK136",
                                                                  ifelse(all_mal_source_btic_down_enrich$indiv_res$tumour == "JK196", "JK142",
                                                                         as.character(all_mal_source_btic_down_enrich$indiv_res$tumour))),
                                                           levels = levels(tis$patient))
all_mal_source_btic_down_enrich$indiv_res$reg_stg <- factor(ifelse(all_mal_source_btic_down_enrich$indiv_res$tumour %in% c("JK202", "JK196"),
                                                                  "rec", as.character(all_mal_source_btic_down_enrich$indiv_res$region)),
                                                           levels = c("reg1", "reg2", "rec"))

#Add term ID to description and make into a factor for plotting
all_mal_source_btic_down_enrich$indiv_res$term_anno <- paste0(all_mal_source_btic_down_enrich$indiv_res$Description, " (",
                                                             all_mal_source_btic_down_enrich$indiv_res$GO, ")")
all_mal_source_btic_down_enrich$indiv_res$term_anno <- factor(all_mal_source_btic_down_enrich$indiv_res$term_anno,
                                                             levels = rev(unique(all_mal_source_btic_down_enrich$indiv_res$term_anno)))

# Plot 
ggplot(all_mal_source_btic_down_enrich$indiv_res[all_mal_source_btic_down_enrich$indiv_res$GROUP_ID %in% c(1:10) & 
                                                  all_mal_source_btic_down_enrich$indiv_res$FirstInGroupByLogP == 1,], 
       aes(x = patient, y = term_anno)) +
  geom_point(aes(size = -Log.q.value., colour = reg_stg), position = position_dodge2(width = 0.8, padding = 0.5)) +
  scale_colour_manual(values = reg_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))