# ---------------------------------------------------------------
## Scripts to get copy number clones from 10X CNV data
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_analysis.R")
source("./funcs_figs.R")


# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

# Subset full dataset to keep JK136 cells
JK136_cnv <- all_cnv[, grepl("JK136|JK202", all_cnv$sample)] #all_cnv object generated in load.R script
JK136_cnv <- JK136_cnv[, JK136_cnv$effective_reads_per_1Mbp >= 50]
colData(JK136_cnv) <- droplevels(colData(JK136_cnv))

# Get coarse-grained copy number profiles
JK136_cnv_crs <- get.coarse.cnv(JK136_cnv)

# Get UMAP
JK136_cnv_crs <- run.umap(JK136_cnv_crs, assay_use = "copy_number", use_red_dim = FALSE,
                          rows_use = intersect(which(!is.na(rowVars(as.matrix(assay(JK136_cnv_crs, "copy_number"))))), 
                                               which(rowVars(as.matrix(assay(JK136_cnv_crs, "copy_number"))) > 0)),
                          metric = "manhattan", k = 15)

# Identify clusters and add to object
JK136_cnv_crs_dbscan <- dbscan::hdbscan(JK136_cnv_crs@reducedDims$UMAP, minPts = 10)
JK136_cnv_crs$dbclust <- as.factor(JK136_cnv_crs_dbscan$cluster)


#Get median copy number at each bin for each cluster
JK136_cnv_crs_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK136_cnv_crs$dbclust), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK136_cnv_crs[, JK136_cnv_crs$dbclust == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK136_cnv_crs$dbclust))

#Remove unmappable bins
JK136_cnv_crs_clustmeds <- JK136_cnv_crs_clustmeds[rowSums(is.na(JK136_cnv_crs_clustmeds)) == 0,]

#Get overlap between medians of each pair of clusters
JK136_cnv_crs_clustmeds_pairs <- expand.grid(clust1 = colnames(JK136_cnv_crs_clustmeds),
                                             clust2 = colnames(JK136_cnv_crs_clustmeds))
JK136_cnv_crs_clustmeds_pairs <- JK136_cnv_crs_clustmeds_pairs[JK136_cnv_crs_clustmeds_pairs$clust1 != JK136_cnv_crs_clustmeds_pairs$clust2,]
JK136_cnv_crs_clustmeds_pairs <- setNames(as.data.frame(unique(t(apply(JK136_cnv_crs_clustmeds_pairs, 1, sort)))),
                                          c("clust1", "clust2"))

JK136_cnv_crs_clustmeds_pairs$overlap <- apply(JK136_cnv_crs_clustmeds_pairs, 1, function(pair) {
  sum(JK136_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK136_cnv_crs_clustmeds)), pair[1]] == 
        JK136_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK136_cnv_crs_clustmeds)), pair[2]]) / 
    nrow(JK136_cnv_crs_clustmeds) * 100
})

#Merge clusters that have >=80% overlap (lower than other tumours b/c of increased ploidy so lower resolution)
JK136_cnv_crs$dbclust_mrg <- as.factor(ifelse(JK136_cnv_crs$dbclust %in% c("1", "2", "3", "4", "5", "6", "7", "11", "12", "13", "14", "15", "17", 
                                                                           "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", 
                                                                           "29", "30", "31", "32", "33"), "a",
                                              ifelse(JK136_cnv_crs$dbclust %in% c("35", "36", "37"),  "b",
                                                     ifelse(JK136_cnv_crs$dbclust %in% c("39", "40", "41", "42", "43", "44", "45", 
                                                                                         "46", "47", "48"), "c",
                                                            ifelse(JK136_cnv_crs$dbclust == "8", "d",
                                                                   ifelse(JK136_cnv_crs$dbclust == "9", "e",
                                                                          ifelse(JK136_cnv_crs$dbclust == "10", "f", 
                                                                                 ifelse(JK136_cnv_crs$dbclust == "16", "g",
                                                                                        ifelse(JK136_cnv_crs$dbclust == "34", "h",
                                                                                               ifelse(JK136_cnv_crs$dbclust == "38", "i", "x"))))))))))

# Plot merged clusters to visually assign to malignant clones or other groups
plot.scDNA.cnv(JK136_cnv_crs[,order(JK136_cnv_crs$dbclust_mrg)], 
               assay_use = "copy_number",
               cell_annots = c("sample", "source", "dbclust_mrg", "is_noisy", "ploidy_confidence"),
               annot_cols = list(is_noisy = c("0" = "black", "1" = "red"),
                                 ploidy_confidence = c("-4" = brewer.pal(11, "PRGn")[2], "-3" = brewer.pal(11, "PRGn")[3], 
                                                       "-2" = brewer.pal(11, "PRGn")[4], "0" = brewer.pal(11, "PuOr")[3], 
                                                       "1" = brewer.pal(11, "PuOr")[4], "2" = brewer.pal(11, "PuOr")[5],
                                                       setNames(colorRampPalette(c(brewer.pal(11, "PRGn")[7], brewer.pal(11, "PRGn")[11]))(length(c(3:max(as.numeric(JK136_cnv_crs$ploidy_confidence))))), 
                                                                c(3:max(as.numeric(JK136_cnv_crs$ploidy_confidence))))),
                                 dbclust_mrg = setNames(iwanthue(length(unique(JK136_cnv_crs$dbclust_mrg))), 
                                                        unique(JK136_cnv_crs$dbclust_mrg)),
                                 sample = setNames(iwanthue(length(unique(JK136_cnv_crs$sample))), unique(JK136_cnv_crs$sample)),
                                 source = c("tis" = "coral3", "org" = "cadetblue"),
                                 chr = setNames(rep(c("grey30", "grey80"), 12), paste0("chr", c(1:22, "X", "Y")))),
               heat_type = "cn", cluster_cells = F,
               out_file = "./figs/copy_number/scDNA/JK136_cnv_coarse_dbclust_dbclust_mrgs.png", height = 22, width = 12,
               gaps_row = cumsum(table(JK136_cnv_crs$dbclust_mrg)))

# Assign clones
JK136_cnv_crs$clone <- ifelse(as.character(JK136_cnv_crs$dbclust_mrg) == "b", "non-malignant",
                              ifelse(as.character(JK136_cnv_crs$dbclust_mrg) %in% c("f", "g", "h"), "s-phase",
                                     ifelse(as.character(JK136_cnv_crs$dbclust_mrg) %in% c("d", "e", "i", "x"), "noisy", 
                                            as.character(JK136_cnv_crs$dbclust_mrg))))

JK136_cnv_crs$clone <- factor(JK136_cnv_crs$clone, 
                              levels = c(unique(JK136_cnv_crs$clone)[!unique(JK136_cnv_crs$clone) %in% c("non-malignant", "s-phase", "noisy")],
                                         "non-malignant", "s-phase", "noisy"))


#Get median copy number at each bin for merged clusters
JK136_cnv_crs_mrg_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK136_cnv_crs$dbclust_mrg), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK136_cnv_crs[, JK136_cnv_crs$dbclust_mrg == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK136_cnv_crs$dbclust_mrg))

# Plot median copy numbers across malignant clones (heatmap shown in Figure 2)
plot.clustmeds(JK136_cnv_crs_mrg_clustmeds[, c("a", "c")], JK136_cnv_crs, plot_type = "heatmap")

# Plot proportions of cells assigned to each malignant clone
ggplot(as.data.frame(colData(JK136_cnv_crs)[JK136_cnv_crs$clone %in% c("a", "c"),]), 
       aes(x = sample, fill = clone)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = setNames(iwanthue(length(levels(JK136_cnv_crs$clone))), levels(JK136_cnv_crs$clone))) +
  labs(x = NULL, y = "percent cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

# Subset full dataset to keep JK142 cells
JK142_cnv <- all_cnv[, grepl("JK142|JK196", all_cnv$sample)] #all_cnv object generated in load.R script
JK142_cnv <- JK142_cnv[, JK142_cnv$effective_reads_per_1Mbp >= 50]
colData(JK142_cnv) <- droplevels(colData(JK142_cnv))

# Get coarse-grained copy number profiles
JK142_cnv_crs <- get.coarse.cnv(JK142_cnv)

# Get UMAP
JK142_cnv_crs <- run.umap(JK142_cnv_crs, assay_use = "copy_number", use_red_dim = FALSE,
                          rows_use = intersect(which(!is.na(rowVars(as.matrix(assay(JK142_cnv_crs, "copy_number"))))), 
                                               which(rowVars(as.matrix(assay(JK142_cnv_crs, "copy_number"))) > 0)),
                          metric = "manhattan", k = 15)

# Identify clusters and add to object
JK142_cnv_crs_dbscan <- dbscan::hdbscan(JK142_cnv_crs@reducedDims$UMAP, minPts = 10)
JK142_cnv_crs$dbclust <- as.factor(JK142_cnv_crs_dbscan$cluster)


#Get median copy number at each bin for each cluster
JK142_cnv_crs_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK142_cnv_crs$dbclust), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK142_cnv_crs[, JK142_cnv_crs$dbclust == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK142_cnv_crs$dbclust))

#Remove unmappable bins
JK142_cnv_crs_clustmeds <- JK142_cnv_crs_clustmeds[rowSums(is.na(JK142_cnv_crs_clustmeds)) == 0,]

#Get overlap between medians of each pair of clusters
JK142_cnv_crs_clustmeds_pairs <- expand.grid(clust1 = colnames(JK142_cnv_crs_clustmeds),
                                             clust2 = colnames(JK142_cnv_crs_clustmeds))
JK142_cnv_crs_clustmeds_pairs <- JK142_cnv_crs_clustmeds_pairs[JK142_cnv_crs_clustmeds_pairs$clust1 != JK142_cnv_crs_clustmeds_pairs$clust2,]
JK142_cnv_crs_clustmeds_pairs <- setNames(as.data.frame(unique(t(apply(JK142_cnv_crs_clustmeds_pairs, 1, sort)))),
                                          c("clust1", "clust2"))

JK142_cnv_crs_clustmeds_pairs$overlap <- apply(JK142_cnv_crs_clustmeds_pairs, 1, function(pair) {
  sum(JK142_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK142_cnv_crs_clustmeds)), pair[1]] == 
        JK142_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK142_cnv_crs_clustmeds)), pair[2]]) / 
    nrow(JK142_cnv_crs_clustmeds) * 100
})

# Merge clusters that have >=90% overlap
JK142_cnv_crs$dbclust_mrg <- as.factor(ifelse(JK142_cnv_crs$dbclust %in% c("1", "2", "3", "4", "5", "6", "7", "8", "24", "25"), "a",
                                              ifelse(JK142_cnv_crs$dbclust %in% c("9", "10", "11", "12", "13"),  "b",
                                                     ifelse(JK142_cnv_crs$dbclust %in% c("22", "31", "32"), "c",
                                                            ifelse(JK142_cnv_crs$dbclust %in% c("27", "28"), "d",
                                                                   ifelse(JK142_cnv_crs$dbclust %in% c("29", "30"), "e",
                                                                          ifelse(JK142_cnv_crs$dbclust == "14", "f", 
                                                                                 ifelse(JK142_cnv_crs$dbclust == "15", "g",
                                                                                        ifelse(JK142_cnv_crs$dbclust == "16", "h",
                                                                                               ifelse(JK142_cnv_crs$dbclust == "17", "i",
                                                                                                      ifelse(JK142_cnv_crs$dbclust == "18", "j",
                                                                                                             ifelse(JK142_cnv_crs$dbclust == "19", "k", 
                                                                                                                    ifelse(JK142_cnv_crs$dbclust == "20", "l",
                                                                                                                           ifelse(JK142_cnv_crs$dbclust == "21", "m",
                                                                                                                                  ifelse(JK142_cnv_crs$dbclust == "23", "n",
                                                                                                                                         ifelse(JK142_cnv_crs$dbclust == "26", "o", "x"))))))))))))))))

# Plot merged clusters to visually assign to malignant clones or other groups
plot.scDNA.cnv(JK142_cnv_crs[,order(JK142_cnv_crs$dbclust_mrg)], 
               assay_use = "copy_number",
               cell_annots = c("sample", "source", "dbclust_mrg", "is_noisy", "ploidy_confidence"),
               annot_cols = list(is_noisy = c("0" = "black", "1" = "red"),
                                 ploidy_confidence = c("-4" = brewer.pal(11, "PRGn")[2], "-3" = brewer.pal(11, "PRGn")[3], 
                                                       "-2" = brewer.pal(11, "PRGn")[4], "0" = brewer.pal(11, "PuOr")[3], 
                                                       "1" = brewer.pal(11, "PuOr")[4], "2" = brewer.pal(11, "PuOr")[5],
                                                       setNames(colorRampPalette(c(brewer.pal(11, "PRGn")[7], brewer.pal(11, "PRGn")[11]))(length(c(3:max(as.numeric(JK142_cnv_crs$ploidy_confidence))))), 
                                                                c(3:max(as.numeric(JK142_cnv_crs$ploidy_confidence))))),
                                 dbclust_mrg = setNames(iwanthue(length(unique(JK142_cnv_crs$dbclust_mrg))), 
                                                        unique(JK142_cnv_crs$dbclust_mrg)),
                                 sample = setNames(iwanthue(length(unique(JK142_cnv_crs$sample))), unique(JK142_cnv_crs$sample)),
                                 source = c("tis" = "coral3", "org" = "cadetblue"),
                                 chr = setNames(rep(c("grey30", "grey80"), 12), paste0("chr", c(1:22, "X", "Y")))),
               heat_type = "cn", cluster_cells = F,
               out_file = "./figs/copy_number/scDNA/JK142_cnv_coarse_dbclust_dbclust_mrgs.png", height = 22, width = 12,
               gaps_row = cumsum(table(JK142_cnv_crs$dbclust_mrg)))

# Assign clones
JK142_cnv_crs$clone <- ifelse(as.character(JK142_cnv_crs$dbclust_mrg) == "c", "non-malignant",
                              ifelse(as.character(JK142_cnv_crs$dbclust_mrg) %in% c("f", "g"), "s-phase",
                                     ifelse(as.character(JK142_cnv_crs$dbclust_mrg) %in% c("h", "i", "l", "o", "x"), "noisy", 
                                            as.character(JK142_cnv_crs$dbclust_mrg))))

JK142_cnv_crs$clone <- factor(JK142_cnv_crs$clone, 
                              levels = c(unique(JK142_cnv_crs$clone)[!unique(JK142_cnv_crs$clone) %in% c("non-malignant", "s-phase", "noisy")],
                                         "non-malignant", "s-phase", "noisy"))


#Get median copy number at each bin for merged clusters
JK142_cnv_crs_mrg_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK142_cnv_crs$dbclust_mrg), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK142_cnv_crs[, JK142_cnv_crs$dbclust_mrg == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK142_cnv_crs$dbclust_mrg))

# Plot median copy numbers across malignant clones (heatmap shown in Figure 2)
plot.clustmeds(JK142_cnv_crs_mrg_clustmeds[ c("d", "e", "a", "m", "b", "j")], JK142_cnv_crs, plot_type = "heatmap")

# Plot proportions of cells assigned to each malignant clone
ggplot(as.data.frame(colData(JK142_cnv_crs)[!JK142_cnv_crs$clone %in% c("k", "n", "non-malignant", "s-phase", "noisy"),]), 
       aes(x = sample, fill = clone)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = setNames(iwanthue(length(levels(JK142_cnv_crs$clone))), levels(JK142_cnv_crs$clone))) +
  labs(x = NULL, y = "percent cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

# Subset full dataset to keep JK153 cells
JK153_cnv <- all_cnv[, grepl("JK153", all_cnv$sample)] #all_cnv object generated in load.R script
JK153_cnv <- JK153_cnv[, JK153_cnv$effective_reads_per_1Mbp >= 50]
colData(JK153_cnv) <- droplevels(colData(JK153_cnv))

# Get coarse-grained copy number profiles
JK153_cnv_crs <- get.coarse.cnv(JK153_cnv)

# Get UMAP
JK153_cnv_crs <- run.umap(JK153_cnv_crs, assay_use = "copy_number", use_red_dim = FALSE,
                          rows_use = intersect(which(!is.na(rowVars(as.matrix(assay(JK153_cnv_crs, "copy_number"))))), 
                                               which(rowVars(as.matrix(assay(JK153_cnv_crs, "copy_number"))) > 0)),
                          metric = "manhattan", k = 15)

# Identify clusters and add to object
JK153_cnv_crs_dbscan <- dbscan::hdbscan(JK153_cnv_crs@reducedDims$UMAP, minPts = 10)
JK153_cnv_crs$dbclust <- as.factor(JK153_cnv_crs_dbscan$cluster)


#Get median copy number at each bin for each cluster
JK153_cnv_crs_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK153_cnv_crs$dbclust), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK153_cnv_crs[, JK153_cnv_crs$dbclust == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK153_cnv_crs$dbclust))

#Remove unmappable bins
JK153_cnv_crs_clustmeds <- JK153_cnv_crs_clustmeds[rowSums(is.na(JK153_cnv_crs_clustmeds)) == 0,]

#Get overlap between medians of each pair of clusters
JK153_cnv_crs_clustmeds_pairs <- expand.grid(clust1 = colnames(JK153_cnv_crs_clustmeds),
                                             clust2 = colnames(JK153_cnv_crs_clustmeds))
JK153_cnv_crs_clustmeds_pairs <- JK153_cnv_crs_clustmeds_pairs[JK153_cnv_crs_clustmeds_pairs$clust1 != JK153_cnv_crs_clustmeds_pairs$clust2,]
JK153_cnv_crs_clustmeds_pairs <- setNames(as.data.frame(unique(t(apply(JK153_cnv_crs_clustmeds_pairs, 1, sort)))),
                                          c("clust1", "clust2"))

JK153_cnv_crs_clustmeds_pairs$overlap <- apply(JK153_cnv_crs_clustmeds_pairs, 1, function(pair) {
  sum(JK153_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK153_cnv_crs_clustmeds)), pair[1]] == 
        JK153_cnv_crs_clustmeds[!grepl("chrX|chrY", rownames(JK153_cnv_crs_clustmeds)), pair[2]]) / 
    nrow(JK153_cnv_crs_clustmeds) * 100
})

# Merge clusters that have >=90% overlap
JK153_cnv_crs$dbclust_mrg <- as.factor(ifelse(JK153_cnv_crs$dbclust %in% c("3", "4"), "a",
                                              ifelse(JK153_cnv_crs$dbclust %in% c("9", "10", "11", "12", "13", "14", "15", "16", "17",
                                                                                  "18"),  "b",
                                                     ifelse(JK153_cnv_crs$dbclust %in% c("24", "25", "26", "27", "28", "29", "30", "31",
                                                                                         "32", "33", "34", "35", "36"), "c",
                                                            ifelse(JK153_cnv_crs$dbclust == "1", "d",
                                                                   ifelse(JK153_cnv_crs$dbclust == "2", "e",
                                                                          ifelse(JK153_cnv_crs$dbclust == "5", "f", 
                                                                                 ifelse(JK153_cnv_crs$dbclust == "6", "g",
                                                                                        ifelse(JK153_cnv_crs$dbclust == "7", "h",
                                                                                               ifelse(JK153_cnv_crs$dbclust == "8", "i", 
                                                                                                      ifelse(JK153_cnv_crs$dbclust == "19", "j",
                                                                                                             ifelse(JK153_cnv_crs$dbclust == "20", "k",
                                                                                                                    ifelse(JK153_cnv_crs$dbclust == "21", "l",
                                                                                                                           ifelse(JK153_cnv_crs$dbclust == "22", "m",
                                                                                                                                  ifelse(JK153_cnv_crs$dbclust == "23", "n", "x")))))))))))))))

# Plot merged clusters to visually assign to malignant clones or other groups
plot.scDNA.cnv(JK153_cnv_crs[,order(JK153_cnv_crs$dbclust_mrg)], 
               assay_use = "copy_number",
               cell_annots = c("sample", "source", "dbclust_mrg", "is_noisy", "ploidy_confidence"),
               annot_cols = list(is_noisy = c("0" = "black", "1" = "red"),
                                 ploidy_confidence = c("-4" = brewer.pal(11, "PRGn")[2], "-3" = brewer.pal(11, "PRGn")[3], 
                                                       "-2" = brewer.pal(11, "PRGn")[4], "0" = brewer.pal(11, "PuOr")[3], 
                                                       "1" = brewer.pal(11, "PuOr")[4], "2" = brewer.pal(11, "PuOr")[5],
                                                       setNames(colorRampPalette(c(brewer.pal(11, "PRGn")[7], brewer.pal(11, "PRGn")[11]))(length(c(3:max(as.numeric(JK153_cnv_crs$ploidy_confidence))))), 
                                                                c(3:max(as.numeric(JK153_cnv_crs$ploidy_confidence))))),
                                 dbclust_mrg = setNames(iwanthue(length(unique(JK153_cnv_crs$dbclust_mrg))), 
                                                        unique(JK153_cnv_crs$dbclust_mrg)),
                                 sample = setNames(iwanthue(length(unique(JK153_cnv_crs$sample))), unique(JK153_cnv_crs$sample)),
                                 source = c("tis" = "coral3", "org" = "cadetblue"),
                                 chr = setNames(rep(c("grey30", "grey80"), 12), paste0("chr", c(1:22, "X", "Y")))),
               heat_type = "cn", cluster_cells = F,
               out_file = "./figs/copy_number/scDNA/JK153_cnv_coarse_dbclust_dbclust_mrgs.png", height = 22, width = 12,
               gaps_row = cumsum(table(JK153_cnv_crs$dbclust_mrg)))

# Assign clones
JK153_cnv_crs$clone <- ifelse(as.character(JK153_cnv_crs$dbclust_mrg) == "a", "non-malignant",
                              ifelse(as.character(JK153_cnv_crs$dbclust_mrg) %in% c("e", "i"), "s-phase",
                                     ifelse(as.character(JK153_cnv_crs$dbclust_mrg) %in% c("f", "g", "h", "j", "k", "l", "m", "n", "x"), "noisy", 
                                            as.character(JK153_cnv_crs$dbclust_mrg))))

JK153_cnv_crs$clone <- factor(JK153_cnv_crs$clone, 
                              levels = c(unique(JK153_cnv_crs$clone)[!unique(JK153_cnv_crs$clone) %in% c("non-malignant", "s-phase", "noisy")],
                                         "non-malignant", "s-phase", "noisy"))


#Get median copy number at each bin for merged clusters
JK153_cnv_crs_mrg_clustmeds <- setNames(do.call(cbind.data.frame, lapply(levels(JK153_cnv_crs$dbclust_mrg), function(clust) {
  #Get CNV matrix
  cnv_mat <- as.matrix(assay(JK153_cnv_crs[, JK153_cnv_crs$dbclust_mrg == clust], "copy_number"))
  
  #Get median copy number for each bin
  data.frame(rowMedians(cnv_mat, na.rm = T),
             row.names = rownames(cnv_mat))
})), levels(JK153_cnv_crs$dbclust_mrg))

# Plot median copy numbers across malignant clones (heatmap shown in Figure 2)
plot.clustmeds(JK153_cnv_crs_mrg_clustmeds[, c("b", "c", "d")], JK153_cnv_crs, plot_type = "heatmap")

# Plot proportions of cells assigned to each malignant clone
ggplot(as.data.frame(colData(JK153_cnv_crs)[JK153_cnv_crs$clone %in% c("b", "c", "d"),]), 
       aes(x = sample, fill = clone)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = setNames(iwanthue(length(levels(JK153_cnv_crs$clone))), levels(JK153_cnv_crs$clone))) +
  labs(x = NULL, y = "percent cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))