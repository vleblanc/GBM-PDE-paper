
tum_cols = c("JK124" = "#BC4765", "JK125" = "#83CCB3", "JK126" = "#493737", "JK134" = "#CEBC52", "JK136" = "#533476", "JK202" = "#7D40D3",
             "JK142" = "#59713C", "JK196" = "#7BD355", "JK152" = "#C55D32", "JK153" = "#7E90C4", "JK156" = "#CAA7A0", "JK163" = "#CD5CBD")

plot.red.dim(tis, dim_use = "UMAP", colour_by = "tumour", col_vals = tum_cols, plot_title = "tis")

plot.red.dim(tis, exprs_use = "logcounts", dim_use = "UMAP", colour_by = "VWF", col_by_exp = TRUE, plot_title = "tis")

singleR.draw.heatmap.custom(tis_snglr$singler[[1]]$SingleR.single.main, top.n = Inf, 
                            clusters = factor(tis_snglr$meta.data$orig.ident, levels = levels(tis$gclust)), 
                            normalize = F, order.by.clusters = T)



### FIGURE 4

tis_reps <- expand.grid(samp1 = levels(tis$sample), samp2 = levels(tis$sample), stringsAsFactors = FALSE)
tis_reps <- tis_reps[tis_reps$samp1 != tis_reps$samp2,]

tis_reps$rep_type <- factor(apply(tis_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(all(grepl("JK136", samp1) & grepl("JK202", samp2)) | all(grepl("JK202", samp1) & grepl("JK136", samp2)) | 
     all(grepl("JK142", samp1) & grepl("JK196", samp2)) | all(grepl("JK196", samp1) & grepl("JK142", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "prim_rec_xchip"
    } else {
      "prim_rec"
    }
  } else if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter"
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "intra_xchip"
    } else {
      "intra"
    }
  } else if(all(grepl("[.]", samp1) & grepl("[.]", samp2))){
    "tech"
  } else if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
    "bio_xchip"
  } else {
    "bio"
  }
}), 
levels = c("same_samp", "tech", "bio", "bio_xchip", "intra", "intra_xchip", "prim_rec", "prim_rec_xchip", "inter"))

rows_remove <- t(apply(tis_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & tis_reps$rep_type != "inter"

tis_reps <- tis_reps[!rows_remove,]

tis_samp_diss <- sc.unifrac.multi(tis, group_by = "sample", clusts_name = "gclust", tree = tis_tree, 
                                  perm_iters = 1000, n_cores = 10)

tis_reps$distance <- apply(tis_reps, 1, function(pair) tis_samp_diss$distance[pair[1], pair[2]])
tis_reps$samp1_tumour <- factor(sapply(strsplit(as.character(tis_reps$samp1), "_"), "[", 1), levels = levels(tis$tumour))


org_reps <- expand.grid(samp1 = levels(org$sample), samp2 = levels(org$sample), stringsAsFactors = FALSE)
org_reps <- org_reps[org_reps$samp1 != org_reps$samp2,]

org_reps$rep_type <- factor(apply(org_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(all(grepl("JK136", samp1) & grepl("JK202", samp2)) | all(grepl("JK202", samp1) & grepl("JK136", samp2)) | 
     all(grepl("JK142", samp1) & grepl("JK196", samp2)) | all(grepl("JK196", samp1) & grepl("JK142", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "prim_rec_xchip"
    } else {
      "prim_rec"
    }
  } else if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter"
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
      "intra_xchip"
    } else {
      "intra"
    }
  } else if(all(grepl("[.]", samp1) & grepl("[.]", samp2))){
    "tech"
  } else if(sum(grepl("_r1|_br", c(samp1, samp2))) == 1) {
    "bio_xchip"
  } else {
    "bio"
  }
}), 
levels = c("same_samp", "tech", "bio", "bio_xchip", "intra", "intra_xchip", "prim_rec", "prim_rec_xchip", "inter"))

rows_remove <- t(apply(org_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & org_reps$rep_type != "inter"

org_reps <- org_reps[!rows_remove,]

org_samp_diss <- sc.unifrac.multi(org, group_by = "sample", clusts_name = "gclust", tree = org_tree, 
                                  perm_iters = 1000, n_cores = 10)

org_reps$distance <- apply(org_reps, 1, function(pair) org_samp_diss$distance[pair[1], pair[2]])
org_reps$samp1_tumour <- factor(sapply(strsplit(org_reps$samp1, "_"), "[", 1), levels = levels(org$tumour))



btic_reps <- expand.grid(samp1 = levels(btic$sample), samp2 = levels(btic$sample), stringsAsFactors = FALSE)
btic_reps <- btic_reps[btic_reps$samp1 != btic_reps$samp2,]

btic_reps$rep_type <- factor(apply(btic_reps, 1, function(pair){
  samp1 <- as.character(pair[1])
  samp2 <- as.character(pair[2])
  
  if(strsplit(samp1, "_")[[1]][1] != strsplit(samp2, "_")[[1]][1]) {
    "inter"
  } else if(all(grepl("reg1", samp1), grepl("reg2", samp2)) |
            all(grepl("reg2", samp1), grepl("reg1", samp2))) {
    "intra"
  }
}), levels = c("intra", "inter"))

rows_remove <- t(apply(btic_reps[,c("samp1", "samp2")], 1, function(pair) {
  pair[order(pair)]
}))
rows_remove <- duplicated(rows_remove) & btic_reps$rep_type != "inter"

btic_reps <- btic_reps[!rows_remove,]

btic_samp_diss <- sc.unifrac.multi(btic, group_by = "sample", clusts_name = "gclust", tree = btic_tree, 
                                   perm_iters = 1000, n_cores = 10)

btic_reps$distance <- apply(btic_reps, 1, function(pair) btic_samp_diss$distance[pair[1], pair[2]])
btic_reps$samp1_tumour <- factor(sapply(strsplit(btic_reps$samp1, "_"), "[", 1), levels = levels(btic$tumour))



all_reps <- rbind(data.frame(tis_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "tis"),
                  data.frame(org_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "org"),
                  data.frame(btic_reps[,c("samp1", "samp2", "rep_type", "distance", "samp1_tumour")], source = "btic"))

all_reps$samp1_patient <- droplevels(factor(ifelse(all_reps$samp1_tumour == "JK202", "JK136",
                                                   ifelse(all_reps$samp1_tumour == "JK196", "JK142", as.character(all_reps$samp1_tumour))),
                                            levels = levels(all_reps$samp1_tumour)))


ggplot(all_reps, aes(x = rep_type, y = distance)) +
  geom_boxplot(aes(fill = source), outlier.shape = NA, position = position_dodge2(preserve = "single")) +
  geom_point(pch = 21, position = position_jitterdodge(), aes(group = source, fill = samp1_tumour)) +
  #stat_summary(fun.y = median, geom = "crossbar", width = 0.5, position = "dodge") +
  scale_fill_manual(values = c(tum_cols, source_cols)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



tis_mal <- tis[, which(tis$cell_type == "malignant")]
colData(tis_mal) <- droplevels(colData(tis_mal))
tis_mal <- tis_mal[Matrix::rowSums(logcounts(tis_mal) > 0) >= 10, ]
tis_mal <- scater::normalize(tis_mal)

tis_mal_reg_de_list <- setNames(lapply(levels(tis_mal$patient), function(pat){
  dat <- tis_mal[!rowData(tis_mal)$is_feature_control, tis_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  
  #Add entrez IDs
  #res$entrez_id <- gconvert(query = rownames(res), target = "ENTREZGENE_ACC", mthreshold = 1, filter_na = FALSE)$target
  res <- as.data.frame(res[order(res$FDR),])
}), levels(tis_mal$patient))

#Remove JK142 (only 1 malignant cell in reg2)
tis_mal_reg_de_list <- tis_mal_reg_de_list[names(tis_mal_reg_de_list) != "JK142"]


org_mal <- org[, which(org$cell_type == "malignant")]
colData(org_mal) <- droplevels(colData(org_mal))
org_mal <- org_mal[Matrix::rowSums(logcounts(org_mal) > 0) >= 10, ]
org_mal <- scater::normalize(org_mal)

org_mal_reg_de_list <- setNames(lapply(levels(org_mal$patient), function(pat){
  dat <- org_mal[!rowData(org_mal)$is_feature_control, org_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  res <- as.data.frame(res[order(res$FDR),])
}), levels(org_mal$patient))

#Remove JK124 (only 1 malignant cell in reg2)
org_mal_reg_de_list <- org_mal_reg_de_list[names(org_mal_reg_de_list) != "JK124"]


btic_mal <- btic

write.table(rownames(btic_mal)[!rowData(btic_mal)$is_feature_control], file = "./results/clustering/btic_mito_dblt_fil_allsamps/btic_mal_background.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

btic_mal_reg_de <- scran::overlapExprs(btic_mal[!rowData(btic_mal)$is_feature_control, ], 
                                       groups = btic_mal$region, 
                                       block = btic_mal$tumour, direction = "any")

btic_mal_reg_de <- as.data.frame(btic_mal_reg_de[[1]])
btic_mal_reg_de$adj.p.val <- p.adjust(btic_mal_reg_de$p.value)

WriteXLS::WriteXLS(btic_mal_reg_de, "./results/clustering/btic_mito_dblt_fil/btic_mal_reg_de.xlsx", row.names = TRUE)


btic_mal_reg_de_list <- setNames(lapply(levels(btic_mal$patient), function(pat){
  dat <- btic_mal[!rowData(btic_mal)$is_feature_control, btic_mal$tumour == pat]
  dat <- dat[Matrix::rowSums(logcounts(dat)) > 0, ]
  res <- scran::pairwiseWilcox(logcounts(dat), clusters = dat$region, BPPARAM = MulticoreParam(5))
  res <- res$statistics[[2]]
  res <- as.data.frame(res[order(res$FDR),])
}), levels(btic_mal$patient))

#Remove JK153 (only one region)
btic_mal_reg_de_list <- btic_mal_reg_de_list[names(btic_mal_reg_de_list) != "JK153"]



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

ggplot(do.call(rbind.data.frame, all_mal_reg_de_res), aes(x = overlap_tis, y = overlap_org)) +
  geom_point(colour = source_cols["org"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ tumour, ncol = 4)

ggplot(do.call(rbind.data.frame, all_mal_reg_de_res), aes(x = overlap_tis, y = overlap_btic)) +
  geom_point(colour = source_cols["btic"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ tumour, ncol = 4)



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

ggplot(do.call(rbind.data.frame, list(do.call(rbind.data.frame, all_mal_primrec_de_res$JK136), all_mal_primrec_de_res$JK142$reg1)), 
       aes(x = overlap_tis, y = overlap_org)) +
  geom_point(colour = source_cols["org"]) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme_bw() +
  facet_wrap(~ patient + region)




#Get multi-list Metascape enrichment results
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

ggplot(all_mal_source_org_up_enrich$indiv_res[all_mal_source_org_up_enrich$indiv_res$GROUP_ID %in% c(1:10) & 
                                                all_mal_source_org_up_enrich$indiv_res$FirstInGroupByLogP == 1,], 
       aes(x = patient, y = term_anno)) +
  #geom_tile(aes(fill = Log.q.value.)) +
  geom_point(aes(size = -Log.q.value., colour = reg_stg), position = position_dodge2(width = 0.8, padding = 0.5)) +
  #scale_alpha_manual(values = c(1, 0.8)) +
  scale_colour_manual(values = reg_cols) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





