# ---------------------------------------------------------------
## Scripts to identify cell types from scRNA-seq data
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

all_final <- readRDS("./data/Robjects/clean/all_190604_190802mitodbltfil_fil.rds") #generated in the pre_processing.R script

# ---------------------------------------------------------------
## Tissue cells
# ---------------------------------------------------------------

## CLUSTER CELLS

# Get tissue cells
tis <- subset.and.scale(full = all_final, subset_by = "source", name = "tis", gene_min_cells = 10, gene_min_umi = 20,
                        scale = FALSE, indiv_scale = FALSE)

# Run clustering analysis
run.full.cluster.analysis(tis, "tis_mito_dblt_fil_allsamps", block = tis$tumour, pca_iters = 100, batch_correct = FALSE, 
                          bio_thresh = 0.01, fdr_thresh = 0.05, k = 30, de_method = "wilcox",
                          save_hvgs = TRUE, get_doublet_scores = FALSE, block_markers = FALSE,
                          n_cores = 15, rand_seed = 12345)

# Load output from clustering analysis
tis <- readRDS("./data/Robjects/scRNA/clean/tis_mito_dblt_fil_allsamps.rds")

# Get number of significant PCs
tis_sig_pcs <- attr(tis@reducedDims$PCA, "npcs")

#Get tree
tis_tree <- clust.tree(tis, method = "pca", weight_vals = attr(tis@reducedDims$PCA, "d")[1:tis_sig_pcs], 
                       pcs_use = tis_sig_pcs, group = "gclust")

# Reorder clusters based on their location on the tree
tis$gclust <- factor(tis$gclust, levels = tis_tree$hclust$order)


## COMPARE TO REFERENCE CELL TYPES

# Compare to reference cell types using SingleR
tis_snglr <- CreateSinglerObject(counts = as.matrix(logcounts(tis)), annot = tis$gclust, clusters = tis$gclust, numCores = 10, project.name = "tis cells")

# Based on cluster similarities to reference cell types (see supplementary figure 3), assign expression-based cell types
colData(tis)$cell_type_exp <- factor(ifelse(tis$gclust %in% c("19", "4", "20"), "immune",
                                            ifelse(tis$gclust == "12", "fibroblast",
                                                   ifelse(tis$gclust == "14", "endothelial",
                                                          ifelse(tis$gclust == "17", "oligodendrocyte",
                                                                 ifelse(tis$gclust == "27", "neuron", "malignant"))))),
                                     levels = c("malignant", "fibroblast", "endothelial", "oligodendrocyte", "neuron", "immune"))


## INFER COPY NUMBER

# Load GTEx references (see gtex_ref.R)
gtex_ref <- readRDS("./data/gtex_ref_matched.rds")
gtex_annot <- readRDS("./data/gtex_ref_annot.rds")

# Make new count matrix with reference data
tis_ref_mat <- merge(as.data.frame(as.matrix(counts(tis))), gtex_ref, by = 0)
rownames(tis_ref_mat) <- tis_ref_mat$Row.names
tis_ref_mat <- tis_ref_mat[,2:ncol(tis_ref_mat)]

#Keep only genes that pass the infercnv cutoff (based on the single cells only, not the reference data)
tis_ref_mat <- tis_ref_mat[Matrix::rowMeans(counts(tis)[rownames(tis_ref_mat),]) > 0.1,]

# Get gene locations
tis_gene_locs <- as.data.frame(rowData(tis)[rownames(tis_ref_mat), c("chr", "start_pos", "end_pos")])
tis_gene_locs <- tis_gene_locs[!is.na(tis_gene_locs$chr),]
tis_gene_locs <- tis_gene_locs[order(tis_gene_locs$chr, tis_gene_locs$start_pos),]

# Get cell type annotations
tis_annot <- rbind(as.data.frame(colData(tis)[,"cell_type_exp", drop = FALSE]), 
                   gtex_annot[colnames(gtex_ref), "cell_type", drop = FALSE])

# Run inferCNV
tis_cnv <- infercnv::CreateInfercnvObject(raw_counts_matrix = tis_ref_mat, 
                                          annotations_file = tis_annot, gene_order_file = tis_gene_locs, chr_exclude = c("X", "Y", "MT"), 
                                          ref_group_names = c("Brain - Cerebellum",
                                                              "Brain - Caudate (basal ganglia)",
                                                              "Brain - Cortex",
                                                              "Brain - Nucleus accumbens (basal ganglia)",
                                                              "Brain - Cerebellar Hemisphere",
                                                              "Brain - Frontal Cortex (BA9)",
                                                              "Brain - Hippocampus"))

tis_cnv <- infercnv::run(tis_cnv, cutoff = 0.1, out_dir = "./results/clustering/tis_mito_dblt_fil_allsamps/infercnv_gtex_ref", 
                         cluster_by_groups = F, cluster_references = F, denoise = T, HMM = F, output_format = NA, num_threads = 10)

# Load tree
tis_cnv_tree <- ape::read.tree("./results/clustering/tis_mito_dblt_fil_allsamps/infercnv_normcell_ref/infercnv.observations_dendrogram.txt")
tis_cnv_tree <- ape::as.hclust.phylo(tis_cnv_tree)

# Cut tree and add cluster info (see supplementary figure 4)
tis_cell_ann <- data.frame(cnv_clust = as.factor(cutree(tis_cnv_tree, k = 2)))

# Add malignant cell type to tissue object
tis$cell_type_cnv <- factor(ifelse(tis_cell_ann[colnames(tis), "cnv_clust"] == "1", "non_malignant",
                                   "malignant"),
                            levels = c("malignant", "non_malignant"))

# Add cell type info based on both expression and inferred CNV (keeping cells that are consistently identified as malignant or non-malignant)
colData(tis)$cell_type <- factor(ifelse((tis$cell_type_exp != "malignant" & tis$cell_type_cnv == "non_malignant") | (tis$cell_type_exp == "malignant" & tis$cell_type_cnv == "malignant"), as.character(tis$cell_type_exp), NA),
                                 levels = c("malignant", "fibroblast", "endothelial", "oligodendrocyte", "neuron", "immune"))




# ---------------------------------------------------------------
## Organoid cells
# ---------------------------------------------------------------

## CLUSTER CELLS

# Get organoid cells
org <- subset.and.scale(full = all_final, subset_by = "source", name = "org", gene_min_cells = 10, gene_min_umi = 20,
                        scale = FALSE, indiv_scale = FALSE)

# Run clustering analysis
run.full.cluster.analysis(org, "org_mito_dblt_fil_allsamps", block = org$tumour, pca_iters = 100, batch_correct = FALSE, 
                          bio_thresh = 0.01, fdr_thresh = 0.05, k = 30, de_method = "wilcox",
                          save_hvgs = TRUE, get_doublet_scores = FALSE, block_markers = FALSE,
                          n_cores = 15, rand_seed = 12345)

# Load output from clustering analysis
org <- readRDS("./data/Robjects/scRNA/clean/org_mito_dblt_fil_allsamps.rds")

# Get number of significant PCs
org_sig_pcs <- attr(org@reducedDims$PCA, "npcs")

#Get tree
org_tree <- clust.tree(org, method = "pca", weight_vals = attr(org@reducedDims$PCA, "d")[1:org_sig_pcs], 
                       pcs_use = org_sig_pcs, group = "gclust")

# Reorder clusters based on their location on the tree
org$gclust <- factor(org$gclust, levels = org_tree$hclust$order)


## COMPARE TO REFERENCE CELL TYPES

# Compare to reference cell types using SingleR
org_snglr <- CreateSinglerObject(counts = as.matrix(logcounts(org)), annot = org$gclust, clusters = org$gclust, numCores = 10, project.name = "org cells")

# Based on cluster similarities to reference cell types (see supplementary figure 3), assign expression-based cell types
colData(org)$cell_type_exp <- factor(ifelse(org$gclust %in% c("21", "27"), "immune",
                                            ifelse(org$gclust == "29", "fibroblast",
                                                   ifelse(org$gclust == "33", "oligodendrocyte", "malignant"))),
                                     levels = c("malignant", "fibroblast", "endothelial", "oligodendrocyte", "neuron", "immune"))


## INFER COPY NUMBER

#Make count matrix with non-malignant tissue cells as a reference
org_tisnormref_mat <- merge(as.data.frame(as.matrix(counts(org))), 
                            as.data.frame(as.matrix(counts(tis)[,which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron"))])), by = 0)
rownames(org_tisnormref_mat) <- org_tisnormref_mat$Row.names
org_tisnormref_mat <- org_tisnormref_mat[,2:ncol(org_tisnormref_mat)]

# Get gene locations
org_gene_locs <- as.data.frame(rowData(org)[rownames(org_tisnormref_mat), c("chr", "start_pos", "end_pos")])
org_gene_locs <- org_gene_locs[!is.na(org_gene_locs$chr),]
org_gene_locs <- org_gene_locs[order(org_gene_locs$chr, org_gene_locs$start_pos),]

# Get cell type annotations
org_annot <- rbind(setNames(as.data.frame(colData(org)[,"cell_type_exp", drop = FALSE]), "cell_type"), 
                   data.frame(cell_type = rep("tis_normal", length(which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron")))),
                              row.names = colnames(tis)[which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron"))]))

# Run inferCNV
org_cnv <- infercnv::CreateInfercnvObject(raw_counts_matrix = org_tisnormref_mat, 
                                          annotations_file = org_annot, gene_order_file = org_gene_locs, chr_exclude = c("X", "Y", "MT"), 
                                          ref_group_names = c("tis_normal"))

org_cnv <- infercnv::run(org_cnv, cutoff = 0.1, out_dir = "./results/clustering/org_mito_dblt_fil_allsamps/infercnv_tisnorm_ref", 
                         cluster_by_groups = F, cluster_references = F, denoise = T, HMM = F, output_format = NA, num_threads = 10)

# Load tree
org_cnv_tree <- ape::read.tree("./results/clustering/org_mito_dblt_fil_allsamps/infercnv_tisnorm_ref/infercnv.observations_dendrogram.txt")
org_cnv_tree <- ape::as.hclust.phylo(org_cnv_tree)

# Cut tree and add cluster info (see supplementary figure 4)
org_cell_ann$cnv_clust <- factor(cutree(org_cnv_tree, h = 20)[rownames(org_cell_ann)])
org_cell_ann$cell_type_cnv <- ifelse(org_cell_ann$cnv_clust %in% c("170", "171", "172", "173", "174", "175", "201", "202", "211"),
                                     "non_malignant", "malignant")

# Add malignant status based on inferred CNVs to organoid object
org$cell_type_cnv <- factor(org_cell_ann[colnames(org), "cell_type_cnv"],
                            levels = c("malignant", "non_malignant"))

# Add cell type info based on both expression and inferred CNV (keeping cells that are consistently identified as malignant or non-malignant)
colData(org)$cell_type <- factor(ifelse((org$cell_type_exp != "malignant" & org$cell_type_cnv == "non_malignant") | (org$cell_type_exp == "malignant" & org$cell_type_cnv == "malignant"), as.character(org$cell_type_exp), NA),
                                 levels = c("malignant", "fibroblast", "endothelial", "oligodendrocyte", "neuron", "immune"))





# ---------------------------------------------------------------
## BTICs
# ---------------------------------------------------------------

## CLUSTER CELLS

# Get BTICs
btic <- subset.and.scale(full = all_final, subset_by = "source", name = "btic", gene_min_cells = 10, gene_min_umi = 20,
                        scale = FALSE, indiv_scale = FALSE)

# Run clustering analysis
run.full.cluster.analysis(btic, "btic_mito_dblt_fil_allsamps", block = btic$tumour, pca_iters = 100, batch_correct = FALSE, 
                          bio_thresh = 0.01, fdr_thresh = 0.05, k = 30, de_method = "wilcox",
                          save_hvgs = TRUE, get_doublet_scores = FALSE, block_markers = FALSE,
                          n_cores = 15, rand_seed = 12345)

# Load output from clustering analysis
btic <- readRDS("./data/Robjects/scRNA/clean/btic_mito_dblt_fil_allsamps.rds")

# Get number of significant PCs
btic_sig_pcs <- attr(btic@reducedDims$PCA, "npcs")

#Get tree
btic_tree <- clust.tree(btic, method = "pca", weight_vals = attr(btic@reducedDims$PCA, "d")[1:btic_sig_pcs], 
                       pcs_use = btic_sig_pcs, group = "gclust")

# Reorder clusters based on their location on the tree
btic$gclust <- factor(btic$gclust, levels = btic_tree$hclust$order)


## COMPARE TO REFERENCE CELL TYPES

# Compare to reference cell types using SingleR
btic_snglr <- CreateSinglerObject(counts = as.matrix(logcounts(btic)), annot = btic$gclust, clusters = btic$gclust, numCores = 10, project.name = "btic cells")


## INFER COPY NUMBER

#Make count matrix with non-malignant tissue cells as a reference
btic_tisnormref_mat <- merge(as.data.frame(as.matrix(counts(btic))), 
                            as.data.frame(as.matrix(counts(tis)[,which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron"))])), by = 0)
rownames(btic_tisnormref_mat) <- btic_tisnormref_mat$Row.names
btic_tisnormref_mat <- btic_tisnormref_mat[,2:ncol(btic_tisnormref_mat)]

# Get gene locations
btic_gene_locs <- as.data.frame(rowData(btic)[rownames(btic_tisnormref_mat), c("chr", "start_pos", "end_pos")])
btic_gene_locs <- btic_gene_locs[!is.na(btic_gene_locs$chr),]
btic_gene_locs <- btic_gene_locs[order(btic_gene_locs$chr, btic_gene_locs$start_pos),]

# Get cell type annotations
btic_annot <- rbind(setNames(as.data.frame(colData(btic)[,"cell_type_exp", drop = FALSE]), "cell_type"), 
                   data.frame(cell_type = rep("tis_normal", length(which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron")))),
                              row.names = colnames(tis)[which(tis$cell_type %in% c("fibroblast", "endothelial", "oligodendrocyte", "neuron"))]))

# Run inferCNV
btic_cnv <- infercnv::CreateInfercnvObject(raw_counts_matrix = btic_tisnormref_mat, 
                                          annotations_file = btic_annot, gene_order_file = btic_gene_locs, chr_exclude = c("X", "Y", "MT"), 
                                          ref_group_names = c("tis_normal"))

btic_cnv <- infercnv::run(btic_cnv, cutoff = 0.1, out_dir = "./results/clustering/btic_mito_dblt_fil_allsamps/infercnv_tisnorm_ref", 
                         cluster_by_groups = F, cluster_references = F, denoise = T, HMM = F, output_format = NA, num_threads = 10)


# Add cell type info (no non-malignant cells in the BTICs)
colData(btic)$cell_type <- factor("malignant", levels = c("malignant", "fibroblast", "endothelial", "oligodendrocyte", "neuron", "immune"))