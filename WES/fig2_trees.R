# ---------------------------------------------------------------
## Scripts to load CITUP and PyClone to plot trees (Figure 2)
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_analysis.R")

library(igraph)
library(rhdf5)
library(ggplot2)

# ---------------------------------------------------------------
## JK124
# ---------------------------------------------------------------

# Load CITUP results
JK124_citup_iter <- h5dump("./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/citup_iter/iter_results.h5")

# Get optimal solution
JK124_citup_iter_sol <- as.character(JK124_citup_iter$results$optimal$index[1])

# Get clone info from optimal solution
JK124_citup_iter_clones <- setNames(reshape2::melt(JK124_citup_iter$trees[[JK124_citup_iter_sol]]$clone_freq$block0_values),
                                    c("clone_id", "sample_id", "clonal_prev"))
JK124_citup_iter_clones$sample_id <- c("JK124_reg1_tis", "JK124_reg1_org", "JK124_reg2_tis")[JK124_citup_iter_clones$sample_id]
JK124_citup_iter_clones$clone_id <- JK124_citup_iter_clones$clone_id - 1

# Load PyClone variants (files generated in prep_pyclone.R)
JK124_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200/pyclone/JK124_reg1_tissue_titan_subclonefil_adfil_noreg2org.tsv"),
                            read.delim("./data/exome_seq/JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200/pyclone/JK124_reg1_organoid_titan_subclonefil_adfil_noreg2org.tsv"),
                            read.delim("./data/exome_seq/JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200/pyclone/JK124_reg2_tissue_titan_subclonefil_adfil_noreg2org.tsv"))
JK124_pyclone_vars <- JK124_pyclone_vars[!duplicated(JK124_pyclone_vars$mutation_id),]

# Load PyClone loci-level results
JK124_pyclone_loci <- read.delim("./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/loci.tsv", stringsAsFactors = FALSE)
JK124_pyclone_loci$sample_id <- factor(JK124_pyclone_loci$sample_id,
                                       levels = c("JK124_reg1_tis", "JK124_reg1_org", "JK124_reg2_tis"))

# Merge PyClone results with variants to get copy number
JK124_pyclone_loci <- merge(JK124_pyclone_loci, JK124_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK124_pyclone_loci$cluster_id <- as.factor(JK124_pyclone_loci$cluster_id)

# Load variants
JK124_vars <- get.var.info.from.table("./results/exome/common_mutations/JK124_noreg2org_strelka_mutect/tp-baseline.table")
JK124_vars$tbl_nosnp$mutation_id <- paste0(JK124_vars$tbl_nosnp$gene, "_", JK124_vars$tbl_nosnp$CHROM, ":", 
                                           as.character(JK124_vars$tbl_nosnp$POS))

JK124_pyclone_loci <- merge(JK124_pyclone_loci, JK124_vars$tbl_nosnp[,c("mutation_id", "TYPE", "gene", "effect", "impact", "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT")],
                            by = "mutation_id", all.x = TRUE)

# Save object
saveRDS(JK124_pyclone_loci, "./results/exome/common_mutations/JK124_pyclone_loci.rds")

#Add clone info
JK124_pyclone_loci$clone_id <- factor(setNames(JK124_citup_iter$trees[[JK124_citup_iter_sol]]$variant_assignment$values, 
                                               c(0,1,2)[JK124_citup_iter$trees[[JK124_citup_iter_sol]]$variant_assignment$index + 1])
                                      [as.character(JK124_pyclone_loci$cluster_id)])

# Get tree
JK124_tree <- graph_from_data_frame(data.frame(source = JK124_citup_iter$trees[[JK124_citup_iter_sol]]$adjacency_list$block0_values[1,],
                                               target = JK124_citup_iter$trees[[JK124_citup_iter_sol]]$adjacency_list$block0_values[2,]))

# Get number of variants per clone
JK124_clone_muts <- table(JK124_pyclone_loci$clone_id) / length(levels(JK124_pyclone_loci$sample_id))

# Build tree with branch lengths proportional to the number of variants
JK124_coords <- layout_(JK124_tree, as_tree())
JK124_coords[,2] <- c(0, -JK124_clone_muts[1], -sum(JK124_clone_muts[1:2]), -sum(JK124_clone_muts[c(1,3)]))
JK124_coords[,1] <- 0

# Plot tree
plot.igraph(JK124_tree, layout = JK124_coords)

# Get variant info for each clone
lapply(levels(JK124_pyclone_loci$clone_id), function(clone) {
  clone_info <- JK124_pyclone_loci[which(JK124_pyclone_loci$clone_id == clone & JK124_pyclone_loci$sample_id == "JK124_reg1_tis"), c(11:16)]
  
  message("genes with high impact variants")
  clone_info[clone_info$impact == "HIGH",]
  
  message("GBM genes")
  clone_info[clone_info$gene %in% brennan_genes,]
  
  message("# of variants")
  nrow(clone_info)
})

# Get proportions of malignant cells assigned to each clone
JK124_citup_iter_clone_freqs <- setNames(as.data.frame(JK124_citup_iter$trees[[JK124_citup_iter_sol]]$clone_freq$block0_values),
                                         levels(JK124_pyclone_loci$sample_id))
JK124_citup_iter_clone_freqs$clone <- as.factor(0:(nrow(JK124_citup_iter_clone_freqs) - 1))

JK124_citup_iter_clone_freqs <- melt(JK124_citup_iter_clone_freqs, id.var = "clone", variable.name = "sample", value.name = "clone_prev")

# Plot
ggplot(JK124_citup_iter_clone_freqs, aes(x = sample, y = clone_prev, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = unname(iwanthue(length(levels(JK124_citup_iter_clone_freqs$clone))))) +
  labs(x = NULL, y = "clonal prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

# Load CITUP results
JK136_citup_iter <- h5dump("./results/exome/pyclone/JK136_subclonefil_adfil/citup_iter/iter_results.h5")

# Get optimal solution
JK136_citup_iter_sol <- as.character(JK136_citup_iter$results$optimal$index[1])

# Get clone info from optimal solution
JK136_citup_iter_clones <- setNames(reshape2::melt(JK136_citup_iter$trees[[JK136_citup_iter_sol]]$clone_freq$block0_values),
                                    c("clone_id", "sample_id", "clonal_prev"))
JK136_citup_iter_clones$sample_id <- c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org", "JK202_tis", "JK202_org")[JK136_citup_iter_clones$sample_id]
JK136_citup_iter_clones$clone_id <- JK136_citup_iter_clones$clone_id - 1

# Load PyClone variants (files generated in prep_pyclone.R)
JK136_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201/pyclone/JK136_reg1_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201/pyclone/JK136_reg1_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201/pyclone/JK136_reg2_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201/pyclone/JK136_reg2_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK202_tissue/hg38_no_alt/EXC/B58180_B58201/pyclone/JK202_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK202_organoid/hg38_no_alt/EXC/B58181_B58201/pyclone/JK202_organoid_titan_subclonefil_adfil.tsv"))
JK136_pyclone_vars <- JK136_pyclone_vars[!duplicated(JK136_pyclone_vars$mutation_id),]

# Load PyClone loci-level results
JK136_pyclone_loci <- read.delim("./results/exome/pyclone/JK136_subclonefil_adfil/tables/loci.tsv", stringsAsFactors = FALSE)
JK136_pyclone_loci$sample_id <- factor(JK136_pyclone_loci$sample_id,
                                       levels = c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org",
                                                  "JK202_tis", "JK202_org"))

# Merge PyClone results with variants to get copy number
JK136_pyclone_loci <- merge(JK136_pyclone_loci, JK136_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK136_pyclone_loci$cluster_id <- as.factor(JK136_pyclone_loci$cluster_id)

# Load variants
JK136_vars <- get.var.info.from.table("./results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.table")
JK136_vars$tbl_nosnp$mutation_id <- paste0(JK136_vars$tbl_nosnp$gene, "_", JK136_vars$tbl_nosnp$CHROM, ":", 
                                           as.character(JK136_vars$tbl_nosnp$POS))

JK136_pyclone_loci <- merge(JK136_pyclone_loci, JK136_vars$tbl_nosnp[,c("mutation_id", "TYPE", "gene", "effect", "impact", "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT")],
                            by = "mutation_id", all.x = TRUE)

# Save object
saveRDS(JK136_pyclone_loci, "./results/exome/common_mutations/JK136_pyclone_loci.rds")

#Add clone info
JK136_pyclone_loci$clone_id <- factor(setNames(JK136_citup_iter$trees[[JK136_citup_iter_sol]]$variant_assignment$values, 
                                               c(0:4)[JK136_citup_iter$trees[[JK136_citup_iter_sol]]$variant_assignment$index + 1])
                                      [as.character(JK136_pyclone_loci$cluster_id)])

# Get tree
JK136_tree <- graph_from_data_frame(data.frame(source = JK136_citup_iter$trees[[JK136_citup_iter_sol]]$adjacency_list$block0_values[1,],
                                               target = JK136_citup_iter$trees[[JK136_citup_iter_sol]]$adjacency_list$block0_values[2,]))

# Get number of variants per clone
JK136_clone_muts <- table(JK136_pyclone_loci$clone_id) / length(levels(JK136_pyclone_loci$sample_id))

# Build tree with branch lengths proportional to the number of variants
JK136_coords <- layout_(JK136_tree, as_tree())
JK136_coords[,2] <- c(0, -JK136_clone_muts[1], -sum(JK136_clone_muts[1:2]), -sum(JK136_clone_muts[1:3]), -sum(JK136_clone_muts[c(1:4)]),
                      -sum(JK136_clone_muts[c(1,2,5)]))
JK136_coords[,1] <- 0

# Plot tree
plot.igraph(JK136_tree, layout = JK136_coords)

# Get variant info for each clone
lapply(levels(JK136_pyclone_loci$clone_id), function(clone) {
  clone_info <- JK136_pyclone_loci[which(JK136_pyclone_loci$clone_id == clone & JK136_pyclone_loci$sample_id == "JK136_reg1_tis"), c(11:16)]
  
  message("genes with high impact variants")
  clone_info[clone_info$impact == "HIGH",]
  
  message("GBM genes")
  clone_info[clone_info$gene %in% brennan_genes,]
  
  message("# of variants")
  nrow(clone_info)
})

# Get proportions of malignant cells assigned to each clone
JK136_citup_iter_clone_freqs <- setNames(as.data.frame(JK136_citup_iter$trees[[JK136_citup_iter_sol]]$clone_freq$block0_values),
                                         levels(JK136_pyclone_loci$sample_id))
JK136_citup_iter_clone_freqs$clone <- as.factor(0:(nrow(JK136_citup_iter_clone_freqs) - 1))

JK136_citup_iter_clone_freqs <- melt(JK136_citup_iter_clone_freqs, id.var = "clone", variable.name = "sample", value.name = "clone_prev")

# Plot
ggplot(JK136_citup_iter_clone_freqs, aes(x = sample, y = clone_prev, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = unname(iwanthue(length(levels(JK136_citup_iter_clone_freqs$clone))))) +
  labs(x = NULL, y = "clonal prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

# Load CITUP results
JK142_citup_iter <- h5dump("./results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/citup_iter/iter_results.h5")

# Get optimal solution
JK142_citup_iter_sol <- as.character(JK142_citup_iter$results$optimal$index[1])

# Get clone info from optimal solution
JK142_citup_iter_clones <- setNames(reshape2::melt(JK142_citup_iter$trees[[JK142_citup_iter_sol]]$clone_freq$block0_values),
                                    c("clone_id", "sample_id", "clonal_prev"))
JK142_citup_iter_clones$sample_id <- c("JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_org", "JK196_tis", "JK196_org")[JK142_citup_iter_clones$sample_id]
JK142_citup_iter_clones$clone_id <- JK142_citup_iter_clones$clone_id - 1

# Load PyClone variants (files generated in prep_pyclone.R)
JK142_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202/pyclone/JK142_reg1_tissue_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202/pyclone/JK142_reg1_organoid_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202/pyclone/JK142_reg2_organoid_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK196_tissue/hg38_no_alt/EXC/B58186_B58202/pyclone/JK196_tissue_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK196_organoid/hg38_no_alt/EXC/B58187_B58202/pyclone/JK196_organoid_titan_subclonefil_adfil_noreg2tis.tsv"))
JK142_pyclone_vars <- JK142_pyclone_vars[!duplicated(JK142_pyclone_vars$mutation_id),]

# Load PyClone loci-level results
JK142_pyclone_loci <- read.delim("./results/exome/pyclone/JK142_subclonefil_adfil_noreg2org/tables/loci.tsv", stringsAsFactors = FALSE)
JK142_pyclone_loci$sample_id <- factor(JK142_pyclone_loci$sample_id,
                                       levels = c("JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_org", "JK196_tis", "JK196_org"))

# Merge PyClone results with variants to get copy number
JK142_pyclone_loci <- merge(JK142_pyclone_loci, JK142_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK142_pyclone_loci$cluster_id <- as.factor(JK142_pyclone_loci$cluster_id)

# Load variants
JK142_vars <- get.var.info.from.table("./results/exome/common_mutations/JK142_JK196_noreg2tis_strelka_mutect/tp-baseline.table")
JK142_vars$tbl_nosnp$mutation_id <- paste0(JK142_vars$tbl_nosnp$gene, "_", JK142_vars$tbl_nosnp$CHROM, ":", 
                                           as.character(JK142_vars$tbl_nosnp$POS))

JK142_pyclone_loci <- merge(JK142_pyclone_loci, JK142_vars$tbl_nosnp[,c("mutation_id", "TYPE", "gene", "effect", "impact", "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT")],
                            by = "mutation_id", all.x = TRUE)

# Save object
saveRDS(JK142_pyclone_loci, "./results/exome/common_mutations/JK142_pyclone_loci.rds")

#Add clone info
JK142_pyclone_loci$clone_id <- factor(setNames(JK142_citup_iter$trees[[JK142_citup_iter_sol]]$variant_assignment$values, 
                                               c(0:10,12,14)[JK142_citup_iter$trees[[JK142_citup_iter_sol]]$variant_assignment$index + 1])
                                      [as.character(JK142_pyclone_loci$cluster_id)])

# Get tree
JK142_tree <- graph_from_data_frame(data.frame(source = JK142_citup_iter$trees[[JK142_citup_iter_sol]]$adjacency_list$block0_values[1,],
                                               target = JK142_citup_iter$trees[[JK142_citup_iter_sol]]$adjacency_list$block0_values[2,]))

# Get number of variants per clone
JK142_clone_muts <- table(JK142_pyclone_loci$clone_id) / length(levels(JK142_pyclone_loci$sample_id))

# Build tree with branch lengths proportional to the number of variants
JK142_coords <- layout_(JK142_tree, as_tree())
JK142_coords[,2] <- c(-JK142_clone_muts["0"], -sum(JK142_clone_muts[c("0", "1")]), 
                      -sum(JK142_clone_muts[c("0", "1", "2")]), -sum(JK142_clone_muts[c("0", "1", "2", "3")]),
                      -sum(JK142_clone_muts[c("0", "1", "5")]), -sum(JK142_clone_muts[c("0", "1", "5", "6")]),
                      -sum(JK142_clone_muts[c("0", "9")]), -sum(JK142_clone_muts[c("0", "9", "10")]),
                      -sum(JK142_clone_muts[c("0", "1", "2", "3", "4")]), 0,
                      -sum(JK142_clone_muts[c("0", "1", "8")]),
                      -sum(JK142_clone_muts[c("0", "9", "10", "11")]),
                      -sum(JK142_clone_muts[c("0", "9", "12")]))
JK142_coords[,1] <- 0

# Plot tree
plot.igraph(JK142_tree, layout = JK142_coords)

# Get variant info for each clone
lapply(levels(JK142_pyclone_loci$clone_id), function(clone) {
  clone_info <- JK142_pyclone_loci[which(JK142_pyclone_loci$clone_id == clone & JK142_pyclone_loci$sample_id == "JK142_reg1_tis"), c(11:16)]

message("genes with high impact variants")
clone_info[clone_info$impact == "HIGH",]

message("GBM genes")
clone_info[clone_info$gene %in% brennan_genes,]

message("# of variants")
nrow(clone_info)
})

# Get proportions of malignant cells assigned to each clone
JK142_citup_iter_clone_freqs <- setNames(as.data.frame(JK142_citup_iter$trees[[JK142_citup_iter_sol]]$clone_freq$block0_values),
                                         levels(JK142_pyclone_loci$sample_id))
JK142_citup_iter_clone_freqs$clone <- as.factor(0:(nrow(JK142_citup_iter_clone_freqs) - 1))

JK142_citup_iter_clone_freqs <- melt(JK142_citup_iter_clone_freqs, id.var = "clone", variable.name = "sample", value.name = "clone_prev")

# Plot
ggplot(JK142_citup_iter_clone_freqs, aes(x = sample, y = clone_prev, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = unname(iwanthue(length(levels(JK142_citup_iter_clone_freqs$clone))))) +
  labs(x = NULL, y = "clonal prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

# Load CITUP results
JK153_citup_iter <- h5dump("./results/exome/pyclone/JK153_subclonefil_adfil/citup_iter/iter_results.h5")

# Get optimal solution
JK153_citup_iter_sol <- as.character(JK153_citup_iter$results$optimal$index[1])

# Get clone info from optimal solution
JK153_citup_iter_clones <- setNames(reshape2::melt(JK153_citup_iter$trees[[JK153_citup_iter_sol]]$clone_freq$block0_values),
                                    c("clone_id", "sample_id", "clonal_prev"))
JK153_citup_iter_clones$sample_id <- c("JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org")[JK153_citup_iter_clones$sample_id]
JK153_citup_iter_clones$clone_id <- JK153_citup_iter_clones$clone_id - 1

# Load PyClone variants (files generated in prep_pyclone.R)
JK153_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203/pyclone/JK153_reg1_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203/pyclone/JK153_reg1_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203/pyclone/JK153_reg2_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203/pyclone/JK153_reg2_organoid_titan_subclonefil_adfil.tsv"))
JK153_pyclone_vars <- JK153_pyclone_vars[!duplicated(JK153_pyclone_vars$mutation_id),]

# Load PyClone loci-level results
JK153_pyclone_loci <- read.delim("./results/exome/pyclone/JK153_subclonefil_adfil_noreg2org/tables/loci.tsv", stringsAsFactors = FALSE)
JK153_pyclone_loci$sample_id <- factor(JK153_pyclone_loci$sample_id,
                                       levels = c("JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org"))

# Merge PyClone results with variants to get copy number
JK153_pyclone_loci <- merge(JK153_pyclone_loci, JK153_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK153_pyclone_loci$cluster_id <- as.factor(JK153_pyclone_loci$cluster_id)

# Load variants
JK153_vars <- get.var.info.from.table("./results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.table")
JK153_vars$tbl_nosnp$mutation_id <- paste0(JK153_vars$tbl_nosnp$gene, "_", JK153_vars$tbl_nosnp$CHROM, ":", 
                                           as.character(JK153_vars$tbl_nosnp$POS))

JK153_pyclone_loci <- merge(JK153_pyclone_loci, JK153_vars$tbl_nosnp[,c("mutation_id", "TYPE", "gene", "effect", "impact", "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT")],
                            by = "mutation_id", all.x = TRUE)

# Save object
saveRDS(JK153_pyclone_loci, "./results/exome/common_mutations/JK153_pyclone_loci.rds")

#Add clone info
JK153_pyclone_loci$clone_id <- factor(setNames(JK153_citup_iter$trees[[JK153_citup_iter_sol]]$variant_assignment$values, 
                                               c(0:4,6,7,10:14)[JK153_citup_iter$trees[[JK153_citup_iter_sol]]$variant_assignment$index + 1])
                                      [as.character(JK153_pyclone_loci$cluster_id)])

# Get tree
JK153_tree <- graph_from_data_frame(data.frame(source = JK153_citup_iter$trees[[JK153_citup_iter_sol]]$adjacency_list$block0_values[1,],
                                               target = JK153_citup_iter$trees[[JK153_citup_iter_sol]]$adjacency_list$block0_values[2,]))

# Get number of variants per clone
JK153_clone_muts <- table(JK153_pyclone_loci$clone_id) / length(levels(JK153_pyclone_loci$sample_id))

# Build tree with branch lengths proportional to the number of variants
JK153_coords <- layout_(JK153_tree, as_tree())
JK153_coords[,2] <- c(0, -JK153_clone_muts[1], -sum(JK153_clone_muts[1:2]), 
                      -sum(JK153_clone_muts[1:3]), -sum(JK153_clone_muts[1:4]), 
                      -sum(JK153_clone_muts[c(1:5)]), 
                      -sum(JK153_clone_muts[c(1,2,9)]),
                      -sum(JK153_clone_muts[c(1:6)]), 
                      -sum(JK153_clone_muts[c(1:5,7)]), -sum(JK153_clone_muts[c(1:5,8)]),
                      -sum(JK153_clone_muts[c(1,2,9,10)]),
                      -sum(JK153_clone_muts[c(1,11)]))
JK153_coords[,1] <- 0

# Plot tree
plot.igraph(JK153_tree, layout = JK153_coords)

# Get variant info for each clone
lapply(levels(JK153_pyclone_loci$clone_id), function(clone) {
  clone_info <- JK153_pyclone_loci[which(JK153_pyclone_loci$clone_id == clone & JK153_pyclone_loci$sample_id == "JK153_reg1_tis"), c(11:16)]
  
  message("genes with high impact variants")
  clone_info[clone_info$impact == "HIGH",]
  
  message("GBM genes")
  clone_info[clone_info$gene %in% brennan_genes,]
  
  message("# of variants")
  nrow(clone_info)
})

# Get proportions of malignant cells assigned to each clone
JK153_citup_iter_clone_freqs <- setNames(as.data.frame(JK153_citup_iter$trees[[JK153_citup_iter_sol]]$clone_freq$block0_values),
                                         levels(JK153_pyclone_loci$sample_id))
JK153_citup_iter_clone_freqs$clone <- as.factor(0:(nrow(JK153_citup_iter_clone_freqs) - 1))

JK153_citup_iter_clone_freqs <- melt(JK153_citup_iter_clone_freqs, id.var = "clone", variable.name = "sample", value.name = "clone_prev")

# Plot
ggplot(JK153_citup_iter_clone_freqs, aes(x = sample, y = clone_prev, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = unname(iwanthue(length(levels(JK153_citup_iter_clone_freqs$clone))))) +
  labs(x = NULL, y = "clonal prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







# ---------------------------------------------------------------
## JK163
# ---------------------------------------------------------------

# Load CITUP results
JK163_citup_iter <- h5dump("./results/exome/pyclone/JK163_subclonefil_adfil50org2/citup_iter/iter_results.h5")

# Get optimal solution
JK163_citup_iter_sol <- as.character(JK163_citup_iter$results$optimal$index[1])

# Get clone info from optimal solution
JK163_citup_iter_clones <- setNames(reshape2::melt(JK163_citup_iter$trees[[JK163_citup_iter_sol]]$clone_freq$block0_values),
                                    c("clone_id", "sample_id", "clonal_prev"))
JK163_citup_iter_clones$sample_id <- c("JK163_reg1_tis", "JK163_reg1_org", "JK163_reg2_tis", "JK163_reg2_org")[JK163_citup_iter_clones$sample_id]
JK163_citup_iter_clones$clone_id <- JK163_citup_iter_clones$clone_id - 1

# Load PyClone variants (files generated in prep_pyclone.R)
JK163_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204/pyclone/JK163_reg1_tissue_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204/pyclone/JK163_reg1_organoid_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204/pyclone/JK163_reg2_tissue_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204/pyclone/JK163_reg2_organoid_titan_subclonefil_adfil50org2.tsv"))
JK163_pyclone_vars <- JK163_pyclone_vars[!duplicated(JK163_pyclone_vars$mutation_id),]

# Load PyClone loci-level results
JK163_pyclone_loci <- read.delim("./results/exome/pyclone/JK163_subclonefil_adfil50org2/tables/loci.tsv", stringsAsFactors = FALSE)
JK163_pyclone_loci$sample_id <- factor(JK163_pyclone_loci$sample_id,
                                       levels = c("JK163_reg1_tis", "JK163_reg1_org", "JK163_reg2_tis", "JK163_reg2_org"))

# Merge PyClone results with variants to get copy number
JK163_pyclone_loci <- merge(JK163_pyclone_loci, JK163_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK163_pyclone_loci$cluster_id <- as.factor(JK163_pyclone_loci$cluster_id)

# Load variants
JK163_vars <- get.var.info.from.table("./results/exome/common_mutations/JK163_strelka_mutect/tp-baseline.table")
JK163_vars$tbl_nosnp$mutation_id <- paste0(JK163_vars$tbl_nosnp$gene, "_", JK163_vars$tbl_nosnp$CHROM, ":", 
                                           as.character(JK163_vars$tbl_nosnp$POS))

JK163_pyclone_loci <- merge(JK163_pyclone_loci, JK163_vars$tbl_nosnp[,c("mutation_id", "TYPE", "gene", "effect", "impact", "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT")],
                            by = "mutation_id", all.x = TRUE)

# Save object
saveRDS(JK163_pyclone_loci, "./results/exome/common_mutations/JK163_pyclone_loci.rds")

#Add clone info
JK163_pyclone_loci$clone_id <- factor(setNames(JK163_citup_iter$trees[[JK163_citup_iter_sol]]$variant_assignment$values, 
                                               c(0,1,4,7,15,34,55)[JK163_citup_iter$trees[[JK163_citup_iter_sol]]$variant_assignment$index + 1])
                                      [as.character(JK163_pyclone_loci$cluster_id)])

# Get tree
JK163_tree <- graph_from_data_frame(data.frame(source = JK163_citup_iter$trees[[JK163_citup_iter_sol]]$adjacency_list$block0_values[1,],
                                               target = JK163_citup_iter$trees[[JK163_citup_iter_sol]]$adjacency_list$block0_values[2,]))

# Get number of variants per clone
JK163_clone_muts <- table(JK163_pyclone_loci$clone_id) / length(levels(JK163_pyclone_loci$sample_id))

# Build tree with branch lengths proportional to the number of variants
JK163_coords <- layout_(JK163_tree, as_tree())
JK163_coords[,2] <- c(0, -JK163_clone_muts[2], -sum(JK163_clone_muts[2:3]), -sum(JK163_clone_muts[2:4]),
                      -sum(JK163_clone_muts[2:5]), -sum(JK163_clone_muts[2:6]), -sum(JK163_clone_muts[c(2:4,7)]))
JK163_coords[,1] <- 0

# Plot tree
plot.igraph(JK163_tree, layout = JK163_coords)

# Get variant info for each clone
lapply(levels(JK163_pyclone_loci$clone_id), function(clone) {
  clone_info <- JK163_pyclone_loci[which(JK163_pyclone_loci$clone_id == clone & JK163_pyclone_loci$sample_id == "JK163_reg1_tis"), c(11:16)]
  
  message("genes with high impact variants")
  clone_info[clone_info$impact == "HIGH",]
  
  message("GBM genes")
  clone_info[clone_info$gene %in% brennan_genes,]
  
  message("# of variants")
  nrow(clone_info)
})

# Get proportions of malignant cells assigned to each clone
JK163_citup_iter_clone_freqs <- setNames(as.data.frame(JK163_citup_iter$trees[[JK163_citup_iter_sol]]$clone_freq$block0_values),
                                         levels(JK163_pyclone_loci$sample_id))
JK163_citup_iter_clone_freqs$clone <- as.factor(0:(nrow(JK163_citup_iter_clone_freqs) - 1))

JK163_citup_iter_clone_freqs <- melt(JK163_citup_iter_clone_freqs, id.var = "clone", variable.name = "sample", value.name = "clone_prev")

# Plot
ggplot(JK163_citup_iter_clone_freqs, aes(x = sample, y = clone_prev, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = unname(iwanthue(length(levels(JK163_citup_iter_clone_freqs$clone))))) +
  labs(x = NULL, y = "clonal prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
