library(data.table)


#JK124

#Get pyclone cluster values
JK124_pyclone_clust <- read.delim("./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/cluster.tsv", stringsAsFactors = FALSE)

#Fix sample names
JK124_pyclone_clust$sample_id <- factor(JK124_pyclone_clust$sample_id,
                                        levels = c("JK124_reg1_tis", "JK124_reg1_org", "JK124_reg2_tis"))

#Get variant info
JK124_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK124_reg1_tissue/hg38_no_alt/EXC/B58192_B58200/pyclone/JK124_reg1_tissue_titan_subclonefil_adfil_noreg2org.tsv"),
                            read.delim("./data/exome_seq/JK124_reg1_organoid/hg38_no_alt/EXC/B58193_B58200/pyclone/JK124_reg1_organoid_titan_subclonefil_adfil_noreg2org.tsv"),
                            read.delim("./data/exome_seq/JK124_reg2_tissue/hg38_no_alt/EXC/B58194_B58200/pyclone/JK124_reg2_tissue_titan_subclonefil_adfil_noreg2org.tsv"))
JK124_pyclone_vars <- JK124_pyclone_vars[!duplicated(JK124_pyclone_vars$mutation_id),]

JK124_pyclone_loci <- read.delim("./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/tables/loci.tsv", stringsAsFactors = FALSE)
JK124_pyclone_loci$sample_id <- factor(JK124_pyclone_loci$sample_id,
                                       levels = c("JK124_reg1_tis", "JK124_reg1_org", "JK124_reg2_tis"))
JK124_pyclone_loci <- merge(JK124_pyclone_loci, JK124_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK124_pyclone_loci$cluster_id <- as.factor(JK124_pyclone_loci$cluster_id)

#Identify clusters to filter (singletons) and check that SNPs and indels don't cluster separately
table(JK124_pyclone_loci$cluster_id) / length(levels(JK124_pyclone_clust$sample_id))
sapply(levels(JK124_pyclone_loci$cluster_id), function(clust) table(JK124_pyclone_loci$type[JK124_pyclone_loci$cluster_id == clust]) / length(levels(JK124_pyclone_clust$sample_id)))

#Get cluster standard deviations across samples
dcast(JK124_pyclone_clust, cluster_id ~ sample_id, value.var = "std")


#CITUP_iter
#Make input table
JK124_citup_iter_input <- dcast(JK124_pyclone_clust, cluster_id ~ sample_id, value.var = "mean")

#Write table
write.table(JK124_citup_iter_input[,2:ncol(JK124_citup_iter_input)], 
            file = "./results/exome/pyclone/JK124_subclonefil_adfil_noreg2org/citup_iter/JK124_subclonefil_adfil_cluster_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")







#JK136

#Get pyclone cluster values
JK136_pyclone_clust <- read.delim("./results/exome/pyclone/JK136_subclonefil_adfil/tables/cluster.tsv", stringsAsFactors = FALSE)

#Fix sample names
JK136_pyclone_clust$sample_id <- factor(JK136_pyclone_clust$sample_id,
                                        levels = c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org",
                                                   "JK202_tis", "JK202_org"))

#Get variant info
JK136_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK136_reg1_tissue/hg38_no_alt/EXC/B58176_B58201/pyclone/JK136_reg1_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg1_organoid/hg38_no_alt/EXC/B58177_B58201/pyclone/JK136_reg1_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg2_tissue/hg38_no_alt/EXC/B58178_B58201/pyclone/JK136_reg2_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK136_reg2_organoid/hg38_no_alt/EXC/B58179_B58201/pyclone/JK136_reg2_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK202_tissue/hg38_no_alt/EXC/B58180_B58201/pyclone/JK202_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK202_organoid/hg38_no_alt/EXC/B58181_B58201/pyclone/JK202_organoid_titan_subclonefil_adfil.tsv"))
JK136_pyclone_vars <- JK136_pyclone_vars[!duplicated(JK136_pyclone_vars$mutation_id),]

JK136_pyclone_loci <- read.delim("./results/exome/pyclone/JK136_subclonefil_adfil/tables/loci.tsv", stringsAsFactors = FALSE)
JK136_pyclone_loci$sample_id <- factor(JK136_pyclone_loci$sample_id,
                                       levels = c("JK136_reg1_tis", "JK136_reg1_org", "JK136_reg2_tis", "JK136_reg2_org",
                                                  "JK202_tis", "JK202_org"))
JK136_pyclone_loci <- merge(JK136_pyclone_loci, JK136_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK136_pyclone_loci$cluster_id <- as.factor(JK136_pyclone_loci$cluster_id)

#Identify clusters to filter (singletons) and check that SNPs and indels don't cluster separately
table(JK136_pyclone_loci$cluster_id) / length(levels(JK136_pyclone_clust$sample_id))
sapply(levels(JK136_pyclone_loci$cluster_id), function(clust) table(JK136_pyclone_loci$type[JK136_pyclone_loci$cluster_id == clust]) / length(levels(JK136_pyclone_clust$sample_id)))

#Get cluster standard deviations across samples
dcast(JK136_pyclone_clust, cluster_id ~ sample_id, value.var = "std")


#CITUP_iter
#Make input table
JK136_citup_iter_input <- dcast(JK136_pyclone_clust[JK136_pyclone_clust$cluster_id %in% c(0:4),], cluster_id ~ sample_id, value.var = "mean")

#Write table
write.table(JK136_citup_iter_input[,2:ncol(JK136_citup_iter_input)], 
            file = "./results/exome/pyclone/JK136_subclonefil_adfil/citup_iter/JK136_subclonefil_adfil_cluster_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")




#Make input table for CITUP_qip
JK136_citup_qip_input <- dcast(JK136_pyclone_loci[JK136_pyclone_loci$cluster_id %in% c(0:3),], 
                               mutation_id + cluster_id ~ sample_id, value.var = "cellular_prevalence")

#Write tables
write.table(JK136_citup_qip_input[,3:ncol(JK136_citup_qip_input)], 
            file = "./results/exome/pyclone/JK136_subclonefil/citup_qip/JK136_subclonefil_var_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(JK136_citup_qip_input[, 2, drop = F], 
            file = "./results/exome/pyclone/JK136_subclonefil/citup_qip/JK136_subclonefil_var_clusters.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")







#JK142

#Get pyclone cluster values
JK142_pyclone_clust <- read.delim("./results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables/cluster.tsv", stringsAsFactors = FALSE)

#Fix sample names
JK142_pyclone_clust$sample_id <- factor(JK142_pyclone_clust$sample_id,
                                        levels = c("JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_org",
                                                   "JK196_tis", "JK196_org"))

#Get variant info
JK142_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK142_reg1_tissue/hg38_no_alt/EXC/B58182_B58202/pyclone/JK142_reg1_tissue_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK142_reg1_organoid/hg38_no_alt/EXC/B58183_B58202/pyclone/JK142_reg1_organoid_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK142_reg2_organoid/hg38_no_alt/EXC/B58185_B58202/pyclone/JK142_reg2_organoid_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK196_tissue/hg38_no_alt/EXC/B58186_B58202/pyclone/JK196_tissue_titan_subclonefil_adfil_noreg2tis.tsv"),
                            read.delim("./data/exome_seq/JK196_organoid/hg38_no_alt/EXC/B58187_B58202/pyclone/JK196_organoid_titan_subclonefil_adfil_noreg2tis.tsv"))
JK142_pyclone_vars <- JK142_pyclone_vars[!duplicated(JK142_pyclone_vars$mutation_id),]

JK142_pyclone_loci <- read.delim("./results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/tables/loci.tsv", stringsAsFactors = FALSE)
JK142_pyclone_loci$sample_id <- factor(JK142_pyclone_loci$sample_id,
                                       levels = c("JK142_reg1_tis", "JK142_reg1_org", "JK142_reg2_org",
                                                  "JK196_tis", "JK196_org"))
JK142_pyclone_loci <- merge(JK142_pyclone_loci, JK142_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK142_pyclone_loci$cluster_id <- as.factor(JK142_pyclone_loci$cluster_id)

#Identify clusters to filter (singletons) and check that SNPs and indels don't cluster separately
table(JK142_pyclone_loci$cluster_id) / length(levels(JK142_pyclone_clust$sample_id))
sapply(levels(JK142_pyclone_loci$cluster_id), function(clust) table(JK142_pyclone_loci$type[JK142_pyclone_loci$cluster_id == clust]) / length(levels(JK142_pyclone_clust$sample_id)))

#Get cluster standard deviations across samples
dcast(JK142_pyclone_clust, cluster_id ~ sample_id, value.var = "std")


#CITUP_iter
#Make input table
JK142_citup_iter_input <- dcast(JK142_pyclone_clust[JK142_pyclone_clust$cluster_id %in% c(0:10,12,14),], cluster_id ~ sample_id, value.var = "mean")

#Write table
write.table(JK142_citup_iter_input[,2:ncol(JK142_citup_iter_input)], 
            file = "./results/exome/pyclone/JK142_subclonefil_adfil_noreg2tis/citup_iter/JK142_subclonefil_adfil_noreg2tis_cluster_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")









#JK153

#Get pyclone cluster values
JK153_pyclone_clust <- read.delim("./results/exome/pyclone/JK153_subclonefil_adfil/tables/cluster.tsv", stringsAsFactors = FALSE)

#Fix sample names
JK153_pyclone_clust$sample_id <- factor(JK153_pyclone_clust$sample_id,
                                        levels = c("JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org"))

#Get variant info
JK153_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK153_reg1_tissue/hg38_no_alt/EXC/B58188_B58203/pyclone/JK153_reg1_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg1_organoid/hg38_no_alt/EXC/B58189_B58203/pyclone/JK153_reg1_organoid_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg2_tissue/hg38_no_alt/EXC/B58190_B58203/pyclone/JK153_reg2_tissue_titan_subclonefil_adfil.tsv"),
                            read.delim("./data/exome_seq/JK153_reg2_organoid/hg38_no_alt/EXC/B58191_B58203//pyclone/JK153_reg2_organoid_titan_subclonefil_adfil.tsv"))
JK153_pyclone_vars <- JK153_pyclone_vars[!duplicated(JK153_pyclone_vars$mutation_id),]

JK153_pyclone_loci <- read.delim("./results/exome/pyclone/JK153_subclonefil_adfil/tables/loci.tsv", stringsAsFactors = FALSE)
JK153_pyclone_loci$sample_id <- factor(JK153_pyclone_loci$sample_id,
                                       levels = c("JK153_reg1_tis", "JK153_reg1_org", "JK153_reg2_tis", "JK153_reg2_org"))
JK153_pyclone_loci <- merge(JK153_pyclone_loci, JK153_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK153_pyclone_loci$cluster_id <- as.factor(JK153_pyclone_loci$cluster_id)

#Identify clusters to filter (singletons) and check that SNPs and indels don't cluster separately
table(JK153_pyclone_loci$cluster_id) / length(levels(JK153_pyclone_clust$sample_id))
sapply(levels(JK153_pyclone_loci$cluster_id), function(clust) table(JK153_pyclone_loci$type[JK153_pyclone_loci$cluster_id == clust]) / length(levels(JK153_pyclone_clust$sample_id)))

#Get cluster standard deviations across samples
dcast(JK153_pyclone_clust, cluster_id ~ sample_id, value.var = "std")


#CITUP_iter
#Make input table
JK153_citup_iter_input <- dcast(JK153_pyclone_clust[JK153_pyclone_clust$cluster_id %in% c(0:4,6,7,10:14),], cluster_id ~ sample_id, value.var = "mean")

#Write table
write.table(JK153_citup_iter_input[,2:ncol(JK153_citup_iter_input)], 
            file = "./results/exome/pyclone/JK153_subclonefil_adfil/citup_iter/JK153_subclonefil_adfil_cluster_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")







#JK163

#Get pyclone cluster values
JK163_pyclone_clust <- read.delim("./results/exome/pyclone/JK163_subclonefil_adfil50org2/tables/cluster.tsv", stringsAsFactors = FALSE)

#Fix sample names
JK163_pyclone_clust$sample_id <- factor(JK163_pyclone_clust$sample_id,
                                        levels = c("JK163_reg1_tis", "JK163_reg1_org", "JK163_reg2_tis", "JK163_reg2_org"))

#Get variant info
JK163_pyclone_vars <- rbind(read.delim("./data/exome_seq/JK163_reg1_tissue/hg38_no_alt/EXC/B58196_B58204/pyclone/JK163_reg1_tissue_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg1_organoid/hg38_no_alt/EXC/B58197_B58204/pyclone/JK163_reg1_organoid_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg2_tissue/hg38_no_alt/EXC/B58198_B58204/pyclone/JK163_reg2_tissue_titan_subclonefil_adfil50org2.tsv"),
                            read.delim("./data/exome_seq/JK163_reg2_organoid/hg38_no_alt/EXC/B58199_B58204//pyclone/JK163_reg2_organoid_titan_subclonefil_adfil50org2.tsv"))
JK163_pyclone_vars <- JK163_pyclone_vars[!duplicated(JK163_pyclone_vars$mutation_id),]

JK163_pyclone_loci <- read.delim("./results/exome/pyclone/JK163_subclonefil_adfil50org2/tables/loci.tsv", stringsAsFactors = FALSE)
JK163_pyclone_loci$sample_id <- factor(JK163_pyclone_loci$sample_id,
                                       levels = c("JK163_reg1_tis", "JK163_reg1_org", "JK163_reg2_tis", "JK163_reg2_org"))
JK163_pyclone_loci <- merge(JK163_pyclone_loci, JK163_pyclone_vars[,c("mutation_id", "type", "clonal_cluster", "cell_prev")],
                            by = "mutation_id")
JK163_pyclone_loci$cluster_id <- as.factor(JK163_pyclone_loci$cluster_id)

#Identify clusters to filter (singletons) and check that SNPs and indels don't cluster separately
table(JK163_pyclone_loci$cluster_id) / length(levels(JK163_pyclone_clust$sample_id))
sapply(levels(JK163_pyclone_loci$cluster_id), function(clust) table(JK163_pyclone_loci$type[JK163_pyclone_loci$cluster_id == clust]) / length(levels(JK163_pyclone_clust$sample_id)))

#Get cluster standard deviations across samples
dcast(JK163_pyclone_clust, cluster_id ~ sample_id, value.var = "std")


#CITUP_iter
#Make input table
JK163_citup_iter_input <- dcast(JK163_pyclone_clust[JK163_pyclone_clust$cluster_id %in% c(0,1,4,7,15,34,55),], cluster_id ~ sample_id, value.var = "mean")

#Write table
write.table(JK163_citup_iter_input[,2:ncol(JK163_citup_iter_input)], 
            file = "./results/exome/pyclone/JK163_subclonefil_adfil50org2/citup_iter/JK163_subclonefil_adfil50org2_cluster_freqs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

