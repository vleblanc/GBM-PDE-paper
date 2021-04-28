# ---------------------------------------------------------------
## Scripts to look at coverage of variants identified in WES data in merged cells from clones identified in scWGS data
## Author: Veronique LeBlanc
# ---------------------------------------------------------------


# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

# Load genotyping results for each clone (generated in mergebams_genotype.sh) and make informative data frame
JK136_cnv_a_muts <- read.vcfR("./data/cellranger_dna/JK136_JK202_mergebams/clone_a.vcf.gz", verbose = FALSE)
JK136_cnv_a_muts <- data.frame(CHROM = JK136_cnv_a_muts@fix[,"CHROM"],
                               POS = JK136_cnv_a_muts@fix[,"POS"],
                               a_ref = JK136_cnv_a_muts@fix[,"REF"],
                               a_alt = JK136_cnv_a_muts@fix[,"ALT"],
                               a_qual = as.numeric(JK136_cnv_a_muts@fix[,"QUAL"]),
                               a_ref_dp = sapply(strsplit(JK136_cnv_a_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               a_alt_dp = sapply(strsplit(JK136_cnv_a_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK136_cnv_a_muts$a_dp <- JK136_cnv_a_muts$a_ref_dp + JK136_cnv_a_muts$a_alt_dp


JK136_cnv_c_muts <- read.vcfR("./data/cellranger_dna/JK136_JK202_mergebams/clone_c.vcf.gz", verbose = FALSE)
JK136_cnv_c_muts <- data.frame(CHROM = JK136_cnv_c_muts@fix[,"CHROM"],
                               POS = JK136_cnv_c_muts@fix[,"POS"],
                               c_ref = JK136_cnv_c_muts@fix[,"REF"],
                               c_alt = JK136_cnv_c_muts@fix[,"ALT"],
                               c_qual = as.numeric(JK136_cnv_c_muts@fix[,"QUAL"]),
                               c_ref_dp = sapply(strsplit(JK136_cnv_c_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               c_alt_dp = sapply(strsplit(JK136_cnv_c_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK136_cnv_c_muts$c_dp <- JK136_cnv_c_muts$c_ref_dp + JK136_cnv_c_muts$c_alt_dp

JK136_cnv_nonmal_muts <- read.vcfR("./data/cellranger_dna/JK136_JK202_mergebams/clone_non-malignant.vcf.gz", verbose = FALSE)
JK136_cnv_nonmal_muts <- data.frame(CHROM = JK136_cnv_nonmal_muts@fix[,"CHROM"],
                                    POS = JK136_cnv_nonmal_muts@fix[,"POS"],
                                    nonmal_ref = JK136_cnv_nonmal_muts@fix[,"REF"],
                                    nonmal_alt = JK136_cnv_nonmal_muts@fix[,"ALT"],
                                    nonmal_qual = as.numeric(JK136_cnv_nonmal_muts@fix[,"QUAL"]),
                                    nonmal_ref_dp = sapply(strsplit(JK136_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                           }),
                                    nonmal_alt_dp = sapply(strsplit(JK136_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                           }),
                                    stringsAsFactors = FALSE)
JK136_cnv_nonmal_muts$nonmal_dp <- JK136_cnv_nonmal_muts$nonmal_ref_dp + JK136_cnv_nonmal_muts$nonmal_alt_dp

# Merge with loci information
JK136_pyclone_loci <- readRDS("./results/exome/common_mutations/JK136_pyclone_loci.rds") #generated in fig2_trees.R

JK136_cnv_muts_ol <- merge(JK136_pyclone_loci[seq(1, nrow(JK136_pyclone_loci), 6), c("mutation_id", "type", "clonal_cluster", "cell_prev", "gene", "effect", "impact", 
                                                                                     "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT", "clone_id")], 
                           JK136_cnv_a_muts, by = c("CHROM", "POS"), all.x = TRUE)

# Qualify overlaps (no coverage, low quality [QUAL < 20], flag if indel for manual review, low depth if <5 reads)
JK136_cnv_muts_ol$a_mut <- ifelse(is.na(JK136_cnv_muts_ol$a_dp), "no_cov", 
                                  ifelse(as.numeric(JK136_cnv_muts_ol$a_qual) < 20, "low_qual",
                                         ifelse(JK136_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK136_cnv_muts_ol$a_alt) & JK136_cnv_muts_ol$a_alt == JK136_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK136_cnv_muts_ol$a_dp < 5, "low_dp", "no")))))

# Repeat for other clones and non-malignant cells
JK136_cnv_muts_ol <- merge(JK136_cnv_muts_ol, JK136_cnv_c_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK136_cnv_muts_ol$c_mut <- ifelse(is.na(JK136_cnv_muts_ol$c_dp), "no_cov", 
                                  ifelse(as.numeric(JK136_cnv_muts_ol$c_qual) < 20, "low_qual",
                                         ifelse(JK136_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK136_cnv_muts_ol$c_alt) & JK136_cnv_muts_ol$c_alt == JK136_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK136_cnv_muts_ol$c_dp < 5, "low_dp", "no")))))

JK136_cnv_muts_ol <- merge(JK136_cnv_muts_ol, JK136_cnv_nonmal_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK136_cnv_muts_ol$nonmal_mut <- ifelse(is.na(JK136_cnv_muts_ol$nonmal_dp), "no_cov", 
                                       ifelse(as.numeric(JK136_cnv_muts_ol$nonmal_qual) < 20, "low_qual",
                                              ifelse(JK136_cnv_muts_ol$type == "INDEL", "flag",
                                                     ifelse(!is.na(JK136_cnv_muts_ol$nonmal_alt) & 
                                                              JK136_cnv_muts_ol$nonmal_alt == JK136_cnv_muts_ol$ALT, "yes",
                                                            ifelse(JK136_cnv_muts_ol$nonmal_dp < 5, "low_dp", "no")))))


#Add mean cellular prevalence for each variant
JK136_cnv_muts_ol <- merge(JK136_cnv_muts_ol, aggregate(JK136_pyclone_loci$cellular_prevalence, by = list(JK136_pyclone_loci$mutation_id), mean),
                           by.x = "mutation_id", by.y = "Group.1")
colnames(JK136_cnv_muts_ol)[colnames(JK136_cnv_muts_ol) == "x"] <- "mean_cell_prev"

write.csv(JK136_cnv_muts_ol, file = "./results/exome/cnv_integration/JK136_cnv_clone_mutations.csv", row.names = FALSE)
#fix indels and duplicates (usually near reads with indels) manually

# Re-load cleaned up data
JK136_cnv_muts_ol <- read.csv("./results/exome/cnv_integration/JK136_cnv_clone_mutations.csv", stringsAsFactors = FALSE)
JK136_cnv_muts_ol$clone_id <- factor(JK136_cnv_muts_ol$clone_id)
JK136_cnv_muts_ol[,grepl("_mut", colnames(JK136_cnv_muts_ol))] <- lapply(JK136_cnv_muts_ol[,grepl("_mut", colnames(JK136_cnv_muts_ol))],
                                                                         function(x) factor(x, levels = c("yes", "no", "low_dp", "low_qual", "no_cov")))

# View overlaps by clone (coloured in Figure 2 if >10% of mutation clone variants are found in copy number clone)
setNames(lapply(colnames(JK136_cnv_muts_ol)[grepl("_mut", colnames(JK136_cnv_muts_ol))], function(cnv_clone) {
  sapply(levels(JK136_cnv_muts_ol$clone_id), function(clone) table(JK136_cnv_muts_ol[which(JK136_cnv_muts_ol$clone_id == clone), cnv_clone]) / 
           length(which(JK136_cnv_muts_ol$clone_id == clone)))
}), colnames(JK136_cnv_muts_ol)[grepl("_mut", colnames(JK136_cnv_muts_ol))])







# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

# Load genotyping results for each clone (generated in mergebams_genotype.sh) and make informative data frame
JK142_cnv_a_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_a.vcf.gz", verbose = FALSE)
JK142_cnv_a_muts <- data.frame(CHROM = JK142_cnv_a_muts@fix[,"CHROM"],
                               POS = JK142_cnv_a_muts@fix[,"POS"],
                               a_ref = JK142_cnv_a_muts@fix[,"REF"],
                               a_alt = JK142_cnv_a_muts@fix[,"ALT"],
                               a_qual = as.numeric(JK142_cnv_a_muts@fix[,"QUAL"]),
                               a_ref_dp = sapply(strsplit(JK142_cnv_a_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               a_alt_dp = sapply(strsplit(JK142_cnv_a_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_a_muts$a_dp <- JK142_cnv_a_muts$a_ref_dp + JK142_cnv_a_muts$a_alt_dp

JK142_cnv_b_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_b.vcf.gz", verbose = FALSE)
JK142_cnv_b_muts <- data.frame(CHROM = JK142_cnv_b_muts@fix[,"CHROM"],
                               POS = JK142_cnv_b_muts@fix[,"POS"],
                               b_ref = JK142_cnv_b_muts@fix[,"REF"],
                               b_alt = JK142_cnv_b_muts@fix[,"ALT"],
                               b_qual = as.numeric(JK142_cnv_b_muts@fix[,"QUAL"]),
                               b_ref_dp = sapply(strsplit(JK142_cnv_b_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               b_alt_dp = sapply(strsplit(JK142_cnv_b_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_b_muts$b_dp <- JK142_cnv_b_muts$b_ref_dp + JK142_cnv_b_muts$b_alt_dp

JK142_cnv_d_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_d.vcf.gz", verbose = FALSE)
JK142_cnv_d_muts <- data.frame(CHROM = JK142_cnv_d_muts@fix[,"CHROM"],
                               POS = JK142_cnv_d_muts@fix[,"POS"],
                               d_ref = JK142_cnv_d_muts@fix[,"REF"],
                               d_alt = JK142_cnv_d_muts@fix[,"ALT"],
                               d_qual = as.numeric(JK142_cnv_d_muts@fix[,"QUAL"]),
                               d_ref_dp = sapply(strsplit(JK142_cnv_d_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               d_alt_dp = sapply(strsplit(JK142_cnv_d_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_d_muts$d_dp <- JK142_cnv_d_muts$d_ref_dp + JK142_cnv_d_muts$d_alt_dp

JK142_cnv_e_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_e.vcf.gz", verbose = FALSE)
JK142_cnv_e_muts <- data.frame(CHROM = JK142_cnv_e_muts@fix[,"CHROM"],
                               POS = JK142_cnv_e_muts@fix[,"POS"],
                               e_ref = JK142_cnv_e_muts@fix[,"REF"],
                               e_alt = JK142_cnv_e_muts@fix[,"ALT"],
                               e_qual = as.numeric(JK142_cnv_e_muts@fix[,"QUAL"]),
                               e_ref_dp = sapply(strsplit(JK142_cnv_e_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               e_alt_dp = sapply(strsplit(JK142_cnv_e_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_e_muts$e_dp <- JK142_cnv_e_muts$e_ref_dp + JK142_cnv_e_muts$e_alt_dp

JK142_cnv_j_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_j.vcf.gz", verbose = FALSE)
JK142_cnv_j_muts <- data.frame(CHROM = JK142_cnv_j_muts@fix[,"CHROM"],
                               POS = JK142_cnv_j_muts@fix[,"POS"],
                               j_ref = JK142_cnv_j_muts@fix[,"REF"],
                               j_alt = JK142_cnv_j_muts@fix[,"ALT"],
                               j_qual = as.numeric(JK142_cnv_j_muts@fix[,"QUAL"]),
                               j_ref_dp = sapply(strsplit(JK142_cnv_j_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               j_alt_dp = sapply(strsplit(JK142_cnv_j_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_j_muts$j_dp <- JK142_cnv_j_muts$j_ref_dp + JK142_cnv_j_muts$j_alt_dp

JK142_cnv_m_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_m.vcf.gz", verbose = FALSE)
JK142_cnv_m_muts <- data.frame(CHROM = JK142_cnv_m_muts@fix[,"CHROM"],
                               POS = JK142_cnv_m_muts@fix[,"POS"],
                               m_ref = JK142_cnv_m_muts@fix[,"REF"],
                               m_alt = JK142_cnv_m_muts@fix[,"ALT"],
                               m_qual = as.numeric(JK142_cnv_m_muts@fix[,"QUAL"]),
                               m_ref_dp = sapply(strsplit(JK142_cnv_m_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               m_alt_dp = sapply(strsplit(JK142_cnv_m_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK142_cnv_m_muts$m_dp <- JK142_cnv_m_muts$m_ref_dp + JK142_cnv_m_muts$m_alt_dp

JK142_cnv_nonmal_muts <- read.vcfR("./data/cellranger_dna/JK142_JK196_mergebams/clone_non-malignant.vcf.gz", verbose = FALSE)
JK142_cnv_nonmal_muts <- data.frame(CHROM = JK142_cnv_nonmal_muts@fix[,"CHROM"],
                                    POS = JK142_cnv_nonmal_muts@fix[,"POS"],
                                    nonmal_ref = JK142_cnv_nonmal_muts@fix[,"REF"],
                                    nonmal_alt = JK142_cnv_nonmal_muts@fix[,"ALT"],
                                    nonmal_qual = as.numeric(JK142_cnv_nonmal_muts@fix[,"QUAL"]),
                                    nonmal_ref_dp = sapply(strsplit(JK142_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                           }),
                                    nonmal_alt_dp = sapply(strsplit(JK142_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                           }),
                                    stringsAsFactors = FALSE)
JK142_cnv_nonmal_muts$nonmal_dp <- JK142_cnv_nonmal_muts$nonmal_ref_dp + JK142_cnv_nonmal_muts$nonmal_alt_dp

# Merge with loci information
JK142_pyclone_loci <- readRDS("./results/exome/common_mutations/JK142_pyclone_loci.rds") #generated in fig2_trees.R

JK142_cnv_muts_ol <- merge(JK142_pyclone_loci[seq(1, nrow(JK142_pyclone_loci), 6), c("mutation_id", "type", "clonal_cluster", "cell_prev", "gene", "effect", "impact", 
                                                                                     "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT", "clone_id")], 
                           JK142_cnv_a_muts, by = c("CHROM", "POS"), all.x = TRUE)

# Qualify overlaps (no coverage, low quality [QUAL < 20], flag if indel for manual review, low depth if <5 reads)
JK142_cnv_muts_ol$a_mut <- ifelse(is.na(JK142_cnv_muts_ol$a_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$a_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$a_alt) & JK142_cnv_muts_ol$a_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$a_dp < 5, "low_dp", "no")))))

# Repeat for other clones and non-malignant cells
JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_b_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$b_mut <- ifelse(is.na(JK142_cnv_muts_ol$b_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$b_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$b_alt) & JK142_cnv_muts_ol$b_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$b_dp < 5, "low_dp", "no")))))

JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_d_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$d_mut <- ifelse(is.na(JK142_cnv_muts_ol$d_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$d_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$d_alt) & JK142_cnv_muts_ol$d_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$d_dp < 5, "low_dp", "no")))))

JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_e_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$e_mut <- ifelse(is.na(JK142_cnv_muts_ol$e_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$e_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$e_alt) & JK142_cnv_muts_ol$e_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$e_dp < 5, "low_dp", "no")))))

JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_j_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$j_mut <- ifelse(is.na(JK142_cnv_muts_ol$j_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$j_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$j_alt) & JK142_cnv_muts_ol$j_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$j_dp < 5, "low_dp", "no")))))

JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_m_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$m_mut <- ifelse(is.na(JK142_cnv_muts_ol$m_dp), "no_cov", 
                                  ifelse(as.numeric(JK142_cnv_muts_ol$m_qual) < 20, "low_qual",
                                         ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK142_cnv_muts_ol$m_alt) & JK142_cnv_muts_ol$m_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK142_cnv_muts_ol$m_dp < 5, "low_dp", "no")))))

JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, JK142_cnv_nonmal_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK142_cnv_muts_ol$nonmal_mut <- ifelse(is.na(JK142_cnv_muts_ol$nonmal_dp), "no_cov", 
                                       ifelse(as.numeric(JK142_cnv_muts_ol$nonmal_qual) < 20, "low_qual",
                                              ifelse(JK142_cnv_muts_ol$type == "INDEL", "flag",
                                                     ifelse(!is.na(JK142_cnv_muts_ol$nonmal_alt) & 
                                                              JK142_cnv_muts_ol$nonmal_alt == JK142_cnv_muts_ol$ALT, "yes",
                                                            ifelse(JK142_cnv_muts_ol$nonmal_dp < 5, "low_dp", "no")))))


#Add mean cellular prevalence for each variant
JK142_cnv_muts_ol <- merge(JK142_cnv_muts_ol, aggregate(JK142_pyclone_loci$cellular_prevalence, by = list(JK142_pyclone_loci$mutation_id), mean),
                           by.x = "mutation_id", by.y = "Group.1")
colnames(JK142_cnv_muts_ol)[colnames(JK142_cnv_muts_ol) == "x"] <- "mean_cell_prev"

write.csv(JK142_cnv_muts_ol, file = "./results/exome/cnv_integration/JK142_cnv_clone_mutations.csv", row.names = FALSE)
#fix indels and duplicates (usually near reads with indels) manually

# Re-load cleaned up data
JK142_cnv_muts_ol <- read.csv("./results/exome/cnv_integration/JK142_cnv_clone_mutations.csv", stringsAsFactors = FALSE)
JK142_cnv_muts_ol$clone_id <- factor(JK142_cnv_muts_ol$clone_id)
JK142_cnv_muts_ol[,grepl("_mut", colnames(JK142_cnv_muts_ol))] <- lapply(JK142_cnv_muts_ol[,grepl("_mut", colnames(JK142_cnv_muts_ol))],
                                                                         function(x) factor(x, levels = c("yes", "no", "low_dp", "low_qual", "no_cov")))

# View overlaps by clone (coloured in Figure 2 if >10% of mutation clone variants are found in copy number clone)
setNames(lapply(colnames(JK142_cnv_muts_ol)[grepl("_mut", colnames(JK142_cnv_muts_ol))], function(cnv_clone) {
  sapply(levels(JK142_cnv_muts_ol$clone_id), function(clone) table(JK142_cnv_muts_ol[which(JK142_cnv_muts_ol$clone_id == clone), cnv_clone]) / 
           length(which(JK142_cnv_muts_ol$clone_id == clone)))
}), colnames(JK142_cnv_muts_ol)[grepl("_mut", colnames(JK142_cnv_muts_ol))])







# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

# Load genotyping results for each clone (generated in mergebams_genotype.sh) and make informative data frame
JK153_cnv_b_muts <- read.vcfR("./data/cellranger_dna/JK153_mergebams/clone_b.vcf.gz", verbose = FALSE)
JK153_cnv_b_muts <- data.frame(CHROM = JK153_cnv_b_muts@fix[,"CHROM"],
                               POS = JK153_cnv_b_muts@fix[,"POS"],
                               b_ref = JK153_cnv_b_muts@fix[,"REF"],
                               b_alt = JK153_cnv_b_muts@fix[,"ALT"],
                               b_qual = as.numeric(JK153_cnv_b_muts@fix[,"QUAL"]),
                               b_ref_dp = sapply(strsplit(JK153_cnv_b_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               b_alt_dp = sapply(strsplit(JK153_cnv_b_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK153_cnv_b_muts$b_dp <- JK153_cnv_b_muts$b_ref_dp + JK153_cnv_b_muts$b_alt_dp


JK153_cnv_c_muts <- read.vcfR("./data/cellranger_dna/JK153_mergebams/clone_c.vcf.gz", verbose = FALSE)
JK153_cnv_c_muts <- data.frame(CHROM = JK153_cnv_c_muts@fix[,"CHROM"],
                               POS = JK153_cnv_c_muts@fix[,"POS"],
                               c_ref = JK153_cnv_c_muts@fix[,"REF"],
                               c_alt = JK153_cnv_c_muts@fix[,"ALT"],
                               c_qual = as.numeric(JK153_cnv_c_muts@fix[,"QUAL"]),
                               c_ref_dp = sapply(strsplit(JK153_cnv_c_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               c_alt_dp = sapply(strsplit(JK153_cnv_c_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK153_cnv_c_muts$c_dp <- JK153_cnv_c_muts$c_ref_dp + JK153_cnv_c_muts$c_alt_dp

JK153_cnv_d_muts <- read.vcfR("./data/cellranger_dna/JK153_mergebams/clone_d.vcf.gz", verbose = FALSE)
JK153_cnv_d_muts <- data.frame(CHROM = JK153_cnv_d_muts@fix[,"CHROM"],
                               POS = JK153_cnv_d_muts@fix[,"POS"],
                               d_ref = JK153_cnv_d_muts@fix[,"REF"],
                               d_alt = JK153_cnv_d_muts@fix[,"ALT"],
                               d_qual = as.numeric(JK153_cnv_d_muts@fix[,"QUAL"]),
                               d_ref_dp = sapply(strsplit(JK153_cnv_d_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                 }),
                               d_alt_dp = sapply(strsplit(JK153_cnv_d_muts@fix[,"INFO"], ";"),
                                                 function(x) {
                                                   dp <- x[grepl("DP4", x)]
                                                   dp_vals <- strsplit(dp, "=")[[1]][2]
                                                   dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                   sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                 }),
                               stringsAsFactors = FALSE)
JK153_cnv_d_muts$d_dp <- JK153_cnv_d_muts$d_ref_dp + JK153_cnv_d_muts$d_alt_dp


JK153_cnv_nonmal_muts <- read.vcfR("./data/cellranger_dna/JK153_mergebams/clone_non-malignant.vcf.gz", verbose = FALSE)
JK153_cnv_nonmal_muts <- data.frame(CHROM = JK153_cnv_nonmal_muts@fix[,"CHROM"],
                                    POS = JK153_cnv_nonmal_muts@fix[,"POS"],
                                    nonmal_ref = JK153_cnv_nonmal_muts@fix[,"REF"],
                                    nonmal_alt = JK153_cnv_nonmal_muts@fix[,"ALT"],
                                    nonmal_qual = as.numeric(JK153_cnv_nonmal_muts@fix[,"QUAL"]),
                                    nonmal_ref_dp = sapply(strsplit(JK153_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[1]), as.numeric(dp_vals[2]))
                                                           }),
                                    nonmal_alt_dp = sapply(strsplit(JK153_cnv_nonmal_muts@fix[,"INFO"], ";"),
                                                           function(x) {
                                                             dp <- x[grepl("DP4", x)]
                                                             dp_vals <- strsplit(dp, "=")[[1]][2]
                                                             dp_vals <- strsplit(dp_vals, ",")[[1]]
                                                             sum(as.numeric(dp_vals[3]), as.numeric(dp_vals[4]))
                                                           }),
                                    stringsAsFactors = FALSE)
JK153_cnv_nonmal_muts$nonmal_dp <- JK153_cnv_nonmal_muts$nonmal_ref_dp + JK153_cnv_nonmal_muts$nonmal_alt_dp


# Merge with loci information
JK153_pyclone_loci <- readRDS("./results/exome/common_mutations/JK153_pyclone_loci.rds") #generated in fig2_trees.R

JK153_cnv_muts_ol <- merge(JK153_pyclone_loci[seq(1, nrow(JK153_pyclone_loci), 6), c("mutation_id", "type", "clonal_cluster", "cell_prev", "gene", "effect", "impact", 
                                                                                     "HVGS.c", "HVGS.p", "CHROM", "POS", "REF", "ALT", "clone_id")], 
                           JK153_cnv_a_muts, by = c("CHROM", "POS"), all.x = TRUE)

# Qualify overlaps (no coverage, low quality [QUAL < 20], flag if indel for manual review, low depth if <5 reads)
JK153_cnv_muts_ol$b_mut <- ifelse(is.na(JK153_cnv_muts_ol$b_dp), "no_cov", 
                                  ifelse(as.numeric(JK153_cnv_muts_ol$b_qual) < 20, "low_qual",
                                         ifelse(JK153_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK153_cnv_muts_ol$b_alt) & JK153_cnv_muts_ol$b_alt == JK153_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK153_cnv_muts_ol$b_dp < 5, "low_dp", "no")))))

JK153_cnv_muts_ol <- merge(JK153_cnv_muts_ol, JK153_cnv_c_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK153_cnv_muts_ol$c_mut <- ifelse(is.na(JK153_cnv_muts_ol$c_dp), "no_cov", 
                                  ifelse(as.numeric(JK153_cnv_muts_ol$c_qual) < 20, "low_qual",
                                         ifelse(JK153_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK153_cnv_muts_ol$c_alt) & JK153_cnv_muts_ol$c_alt == JK153_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK153_cnv_muts_ol$c_dp < 5, "low_dp", "no")))))

JK153_cnv_muts_ol <- merge(JK153_cnv_muts_ol, JK153_cnv_d_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK153_cnv_muts_ol$d_mut <- ifelse(is.na(JK153_cnv_muts_ol$d_dp), "no_cov", 
                                  ifelse(as.numeric(JK153_cnv_muts_ol$d_qual) < 20, "low_qual",
                                         ifelse(JK153_cnv_muts_ol$type == "INDEL", "flag",
                                                ifelse(!is.na(JK153_cnv_muts_ol$d_alt) & JK153_cnv_muts_ol$d_alt == JK153_cnv_muts_ol$ALT, "yes",
                                                       ifelse(JK153_cnv_muts_ol$d_dp < 5, "low_dp", "no")))))

JK153_cnv_muts_ol <- merge(JK153_cnv_muts_ol, JK153_cnv_nonmal_muts, by = c("CHROM", "POS"), all.x = TRUE)
JK153_cnv_muts_ol$nonmal_mut <- ifelse(is.na(JK153_cnv_muts_ol$nonmal_dp), "no_cov", 
                                       ifelse(as.numeric(JK153_cnv_muts_ol$nonmal_qual) < 20, "low_qual",
                                              ifelse(JK153_cnv_muts_ol$type == "INDEL", "flag",
                                                     ifelse(!is.na(JK153_cnv_muts_ol$nonmal_alt) & 
                                                              JK153_cnv_muts_ol$nonmal_alt == JK153_cnv_muts_ol$ALT, "yes",
                                                            ifelse(JK153_cnv_muts_ol$nonmal_dp < 5, "low_dp", "no")))))


#Add mean cellular prevalence for each variant
JK153_cnv_muts_ol <- merge(JK153_cnv_muts_ol, aggregate(JK153_pyclone_loci$cellular_prevalence, by = list(JK153_pyclone_loci$mutation_id), mean),
                           by.x = "mutation_id", by.y = "Group.1")
colnames(JK153_cnv_muts_ol)[colnames(JK153_cnv_muts_ol) == "x"] <- "mean_cell_prev"

write.csv(JK153_cnv_muts_ol, file = "./results/exome/cnv_integration/JK153_cnv_clone_mutations.csv", row.names = FALSE)
#fix indels and duplicates (usually near reads with indels) manually

# Re-load cleaned up data
JK153_cnv_muts_ol <- read.csv("./results/exome/cnv_integration/JK153_cnv_clone_mutations.csv", stringsAsFactors = FALSE)
JK153_cnv_muts_ol$clone_id <- factor(JK153_cnv_muts_ol$clone_id)
JK153_cnv_muts_ol[,grepl("_mut", colnames(JK153_cnv_muts_ol))] <- lapply(JK153_cnv_muts_ol[,grepl("_mut", colnames(JK153_cnv_muts_ol))],
                                                                         function(x) factor(x, levels = c("yes", "no", "low_dp", "low_qual", "no_cov")))

# View overlaps by clone (coloured in Figure 2 if >10% of mutation clone variants are found in copy number clone)
setNames(lapply(colnames(JK153_cnv_muts_ol)[grepl("_mut", colnames(JK153_cnv_muts_ol))], function(cnv_clone) {
  sapply(levels(JK153_cnv_muts_ol$clone_id), function(clone) table(JK153_cnv_muts_ol[which(JK153_cnv_muts_ol$clone_id == clone), cnv_clone]) / 
           length(which(JK153_cnv_muts_ol$clone_id == clone)))
}), colnames(JK153_cnv_muts_ol)[grepl("_mut", colnames(JK153_cnv_muts_ol))])