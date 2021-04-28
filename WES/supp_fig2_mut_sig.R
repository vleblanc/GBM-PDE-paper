# ---------------------------------------------------------------
## Scripts to run mutation signature analyses
## Author: Veronique LeBlanc
# ---------------------------------------------------------------

source("./funcs_analysis.R")

library(SignIT)
library(ggraph) #needed for some of the SignIT plots
library(BSgenome.Hsapiens.UCSC.hg38)

# ---------------------------------------------------------------
## JK124
# ---------------------------------------------------------------

# Load variants
JK124_vars <- get.var.info.from.table("./results/exome/common_mutations/JK124_noreg2org_strelka_mutect/tp-baseline.table")

# Split per sample and keep variants with at least 20 reads at least 10 alt reads
JK124_samp_vars <- setNames(lapply(c("JK124_reg1_tissue", "JK124_reg1_organoid", "JK124_reg2_tissue"),
                                   function(samp) {
                                     samp_vars <- JK124_vars$tbl_nosnp[JK124_vars$tbl_nosnp[,paste0(samp, ".DP")] >= 20 &
                                                                         JK124_vars$tbl_nosnp[,paste0(samp, "_alt_ad")] >= 10 &
                                                                         JK124_vars$tbl_nosnp$TYPE == "SNP",]
                                     samp_vars_df <- data.frame(chr = samp_vars$CHROM, pos = samp_vars$POS,
                                                                ref = samp_vars$REF, alt = samp_vars$ALT)
                                   }), c("JK124_reg1_tissue", "JK124_reg1_organoid", "JK124_reg2_tissue"))

# REG1 TISSUE

# Convert to catalog
JK124_reg1_tis_signit_cat <- mutations_to_catalog(chr = JK124_samp_vars$JK124_reg1_tissue$chr,
                                                  pos = JK124_samp_vars$JK124_reg1_tissue$pos,
                                                  ref = JK124_samp_vars$JK124_reg1_tissue$ref,
                                                  alt = JK124_samp_vars$JK124_reg1_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK124_reg1_tis_signit <- get_exposures(JK124_reg1_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK124_reg1_tis_signit, units = "fraction")


# REG1 PDO

# Convert to catalog
JK124_reg1_org_signit_cat <- mutations_to_catalog(chr = JK124_samp_vars$JK124_reg1_organoid$chr,
                                                  pos = JK124_samp_vars$JK124_reg1_organoid$pos,
                                                  ref = JK124_samp_vars$JK124_reg1_organoid$ref,
                                                  alt = JK124_samp_vars$JK124_reg1_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK124_reg1_org_signit <- get_exposures(JK124_reg1_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK124_reg1_org_signit, units = "fraction")


# REG2 TISSUE

# Convert to catalog
JK124_reg2_tis_signit_cat <- mutations_to_catalog(chr = JK124_samp_vars$JK124_reg2_tissue$chr,
                                                  pos = JK124_samp_vars$JK124_reg2_tissue$pos,
                                                  ref = JK124_samp_vars$JK124_reg2_tissue$ref,
                                                  alt = JK124_samp_vars$JK124_reg2_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK124_reg2_tis_signit <- get_exposures(JK124_reg2_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK124_reg2_tis_signit, units = "fraction")






# ---------------------------------------------------------------
## JK136
# ---------------------------------------------------------------

# Load variants
JK136_vars <- get.var.info.from.table("./results/exome/common_mutations/JK136_JK202_strelka_mutect/tp-baseline.table")

# Split per sample and keep variants with at least 20 reads at least 10 alt reads
JK136_samp_vars <- setNames(lapply(c("JK136_reg1_tissue", "JK136_reg1_organoid", "JK136_reg2_tissue", "JK136_reg2_organoid",
                                     "JK202_tissue", "JK202_organoid"),
                                   function(samp) {
                                     samp_vars <- JK136_vars$tbl_nosnp[JK136_vars$tbl_nosnp[,paste0(samp, ".DP")] >= 20 &
                                                                         JK136_vars$tbl_nosnp[,paste0(samp, "_alt_ad")] >= 10 &
                                                                         JK136_vars$tbl_nosnp$TYPE == "SNP",]
                                     samp_vars_df <- data.frame(chr = samp_vars$CHROM, pos = samp_vars$POS,
                                                                ref = samp_vars$REF, alt = samp_vars$ALT)
                                   }), c("JK136_reg1_tissue", "JK136_reg1_organoid", "JK136_reg2_tissue", "JK136_reg2_organoid",
                                         "JK202_tissue", "JK202_organoid"))

# REG1 TISSUE

# Convert to catalog
JK136_reg1_tis_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK136_reg1_tissue$chr,
                                                  pos = JK136_samp_vars$JK136_reg1_tissue$pos,
                                                  ref = JK136_samp_vars$JK136_reg1_tissue$ref,
                                                  alt = JK136_samp_vars$JK136_reg1_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK136_reg1_tis_signit <- get_exposures(JK136_reg1_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_reg1_tis_signit, units = "fraction")


# REG1 PDO

# Convert to catalog
JK136_reg1_org_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK136_reg1_organoid$chr,
                                                  pos = JK136_samp_vars$JK136_reg1_organoid$pos,
                                                  ref = JK136_samp_vars$JK136_reg1_organoid$ref,
                                                  alt = JK136_samp_vars$JK136_reg1_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK136_reg1_org_signit <- get_exposures(JK136_reg1_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_reg1_org_signit, units = "fraction")


# REG2 TISSUE

# Convert to catalog
JK136_reg2_tis_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK136_reg2_tissue$chr,
                                                  pos = JK136_samp_vars$JK136_reg2_tissue$pos,
                                                  ref = JK136_samp_vars$JK136_reg2_tissue$ref,
                                                  alt = JK136_samp_vars$JK136_reg2_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK136_reg2_tis_signit <- get_exposures(JK136_reg2_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_reg2_tis_signit, units = "fraction")


# REG2 PDO

# Convert to catalog
JK136_reg2_org_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK136_reg2_organoid$chr,
                                                  pos = JK136_samp_vars$JK136_reg2_organoid$pos,
                                                  ref = JK136_samp_vars$JK136_reg2_organoid$ref,
                                                  alt = JK136_samp_vars$JK136_reg2_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK136_reg2_org_signit <- get_exposures(JK136_reg2_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_reg2_org_signit, units = "fraction")


# RECURRENT TISSUE

# Convert to catalog
JK136_rec_tis_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK202_tissue$chr,
                                                 pos = JK136_samp_vars$JK202_tissue$pos,
                                                 ref = JK136_samp_vars$JK202_tissue$ref,
                                                 alt = JK136_samp_vars$JK202_tissue$alt,
                                                 genome = BSgenome.Hsapiens.UCSC.hg38)


# Calculate exposures 
JK136_rec_tis_signit <- get_exposures(JK136_rec_tis_signit_cat, 
                                      reference_signatures = cosmic_v3, 
                                      n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_rec_tis_signit, unit = "fraction")


# RECURRENT PDO

# Convert to catalog
JK136_rec_org_signit_cat <- mutations_to_catalog(chr = JK136_samp_vars$JK202_organoid$chr,
                                                 pos = JK136_samp_vars$JK202_organoid$pos,
                                                 ref = JK136_samp_vars$JK202_organoid$ref,
                                                 alt = JK136_samp_vars$JK202_organoid$alt,
                                                 genome = BSgenome.Hsapiens.UCSC.hg38)


# Calculate exposures 
JK136_rec_org_signit <- get_exposures(JK136_rec_org_signit_cat, 
                                      reference_signatures = cosmic_v3, 
                                      n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK136_rec_org_signit, unit = "fraction")






# ---------------------------------------------------------------
## JK142
# ---------------------------------------------------------------

# Load variants
JK142_vars <- get.var.info.from.table("./results/exome/common_mutations/JK142_JK196_noreg2tis_strelka_mutect/tp-baseline.table")

# Split per sample and keep variants with at least 20 reads at least 10 alt reads
JK142_samp_vars <- setNames(lapply(c("JK142_reg1_tissue", "JK142_reg1_organoid", "JK142_reg2_organoid",
                                     "JK196_tissue", "JK196_organoid"),
                                   function(samp) {
                                     samp_vars <- JK142_vars$tbl_nosnp[JK142_vars$tbl_nosnp[,paste0(samp, ".DP")] >= 20 &
                                                                         JK142_vars$tbl_nosnp[,paste0(samp, "_alt_ad")] >= 10 &
                                                                         JK142_vars$tbl_nosnp$TYPE == "SNP",]
                                     samp_vars_df <- data.frame(chr = samp_vars$CHROM, pos = samp_vars$POS,
                                                                ref = samp_vars$REF, alt = samp_vars$ALT)
                                   }), c("JK142_reg1_tissue", "JK142_reg1_organoid", "JK142_reg2_organoid", "JK196_tissue", "JK196_organoid"))

# REG1 TISSUE

# Convert to catalog
JK142_reg1_tis_signit_cat <- mutations_to_catalog(chr = JK142_samp_vars$JK142_reg1_tissue$chr,
                                                  pos = JK142_samp_vars$JK142_reg1_tissue$pos,
                                                  ref = JK142_samp_vars$JK142_reg1_tissue$ref,
                                                  alt = JK142_samp_vars$JK142_reg1_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK142_reg1_tis_signit <- get_exposures(JK142_reg1_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK142_reg1_tis_signit, units = "fraction")


# REG1 PDO

# Convert to catalog
JK142_reg1_org_signit_cat <- mutations_to_catalog(chr = JK142_samp_vars$JK142_reg1_organoid$chr,
                                                  pos = JK142_samp_vars$JK142_reg1_organoid$pos,
                                                  ref = JK142_samp_vars$JK142_reg1_organoid$ref,
                                                  alt = JK142_samp_vars$JK142_reg1_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK142_reg1_org_signit <- get_exposures(JK142_reg1_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK142_reg1_org_signit, units = "fraction")


# REG2 PDO

# Convert to catalog
JK142_reg2_org_signit_cat <- mutations_to_catalog(chr = JK142_samp_vars$JK142_reg2_organoid$chr,
                                                  pos = JK142_samp_vars$JK142_reg2_organoid$pos,
                                                  ref = JK142_samp_vars$JK142_reg2_organoid$ref,
                                                  alt = JK142_samp_vars$JK142_reg2_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK142_reg2_org_signit <- get_exposures(JK142_reg2_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK142_reg2_org_signit, units = "fraction")


# RECURRENT TISSUE

# Convert to catalog
JK142_rec_tis_signit_cat <- mutations_to_catalog(chr = JK142_samp_vars$JK196_tissue$chr,
                                                 pos = JK142_samp_vars$JK196_tissue$pos,
                                                 ref = JK142_samp_vars$JK196_tissue$ref,
                                                 alt = JK142_samp_vars$JK196_tissue$alt,
                                                 genome = BSgenome.Hsapiens.UCSC.hg38)


# Calculate exposures 
JK142_rec_tis_signit <- get_exposures(JK142_rec_tis_signit_cat, 
                                      reference_signatures = cosmic_v3, 
                                      n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK142_rec_tis_signit, unit = "fraction")


# RECURRENT PDO

# Convert to catalog
JK142_rec_org_signit_cat <- mutations_to_catalog(chr = JK142_samp_vars$JK196_organoid$chr,
                                                 pos = JK142_samp_vars$JK196_organoid$pos,
                                                 ref = JK142_samp_vars$JK196_organoid$ref,
                                                 alt = JK142_samp_vars$JK196_organoid$alt,
                                                 genome = BSgenome.Hsapiens.UCSC.hg38)


# Calculate exposures 
JK142_rec_org_signit <- get_exposures(JK142_rec_org_signit_cat, 
                                      reference_signatures = cosmic_v3, 
                                      n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK142_rec_org_signit, unit = "fraction")








# ---------------------------------------------------------------
## JK153
# ---------------------------------------------------------------

# Load variants
JK153_vars <- get.var.info.from.table("./results/exome/common_mutations/JK153_strelka_mutect/tp-baseline.table")

# Split per sample and keep variants with at least 20 reads at least 10 alt reads
JK153_samp_vars <- setNames(lapply(c("JK153_reg1_tissue", "JK153_reg1_organoid", "JK153_reg2_tissue", "JK153_reg2_organoid"),
                                   function(samp) {
                                     samp_vars <- JK153_vars$tbl_nosnp[JK153_vars$tbl_nosnp[,paste0(samp, ".DP")] >= 20 &
                                                                         JK153_vars$tbl_nosnp[,paste0(samp, "_alt_ad")] >= 10 &
                                                                         JK153_vars$tbl_nosnp$TYPE == "SNP",]
                                     samp_vars_df <- data.frame(chr = samp_vars$CHROM, pos = samp_vars$POS,
                                                                ref = samp_vars$REF, alt = samp_vars$ALT)
                                   }), c("JK153_reg1_tissue", "JK153_reg1_organoid", "JK153_reg2_tissue", "JK153_reg2_organoid"))

# REG1 TISSUE

# Convert to catalog
JK153_reg1_tis_signit_cat <- mutations_to_catalog(chr = JK153_samp_vars$JK153_reg1_tissue$chr,
                                                  pos = JK153_samp_vars$JK153_reg1_tissue$pos,
                                                  ref = JK153_samp_vars$JK153_reg1_tissue$ref,
                                                  alt = JK153_samp_vars$JK153_reg1_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK153_reg1_tis_signit <- get_exposures(JK153_reg1_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK153_reg1_tis_signit, units = "fraction")


# REG1 PDO

# Convert to catalog
JK153_reg1_org_signit_cat <- mutations_to_catalog(chr = JK153_samp_vars$JK153_reg1_organoid$chr,
                                                  pos = JK153_samp_vars$JK153_reg1_organoid$pos,
                                                  ref = JK153_samp_vars$JK153_reg1_organoid$ref,
                                                  alt = JK153_samp_vars$JK153_reg1_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK153_reg1_org_signit <- get_exposures(JK153_reg1_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK153_reg1_org_signit, units = "fraction")


# REG2 TISSUE

# Convert to catalog
JK153_reg2_tis_signit_cat <- mutations_to_catalog(chr = JK153_samp_vars$JK153_reg2_tissue$chr,
                                                  pos = JK153_samp_vars$JK153_reg2_tissue$pos,
                                                  ref = JK153_samp_vars$JK153_reg2_tissue$ref,
                                                  alt = JK153_samp_vars$JK153_reg2_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK153_reg2_tis_signit <- get_exposures(JK153_reg2_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK153_reg2_tis_signit, units = "fraction")


# REG2 PDO

# Convert to catalog
JK153_reg2_org_signit_cat <- mutations_to_catalog(chr = JK153_samp_vars$JK153_reg2_organoid$chr,
                                                  pos = JK153_samp_vars$JK153_reg2_organoid$pos,
                                                  ref = JK153_samp_vars$JK153_reg2_organoid$ref,
                                                  alt = JK153_samp_vars$JK153_reg2_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK153_reg2_org_signit <- get_exposures(JK153_reg2_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK153_reg2_org_signit, units = "fraction")







# ---------------------------------------------------------------
## JK163
# ---------------------------------------------------------------

# Load variants
JK163_vars <- get.var.info.from.table("./results/exome/common_mutations/JK163_strelka_mutect/tp-baseline.table")

# Split per sample and keep variants with at least 20 reads at least 10 alt reads
JK163_samp_vars <- setNames(lapply(c("JK163_reg1_tissue", "JK163_reg1_organoid", "JK163_reg2_tissue", "JK163_reg2_organoid"),
                                   function(samp) {
                                     samp_vars <- JK163_vars$tbl_nosnp[JK163_vars$tbl_nosnp[,paste0(samp, ".DP")] >= 20 &
                                                                         JK163_vars$tbl_nosnp[,paste0(samp, "_alt_ad")] >= 10 &
                                                                         JK163_vars$tbl_nosnp$TYPE == "SNP",]
                                     samp_vars_df <- data.frame(chr = samp_vars$CHROM, pos = samp_vars$POS,
                                                                ref = samp_vars$REF, alt = samp_vars$ALT)
                                   }), c("JK163_reg1_tissue", "JK163_reg1_organoid", "JK163_reg2_tissue", "JK163_reg2_organoid"))

# REG1 TISSUE

# Convert to catalog
JK163_reg1_tis_signit_cat <- mutations_to_catalog(chr = JK163_samp_vars$JK163_reg1_tissue$chr,
                                                  pos = JK163_samp_vars$JK163_reg1_tissue$pos,
                                                  ref = JK163_samp_vars$JK163_reg1_tissue$ref,
                                                  alt = JK163_samp_vars$JK163_reg1_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK163_reg1_tis_signit <- get_exposures(JK163_reg1_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK163_reg1_tis_signit, units = "fraction")


# REG1 PDO

# Convert to catalog
JK163_reg1_org_signit_cat <- mutations_to_catalog(chr = JK163_samp_vars$JK163_reg1_organoid$chr,
                                                  pos = JK163_samp_vars$JK163_reg1_organoid$pos,
                                                  ref = JK163_samp_vars$JK163_reg1_organoid$ref,
                                                  alt = JK163_samp_vars$JK163_reg1_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK163_reg1_org_signit <- get_exposures(JK163_reg1_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK163_reg1_org_signit, units = "fraction")


# REG2 TISSUE

# Convert to catalog
JK163_reg2_tis_signit_cat <- mutations_to_catalog(chr = JK163_samp_vars$JK163_reg2_tissue$chr,
                                                  pos = JK163_samp_vars$JK163_reg2_tissue$pos,
                                                  ref = JK163_samp_vars$JK163_reg2_tissue$ref,
                                                  alt = JK163_samp_vars$JK163_reg2_tissue$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK163_reg2_tis_signit <- get_exposures(JK163_reg2_tis_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK163_reg2_tis_signit, units = "fraction")


# REG2 PDO

# Convert to catalog
JK163_reg2_org_signit_cat <- mutations_to_catalog(chr = JK163_samp_vars$JK163_reg2_organoid$chr,
                                                  pos = JK163_samp_vars$JK163_reg2_organoid$pos,
                                                  ref = JK163_samp_vars$JK163_reg2_organoid$ref,
                                                  alt = JK163_samp_vars$JK163_reg2_organoid$alt,
                                                  genome = BSgenome.Hsapiens.UCSC.hg38)

# Calculate exposures 
JK163_reg2_org_signit <- get_exposures(JK163_reg2_org_signit_cat, 
                                       reference_signatures = cosmic_v3, 
                                       n_iter = 1000, n_cores = 10)

# Plot
plot_exposure_posteriors(JK163_reg2_org_signit, units = "fraction")







# ---------------------------------------------------------------
## All samples
# ---------------------------------------------------------------

all_signit <- data.frame(signature = levels(JK124_reg1_tis_signit$exposure_chain$signature),
                         JK124_reg1_tis = aggregate(JK124_reg1_tis_signit$exposure_chain$exposure / JK124_reg1_tis_signit$n_mutations, 
                                                    by = list(JK124_reg1_tis_signit$exposure_chain$signature), median)$x,
                         JK124_reg1_org = aggregate(JK124_reg1_org_signit$exposure_chain$exposure / JK124_reg1_org_signit$n_mutations, 
                                                    by = list(JK124_reg1_org_signit$exposure_chain$signature), median)$x,
                         JK124_reg2_tis = aggregate(JK124_reg2_tis_signit$exposure_chain$exposure / JK124_reg2_tis_signit$n_mutations, 
                                                    by = list(JK136_reg2_tis_signit$exposure_chain$signature), median)$x,
                         JK136_reg1_tis = aggregate(JK136_reg1_tis_signit$exposure_chain$exposure / JK136_reg1_tis_signit$n_mutations, 
                                                    by = list(JK136_reg1_tis_signit$exposure_chain$signature), median)$x,
                         JK136_reg1_org = aggregate(JK136_reg1_org_signit$exposure_chain$exposure / JK136_reg1_org_signit$n_mutations, 
                                                    by = list(JK136_reg1_org_signit$exposure_chain$signature), median)$x,
                         JK136_reg2_tis = aggregate(JK136_reg2_tis_signit$exposure_chain$exposure / JK136_reg2_tis_signit$n_mutations, 
                                                    by = list(JK136_reg2_tis_signit$exposure_chain$signature), median)$x,
                         JK136_reg2_org = aggregate(JK136_reg2_org_signit$exposure_chain$exposure / JK136_reg2_org_signit$n_mutations, 
                                                    by = list(JK136_reg2_org_signit$exposure_chain$signature), median)$x,
                         JK136_rec_tis = aggregate(JK136_rec_tis_signit$exposure_chain$exposure / JK136_rec_tis_signit$n_mutations, 
                                                   by = list(JK136_rec_tis_signit$exposure_chain$signature), median)$x,
                         JK136_rec_org = aggregate(JK136_rec_org_signit$exposure_chain$exposure / JK136_rec_org_signit$n_mutations, 
                                                   by = list(JK136_rec_org_signit$exposure_chain$signature), median)$x,
                         JK142_reg1_tis = aggregate(JK142_reg1_tis_signit$exposure_chain$exposure / JK142_reg1_tis_signit$n_mutations, 
                                                    by = list(JK142_reg1_tis_signit$exposure_chain$signature), median)$x,
                         JK142_reg1_org = aggregate(JK142_reg1_org_signit$exposure_chain$exposure / JK142_reg1_org_signit$n_mutations, 
                                                    by = list(JK142_reg1_org_signit$exposure_chain$signature), median)$x,
                         JK142_reg2_org = aggregate(JK142_reg2_org_signit$exposure_chain$exposure / JK142_reg2_org_signit$n_mutations, 
                                                    by = list(JK142_reg2_org_signit$exposure_chain$signature), median)$x,
                         JK142_rec_tis = aggregate(JK142_rec_tis_signit$exposure_chain$exposure / JK142_rec_tis_signit$n_mutations, 
                                                   by = list(JK142_rec_tis_signit$exposure_chain$signature), median)$x,
                         JK142_rec_org = aggregate(JK142_rec_org_signit$exposure_chain$exposure / JK142_rec_org_signit$n_mutations, 
                                                   by = list(JK142_rec_org_signit$exposure_chain$signature), median)$x,
                         JK153_reg1_tis = aggregate(JK153_reg1_tis_signit$exposure_chain$exposure / JK153_reg1_tis_signit$n_mutations, 
                                                    by = list(JK153_reg1_tis_signit$exposure_chain$signature), median)$x,
                         JK153_reg1_org = aggregate(JK153_reg1_org_signit$exposure_chain$exposure / JK153_reg1_org_signit$n_mutations, 
                                                    by = list(JK153_reg1_org_signit$exposure_chain$signature), median)$x,
                         JK153_reg2_tis = aggregate(JK153_reg2_tis_signit$exposure_chain$exposure / JK153_reg2_tis_signit$n_mutations, 
                                                    by = list(JK153_reg2_tis_signit$exposure_chain$signature), median)$x,
                         JK153_reg2_org = aggregate(JK153_reg2_org_signit$exposure_chain$exposure / JK153_reg2_org_signit$n_mutations, 
                                                    by = list(JK153_reg2_org_signit$exposure_chain$signature), median)$x,
                         JK163_reg1_tis = aggregate(JK163_reg1_tis_signit$exposure_chain$exposure / JK163_reg1_tis_signit$n_mutations, 
                                                    by = list(JK163_reg1_tis_signit$exposure_chain$signature), median)$x,
                         JK163_reg1_org = aggregate(JK163_reg1_org_signit$exposure_chain$exposure / JK163_reg1_org_signit$n_mutations, 
                                                    by = list(JK163_reg1_org_signit$exposure_chain$signature), median)$x,
                         JK163_reg2_tis = aggregate(JK163_reg2_tis_signit$exposure_chain$exposure / JK163_reg2_tis_signit$n_mutations, 
                                                    by = list(JK163_reg2_tis_signit$exposure_chain$signature), median)$x,
                         JK163_reg2_org = aggregate(JK163_reg2_org_signit$exposure_chain$exposure / JK163_reg2_org_signit$n_mutations, 
                                                    by = list(JK163_reg2_org_signit$exposure_chain$signature), median)$x)

WriteXLS::WriteXLS(all_signit, "./results/exome/all_samps_signit_median_fractions.xlsx", row.names = FALSE) #Supplementary Table 3
